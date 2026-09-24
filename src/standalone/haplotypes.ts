/**
 * The two haplotypes of a reads window as consensus sequences, from either source of phasing:
 *
 * - **the file's haplotags**: phasing tools write, on every read they can place, the haplotype (HP:i:1 / 2) and
 *   the phase set it belongs to (PS:i), sometimes a Phred-scaled confidence (PC:i): WhatsHap and LongPhase
 *   `haplotag` (the ONT wf-human-variation outputs), PacBio HiPhase, Illumina DRAGEN (TruPath, where HP is a
 *   "copy label"). That phasing comes from the whole genome's variants and reaches across what a window's reads
 *   cannot link;
 * - **the in-page read-based phasing** (phasing.ts): each fragment goes to the haplotype of its block it matches
 *   at more sites; the blocks play the part of the phase sets.
 *
 * Each (phase set, haplotype) then gets its consensus: callSites on its reads with the allele fraction taken
 * within the haplotype and at least one half, so a heterozygous variant lands on one haplotype and a homozygous
 * one on both; deletions of any length included, the long-read minimum indel size applied. Covered stretches
 * count deletions as covered (the haplotype is known there: it lacks those bases).
 *
 * Checks: at the heterozygous sites of the window (allele fraction 25–75 % over all reads), the two haplotypes of
 * a set should split the alleles (at least SPLIT_MAJOR of one carrying the alt, at most 1 − SPLIT_MAJOR of the
 * other). Sites that do not split are listed (a mis-phased or mosaic site, a third haplotype, a collapsed
 * duplication), and reads whose allele contradicts their haplotype's are counted (tagging errors, chimeras).
 */
import type { AlignedRead, AllelicJunction, HapJunction, HaplotypeConsensus, HaplotypeSet, HaplotypeView, PhaseResult, VariantSite } from '../components/sashimi/types';
import { fisher } from '../components/sashimi/siteQuality';
import { callSites, jKey, junctionSnap, readJunctions } from './collapse';
import { HET_MAX, HET_MIN, alleleAt, fragmentsOf } from './phasing';

/** Reads of one haplotype needed for a stretch to count as covered, and for a site to be judged. */
export const HAP_MIN_DEPTH = 3;
/** Share of a haplotype's reads that carry an allele for it to be in the consensus. */
const CONSENSUS_FRACTION = 0.5;
/** A heterozygous site is split by the haplotypes when one carries the alt in at least this share of its reads and the other in at most 1 − it. */
const SPLIT_MAJOR = 0.8;

/** Reads of the window carrying a haplotag, and the phase sets among them. */
export function haplotagCounts(reads: AlignedRead[]): { tagged: number; sets: number } {
  let tagged = 0;
  const sets = new Set<number>();
  for (const r of reads) if (r.hp) { tagged++; if (r.ps != null) sets.add(r.ps); }
  return { tagged, sets: sets.size };
}

interface Group { key: string; ps: number | null; hap: number; frags: AlignedRead[][] }
interface Options { start: number; end: number; ref: string | null; refStart: number; minBq: number; minIndel: number; minAlt: number; minVaf: number; /** the window's sites over all reads, when the caller has them */ sites?: VariantSite[]; /** long reads: junctions a few bases off taken for the common one (collapse.ts, junctionSnap) */ longReads?: boolean }

/** Haplotypes from the HP / PS tags: one group per (phase set, haplotype); a fragment takes the tag of the mate that has one. */
export function haplotypesFromTags(reads: AlignedRead[], o: Options): HaplotypeView {
  const groups = new Map<string, Group>();
  let unassigned = 0, assigned = 0;
  for (const fr of fragmentsOf(reads)) {
    const tagged = fr.find(r => r.hp);
    if (!tagged) { unassigned += fr.length; continue; }
    const ps = tagged.ps ?? null, key = `${ps ?? 'none'}:${tagged.hp}`;
    let g = groups.get(key);
    if (!g) groups.set(key, g = { key, ps, hap: tagged.hp!, frags: [] });
    g.frags.push(fr);
    assigned += fr.length;
  }
  const bySet = new Map<string, Group[]>();
  for (const g of groups.values()) { const k = `${g.ps ?? 'none'}`; const l = bySet.get(k); if (l) l.push(g); else bySet.set(k, [g]); }
  const sets: { id: string; ps: number | null; groups: Group[] }[] = [...bySet.entries()].map(([k, gs]) => ({
    id: k === 'none' ? 'no phase set' : `PS ${Number(k).toLocaleString('en-US')}`, ps: gs[0].ps, groups: gs.sort((a, b) => a.hap - b.hap),
  }));
  return build('tags', sets, reads, assigned, unassigned, o, true);
}

/** Haplotypes from the in-page phasing: each fragment to the haplotype of the block it matches at more sites (ties unassigned). */
export function haplotypesFromPhase(reads: AlignedRead[], phase: PhaseResult, o: Options): HaplotypeView {
  const frags = fragmentsOf(reads);
  const sets = phase.blocks.map(b => ({ id: b.id, ps: null as number | null, groups: [1, 2].map(hap => ({ key: `${b.id}:${hap}`, ps: null, hap, frags: [] as AlignedRead[][] })) }));
  // the phased sites in position order (block, index in it): each fragment looks only at those inside its span, not
  // at every site of every block
  const phased = phase.blocks.flatMap((b, bi) => b.sites.map((si, k) => ({ pos: phase.sites[si].pos, bi, k, s: phase.sites[si] }))).sort((x, y) => x.pos - y.pos);
  const firstAt = (pos: number) => { let lo = 0, hi = phased.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (phased[mid].pos < pos) lo = mid + 1; else hi = mid; } return lo; };
  let assigned = 0, unassigned = 0;
  for (const fr of frags) {
    let best = -1, bestMargin = 0, bestHap = 0;
    let lo = Infinity, hi = -Infinity;
    for (const r of fr) { if (r.s < lo) lo = r.s; if (r.e > hi) hi = r.e; }
    const tally = new Map<number, [number, number]>();
    for (let x = firstAt(lo); x < phased.length && phased[x].pos < hi; x++) {
      const { bi, k, s } = phased[x];
      let a = -1;
      for (const r of fr) { const y = alleleAt(r, s, o.minBq); if (y < 0) continue; if (a < 0) a = y; else if (a !== y) { a = -1; break; } }
      if (a < 0) continue;
      let t = tally.get(bi);
      if (!t) tally.set(bi, t = [0, 0]);
      if ((a === 1) === (phase.blocks[bi].h1[k] === 'alt')) t[0]++; else t[1]++;
    }
    // the block where the fragment is the most clearly on one side (the first one on a tie, as before)
    for (const [bi, [m1, m2]] of [...tally].sort((x, y) => x[0] - y[0])) {
      const margin = Math.abs(m1 - m2);
      if (margin > bestMargin) { best = bi; bestMargin = margin; bestHap = m1 > m2 ? 1 : 2; }
    }
    if (best < 0) { unassigned += fr.length; continue; }
    sets[best].groups[bestHap - 1].frags.push(fr);
    assigned += fr.length;
  }
  return build('reads', sets, reads, assigned, unassigned, o, false);
}

function build(source: HaplotypeView['source'], sets: { id: string; ps: number | null; groups: Group[] }[], reads: AlignedRead[],
  assigned: number, unassigned: number, o: Options, countReads: boolean): HaplotypeView {
  const out: HaplotypeSet[] = [];
  const notSplit: HaplotypeView['notSplit'] = [];
  let checked = 0;
  const conflicting = new Set<AlignedRead>();
  // heterozygous sites of the window, over all reads, with the sample's thresholds
  const all = o.sites ?? callSites(reads, o.start, o.end, o.ref, o.refStart, o.minAlt, o.minVaf, o.minBq, o.minIndel);
  const het = all.filter(s => s.vaf >= HET_MIN && s.vaf <= HET_MAX).sort((a, b) => a.pos - b.pos);
  const hetFirst = (pos: number) => { let lo = 0, hi = het.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (het[mid].pos < pos) lo = mid + 1; else hi = mid; } return lo; };
  const windowKeys = new Set(all.map(siteKey));
  const snap = o.longReads ? junctionSnap(reads) : null;
  for (const set of sets) {
    const setReads = set.groups.flatMap(g => g.frags.flat());
    if (!setReads.length) continue;
    let lo = Infinity, hi = -Infinity;
    for (const r of setReads) { if (r.s < lo) lo = r.s; if (r.e > hi) hi = r.e; }
    const s0 = Math.max(o.start, lo), e0 = Math.min(o.end, hi);
    if (e0 <= s0) continue;
    const haps: HaplotypeConsensus[] = set.groups.filter(g => g.frags.length).map(g => {
      const rs = g.frags.flat();
      const covered = coveredStretches(rs, s0, e0);
      // a variant counts only where the haplotype is covered by HAP_MIN_DEPTH reads, as the row is drawn
      const inCovered = (p: number) => covered.some(([a, b]) => a <= p && p < b);
      // a site of one haplotype must be a site of the window over all its reads too: a real allele always is, while two
      // sequencing errors of the same base among the 3–4 reads of a haplotype at the window's edge otherwise make one
      const sites = callSites(rs, s0, e0, o.ref, o.refStart, Math.min(o.minAlt, 2), CONSENSUS_FRACTION, o.minBq, o.minIndel)
        .filter(x => inCovered(x.pos) && windowKeys.has(siteKey(x)));
      const pcs = rs.map(r => r.pc).filter((x): x is number => x != null).sort((a, b) => a - b);
      return { hap: g.hap, reads: countReads ? rs.length : g.frags.length, covered, sites, ...(pcs.length ? { pc: pcs[pcs.length >> 1] } : {}) };
    });
    const withReads = set.groups.filter(g => g.frags.length);
    const { per, allelic } = haplotypeJunctions(withReads, snap);
    per.forEach((js, k) => { if (js.length) haps[k].junctions = js; });
    out.push({ id: set.id, ps: set.ps, start: s0, end: e0, haps, ...(allelic.length ? { allelic } : {}) });
    // the split check needs two haplotypes
    const [g1, g2] = set.groups;
    if (!g1 || !g2 || !g1.frags.length || !g2.frags.length) continue;
    // each fragment's alleles at the heterozygous sites inside its span (and the set's): not every site × every fragment
    const seen = [g1, g2].map(g => g.frags.map(fr => {
      let lo = Infinity, hi = -Infinity;
      for (const r of fr) { if (r.s < lo) lo = r.s; if (r.e > hi) hi = r.e; }
      const out: [number, number][] = [];
      for (let x = hetFirst(Math.max(lo, s0)); x < het.length && het[x].pos < Math.min(hi, e0); x++) { const a = fragAllele(fr, het[x], o.minBq); if (a >= 0) out.push([x, a]); }
      return out;
    }));
    const counts = seen.map(fs => { const m = new Map<number, [number, number]>(); for (const f of fs) for (const [x, a] of f) { const c = m.get(x) ?? [0, 0]; c[1]++; if (a === 1) c[0]++; m.set(x, c); } return m; });
    const altOnSite = new Map<number, number>();
    for (const x of [...counts[0].keys()].sort((a, b) => a - b)) {
      const c = [counts[0].get(x), counts[1].get(x)];
      if (!c[0] || !c[1] || c[0][1] < HAP_MIN_DEPTH || c[1][1] < HAP_MIN_DEPTH) continue;
      const frac = c.map(t => t![0] / t![1]);
      checked++;
      const altOn = frac[0] >= SPLIT_MAJOR && frac[1] <= 1 - SPLIT_MAJOR ? 0 : frac[1] >= SPLIT_MAJOR && frac[0] <= 1 - SPLIT_MAJOR ? 1 : -1;
      const s = het[x];
      if (altOn < 0) notSplit.push({ pos: s.pos, kind: s.kind, alt: s.alt, fractions: frac, set: set.id });
      else altOnSite.set(x, altOn);
    }
    [g1, g2].forEach((g, gi) => g.frags.forEach((fr, fi) => {
      for (const [x, a] of seen[gi][fi]) { const on = altOnSite.get(x); if (on !== undefined && a !== (gi === on ? 1 : 0)) { for (const r of fr) conflicting.add(r); break; } }
    }));
  }
  out.sort((a, b) => a.start - b.start || a.end - b.end);
  return { source, sets: out, assigned, unassigned, notSplit, checked, conflicting: countReads ? conflicting.size : new Set([...conflicting].map(r => r.n + r.s)).size };
}

/** Fisher p-value and difference of the two shares for a junction to count as used differently by the haplotypes. */
const ALLELIC_P = 1e-3, ALLELIC_DELTA = 0.2;

/**
 * Splice junctions of a set's haplotypes (RNA): the junctions carried by at least HAP_MIN_DEPTH fragments of the set,
 * and for each haplotype the fragments that carry one and those that go another way at one of its ends (another
 * junction from its donor or to its acceptor, an exon skip included, or aligned bases across the exon–intron
 * boundary). Fragments, not reads: two mates over one junction are one molecule. Junctions whose share differs
 * between the two haplotypes are tested (Fisher's exact test, two-sided) and listed.
 */
function haplotypeJunctions(groups: Group[], snap: Map<number, [number, number]> | null): { per: HapJunction[][]; allelic: AllelicJunction[] } {
  const none = { per: groups.map(() => []), allelic: [] };
  const units = groups.map(g => g.frags.map(fr => {
    const u: { fr: AlignedRead[]; js: Map<number, [number, number]> | null } = { fr, js: null };
    for (const r of fr) readJunctions(r, (_, a, b) => { const j = snap?.get(jKey(a, b)) ?? [a, b]; (u.js ??= new Map()).set(jKey(j[0], j[1]), j); });
    return u;
  }));
  const total = new Map<number, { s: number; e: number; n: number }>();
  for (const us of units) for (const u of us) if (u.js) for (const [k, [a, b]] of u.js) { const t = total.get(k); if (t) t.n++; else total.set(k, { s: a, e: b, n: 1 }); }
  const cand = [...total.entries()].filter(([, t]) => t.n >= HAP_MIN_DEPTH).sort((x, y) => x[1].s - y[1].s || x[1].e - y[1].e);
  if (!cand.length) return none;
  // bases p − 1 and p aligned in one block: the exon goes on into the intron (or the intron into the exon) at p
  const across = (fr: AlignedRead[], p: number) => fr.some(r => r.s < p && p < r.e && r.b.some(([a, b]) => a < p && p < b));
  const per = units.map(us => cand.map(([k, t]): HapJunction => {
    let n = 0, other = 0;
    for (const u of us) {
      if (u.js?.has(k)) { n++; continue; }
      let alt = false;
      if (u.js) for (const [a, b] of u.js.values()) if (a === t.s || b === t.e) { alt = true; break; }
      if (alt || across(u.fr, t.s) || across(u.fr, t.e)) other++;
    }
    return { start: t.s, end: t.e, n, other, psi: n + other ? n / (n + other) : 0 };
  }));
  const allelic: AllelicJunction[] = [];
  if (per.length >= 2) cand.forEach((_, i) => {
    const a = per[0][i], b = per[1][i];
    if (a.n + a.other < HAP_MIN_DEPTH || b.n + b.other < HAP_MIN_DEPTH) return;
    const p = fisher(a.n, a.other, b.n, b.other, 'two');
    if (p < ALLELIC_P && Math.abs(a.psi - b.psi) >= ALLELIC_DELTA) allelic.push({ start: a.start, end: a.end, n: [a.n, b.n], other: [a.other, b.other], psi: [a.psi, b.psi], p });
  });
  return { per: per.map(l => l.filter(j => j.n > 0)), allelic };
}

/** A fragment's allele at a site: its reads agreeing, unknown otherwise. */
function fragAllele(fr: AlignedRead[], s: VariantSite, minBq: number): number {
  let a = -1;
  for (const r of fr) { const x = alleleAt(r, s, minBq); if (x < 0) continue; if (a < 0) a = x; else if (a !== x) return -1; }
  return a;
}

const siteKey = (s: VariantSite) => `${s.pos}:${s.kind}:${s.alt}`;

/** Stretches of [start, end) where at least HAP_MIN_DEPTH reads align or delete bases. */
function coveredStretches(reads: AlignedRead[], start: number, end: number): [number, number][] {
  const w = end - start, diff = new Int32Array(w + 1);
  const add = (a: number, b: number) => { const x = Math.max(a, start) - start, y = Math.min(b, end) - start; if (y > x) { diff[x]++; diff[y]--; } };
  for (const r of reads) { for (const [a, b] of r.b) add(a, b); for (const [a, b] of r.d) add(a, b); }
  const out: [number, number][] = [];
  let d = 0, from = -1;
  for (let i = 0; i < w; i++) {
    d += diff[i];
    if (d >= HAP_MIN_DEPTH) { if (from < 0) from = i; }
    else if (from >= 0) { out.push([start + from, start + i]); from = -1; }
  }
  if (from >= 0) out.push([start + from, end]);
  return out;
}

/**
 * The two-haplotype answer of a reads window: the file's haplotags when the window has tagged reads (unless the
 * caller asks for the reads' own phasing), else the in-page read-based phasing, whose blocks are also returned.
 */
export function windowHaplotypes(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minVaf: number, minIndel: number, source: 'auto' | 'reads', phaseReads: () => PhaseResult, sites?: VariantSite[], longReads = false): { phase?: PhaseResult; haplotypes: HaplotypeView } {
  const o: Options = { start, end, ref, refStart, minBq: 20, minIndel, minAlt: 3, minVaf, sites, longReads };
  if (source !== 'reads' && reads.some(r => r.hp)) return { haplotypes: haplotypesFromTags(reads, o) };
  const phase = phaseReads();
  return { phase, haplotypes: haplotypesFromPhase(reads, phase, o) };
}
