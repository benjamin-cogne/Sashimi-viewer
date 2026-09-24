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
import type { AlignedRead, HaplotypeConsensus, HaplotypeSet, HaplotypeView, PhaseResult, VariantSite } from '../components/sashimi/types';
import { callSites } from './collapse';
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
interface Options { start: number; end: number; ref: string | null; refStart: number; minBq: number; minIndel: number; minAlt: number; minVaf: number; /** the window's sites over all reads, when the caller has them */ sites?: VariantSite[] }

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
    id: k === 'none' ? 'no phase set' : `PS ${Number(k).toLocaleString()}`, ps: gs[0].ps, groups: gs.sort((a, b) => a.hap - b.hap),
  }));
  return build('tags', sets, reads, assigned, unassigned, o, true);
}

/** Haplotypes from the in-page phasing: each fragment to the haplotype of the block it matches at more sites (ties unassigned). */
export function haplotypesFromPhase(reads: AlignedRead[], phase: PhaseResult, o: Options): HaplotypeView {
  const frags = fragmentsOf(reads);
  const sets = phase.blocks.map(b => ({ id: b.id, ps: null as number | null, groups: [1, 2].map(hap => ({ key: `${b.id}:${hap}`, ps: null, hap, frags: [] as AlignedRead[][] })) }));
  let assigned = 0, unassigned = 0;
  for (const fr of frags) {
    let best = -1, bestMargin = 0, bestHap = 0;
    phase.blocks.forEach((b, bi) => {
      let m1 = 0, m2 = 0;
      b.sites.forEach((si, k) => {
        const s = phase.sites[si];
        let a = -1;
        for (const r of fr) { const x = alleleAt(r, s, o.minBq); if (x < 0) continue; if (a < 0) a = x; else if (a !== x) { a = -1; break; } }
        if (a < 0) return;
        if ((a === 1) === (b.h1[k] === 'alt')) m1++; else m2++;
      });
      const margin = Math.abs(m1 - m2);
      if (margin > bestMargin) { best = bi; bestMargin = margin; bestHap = m1 > m2 ? 1 : 2; }
    });
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
  const het = all.filter(s => s.vaf >= HET_MIN && s.vaf <= HET_MAX);
  const windowKeys = new Set(all.map(siteKey));
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
    out.push({ id: set.id, ps: set.ps, start: s0, end: e0, haps });
    // the split check needs two haplotypes
    const [g1, g2] = set.groups;
    if (!g1 || !g2 || !g1.frags.length || !g2.frags.length) continue;
    for (const s of het) {
      if (s.pos < s0 || s.pos >= e0) continue;
      const frac = [g1, g2].map(g => {
        let alt = 0, n = 0;
        for (const fr of g.frags) { const a = fragAllele(fr, s, o.minBq); if (a < 0) continue; n++; if (a === 1) alt++; }
        return n >= HAP_MIN_DEPTH ? alt / n : NaN;
      });
      if (frac.some(Number.isNaN)) continue;
      checked++;
      const altOn = frac[0] >= SPLIT_MAJOR && frac[1] <= 1 - SPLIT_MAJOR ? 0 : frac[1] >= SPLIT_MAJOR && frac[0] <= 1 - SPLIT_MAJOR ? 1 : -1;
      if (altOn < 0) { notSplit.push({ pos: s.pos, kind: s.kind, alt: s.alt, fractions: frac, set: set.id }); continue; }
      [g1, g2].forEach((g, gi) => {
        const expect = gi === altOn ? 1 : 0;
        for (const fr of g.frags) { const a = fragAllele(fr, s, o.minBq); if (a >= 0 && a !== expect) for (const r of fr) conflicting.add(r); }
      });
    }
  }
  out.sort((a, b) => a.start - b.start || a.end - b.end);
  return { source, sets: out, assigned, unassigned, notSplit, checked, conflicting: countReads ? conflicting.size : new Set([...conflicting].map(r => r.n + r.s)).size };
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
  minVaf: number, minIndel: number, source: 'auto' | 'reads', phaseReads: () => PhaseResult, sites?: VariantSite[]): { phase?: PhaseResult; haplotypes: HaplotypeView } {
  const o: Options = { start, end, ref, refStart, minBq: 20, minIndel, minAlt: 3, minVaf, sites };
  if (source !== 'reads' && reads.some(r => r.hp)) return { haplotypes: haplotypesFromTags(reads, o) };
  const phase = phaseReads();
  return { phase, haplotypes: haplotypesFromPhase(reads, phase, o) };
}
