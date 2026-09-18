/**
 * Read-based phasing of the variant sites of a window: two haplotypes per block, linked through the reads and
 * their mates, in the spirit of WhatsHap / HapCUT2 in a greedy form that runs in milliseconds on a reads window.
 *
 *  1. Sites: every site called above the thresholds; those with an allele fraction in the heterozygous window
 *     (HET_MIN..HET_MAX) are the ones to phase. Homozygous sites (above HET_MAX) sit on both haplotypes; sites below
 *     HET_MIN (mosaic, subclonal, errors) are reported unphased.
 *  2. Fragments: a read and its mate are one fragment (mates identified by their positions, as the reads track
 *     does); the fragment's allele at every heterozygous site it covers, unknown when the base is of low quality,
 *     a third allele, or the two mates disagree.
 *  3. Links: for every pair of sites covered by the same fragments, the number of fragments seeing the same phase
 *     (ref–ref or alt–alt) and the opposite one. A link is trusted with at least MIN_LINK fragments and at most
 *     MAX_CONFLICT of them disagreeing.
 *  4. Blocks: the sites in order; each joins the current block when trusted links to sites already in it agree on
 *     its phase; else a new block starts (no link, or contradicting links: the reason is kept).
 *  5. Fragments are assigned to the haplotype they match best inside each block; support and disagreement counted.
 */
import type { AlignedRead, VariantSite, PhaseBlock, PhaseResult, UnphasedSite } from '../components/sashimi/types';
import { callSites } from './collapse';

export const HET_MIN = 0.25, HET_MAX = 0.75;
export const MIN_LINK = 2, MAX_CONFLICT = 0.2;

type Allele = 0 | 1 | -1;   // 0 ref, 1 alt, -1 unknown

/** Allele of one read at one site (unknown when not covered, low quality, or another allele). */
function alleleAt(r: AlignedRead, s: VariantSite, minBq: number): Allele {
  if (s.kind === 'snv') {
    if (!r.b.some(([bs, be]) => bs <= s.pos && s.pos < be)) return -1;
    for (const [pos, base, qual] of r.m) if (pos === s.pos) return qual < minBq ? -1 : base === s.alt ? 1 : -1;
    return 0;
  }
  if (s.kind === 'ins') {
    if (r.i.some(([p, l]) => p === s.pos && l === s.length)) return 1;
    return r.b.some(([bs, be]) => bs < s.pos && s.pos < be) ? 0 : -1;
  }
  if (r.d.some(([a, b]) => a === s.pos && b === s.pos + s.length)) return 1;
  return r.b.some(([bs, be]) => bs <= s.pos && s.pos < be) ? 0 : -1;
}

/** Reads joined into fragments: a read and its mate (both in the set) count once. */
export function fragmentsOf(reads: AlignedRead[]): AlignedRead[][] {
  const byStart = new Map<number, number[]>();
  reads.forEach((r, i) => { const l = byStart.get(r.s); if (l) l.push(i); else byStart.set(r.s, [i]); });
  const mate = new Int32Array(reads.length).fill(-1);
  reads.forEach((r, i) => {
    if (mate[i] >= 0 || r.mp == null || r.mc) return;
    for (const j of byStart.get(r.mp) ?? []) {
      if (j !== i && mate[j] < 0 && reads[j].mp === r.s && (reads[j].f & 192) !== (r.f & 192)) { mate[i] = j; mate[j] = i; break; }
    }
  });
  const out: AlignedRead[][] = [];
  reads.forEach((r, i) => { if (mate[i] >= 0 && mate[i] < i) return; out.push(mate[i] >= 0 ? [r, reads[mate[i]]] : [r]); });
  return out;
}

export function phaseReads(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20, minIndel = 1): PhaseResult {
  const sites = callSites(reads, start, end, ref, refStart, minAlt, minVaf, minBq, minIndel);
  const het = sites.map((s, i) => i).filter(i => sites[i].vaf >= HET_MIN && sites[i].vaf <= HET_MAX);
  const hom = sites.map((s, i) => i).filter(i => sites[i].vaf > HET_MAX);
  const unphased: UnphasedSite[] = sites.map((s, i) => i).filter(i => sites[i].vaf < HET_MIN).map(i => ({ site: i, reason: 'low' as const }));
  const frags = fragmentsOf(reads);
  // allele matrix: fragment × heterozygous site, evaluated only at the sites inside the fragment's span
  const hetPos = het.map(si => sites[si].pos);
  const lowerBound = (pos: number) => { let lo = 0, hi = hetPos.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (hetPos[mid] < pos) lo = mid + 1; else hi = mid; } return lo; };
  const alleles: Allele[][] = frags.map(fr => {
    const row: Allele[] = new Array(het.length).fill(-1);
    let lo = Infinity, hi = -Infinity;
    for (const r of fr) { if (r.s < lo) lo = r.s; if (r.e > hi) hi = r.e; }
    for (let k = lowerBound(lo); k < het.length && hetPos[k] < hi; k++) {
      let a: Allele = -1;
      for (const r of fr) { const x = alleleAt(r, sites[het[k]], minBq); if (x < 0) continue; if (a < 0) a = x; else if (a !== x) { a = -1; break; } }
      row[k] = a;
    }
    return row;
  });
  // links between pairs of heterozygous sites (indices into `het`)
  const same = new Map<string, number>(), diff = new Map<string, number>();
  const key = (a: number, b: number) => `${a},${b}`;
  const bump = (m: Map<string, number>, k: string) => m.set(k, (m.get(k) ?? 0) + 1);
  for (const row of alleles) {
    const covered: number[] = []; row.forEach((a, k) => { if (a >= 0) covered.push(k); });
    for (let x = 0; x < covered.length; x++) for (let y = x + 1; y < covered.length; y++) {
      const i = covered[x], j = covered[y];
      bump(row[i] === row[j] ? same : diff, key(i, j));
    }
  }
  // greedy blocks; phase[k] = 0 when the alt allele of het site k is on H1, 1 when on H2
  const phase = new Int8Array(het.length).fill(0);
  const blocks: PhaseBlock[] = [];
  let cur: number[] = [];
  let nextReason: PhaseBlock['breakBefore'] = undefined;
  const linksOf = (j: number) => cur.map(i => ({ i, s: same.get(key(i, j)) ?? 0, d: diff.get(key(i, j)) ?? 0 })).filter(l => l.s + l.d >= MIN_LINK && Math.min(l.s, l.d) / (l.s + l.d) <= MAX_CONFLICT);
  const close = () => {
    if (!cur.length) return;
    const first = het[cur[0]], last = het[cur[cur.length - 1]];
    blocks.push({ id: `block ${blocks.length + 1}`, start: sites[first].pos, end: sites[last].pos + Math.max(1, sites[last].length), sites: cur.map(k => het[k]),
      h1: cur.map(k => (phase[k] === 0 ? 'alt' : 'ref')), h2: cur.map(k => (phase[k] === 0 ? 'ref' : 'alt')),
      support: [0, 0], ambiguous: 0, conflicting: 0, links: [], breakBefore: nextReason });
    cur = [];
  };
  for (let j = 0; j < het.length; j++) {
    if (!cur.length) { cur.push(j); continue; }
    const links = linksOf(j);
    if (!links.length) { close(); nextReason = 'no link'; cur.push(j); continue; }
    // vote: each trusted link says whether j is in phase with i (same) or opposite, weighted by its fragments
    let vote0 = 0, vote1 = 0;
    for (const l of links) { const inPhase = l.s >= l.d; const p = inPhase ? phase[l.i] : 1 - phase[l.i]; if (p === 0) vote0 += l.s + l.d; else vote1 += l.s + l.d; }
    if (Math.min(vote0, vote1) / (vote0 + vote1) > MAX_CONFLICT) { close(); nextReason = 'conflict'; cur.push(j); continue; }
    phase[j] = vote0 >= vote1 ? 0 : 1;
    cur.push(j);
  }
  close();
  // fragments on haplotypes, per block; links kept for the tooltips
  const hetIndex = new Map(het.map((si, k) => [si, k]));
  for (const b of blocks) {
    const ks = b.sites.map(si => hetIndex.get(si)!);
    for (const row of alleles) {
      let m1 = 0, m2 = 0, n = 0;
      for (const k of ks) { const a = row[k]; if (a < 0) continue; n++; const alt1 = phase[k] === 0; if ((a === 1) === alt1) m1++; else m2++; }
      if (!n) continue;
      if (m1 === m2) b.ambiguous++;
      else { b.support[m1 > m2 ? 0 : 1]++; if (Math.min(m1, m2) > 0) b.conflicting++; }
    }
    for (let x = 0; x < ks.length; x++) for (let y = x + 1; y < ks.length; y++) {
      const s = same.get(key(ks[x], ks[y])) ?? 0, d = diff.get(key(ks[x], ks[y])) ?? 0;
      if (s + d) b.links.push({ a: b.sites[x], b: b.sites[y], same: s, diff: d });
    }
  }
  // single-site blocks are heterozygous sites no fragment ties to another: unphased, not haplotypes
  const phased = blocks.filter(b => b.sites.length >= 2).map((b, i) => ({ ...b, id: `block ${i + 1}` }));
  for (const b of blocks) if (b.sites.length < 2) unphased.push({ site: b.sites[0], reason: b.breakBefore === 'conflict' ? 'conflict' : 'unlinked' });
  unphased.sort((a, b) => sites[a.site].pos - sites[b.site].pos);
  return { sites, hom, blocks: phased, unphased, fragments: frags.length, het: het.length };
}
