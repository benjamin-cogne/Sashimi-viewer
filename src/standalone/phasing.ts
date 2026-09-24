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
import { callSites, fragmentsOf, inBlocks, mismatchAt } from './collapse';

export const HET_MIN = 0.25, HET_MAX = 0.75;
export const MIN_LINK = 2, MAX_CONFLICT = 0.2;
/** Sites after each one (in a fragment) that it is linked to. */
const LINK_AHEAD = 32;

export type Allele = 0 | 1 | -1;   // 0 ref, 1 alt, -1 unknown

/** Allele of one read at one site (unknown when not covered, low quality, or another allele); lookups by binary search. */
export function alleleAt(r: AlignedRead, s: VariantSite, minBq: number): Allele {
  if (s.kind === 'snv') {
    if (!inBlocks(r, s.pos)) return -1;
    const m = mismatchAt(r, s.pos);
    if (!m) return 0;
    return m[1] < minBq ? -1 : m[0] === s.alt ? 1 : -1;
  }
  if (s.kind === 'ins') {
    if (r.i.some(([p, l]) => p === s.pos && l === s.length)) return 1;
    return inBlocks(r, s.pos, true) ? 0 : -1;
  }
  if (r.d.some(([a, b]) => a === s.pos && b === s.pos + s.length)) return 1;
  return inBlocks(r, s.pos) ? 0 : -1;
}

/** Reads joined into fragments: a read and its mate (both in the set) count once (shared with the consensus groups). */
export { fragmentsOf };

export function phaseReads(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20, minIndel = 1, called?: VariantSite[]): PhaseResult {
  // the window's sites, called by the caller already or here
  const sites = called ?? callSites(reads, start, end, ref, refStart, minAlt, minVaf, minBq, minIndel);
  const het = sites.map((s, i) => i).filter(i => sites[i].vaf >= HET_MIN && sites[i].vaf <= HET_MAX);
  const hom = sites.map((s, i) => i).filter(i => sites[i].vaf > HET_MAX);
  const unphased: UnphasedSite[] = sites.map((s, i) => i).filter(i => sites[i].vaf < HET_MIN).map(i => ({ site: i, reason: 'low' as const }));
  const frags = fragmentsOf(reads);
  // each fragment's known alleles at the heterozygous sites inside its span, sparse (indices into `het`, in order):
  // a dense fragment × site matrix took gigabytes on a window with thousands of sites
  const hetPos = het.map(si => sites[si].pos);
  const lowerBound = (pos: number) => { let lo = 0, hi = hetPos.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (hetPos[mid] < pos) lo = mid + 1; else hi = mid; } return lo; };
  const rows: { k: number[]; a: Allele[] }[] = frags.map(fr => {
    const row = { k: [] as number[], a: [] as Allele[] };
    let lo = Infinity, hi = -Infinity;
    for (const r of fr) { if (r.s < lo) lo = r.s; if (r.e > hi) hi = r.e; }
    for (let k = lowerBound(lo); k < het.length && hetPos[k] < hi; k++) {
      let a: Allele = -1;
      for (const r of fr) { const x = alleleAt(r, sites[het[k]], minBq); if (x < 0) continue; if (a < 0) a = x; else if (a !== x) { a = -1; break; } }
      if (a >= 0) { row.k.push(k); row.a.push(a); }
    }
    return row;
  });
  // links between pairs of heterozygous sites (indices into `het`)
  // (numeric keys; each site linked to the LINK_AHEAD next ones a fragment covers: a long read over hundreds of sites
  // made tens of thousands of pairs, and the blocks are built from neighbouring sites anyway)
  const same = new Map<number, number>(), diff = new Map<number, number>();
  const H = het.length;
  const key = (a: number, b: number) => a * H + b;
  const bump = (m: Map<number, number>, k: number) => m.set(k, (m.get(k) ?? 0) + 1);
  for (const row of rows) {
    for (let x = 0; x < row.k.length; x++) for (let y = x + 1; y < row.k.length && y <= x + LINK_AHEAD; y++) bump(row.a[x] === row.a[y] ? same : diff, key(row.k[x], row.k[y]));
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
  /** site j against the current block: its phase set when its trusted links agree, else why it does not fit */
  const tryJoin = (j: number): 'ok' | 'no link' | 'conflict' => {
    const links = linksOf(j);
    if (!links.length) return 'no link';
    // vote: each trusted link says whether j is in phase with i (same) or opposite, weighted by its fragments
    let vote0 = 0, vote1 = 0;
    for (const l of links) { const inPhase = l.s >= l.d; const p = inPhase ? phase[l.i] : 1 - phase[l.i]; if (p === 0) vote0 += l.s + l.d; else vote1 += l.s + l.d; }
    if (Math.min(vote0, vote1) / (vote0 + vote1) > MAX_CONFLICT) return 'conflict';
    phase[j] = vote0 >= vote1 ? 0 : 1;
    return 'ok';
  };
  // a site that does not fit while the next one does is an outlier (an error-made or mistyped site, a homozygous one
  // read as heterozygous): it is left unphased and the block goes on, instead of ending there
  const skipped: UnphasedSite[] = [];
  for (let j = 0; j < het.length; j++) {
    if (!cur.length) { cur.push(j); continue; }
    const r = tryJoin(j);
    if (r === 'ok') { cur.push(j); continue; }
    if (j + 1 < het.length && tryJoin(j + 1) === 'ok') { skipped.push({ site: het[j], reason: r === 'conflict' ? 'conflict' : 'unlinked' }); cur.push(j + 1); j++; continue; }
    close(); nextReason = r; cur.push(j);
  }
  close();
  unphased.push(...skipped);
  // fragments on haplotypes, per block (each fragment through its own sites only); links kept for the tooltips
  const hetIndex = new Map(het.map((si, k) => [si, k]));
  const blockOf = new Int32Array(het.length).fill(-1);
  blocks.forEach((b, bi) => { for (const si of b.sites) blockOf[hetIndex.get(si)!] = bi; });
  for (const row of rows) {
    const tally = new Map<number, [number, number]>();
    row.k.forEach((k, x) => {
      const bi = blockOf[k];
      if (bi < 0) return;
      let t = tally.get(bi);
      if (!t) tally.set(bi, t = [0, 0]);
      if ((row.a[x] === 1) === (phase[k] === 0)) t[0]++; else t[1]++;
    });
    for (const [bi, [m1, m2]] of tally) {
      const b = blocks[bi];
      if (m1 === m2) b.ambiguous++;
      else { b.support[m1 > m2 ? 0 : 1]++; if (Math.min(m1, m2) > 0) b.conflicting++; }
    }
  }
  for (const b of blocks) {
    const ks = b.sites.map(si => hetIndex.get(si)!);
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
