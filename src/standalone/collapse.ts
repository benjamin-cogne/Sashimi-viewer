/**
 * Variant sites of a reads window, and its consensus groups (Collapse, Haplotypes: any): one group per local
 * haplotype × splice pattern that the reads support. All coordinates are 0-based half-open.
 *
 * - Units are fragments: a read and its mate are one molecule, so their sites and junctions are seen together.
 * - Each fragment has a pattern: its splice junctions (the gaps between aligned blocks that are not deletions) and its
 *   allele at every called site it covers.
 * - Unspliced short reads (DNA) are grouped by exact patterns: seeds of at least *Min reads* fragments, most
 *   specific first, and every pattern placed on the one seed it fits (ambiguous when several, minor when none).
 * - Long reads, and spliced short reads (RNA-seq), are clustered with a tolerance (clusterTolerant): pieces of one
 *   isoform of one haplotype seen by different reads join through the sites and junctions they share, and through
 *   the fragments linking them; splicing never disagrees within a group.
 */
import type { AlignedRead, ReadGroup, VariantSite } from '../components/sashimi/types';
import { depthArray } from './alignments';

const BASE_CODE: Record<string, number> = { A: 0, C: 1, G: 2, T: 3, N: 4 };
const CODE_BASE = 'ACGTN';
/** position × INDEL_KEY + length, for insertions and deletions (positions up to 1e9 stay exact) */
const INDEL_KEY = 8_388_608;

/**
 * Variant sites of a set of reads over [start, end): an allele needs `minAlt` reads and `minVaf` of the depth. Every
 * alternate base counts, whatever its base quality, as in the full scan (alleles.ts): leaving the low-quality ones out
 * of the count while their reads stay in the depth under-read every allele fraction, by a third on ONT data, where it
 * took homozygous sites (~0.67) for heterozygous ones and dropped heterozygous ones (~0.33) under the phasing window.
 * `minBq` is kept for the callers' signature; the phasing (alleleAt) still trusts only alleles of that quality.
 * Keys are numbers, not strings: this runs several times per collapse over every mismatch of every read.
 */
export function callSites(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, _minBq = 20, minIndel = 1): VariantSite[] {
  const depth = depthArray(reads, start, end), w = Math.max(0, end - start);
  // SNV counts in one typed array (position × 5 bases): no Map, no key, nothing for the collector
  const snv = new Int32Array(w * 5), otherSnv = new Map<string, number>(), ins = new Map<number, number>(), del = new Map<number, number>();
  const bump = <K,>(m: Map<K, number>, k: K) => m.set(k, (m.get(k) ?? 0) + 1);
  for (const r of reads) {
    for (const m of r.m) {
      const pos = m[0];
      if (pos < start || pos >= end) continue;
      const c = BASE_CODE[m[1]];
      if (c === undefined) bump(otherSnv, `${pos}\t${m[1]}`); else snv[(pos - start) * 5 + c]++;
    }
    for (const [pos, len] of r.i) if (len >= minIndel && pos >= start && pos < end) bump(ins, pos * INDEL_KEY + Math.min(len, INDEL_KEY - 1));
    for (const [ds, de] of r.d) if (de - ds >= minIndel && ds >= start && ds < end) bump(del, ds * INDEL_KEY + Math.min(de - ds, INDEL_KEY - 1));
  }
  const sites: VariantSite[] = [];
  const addSnv = (pos: number, base: string, n: number) => {
    const d = depth[pos - start];
    if (n >= minAlt && d && n / d >= minVaf) {
      const rb = ref && pos - refStart >= 0 && pos - refStart < ref.length ? ref[pos - refStart] : '?';
      sites.push({ pos, kind: 'snv', ref: rb, alt: base, length: 0, alt_count: n, depth: d, vaf: n / d });
    }
  };
  for (let i = 0; i < snv.length; i++) if (snv[i] >= minAlt) addSnv(start + Math.floor(i / 5), CODE_BASE[i % 5], snv[i]);
  for (const [k, n] of otherSnv) { const [p, base] = k.split('\t'); addSnv(parseInt(p), base, n); }
  for (const [k, n] of ins) {
    const pos = Math.floor(k / INDEL_KEY), len = k % INDEL_KEY;
    const d = pos - start < depth.length ? depth[pos - start] : 0;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos, kind: 'ins', ref: '', alt: `+${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  for (const [k, n] of del) {
    const ds = Math.floor(k / INDEL_KEY), len = k % INDEL_KEY;
    const d = depth[ds - start] + n;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos: ds, kind: 'del', ref: '', alt: `-${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  return sites.sort((x, y) => x.pos - y.pos || x.kind.localeCompare(y.kind) || (x.alt < y.alt ? -1 : x.alt > y.alt ? 1 : 0));
}

/**
 * A read's mismatches by position, for lookups in O(log n): alleleAt and readAlleles used to scan them all for every
 * site, which on a long read (hundreds of sequencing errors) times every site of the window was most of a collapse.
 * Built once per read, kept with it.
 */
const mismatchIndex = new WeakMap<AlignedRead, { pos: Int32Array; at: Int32Array }>();
/** The mismatch of `r` at `pos` (base, quality), or null. */
export function mismatchAt(r: AlignedRead, pos: number): [string, number] | null {
  let ix = mismatchIndex.get(r);
  if (!ix) {
    const order = r.m.map((_, i) => i).sort((i, j) => r.m[i][0] - r.m[j][0]);
    ix = { pos: Int32Array.from(order, i => r.m[i][0]), at: Int32Array.from(order) };
    mismatchIndex.set(r, ix);
  }
  let lo = 0, hi = ix.pos.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (ix.pos[mid] < pos) lo = mid + 1; else hi = mid; }
  if (lo >= ix.pos.length || ix.pos[lo] !== pos) return null;
  const m = r.m[ix.at[lo]];
  return [m[1], m[2]];
}
/** Whether `pos` lies in one of the read's aligned blocks (sorted, disjoint): bs <= pos < be, or bs < pos < be when `inner`. */
export function inBlocks(r: AlignedRead, pos: number, inner = false): boolean {
  const b = r.b;
  let lo = 0, hi = b.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (b[mid][1] <= pos) lo = mid + 1; else hi = mid; }
  return lo < b.length && (inner ? b[lo][0] < pos : b[lo][0] <= pos) && pos < b[lo][1];
}

/** Index of the first site at or after `pos` (sites sorted by position). */
function lowerBound(pos: Int32Array, x: number): number {
  let lo = 0, hi = pos.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (pos[mid] < x) lo = mid + 1; else hi = mid; }
  return lo;
}

/**
 * A read's alleles at the sites it says something of, into `out` (site index → 'ref', 'alt' or another base). Only
 * the sites under the read are looked at (binary search on their positions): a short read covers one or two of a
 * window's hundreds of sites, and all of them were tested for every read. Sparse: a window with thousands of sites (a
 * wrong reference, a hypermutated stretch) held one array of every site per read, gigabytes on 40,000 reads.
 */
function readAlleles(r: AlignedRead, sites: VariantSite[], sitePos: Int32Array, minBq: number, out: Map<number, string>): void {
  for (let k = lowerBound(sitePos, r.s), hi = lowerBound(sitePos, r.e + 1); k < hi; k++) {
    const s = sites[k];
    let a: string;
    if (s.kind === 'snv') {
      if (!inBlocks(r, s.pos)) a = '.';
      else { const m = mismatchAt(r, s.pos); a = m ? (m[1] < minBq ? '.' : (m[0] === s.alt ? 'alt' : m[0])) : 'ref'; }
    } else if (s.kind === 'ins') {
      const alt = r.i.some(([p, l]) => p === s.pos && l === s.length);
      a = alt ? 'alt' : (inBlocks(r, s.pos, true) || r.i.some(([p]) => p === s.pos)) ? 'ref' : '.';
    } else a = r.d.some(([x, y]) => x === s.pos && y === s.pos + s.length) ? 'alt' : inBlocks(r, s.pos) ? 'ref' : '.';
    // the two mates of a fragment over the same site: kept when they agree, unknown when they do not
    if (a === '.') continue;
    const was = out.get(k);
    out.set(k, was === undefined || was === a ? a : '?');
  }
}

/**
 * Pairs a read with its mate (same start as the mate's recorded position and the reverse, first/second of pair
 * flags differing): a fragment is one molecule, whose two mates together link sites and junctions further apart
 * than either read does.
 */
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

/** Numeric key of a junction [start, end): start × INDEL_KEY + length (introns of a reads window are far under 8 Mb). */
export const jKey = (s: number, e: number) => s * INDEL_KEY + Math.min(e - s, INDEL_KEY - 1);

/** The gaps between a read's aligned blocks that are introns (N), not deletions nor insertion splits. */
export function readJunctions(r: AlignedRead, visit: (k: number, be: number, ns: number) => void): void {
  // blocks and deletions are both in reference order: one pointer walks the deletions (a long read has hundreds)
  let d = 0;
  for (let k = 0; k + 1 < r.b.length; k++) {
    const be = r.b[k][1], ns = r.b[k + 1][0];
    if (ns <= be) continue;
    while (d < r.d.length && r.d[d][0] < be) d++;
    if (d < r.d.length && r.d[d][0] === be && r.d[d][1] === ns) continue;
    visit(k, be, ns);
  }
}

/**
 * Long reads place a junction a few bases off where the bases next to it carry errors (ONT especially), so the reads
 * of one isoform disagree on it and no two would share a splice pattern. A junction seen in few reads, within
 * JUNCTION_SNAP_BP at both ends of one seen at least 1 / JUNCTION_SNAP_RATIO times as often, is taken for that one:
 * IsoQuant corrects to the annotation within 6 bp on ONT data (its `delta`, 4 on PacBio), FLAIR within 15. Two real
 * splice sites as close as that (NAGNAG acceptors, 3 bp) stay apart when both are common.
 */
const JUNCTION_SNAP_BP = 6, JUNCTION_SNAP_RATIO = 4;
export function junctionSnap(reads: AlignedRead[]): Map<number, [number, number]> {
  const count = new Map<number, { s: number; e: number; n: number }>();
  for (const r of reads) readJunctions(r, (_, s, e) => { const k = jKey(s, e); const c = count.get(k); if (c) c.n++; else count.set(k, { s, e, n: 1 }); });
  const all = [...count.values()].sort((a, b) => b.n - a.n || a.s - b.s);
  const byStart = [...all].sort((a, b) => a.s - b.s);
  const snap = new Map<number, [number, number]>();
  const done = new Set<number>();
  for (const c of all) {
    const ck = jKey(c.s, c.e);
    if (done.has(ck)) continue;
    done.add(ck);
    let lo = 0, hi = byStart.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (byStart[mid].s < c.s - JUNCTION_SNAP_BP) lo = mid + 1; else hi = mid; }
    for (let i = lo; i < byStart.length && byStart[i].s <= c.s + JUNCTION_SNAP_BP; i++) {
      const o = byStart[i], ok = jKey(o.s, o.e);
      if (done.has(ok) || Math.abs(o.e - c.e) > JUNCTION_SNAP_BP || o.n * JUNCTION_SNAP_RATIO > c.n) continue;
      done.add(ok); snap.set(ok, [c.s, c.e]);
    }
  }
  return snap;
}

function mergeBlocks(blocks: [number, number][]): [number, number][] {
  const out: [number, number][] = [];
  for (const [s, e] of [...blocks].sort((a, b) => a[0] - b[0] || a[1] - b[1])) {
    if (out.length && s <= out[out.length - 1][1]) out[out.length - 1][1] = Math.max(out[out.length - 1][1], e);
    else out.push([s, e]);
  }
  return out;
}

/**
 * A set of fragments sharing one pattern (splice junctions and alleles), or a group being built from several.
 * `blocks` are merged (sorted, disjoint), with the junctions' ends snapped as the chain's are; `n` counts reads,
 * `f` fragments (the support: two mates are one molecule).
 */
interface Group {
  chain: [number, number][]; jkeys: Set<number>; n: number; f: number; reads: AlignedRead[]; blocks: [number, number][];
  /** alleles at the sites it covers (site index → 'ref', 'alt' or another base), none where it says nothing */
  calls: Map<number, string>;
  start: number; end: number; spec: number;
}
const groupSpec = (g: Group) => g.chain.length + g.calls.size;
function finish(g: Group): Group {
  g.blocks = mergeBlocks(g.blocks);
  g.start = g.blocks.length ? g.blocks[0][0] : 0; g.end = g.blocks.length ? g.blocks[g.blocks.length - 1][1] : 0;
  g.spec = groupSpec(g);
  return g;
}

/** Whether any of the sorted, disjoint `blocks` overlaps [s, e). */
function blocksOverlap(blocks: [number, number][], s: number, e: number): boolean {
  let lo = 0, hi = blocks.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (blocks[mid][1] <= s) lo = mid + 1; else hi = mid; }
  return lo < blocks.length && blocks[lo][0] < e;
}

/**
 * Splicing of two groups where both are seen: no junction of one crossed by a block of the other (an exon the
 * other skips, a retained intron) and no two different junctions overlapping (another donor or acceptor).
 */
function structureAgrees(a: Group, b: Group): boolean {
  if (!a.chain.length && !b.chain.length) return true;
  for (const [x, y] of [[a, b], [b, a]] as const) {
    for (const [js, je] of x.chain) {
      if (y.end <= js || je <= y.start) continue;
      if (blocksOverlap(y.blocks, js, je)) return false;
      if (!y.jkeys.has(jKey(js, je)) && y.chain.some(([ks, ke]) => ks < je && js < ke)) return false;
    }
  }
  return true;
}

/**
 * Short reads: a group joins an anchor that holds every junction it has, has none of its blocks across the anchor's
 * junctions, and calls every site it calls with the same allele. A group without any call or junction joins what it
 * overlaps.
 */
function compatible(g: Group, anchor: Group): boolean {
  if (g.end <= anchor.start || anchor.end <= g.start) return false;
  for (const k of g.jkeys) if (!anchor.jkeys.has(k)) return false;
  for (const [js, je] of anchor.chain) if (js < g.end && g.start < je && blocksOverlap(g.blocks, js, je)) return false;
  for (const [i, a] of g.calls) if (anchor.calls.get(i) !== a) return false;
  return g.chain.length + g.calls.size > 0 || g.spec === 0;
}

/**
 * Long reads, clustered with a tolerance: a read (or a set of identical ones) joins the group whose consensus it
 * agrees with at their shared sites, up to TOLERANT_MAX_DISAGREE of them disagreeing (a sequencing error, a
 * mistyped site), and any read may seed a group. Exact patterns cannot group long reads: each covers its own run of
 * sites and carries its own errors, so no two are identical and nearly all fell into "minor" (301 of 324 on a
 * synthetic 30x ONT window, two haplotypes). A read agreeing with two groups goes to the one it agrees with by
 * TOLERANT_MARGIN more sites, else stays ambiguous between them. Splicing (RNA) is a hard constraint, checked
 * against the whole group (its junctions and blocks so far, not its first read's): a read of the same isoform that
 * starts or ends elsewhere joins, one that skips an exon the group holds does not; shared junctions count as
 * agreement, never towards the disagreeing share.
 */
const TOLERANT_MAX_DISAGREE = 0.2, TOLERANT_MARGIN = 2;

interface Cluster { g: Group; members: Group[]; tally: Map<number, Map<string, number>>; n: number; calls: number }
interface Score { i: number; score: number; sites: number }

/**
 * Agreement of a group with a cluster, or null when their splicing, or more than TOLERANT_MAX_DISAGREE of their
 * shared sites, disagree. `sites`: shared sites with the same allele; `score`: those plus the shared junctions, less
 * the disagreeing sites. A group with calls meets a cluster with calls only through shared sites: junctions both
 * haplotypes splice say nothing of the phase, and would join the alleles of a site to either one.
 */
function agreement(g: Group, calls: [number, string][], c: Cluster): Omit<Score, 'i'> | null {
  if (g.end <= c.g.start || c.g.end <= g.start || !structureAgrees(g, c.g)) return null;
  let agree = 0, disagree = 0;
  for (const [i, a] of calls) { const b = c.g.calls.get(i); if (b === undefined) continue; if (b === a) agree++; else disagree++; }
  if (disagree > TOLERANT_MAX_DISAGREE * (agree + disagree)) return null;
  if (!calls.length && !g.chain.length) return { score: 0, sites: 0 };
  if (calls.length && c.calls && !(agree + disagree)) return null;
  let shared = agree;
  for (const k of g.jkeys) if (c.g.jkeys.has(k)) shared++;
  return shared >= 1 ? { score: shared - disagree, sites: agree } : null;
}

/** Fragments that agree with two groups at a site of each, needed to join groups that share no site (short reads). */
const MIN_BRIDGES = 2;
/** Share of the patterns seeding a cluster of their own past which the reads are taken to agree on nothing. */
const CLUSTER_GIVE_UP = 0.3;
/**
 * Sites where nearly every read carries the alternate allele (homozygous) are no evidence of which group a read
 * belongs to: every read agrees with every group there, which joined the haplotypes through them. They stay in the
 * consensus shown, not in the clustering. An allele-specific expression of 80–90 % still counts (RNA).
 */
const HOM_VAF = 0.9;

/**
 * Fills anchors / members / ambiguous / minor from the patterns, clustering them with a tolerance (see above): long
 * reads, and spliced short reads (RNA-seq), whose pieces along a transcript (exons 2–3, 3–4, …) differ in the sites
 * and junctions they see, so that exact patterns split one isoform of one haplotype into dozens of rows.
 */
function clusterTolerant(ordered: Group[], sites: VariantSite[], minSupport: number, anchors: Group[], members: Group[][],
  ambiguous: Map<string, { hits: number[]; groups: Group[] }>, minor: Group[]): boolean {
  const informative = sites.map(s => s.vaf <= HOM_VAF);
  const clusters: Cluster[] = [];
  const join = (c: Cluster, g: Group) => {
    for (const [i, a] of g.calls) {
      let t = c.tally.get(i);
      if (!t) c.tally.set(i, t = new Map());
      t.set(a, (t.get(a) ?? 0) + g.n);
      let best = '', bn = 0, tie = false;
      for (const [k, v] of t) { if (v > bn) { best = k; bn = v; tie = false; } else if (v === bn) tie = true; }
      const had = c.g.calls.has(i);
      if (tie) c.g.calls.delete(i); else c.g.calls.set(i, best);
      if (informative[i]) c.calls += (tie ? 0 : 1) - (had ? 1 : 0);
    }
    let added = false;
    for (const j of g.chain) { const k = jKey(j[0], j[1]); if (!c.g.jkeys.has(k)) { c.g.jkeys.add(k); c.g.chain.push(j); added = true; } }
    if (added) c.g.chain.sort((x, y) => x[0] - y[0]);
    c.g.blocks = mergeBlocks([...c.g.blocks, ...g.blocks]);
    c.g.start = Math.min(c.g.start, g.start); c.g.end = Math.max(c.g.end, g.end);
    c.n += g.n;
  };
  const callsOf = (g: Group) => { const out: [number, string][] = []; for (const e of g.calls) if (informative[e[0]]) out.push(e); return out; };
  // a read with calls goes by them: groups sharing none of its sites (junctions only) compete only when none does
  const score = (g: Group, pool: number[]) => {
    const calls = callsOf(g);
    let scored: Score[] = [];
    for (const ci of pool) { const s = agreement(g, calls, clusters[ci]); if (s) scored.push({ i: ci, ...s }); }
    if (scored.some(s => s.sites > 0)) scored = scored.filter(s => s.sites > 0);
    return scored.sort((a, b) => b.score - a.score);
  };
  const decided = (s: Score[]) => s.length === 1 || (s.length > 1 && s[0].score - s[1].score >= TOLERANT_MARGIN);
  // left to right, so that each read meets the group of the reads just before it over the sites they share (by how
  // informative they are, groups were seeded all along the window and never met)
  const byPos = [...ordered].sort((a, b) => a.start - b.start || b.spec - a.spec);
  /**
   * The clusters overlapping [s, e) among `pool`: sorted by start, with the running maximum of their ends, the scan
   * goes back from the last one starting before e and stops where no earlier one reaches s. Scoring every group
   * against every cluster was the whole cost of a window with thousands of sites (a wrong reference: 158 s).
   */
  const indexOf = (pool: number[]) => {
    const sorted = [...pool].sort((a, b) => clusters[a].g.start - clusters[b].g.start);
    const maxEnd = new Float64Array(sorted.length);
    sorted.forEach((c, k) => { maxEnd[k] = Math.max(k ? maxEnd[k - 1] : -Infinity, clusters[c].g.end); });
    return (st: number, en: number) => {
      let lo = 0, hi = sorted.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (clusters[sorted[mid]].g.start < en) lo = mid + 1; else hi = mid; }
      const out: number[] = [];
      for (let k = lo - 1; k >= 0 && maxEnd[k] > st; k--) if (clusters[sorted[k]].g.end > st) out.push(sorted[k]);
      return out;
    };
  };
  // first pass: a cluster ending before the group starts can take no later group (they start further on)
  let open: number[] = [];
  let seen = 0;
  for (const g of byPos) {
    // reads that agree with nobody (a wrong reference, a window of junk: every read its own "site" pattern) have no
    // groups to find, and scoring each against all the others took minutes: the exact patterns do instead
    if (++seen >= 300 && clusters.length > CLUSTER_GIVE_UP * seen) return false;
    open = open.filter(ci => clusters[ci].g.end > g.start);
    const scored = score(g, open);
    if (decided(scored)) { join(clusters[scored[0].i], g); clusters[scored[0].i].members.push(g); continue; }
    if (!scored.length) {
      if (g.spec > 0) {
        const c: Cluster = { g: { ...g, chain: [], jkeys: new Set(), calls: new Map(), blocks: [], reads: [] }, members: [g], tally: new Map(), n: 0, calls: 0 };
        clusters.push(c); join(c, g); open.push(clusters.length - 1);
      }
      continue;
    }
  }
  // links: fragments agreeing with two groups at a site of each (counted once the groups exist: a fragment seen
  // before the second group was seeded links it too)
  const bridges = clusters.map(() => new Map<number, number>());
  const link = (x: number, y: number, n: number) => { bridges[x].set(y, (bridges[x].get(y) ?? 0) + n); bridges[y].set(x, (bridges[y].get(x) ?? 0) + n); };
  const near = indexOf(clusters.map((_, i) => i));
  for (const g of byPos) {
    if (!callsOf(g).length) continue;
    const linked = score(g, near(g.start, g.end)).filter(s => s.sites > 0);
    for (let x = 0; x < linked.length; x++) for (let y = x + 1; y < linked.length; y++) link(linked[x].i, linked[y].i, g.f);
  }
  // two groups that overlap and differ (splicing, or an allele): a third one agreeing with both is not either's
  const conflict = (a: Group, b: Group) => {
    if (a.end <= b.start || b.end <= a.start) return false;
    if (!structureAgrees(a, b)) return true;
    const [x, y] = a.calls.size <= b.calls.size ? [a, b] : [b, a];
    for (const [i, v] of x.calls) { const w = y.calls.get(i); if (informative[i] && w !== undefined && w !== v) return true; }
    return false;
  };
  // groups that are one are merged: at least TOLERANT_MARGIN shared sites agreeing, or MIN_BRIDGES fragments linking
  // them (short reads, where two groups of one haplotype often share no site), or, for two groups without any call,
  // TOLERANT_MARGIN shared junctions with no other group as close; at most TOLERANT_MAX_DISAGREE of their shared
  // sites disagreeing, and their splicing agreeing
  // A worklist, not rounds over every pair: a group is looked at again only when it grew (it may now meet others),
  // or when it was held back by two differing candidates and a merge nearby changed them.
  const alive = clusters.map(() => true);
  const queue = clusters.map((_, i) => i), queued = new Uint8Array(clusters.length).fill(1), held = new Set<number>();
  const push = (i: number) => { if (!queued[i] && alive[i]) { queued[i] = 1; queue.push(i); } };
  for (let head = 0; head < queue.length; head++) {
    const x = queue[head];
    queued[x] = 0;
    if (!alive[x]) continue;
    {
      const cx = clusters[x].g, calls = callsOf(cx);
      const ok: Score[] = [];
      for (let y = 0; y < clusters.length; y++) {
        if (y === x || !alive[y]) continue;
        const cy = clusters[y].g;
        if (cy.end <= cx.start || cx.end <= cy.start) continue;
        const s = agreement(cx, calls, clusters[y]);
        if (!s) continue;
        if (s.sites >= TOLERANT_MARGIN || (bridges[x].get(y) ?? 0) >= MIN_BRIDGES || (!clusters[x].calls && !clusters[y].calls && s.score >= TOLERANT_MARGIN)) ok.push({ i: y, ...s });
      }
      if (!ok.length) continue;
      ok.sort((a, b) => b.score - a.score);
      // agreeing with two groups that differ from each other (two isoforms, two haplotypes where it sees no site of
      // either), the group stays apart: joining it to either would claim what the reads do not show
      if (ok.some((a, k) => ok.some((b, l) => l > k && conflict(clusters[a.i].g, clusters[b.i].g)))) { held.add(x); continue; }
      const best = ok[0].i;
      const [keep, drop] = clusters[x].n >= clusters[best].n ? [x, best] : [best, x];
      for (const g of clusters[drop].members) { join(clusters[keep], g); clusters[keep].members.push(g); }
      alive[drop] = false; held.delete(drop);
      push(keep);
      const kg = clusters[keep].g;
      for (const h of held) { const hg = clusters[h].g; if (hg.end > kg.start && kg.end > hg.start) { held.delete(h); push(h); } }
      // the dropped group's links are the kept one's
      for (const [o, v] of bridges[drop]) { bridges[o].delete(drop); if (o !== keep) link(keep, o, v); }
      bridges[drop].clear();
    }
  }
  // every read placed again against the final groups: the first reads of a group joined it before its consensus was
  // known. Groups left with fewer than minSupport fragments are dropped and the placement redone: their reads are
  // minor patterns, or go to the groups that remain; ambiguity is between surviving groups only
  let live = clusters.map((_, i) => i).filter(i => alive[i] && clusters[i].members.reduce((s, g) => s + g.f, 0) >= minSupport);
  for (;;) {
    const got = live.map(() => [] as Group[]), index = new Map(live.map((ci, k) => [ci, k])), nearLive = indexOf(live);
    ambiguous.clear(); minor.length = 0;
    for (const g of byPos) {
      const scored = score(g, nearLive(g.start, g.end));
      if (decided(scored)) got[index.get(scored[0].i)!].push(g);
      else if (scored.length > 1) {
        const hits = scored.map(s => index.get(s.i)!).sort((a, b) => a - b), k = hits.join(',');
        const e = ambiguous.get(k); if (e) e.groups.push(g); else ambiguous.set(k, { hits, groups: [g] });
      } else minor.push(g);
    }
    const weak = new Set(live.filter((_, k) => got[k].reduce((s, g) => s + g.f, 0) < minSupport));
    if (!weak.size) {
      live.forEach((ci, k) => { anchors.push(finish({ ...clusters[ci].g, n: 0, f: 0, reads: [] })); members.push(got[k]); });
      return true;
    }
    live = live.filter(ci => !weak.has(ci));
  }
}

/** Stretches covered by at least `thr` of the reads (aligned blocks), from sorted block ends: no array over the span. */
function denseStretches(reads: AlignedRead[], thr: number): [number, number][] {
  const ev: number[] = [];
  for (const r of reads) for (const [a, b] of r.b) { ev.push(a * 2 + 1, b * 2); }
  ev.sort((x, y) => x - y);   // at one position, ends (even) before starts (odd)
  const out: [number, number][] = [];
  let d = 0, from = -1;
  for (const v of ev) {
    const pos = v >> 1;
    d += v & 1 ? 1 : -1;
    if (d >= thr && from < 0) from = pos;
    else if (d < thr && from >= 0) { if (pos > from) { const last = out[out.length - 1]; if (last && last[1] === from) last[1] = pos; else out.push([from, pos]); } from = -1; }
  }
  return out;
}

export function collapseReads(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20, minSupport = 3, minIndel = 1, tolerant = false): { sites: VariantSite[]; groups: ReadGroup[]; total: number } {
  const sites = callSites(reads, start, end, ref, refStart, minAlt, minVaf, minBq, minIndel);
  const sitePos = Int32Array.from(sites, s => s.pos);
  const snap = tolerant ? junctionSnap(reads) : null;
  const code: Record<string, string> = { ref: 'r', alt: 'a' };
  // fragments with one pattern (junctions, alleles) make one group
  const exact = new Map<string, Group>();
  for (const fr of fragmentsOf(reads)) {
    const calls = new Map<number, string>();
    const chain = new Map<number, [number, number]>();
    const blocks: [number, number][] = [];
    for (const r of fr) {
      readAlleles(r, sites, sitePos, minBq, calls);
      const rb = r.b.map(b => [b[0], b[1]] as [number, number]);
      readJunctions(r, (k, be, ns) => {
        const j = snap?.get(jKey(be, ns)) ?? [be, ns];
        chain.set(jKey(j[0], j[1]), j);
        // the blocks follow the junction's snapped ends (kept non-empty)
        if (j[0] > rb[k][0]) rb[k][1] = j[0];
        if (j[1] < rb[k + 1][1]) rb[k + 1][0] = j[1];
      });
      blocks.push(...rb);
    }
    // the mates disagreeing at a site: unknown there
    for (const [i, a] of calls) if (a === '?') calls.delete(i);
    const sorted = [...chain.values()].sort((x, y) => x[0] - y[0] || x[1] - y[1]);
    const key = `${[...chain.keys()].sort((x, y) => x - y).join(',')}|${[...calls].sort((x, y) => x[0] - y[0]).map(([i, a]) => `${i}${code[a] ?? a}`).join(',')}`;
    let g = exact.get(key);
    if (!g) { g = { chain: sorted, jkeys: new Set(chain.keys()), calls, n: 0, f: 0, reads: [], blocks: [], start: 0, end: 0, spec: 0 }; exact.set(key, g); }
    g.n += fr.length; g.f++; g.reads.push(...fr); g.blocks.push(...blocks);
  }
  const ordered = [...exact.values()].map(finish).sort((a, b) => b.spec - a.spec || b.f - a.f);
  // spliced reads (RNA-seq) are clustered like long reads (see clusterTolerant)
  let cluster = tolerant || ordered.some(g => g.chain.length);
  let anchors: Group[] = [], members: Group[][] = [];
  const ambiguous = new Map<string, { hits: number[]; groups: Group[] }>(), minor: Group[] = [];
  if (cluster && !clusterTolerant(ordered, sites, minSupport, anchors, members, ambiguous, minor)) cluster = false;
  if (!cluster) {
    // 1. seeds: the patterns of at least minSupport fragments, most specific first, that fit no seed already taken
    const seeds: Group[] = [];
    for (const g of ordered) if (g.f >= minSupport && !seeds.some(a => compatible(g, a))) seeds.push(g);
    // 2. every pattern placed against the seeds: one compatible seed takes it, several make it ambiguous, none
    //    minor (a pattern met before a less specific seed it also fits is ambiguous too); a seed left with fewer than
    //    minSupport fragments is dropped and the placement redone
    let keep = seeds;
    for (;;) {
      members = keep.map(() => []); ambiguous.clear(); minor.length = 0;
      for (const g of ordered) {
        const hits: number[] = [];
        for (let i = 0; i < keep.length; i++) if (compatible(g, keep[i])) hits.push(i);
        if (hits.length === 1) members[hits[0]].push(g);
        else if (hits.length > 1) { const k = hits.join(','); const e = ambiguous.get(k); if (e) e.groups.push(g); else ambiguous.set(k, { hits, groups: [g] }); }
        else minor.push(g);
      }
      const weak = new Set(keep.filter((_, i) => members[i].reduce((n, g) => n + g.f, 0) < minSupport));
      if (!weak.size) break;
      keep = keep.filter(k => !weak.has(k));
    }
    anchors = keep.map(k => ({ ...k, n: 0, f: 0, reads: [] }));
  }
  const total = reads.length;
  const envelope = (gs: Group[]) => {
    const all = gs.flatMap(g => g.reads);
    const n = all.length;
    return { blocks: mergeBlocks(all.flatMap(r => r.b.map(b => [b[0], b[1]] as [number, number]))), dense: denseStretches(all, Math.max(1, n * 0.25)), n };
  };
  let groups: ReadGroup[] = anchors.map((a, i) => {
    const { blocks, dense, n } = envelope([a, ...members[i]]);
    const alleles = sites.map(() => '.');
    for (const [k, v] of a.calls) alleles[k] = v;
    return { id: `H${i + 1}`, kind: 'consensus' as const, n, frac: total ? n / total : 0, chain: a.chain, alleles, blocks, dense, absorbed: n - Math.max(a.n, ...members[i].map(g => g.n)) };
  });
  const order = groups.map((_, i) => i).sort((i, j) => groups[j].n - groups[i].n);
  const rename = new Map(order.map((i, rank) => [`H${i + 1}`, `H${rank + 1}`]));
  groups = order.map(i => ({ ...groups[i], id: rename.get(groups[i].id)! }));
  const sizeOf = (gs: Group[]) => gs.reduce((s, g) => s + g.n, 0);
  for (const { hits, groups: gs } of [...ambiguous.values()].sort((a, b) => sizeOf(b.groups) - sizeOf(a.groups))) {
    const { blocks, dense, n } = envelope(gs);
    const ids = hits.map(h => rename.get(`H${h + 1}`)!);
    groups.push({ id: ids.join('|'), kind: 'ambiguous', n, frac: total ? n / total : 0, compatible: ids, chain: [], alleles: [], blocks, dense, absorbed: 0 });
  }
  if (minor.length) {
    const n = sizeOf(minor);
    groups.push({ id: 'minor', kind: 'minor', n, frac: total ? n / total : 0, patterns: minor.length, chain: [], alleles: [], blocks: mergeBlocks(minor.flatMap(g => g.blocks)), dense: [], absorbed: 0 });
  }
  return { sites, groups, total };
}
