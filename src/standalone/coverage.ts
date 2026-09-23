/**
 * Exact coverage streamed out of the alignment records, and kept.
 *
 * A coverage request used to decode the reads of its window into objects, then walk them four times
 * (depth, junctions, boundary-spanning reads, spliced fraction). Past a read budget it kept one read
 * in k and scaled the counts back. Measured on deep BAMs, 75–92 % of the time goes to inflating the
 * BGZF blocks and a few percent to the records themselves, so the budget bounded the memory but not
 * the time — and it cost exactness.
 *
 * Here every aligned block of every read goes straight into two sparse histograms: where blocks start
 * and where they end. Everything the viewer draws is read back from those two:
 *
 * - depth at x = (block starts ≤ x) − (block ends ≤ x);
 * - reads spanning an exon–intron boundary with anchors [A, B) = blocks with start ≤ A and end ≥ B,
 *   which is starts(≤ A) − ends(< B) + (blocks strictly inside (A, B)); the last term only involves
 *   blocks shorter than the 16-base anchor window, kept apart in a short list;
 * - junctions (the CIGAR N operations), counted as the reads go by.
 *
 * so any window, and any list of boundaries (the ones of a gene model chosen later included), is
 * answered from the same pass. Reads that are not uniquely mapped go into a second, usually tiny,
 * layer: "unique only" is the first minus the second, so toggling it costs nothing either.
 *
 * Nothing here decodes: the caller hands over each record's start, packed CIGAR and flags.
 */
import type { BoundaryHint, BoundarySpanning, CoverageRun, JunctionArc } from '../components/sashimi/types';
import { SPAN_EXON_ANCHOR, SPAN_INTRON_ANCHOR, type RawRead } from './alignments';
import { DepthIndexBuilder, registerDepthIndex } from '../components/sashimi/geometry';

const CHUNK_BITS = 12, CHUNK = 1 << CHUNK_BITS, CHUNK_MASK = CHUNK - 1;
/** Blocks this short or shorter can sit strictly inside a boundary's anchor window (16 bases) and are kept in a list. */
const SHORT_BLOCK = SPAN_EXON_ANCHOR + SPAN_INTRON_ANCHOR;
/** Junction key: start and length packed in one safe integer (introns up to 8.4 Mb, starts up to 2^30). */
const JUNC_LEN = 8_388_608;

/** BAM CIGAR operation codes, as packed in NUMERIC_CIGAR (length << 4 | op). */
export const CIGAR_OP = { M: 0, I: 1, D: 2, N: 3, S: 4, H: 5, P: 6, EQ: 7, X: 8 } as const;
const OP_CODE: Record<string, number> = { M: 0, I: 1, D: 2, N: 3, S: 4, H: 5, P: 6, '=': 7, X: 8 };
/** A CIGAR string packed the way BAM stores it, for records (CRAM) that only give the string. */
export function packCigar(cigar: string): number[] {
  const out: number[] = [];
  let n = 0;
  for (let i = 0; i < cigar.length; i++) {
    const c = cigar.charCodeAt(i);
    if (c >= 48 && c <= 57) n = n * 10 + (c - 48);
    else { const op = OP_CODE[cigar[i]]; if (op != null) out.push(n * 16 + op); n = 0; }
  }
  return out;
}

/** Counts per genomic position, stored by 4 kb chunk only where something landed: an intron, or the gap between two panel targets, costs nothing. */
export class Hist {
  private chunks = new Map<number, Int32Array>();
  private lastIdx = -1;
  private last: Int32Array | null = null;
  private version = 0;
  /** chunk indices in order, and the total of the chunks before each: rebuilt after new counts */
  private order: Int32Array | null = null;
  private before: Float64Array | null = null;
  private orderVersion = -1;

  add(pos: number): void {
    const idx = pos >> CHUNK_BITS;
    let a = idx === this.lastIdx ? this.last : this.chunks.get(idx);
    if (!a) { a = new Int32Array(CHUNK); this.chunks.set(idx, a); }
    this.last = a; this.lastIdx = idx;
    a[pos & CHUNK_MASK]++;
    this.version++;
  }
  get bytes(): number { return this.chunks.size * CHUNK * 4; }
  /** The chunk indices in order (cached until the next count). */
  private sorted(): Int32Array {
    if (this.orderVersion !== this.version || !this.order) {
      const order = Int32Array.from(this.chunks.keys()).sort();
      const before = new Float64Array(order.length + 1);
      for (let k = 0; k < order.length; k++) {
        const a = this.chunks.get(order[k])!;
        let t = 0;
        for (let i = 0; i < CHUNK; i++) t += a[i];
        before[k + 1] = before[k] + t;
      }
      this.order = order; this.before = before; this.orderVersion = this.version;
    }
    return this.order;
  }
  /** Sum of the counts at positions < pos: the chunks before pos's by a binary search, then part of pos's own. */
  below(pos: number): number {
    const order = this.sorted(), ci = pos >> CHUNK_BITS;
    let l = 0, h = order.length;
    while (l < h) { const m = (l + h) >> 1; if (order[m] < ci) l = m + 1; else h = m; }
    let s = this.before![l];
    const a = l < order.length && order[l] === ci ? this.chunks.get(ci)! : null;
    if (a) for (let i = 0, n = pos & CHUNK_MASK; i < n; i++) s += a[i];
    return s;
  }
  /** Adds the indices of the chunks holding counts of [start, end) to `into`. */
  chunksIn(start: number, end: number, into: Set<number>): void {
    const order = this.sorted(), c0 = start >> CHUNK_BITS, c1 = (end - 1) >> CHUNK_BITS;
    let l = 0, h = order.length;
    while (l < h) { const m = (l + h) >> 1; if (order[m] < c0) l = m + 1; else h = m; }
    for (let k = l; k < order.length && order[k] <= c1; k++) into.add(order[k]);
  }
  /** Adds `sign` × the counts of chunk `ci` into `out` (CHUNK long). */
  addChunk(ci: number, out: Int32Array, sign: number): void {
    const a = this.chunks.get(ci);
    if (a) for (let i = 0; i < CHUNK; i++) out[i] += sign * a[i];
  }
  /** Adds `sign` × the counts of [start, start + out.length) into `out`. */
  addInto(out: Int32Array, start: number, sign: number): void {
    const end = start + out.length;
    for (let ci = start >> CHUNK_BITS; ci <= (end - 1) >> CHUNK_BITS; ci++) {
      const a = this.chunks.get(ci);
      if (!a) continue;
      const base = ci << CHUNK_BITS;
      const lo = Math.max(start, base), hi = Math.min(end, base + CHUNK);
      for (let p = lo; p < hi; p++) out[p - start] += sign * a[p - base];
    }
  }
}

/** Everything counted for one set of reads. */
class Counts {
  S = new Hist();
  E = new Hist();
  /** blocks of SHORT_BLOCK bases or fewer, start and end */
  shortS: number[] = [];
  shortE: number[] = [];
  junc = new Map<number, number>();
  reads = 0;
  spliced = 0;
  /** |template length| of proper pairs, for the median insert the structural evidence measures discordance against */
  inserts = new Map<number, number>();
  /** the short blocks as start × 32 + length, sorted; rebuilt after new reads (a deep RNA gene holds hundreds of thousands) */
  private sortedShort: Float64Array | null = null;

  get bytes(): number { return this.S.bytes + this.E.bytes + this.shortS.length * 24 + this.junc.size * 32 + this.inserts.size * 32; }

  shortKeys(): Float64Array {
    if (!this.sortedShort) {
      const k = new Float64Array(this.shortS.length);
      for (let i = 0; i < k.length; i++) k[i] = this.shortS[i] * 32 + (this.shortE[i] - this.shortS[i]);
      this.sortedShort = k.sort();
    }
    return this.sortedShort;
  }

  add(start: number, ops: ArrayLike<number>, insert: number): void {
    let pos = start, spliced = false;
    for (let k = 0; k < ops.length; k++) {
      const v = ops[k], len = v >>> 4, op = v & 15;
      if (op === 0 || op === 7 || op === 8) {
        this.S.add(pos); this.E.add(pos + len);
        if (len <= SHORT_BLOCK) { this.shortS.push(pos); this.shortE.push(pos + len); this.sortedShort = null; }
        pos += len;
      } else if (op === 3) {
        // a junction is the N operation itself, as STAR's SJ.out.tab, regtools and pysam count it. The gap between two
        // aligned blocks, which the previous counting used, also took in a deletion next to the intron (the junction
        // then reported bases off its true donor) and made a junction of two deletions in a row.
        const key = pos * JUNC_LEN + Math.min(len, JUNC_LEN - 1);
        this.junc.set(key, (this.junc.get(key) ?? 0) + 1);
        pos += len; spliced = true;
      } else if (op === 2) pos += len;
    }
    this.reads++;
    if (spliced) this.spliced++;
    if (insert > 0) this.inserts.set(insert, (this.inserts.get(insert) ?? 0) + 1);
  }
}

/** Short blocks strictly inside (a, b): start > a and end < b. */
function shortInside(c: Counts, a: number, b: number): number {
  const keys = c.shortKeys();
  // first key whose start is > a, i.e. key ≥ (a + 1) × 32
  const lo = (a + 1) * 32;
  let l = 0, h = keys.length;
  while (l < h) { const m = (l + h) >> 1; if (keys[m] < lo) l = m + 1; else h = m; }
  let n = 0;
  for (let i = l; i < keys.length; i++) {
    const s = Math.floor(keys[i] / 32);
    if (s >= b) break;
    if (s + (keys[i] - s * 32) < b) n++;
  }
  return n;
}

/** The reads of one part of a scan: every read goes to `all`, a read that is not uniquely mapped to `multi` as well. */
export class Layer {
  all = new Counts();
  multi = new Counts();
  /** records the structural evidence needs (split, clipped, deleted, discordant), when it was asked for */
  sv: { r: RawRead; unique: boolean; end: number }[] = [];

  get bytes(): number { return this.all.bytes + this.multi.bytes + this.sv.length * 400; }

  /**
   * One alignment. `insert` is the |template length| of a proper pair (0 otherwise). Counting
   * starts at the record's own start, wherever the caller's window is: the caller decides
   * which layer owns which reads so that each is counted once.
   */
  add(start: number, ops: ArrayLike<number>, unique: boolean, insert: number): void {
    this.all.add(start, ops, insert);
    if (!unique) this.multi.add(start, ops, insert);
  }
}

export interface CoverageSlice {
  coverage: CoverageRun[];
  junctions: JunctionArc[];
  spanning: BoundarySpanning;
  /** over every read of the layers, not only the window's: library-type evidence, the more reads the better */
  spliced: { reads: number; fraction: number };
  /** the structural records of the reads counted, the insert median of their proper pairs */
  sv: RawRead[];
  insertMedian: number | null;
}

/** A depth sweep of [start, end) under way: the chunks holding counts, the next one, and the run being extended. */
interface Sweep {
  parts: [Counts, number][]; start: number; end: number; ids: number[]; next: number;
  d: number; runStart: number; coverage: CoverageRun[]; pyr: DepthIndexBuilder; delta: Int32Array;
}
/**
 * Sweeps chunks until they are done (false) or `budgetMs` has passed (true: call again). A plain function, not the
 * async caller: the hot loop then gets the engine's full optimisation.
 */
function sweepChunks(sw: Sweep, budgetMs: number): boolean {
  const { parts, start, end, ids, coverage, pyr, delta } = sw;
  const t0 = budgetMs < Infinity ? performance.now() : 0;
  let d = sw.d, runStart = sw.runStart;
  while (sw.next < ids.length) {
    if (budgetMs < Infinity && performance.now() - t0 > budgetMs) { sw.d = d; sw.runStart = runStart; return true; }
    const ci = ids[sw.next++];
    delta.fill(0);
    for (const [c, sign] of parts) { c.S.addChunk(ci, delta, sign); c.E.addChunk(ci, delta, -sign); }
    const base = ci << CHUNK_BITS, lo = Math.max(start, base), hi = Math.min(end, base + CHUNK);
    for (let p = lo; p < hi; p++) {
      const x = delta[p - base];
      if (x === 0) continue;
      if (p > runStart) { coverage.push({ start: runStart, end: p, depth: d }); pyr.add(runStart, p, d); runStart = p; }
      d += x;
    }
  }
  sw.d = d; sw.runStart = runStart;
  return false;
}

/** Work done between two hand-backs to the browser while a slice is read (a wide, deep window holds ~10⁶ runs). */
const SLICE_STEP_MS = 8;

/**
 * Coverage, junctions and boundary-spanning counts of [start, end) from the layers that hold the
 * reads overlapping it, each read in exactly one of them (the caller's invariant). `pause`, when
 * given, is awaited every SLICE_STEP_MS of the depth sweep, so a wide window never holds the page;
 * the layers must not change meanwhile (the caller holds the sample's lock).
 */
export async function readSlice(layers: Layer[], uniqueOnly: boolean, start: number, end: number, boundaries?: BoundaryHint,
  pause?: () => Promise<void>): Promise<CoverageSlice> {
  const parts: [Counts, number][] = [];
  for (const l of layers) { parts.push([l.all, 1]); if (uniqueOnly) parts.push([l.multi, -1]); }
  /** block starts ≤ x, block ends ≤ x */
  const startsTo = (x: number) => { let n = 0; for (const [c, sign] of parts) n += sign * c.S.below(x + 1); return n; };
  const endsTo = (x: number) => { let n = 0; for (const [c, sign] of parts) n += sign * c.E.below(x + 1); return n; };

  // depth over [start, end) as runs, one per stretch of equal depth: depth(x) = depth(x − 1) + starts at x − ends at x,
  // swept over the 4 kb chunks that hold counts only. A stretch without any (an intron, the gap between two panel
  // targets, everything around a gene on a zoomed-out view) is one run whatever its length, so a slice costs what the
  // window holds, not how wide it is. The min/max pyramid the viewer draws from is built from the runs, flat stretches
  // included in one step, and attached to them.
  const coverage: CoverageRun[] = [];
  if (end > start) {
    const ids = new Set<number>();
    for (const [c] of parts) { c.S.chunksIn(start, end, ids); c.E.chunksIn(start, end, ids); }
    // an integer from the start (the chunk totals are doubles): runs then hold small integers, not boxed numbers
    const sw: Sweep = { parts, start, end, ids: [...ids].sort((a, b) => a - b), next: 0, d: (startsTo(start - 1) - endsTo(start - 1)) | 0,
      runStart: start, coverage, pyr: new DepthIndexBuilder(start, end), delta: new Int32Array(CHUNK) };
    while (sweepChunks(sw, pause ? SLICE_STEP_MS : Infinity)) await pause!();
    coverage.push({ start: sw.runStart, end, depth: sw.d });
    sw.pyr.add(sw.runStart, end, sw.d);
    registerDepthIndex(coverage, sw.pyr.finish());
  }

  // junctions overlapping the window
  const junc = new Map<number, number>();
  for (const [c, sign] of parts) for (const [k, v] of c.junc) junc.set(k, (junc.get(k) ?? 0) + sign * v);
  const junctions: JunctionArc[] = [];
  for (const [k, count] of junc) {
    if (count <= 0) continue;
    const js = Math.floor(k / JUNC_LEN), je = js + (k % JUNC_LEN);
    if (je > start && js < end) junctions.push({ start: js, end: je, count });
  }
  junctions.sort((a, b) => a.start - b.start || a.end - b.end);

  /** blocks with start ≤ a and end ≥ b */
  const covering = (a: number, b: number): number => {
    let inside = 0;
    for (const [c, sign] of parts) inside += sign * shortInside(c, a, b);
    return startsTo(a) - endsTo(b - 1) + inside;
  };
  const starts = [...new Set([...junctions.map(j => j.start), ...(boundaries?.intronStarts ?? [])].filter(p => p >= start && p < end))].sort((a, b) => a - b);
  const ends = [...new Set([...junctions.map(j => j.end), ...(boundaries?.intronEnds ?? [])].filter(p => p > start && p <= end))].sort((a, b) => a - b);
  const spanning: BoundarySpanning = { intronStart: {}, intronEnd: {} };
  // intron start p: a block over [p − exon anchor, p + intron anchor); intron end q: over [q − intron anchor, q + exon anchor)
  for (const p of starts) spanning.intronStart[p] = covering(p - SPAN_EXON_ANCHOR, p + SPAN_INTRON_ANCHOR);
  for (const q of ends) spanning.intronEnd[q] = covering(q - SPAN_INTRON_ANCHOR, q + SPAN_EXON_ANCHOR);

  let reads = 0, spliced = 0;
  for (const [c, sign] of parts) { reads += sign * c.reads; spliced += sign * c.spliced; }

  // insert median, at the same rank structuralEvidence takes it (sorted list, index length >> 1)
  const inserts = new Map<number, number>();
  for (const [c, sign] of parts) for (const [k, v] of c.inserts) inserts.set(k, (inserts.get(k) ?? 0) + sign * v);
  let total = 0;
  for (const v of inserts.values()) total += v;
  let insertMedian: number | null = null;
  if (total > 0) {
    let rank = total >> 1;
    for (const k of [...inserts.keys()].sort((a, b) => a - b)) { const v = inserts.get(k)!; if (rank < v) { insertMedian = k; break; } rank -= v; }
  }
  const sv: RawRead[] = [];
  // the records of the reads overlapping the window, as a scan of the window would have handed them over
  for (const l of layers) for (const x of l.sv) if ((!uniqueOnly || x.unique) && x.r.start < end && x.end > start) sv.push(x.r);

  return { coverage, junctions, spanning, spliced: { reads, fraction: reads ? spliced / reads : 0 }, sv, insertMedian };
}

/**
 * What one sample has counted on one chromosome: every read whose start lies in [ps, pe) in `owned`,
 * and the reads that start before ps but reach into it in `spill`. Coverage, junctions and
 * boundary counts are then exact over [ps, pe), whatever length the reads or their introns have.
 *
 * Growing the range to the right only adds owned reads. Growing it to the left moves ps: the
 * reads of the old spill are counted again as they come — owned when they start in the new tiles,
 * in the new spill when they start before them — so the old spill is dropped whole.
 */
export class CoverageState {
  owned = new Layer();
  spill = new Layer();
  ps = 0;
  pe = 0;
  lastUsed = 0;
  constructor(readonly chrom: string, readonly structural: boolean) {}
  get empty(): boolean { return this.pe <= this.ps; }
  get bytes(): number { return this.owned.bytes + this.spill.bytes; }
  covers(start: number, end: number): boolean { return !this.empty && this.ps <= start && this.pe >= end; }
  slice(uniqueOnly: boolean, start: number, end: number, boundaries?: BoundaryHint, pause?: () => Promise<void>): Promise<CoverageSlice> {
    return readSlice([this.owned, this.spill], uniqueOnly, start, end, boundaries, pause);
  }
}
