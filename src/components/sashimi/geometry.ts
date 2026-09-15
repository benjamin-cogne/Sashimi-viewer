/**
 * Pure geometry helpers for the Sashimi viewer.
 *
 * Conventions
 * -----------
 * - All genomic coordinates in this module are 0-based, half-open (BAM convention).
 *   Ensembl transcripts (1-based inclusive) are converted once by `toTxModel`.
 * - A "virtual" axis is a monotonic, piecewise-linear transform of genomic
 *   position. The linear axis is the identity; the equal-intron axis keeps exons
 *   and flanks at their genomic length but gives every intron the same virtual
 *   length (the display convention of MISO's sashimi_plot `intron_scale` and
 *   ggsashimi `--shrink`, made strictly equal).
 * - A `Scale` maps genomic → pixel through the virtual axis, honouring the strand
 *   flip used to draw reverse-strand genes 5'→3' left to right.
 */
import type { CoverageRun, JunctionArc, TranscriptData, ExonDepth } from './types';

// ======================== Transcript model ========================

export interface Exon0 { start: number; end: number; rank: number }

export interface TxModel {
  geneName: string;
  transcriptId: string;
  translationId?: string | null;
  isMane: boolean;
  modelKind: 'mane' | 'canonical' | 'longest';
  biotype?: string;
  chrom: string;
  strand: 1 | -1;
  start: number;
  end: number;
  /** Sorted by genomic start, 0-based half-open. */
  exons: Exon0[];
  /** 0-based half-open CDS bounds, or null for non-coding transcripts. */
  cdsStart: number | null;
  cdsEnd: number | null;
}

export function toTxModel(t: TranscriptData): TxModel {
  const exons = [...t.exons]
    .map(e => ({ start: e.start - 1, end: e.end, rank: e.rank }))
    .sort((a, b) => a.start - b.start);
  const hasCds = t.cds_start != null && t.cds_end != null;
  return {
    geneName: t.gene_name,
    transcriptId: t.transcript_id,
    translationId: t.translation_id ?? null,
    isMane: !!t.is_mane_select,
    modelKind: t.model_kind ?? (t.is_mane_select ? 'mane' : 'canonical'),
    biotype: t.biotype,
    chrom: t.chrom,
    strand: t.strand < 0 ? -1 : 1,
    start: t.start - 1,
    end: t.end,
    exons,
    cdsStart: hasCds ? (t.cds_start as number) - 1 : null,
    cdsEnd: hasCds ? (t.cds_end as number) : null,
  };
}

/** Introns of a transcript model, genomic order. */
export function intronsOf(tx: TxModel): { start: number; end: number }[] {
  const out: { start: number; end: number }[] = [];
  for (let i = 0; i + 1 < tx.exons.length; i++) {
    const s = tx.exons[i].end, e = tx.exons[i + 1].start;
    if (e > s) out.push({ start: s, end: e });
  }
  return out;
}

// ======================== Virtual axis ========================

interface Segment { gStart: number; vStart: number; factor: number }

export interface VirtualAxis {
  kind: 'linear' | 'equal-intron';
  toV(pos: number): number;
  fromV(v: number): number;
  /** Virtual length given to every intron (equal-intron axis only). */
  intronV: number;
}

export const LINEAR_AXIS: VirtualAxis = {
  kind: 'linear',
  toV: p => p,
  fromV: v => v,
  intronV: 0,
};

export function median(xs: number[]): number {
  if (xs.length === 0) return 0;
  const s = [...xs].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

/**
 * Build the equal-intron axis. Every intron is drawn with the same virtual
 * length (median exon length, clamped to [80, 300] bp-equivalents) so it reads
 * like one more exon-sized block; the flanks outside the transcript are
 * compressed at the median intron's factor so they do not dwarf the gene.
 */
export function equalIntronAxis(tx: TxModel): VirtualAxis {
  const introns = intronsOf(tx);
  if (introns.length === 0) return LINEAR_AXIS;
  const intronV = Math.min(300, Math.max(80, Math.round(median(tx.exons.map(e => e.end - e.start)))));
  const flankFactor = Math.min(1, intronV / median(introns.map(i => i.end - i.start)));

  // Segments in genomic order; each maps [gStart, next.gStart) linearly with `factor`.
  const segs: Segment[] = [];
  const first = tx.exons[0];
  let v = 0;
  segs.push({ gStart: -Infinity, vStart: -Infinity, factor: flankFactor }); // left flank (anchored at first exon)
  for (let i = 0; i < tx.exons.length; i++) {
    const ex = tx.exons[i];
    segs.push({ gStart: ex.start, vStart: v, factor: 1 });
    v += ex.end - ex.start;
    if (i + 1 < tx.exons.length) {
      const next = tx.exons[i + 1];
      const len = next.start - ex.end;
      if (len > 0) {
        segs.push({ gStart: ex.end, vStart: v, factor: intronV / len });
        v += intronV;
      }
    }
  }
  const last = tx.exons[tx.exons.length - 1];
  segs.push({ gStart: last.end, vStart: v, factor: flankFactor }); // right flank

  const toV = (pos: number): number => {
    if (pos < first.start) return (pos - first.start) * flankFactor;
    // binary search for the last segment with gStart <= pos (skip the sentinel at index 0)
    let lo = 1, hi = segs.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (segs[mid].gStart <= pos) lo = mid; else hi = mid - 1;
    }
    const s = segs[lo];
    return s.vStart + (pos - s.gStart) * s.factor;
  };
  const fromV = (val: number): number => {
    if (val < 0) return first.start + val / flankFactor;
    let lo = 1, hi = segs.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (segs[mid].vStart <= val) lo = mid; else hi = mid - 1;
    }
    const s = segs[lo];
    return s.gStart + (val - s.vStart) / s.factor;
  };
  return { kind: 'equal-intron', toV, fromV, intronV };
}

// ======================== Pixel scale ========================

export interface Scale {
  axis: VirtualAxis;
  reverse: boolean;
  plotLeft: number;
  plotRight: number;
  plotWidth: number;
  vStart: number;
  vEnd: number;
  /** genomic → px */
  x(pos: number): number;
  /** px → genomic */
  invert(px: number): number;
  /** px → virtual */
  pxToV(px: number): number;
  /** virtual → px */
  vToPx(v: number): number;
}

export function makeScale(
  axis: VirtualAxis, viewStart: number, viewEnd: number,
  plotLeft: number, plotWidth: number, reverse: boolean,
): Scale {
  const vStart = axis.toV(viewStart);
  const vEnd = axis.toV(viewEnd);
  const vSpan = Math.max(1e-9, vEnd - vStart);
  const vToPx = (v: number) => {
    const frac = (v - vStart) / vSpan;
    return plotLeft + (reverse ? 1 - frac : frac) * plotWidth;
  };
  const pxToV = (px: number) => {
    let frac = (px - plotLeft) / plotWidth;
    if (reverse) frac = 1 - frac;
    return vStart + frac * vSpan;
  };
  return {
    axis, reverse, plotLeft, plotRight: plotLeft + plotWidth, plotWidth, vStart, vEnd,
    x: pos => vToPx(axis.toV(pos)),
    invert: px => axis.fromV(pxToV(px)),
    pxToV, vToPx,
  };
}

// ======================== Coverage ========================

/** Depth at a genomic position (0 outside the fetched runs). Binary search on sorted runs. */
export function depthAt(runs: CoverageRun[], pos: number): number {
  let lo = 0, hi = runs.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    const r = runs[mid];
    if (pos < r.start) hi = mid - 1;
    else if (pos >= r.end) lo = mid + 1;
    else return r.depth;
  }
  return 0;
}

/** Index of the first run whose end is > pos (i.e. first run overlapping [pos, ∞)). */
function firstRunFrom(runs: CoverageRun[], pos: number): number {
  let lo = 0, hi = runs.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (runs[mid].end <= pos) lo = mid + 1; else hi = mid;
  }
  return lo;
}

export function maxDepthIn(runs: CoverageRun[], start: number, end: number): number {
  let m = 0;
  for (let i = firstRunFrom(runs, start); i < runs.length && runs[i].start < end; i++) {
    if (runs[i].depth > m) m = runs[i].depth;
  }
  return m;
}

export interface CoveragePaths { fill: string; stroke: string }

/**
 * Step-function coverage paths, decimated to pixel resolution: runs narrower than
 * one pixel are pooled into a per-column maximum, so the path length is bounded
 * by the plot width whatever the data size, and narrow high peaks are never lost.
 */
export function buildCoveragePaths(
  runs: CoverageRun[], scale: Scale, viewStart: number, viewEnd: number,
  baseline: number, depthToY: (d: number) => number,
): CoveragePaths {
  if (runs.length === 0) return { fill: '', stroke: '' };
  // Traversal is in genomic order; on the reverse strand that is right→left in pixels.
  // A run that does not start where the previous one ended (a source that omits zero-depth
  // stretches) is a gap: the profile drops to the baseline there instead of ramping across.
  const segs: { a: number; b: number; d: number; gap: boolean }[] = [];
  let bucketCol: number | null = null, bucketMax = 0, bucketA = 0, bucketB = 0, bucketGap = false;
  const flush = () => {
    if (bucketCol !== null) { segs.push({ a: bucketA, b: bucketB, d: bucketMax, gap: bucketGap }); bucketCol = null; }
  };
  let prevEnd: number | null = null;
  for (let i = firstRunFrom(runs, viewStart); i < runs.length && runs[i].start < viewEnd; i++) {
    const r = runs[i];
    const gs = Math.max(r.start, viewStart), ge = Math.min(r.end, viewEnd);
    const gap = prevEnd !== null && r.start > prevEnd;
    prevEnd = r.end;
    const pa = scale.x(gs), pb = scale.x(ge);
    const lo = Math.min(pa, pb), hi = Math.max(pa, pb);
    if (hi - lo >= 1) {
      flush();
      segs.push({ a: pa, b: pb, d: r.depth, gap });
    } else {
      const col = Math.floor((lo + hi) / 2);
      if (col !== bucketCol) {
        flush();
        bucketCol = col; bucketMax = r.depth; bucketGap = gap;
        bucketA = scale.reverse ? col + 1 : col;
        bucketB = scale.reverse ? col : col + 1;
      } else if (r.depth > bucketMax) bucketMax = r.depth;
    }
  }
  flush();
  if (segs.length === 0) return { fill: '', stroke: '' };
  const f = (n: number) => n.toFixed(1);
  let fill = `M${f(segs[0].a)},${f(baseline)}`;
  let stroke = `M${f(segs[0].a)},${f(depthToY(segs[0].d))}`;
  let lastY = NaN, lastB = segs[0].a;
  for (const s of segs) {
    const y = depthToY(s.d);
    if (s.gap && lastY !== baseline) {
      fill += `L${f(lastB)},${f(baseline)}L${f(s.a)},${f(baseline)}`;
      stroke += `L${f(lastB)},${f(baseline)}L${f(s.a)},${f(baseline)}`;
      lastY = baseline;
    }
    if (y !== lastY) { fill += `L${f(s.a)},${f(y)}`; stroke += `L${f(s.a)},${f(y)}`; }
    fill += `L${f(s.b)},${f(y)}`; stroke += `L${f(s.b)},${f(y)}`;
    lastY = y; lastB = s.b;
  }
  fill += `L${f(segs[segs.length - 1].b)},${f(baseline)}Z`;
  return { fill, stroke };
}

// ======================== Junctions ========================

export type JunctionClass = 'canonical' | 'exon_skipping' | 'novel_donor' | 'novel_acceptor' | 'novel';

export interface JunctionInfo {
  cls: JunctionClass;
  /** Human label, e.g. "Exon 3 → Exon 5 · exon skipping (1 exon)". */
  label: string;
  /** Display rank of the exon at the junction's left / right end, when matched. */
  leftExon: number | null;
  rightExon: number | null;
}

export const JUNCTION_CLASS_LABEL: Record<JunctionClass, string> = {
  canonical: 'canonical (consecutive exons)',
  exon_skipping: 'exon skipping',
  novel_donor: 'novel donor (5′ splice site)',
  novel_acceptor: 'novel acceptor (3′ splice site)',
  novel: 'novel junction (no annotated boundary)',
};

/**
 * Classify a junction against the transcript model. The intron [s, e) is
 * annotated when `s` equals an exon end and/or `e` equals an exon start.
 * Donor/acceptor naming follows transcription direction: on the minus strand the
 * left boundary of an intron is the acceptor side.
 */
export function classifyJunction(j: JunctionArc, tx: TxModel | null): JunctionInfo {
  if (!tx || tx.exons.length === 0) return { cls: 'novel', label: 'unannotated transcript', leftExon: null, rightExon: null };
  let li = -1, ri = -1;
  for (let i = 0; i < tx.exons.length; i++) {
    if (tx.exons[i].end === j.start) li = i;
    if (tx.exons[i].start === j.end) ri = i;
  }
  const plus = tx.strand > 0;
  const leftExon = li >= 0 ? tx.exons[li].rank : null;
  const rightExon = ri >= 0 ? tx.exons[ri].rank : null;
  const fromTo = (a: number | null, b: number | null) => plus ? [a, b] : [b, a];
  const exonTxt = (n: number | null) => n == null ? '?' : `Exon ${n}`;
  if (li >= 0 && ri >= 0) {
    const [from, to] = fromTo(leftExon, rightExon);
    if (ri === li + 1) return { cls: 'canonical', label: `${exonTxt(from)} → ${exonTxt(to)} · canonical`, leftExon, rightExon };
    if (ri > li) {
      const skipped = ri - li - 1;
      return { cls: 'exon_skipping', label: `${exonTxt(from)} → ${exonTxt(to)} · skips ${skipped} exon${skipped > 1 ? 's' : ''}`, leftExon, rightExon };
    }
    return { cls: 'novel', label: `${exonTxt(from)} → ${exonTxt(to)} · unexpected order`, leftExon, rightExon };
  }
  if (li >= 0) {
    // left boundary annotated, right novel: on + the right end is the acceptor
    const cls: JunctionClass = plus ? 'novel_acceptor' : 'novel_donor';
    const [from, to] = fromTo(leftExon, null);
    return { cls, label: `${exonTxt(from)} → ${to == null ? 'novel site' : exonTxt(to)} · ${JUNCTION_CLASS_LABEL[cls]}`, leftExon, rightExon };
  }
  if (ri >= 0) {
    const cls: JunctionClass = plus ? 'novel_donor' : 'novel_acceptor';
    const [from, to] = fromTo(null, rightExon);
    return { cls, label: `${from == null ? 'novel site' : exonTxt(from)} → ${exonTxt(to)} · ${JUNCTION_CLASS_LABEL[cls]}`, leftExon, rightExon };
  }
  return { cls: 'novel', label: JUNCTION_CLASS_LABEL.novel, leftExon, rightExon };
}

export const junctionKey = (j: { start: number; end: number }) => `${j.start}-${j.end}`;

/**
 * Nesting level for each junction: 1 for arcs containing no other arc, otherwise
 * one more than the highest level among the narrower arcs they overlap. Wide arcs
 * therefore always sit above the arcs they contain and never cross them.
 */
export function layerJunctions(js: JunctionArc[]): Map<string, number> {
  const sorted = [...js].sort((a, b) => (a.end - a.start) - (b.end - b.start) || a.start - b.start);
  const levels = new Map<string, number>();
  const done: JunctionArc[] = [];
  for (const j of sorted) {
    let lvl = 1;
    for (const k of done) {
      if (k.start < j.end && j.start < k.end) lvl = Math.max(lvl, (levels.get(junctionKey(k)) || 1) + 1);
    }
    levels.set(junctionKey(j), lvl);
    done.push(j);
  }
  return levels;
}

// ======================== Arcs (cubic Bézier with flat control points) ========================

export interface ArcGeom { x1: number; y1: number; x2: number; y2: number; cy: number; d: string }

/**
 * Arc from (x1,y1) to (x2,y2) whose apex sits exactly `apexHeight` above the
 * higher endpoint. With both control points at height cy the curve's midpoint is
 * 0.125·(y1+y2) + 0.75·cy, which is solved for cy here.
 */
export function arcGeom(x1: number, y1: number, x2: number, y2: number, apexHeight: number): ArcGeom {
  const apexY = Math.min(y1, y2) - apexHeight;
  const cy = (apexY - 0.125 * (y1 + y2)) / 0.75;
  const f = (n: number) => n.toFixed(1);
  return { x1, y1, x2, y2, cy, d: `M${f(x1)},${f(y1)}C${f(x1)},${f(cy)} ${f(x2)},${f(cy)} ${f(x2)},${f(y2)}` };
}

/** y of the arc at pixel x (bisection on the monotone x(t)). */
export function arcYAtX(g: ArcGeom, x: number): number {
  const { x1, y1, x2, y2, cy } = g;
  if (x2 === x1) return Math.min(y1, y2);
  let lo = 0, hi = 1;
  const xt = (t: number) => x1 * (1 - t) * (1 - t) * (1 + 2 * t) + x2 * t * t * (3 - 2 * t);
  const inc = x2 > x1;
  for (let i = 0; i < 24; i++) {
    const mid = (lo + hi) / 2;
    if ((xt(mid) < x) === inc) lo = mid; else hi = mid;
  }
  const t = (lo + hi) / 2, u = 1 - t;
  return u * u * u * y1 + 3 * u * u * t * cy + 3 * u * t * t * cy + t * t * t * y2;
}

// ======================== Axis ticks ========================

/** "Nice" tick step: 1, 2 or 5 × 10^k, aiming for roughly `target` ticks. */
export function niceStep(span: number, target: number): number {
  const raw = span / Math.max(1, target);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const r = raw / mag;
  const m = r < 1.5 ? 1 : r < 3.5 ? 2 : r < 7.5 ? 5 : 10;
  return Math.max(1, m * mag);
}

export function niceTicks(start: number, end: number, target: number): number[] {
  const step = niceStep(end - start, target);
  const out: number[] = [];
  for (let t = Math.ceil(start / step) * step; t <= end; t += step) out.push(t);
  return out;
}

/** Round a depth up to a nice axis maximum (1, 2, 2.5, 5 × 10^k). */
export function niceMax(v: number, floor = 10): number {
  if (v <= floor) return floor;
  if (v <= 0) return floor;
  const mag = Math.pow(10, Math.floor(Math.log10(v)));
  const r = v / mag;
  const m = r <= 1 ? 1 : r <= 2 ? 2 : r <= 2.5 ? 2.5 : r <= 5 ? 5 : 10;
  return m * mag;
}

/**
 * Genomic coordinates typed in the search box: `chr17:43,094,464`, `17:43094464-43094500`,
 * `chrX:31.1-31.2 Mb` → 1-based inclusive locus (a single position has start = end), `chr`-prefixed;
 * null when the text is not coordinates (a gene symbol, an accession).
 */
export function parseLocus(text: string): { chrom: string; start: number; end: number } | null {
  const num = "\\d[\\d,_' ]*(?:\\.\\d+)?";
  const m = text.trim().match(new RegExp(`^(?:chr)?([0-9]{1,2}|X|Y|MT|M)\\s*:\\s*(${num})\\s*(bp|kb|mb)?\\s*(?:(?:-|–|—|\\.\\.)\\s*(${num})\\s*(bp|kb|mb)?)?$`, 'i'));
  if (!m) return null;
  let c = m[1].toUpperCase();
  if (c === 'MT') c = 'M';
  const scale: Record<string, number> = { bp: 1, kb: 1e3, mb: 1e6 };
  const unit = (m[3] || m[5] || '').toLowerCase();
  const toBp = (raw: string, u: string) => Math.round(Number(raw.replace(/[,_' ]/g, '')) * (scale[(u || unit).toLowerCase()] ?? 1));
  const a = toBp(m[2], m[3] || ''), b = m[4] ? toBp(m[4], m[5] || '') : a;
  if (!Number.isFinite(a) || !Number.isFinite(b) || a < 1 || b < 1) return null;
  return { chrom: `chr${c}`, start: Math.min(a, b), end: Math.max(a, b) };
}

export function formatBp(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(n % 1_000_000 === 0 ? 0 : 2)} Mb`;
  if (n >= 1000) return `${(n / 1000).toFixed(n % 1000 === 0 ? 0 : 1)} kb`;
  return `${n} bp`;
}

// ======================== Read packing (alignment track) ========================

/**
 * Greedy first-fit row packing of intervals sorted by start: each read goes to the
 * first row whose last read ends before it (with a 1 bp gap). Reads beyond
 * `maxRows` rows are not placed and counted in `hidden`. Overlap is preserved by
 * any monotone axis, so packing in genomic space is valid in equal-intron mode.
 */
export function packReads(reads: { s: number; e: number }[], maxRows: number): { rows: Int32Array; nRows: number; hidden: number } {
  const order = reads.map((_, i) => i).sort((a, b) => reads[a].s - reads[b].s || reads[a].e - reads[b].e);
  const rows = new Int32Array(reads.length).fill(-1);
  const rowEnd: number[] = [];
  let hidden = 0;
  for (const i of order) {
    const r = reads[i];
    let placed = false;
    for (let k = 0; k < rowEnd.length; k++) {
      if (rowEnd[k] < r.s) { rowEnd[k] = r.e; rows[i] = k; placed = true; break; }
    }
    if (!placed) {
      if (rowEnd.length < maxRows) { rowEnd.push(r.e); rows[i] = rowEnd.length - 1; }
      else hidden++;
    }
  }
  return { rows, nRows: rowEnd.length, hidden };
}

// ======================== HGVS cDNA coordinates (MANE Select) ========================

export interface CdnaPos {
  /** e.g. "c.123", "c.-45", "c.*12", "c.123+45", "c.124-12", or "n.123" for non-coding transcripts. */
  label: string;
  /** Numbering without the prefix, as used in r. descriptions (e.g. "123+45"). */
  bare: string;
  kind: 'cds' | 'utr5' | 'utr3' | 'intron' | 'upstream' | 'downstream' | 'noncoding';
}

interface TxLayout {
  /** exons in transcription order with the 1-based transcript coordinate of their first transcribed base */
  order: { start: number; end: number; rank: number; tStart: number }[];
  plus: boolean;
  tCdsStart: number | null; // transcript coordinate of the A of ATG
  tCdsEnd: number | null;   // transcript coordinate of the last base of the stop codon
  tLength: number;
}

function txLayout(tx: TxModel): TxLayout {
  const plus = tx.strand > 0;
  const order = (plus ? [...tx.exons] : [...tx.exons].reverse()).map(e => ({ ...e, tStart: 0 }));
  let t = 1;
  for (const e of order) { e.tStart = t; t += e.end - e.start; }
  const tPos = (pos: number): number | null => {
    for (const e of order) if (pos >= e.start && pos < e.end) return e.tStart + (plus ? pos - e.start : e.end - 1 - pos);
    return null;
  };
  let tCdsStart: number | null = null, tCdsEnd: number | null = null;
  if (tx.cdsStart != null && tx.cdsEnd != null) {
    tCdsStart = tPos(plus ? tx.cdsStart : tx.cdsEnd - 1);
    tCdsEnd = tPos(plus ? tx.cdsEnd - 1 : tx.cdsStart);
  }
  return { order, plus, tCdsStart, tCdsEnd, tLength: t - 1 };
}

function cLabel(t: number, L: TxLayout): { label: string; bare: string; kind: CdnaPos['kind'] } {
  if (L.tCdsStart == null || L.tCdsEnd == null) return { label: `n.${t}`, bare: String(t), kind: 'noncoding' };
  if (t < L.tCdsStart) { const b = `-${L.tCdsStart - t}`; return { label: `c.${b}`, bare: b, kind: 'utr5' }; }
  if (t > L.tCdsEnd) { const b = `*${t - L.tCdsEnd}`; return { label: `c.${b}`, bare: b, kind: 'utr3' }; }
  const b = String(t - L.tCdsStart + 1);
  return { label: `c.${b}`, bare: b, kind: 'cds' };
}

/**
 * HGVS cDNA coordinate of a 0-based genomic position on the transcript model
 * (https://hgvs-nomenclature.org): c.N in the CDS, c.-N / c.*N in the UTRs, intronic
 * positions as c.N+M from the closer donor or c.N-M from the closer acceptor, and
 * n.N for non-coding transcripts. Positions outside the transcript continue the
 * UTR numbering.
 */
export function cdnaPosition(pos: number, tx: TxModel): CdnaPos {
  const L = txLayout(tx);
  const { order, plus } = L;
  for (const e of order) {
    if (pos >= e.start && pos < e.end) return cLabel(e.tStart + (plus ? pos - e.start : e.end - 1 - pos), L);
  }
  // transcription-order position helpers
  const firstBase = (e: TxLayout['order'][number]) => (plus ? e.start : e.end - 1);
  const lastBase = (e: TxLayout['order'][number]) => (plus ? e.end - 1 : e.start);
  const before = (a: number, b: number) => (plus ? a < b : a > b); // a is upstream of b
  const first = order[0], last = order[order.length - 1];
  if (before(pos, firstBase(first))) {
    const d = Math.abs(firstBase(first) - pos);
    const r = cLabel(1 - d, L);
    return { ...r, kind: 'upstream' };
  }
  if (before(lastBase(last), pos)) {
    const d = Math.abs(pos - lastBase(last));
    const r = cLabel(L.tLength + d, L);
    return { ...r, kind: 'downstream' };
  }
  for (let i = 0; i + 1 < order.length; i++) {
    const a = order[i], b = order[i + 1];
    if (before(lastBase(a), pos) && before(pos, firstBase(b))) {
      const d1 = Math.abs(pos - lastBase(a)), d2 = Math.abs(firstBase(b) - pos);
      if (d1 <= d2) { const c = cLabel(a.tStart + (a.end - a.start) - 1, L); return { label: `${c.label}+${d1}`, bare: `${c.bare}+${d1}`, kind: 'intron' }; }
      const c = cLabel(b.tStart, L);
      return { label: `${c.label}-${d2}`, bare: `${c.bare}-${d2}`, kind: 'intron' };
    }
  }
  return { label: '?', bare: '?', kind: 'intron' };
}

export interface JunctionHgvs {
  /** c. position of the last exonic base before the intron on the donor (5′) side */
  donor: string;
  /** c. position of the first exonic base after the intron on the acceptor (3′) side */
  acceptor: string;
  /** Predicted RNA-level description (r.), or a plain statement when none applies */
  effect: string;
  /** Short human label, e.g. "exon 4 skipped", "cryptic donor 35 nt into intron 3" */
  summary: string;
}

/**
 * HGVS description of a splice junction against the MANE Select model. Skipped exons
 * and cryptic exonic sites give r.X_Ydel; cryptic intronic sites give
 * r.X_YinsX+1_X+n; anything else is described as the join of its two c. positions.
 */
export function junctionHgvs(j: { start: number; end: number }, tx: TxModel, others: { start: number; end: number }[] = []): JunctionHgvs {
  const plus = tx.strand > 0;
  const donorPos = plus ? j.start - 1 : j.end;       // last exonic base upstream of the intron
  const acceptorPos = plus ? j.end : j.start - 1;    // first exonic base downstream of the intron
  const donor = cdnaPosition(donorPos, tx), acceptor = cdnaPosition(acceptorPos, tx);
  // A cryptic exon is described from both of its flanking junctions: r.X_X+1insX+a_X+b (pseudo-exon bases in intronic offsets)
  const pe = pseudoExonOf(j, tx, others);
  if (pe) {
    const first = plus ? pe.upExon : pe.downExon, second = plus ? pe.downExon : pe.upExon;
    const lastUp = plus ? first.end - 1 : first.start, firstDown = plus ? second.start : second.end - 1;
    const peFirst = plus ? pe.start : pe.end - 1, peLast = plus ? pe.end - 1 : pe.start;
    const b = (pos: number) => cdnaPosition(pos, tx).bare;
    return { donor: donor.label, acceptor: acceptor.label,
      effect: `r.${b(lastUp)}_${b(firstDown)}ins${b(peFirst)}_${b(peLast)}`,
      summary: `${pe.end - pe.start}-nt cryptic exon ${tx.chrom}:${(pe.start + 1).toLocaleString()}-${pe.end.toLocaleString()} between exons ${first.rank} and ${second.rank}` };
  }
  const order = plus ? [...tx.exons] : [...tx.exons].reverse();  // transcription order
  const firstBase = (e: Exon0) => (plus ? e.start : e.end - 1);
  const lastBase = (e: Exon0) => (plus ? e.end - 1 : e.start);
  const di = order.findIndex(e => lastBase(e) === donorPos);      // exon whose canonical donor is used
  const ai = order.findIndex(e => firstBase(e) === acceptorPos);  // exon whose canonical acceptor is used
  const c = (pos: number) => cdnaPosition(pos, tx).bare;
  const step = plus ? 1 : -1;
  const generic = { donor: donor.label, acceptor: acceptor.label };

  if (di >= 0 && ai >= 0) {
    if (ai === di + 1) return { ...generic, effect: 'r.(=) canonical junction', summary: `exon ${order[di].rank} → exon ${order[ai].rank}, canonical` };
    if (ai > di + 1) {
      const skipped = order.slice(di + 1, ai);
      const from = c(firstBase(skipped[0])), to = c(lastBase(skipped[skipped.length - 1]));
      const ranks = skipped.map(e => e.rank);
      return { ...generic, effect: `r.${from}_${to}del`, summary: `exon${ranks.length > 1 ? 's' : ''} ${ranks[0]}${ranks.length > 1 ? `–${ranks[ranks.length - 1]}` : ''} skipped` };
    }
    return { ...generic, effect: `r.${donor.bare}_${acceptor.bare} join (back-splice or unexpected order)`, summary: 'unexpected exon order' };
  }
  if (ai >= 0 && di < 0) {
    // canonical acceptor, novel donor: compare with the canonical donor of the upstream exon
    const up = order[ai - 1];
    if (up) {
      const canLast = lastBase(up);
      const inExon = plus ? donorPos >= up.start && donorPos < canLast : donorPos <= up.end - 1 && donorPos > canLast;
      if (inExon) {
        const n = Math.abs(canLast - donorPos);
        return { ...generic, effect: `r.${c(donorPos + step)}_${c(canLast)}del`, summary: `cryptic donor ${n} nt inside exon ${up.rank} (${n} nt lost)` };
      }
      const inIntron = plus ? donorPos > canLast && donorPos < acceptorPos : donorPos < canLast && donorPos > acceptorPos;
      if (inIntron) {
        const n = Math.abs(donorPos - canLast);
        return { ...generic, effect: `r.${c(canLast)}_${c(firstBase(order[ai]))}ins${c(canLast)}+1_${c(canLast)}+${n}`, summary: `cryptic donor ${n} nt into intron ${up.rank} (${n} nt retained)` };
      }
    }
  }
  if (di >= 0 && ai < 0) {
    const down = order[di + 1];
    if (down) {
      const canFirst = firstBase(down);
      const inExon = plus ? acceptorPos > canFirst && acceptorPos < down.end : acceptorPos < canFirst && acceptorPos >= down.start;
      if (inExon) {
        const n = Math.abs(acceptorPos - canFirst);
        return { ...generic, effect: `r.${c(canFirst)}_${c(acceptorPos - step)}del`, summary: `cryptic acceptor ${n} nt inside exon ${down.rank} (${n} nt lost)` };
      }
      const inIntron = plus ? acceptorPos < canFirst && acceptorPos > donorPos : acceptorPos > canFirst && acceptorPos < donorPos;
      if (inIntron) {
        const n = Math.abs(canFirst - acceptorPos);
        return { ...generic, effect: `r.${c(lastBase(order[di]))}_${c(canFirst)}ins${c(canFirst)}-${n}_${c(canFirst)}-1`, summary: `cryptic acceptor ${n} nt into intron ${order[di].rank} (${n} nt retained)` };
      }
    }
  }
  const sameIntron = donor.kind === 'intron' && acceptor.kind === 'intron';
  return {
    ...generic,
    effect: `r.${donor.bare}_${acceptor.bare} join`,
    summary: sameIntron ? 'both ends intronic: pseudo-exon flank, combine with the partner junction' : 'novel junction',
  };
}

// ======================== Percent spliced-in (junction-based) ========================

export interface ExonPsi { inclusionUp: number; inclusionDown: number; exclusion: number; psi: number | null }

/**
 * Exon inclusion from junction reads: inclusion = junctions ending at the exon start or
 * starting at its end (mean of both sides), exclusion = junctions spanning the whole exon.
 * ψ = inc / (inc + exc), null when no read informs the exon.
 */
export function exonPsi(exon: { start: number; end: number }, junctions: JunctionArc[], strand: number): ExonPsi {
  let left = 0, right = 0, exclusion = 0;
  for (const j of junctions) {
    if (j.end === exon.start) left += j.count;
    else if (j.start === exon.end) right += j.count;
    else if (j.start < exon.start && j.end > exon.end) exclusion += j.count;
  }
  const [inclusionUp, inclusionDown] = strand > 0 ? [left, right] : [right, left];
  const inc = (inclusionUp + inclusionDown) / 2;
  return { inclusionUp, inclusionDown, exclusion, psi: inc + exclusion > 0 ? inc / (inc + exclusion) : null };
}

export interface ExonSiteUsage {
  /** mean of the two boundary usages (one when only one is informative); null when no junction informs the exon */
  usage: number | null;
  /** usage of the left (genomic) boundary: junctions ending exactly there / (those + junctions spanning it) */
  left: number | null;
  right: number | null;
  /** junction reads supporting inclusion at each boundary and reads bypassing it (skipping the exon or using a cryptic site) */
  leftIn: number; leftOut: number; rightIn: number; rightOut: number;
}

/**
 * Exon usage from junction counts alone (the only information GTEx exposes): at each annotated
 * boundary, the reads of the junctions that end exactly there divided by those plus the reads of
 * the junctions that span that boundary, i.e. skip the exon or land at a cryptic site inside it.
 * A constitutive exon scores 1; an exon skipped by 20 % of the transcripts scores 0.8; a cryptic
 * acceptor used in half the transcripts halves the left-side usage.
 */
export function exonSiteUsage(exon: { start: number; end: number }, junctions: JunctionArc[]): ExonSiteUsage {
  let leftIn = 0, leftOut = 0, rightIn = 0, rightOut = 0;
  for (const j of junctions) {
    if (j.end === exon.start) leftIn += j.count;
    else if (j.start < exon.start && j.end > exon.start) leftOut += j.count;
    if (j.start === exon.end) rightIn += j.count;
    else if (j.start < exon.end && j.end > exon.end) rightOut += j.count;
  }
  const left = leftIn + leftOut > 0 ? leftIn / (leftIn + leftOut) : null;
  const right = rightIn + rightOut > 0 ? rightIn / (rightIn + rightOut) : null;
  const sides = [left, right].filter((v): v is number => v != null);
  return { usage: sides.length ? sides.reduce((a, b) => a + b, 0) / sides.length : null, left, right, leftIn, leftOut, rightIn, rightOut };
}

export interface JunctionPsi { psi5: number | null; psi3: number | null; donorTotal: number; acceptorTotal: number }

/**
 * FRASER-style ψ5 and ψ3 of a junction: its reads over all junctions sharing its donor
 * (ψ5) or its acceptor (ψ3), donor/acceptor taken in transcription direction.
 */
export function junctionPsi(j: JunctionArc, junctions: JunctionArc[], strand: number): JunctionPsi {
  const donorKey = strand > 0 ? 'start' : 'end', acceptorKey = strand > 0 ? 'end' : 'start';
  let donorTotal = 0, acceptorTotal = 0;
  for (const k of junctions) {
    if (k[donorKey] === j[donorKey]) donorTotal += k.count;
    if (k[acceptorKey] === j[acceptorKey]) acceptorTotal += k.count;
  }
  return { psi5: donorTotal ? j.count / donorTotal : null, psi3: acceptorTotal ? j.count / acceptorTotal : null, donorTotal, acceptorTotal };
}

// ======================== Reading frame of a splicing anomaly ========================

export type FrameEffect = 'in' | 'out' | 'utr' | 'unknown';

export interface FrameInfo {
  frame: FrameEffect;
  /** Signed change in transcript length (bases), when derivable. */
  delta: number | null;
  /** Bases of the change that fall in the CDS. */
  cdsBases: number;
  /** Human-readable explanation. */
  text: string;
}

/**
 * Reading-frame consequence of a junction, from the same rules as `junctionHgvs`:
 * skipped exons and cryptic exonic sites delete transcript bases, cryptic intronic
 * sites insert bases. The frame is kept when the number of *coding* bases changed is
 * a multiple of three and neither the start nor the stop codon is touched; changes
 * confined to the UTRs do not affect the frame; pseudo-exon flanks and novel joins
 * cannot be judged from one junction alone.
 */
export function junctionFrame(j: { start: number; end: number }, tx: TxModel, others: { start: number; end: number }[] = []): FrameInfo {
  const unknown: FrameInfo = { frame: 'unknown', delta: null, cdsBases: 0, text: 'frame cannot be judged from this junction alone' };
  if (tx.cdsStart == null || tx.cdsEnd == null) return { frame: 'unknown', delta: null, cdsBases: 0, text: 'non-coding transcript' };
  const L = txLayout(tx);
  if (L.tCdsStart == null || L.tCdsEnd == null) return unknown;
  const pseudo = pseudoExonOf(j, tx, others);
  if (pseudo) return pseudoExonFrame(pseudo, tx, L);
  const plus = tx.strand > 0;
  const donorPos = plus ? j.start - 1 : j.end, acceptorPos = plus ? j.end : j.start - 1;
  const order = plus ? [...tx.exons] : [...tx.exons].reverse();
  const firstBase = (e: Exon0) => (plus ? e.start : e.end - 1);
  const lastBase = (e: Exon0) => (plus ? e.end - 1 : e.start);
  const tPos = (pos: number): number | null => {
    for (const e of L.order) if (pos >= e.start && pos < e.end) return e.tStart + (plus ? pos - e.start : e.end - 1 - pos);
    return null;
  };
  const di = order.findIndex(e => lastBase(e) === donorPos), ai = order.findIndex(e => firstBase(e) === acceptorPos);
  // transcript interval [tA, tB] of the bases removed (deletion) or the CDS-coordinate anchor of an insertion
  let delta: number | null = null, tA: number | null = null, tB: number | null = null, kind: 'del' | 'ins' | null = null;
  if (di >= 0 && ai >= 0 && ai > di + 1) {
    tA = L.order[di + 1].tStart; tB = L.order[ai - 1].tStart + (L.order[ai - 1].end - L.order[ai - 1].start) - 1; kind = 'del';
  } else if (ai >= 0 && di < 0 && order[ai - 1]) {
    const up = order[ai - 1], canLast = lastBase(up);
    const inExon = plus ? donorPos >= up.start && donorPos < canLast : donorPos <= up.end - 1 && donorPos > canLast;
    const inIntron = plus ? donorPos > canLast && donorPos < acceptorPos : donorPos < canLast && donorPos > acceptorPos;
    if (inExon) { tA = tPos(donorPos)! + 1; tB = tPos(canLast)!; kind = 'del'; }
    else if (inIntron) { delta = Math.abs(donorPos - canLast); tA = tPos(canLast); tB = tA; kind = 'ins'; }
  } else if (di >= 0 && ai < 0 && order[di + 1]) {
    const down = order[di + 1], canFirst = firstBase(down);
    const inExon = plus ? acceptorPos > canFirst && acceptorPos < down.end : acceptorPos < canFirst && acceptorPos >= down.start;
    const inIntron = plus ? acceptorPos < canFirst && acceptorPos > donorPos : acceptorPos > canFirst && acceptorPos < donorPos;
    if (inExon) { tA = tPos(canFirst)!; tB = tPos(acceptorPos)! - 1; kind = 'del'; }
    else if (inIntron) { delta = Math.abs(canFirst - acceptorPos); tA = tPos(canFirst)! - 1; tB = tA; kind = 'ins'; }
  } else if (di >= 0 && ai >= 0 && ai === di + 1) {
    return { frame: 'in', delta: 0, cdsBases: 0, text: 'canonical junction' };
  }
  if (kind == null || tA == null || tB == null) return unknown;
  const cs = L.tCdsStart, ce = L.tCdsEnd;
  if (kind === 'del') {
    delta = -(tB - tA + 1);
    const cdsBases = Math.max(0, Math.min(tB, ce) - Math.max(tA, cs) + 1);
    if (cdsBases === 0) return { frame: 'utr', delta, cdsBases, text: `${-delta} nt deleted in the UTR, coding frame unaffected` };
    const hitsStart = tA <= cs + 2 && tB >= cs, hitsStop = tA <= ce && tB >= ce - 2;
    if (hitsStart) return { frame: 'out', delta, cdsBases, text: `${cdsBases} coding nt deleted including the start codon` };
    if (hitsStop) return { frame: 'out', delta, cdsBases, text: `${cdsBases} coding nt deleted including the stop codon` };
    return cdsBases % 3 === 0
      ? { frame: 'in', delta, cdsBases, text: `${cdsBases} coding nt deleted, in frame (−${cdsBases / 3} aa)` }
      : { frame: 'out', delta, cdsBases, text: `${cdsBases} coding nt deleted, frameshift` };
  }
  // insertion between transcript bases tA and tA+1
  const inCds = tA >= cs && tA < ce;
  if (!inCds) return { frame: 'utr', delta, cdsBases: 0, text: `${delta} nt of intron retained in the UTR, coding frame unaffected` };
  return delta! % 3 === 0
    ? { frame: 'in', delta, cdsBases: delta!, text: `${delta} intronic nt inserted in the CDS, in frame (+${delta! / 3} aa; a premature stop is possible)` }
    : { frame: 'out', delta, cdsBases: delta!, text: `${delta} intronic nt inserted in the CDS, frameshift` };
}

/**
 * A cryptic (pseudo-)exon shows as two junctions inside one annotated intron: one leaves the
 * upstream exon for a novel acceptor, the other leaves a novel donor for the downstream exon.
 * Given one of them, find the partner among `others` and return the pseudo-exon [start, end)
 * with its flanking annotated exons (genomic order, so strand-independent).
 */
export function pseudoExonOf(j: { start: number; end: number }, tx: TxModel, others: { start: number; end: number }[]): { start: number; end: number; upExon: Exon0; downExon: Exon0 } | null {
  if (!others.length) return null;
  const ex = tx.exons;
  const endIdx = ex.findIndex(e => e.end === j.start), startIdx = ex.findIndex(e => e.start === j.end);
  if (endIdx >= 0 && startIdx >= 0) return null;                       // both ends annotated: not a pseudo-exon flank
  if (endIdx < 0 && startIdx < 0) return null;
  if (endIdx >= 0) {
    // j leaves annotated exon `endIdx` for a novel site at j.end: the partner starts after it and lands on the next exon
    const down = ex[endIdx + 1];
    if (!down || j.end >= down.start) return null;
    const k = others.find(o => o.start > j.end && o.end === down.start && o.start < down.start);
    return k ? { start: j.end, end: k.start, upExon: ex[endIdx], downExon: down } : null;
  }
  const up = ex[startIdx - 1];
  if (!up || j.start <= up.end) return null;
  const k = others.find(o => o.end < j.start && o.start === up.end && o.end > up.end);
  return k ? { start: k.end, end: j.start, upExon: up, downExon: ex[startIdx] } : null;
}

function pseudoExonFrame(pe: { start: number; end: number; upExon: Exon0; downExon: Exon0 }, tx: TxModel, L: ReturnType<typeof txLayout>): FrameInfo {
  const len = pe.end - pe.start;
  const cs = L.tCdsStart as number, ce = L.tCdsEnd as number;
  // insertion point in transcript coordinates: last base of the exon transcribed first
  const plus = tx.strand > 0;
  const first = plus ? pe.upExon : pe.downExon;
  const lay = L.order.find(e => e.start === first.start);
  const tA = lay ? lay.tStart + (lay.end - lay.start) - 1 : null;
  if (tA == null) return { frame: 'unknown', delta: len, cdsBases: 0, text: `${len}-nt cryptic exon, position unknown` };
  if (!(tA >= cs && tA < ce)) return { frame: 'utr', delta: len, cdsBases: 0, text: `${len}-nt cryptic exon inserted in the UTR, coding frame unaffected` };
  return len % 3 === 0
    ? { frame: 'in', delta: len, cdsBases: len, text: `${len}-nt cryptic exon inserted in the CDS, in frame (+${len / 3} aa; a premature stop is possible)` }
    : { frame: 'out', delta: len, cdsBases: len, text: `${len}-nt cryptic exon inserted in the CDS, frameshift` };
}

// ======================== Translation of the MANE CDS ========================

const CODON_TABLE: Record<string, string> = (() => {
  const bases = 'TCAG', aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const t: Record<string, string> = {};
  let i = 0;
  for (const a of bases) for (const b of bases) for (const c of bases) t[a + b + c] = aas[i++];
  return t;
})();
const COMPLEMENT: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };

export function translateCodon(codon: string): string {
  return CODON_TABLE[codon.toUpperCase()] ?? 'X';
}

export interface CodonBox {
  /** 1-based codon index in the protein */
  index: number;
  aa: string;
  /** genomic segments [start, end) covered by the codon (two when it spans a splice junction) */
  segments: [number, number][];
}

/**
 * Codons of the CDS whose three bases all fall inside the reference window
 * [refStart, refStart + ref.length) and that overlap [from, to). Bases are read in
 * transcription order (reverse-complemented on the minus strand).
 */
export function codonsInWindow(tx: TxModel, ref: string, refStart: number, from: number, to: number): CodonBox[] {
  if (tx.cdsStart == null || tx.cdsEnd == null) return [];
  const L = txLayout(tx);
  if (L.tCdsStart == null || L.tCdsEnd == null) return [];
  const plus = tx.strand > 0;
  // transcript coordinate → genomic position
  const genomicOf = (t: number): number | null => {
    for (const e of L.order) {
      const len = e.end - e.start;
      if (t >= e.tStart && t < e.tStart + len) return plus ? e.start + (t - e.tStart) : e.end - 1 - (t - e.tStart);
    }
    return null;
  };
  const baseAt = (g: number): string | null => {
    if (g < refStart || g >= refStart + ref.length) return null;
    const b = ref[g - refStart].toUpperCase();
    return plus ? b : (COMPLEMENT[b] ?? 'N');
  };
  const out: CodonBox[] = [];
  const cs = L.tCdsStart, ce = L.tCdsEnd;
  const nCodons = Math.floor((ce - cs + 1) / 3);
  for (let c = 0; c < nCodons; c++) {
    const gs = [0, 1, 2].map(k => genomicOf(cs + 3 * c + k));
    if (gs.some(g => g == null)) continue;
    const g3 = gs as number[];
    if (Math.max(...g3) < from || Math.min(...g3) >= to) continue;
    const bs = g3.map(baseAt);
    if (bs.some(b => b == null)) continue;
    const segments: [number, number][] = [];
    for (const g of [...g3].sort((a, b) => a - b)) {
      const last = segments[segments.length - 1];
      if (last && last[1] === g) last[1] = g + 1; else segments.push([g, g + 1]);
    }
    out.push({ index: c + 1, aa: translateCodon(bs.join('')), segments });
  }
  return out;
}

// ======================== Depth-based exon usage ========================

/** Exon intervals evaluated for usage: the coding part of each exon of a coding model (UTR halves of the
 *  terminal exons are left out: 5′ drop-off and 3′ poly(A) bias), the whole exon otherwise. */
export function usageIntervals(tx: TxModel): { rank: number; start: number; end: number; coding: boolean }[] {
  return tx.exons.map(ex => {
    if (tx.cdsStart == null || tx.cdsEnd == null) return { rank: ex.rank, start: ex.start, end: ex.end, coding: false };
    const s = Math.max(ex.start, tx.cdsStart), e = Math.min(ex.end, tx.cdsEnd);
    return e > s ? { rank: ex.rank, start: s, end: e, coding: true } : { rank: ex.rank, start: ex.start, end: ex.end, coding: false };
  });
}

/** Indices of the exons that define the gene depth for exon `idx`: the other coding exons (all others when
 *  fewer than two coding ones remain, e.g. non-coding transcripts). */
export function referenceExons(intervals: { coding: boolean }[], idx: number): number[] {
  const coding = intervals.map((iv, i) => (iv.coding && i !== idx ? i : -1)).filter(i => i >= 0);
  if (coding.length >= 2) return coding;
  return intervals.map((_, i) => i).filter(i => i !== idx);
}

export interface ExonUsage {
  /** depth(exon) / median depth of the reference exons; null when the gene has no depth */
  usage: number | null;
  /** delta-method standard deviation from the read counts: usage · √(1/n_exon + 1/n_ref) */
  sd: number | null;
  depth: number;
  refDepth: number;
  reads: number;
  refReads: number;
}

export function exonUsage(depths: ExonDepth[], idx: number, refIdx: number[]): ExonUsage {
  const d = depths[idx];
  const refs = refIdx.map(i => depths[i]).filter(Boolean);
  const refDepth = refs.length ? median(refs.map(r => r.median)) : 0;
  const refReads = refs.reduce((a, r) => a + r.reads, 0);
  if (!d || refDepth <= 0) return { usage: null, sd: null, depth: d?.median ?? 0, refDepth, reads: d?.reads ?? 0, refReads };
  const usage = d.median / refDepth;
  const sd = d.reads > 0 && refReads > 0 ? usage * Math.sqrt(1 / d.reads + 1 / refReads) : null;
  return { usage, sd, depth: d.median, refDepth, reads: d.reads, refReads };
}

export interface Cohort { n: number; median: number | null; mad: number | null }

/** Median and MAD (scaled to σ for a normal distribution) of the usage values of a cohort. */
export function usageCohort(values: number[]): Cohort {
  const v = values.filter(x => Number.isFinite(x));
  if (!v.length) return { n: 0, median: null, mad: null };
  const m = median(v);
  return { n: v.length, median: m, mad: 1.4826 * median(v.map(x => Math.abs(x - m))) };
}

/** Robust z-score of a usage against the cohort; null below five controls or when the cohort has no spread. */
export function usageZ(u: number | null, c: Cohort): number | null {
  if (u == null || c.median == null || c.mad == null || c.n < 5) return null;
  const mad = Math.max(c.mad, 0.02);   // floor: 2 % usage, so a perfectly flat cohort does not blow up
  return (u - c.median) / mad;
}

// ======================== Junction share vs the canonical alternative (rMATS-style) ========================

export interface JunctionAlternative {
  /** reads of the canonical alternative(s) this junction competes with (mean of both inclusion junctions for skipping) */
  canonical: number | null;
  /** count / (count + canonical) */
  share: number | null;
  /** what `canonical` is */
  label: string;
}

/**
 * How much a non-canonical junction is used against its canonical alternative, from junction reads only:
 * exon skipping → skipping reads vs the mean of the two inclusion junctions of the skipped block (rMATS SE);
 * cryptic donor/acceptor → the novel junction vs the canonical junction of the annotated exon on that side
 * (rMATS A5SS/A3SS). Canonical junctions have no single alternative and return nulls.
 */
export function junctionAlternative(j: JunctionArc, junctions: JunctionArc[], tx: TxModel): JunctionAlternative {
  const ex = tx.exons;
  const li = ex.findIndex(e => e.end === j.start), ri = ex.findIndex(e => e.start === j.end);
  const count = (s: number, e: number) => junctions.find(k => k.start === s && k.end === e)?.count ?? 0;
  const done = (canonical: number | null, label: string): JunctionAlternative => ({
    canonical, share: canonical == null ? null : j.count + canonical > 0 ? j.count / (j.count + canonical) : null, label,
  });
  if (li >= 0 && ri >= 0) {
    if (ri === li + 1) return done(null, 'canonical junction');
    if (ri > li + 1) {
      const inc = (count(ex[li].end, ex[li + 1].start) + count(ex[ri - 1].end, ex[ri].start)) / 2;
      return done(inc, 'mean of the two inclusion junctions of the skipped block');
    }
    return done(null, 'unexpected order');
  }
  const pair = (a: number, b: number) => (tx.strand > 0 ? `exon ${ex[a].rank} → exon ${ex[b].rank}` : `exon ${ex[b].rank} → exon ${ex[a].rank}`);
  if (li >= 0 && ex[li + 1]) return done(count(ex[li].end, ex[li + 1].start), `canonical junction ${pair(li, li + 1)}`);
  if (ri >= 0 && ex[ri - 1]) return done(count(ex[ri - 1].end, ex[ri].start), `canonical junction ${pair(ri - 1, ri)}`);
  return done(null, 'no annotated boundary');
}
