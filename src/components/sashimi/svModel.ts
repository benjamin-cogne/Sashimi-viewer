/**
 * Structural-variant reconstruction: the rearranged chromosome implied by a set of junctions, and whether the read
 * depth agrees with it.
 *
 * A junction joins two *breakends*. A breakend is a position and the side of it the junction keeps: `tail` keeps the
 * sequence to the left of the position (the end of a segment), `head` the sequence to the right (the start of one).
 * Every junction the viewer finds has this form: a deletion-type split read joins tail(start) to head(end); a tandem
 * duplication's junction joins tail(end) to head(start); an inversion shows two junctions, tail(start)–tail(end) (the
 * reads that leave the + strand, `+/-`) and head(start)–head(end) (`-/+`).
 *
 * The chosen breakends cut the chromosome into segments A, B, C… Each rebuilt chromosome is then walked from its
 * start: at the end of a segment an unused junction leads on (entering the next segment forward at its head, or
 * backward at its tail); without one the reference continues. Segments the walk never reaches are lost from the
 * rearranged homologue. The model then expects a copy number per segment (the intact homologues plus every visit) and
 * compares it with the read depth, relative to the depth of the flanks (or to comparison samples), with a confidence
 * interval that grows as segments get short.
 */
import type { CoverageRun, StructuralEvidence, SvArc, DiscordantArc } from './types';
import type { TxModel } from './geometry';

export type BreakendSide = 'tail' | 'head';
export interface Breakend { chrom: string; pos: number; side: BreakendSide }
export type SvJunctionKind = 'deletion' | 'duplication' | 'inversion' | 'translocation';

export interface SvJunction {
  /** short label (J1, J2…) */
  id: string;
  a: Breakend; b: Breakend;
  reads: number;
  kind: SvJunctionKind;
  /** the evidence behind the reads: "split 12 · pairs 9 · CIGAR 3" */
  evidence: string;
  /** how precisely the breakpoints are placed, in bp (5 for split reads, about the insert size for pairs alone) */
  precision: number;
  /** a note when the junction's orientation was assumed rather than read (pairs only, no split read) */
  assumed?: string;
  /** the viewer's arc key it came from */
  key?: string;
}

// ======================== Thresholds ========================
// Every number here changes what the figure says; they are listed in the cartoon's "method" note.

/** Bases left out on each side of a breakpoint when the depth of a segment is measured: split and clipped reads lower
 *  the depth right at a breakpoint, and the breakpoint itself is placed within `precision`. */
export const SV_DEPTH_MARGIN_BP = 300;
/** Target number of depth bins per segment: a segment is cut into at least this many bins (down to BIN_MIN_BP). */
export const SV_DEPTH_BINS = 12;
export const SV_DEPTH_BIN_MAX_BP = 1000, SV_DEPTH_BIN_MIN_BP = 100;
/** Bins a segment needs for a verdict; fewer: "too short to assess". */
export const SV_DEPTH_MIN_BINS = 3;
/** A segment fits when |observed − expected| ≤ max(SV_FIT_MIN_COPIES, SV_FIT_Z × SE). */
export const SV_FIT_MIN_COPIES = 0.25, SV_FIT_Z = 3;
/** Mosaic fractions reported when the fitted fraction of cells lies in this range (outside, it is a plain mismatch). */
export const SV_MOSAIC_MIN = 0.1, SV_MOSAIC_MAX = 0.9;
/** Bins the flanks must hold together for the reference level to be taken from them; fewer: the whole window. */
export const SV_FLANK_MIN_BINS = 10;

// ======================== Reconstruction ========================

export interface SvSegment {
  id: number; chrom: string; letter: string;
  /** 0-based half-open; the first and last segment of a chromosome reach the window's ends (the chromosome goes on) */
  start: number; end: number;
  first: boolean; last: boolean;
  /** index along its chromosome */
  i: number;
}
export interface SvPiece { seg: SvSegment; reversed: boolean; /** how the piece was reached: the reference order, or a junction */ via: SvJunction | 'ref' | null; /** 0 for the first visit of the segment, 1 for a second copy… */ copy: number }
export interface SvDerivative { chrom: string; name: string; pieces: SvPiece[]; /** "22pter → 9qter" when it ends on another chromosome */ ends: string; rearranged: boolean }
export interface SvReconstruction {
  segs: SvSegment[]; byChr: Record<string, SvSegment[]>; chrs: string[];
  ders: SvDerivative[]; lost: SvSegment[]; warnings: string[]; junctions: SvJunction[];
}

const chrNum = (c: string) => c.replace(/^chr/i, '');

export function reconstruct(windows: { chrom: string; start: number; end: number }[], junctions: SvJunction[]): SvReconstruction {
  const chrs = windows.map(w => w.chrom);
  const segs: SvSegment[] = [];
  const byChr: Record<string, SvSegment[]> = {};
  for (const w of windows) {
    const cuts = [...new Set(junctions.flatMap(j => [j.a, j.b]).filter(e => e.chrom === w.chrom).map(e => e.pos))].filter(p => p > w.start && p < w.end).sort((x, y) => x - y);
    const bounds = [w.start, ...cuts, w.end];
    byChr[w.chrom] = [];
    for (let i = 0; i < bounds.length - 1; i++) {
      const s: SvSegment = { id: segs.length, chrom: w.chrom, letter: '', start: bounds[i], end: bounds[i + 1], first: i === 0, last: i === bounds.length - 2, i };
      segs.push(s); byChr[w.chrom].push(s);
    }
  }
  const letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  segs.forEach((s, k) => { s.letter = k < letters.length ? letters[k] : `S${k + 1}`; });
  const keyOf = (e: Breakend): string | null => {
    const list = byChr[e.chrom];
    if (!list) return null;
    const s = e.side === 'tail' ? list.find(x => x.end === e.pos) : list.find(x => x.start === e.pos);
    return s ? `${s.id}:${e.side === 'tail' ? 't' : 'h'}` : null;
  };
  const at = new Map<string, { j: SvJunction; other: string }[]>();
  const warnings: string[] = [];
  for (const j of junctions) {
    const ka = keyOf(j.a), kb = keyOf(j.b);
    if (!ka || !kb) { warnings.push(`${j.id} reaches outside the window and was left out.`); continue; }
    if (!at.has(ka)) at.set(ka, []); at.get(ka)!.push({ j, other: kb });
    if (!at.has(kb)) at.set(kb, []); at.get(kb)!.push({ j, other: ka });
  }
  for (const [k, list] of at) if (list.length > 1) warnings.push(`${list.map(x => x.j.id).join(' and ')} join the same breakend (${k.endsWith('t') ? 'end' : 'start'} of ${segs[+k.split(':')[0]].letter}): only one can be on a chromosome; the first is taken.`);
  const used = new Set<string>(), seen = new Set<number>();
  const ders: SvDerivative[] = [];
  for (const chr of chrs) {
    const s0 = byChr[chr][0];
    if (!s0 || seen.has(s0.id)) continue;
    const pieces: SvPiece[] = [];
    const visits = new Map<number, number>();
    let key = `${s0.id}:h`, via: SvPiece['via'] = null, folded: SvSegment | null = null;
    for (let guard = 0; guard < 64; guard++) {
      const [sid, end] = key.split(':');
      const seg = segs[+sid], reversed = end === 't';
      if (via === 'ref' && seen.has(seg.id)) { folded = seg; break; }
      seen.add(seg.id);
      const copy = visits.get(seg.id) ?? 0; visits.set(seg.id, copy + 1);
      pieces.push({ seg, reversed, via, copy });
      const exit = `${seg.id}:${reversed ? 'h' : 't'}`;
      const next = (at.get(exit) ?? []).find(x => !used.has(x.j.id));
      if (next) { used.add(next.j.id); via = next.j; key = next.other; continue; }
      const list = byChr[seg.chrom];
      if (!reversed) { if (seg.i + 1 < list.length) { key = `${list[seg.i + 1].id}:h`; via = 'ref'; continue; } }
      else if (seg.i > 0) { key = `${list[seg.i - 1].id}:t`; via = 'ref'; continue; }
      break;
    }
    if (folded) warnings.push(`The walk comes back to ${folded.letter} through the reference: a junction is probably missing from the selection (the other side of an inversion, say).`);
    const rearranged = pieces.some(p => p.via && p.via !== 'ref');
    const last = pieces[pieces.length - 1];
    ders.push({ chrom: chr, pieces, rearranged, name: rearranged ? `der(${chrNum(chr)})` : chrNum(chr), ends: last && last.seg.chrom !== chr ? `${chrNum(chr)}pter → ${chrNum(last.seg.chrom)}qter` : '' });
  }
  for (const j of junctions) if (!used.has(j.id) && !warnings.some(w => w.startsWith(j.id))) warnings.push(`${j.id} is on no rebuilt chromosome: its breakends were already joined by another junction.`);
  const lost = segs.filter(s => !seen.has(s.id));
  return { segs, byChr, chrs, ders, lost, warnings, junctions };
}

// ======================== Depth ========================

export interface DepthSample { name: string; coverage: CoverageRun[] }

export interface SegmentDepth {
  seg: SvSegment;
  /** visits of the segment on the rebuilt chromosomes */
  visits: number;
  /** copies expected under the chosen zygosity */
  expected: number;
  /** copies the depth shows (relative to the reference level), with its standard error; null when too short */
  observed: number | null; se: number | null;
  /** bins measured, and the effective number after the correlation of neighbouring bins */
  bins: number; nEff: number;
  fits: boolean | null;
  /** when the depth lies between the reference copies and the expected ones: the fraction of cells carrying the event that would give it */
  mosaic: { fraction: number; lo: number; hi: number } | null;
  /** the integer copy number closest to the observed depth */
  nearest: number | null;
}
export interface DepthCheck {
  segments: SegmentDepth[];
  /** 'explained': every assessable segment fits; 'unexplained': at least one does not; 'balanced': no segment's copy number is expected to change; 'none': nothing could be measured */
  verdict: 'explained' | 'unexplained' | 'balanced' | 'none';
  /** where the reference level came from */
  reference: { how: 'controls' | 'flanks' | 'window'; depth: number; controls: string[]; note: string };
  flankCopies: number; zygosity: 'het' | 'hom';
  binBp: number;
}

/** Mean depth of [start, end) from coverage runs (0-based half-open), NaN when uncovered. */
function meanDepth(runs: CoverageRun[], start: number, end: number): number {
  let t = 0, w = 0;
  // runs are sorted; a binary search for the first run reaching `start`
  let lo = 0, hi = runs.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (runs[m].end <= start) lo = m + 1; else hi = m; }
  for (let i = lo; i < runs.length && runs[i].start < end; i++) {
    const o = Math.min(runs[i].end, end) - Math.max(runs[i].start, start);
    if (o > 0) { t += runs[i].depth * o; w += o; }
  }
  // the coverage is run-length encoded with the gaps at depth 0 left out: bases no run covers count as 0
  const span = end - start;
  return span > 0 && w >= 0 ? t / span : NaN;
}
const median = (v: number[]): number => { const s = v.filter(x => Number.isFinite(x)).sort((a, b) => a - b); return s.length ? s[s.length >> 1] : NaN; };

/** The bins of a segment, inside its margins. */
function binsOf(seg: SvSegment, margin: number, binBp: number): [number, number][] {
  const a = seg.first ? seg.start : seg.start + margin, b = seg.last ? seg.end : seg.end - margin;
  const out: [number, number][] = [];
  for (let x = a; x + binBp <= b; x += binBp) out.push([x, x + binBp]);
  return out;
}

/**
 * Copy number of every segment from the read depth. Each bin's depth is divided by the sample's reference level (the
 * median bin depth of the flanks, outside every breakpoint; with comparison samples, the ratio sample / median of the
 * controls, normalised on the flanks the same way), and multiplied by the copies of the flanks. The standard error of a
 * segment's mean follows the spread of its bins and their lag-1 correlation.
 */
export function checkDepth(rec: SvReconstruction, sample: DepthSample, controls: DepthSample[], opts: { zygosity: 'het' | 'hom'; flankCopies: number; margin?: number }): DepthCheck {
  const margin = opts.margin ?? SV_DEPTH_MARGIN_BP;
  const { zygosity, flankCopies } = opts;
  const visits = new Map<number, number>(rec.segs.map(s => [s.id, 0]));
  for (const d of rec.ders) for (const p of d.pieces) visits.set(p.seg.id, (visits.get(p.seg.id) ?? 0) + 1);
  // one bin size for every segment: the shortest inner segment gets SV_DEPTH_BINS bins
  const inner = rec.segs.filter(s => !s.first && !s.last);
  const shortest = Math.min(...(inner.length ? inner : rec.segs).map(s => s.end - s.start - (s.first || s.last ? 0 : 2 * margin)));
  const binBp = Math.max(SV_DEPTH_BIN_MIN_BP, Math.min(SV_DEPTH_BIN_MAX_BP, Math.floor(shortest / SV_DEPTH_BINS)));
  const ctrlRuns = controls.filter(c => c.coverage.length);
  const ratioOf = (bin: [number, number]): number => {
    const s = meanDepth(sample.coverage, bin[0], bin[1]);
    if (!ctrlRuns.length) return s;
    const c = median(ctrlRuns.map(k => meanDepth(k.coverage, bin[0], bin[1])));
    return c > 0 ? s / c : NaN;
  };
  const perSeg = rec.segs.map(seg => binsOf(seg, margin, binBp).map(ratioOf).filter(Number.isFinite));
  // the reference level: the flanks when they hold enough bins, else every bin of the window
  const flankIdx = rec.segs.map((s, i) => (s.first || s.last ? i : -1)).filter(i => i >= 0);
  const flankVals = flankIdx.flatMap(i => perSeg[i]);
  const how: DepthCheck['reference']['how'] = ctrlRuns.length ? 'controls' : flankVals.length >= SV_FLANK_MIN_BINS ? 'flanks' : 'window';
  const level = median(flankVals.length >= SV_FLANK_MIN_BINS ? flankVals : perSeg.flat());
  const normals = zygosity === 'het' ? Math.max(0, flankCopies - 1) : 0;
  const segments: SegmentDepth[] = rec.segs.map((seg, i) => {
    const v = visits.get(seg.id) ?? 0;
    const expected = normals + v * (zygosity === 'het' ? 1 : flankCopies);
    const vals = perSeg[i];
    if (!(level > 0) || vals.length < SV_DEPTH_MIN_BINS) return { seg, visits: v, expected, observed: null, se: null, bins: vals.length, nEff: 0, fits: null, mosaic: null, nearest: null };
    const r = vals.map(x => (x / level) * flankCopies);
    const mean = r.reduce((a, b) => a + b, 0) / r.length;
    const sd = Math.sqrt(r.reduce((a, b) => a + (b - mean) ** 2, 0) / Math.max(1, r.length - 1));
    // lag-1 autocorrelation of neighbouring bins (GC waves, mappability): fewer effective bins
    let num = 0, den = 0;
    for (let k = 0; k < r.length; k++) { den += (r[k] - mean) ** 2; if (k) num += (r[k] - mean) * (r[k - 1] - mean); }
    const rho = den > 0 ? Math.max(0, Math.min(0.95, num / den)) : 0;
    const nEff = Math.max(1, (r.length * (1 - rho)) / (1 + rho));
    const se = sd / Math.sqrt(nEff);
    const tol = Math.max(SV_FIT_MIN_COPIES, SV_FIT_Z * se);
    const fits = Math.abs(mean - expected) <= tol;
    let mosaic: SegmentDepth['mosaic'] = null;
    if (!fits && expected !== flankCopies) {
      const f = (mean - flankCopies) / (expected - flankCopies);
      const df = Math.abs(se / (expected - flankCopies));
      if (f >= SV_MOSAIC_MIN && f <= SV_MOSAIC_MAX) mosaic = { fraction: f, lo: Math.max(0, f - 2 * df), hi: Math.min(1, f + 2 * df) };
    }
    return { seg, visits: v, expected, observed: mean, se, bins: r.length, nEff, fits, mosaic, nearest: Math.max(0, Math.round(mean)) };
  });
  const assessed = segments.filter(s => s.fits != null);
  const verdict: DepthCheck['verdict'] = !assessed.length ? 'none' : assessed.some(s => !s.fits) ? 'unexplained' : assessed.every(s => s.expected === flankCopies) ? 'balanced' : 'explained';
  const note = how === 'controls' ? `ratio to the median depth of ${ctrlRuns.length} comparison sample${ctrlRuns.length > 1 ? 's' : ''} (${ctrlRuns.map(c => c.name).join(', ')}), set to ${flankCopies} copies on the flanks`
    : how === 'flanks' ? `median depth of the flanks (outside every selected breakpoint), taken as ${flankCopies} copies`
    : `median depth of the whole window (the flanks hold fewer than ${SV_FLANK_MIN_BINS} bins), taken as ${flankCopies} copies`;
  return { segments, verdict, reference: { how, depth: level, controls: ctrlRuns.map(c => c.name), note }, flankCopies, zygosity, binBp };
}

// ======================== Annotation ========================

export interface SegmentContent { genes: { name: string; whole: number[]; cut: number[]; codingBases: number; frame: 'kept' | 'shifted' | 'none' }[] }
/** The exons of the gene models inside a segment (whole, or cut by its ends) and whether removing or repeating it keeps the frame. */
export function segmentContent(seg: SvSegment, models: TxModel[]): SegmentContent {
  const genes: SegmentContent['genes'] = [];
  for (const m of models) {
    if (m.end <= seg.start || m.start >= seg.end) continue;
    const whole = m.exons.filter(e => e.start >= seg.start && e.end <= seg.end).map(e => e.rank).sort((a, b) => a - b);
    const cut = m.exons.filter(e => e.start < seg.end && e.end > seg.start && !(e.start >= seg.start && e.end <= seg.end)).map(e => e.rank).sort((a, b) => a - b);
    let coding = 0;
    if (m.cdsStart != null && m.cdsEnd != null) for (const e of m.exons) if (e.start >= seg.start && e.end <= seg.end) coding += Math.max(0, Math.min(e.end, m.cdsEnd) - Math.max(e.start, m.cdsStart));
    genes.push({ name: m.geneName, whole, cut, codingBases: coding, frame: !coding ? 'none' : coding % 3 ? 'shifted' : 'kept' });
  }
  return { genes };
}
/** Where a breakpoint falls: "intron 3 of BRCA1", "exon 5 of BRCA1", or "outside the gene models". */
export function breakpointSite(pos: number, models: TxModel[]): string {
  for (const m of models) {
    if (pos < m.start || pos >= m.end) continue;
    const ex = m.exons.find(e => pos >= e.start && pos < e.end);
    if (ex) return `exon ${ex.rank} of ${m.geneName}`;
    const sorted = [...m.exons].sort((a, b) => a.start - b.start);
    for (let i = 0; i + 1 < sorted.length; i++) if (pos >= sorted[i].end && pos < sorted[i + 1].start) {
      const n = m.strand > 0 ? sorted[i].rank : sorted[i + 1].rank;
      return `intron ${n} of ${m.geneName}`;
    }
    return `${m.geneName}`;
  }
  return 'outside the gene models';
}

/** Short names of the exon runs: "exons 3–5", "exon 7". */
export const exonRun = (ranks: number[]): string => !ranks.length ? '' : ranks.length > 1 && ranks[ranks.length - 1] - ranks[0] === ranks.length - 1 ? `exons ${ranks[0]}–${ranks[ranks.length - 1]}` : `exon${ranks.length > 1 ? 's' : ''} ${ranks.join(', ')}`;

// ======================== Notation ========================

const g = (chrom: string, a: number, b: number) => `${chrom}:g.${a + 1}_${b}`;
/**
 * A draft HGVS / ISCN description of the simple events (one deletion, one tandem duplication, one inversion made of its
 * two junctions, one reciprocal translocation); a complex event is described by its rebuilt chromosomes. Uncertain
 * breakpoints (pairs only) are written with their uncertainty.
 */
export function draftNotation(rec: SvReconstruction, build: string): { text: string; kind: 'hgvs' | 'iscn' | 'none'; note: string } {
  const J = rec.junctions;
  const unc = (j: SvJunction) => (j.precision > 10 ? ` (±${j.precision} bp)` : '');
  if (J.length === 1 && J[0].kind === 'deletion') { const [a, b] = [J[0].a.pos, J[0].b.pos].sort((x, y) => x - y); return { text: `${g(J[0].a.chrom, a, b)}del${unc(J[0])}`, kind: 'hgvs', note: `${build ? `${build} · ` : ''}${(b - a).toLocaleString('en-US')} bp deleted` }; }
  if (J.length === 1 && J[0].kind === 'duplication') { const [a, b] = [J[0].a.pos, J[0].b.pos].sort((x, y) => x - y); return { text: `${g(J[0].a.chrom, a, b)}dup${unc(J[0])}`, kind: 'hgvs', note: `${build ? `${build} · ` : ''}${(b - a).toLocaleString('en-US')} bp in tandem` }; }
  if (J.length === 2 && J.every(j => j.kind === 'inversion') && !rec.lost.length) {
    const p = [...new Set(J.flatMap(j => [j.a.pos, j.b.pos]))].sort((x, y) => x - y);
    if (p.length === 2) return { text: `${g(J[0].a.chrom, p[0], p[1])}inv${unc(J[0])}`, kind: 'hgvs', note: `${build ? `${build} · ` : ''}${(p[1] - p[0]).toLocaleString('en-US')} bp inverted, both junctions seen` };
  }
  if (J.length === 2 && J.every(j => j.kind === 'translocation') && !rec.lost.length) {
    const c = [...new Set(J.flatMap(j => [j.a.chrom, j.b.chrom]))].map(chrNum).sort((a, b) => (+a || 99) - (+b || 99));
    return { text: `t(${c.join(';')})`, kind: 'iscn', note: 'reciprocal and balanced on the junctions; bands to be added by hand' };
  }
  if (!J.length) return { text: '', kind: 'none', note: '' };
  return { text: rec.ders.filter(d => d.rearranged).map(d => `${d.name}: ${d.pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ')}`).join('; '), kind: 'none', note: 'complex event: the rebuilt order of the segments; the ISCN / HGVS description is to be written by hand from the breakpoints' };
}

// ======================== From the viewer's evidence ========================

const SRC_LABEL: Record<string, string> = { cigar: 'CIGAR', split: 'split', '+/-': 'split', '-/+': 'split', clip: 'clipped', rescued: 'rescued', realigned: 'realigned', pair: 'pairs', insertion: 'tandem insertion' };
/** "split 12 · pairs 9 · CIGAR 3" from an arc's sources. */
export function evidenceText(sources: Record<string, number> | undefined, fallback: string, count: number): string {
  if (!sources) return `${fallback} ${count}`;
  const by = new Map<string, number>();
  for (const [k, n] of Object.entries(sources)) if (n > 0) { const l = SRC_LABEL[k] ?? k; by.set(l, (by.get(l) ?? 0) + n); }
  return [...by].map(([l, n]) => `${l} ${n}`).join(' · ') || `${fallback} ${count}`;
}
/** Breakpoints placed by aligned bases (split, clipped, CIGAR): to 5 bp; by pairs alone: within an insert size. */
const BASE_SOURCES = new Set(['cigar', 'split', '+/-', '-/+', 'clip', 'rescued', 'realigned', 'insertion']);
export const SV_SPLIT_PRECISION_BP = 5, SV_PAIR_PRECISION_FALLBACK_BP = 500;

/**
 * The junctions of a DNA sample's structural evidence in one window: deletions (CIGAR and deletion-type splits),
 * duplications, inversions (two junctions each, when both orientations were read) and discordant pairs by class.
 * Ordered along the chromosome and named J1, J2…; `key` is the viewer's arc key (`kind:start-end`), shared by the two
 * junctions of an inversion arc.
 */
export function junctionsOfEvidence(sv: StructuralEvidence, chrom: string, minReads: number, hidden: (j: { start: number; end: number }) => boolean): SvJunction[] {
  const out: Omit<SvJunction, 'id'>[] = [];
  const pairPrecision = sv.insertMedian ?? SV_PAIR_PRECISION_FALLBACK_BP;
  const precisionOf = (src?: Record<string, number>) => (src && Object.entries(src).some(([k, n]) => n > 0 && BASE_SOURCES.has(k)) ? SV_SPLIT_PRECISION_BP : pairPrecision);
  const be = (pos: number, side: BreakendSide): Breakend => ({ chrom, pos, side });
  const keep = (j: { start: number; end: number; count: number }) => j.count >= minReads && j.end > j.start && !hidden(j);
  const push = (kind: 'deletion' | 'split' | 'duplication' | 'inversion' | 'discordant', j: SvArc | DiscordantArc, k: SvJunctionKind, a: Breakend, b: Breakend, reads: number, assumed?: string) =>
    out.push({ a, b, reads, kind: k, evidence: evidenceText((j as SvArc).sources, kind === 'discordant' ? 'pairs' : 'reads', j.count), precision: kind === 'discordant' ? pairPrecision : precisionOf((j as SvArc).sources), assumed, key: `${kind}:${j.start}-${j.end}` });
  for (const j of sv.deletions) if (keep(j)) push('deletion', j, 'deletion', be(j.start, 'tail'), be(j.end, 'head'), j.count);
  for (const j of sv.splits) if (keep(j)) push('split', j, 'deletion', be(j.start, 'tail'), be(j.end, 'head'), j.count);
  for (const j of sv.duplications ?? []) if (keep(j)) push('duplication', j, 'duplication', be(j.end, 'tail'), be(j.start, 'head'), j.count);
  for (const j of sv.inversions ?? []) {
    if (!keep(j)) continue;
    const pm = j.sources?.['+/-'] ?? 0, mp = j.sources?.['-/+'] ?? 0;
    if (pm > 0) push('inversion', j, 'inversion', be(j.start, 'tail'), be(j.end, 'tail'), pm);
    if (mp > 0) push('inversion', j, 'inversion', be(j.start, 'head'), be(j.end, 'head'), mp);
    if (!pm && !mp) {
      const note = 'orientation not read (no split read): both junctions of a simple inversion assumed';
      push('inversion', j, 'inversion', be(j.start, 'tail'), be(j.end, 'tail'), j.count, note);
      push('inversion', j, 'inversion', be(j.start, 'head'), be(j.end, 'head'), j.count, note);
    }
  }
  for (const j of sv.discordant) {
    if (!keep(j)) continue;
    const k = j.kind ?? 'deletion';
    if (k === 'deletion') push('discordant', j, 'deletion', be(j.start, 'tail'), be(j.end, 'head'), j.count, 'from pairs: breakpoints within an insert size');
    else if (k === 'duplication') push('discordant', j, 'duplication', be(j.end, 'tail'), be(j.start, 'head'), j.count, 'from pairs: breakpoints within an insert size');
    else { const note = 'from pairs: both junctions of a simple inversion assumed, breakpoints within an insert size'; push('discordant', j, 'inversion', be(j.start, 'tail'), be(j.end, 'tail'), j.count, note); push('discordant', j, 'inversion', be(j.start, 'head'), be(j.end, 'head'), j.count, note); }
  }
  out.sort((x, y) => Math.min(x.a.pos, x.b.pos) - Math.min(y.a.pos, y.b.pos) || Math.max(x.a.pos, x.b.pos) - Math.max(y.a.pos, y.b.pos) || (x.a.side === 'tail' ? 0 : 1) - (y.a.side === 'tail' ? 0 : 1));
  return out.map((j, i) => ({ ...j, id: `J${i + 1}` }));
}
