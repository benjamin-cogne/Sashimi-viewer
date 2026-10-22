/**
 * Tandem-repeat sizing from the reads of a window, for the repeat inspector (RepeatInspector.tsx).
 *
 * Long reads of an expansion rarely align through it: the aligner soft-clips the repeat, or the read stops inside it.
 * Each read is therefore rebuilt over the locus from what it carries (its soft clips, its aligned bases with their
 * mismatches and insertions) and the tract is measured between the two unique flanks of the reference (FLANK bp on
 * each side), wherever they lie in the read:
 *
 * - **spanning**: both flanks found; the repeat is what lies between them, its size in units its length / the motif
 *   length (rounded), so that sequencing indels inside the tract average out rather than break the count;
 * - **truncated**: one flank found and the read ends inside the repeat (or a read that is all repeat): a lower bound,
 *   the motif-rich stretch next to that flank.
 *
 * The tract's composition is read unit by unit in the phase of the motif: a unit of the pathogenic motif, of a benign or
 * reference motif (when they differ, as at RFC1), a known interruption (AGG at FMR1, CAA at HTT), or another unit; a
 * base or two that break the phase (a sequencing indel) are skipped to the next motif unit.
 *
 * Alleles are the modes of the sizes (callAlleles): a kernel density on the log of the size, so that a peak's width
 * grows with its size as sequencing and PCR stutter do, up to two peaks kept when each stands out of the valley between
 * them. A peak whose reads spread widely (P10–P90 over a quarter of its median) is reported as broad: somatic mosaicism,
 * or stutter of the amplification or the sequencing.
 *
 * Coordinates are 0-based half-open, as the reads'.
 */
import type { AlignedRead } from '../components/sashimi/types';
import type { StrLocus } from './strCatalog';
import { STR_FLANKS } from './strFlanks';

/** anchor length: a locus-specific stretch on each side of the tract, found in a read to place the tract's ends */
export const ANCHOR = 40;
/** how far from the tract an anchor may sit (reference bases), when the stretches next to it look like the repeat */
const ANCHOR_SEARCH = 150;
/** reference bases around the tract a locus needs (the anchors' room) */
export const FLANK = ANCHOR_SEARCH + ANCHOR;
/** reference bases shown next to a read's tract, on each side */
export const FLANK_SHOWN = 120;
/** an anchor must differ from the repeat by this share of its bases at least (FMR1's 40 bp next to the CGG tract differ by 28 %) */
const ANCHOR_MIN_FAR = 0.35;
/** edits allowed between an anchor and the read: 15 %, above nanopore error rates and well below ANCHOR_MIN_FAR */
const ANCHOR_MAX_EDITS = 6;
/** a tract whose units are this share of the motif or a known interruption at least counts; under it, impure (chimera, noise) */
export const PURITY_MIN = 0.8;
/** seeds taken in an anchor: exact matches of this length, every SEED_STEP bases */
const SEED = 10, SEED_STEP = 2;
/** spanning reads needed to size a sample from them alone (below it, the reads stopping inside count as lower bounds) */
export const SPAN_MIN = 20;

export type CategoryTone = 'normal' | 'intermediate' | 'premutation' | 'reduced' | 'pathogenic';
export interface Category { name: string; min: number; max: number; tone: CategoryTone }

export interface RepeatLocus {
  chrom: string;
  /** the tract in the reference, 0-based half-open */
  start: number; end: number;
  /** motif length and the motif as it reads at the tract's start (reference orientation) */
  k: number; motif: string;
  /** motifs in the reference orientation: pathogenic, benign (reference ones not pathogenic), known interruptions */
  pathogenic: string[]; benign: string[]; interruptions: string[];
  refUnits: number;
  /** the reference tract in units: "(CGG)10 AGG (CGG)9" */
  refStructure: string;
  /** the stretches found in the reads to place the tract's ends, and where they come from */
  anchorL: Anchor; anchorR: Anchor;
  anchorSource: string;
  /** the reference next to the tract, FLANK_SHOWN bases on each side: what a read's flanks are compared with */
  refL: string; refR: string;
  /** "FMR1 CGG" or "CAG repeat" */
  label: string;
  catalog?: StrLocus;
  categories: Category[];
  /** where the categories come from */
  categorySource?: string;
}

/** An anchor: its sequence (reference orientation), the bases between it and the tract, its share of edits from the repeat. */
export interface Anchor { seq: string; gap: number; far: number }

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
export const revComp = (s: string) => { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] ?? 'N'; return o; };
const rotations = (m: string) => Array.from({ length: m.length }, (_, i) => m.slice(i) + m.slice(0, i));
const reverse = (s: string) => s.split('').reverse().join('');
const hamming = (a: string, b: string) => { let n = 0; for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) n++; return n; };

/**
 * The tract of a motif around [hintStart, hintEnd) in the reference `ref` (starting at `refStart`): from the seed that
 * gives the longest one, units of the motif in its phase there, extended each way while no two units in a row differ
 * from it by more than one base (an interruption such as AGG in CGG); the tract ends on exact units.
 */
export function findTract(ref: string, refStart: number, hintStart: number, hintEnd: number, motifs: string[]): { start: number; end: number; unit: string } | null {
  const k = motifs[0]?.length ?? 0;
  if (!k) return null;
  // seeds on the motifs as written (CGG, not GCG): the tract is then read in their phase, where the interruptions are named
  const units = new Set(motifs);
  let best: { start: number; end: number; unit: string } | null = null;
  const lo = Math.max(0, hintStart - refStart - 2 * k), hi = Math.min(ref.length - k, hintEnd - refStart + 2 * k);
  const tried = new Set<number>();
  for (let p = lo; p <= hi; p++) {
    const u = ref.substr(p, k);
    if (!units.has(u) || tried.has(p)) continue;
    // walk right then left in steps of k
    let e = p + k, q = p + k, miss = 0;
    while (q + k <= ref.length && miss < 2) { const s = ref.substr(q, k); if (s === u) { e = q + k; miss = 0; tried.add(q); } else if (hamming(s, u) > 1) miss++; q += k; }
    let s0 = p; q = p - k; miss = 0;
    while (q >= 0 && miss < 2) { const s = ref.substr(q, k); if (s === u) { s0 = q; miss = 0; } else if (hamming(s, u) > 1) miss++; q -= k; }
    if (!best || e - s0 > best.end - best.start) best = { start: s0 + refStart, end: e + refStart, unit: ref.substr(s0, k) };
  }
  return best && best.end - best.start >= 3 * k ? best : null;
}

/** The motif of a stretch: the shortest unit (1–6 bp) whose rotations cover most of it. */
export function detectMotif(seq: string): { motif: string; coverage: number } | null {
  let best: { motif: string; coverage: number; k: number } | null = null;
  for (let k = 1; k <= 6; k++) {
    if (seq.length < 3 * k) break;
    const counts = new Map<string, number>();
    for (let i = 0; i + k <= seq.length; i++) {
      const s = seq.substr(i, k);
      if (s.includes('N')) continue;
      const key = rotations(s).sort()[0];
      counts.set(key, (counts.get(key) ?? 0) + 1);
    }
    let top = '', n = 0;
    for (const [m, c] of counts) if (c > n) { top = m; n = c; }
    // a unit made of a shorter one (CGGCGG) covers as well as it: keep the shorter
    if (k > 1 && /^(.+)\1+$/.test(top)) continue;
    const coverage = n / (seq.length - k + 1);
    if (!best || coverage > best.coverage + 0.05) best = { motif: top, coverage, k };
  }
  return best && best.coverage >= 0.5 ? { motif: best.motif, coverage: best.coverage } : null;
}

/** The reference tract in units, runs compressed: "(CGG)10 AGG (CGG)9". */
export function tractStructure(seq: string, unit: string): string {
  const k = unit.length, out: [string, number][] = [];
  for (let i = 0; i + k <= seq.length; i += k) {
    const s = seq.substr(i, k);
    if (out.length && out[out.length - 1][0] === s) out[out.length - 1][1]++; else out.push([s, 1]);
  }
  return out.map(([s, n]) => (n > 1 ? `(${s})${n}` : s)).join(' ');
}

/** A locus of the catalogue, or a stretch picked by hand, made a locus from the reference around it. */
export function makeLocus(ref: string, refStart: number, chrom: string, hint: { start: number; end: number }, entry?: StrLocus, motif?: string): RepeatLocus | null {
  const motifs = entry ? [...new Set([...entry.ref, ...entry.path, ...entry.benign])].filter(m => m.length === entry.k) : motif ? [motif] : [];
  if (!motifs.length) return null;
  const t = findTract(ref, refStart, hint.start, hint.end, motifs);
  if (!t) return null;
  const a = t.start - refStart, b = t.end - refStart;
  if (a < ANCHOR + 10 || b + ANCHOR + 10 > ref.length) return null;
  const k = t.unit.length;
  const pathogenic = entry?.path.length ? entry.path : motifs;
  const benign = entry ? [...new Set([...entry.benign, ...entry.ref])].filter(m => !pathogenic.includes(m)) : [];
  const tract = ref.slice(a, b);
  // interruptions: the catalogue's, else the units one base off the motif found in the reference tract
  const interruptions = entry?.intr.length ? entry.intr.filter(m => m.length === k) : [...new Set(Array.from({ length: Math.floor(tract.length / k) }, (_, i) => tract.substr(i * k, k)).filter(s => s !== t.unit && hamming(s, t.unit) === 1))];
  const { categories, source } = categoriesOf(entry);
  // the stretches next to the tract: published flanks where a tool gives them (STRique), else the reference's
  const pub = entry ? STR_FLANKS[entry.id] : undefined;
  const exact = [...new Set([...pathogenic, ...benign, t.unit])];
  const before = pub ? trimUnits(pub.prefix, exact, true) : '';
  const after = pub ? trimUnits(pub.suffix, exact, false) : '';
  const reps = [...new Set([...pathogenic, ...benign])].map(m => m.repeat(Math.ceil((3 * ANCHOR) / m.length)));
  const refBefore = ref.slice(Math.max(0, a - ANCHOR_SEARCH - ANCHOR), a), refAfter = ref.slice(b, b + ANCHOR_SEARCH + ANCHOR);
  // an anchor from a published flank is placed in the reference next to this tract (whose ends may differ from the
  // tool's by a phase): its gap is measured there; not found there, the reference's own stretches are used
  const placed = (anc: Anchor | null, stretch: string, isBefore: boolean): Anchor | null => {
    if (!anc) return null;
    const hits = anchorHits(stretch, anc.seq, isBefore);
    if (hits.length !== 1) return null;
    return { ...anc, gap: isBefore ? stretch.length - hits[0].pos : hits[0].pos };
  };
  let anchorL = pub ? placed(chooseAnchor(before, reps, true), refBefore, true) : null;
  let anchorR = pub ? placed(chooseAnchor(after, reps, false), refAfter, false) : null;
  const fromPub = !!(anchorL && anchorR);
  if (!fromPub) { anchorL = chooseAnchor(refBefore, reps, true); anchorR = chooseAnchor(refAfter, reps, false); }
  if (!anchorL || !anchorR) return null;
  return {
    chrom, start: t.start, end: t.end, k, motif: t.unit, pathogenic, benign, interruptions,
    refUnits: Math.round(tract.length / k), refStructure: tractStructure(tract, t.unit),
    refL: ref.slice(Math.max(0, a - FLANK_SHOWN), a), refR: ref.slice(b, b + FLANK_SHOWN),
    anchorL, anchorR, anchorSource: fromPub ? `${pub!.tool} flanks (${pub!.name}), Giesselmann et al. 2019` : 'reference flanks',
    label: entry ? `${entry.gene} ${entry.path[0] ?? t.unit}` : `${t.unit} repeat`, catalog: entry, categories, categorySource: source,
  };
}

/**
 * A flank without the repeat units at its tract end (a tool's flank may keep some outside its tract): cut after the last
 * exact motif unit reached through units one base off at most, no two misses in a row, as findTract ends the tract.
 */
function trimUnits(flank: string, motifs: string[], atEnd: boolean): string {
  const k = motifs[0]?.length ?? 3;
  const unitAt = (i: number) => (atEnd ? flank.slice(flank.length - (i + 1) * k, flank.length - i * k) : flank.slice(i * k, (i + 1) * k));
  let cut = 0, miss = 0;
  for (let i = 0; (i + 1) * k <= flank.length && miss < 2; i++) {
    const u = unitAt(i);
    if (motifs.includes(u)) { cut = i + 1; miss = 0; } else if (!motifs.some(m => hamming(m, u) <= 1)) miss++;
  }
  return atEnd ? flank.slice(0, flank.length - cut * k) : flank.slice(cut * k);
}

/**
 * The anchor of one side: the ANCHOR-bp window nearest the tract that differs from the repeat (`reps`: motifs written
 * out) by ANCHOR_MIN_FAR of its bases, else the one that differs most. `before`: the stretch ends at the tract;
 * otherwise it starts there.
 */
function chooseAnchor(stretch: string, reps: string[], before: boolean): Anchor | null {
  if (stretch.length < ANCHOR) return null;
  let best: Anchor | null = null;
  for (let gap = 0; gap + ANCHOR <= stretch.length; gap++) {
    const seq = before ? stretch.slice(stretch.length - gap - ANCHOR, stretch.length - gap) : stretch.slice(gap, gap + ANCHOR);
    if (seq.includes('N')) continue;
    const far = Math.min(...reps.map(r => fitPattern(seq, r).ed)) / ANCHOR;
    if (far >= ANCHOR_MIN_FAR) return { seq, gap, far };
    if (!best || far > best.far + 0.01) best = { seq, gap, far };
  }
  return best;
}

/**
 * The size categories of a locus, in units. FMR1 and HTT follow their laboratory standards (four classes); the others
 * STRchive's benign / intermediate / pathogenic ranges.
 */
export function categoriesOf(entry?: StrLocus): { categories: Category[]; source?: string } {
  if (!entry) return { categories: [] };
  if (entry.gene === 'FMR1') return {
    source: 'ACMG technical standard (Spector et al. 2021)',
    categories: [
      { name: 'Normal', min: 0, max: 44, tone: 'normal' }, { name: 'Intermediate', min: 45, max: 54, tone: 'intermediate' },
      { name: 'Premutation', min: 55, max: 200, tone: 'premutation' }, { name: 'Full mutation', min: 201, max: Infinity, tone: 'pathogenic' },
    ],
  };
  if (entry.gene === 'HTT') return {
    source: 'ACMG/ASHG Huntington disease testing guidelines (1998)',
    categories: [
      { name: 'Normal', min: 0, max: 26, tone: 'normal' }, { name: 'Intermediate', min: 27, max: 35, tone: 'intermediate' },
      { name: 'Reduced penetrance', min: 36, max: 39, tone: 'reduced' }, { name: 'Full penetrance', min: 40, max: Infinity, tone: 'pathogenic' },
    ],
  };
  const [, bMax, iMin, iMax, pMin] = entry.ranges;
  const out: Category[] = [];
  if (bMax != null) out.push({ name: 'Benign', min: 0, max: bMax, tone: 'normal' });
  if (iMin != null && iMax != null) out.push({ name: 'Intermediate', min: iMin, max: iMax, tone: 'intermediate' });
  if (pMin != null) out.push({ name: 'Pathogenic', min: pMin, max: Infinity, tone: 'pathogenic' });
  return { categories: out, source: 'STRchive' };
}
export const categoryOf = (cats: Category[], units: number) => cats.find(c => units >= c.min && units <= c.max) ?? (cats.length && units > cats[cats.length - 1].max ? cats[cats.length - 1] : null);

// ---------------------------------------------------------------------------------------------------------------
// Reads

/** A read's sequence on the reference strand: its soft clips, its aligned bases (the reference with its mismatches) and its insertions. */
export function readSequence(r: AlignedRead, ref: string, refStart: number): string {
  const mm = new Map<number, string>(r.m.map(([p, b]) => [p, b]));
  const ins = new Map<number, string>(r.i.map(([p, len], j) => [p, r.is?.[j] ?? 'N'.repeat(len)]));
  const parts: string[] = [];
  parts.push(r.c[0] ? (r.cs?.[0] ?? 'N'.repeat(r.c[0])) : '');
  // an insertion sits before the base at its position, or at the end of a block (before a deletion or the read's end):
  // each one is written once
  const put = (p: number) => { const s = ins.get(p); if (s) { parts.push(s); ins.delete(p); } };
  for (const [a, b] of r.b) {
    for (let p = a; p < b; p++) {
      put(p);
      parts.push(mm.get(p) ?? ref[p - refStart] ?? 'N');
    }
    put(b);
  }
  parts.push(r.c[1] ? (r.cs?.[1] ?? 'N'.repeat(r.c[1])) : '');
  return parts.join('');
}

/**
 * Semi-global edit distance of `pat` (used whole) in `text` (free ends): the fewest edits, and where the best alignment
 * ends in `text` (exclusive).
 */
function fitPattern(pat: string, text: string): { ed: number; end: number } {
  const n = text.length;
  let prev = new Array<number>(n + 1).fill(0), cur = new Array<number>(n + 1);
  for (let i = 1; i <= pat.length; i++) {
    cur[0] = i;
    const c = pat[i - 1];
    for (let j = 1; j <= n; j++) cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (text[j - 1] === c ? 0 : 1));
    [prev, cur] = [cur, prev];
  }
  let ed = Infinity, end = 0;
  for (let j = 0; j <= n; j++) if (prev[j] < ed) { ed = prev[j]; end = j; }
  return { ed, end };
}

/**
 * Every place an anchor sits in a read: where it ends (`before` the tract) or starts (after it), with its edits. Exact
 * seeds of the anchor propose places; each is checked by aligning the whole anchor there, ANCHOR_MAX_EDITS at most: a
 * GC-rich seed also matches inside a GC-rich repeat, where the anchor itself is not.
 */
function anchorHits(seq: string, anchor: string, before: boolean): { pos: number; ed: number }[] {
  const at: number[] = [];
  for (let o = 0; o + SEED <= anchor.length; o += SEED_STEP) {
    const key = anchor.substr(o, SEED);
    for (let j = seq.indexOf(key); j >= 0; j = seq.indexOf(key, j + 1)) at.push(j - o);
  }
  at.sort((a, b) => a - b);
  const slack = 8, out: { pos: number; ed: number }[] = [];
  for (let i = 0; i < at.length;) {
    let j = i; while (j + 1 < at.length && at[j + 1] - at[i] <= slack) j++;
    const start = Math.max(0, at[(i + j) >> 1] - slack), text = seq.slice(start, start + anchor.length + 2 * slack);
    const f = before ? fitPattern(anchor, text) : fitPattern(reverse(anchor), reverse(text));
    const pos = before ? start + f.end : start + text.length - f.end;
    if (f.ed <= ANCHOR_MAX_EDITS && !out.some(h => Math.abs(h.pos - pos) < anchor.length / 2)) out.push({ pos, ed: f.ed });
    i = j + 1;
  }
  return out;
}

/**
 * The units of a tract in the motif's phase: P pathogenic motif, B benign or reference motif, I known interruption,
 * o another unit; a base or two breaking the phase (an indel) are skipped to the next motif unit.
 */
export function tokenize(tract: string, locus: RepeatLocus): string { return units(tract, locus).tok; }

/** The units of a tract (tokenize) with where each one ends in it. */
function units(tract: string, locus: RepeatLocus): { tok: string; ends: number[] } {
  const k = locus.k;
  // in the motif's phase only: a unit of another phase is a shifted read, brought back by the resync below
  const P = new Set(locus.pathogenic), B = new Set(locus.benign), I = new Set(locus.interruptions);
  let tok = '';
  const ends: number[] = [];
  for (let i = 0; i + k <= tract.length;) {
    const s = tract.substr(i, k);
    const c = P.has(s) ? 'P' : B.has(s) ? 'B' : I.has(s) ? 'I' : '';
    if (c) { tok += c; i += k; ends.push(i); continue; }
    let j = 1;
    while (j < k && !P.has(tract.substr(i + j, k)) && !B.has(tract.substr(i + j, k))) j++;
    if (j < k) { i += j; continue; }
    tok += 'o'; i += k; ends.push(i);
  }
  return { tok, ends };
}

/** units of the window over which a run's motif share is judged, and the share it must keep */
const RUN_WINDOW = 10, RUN_SHARE = 0.6;
/**
 * The motif-rich stretch at the start of `s`: its units and bases, up to the last motif unit before the share of motif
 * units over the last RUN_WINDOW drops under RUN_SHARE (sequencing errors in a long GC-rich repeat break a stricter run).
 */
function motifRun(s: string, locus: RepeatLocus): { units: number; bases: number } {
  const { tok, ends } = units(s, locus);
  let last = -1, inWindow = 0;
  for (let i = 0; i < tok.length; i++) {
    const good = tok[i] !== 'o';
    if (good) inWindow++;
    if (i >= RUN_WINDOW && tok[i - RUN_WINDOW] !== 'o') inWindow--;
    if (i + 1 >= Math.min(RUN_WINDOW, 3) && inWindow < RUN_SHARE * Math.min(RUN_WINDOW, i + 1)) break;
    if (good) last = i;
  }
  return { units: last + 1, bases: last >= 0 ? ends[last] : 0 };
}
/** The same locus read backwards (its motifs reversed): the run before a right flank, measured from the flank. */
const backwards = (locus: RepeatLocus): RepeatLocus => ({ ...locus, pathogenic: locus.pathogenic.map(reverse), benign: locus.benign.map(reverse), interruptions: locus.interruptions.map(reverse) });

export interface ReadRepeat {
  name: string;
  /** size in units (a lower bound when truncated) */
  units: number;
  /** the units' classes, from the 5′ flank of the reference (P B I o) */
  tokens: string;
  truncated: boolean;
  /** truncated: the anchor it was measured from */
  from?: 'left' | 'right' | 'none';
  reverse: boolean;
  /** the anchors' edits (absent: not found) and the share of motif or interruption units of the tract */
  edL?: number; edR?: number; purity: number;
  /** the read's sequence as it was measured (reference orientation; reverse-complemented when it read the locus backwards) */
  seq: string;
  /** the tract in `seq` (a lower bound: the stretch measured) */
  tract: [number, number];
  /** every place each anchor was found in `seq`: where the 5′ one ends, where the 3′ one starts */
  hitsL: number[]; hitsR: number[];
  /** a read set apart: why */
  note?: string;
}
export interface RepeatReads {
  /** both anchors, once each and in order, and a tract of PURITY_MIN or more */
  spanning: ReadRepeat[];
  /** one anchor: a lower bound */
  truncated: ReadRepeat[];
  /** both anchors but an impure tract (a chimera, or a read too noisy to size) */
  impure: ReadRepeat[];
  /** an anchor found twice, or the 3′ one before the 5′ one: concatemers and fold-back reads */
  chimeric: ReadRepeat[];
  /** no anchor and not all repeat: the read does not reach the locus */
  skipped: number;
  total: number;
}

const purityOf = (tok: string) => (tok.length ? [...tok].filter(c => c !== 'o').length / tok.length : 1);

/**
 * The repeat of each read (primary alignments only: a supplementary part is the same read again), from the anchors: the
 * tract lies between the 5′ anchor's end plus its gap and the 3′ anchor's start minus its gap.
 */
export function measureReads(reads: AlignedRead[], ref: string, refStart: number, locus: RepeatLocus): RepeatReads {
  const out: RepeatReads = { spanning: [], truncated: [], impure: [], chimeric: [], skipped: 0, total: 0 };
  const k = locus.k, back = backwards(locus), gL = locus.anchorL.gap, gR = locus.anchorR.gap;
  for (const r of reads) {
    if (r.f & 0x900) continue;
    out.total++;
    const seq = readSequence(r, ref, refStart);
    const base = { name: r.n, reverse: r.r === 1 };
    let done = false;
    // the read as stored, then reverse-complemented (a read folded back on itself reads the locus the other way)
    for (const s of [seq, revComp(seq)]) {
      const Ls = anchorHits(s, locus.anchorL.seq, true), Rs = anchorHits(s, locus.anchorR.seq, false);
      if (!Ls.length && !Rs.length) continue;
      done = true;
      const at = { seq: s, hitsL: Ls.map(h => h.pos), hitsR: Rs.map(h => h.pos) };
      if (Ls.length > 1 || Rs.length > 1 || (Ls.length && Rs.length && Rs[0].pos < Ls[0].pos)) {
        out.chimeric.push({ ...base, ...at, units: 0, tokens: '', truncated: false, purity: 0, tract: [0, 0], note: Ls.length > 1 || Rs.length > 1 ? 'an anchor found twice' : "the 3′ anchor before the 5′ one" });
        break;
      }
      const L = Ls[0], R = Rs[0];
      if (L && R) {
        const a = L.pos + gL, b = Math.max(a, R.pos - gR), tract = s.slice(a, b), tokens = tokenize(tract, locus), purity = purityOf(tokens);
        const rr: ReadRepeat = { ...base, ...at, tract: [a, b], units: Math.round(tract.length / k), tokens, truncated: false, edL: L.ed, edR: R.ed, purity };
        (purity >= PURITY_MIN ? out.spanning : out.impure).push(rr);
        break;
      }
      // one anchor and not the other: the repeat is at least the motif-rich stretch next to it, whatever follows (the
      // read's end, or sequence that is not the other anchor)
      if (L) {
        const a = L.pos + gL, run = motifRun(s.slice(a), locus), tokens = tokenize(s.slice(a, a + run.bases), locus);
        if (run.units >= 3) out.truncated.push({ ...base, ...at, tract: [a, a + run.bases], units: Math.round(run.bases / k), tokens, truncated: true, from: 'left', edL: L.ed, purity: purityOf(tokens) });
        else out.skipped++;
      } else {
        const b = R.pos - gR, run = motifRun(reverse(s.slice(0, Math.max(0, b))), back), tokens = tokenize(s.slice(b - run.bases, b), locus);
        if (run.units >= 3) out.truncated.push({ ...base, ...at, tract: [b - run.bases, b], units: Math.round(run.bases / k), tokens, truncated: true, from: 'right', edR: R.ed, purity: purityOf(tokens) });
        else out.skipped++;
      }
      break;
    }
    if (!done) {
      // a read that is all repeat (inside a long expansion)
      const run = motifRun(seq, locus);
      if (run.units >= 10 && run.bases >= 0.8 * seq.length) { const tokens = tokenize(seq.slice(0, run.bases), locus); out.truncated.push({ ...base, seq, hitsL: [], hitsR: [], tract: [0, run.bases], units: Math.round(run.bases / k), tokens, truncated: true, from: 'none', purity: purityOf(tokens) }); }
      else out.skipped++;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Showing a read

/** A tract as runs of units, as read in the motif's phase: `cls` P B I o, or x for bases skipped to regain the phase. */
export function unitRuns(tract: string, locus: RepeatLocus): { unit: string; cls: 'P' | 'B' | 'I' | 'o' | 'x'; n: number }[] {
  const k = locus.k;
  const P = new Set(locus.pathogenic), B = new Set(locus.benign), I = new Set(locus.interruptions);
  const out: { unit: string; cls: 'P' | 'B' | 'I' | 'o' | 'x'; n: number }[] = [];
  const push = (unit: string, cls: 'P' | 'B' | 'I' | 'o' | 'x') => {
    const last = out[out.length - 1];
    if (last && last.unit === unit && last.cls === cls && cls !== 'x') last.n++; else out.push({ unit, cls, n: 1 });
  };
  let i = 0;
  for (; i + k <= tract.length;) {
    const u = tract.substr(i, k);
    const c = P.has(u) ? 'P' : B.has(u) ? 'B' : I.has(u) ? 'I' : '';
    if (c) { push(u, c); i += k; continue; }
    let j = 1;
    while (j < k && !P.has(tract.substr(i + j, k)) && !B.has(tract.substr(i + j, k))) j++;
    if (j < k) { push(tract.substr(i, j), 'x'); i += j; continue; }
    push(u, 'o'); i += k;
  }
  if (i < tract.length) push(tract.slice(i), 'x');
  return out;
}

/**
 * A read's flank against the reference's, base to base: an alignment fixed at the tract (`atEnd`: both strings end
 * there; else both start there) and free at the other end, where the read may stop early or go on. Returns the two
 * rows with '-' for gaps, the tract side last (atEnd) or first.
 */
export function alignFlank(read: string, refFlank: string, atEnd: boolean): { ref: string; read: string } {
  const a = atEnd ? reverse(refFlank) : refFlank, b = atEnd ? reverse(read) : read;
  const n = a.length, m = b.length;
  const D: Uint16Array[] = Array.from({ length: n + 1 }, () => new Uint16Array(m + 1));
  for (let i = 0; i <= n; i++) D[i][0] = i;
  for (let j = 0; j <= m; j++) D[0][j] = j;
  for (let i = 1; i <= n; i++) for (let j = 1; j <= m; j++) D[i][j] = Math.min(D[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1), D[i - 1][j] + 1, D[i][j - 1] + 1);
  // the free end: the best cell of the last row or column
  let bi = n, bj = m, best = D[n][m];
  for (let j = 0; j <= m; j++) if (D[n][j] < best) { best = D[n][j]; bi = n; bj = j; }
  for (let i = 0; i <= n; i++) if (D[i][m] < best) { best = D[i][m]; bi = i; bj = m; }
  let ra = '', rb = '';
  let i = bi, j = bj;
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && D[i][j] === D[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)) { ra = a[i - 1] + ra; rb = b[j - 1] + rb; i--; j--; }
    else if (i > 0 && D[i][j] === D[i - 1][j] + 1) { ra = a[i - 1] + ra; rb = '-' + rb; i--; }
    else { ra = '-' + ra; rb = b[j - 1] + rb; j--; }
  }
  // beyond the free end: the rest of the longer one, against nothing
  ra += a.slice(bi) + ' '.repeat(Math.max(0, (m - bj) - (n - bi)));
  rb += b.slice(bj) + ' '.repeat(Math.max(0, (n - bi) - (m - bj)));
  return atEnd ? { ref: reverse(ra), read: reverse(rb) } : { ref: ra, read: rb };
}

// ---------------------------------------------------------------------------------------------------------------
// Alleles

export interface Allele {
  /** the peak, in units */
  mode: number;
  n: number; median: number; p5: number; p95: number; p10: number; p90: number;
  /** spread wide for its size: P10–P90 over a quarter of the median */
  broad: boolean;
  /** the reads' sizes are in [lo, hi] */
  lo: number; hi: number;
}
const quantile = (sorted: number[], p: number) => sorted.length ? sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(p * (sorted.length - 1))))] : NaN;
export { quantile };

/**
 * Up to `max` alleles from the sizes: modes of a Gaussian kernel density on log(size) (Silverman's bandwidth, 0.04 at
 * least), a peak kept when it reaches 15 % of the highest and the valley between it and a higher one drops under 70 %
 * of it; the reads split at the valleys.
 */
export function callAlleles(units: number[], max = 2): Allele[] {
  const v = units.filter(u => u > 0).sort((a, b) => a - b);
  if (v.length < 3) return [];
  const xs = v.map(u => Math.log(u));
  const mean = xs.reduce((a, b) => a + b, 0) / xs.length;
  const sd = Math.sqrt(xs.reduce((a, b) => a + (b - mean) ** 2, 0) / xs.length);
  const iqr = Math.log(quantile(v, 0.75)) - Math.log(quantile(v, 0.25));
  const h = Math.max(0.04, 0.9 * Math.min(sd, iqr / 1.34 || sd) * Math.pow(xs.length, -0.2));
  const g0 = Math.log(Math.max(1, v[0] * 0.8)), g1 = Math.log(v[v.length - 1] * 1.25), step = 0.005;
  const grid: number[] = [], dens: number[] = [];
  for (let g = g0; g <= g1; g += step) grid.push(g);
  for (const g of grid) { let d = 0; for (const x of xs) { const z = (g - x) / h; if (z > -4 && z < 4) d += Math.exp(-0.5 * z * z); } dens.push(d); }
  const top = Math.max(...dens);
  const peaks: number[] = [];
  for (let i = 1; i < dens.length - 1; i++) if (dens[i] >= dens[i - 1] && dens[i] > dens[i + 1] && dens[i] >= 0.15 * top) peaks.push(i);
  // prominence: the valley towards a higher peak
  const kept = peaks.filter(i => {
    const higher = peaks.filter(j => dens[j] > dens[i]);
    if (!higher.length) return true;
    return higher.every(j => { const [a, b] = i < j ? [i, j] : [j, i]; let m = Infinity; for (let t = a; t <= b; t++) m = Math.min(m, dens[t]); return m < 0.7 * dens[i]; });
  }).sort((a, b) => dens[b] - dens[a]).slice(0, max).sort((a, b) => a - b);
  // split at the lowest point between consecutive peaks
  const cuts: number[] = [];
  for (let p = 0; p + 1 < kept.length; p++) { let m = kept[p]; for (let t = kept[p]; t <= kept[p + 1]; t++) if (dens[t] < dens[m]) m = t; cuts.push(Math.exp(grid[m])); }
  return kept.map((i, p) => {
    const lo = p === 0 ? 0 : cuts[p - 1], hi = p === kept.length - 1 ? Infinity : cuts[p];
    const mine = v.filter(u => u > lo && u <= hi);
    const median = quantile(mine, 0.5), p10 = quantile(mine, 0.1), p90 = quantile(mine, 0.9);
    return { mode: Math.round(Math.exp(grid[i])), n: mine.length, median, p5: quantile(mine, 0.05), p95: quantile(mine, 0.95), p10, p90, broad: mine.length >= 10 && (p90 - p10) / median > 0.25, lo, hi };
  });
}

/** The interruptions of a read's first units, as positions (1-based units): "10,20". */
export function interruptionPattern(tokens: string, within = 60): string {
  const at: number[] = [];
  for (let i = 0; i < Math.min(tokens.length, within); i++) if (tokens[i] === 'I') at.push(i + 1);
  return at.join(',');
}
