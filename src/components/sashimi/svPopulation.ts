/**
 * Copy-number variants of the general population against the losses and gains of a rebuilt SV: the DGV Gold Standard
 * (UCSC track dgvGold), DGV's curated set of CNV regions found in at least two studies and two samples, the variants of
 * a region (≥ 70 % reciprocal overlap) merged into one record. A record has an outer extent (the widest of its
 * variants: chromStart-chromEnd) and an inner one (the narrowest: thickStart-thickEnd), drawn by UCSC as a box plot.
 *
 * Only the geometry (BED columns) is relied upon. Type, frequency, number of samples and of studies are read from
 * whichever field carries them, and every field is kept as given for the card: the track's extra columns are not
 * assumed to have fixed names.
 */
import type { SvReconstruction } from './svModel';

export type PopType = 'loss' | 'gain' | 'both' | 'unknown';
export interface PopulationSv {
  id: string; type: PopType;
  /** 0-based half-open: the widest and the narrowest extent of the merged variants */
  outer: [number, number]; inner: [number, number];
  /** the frequency as the source gives it, and as a fraction for sorting (null when absent) */
  freqText: string | null; freq: number | null;
  samples: number | null; studies: number | null;
  /** every other field, as given */
  fields: [string, string][];
}

const GEOMETRY = new Set(['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'itemRgb', 'blockCount', 'blockSizes', 'chromStarts']);
const num = (v: unknown): number | null => { const n = typeof v === 'number' ? v : Number(String(v ?? '').replace(/[,%\s]/g, '')); return Number.isFinite(n) && String(v ?? '').trim() !== '' ? n : null; };

/** The type a record states: a field naming gain / loss (or duplication / deletion), else the colour UCSC gives it. */
function typeOf(it: Record<string, unknown>): PopType {
  const text = Object.entries(it).filter(([k, v]) => !GEOMETRY.has(k) && typeof v === 'string' && /type|class|sub/i.test(k)).map(([, v]) => String(v)).join(' ').toLowerCase()
    || String(it.name ?? '').toLowerCase();
  const loss = /loss|deletion|\bdel\b/.test(text), gain = /gain|duplication|\bdup\b|insertion/.test(text);
  if (loss && gain) return 'both';
  if (loss) return 'loss';
  if (gain) return 'gain';
  // the colour of the item: red-ish losses, blue-ish gains (DGV's convention)
  const rgb = String(it.itemRgb ?? it.reserved ?? '').split(',').map(Number);
  if (rgb.length === 3 && rgb.every(Number.isFinite)) {
    const [r, , b] = rgb;
    if (r > b + 60) return 'loss';
    if (b > r + 60) return 'gain';
    if (r > 100 && b > 100) return 'both';
  }
  return 'unknown';
}

/** One record of the dgvGold track (a getData/track item); null when it has no usable coordinates. */
export function parseDgvItem(it: Record<string, unknown>): PopulationSv | null {
  const s = num(it.chromStart), e = num(it.chromEnd);
  if (s == null || e == null || e <= s) return null;
  let ts = num(it.thickStart), te = num(it.thickEnd);
  if (ts == null || te == null || te <= ts || ts < s || te > e) { ts = s; te = e; }
  const fields: [string, string][] = Object.entries(it).filter(([k]) => !GEOMETRY.has(k)).map(([k, v]) => [k, String(v)]);
  const pick = (re: RegExp) => fields.find(([k]) => re.test(k));
  const f = pick(/freq/i);
  let freq = f ? num(f[1]) : null;
  // a frequency written with % (or above 1) is a percentage
  if (freq != null && (/%/.test(f![1]) || freq > 1)) freq /= 100;
  const samples = pick(/unique.*sample|num.*sample|sample.*(count|num)/i), studies = pick(/num.*stud|stud.*(count|num)/i);
  return {
    id: String(it.name ?? `${s}-${e}`), type: typeOf(it), outer: [s, e], inner: [ts, te],
    freqText: f ? f[1] : null, freq, samples: samples ? num(samples[1]) : null, studies: studies ? num(studies[1]) : null, fields,
  };
}

/** A loss or a gain of the rebuilt SV: a run of adjacent segments with no copy (loss) or more than one (gain) on the derivatives. */
export interface CnvEvent { kind: 'loss' | 'gain'; start: number; end: number; letters: string; copies: number }

export function cnvEvents(rec: SvReconstruction, chrom: string): CnvEvent[] {
  const count = new Map<number, number>();
  for (const d of rec.ders) for (const p of d.pieces) count.set(p.seg.id, (count.get(p.seg.id) ?? 0) + 1);
  const segs = (rec.byChr[chrom] ?? rec.segs.filter(s => s.chrom === chrom)).slice().sort((a, b) => a.start - b.start);
  const out: CnvEvent[] = [];
  for (const s of segs) {
    const n = count.get(s.id) ?? 0;
    const kind = n === 0 ? 'loss' : n > 1 ? 'gain' : null;
    const last = out[out.length - 1];
    if (!kind) continue;
    if (last && last.kind === kind && last.end === s.start) { last.end = s.end; last.letters += s.letter; last.copies = Math.max(last.copies, n); }
    else out.push({ kind, start: s.start, end: s.end, letters: s.letter, copies: n });
  }
  return out;
}

const overlap = (a: [number, number], b: [number, number]) => Math.max(0, Math.min(a[1], b[1]) - Math.max(a[0], b[0]));
/** Reciprocal overlap: the overlap over the longer of the two (both must be covered by that fraction). */
export const reciprocal = (a: [number, number], b: [number, number]) => overlap(a, b) / Math.max(a[1] - a[0], b[1] - b[0], 1);

export interface PopMatch { sv: PopulationSv; ro: number; roInner: number; roOuter: number; covers: boolean; sameType: boolean }
/** The share of reciprocal overlap from which a population CNV counts as the same event. */
export const RO_MATCH = 0.5;

/** The population CNVs overlapping an event, best first: reciprocal overlap with the inner or the outer extent (the better), same type first. */
export function matchesOf(ev: CnvEvent, svs: PopulationSv[]): PopMatch[] {
  const e: [number, number] = [ev.start, ev.end];
  return svs.filter(v => overlap(e, v.outer) > 0).map(v => {
    const roInner = reciprocal(e, v.inner), roOuter = reciprocal(e, v.outer);
    const sameType = v.type === ev.kind || v.type === 'both' || v.type === 'unknown';
    // the narrowest extent of the population variant holds the whole event
    const covers = v.inner[0] <= ev.start && v.inner[1] >= ev.end;
    return { sv: v, ro: Math.max(roInner, roOuter), roInner, roOuter, covers, sameType };
  }).sort((a, b) => Number(b.sameType) - Number(a.sameType) || b.ro - a.ro || (b.sv.freq ?? 0) - (a.sv.freq ?? 0));
}
