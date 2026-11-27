/**
 * The array format of a patient: the design id as the files give it (CytoGenomics' probe file `AMADID <id>` line, the
 * interval report's `Array Design Info: <sample>[<design>_<date>]`), and, from the probe file, the nearest standard
 * format by its number of probes (8x60K, 4x180K, 2x400K, 1x1M) with the median spacing of the autosomes' probes. A design
 * id alone names no format here (no id → product mapping is embedded: none is certain for custom designs), so without a
 * probe file no format is deduced.
 */
export type ArrayChip = '60K' | '180K' | '400K' | '1M' | 'other';
/** the standard formats by their nominal number of probes per sample */
export const ARRAY_FORMATS: { chip: Exclude<ArrayChip, 'other'>; name: string; n: number }[] = [
  { chip: '60K', name: '8x60K', n: 60_000 }, { chip: '180K', name: '4x180K', n: 180_000 }, { chip: '400K', name: '2x400K', n: 400_000 }, { chip: '1M', name: '1x1M', n: 1_000_000 },
];
/** a file's probe count within this factor of a format's nominal count reads as that format; else `other` */
export const FORMAT_FACTOR = 1.4;
/** the format nearest a probe count (on a log scale), or `other` when none is within FORMAT_FACTOR */
export function nearestFormat(n: number): ArrayChip {
  if (!(n > 0)) return 'other';
  const best = ARRAY_FORMATS.reduce((a, f) => (Math.abs(Math.log(n / f.n)) < Math.abs(Math.log(n / a.n)) ? f : a));
  return Math.abs(Math.log(n / best.n)) <= Math.log(FORMAT_FACTOR) ? best.chip : 'other';
}
const AUTOSOME = /^chr([1-9]|1\d|2[0-2])$/;
/** the median distance between consecutive probe starts on chr1-22 (bp), or null under 2 probes */
export function medianSpacing(byChrom: Map<string, { start: ArrayLike<number> }>): number | null {
  const d: number[] = [];
  for (const [c, x] of byChrom) {
    if (!AUTOSOME.test(c)) continue;
    for (let i = 1; i < x.start.length; i++) d.push(x.start[i] - x.start[i - 1]);
  }
  if (!d.length) return null;
  d.sort((a, b) => a - b);
  const m = d.length >> 1;
  return d.length % 2 ? d[m] : (d[m - 1] + d[m]) / 2;
}
/** the design id in a file's header lines: `Array Design Info: sample[031746_20200303]` → `031746_20200303`; an AMADID / Design ID line */
export function designOfMeta(meta: [string, string][]): string | null {
  for (const [k, v] of meta) {
    if (/^array design info$/i.test(k.trim())) { const m = v.match(/\[([^\]]+)\]/); if (m) return m[1].trim(); }
    if (/^(amadid|design\s*id|array\s*design)$/i.test(k.trim()) && v.trim()) return v.trim().split(/\s+/)[0];
  }
  return null;
}
const SPACING = (bp: number) => (bp >= 1000 ? `~${Math.round(bp / 1000)} kb` : `~${Math.round(bp)} bp`);
export interface ArrayInfo {
  /** the tab's chip: a standard format, `other` (a probe file of no standard count), or null (no probe file: nothing shown) */
  chip: ArrayChip | null;
  /** the design ids as given, with where each was read */
  designs: { id: string; from: string }[];
  probes: number | null; spacing: number | null;
  /** one line: `≈ 60K array (61,246 probes, median spacing ~41 kb) · design 031746` */
  text: string;
  /** how it was deduced (the chip's tooltip) */
  how: string;
}
/** what is known of a patient's array: its design ids (report, probe file) and its probes (count incl. NaN, positions) */
export function arrayInfo(o: { designs?: { id: string | null | undefined; from: string }[]; probes?: { n: number; byChrom: Map<string, { start: ArrayLike<number> }> } | null }): ArrayInfo | null {
  const designs: { id: string; from: string }[] = [];
  // one design given twice (the report's `031746_20200303`, the probe file's `031746`): written once, where it was read first
  const same = (a: string, b: string) => a === b || a.startsWith(`${b}_`) || b.startsWith(`${a}_`);
  for (const d of o.designs ?? []) {
    if (!d.id) continue;
    const had = designs.find(x => same(x.id, d.id!));
    if (had) { if (!had.from.includes(d.from)) had.from += ` and ${d.from}`; } else designs.push({ id: d.id, from: d.from });
  }
  const des = designs.length ? ` · design ${designs.map(d => d.id).join(' / ')}` : '';
  const desHow = designs.map(d => `Design id ${d.id}: as given in ${d.from}.`).join('\n');
  if (o.probes && o.probes.n > 0) {
    const n = o.probes.n, sp = medianSpacing(o.probes.byChrom), chip = nearestFormat(n), f = ARRAY_FORMATS.find(x => x.chip === chip);
    const nums = `${n.toLocaleString('en-US')} probes${sp != null ? `, median spacing ${SPACING(sp)}` : ''}`;
    const text = `${chip === 'other' ? 'array of no standard format' : `≈ ${chip} array`} (${nums})${des}`;
    const how = [text, `Deduced from the probe file: ${n.toLocaleString('en-US')} probe rows (header left out, NaN counted)${sp != null ? `, median distance between consecutive probes on chr1-22 ${SPACING(sp)}` : ''}; ${
      chip === 'other' ? `no standard format (8x60K, 4x180K, 2x400K, 1x1M) within ×${FORMAT_FACTOR} of the count` : `nearest standard format by probe count: ${f!.name} (~${f!.n.toLocaleString('en-US')}, within ×${FORMAT_FACTOR})`}.`, desHow].filter(Boolean).join('\n');
    return { chip, designs, probes: n, spacing: sp, text, how };
  }
  if (!designs.length) return null;
  const text = `design ${designs.map(d => d.id).join(' / ')} (no probe file: format not deduced)`;
  return { chip: null, designs, probes: null, spacing: null, text, how: `${text}\n${desHow}\nA design id alone names no format here (no id → product list embedded).` };
}
/** the chip's colours: green 60K, orange 180K, blue 400K, purple 1M, grey other */
export const CHIP_STYLE: Record<ArrayChip, string> = {
  '60K': 'bg-emerald-100 text-emerald-800 border-emerald-300', '180K': 'bg-orange-100 text-orange-800 border-orange-300',
  '400K': 'bg-sky-100 text-sky-800 border-sky-300', '1M': 'bg-fuchsia-100 text-fuchsia-800 border-fuchsia-300', other: 'bg-slate-100 text-slate-600 border-slate-300',
};
