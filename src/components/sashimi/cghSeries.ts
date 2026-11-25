/**
 * Recurrence of the calls in the loaded series: for each call, the OTHER patients having a call of the same kind on
 * the same chromosome with a reciprocal overlap of at least SERIES_OVERLAP. Every call of every patient counts, those
 * under the probe or size minimum too (a recurrent call is worth knowing even when not analysed in the other patient).
 */
import type { CghFile, CghCnv, CghSex } from './cghIntervals';

/** a call hidden because it is the normal one for the patient's sex (a chrX row typed by copy number in a file without the sex: x1 male, x2 female) */
export const hiddenBySex = (c: CghCnv, sex: CghSex | undefined) => !!c.cnSex && !!sex && c.copies === (sex === 'male' ? 1 : 2);

/** reciprocal overlap: the shared part is at least this share of each of the two calls */
export const SERIES_OVERLAP = 0.7;
export interface SeriesHit { name: string; key: string; pct: number }
export interface SeriesInfo { hits: SeriesHit[]; /** the other patients loaded */ others: number }

/** a sample name compared: trimmed, case-insensitive */
export const nameKey = (n: string) => n.trim().toLowerCase();
/** the sample names loaded more than once, with the files holding them: "26X0001 loaded twice (file A, file B)" */
export function duplicateNames(files: CghFile[]): string[] {
  const by = new Map<string, { name: string; files: string[] }>();
  for (const f of files) for (const x of f.samples) { const k = nameKey(x.name), e = by.get(k) ?? { name: x.name.trim(), files: [] }; e.files.push(f.name); by.set(k, e); }
  return [...by.values()].filter(e => e.files.length > 1).map(e => `${e.name} loaded ${e.files.length === 2 ? 'twice' : `${e.files.length} times`} (${e.files.join(', ')}): not counted as another patient in Series`);
}
let memoFiles: CghFile[] | null = null, memoSex: Record<string, CghSex> | null = null, memo = new Map<string, SeriesInfo>();
/** by call id, computed once per list of files and sexes; a call hidden for its patient's sex neither counts nor gets a value */
export function seriesOf(files: CghFile[], sex: Record<string, CghSex> = {}): Map<string, SeriesInfo> {
  if (files === memoFiles && sex === memoSex) return memo;
  const samples = files.flatMap(f => f.samples);
  // a patient loaded twice (same name) is one patient: never its own recurrence
  const names = new Set(samples.map(x => nameKey(x.name)));
  const groups = new Map<string, { start: number; end: number; key: string; name: string }[]>();
  let maxLen = 0;
  for (const x of samples) for (const c of x.cnvs) {
    if (hiddenBySex(c, sex[x.key])) continue;
    const g = `${c.chrom}|${c.kind}`; let l = groups.get(g); if (!l) groups.set(g, l = []);
    l.push({ start: c.start, end: c.end, key: x.key, name: x.name }); maxLen = Math.max(maxLen, c.end - c.start + 1);
  }
  for (const l of groups.values()) l.sort((a, b) => a.start - b.start);
  const out = new Map<string, SeriesInfo>();
  for (const x of samples) for (const c of x.cnvs) {
    if (hiddenBySex(c, sex[x.key])) continue;
    const l = groups.get(`${c.chrom}|${c.kind}`) ?? [], len = c.end - c.start + 1, best = new Map<string, SeriesHit>(); // by sample name: one patient loaded twice counts once
    // first call starting late enough to reach the call (start >= c.start - maxLen), by bisection
    let lo = 0, hi = l.length; while (lo < hi) { const m = (lo + hi) >> 1; if (l[m].start < c.start - maxLen) lo = m + 1; else hi = m; }
    for (let i = lo; i < l.length && l[i].start <= c.end; i++) {
      const o = l[i]; if (o.key === x.key || nameKey(o.name) === nameKey(x.name)) continue;
      const ov = Math.min(c.end, o.end) - Math.max(c.start, o.start) + 1;
      if (ov <= 0) continue;
      const pct = ov / Math.max(len, o.end - o.start + 1);
      if (pct < SERIES_OVERLAP) continue;
      const nk = nameKey(o.name), p = best.get(nk); if (!p || pct > p.pct) best.set(nk, { name: o.name, key: o.key, pct });
    }
    out.set(c.id, { hits: [...best.values()].sort((a, b) => b.pct - a.pct), others: names.size - (names.has(nameKey(x.name)) ? 1 : 0) });
  }
  memoFiles = files; memoSex = sex; memo = out;
  return out;
}
