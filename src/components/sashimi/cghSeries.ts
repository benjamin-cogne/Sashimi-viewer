/**
 * Recurrence of the calls in the loaded series: for each call, the OTHER patients having a call of the same kind on
 * the same chromosome with a reciprocal overlap of at least SERIES_OVERLAP. Every call of every patient counts, those
 * under the probe or size minimum too (a recurrent call is worth knowing even when not analysed in the other patient).
 */
import type { CghFile } from './cghIntervals';

/** reciprocal overlap: the shared part is at least this share of each of the two calls */
export const SERIES_OVERLAP = 0.7;
export interface SeriesHit { name: string; key: string; pct: number }
export interface SeriesInfo { hits: SeriesHit[]; /** the other patients loaded */ others: number }

let memoFiles: CghFile[] | null = null, memo = new Map<string, SeriesInfo>();
/** by call id, computed once per list of files */
export function seriesOf(files: CghFile[]): Map<string, SeriesInfo> {
  if (files === memoFiles) return memo;
  const samples = files.flatMap(f => f.samples);
  const total = samples.length;
  const groups = new Map<string, { start: number; end: number; key: string; name: string }[]>();
  let maxLen = 0;
  for (const x of samples) for (const c of x.cnvs) {
    const g = `${c.chrom}|${c.kind}`; let l = groups.get(g); if (!l) groups.set(g, l = []);
    l.push({ start: c.start, end: c.end, key: x.key, name: x.name }); maxLen = Math.max(maxLen, c.end - c.start + 1);
  }
  for (const l of groups.values()) l.sort((a, b) => a.start - b.start);
  const out = new Map<string, SeriesInfo>();
  for (const x of samples) for (const c of x.cnvs) {
    const l = groups.get(`${c.chrom}|${c.kind}`) ?? [], len = c.end - c.start + 1, best = new Map<string, SeriesHit>();
    // first call starting late enough to reach the call (start >= c.start - maxLen), by bisection
    let lo = 0, hi = l.length; while (lo < hi) { const m = (lo + hi) >> 1; if (l[m].start < c.start - maxLen) lo = m + 1; else hi = m; }
    for (let i = lo; i < l.length && l[i].start <= c.end; i++) {
      const o = l[i]; if (o.key === x.key) continue;
      const ov = Math.min(c.end, o.end) - Math.max(c.start, o.start) + 1;
      if (ov <= 0) continue;
      const pct = ov / Math.max(len, o.end - o.start + 1);
      if (pct < SERIES_OVERLAP) continue;
      const p = best.get(o.key); if (!p || pct > p.pct) best.set(o.key, { name: o.name, key: o.key, pct });
    }
    out.set(c.id, { hits: [...best.values()].sort((a, b) => b.pct - a.pct), others: total - 1 });
  }
  memoFiles = files; memo = out;
  return out;
}
