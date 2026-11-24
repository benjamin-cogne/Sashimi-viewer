/**
 * The probes' log ratios of an array (CytoGenomics' CGH_Log_Ratio.txt: `AMADID <design> hg19`, then Probe Name, ChrName,
 * Start, Stop, FeatureNum, Description, Name of Gene, Normalized log ratio(<sample>); NaN for a probe left out). Kept
 * as the file's bytes and read on first use into typed arrays per chromosome (about 12 bytes a probe), for a probe
 * the probe panel above the browser (CghProbePanel.tsx) and the DLRS.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { decodeText } from './cghFormats';
import { cghBuild, cghChrom } from './cghIntervals';

export interface CghProbeChrom { start: Uint32Array; end: Uint32Array; lr: Float32Array; name: string[] }
export interface CghProbes {
  sample: string; build: GenomeBuild | null;
  /** probes read (NaN included) */
  n: number;
  /** chr1…chrY, each sorted by start */
  byChrom: Map<string, CghProbeChrom>;
  /** derivative log ratio spread over all chromosomes (null: under 10 probes) */
  dlrs: number | null;
}

const LR_HEAD = /normalized\s+log\s*ratio\s*\(([^)]*)\)/i;
/** the sample a probe log-ratio file is of (its `Normalized log ratio(<sample>)` header), or null when the bytes are not one */
export function probeSampleOf(bytes: Uint8Array): string | null {
  const head = decodeText(bytes.subarray(0, 8192));
  for (const l of head.split(/\r\n|\n|\r/).slice(0, 30)) if (/^probe\s*name\t/i.test(l)) return l.match(LR_HEAD)?.[1].trim() ?? null;
  return null;
}

/** the file read: probes per chromosome, sorted by start */
export function parseProbes(bytes: Uint8Array): CghProbes {
  const lines = decodeText(bytes).split(/\r\n|\n|\r/);
  const hi = lines.findIndex(l => /^probe\s*name\t/i.test(l));
  if (hi < 0) throw new Error('not a probe log-ratio file');
  const h = lines[hi].split('\t').map(x => x.trim().toLowerCase());
  const cChr = h.findIndex(x => x === 'chrname' || x === 'chr' || x === 'chromosome'), cS = h.indexOf('start'), cE = h.indexOf('stop'), cL = h.findIndex(x => LR_HEAD.test(x));
  const sample = lines[hi].match(LR_HEAD)?.[1].trim() ?? '';
  const build = cghBuild(lines.slice(0, hi).join(' '));
  const acc = new Map<string, { s: number[]; e: number[]; v: number[]; p: string[] }>();
  const cName = Math.max(0, h.findIndex(x => x === 'probe name' || x === 'probename'));
  let n = 0;
  for (let i = hi + 1; i < lines.length; i++) {
    const f = lines[i].split('\t');
    const chrom = cghChrom(f[cChr] ?? ''), s = Number(f[cS]), e = Number(f[cE]);
    if (!chrom || chrom === 'chrM' || !Number.isFinite(s) || !Number.isFinite(e)) continue;
    let a = acc.get(chrom); if (!a) acc.set(chrom, a = { s: [], e: [], v: [], p: [] });
    const v = Number(f[cL]);
    a.s.push(s); a.e.push(e); a.v.push(Number.isFinite(v) ? v : NaN); a.p.push(f[cName] ?? ''); n++;
  }
  const byChrom = new Map<string, CghProbeChrom>();
  for (const [chrom, a] of acc) {
    const order = a.s.map((_, i) => i).sort((x, y) => a.s[x] - a.s[y]);
    byChrom.set(chrom, { start: Uint32Array.from(order, i => a.s[i]), end: Uint32Array.from(order, i => a.e[i]), lr: Float32Array.from(order, i => a.v[i]), name: order.map(i => a.p[i]) });
  }
  return { sample, build, n, byChrom, dlrs: dlrsOf(byChrom) };
}

/**
 * DLRS (derivative log ratio spread, Agilent's QC metric): the standard deviation of the differences between consecutive
 * probes (by position, within a chromosome, NaN skipped), divided by √2 (the difference of two noises doubles the variance).
 */
export function dlrsOf(byChrom: Map<string, CghProbeChrom>): number | null {
  let n = 0, sum = 0, sq = 0;
  for (const { lr } of byChrom.values()) {
    let prev = NaN;
    for (let i = 0; i < lr.length; i++) {
      const v = lr[i]; if (Number.isNaN(v)) continue;
      if (!Number.isNaN(prev)) { const d = v - prev; n++; sum += d; sq += d * d; }
      prev = v;
    }
  }
  return n < 10 ? null : Math.sqrt(Math.max(0, (sq - sum * sum / n) / (n - 1))) / Math.SQRT2;
}
/** Agilent guidance on the DLRS */
export const dlrsWord = (d: number) => (d < 0.2 ? 'excellent' : d <= 0.3 ? 'good' : 'poor');
export const DLRS_TIP = 'DLRS (derivative log ratio spread): the standard deviation of the log-ratio differences between consecutive probes, divided by √2; the noise of the array. Agilent guidance: below 0.2 excellent, 0.2 to 0.3 good, above 0.3 poor.';

/** the first probe whose start is ≥ pos (binary search) */
export function lowerBound(a: Uint32Array, pos: number): number {
  let lo = 0, hi = a.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < pos) lo = m + 1; else hi = m; }
  return lo;
}
/** the mean log ratio of the probes in [from, to] (NaN skipped), and their number */
export function meanIn(p: CghProbeChrom, from: number, to: number): { mean: number | null; n: number } {
  let n = 0, sum = 0;
  for (let i = lowerBound(p.start, from); i < p.start.length && p.start[i] <= to; i++) { const v = p.lr[i]; if (!Number.isNaN(v)) { n++; sum += v; } }
  return { mean: n ? sum / n : null, n };
}
