/**
 * Copy-number events detected from a sample's probes, a second opinion beside the file's calls. Segmentation: Circular
 * Binary Segmentation (Olshen AB et al., Biostatistics 2004;5:557-72; DNAcopy), with a fixed threshold on the statistic
 * instead of its permutation p-value (no permutations: under a second for 60k probes), the noise taken once per
 * sample. Segments are then called by their mean log ratio and number of probes; adjacent same-direction ones merged.
 * Pure: no React, no store (unit-tested on invented data, scripts/cgh-probe-fixture.mjs).
 */
import type { CghProbeChrom, CghProbes } from './cghProbes';
import { PAR } from './cghMosaic';

export interface SegOpts {
  /** fewest probes of an event */
  minProbes: number;
  /** smallest |mean log ratio| of an event */
  minLr: number;
  /** the split threshold: the segment-vs-rest difference in noise units (higher: fewer, surer segments) */
  t: number;
}
/** sensitivity → split threshold. 8 calibrated locally on real arrays against ADM-2 calls (ARRAY-CGH.md) */
export const SENSITIVITY = { high: 6, normal: 8, low: 10 } as const;
export type Sensitivity = keyof typeof SENSITIVITY;
export const SEG_DEFAULTS: SegOpts = { minProbes: 3, minLr: 0.25, t: SENSITIVITY.normal };

/** a segment: probes i0..i1 of the chromosome (indices into its arrays, NaN ones inside not counted) */
export interface ProbeSeg { chrom: string; i0: number; i1: number; start: number; end: number; n: number; mean: number }
export interface ProbeEvent extends ProbeSeg { kind: 'gain' | 'loss' }
export interface ProbeCalls { sigma: number; segs: Map<string, ProbeSeg[]>; events: ProbeEvent[]; ms: number }

/** the noise per probe, robust to the events and to slow waves: MAD of consecutive differences / √2 */
export function noiseOf(byChrom: Map<string, CghProbeChrom>): number {
  const d: number[] = [];
  for (const { lr } of byChrom.values()) {
    let prev = NaN;
    for (let i = 0; i < lr.length; i++) { const v = lr[i]; if (Number.isNaN(v)) continue; if (!Number.isNaN(prev)) d.push(Math.abs(v - prev)); prev = v; }
  }
  if (!d.length) return NaN;
  d.sort((a, b) => a - b);
  return (d[d.length >> 1] * 1.4826) / Math.SQRT2;
}

/**
 * CBS of x[l, r): the arc (i, j] whose mean differs most from the rest, as Z = (S_j - S_i - k·T/n)·√(n / (k(n-k))) / σ
 * (k = j - i); split there when |Z| ≥ t, and recurse on the up-to-three pieces. Returns the breakpoints.
 */
function cbs(x: Float64Array, l: number, r: number, sigma: number, t: number, out: number[]) {
  const n = r - l;
  if (n >= 2) {
    const S = new Float64Array(n + 1);
    for (let i = 0; i < n; i++) S[i + 1] = S[i] + x[l + i];
    const T = S[n], w = new Float64Array(n);
    for (let k = 1; k < n; k++) w[k] = Math.sqrt(n / (k * (n - k))) / sigma;
    const m = T / n;
    let best = 0, bi = 0, bj = 0;
    for (let i = 0; i < n; i++) {
      const si = S[i];
      for (let j = i + 1; j <= n; j++) {
        const k = j - i; if (k === n) break;
        const z = Math.abs(S[j] - si - k * m) * w[k];
        if (z > best) { best = z; bi = i; bj = j; }
      }
    }
    if (best >= t) {
      const cuts = [0, bi, bj, n].filter((v, q, a) => q === 0 || v !== a[q - 1]);
      for (let q = 0; q + 1 < cuts.length; q++) {
        cbs(x, l + cuts[q], l + cuts[q + 1], sigma, t, out);
        if (q + 2 < cuts.length) out.push(l + cuts[q + 1]);
      }
      return;
    }
  }
}

/** one chromosome's segments (NaN probes skipped) */
export function segmentChrom(chrom: string, p: CghProbeChrom, sigma: number, t: number): ProbeSeg[] {
  const idx: number[] = [];
  for (let i = 0; i < p.lr.length; i++) if (!Number.isNaN(p.lr[i])) idx.push(i);
  if (!idx.length) return [];
  const x = Float64Array.from(idx, i => p.lr[i]), cuts: number[] = [];
  cbs(x, 0, x.length, sigma, t, cuts);
  const bounds = [0, ...cuts.sort((a, b) => a - b), x.length], segs: ProbeSeg[] = [];
  for (let q = 0; q + 1 < bounds.length; q++) {
    const a = bounds[q], b = bounds[q + 1];
    let s = 0; for (let i = a; i < b; i++) s += x[i];
    const i0 = idx[a], i1 = idx[b - 1];
    segs.push({ chrom, i0, i1, start: p.start[i0], end: p.end[i1], n: b - a, mean: s / (b - a) });
  }
  return segs;
}

/** the events of a chromosome's segments: adjacent ones of the same direction merged, then kept by size and mean */
export function callSegs(segs: ProbeSeg[], o: SegOpts): ProbeEvent[] {
  const dirOf = (s: ProbeSeg) => (s.mean >= o.minLr ? 1 : s.mean <= -o.minLr ? -1 : 0);
  const ds = segs.map(dirOf);
  // a segment of under minProbes probes between two of one direction is noise inside one event (one stray probe split a gain in two)
  for (let k = 1; k + 1 < segs.length; k++) if (segs[k].n < o.minProbes && ds[k - 1] && ds[k - 1] === ds[k + 1]) ds[k] = ds[k - 1];
  const merged: ProbeSeg[] = [], md: number[] = [];
  segs.forEach((s, k) => {
    const last = merged[merged.length - 1], d = ds[k];
    if (d && last && md[md.length - 1] === d) {
      const n = last.n + s.n;
      merged[merged.length - 1] = { ...last, i1: s.i1, end: s.end, n, mean: (last.mean * last.n + s.mean * s.n) / n };
    } else { merged.push({ ...s }); md.push(d); }
  });
  const ev: ProbeEvent[] = [];
  for (const [k, s] of merged.entries()) {
    const d = md[k];
    if (d && Math.abs(s.mean) >= o.minLr && s.n >= o.minProbes) ev.push({ ...s, kind: d > 0 ? 'gain' : 'loss' });
  }
  return ev;
}

const inPar = (chrom: string, pos: number) => (PAR[chrom] ?? []).some(([a, b]) => pos >= a && pos <= b);
/**
 * chrY: its PAR probes are the chrX PAR ones again (same values), so left out (else each PAR event twice). Its other
 * probes are cross-hybridisation noise in a female (many false events), so chrY is read
 * only for a patient said male.
 */
const chrYPart = (p: CghProbeChrom, male: boolean): CghProbeChrom =>
  ({ ...p, lr: male ? p.lr.map((v, i) => (inPar('chrY', p.start[i]) ? NaN : v)) : p.lr.map(() => NaN) });

/** a sample's segments (kept for drawing) and events */
export function detect(pr: CghProbes, o: SegOpts = SEG_DEFAULTS, male = false): ProbeCalls {
  const t0 = performance.now();
  const sigma = noiseOf(pr.byChrom), segs = new Map<string, ProbeSeg[]>(), events: ProbeEvent[] = [];
  if (Number.isFinite(sigma) && sigma > 0) for (const [chrom, p] of pr.byChrom) {
    const s = segmentChrom(chrom, chrom === 'chrY' ? chrYPart(p, male) : p, sigma, o.t);
    segs.set(chrom, s); events.push(...callSegs(s, o));
  }
  return { sigma, segs, events, ms: performance.now() - t0 };
}

/** reciprocal overlap of two intervals (the overlap over the longer one) */
export const recip = (a: { start: number; end: number }, b: { start: number; end: number }) => {
  const ov = Math.min(a.end, b.end) - Math.max(a.start, b.start) + 1;
  return ov <= 0 ? 0 : ov / Math.max(a.end - a.start + 1, b.end - b.start + 1);
};
/** two calls agree: same chromosome, same direction, reciprocal overlap ≥ 50 % */
export const RECIP_MIN = 0.5;
type Iv = { chrom: string; kind: string; start: number; end: number };
/** a chrY PAR interval as its chrX twin (GRCh37 PARs have equal lengths): a file's chrY PAR call meets the chrX event */
const asX = (c: Iv): Iv => {
  const k = PAR.chrY.findIndex(([a, b]) => c.start >= a && c.end <= b);
  if (c.chrom !== 'chrY' || k < 0) return c;
  const d = PAR.chrX[k][0] - PAR.chrY[k][0];
  return { ...c, chrom: 'chrX', start: c.start + d, end: c.end + d };
};
export const agrees = (a: Iv, b: Iv) => { const x = asX(a), y = asX(b); return x.chrom === y.chrom && x.kind === y.kind && recip(x, y) >= RECIP_MIN; };
