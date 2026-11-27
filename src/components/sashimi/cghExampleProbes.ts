/**
 * Invented probe-level log ratios for the built-in examples (cghExamples.ts), written as a CytoGenomics CGH_Log_Ratio.txt
 * at load time: no embedded data, a seeded random generator (same file at every load). About 60k probes over the GRCh37
 * chromosomes (like a 60K array: none in the centromeres, the acrocentric p arms nor chrY's heterochromatin), Gaussian
 * noise (SD 0.15). Hybridised against a sex-matched reference: a male's X and Y lie around 0 like the autosomes; a female has
 * no chrY probe written (no Y DNA: CytoGenomics' ratios there are noise, not data). Each event gets exactly its call's number of probes, evenly spaced from the call's start to its stop,
 * their mean set to the call's mean log ratio, so the file's calls and the probes agree. Invented data only.
 */
import { chromLength } from '../../standalone/assembly';

/** an event in the probes: [chrom, start, stop, probes, mean log ratio] */
export type ProbeEvent = [chrom: string, start: number, stop: number, n: number, lr: number];

const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
/** no probe there (Mb, GRCh37, rounded): centromeres and heterochromatin, the acrocentric p arms, chrY outside its euchromatin */
const GAPS: Record<string, [number, number][]> = {
  1: [[121, 144]], 2: [[91.8, 95.3]], 3: [[89.5, 93.6]], 4: [[48.9, 52]], 5: [[46.9, 49.9]], 6: [[59.5, 62.5]], 7: [[58.4, 61.4]],
  8: [[44.1, 47.1]], 9: [[38.7, 71]], 10: [[38.7, 42.4]], 11: [[51.6, 55]], 12: [[34.3, 38.3]], 13: [[0, 19.2]], 14: [[0, 19.2]],
  15: [[0, 20.1]], 16: [[35.3, 46.5]], 17: [[22.3, 25.3]], 18: [[15.4, 18.6]], 19: [[24.6, 28.6]], 20: [[26.3, 29.8]],
  21: [[0, 14.4]], 22: [[0, 16.2]], X: [[58.6, 61.7]], Y: [[0, 2.7], [28.8, 60]],
};
/** probes of the backbone, events aside */
const BACKBONE = 59_600, PROBE = 60, EDGE = 200_000, SD = 0.15;

/** mulberry32, seeded from the sample's name */
function rng(name: string) {
  let seed = 0; for (const ch of name) seed = (Math.imul(seed, 31) + ch.charCodeAt(0)) | 0;
  const rnd = () => { seed = (seed + 0x6d2b79f5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  return { rnd, gauss: () => Math.sqrt(-2 * Math.log(rnd() || 1e-12)) * Math.cos(2 * Math.PI * rnd()) };
}
const inGap = (c: string, p: number) => (GAPS[c] ?? []).some(([a, b]) => p >= a * 1e6 && p <= b * 1e6);

/** the CGH_Log_Ratio.txt of an invented sample: its events, its sex (a sex-matched reference; no chrY probe for a female) */
export function exampleProbeFile(sample: string, events: ProbeEvent[], male: boolean): string {
  const { rnd, gauss } = rng(sample);
  const len = (c: string) => chromLength('GRCh37', c) ?? 0;
  // the span the backbone covers, for its step
  let covered = 0;
  for (const c of CHROMS) { const L = len(c) - 2 * EDGE; covered += L - (GAPS[c] ?? []).reduce((k, [a, b]) => k + Math.max(0, Math.min(b * 1e6, len(c)) - a * 1e6), 0); }
  const step = covered / BACKBONE;
  const out = ['AMADID\t000000\thg19', `Probe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(${sample})`];
  let k = 0;
  for (const c of CHROMS) {
    if (c === 'Y' && !male) continue;
    const evs = events.filter(e => e[0] === 'chr' + c || e[0] === c);
    // an event within 1 Mb of an end is terminal: no backbone probe beyond it
    const lo = Math.min(Infinity, ...evs.filter(e => e[1] < 1e6).map(e => e[1])), hi = Math.max(-Infinity, ...evs.filter(e => e[2] > len(c) - 1e6).map(e => e[2]));
    const ps: { s: number; v: number }[] = [];
    for (let p = EDGE + rnd() * step; p < len(c) - EDGE; p += step * (0.6 + 0.8 * rnd())) {
      const s = Math.round(p);
      if (inGap(c, s) || (lo < Infinity && s < lo) || (hi > -Infinity && s > hi) || evs.some(e => s + PROBE - 1 >= e[1] && s <= e[2])) continue;
      ps.push({ s, v: rnd() < 0.003 ? NaN : gauss() * SD });
    }
    // each event: exactly its probes, first at its start, last ending at its stop, their mean the call's
    for (const [, a, b, n, lr] of evs) {
      const vs = Array.from({ length: n }, () => gauss() * SD), m = vs.reduce((x, y) => x + y, 0) / n;
      for (let i = 0; i < n; i++) ps.push({ s: n === 1 ? a : Math.round(a + (i * (b - PROBE + 1 - a)) / (n - 1)), v: vs[i] - m + lr });
    }
    ps.sort((x, y) => x.s - y.s);
    for (const { s, v } of ps) { k++; out.push(`A_00_P${String(k).padStart(6, '0')}\t${c}\t${s}\t${s + PROBE - 1}\t${k}\t\t\t${Number.isNaN(v) ? 'NaN' : v.toFixed(4)}`); }
  }
  return out.join('\n') + '\n';
}
