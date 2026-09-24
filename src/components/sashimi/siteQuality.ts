/**
 * Quality verdicts of a called variant site, from its evidence (VariantSite.q, counted by the full variant scan): four
 * checks in the spirit of GATK's hard filters and of what a reviewer looks at in IGV before trusting a call.
 *
 * - **BQ** (SNV): share of the alternate bases seen with a base quality under 20. They are left out of the allele
 *   fraction; when they are most of what carries the allele, the call leans on its few good bases (sequencing errors
 *   of a noisy cycle or homopolymer tail). Warn from 25 %, bad from 50 %.
 * - **HP** (indel, in place of BQ): length of the reference homopolymer at the site, where polymerase slippage (short
 *   reads) and basecalling (ONT, less PacBio HiFi) make indel errors. Short reads: warn from 6, bad from 10; long
 *   reads: warn from 4, bad from 7.
 * - **MQ**: share of the alternate reads with a mapping quality under 20, against the other reads over the site
 *   (GATK's MQ / MQRankSum idea). Bad when at least half of them map poorly while the other reads do not (+30 points),
 *   warn from 20 % (+15 points), or when the region maps poorly for every read (other reads ≥ 50 %).
 * - **SB** (strand bias): the alternate reads' + strand share tested against the other reads' (binomial, two-sided;
 *   GATK's FisherStrand compares the same two strand splits). Warn p < 0.01, bad p < 0.001: an allele carried on one
 *   strand only is the signature of oxidative damage (8-oxoG, G>T), of a PCR or library artefact.
 * - **END** (read-position bias): the share of alternate calls within 10 bases of an alignment end, tested one-sided
 *   against the other reads' share (GATK's ReadPosRankSum idea): misaligned ends near an indel, adapter or clip
 *   artefacts. Warn p < 0.01, bad p < 0.001.
 *
 * Each test needs a few alternate reads to say anything; with fewer, the check is not judged.
 */
import type { AlignedRead, VariantSite } from './types';

export type QLevel = 'good' | 'warn' | 'bad' | 'na';
export interface QCheck { key: 'BQ' | 'HP' | 'MQ' | 'SB' | 'END'; level: QLevel; value: string; detail: string }

export const Q_COLORS: Record<QLevel, string> = { good: '#16a34a', warn: '#f59e0b', bad: '#dc2626', na: '#d1d5db' };
const RANK: Record<QLevel, number> = { na: 0, good: 1, warn: 2, bad: 3 };
export const worstLevel = (checks: QCheck[]): QLevel => checks.reduce<QLevel>((w, c) => (RANK[c.level] > RANK[w] ? c.level : w), 'na');

/** ln Γ(x), Lanczos approximation (g = 7, 9 terms; relative error under 1e-13). */
function lnGamma(x: number): number {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lnGamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + 7.5;
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}
/** P(X ≥ k) for X ~ Binomial(n, p): the first term from log-gamma, the next ones by recurrence, O(n − k). */
export function binomUpper(n: number, k: number, p: number): number {
  if (k <= 0) return 1;
  if (k > n) return 0;
  let term = Math.exp(lnGamma(n + 1) - lnGamma(k + 1) - lnGamma(n - k + 1) + k * Math.log(p) + (n - k) * Math.log(1 - p)), s = 0;
  const r = p / (1 - p);
  for (let i = k; i <= n; i++) { s += term; term *= ((n - i) / (i + 1)) * r; if (term < s * 1e-16) break; }
  return Math.min(1, s);
}
/** Two-sided binomial p-value (twice the smaller tail). */
function binomTwoSided(n: number, k: number, p: number): number {
  const upper = binomUpper(n, k, p), lower = 1 - binomUpper(n, k + 1, p);
  return Math.min(1, 2 * Math.min(upper, lower));
}
const pct = (x: number) => `${Math.round(x * 100)} %`;
const pText = (p: number) => (p < 1e-4 ? p.toExponential(0) : p < 0.01 ? p.toFixed(4) : p.toFixed(2));
/** Alternate reads needed before a strand or read-position test is judged. */
const MIN_TEST_READS = 5;

export function siteChecks(s: VariantSite, longReads: boolean): QCheck[] {
  const q = s.q;
  if (!q) return [];
  const n = s.alt_count, out: QCheck[] = [];
  if (s.kind === 'snv') {
    const lb = q.lowBq ?? 0;
    out.push({ key: 'BQ', level: lb >= 0.5 ? 'bad' : lb >= 0.25 ? 'warn' : 'good', value: pct(lb),
      detail: `${pct(lb)} of the alternate bases have a base quality under 20 (left out of the allele fraction)` });
  } else {
    const hp = q.hp ?? 0, [w, b] = longReads ? [4, 7] : [6, 10];
    out.push({ key: 'HP', level: hp >= b ? 'bad' : hp >= w ? 'warn' : 'good', value: `${hp}`,
      detail: hp > 1 ? `reference homopolymer of ${hp} at the site (indel errors grow with its length${longReads ? ', long reads especially' : ''})` : 'no homopolymer at the site' });
  }
  const lm = q.lowMq, lr = q.lowMqRef;
  out.push({ key: 'MQ', level: (lm >= 0.5 && lm - lr >= 0.3) ? 'bad' : ((lm >= 0.2 && lm - lr >= 0.15) || lr >= 0.5) ? 'warn' : 'good', value: pct(lm),
    detail: `${pct(lm)} of the alternate reads map with MAPQ < 20 (other reads: ${pct(lr)})${lr >= 0.5 ? '; the region maps poorly for every read' : ''}` });
  if (q.fwdRef == null || n < MIN_TEST_READS) out.push({ key: 'SB', level: 'na', value: pct(q.fwd), detail: `${pct(q.fwd)} of the alternate reads on the + strand (too few reads to test)` });
  else {
    const p0 = Math.min(0.95, Math.max(0.05, q.fwdRef)), pv = binomTwoSided(n, Math.round(q.fwd * n), p0);
    out.push({ key: 'SB', level: pv < 1e-3 ? 'bad' : pv < 1e-2 ? 'warn' : 'good', value: pct(q.fwd),
      detail: `${pct(q.fwd)} of the alternate reads on the + strand, ${pct(q.fwdRef)} of the other reads (binomial p = ${pText(pv)})` });
  }
  if (q.endRef == null || n < MIN_TEST_READS) out.push({ key: 'END', level: 'na', value: pct(q.end), detail: `${pct(q.end)} of the alternate calls within 10 bases of a read end (too few reads to test)` });
  else {
    const p0 = Math.min(0.95, Math.max(0.01, q.endRef)), pv = binomUpper(n, Math.round(q.end * n), p0);
    out.push({ key: 'END', level: pv < 1e-3 ? 'bad' : pv < 1e-2 ? 'warn' : 'good', value: pct(q.end),
      detail: `${pct(q.end)} of the alternate calls within 10 bases of a read end, ${pct(q.endRef)} of the other reads' bases (binomial p = ${pText(pv)})` });
  }
  return out;
}

/** One read over a site, for the detail panel: its mapping quality, strand, base quality at an SNV, distance to its nearer alignment end. */
export interface SiteObs { mq: number; fwd: boolean; bq: number | null; end: number }
export interface ReadEvidence { alt: SiteObs[]; ref: SiteObs[]; other: number }
/**
 * The reads over a site split by allele, exact (the reads decoded around it): what the scan's shares summarise, as
 * distributions. The reference side has no base quality here (only the mismatching bases' are kept with a read).
 */
export function readEvidence(reads: AlignedRead[], s: VariantSite): ReadEvidence {
  const out: ReadEvidence = { alt: [], ref: [], other: 0 };
  const endLen = s.kind === 'del' ? s.length : 1;
  for (const r of reads) {
    const inBlock = (p: number) => r.b.some(([a, b]) => a <= p && p < b);
    const obs = (bq: number | null): SiteObs => ({ mq: r.q, fwd: r.r === 0, bq, end: Math.max(0, Math.min(s.pos - r.s, r.e - (s.pos + endLen))) });
    if (s.kind === 'snv') {
      if (!inBlock(s.pos)) continue;
      const m = r.m.find(x => x[0] === s.pos);
      if (!m) out.ref.push(obs(null));
      else if (m[1] === s.alt) out.alt.push(obs(m[2]));
      else out.other++;
    } else if (s.kind === 'ins') {
      const alt = r.i.some(([p, l]) => p === s.pos && l === s.length);
      if (alt) out.alt.push(obs(null));
      else if (r.i.some(([p]) => p === s.pos)) out.other++;
      else if (r.b.some(([a, b]) => a < s.pos && s.pos < b)) out.ref.push(obs(null));
    } else {
      const alt = r.d.some(([a, b]) => a === s.pos && b === s.pos + s.length);
      if (alt) out.alt.push(obs(null));
      else if (r.d.some(([a, b]) => a <= s.pos && s.pos < b)) out.other++;
      else if (inBlock(s.pos)) out.ref.push(obs(null));
    }
  }
  return out;
}
