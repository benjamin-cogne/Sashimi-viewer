/**
 * Quality verdicts of a called variant site, from its evidence (VariantSite.q, counted by the full variant scan): four
 * checks in the spirit of GATK's hard filters and of what a reviewer looks at in IGV before trusting a call.
 *
 * - **BQ** (SNV): share of the alternate bases seen with a base quality under 20. They count in the allele fraction
 *   like the others; when they are most of what carries the allele, the call may be sequencing errors (a noisy cycle,
 *   a homopolymer tail) rather than an allele. Warn from 25 %, bad from 50 %.
 * - **HP** (indel, in place of BQ): length of the reference homopolymer at the site, where polymerase slippage (short
 *   reads) and basecalling (ONT, less PacBio HiFi) make indel errors. Short reads: warn from 6, bad from 10; long
 *   reads: warn from 4, bad from 7.
 * - **MQ**: share of the alternate reads with a mapping quality under 20, against the other reads over the site
 *   (GATK's MQ / MQRankSum idea). Bad when at least half of them map poorly while the other reads do not (+30 points),
 *   warn from 20 % (+15 points), or when the region maps poorly for every read (other reads ≥ 50 %).
 * - **SB** (strand bias): the alternate reads' strand split against the other reads' (Fisher's exact test, two-sided,
 *   as GATK's FisherStrand). Warn p < 0.001, bad p < 1e-4: an allele carried on one strand only is the signature of
 *   oxidative damage (8-oxoG, G>T), of a PCR or library artefact.
 * - **END** (read-position bias): the alternate calls within 10 bases of an alignment end against the other reads'
 *   (Fisher, one-sided; GATK's ReadPosRankSum idea): misaligned ends near an indel, adapter or clip artefacts.
 *   Warn p < 0.001, bad p < 1e-4.
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
/** Probability of one 2×2 table with the given margins (hypergeometric): a successes among r draws from a + c successes in n. */
const lnChoose = (n: number, k: number) => lnGamma(n + 1) - lnGamma(k + 1) - lnGamma(n - k + 1);
/**
 * Fisher's exact test of a 2×2 table [[a, b], [c, d]] (rows: alternate / reference reads; columns: yes / no), as
 * GATK's FisherStrand does for the strand split. `side`: 'two' for any difference, 'greater' for a larger share among
 * the alternate reads. Both rows are counts, so a reference side of few reads weighs as little as it should.
 */
export function fisher(a: number, b: number, c: number, d: number, side: 'two' | 'greater'): number {
  const r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
  if (!r1 || !r2 || !c1 || c1 === n) return 1;
  const lo = Math.max(0, c1 - r2), hi = Math.min(r1, c1), base = lnChoose(n, c1);
  const p = (x: number) => Math.exp(lnChoose(r1, x) + lnChoose(r2, c1 - x) - base);
  const pa = p(a);
  let s = 0;
  if (side === 'greater') for (let x = a; x <= hi; x++) s += p(x);
  else for (let x = lo; x <= hi; x++) { const px = p(x); if (px <= pa * (1 + 1e-7)) s += px; }
  return Math.min(1, s);
}
const pct = (x: number) => `${Math.round(x * 100)} %`;
const pText = (p: number) => (p < 1e-4 ? p.toExponential(0) : p < 0.01 ? p.toFixed(4) : p.toFixed(2));
/** Alternate reads needed before a strand or read-position test is judged. */
const MIN_TEST_READS = 5;
/**
 * Fisher p-values for amber and red. A window holds hundreds of sites, so a per-site cut-off at 0.01 would colour
 * a few of them by chance alone; GATK's hard filter for SNVs (FS > 60) sits near p = 1e-6.
 */
const P_WARN = 1e-3, P_BAD = 1e-4;

export function siteChecks(s: VariantSite, longReads: boolean): QCheck[] {
  const q = s.q;
  if (!q) return [];
  const n = s.alt_count, out: QCheck[] = [];
  if (s.kind === 'snv') {
    const lb = q.lowBq ?? 0;
    out.push({ key: 'BQ', level: lb >= 0.5 ? 'bad' : lb >= 0.25 ? 'warn' : 'good', value: pct(lb),
      detail: `${pct(lb)} of the alternate bases have a base quality under 20 (all are counted in the allele fraction)` });
  } else {
    const hp = q.hp ?? 0, [w, b] = longReads ? [4, 7] : [6, 10];
    out.push({ key: 'HP', level: hp >= b ? 'bad' : hp >= w ? 'warn' : 'good', value: `${hp}`,
      detail: hp > 1 ? `reference homopolymer of ${hp} at the site (indel errors grow with its length${longReads ? ', long reads especially' : ''})` : 'no homopolymer at the site' });
  }
  const lm = q.lowMq, lr = q.lowMqRef;
  out.push({ key: 'MQ', level: (lm >= 0.5 && lm - lr >= 0.3) ? 'bad' : ((lm >= 0.2 && lm - lr >= 0.15) || lr >= 0.5) ? 'warn' : 'good', value: pct(lm),
    detail: `${pct(lm)} of the alternate reads map with MAPQ < 20 (other reads: ${pct(lr)})${lr >= 0.5 ? '; the region maps poorly for every read' : ''}` });
  const nRef = q.nRef ?? 0;
  if (q.fwdRef == null || n < MIN_TEST_READS || nRef < MIN_TEST_READS) out.push({ key: 'SB', level: 'na', value: pct(q.fwd), detail: `${pct(q.fwd)} of the alternate reads on the + strand (too few reads to test)` });
  else {
    const af = Math.round(q.fwd * n), rf = Math.round(q.fwdRef * nRef), pv = fisher(af, n - af, rf, nRef - rf, 'two');
    out.push({ key: 'SB', level: pv < P_BAD ? 'bad' : pv < P_WARN ? 'warn' : 'good', value: pct(q.fwd),
      detail: `${pct(q.fwd)} of the ${n} alternate reads on the + strand, ${pct(q.fwdRef)} of the ${nRef} other reads (Fisher p = ${pText(pv)})` });
  }
  if (q.endRef == null || n < MIN_TEST_READS || nRef < MIN_TEST_READS) out.push({ key: 'END', level: 'na', value: pct(q.end), detail: `${pct(q.end)} of the alternate calls within 10 bases of a read end (too few reads to test)` });
  else {
    const ae = Math.round(q.end * n), re = Math.round(q.endRef * nRef), pv = fisher(ae, n - ae, re, nRef - re, 'greater');
    out.push({ key: 'END', level: pv < P_BAD ? 'bad' : pv < P_WARN ? 'warn' : 'good', value: pct(q.end),
      detail: `${pct(q.end)} of the ${n} alternate calls within 10 bases of a read end, ${pct(q.endRef)} of the ${nRef} other reads (Fisher p = ${pText(pv)}, one-sided)` });
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
