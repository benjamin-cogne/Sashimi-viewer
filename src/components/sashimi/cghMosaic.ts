/**
 * Fraction of cells carrying a one-copy change, from the mean log2 ratio of the call: a port of the AChro-Puce
 * mosaicism calculator (https://acpa-achropuce.com/calcul-du-mosaicisme/, tool of July 2025), formulas unchanged.
 * Indicative only: the calculator ignores the biological variability of the ratio.
 */
import type { CghCnv, CghSex } from './cghIntervals';
import type { GenomeBuild } from '../../standalone/ensembl';

/** below this fraction a call is tagged mosaic in the table */
export const MOSAIC_TAG_BELOW = 0.8;

/** the calculator's own function: chrom 'autosome' | 'X' | 'Y', sex 'F' (anything else: male), event 'del' | 'dup' */
function computeMosaicism(lr: number, event: 'del' | 'dup', chrom: 'autosome' | 'X' | 'Y', sexe: 'F' | 'M'): number {
  let f = 0;
  const hemi = (chrom === 'X' || chrom === 'Y') && sexe !== 'F';
  if (event === 'del') f = hemi ? 1 - Math.pow(2, lr) : 2 - 2 * Math.pow(2, lr);
  else f = hemi ? Math.pow(2, lr) - 1 : Math.pow(2, lr) * 2 - 2;
  return Math.min(1, Math.max(0, f));
}

/**
 * Pseudo-autosomal regions per build (1-based inclusive; GRC's PAR1 / PAR2). Checked on the UCSC API (2026-11-25): the
 * chrX and chrY sequences are identical inside each region and differ just outside its inner ends.
 */
export const PARS: Record<GenomeBuild, Record<string, [number, number][]>> = {
  GRCh37: { chrX: [[60001, 2699520], [154931044, 155260560]], chrY: [[10001, 2649520], [59034050, 59363566]] },
  GRCh38: { chrX: [[10001, 2781479], [155701383, 156030895]], chrY: [[10001, 2781479], [56887903, 57217415]] },
};
/** the PARs of a build (GRCh37 when none is stated) */
export const parOf = (build: GenomeBuild | null | undefined) => PARS[build ?? 'GRCh37'] ?? PARS.GRCh37;
/** GRCh37's, kept for callers without a build */
export const PAR = PARS.GRCh37;
const inPar = (c: CghCnv, build?: GenomeBuild | null) => (parOf(build)[c.chrom] ?? []).some(([a, b]) => c.start >= a && c.end <= b);

export type Mosaic =
  | { ok: true; f: number; expected: number; hemi: boolean; why?: undefined }
  | { ok: false; why: 'sex' };

/** the estimate for a call, or null when none applies (no log ratio, not a one-copy change); `sex needed` for chrX / chrY outside the PARs */
export function mosaicOf(c: CghCnv, sex: CghSex | undefined, build?: GenomeBuild | null): Mosaic | null {
  if (c.logRatio == null || c.amp || c.copiesHi != null) return null;
  const sexChrom = (c.chrom === 'chrX' || c.chrom === 'chrY') && !inPar(c, build);
  if (sexChrom && !sex) return { ok: false, why: 'sex' };
  const hemi = sexChrom && sex === 'male';
  // a male's single X / Y: the one-copy change is x0 (loss) or x2 (gain); elsewhere x1 or x3
  const base = hemi ? 1 : 2;
  if (c.copies != null && c.copies !== (c.kind === 'loss' ? base - 1 : base + 1)) return null;
  const chrom = sexChrom ? (c.chrom === 'chrX' ? 'X' : 'Y') : 'autosome';
  const ev = c.kind === 'loss' ? 'del' : 'dup', sx = sex === 'male' ? 'M' : 'F';
  const f = computeMosaicism(c.logRatio, ev, chrom, sx);
  return { ok: true, f, expected: Math.log2((c.kind === 'loss' ? base - 1 : base + 1) / base), hemi };
}

export interface EstCopies { n: number; expected: number; hemi: boolean; mosaic: number | null; text: string }
/**
 * A copy number estimated from the log ratio, for a call whose file gives none (an interval-based report): the n whose
 * log2(n / normal) is nearest (x0 taken at -2, x1 -1, x3 +0.58, x4 +1…; a male's X / Y: normal 1), written as any copy
 * number (also in the ISCN). `mosaic`: the fraction of cells (AChro-Puce) when a one-copy change is below 80 %.
 */
export function estCopies(c: CghCnv, sex: CghSex | undefined, build?: GenomeBuild | null): EstCopies | null {
  if (c.logRatio == null || c.logRatio === 0 || c.copies != null || c.amp || c.iscn) return null;
  const hemi = (c.chrom === 'chrX' || c.chrom === 'chrY') && !inPar(c, build) && (sex === 'male' || c.chrom === 'chrY');
  const base = hemi ? 1 : 2, lr = c.logRatio;
  const cand = c.kind === 'loss' ? Array.from({ length: base }, (_, n) => n) : [base + 1, base + 2, base + 3, base + 4];
  const exp = (n: number) => (n === 0 ? -2 : Math.log2(n / base));
  const n = cand.reduce((b, k) => (Math.abs(lr - exp(k)) < Math.abs(lr - exp(b)) ? k : b), cand[0]);
  const m = Math.abs(n - base) === 1 ? mosaicOf(c, sex, build) : null;
  const mosaic = m?.ok && m.f < MOSAIC_TAG_BELOW ? m.f : null;
  const top = c.kind === 'gain' && n === base + 4 ? '+' : '';
  return { n, expected: exp(n), hemi, mosaic, text: `x${n}${top}${mosaic != null ? ` mosaic ~${Math.round(mosaic * 100)}%` : ''}` };
}
/** the tooltip of an estimate */
export const estTip = (e: EstCopies, lr: number) => `Copy number from the log ratio ${lr.toFixed(2)} (the file gives none): nearest expected log2(n/${e.hemi ? 1 : 2}) is ${e.expected.toFixed(2)} for x${e.n}${e.mosaic != null ? `; a one-copy change in ~${Math.round(e.mosaic * 100)} % of cells (AChro-Puce: mosaic below 80 %)` : ''}.`;
