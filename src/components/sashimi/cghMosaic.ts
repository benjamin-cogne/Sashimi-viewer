/**
 * Fraction of cells carrying a one-copy change, from the mean log2 ratio of the call: a port of the AChro-Puce
 * mosaicism calculator (https://acpa-achropuce.com/calcul-du-mosaicisme/, tool of July 2025), formulas unchanged.
 * Indicative only: the calculator ignores the biological variability of the ratio.
 */
import type { CghCnv, CghSex } from './cghIntervals';

/** below this fraction a call is tagged mosaic in the table */
export const MOSAIC_TAG_BELOW = 0.8;

/** the calculator's own function: chrom 'autosome' | 'X' | 'Y', sex 'F' (anything else: male), event 'del' | 'dup' */
export function computeMosaicism(lr: number, event: 'del' | 'dup', chrom: 'autosome' | 'X' | 'Y', sexe: 'F' | 'M'): number {
  let f = 0;
  const hemi = (chrom === 'X' || chrom === 'Y') && sexe !== 'F';
  if (event === 'del') f = hemi ? 1 - Math.pow(2, lr) : 2 - 2 * Math.pow(2, lr);
  else f = hemi ? Math.pow(2, lr) - 1 : Math.pow(2, lr) * 2 - 2;
  return Math.min(1, Math.max(0, f));
}

/** pseudo-autosomal regions, GRCh37 (1-based inclusive) */
const PAR: Record<string, [number, number][]> = {
  chrX: [[60001, 2699520], [154931044, 155260560]],
  chrY: [[10001, 2649520], [59034050, 59363566]],
};
const inPar = (c: CghCnv) => (PAR[c.chrom] ?? []).some(([a, b]) => c.start >= a && c.end <= b);

export type Mosaic =
  | { ok: true; f: number; lo: number | null; hi: number | null; expected: number; hemi: boolean; why?: undefined }
  | { ok: false; why: 'sex' };

/** the estimate for a call, or null when none applies (no log ratio, not a one-copy change); `sex needed` for chrX / chrY outside the PARs */
export function mosaicOf(c: CghCnv, sex: CghSex | undefined, dlrs?: number): Mosaic | null {
  if (c.logRatio == null || c.amp || c.copiesHi != null) return null;
  const sexChrom = (c.chrom === 'chrX' || c.chrom === 'chrY') && !inPar(c);
  if (sexChrom && !sex) return { ok: false, why: 'sex' };
  const hemi = sexChrom && sex === 'male';
  // a male's single X / Y: the one-copy change is x0 (loss) or x2 (gain); elsewhere x1 or x3
  const base = hemi ? 1 : 2;
  if (c.copies != null && c.copies !== (c.kind === 'loss' ? base - 1 : base + 1)) return null;
  const chrom = sexChrom ? (c.chrom === 'chrX' ? 'X' : 'Y') : 'autosome';
  const ev = c.kind === 'loss' ? 'del' : 'dup', sx = sex === 'male' ? 'M' : 'F';
  const f = computeMosaicism(c.logRatio, ev, chrom, sx);
  const d = dlrs != null && dlrs > 0 ? dlrs : null;
  const a = d != null ? computeMosaicism(c.logRatio - d, ev, chrom, sx) : null, b = d != null ? computeMosaicism(c.logRatio + d, ev, chrom, sx) : null;
  return { ok: true, f, lo: a != null && b != null ? Math.min(a, b) : null, hi: a != null && b != null ? Math.max(a, b) : null,
    expected: Math.log2((c.kind === 'loss' ? base - 1 : base + 1) / base), hemi };
}
