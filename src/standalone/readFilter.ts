/**
 * Read filters of a reads request (ReadsOptions), applied by every source before its cap, so that the reads kept are
 * drawn from all those that pass rather than from a sample of the window:
 *
 * - **read length** (`minReadLength`): the read's own length, as sequenced: aligned bases, inserted, soft- and
 *   hard-clipped (the bases a supplementary record leaves to its primary). The same from a record's CIGAR and from a
 *   decoded read; introns (N) and deletions are not read bases, so a spliced RNA read counts its bases only.
 * - **supporting reads** (`support`): an arc's (arcSupport.ts), or a variant's (kind 'variant'): the reads carrying its
 *   alternate allele as the site counts them (phasing.ts alleleAt): an SNV base of quality ≥ MIN_BQ, an insertion or a
 *   deletion of the site's position and length.
 */
import type { AlignedRead, VariantSite } from '../components/sashimi/types';
import type { ReadsOptions } from '../components/sashimi/datasource';
import type { RawRead } from './alignments';
import { arcReadFromAligned, arcReadFromCigar, supportsArc, type ArcSupport } from './arcSupport';
import { alleleAt } from './phasing';
import { MIN_BQ } from './alleles';
import { plusStrand } from './coverage';

/** Length of a decoded read, as sequenced. */
export function readLength(r: AlignedRead): number {
  let n = r.c[0] + r.c[1] + (r.h?.[0] ?? 0) + (r.h?.[1] ?? 0);
  for (const [a, b] of r.b) n += b - a;
  for (const [, l] of r.i) n += l;
  return n;
}

/** Length of a record from its packed CIGAR (length << 4 | op, as BAM stores it): M I S H = X. */
export function opsReadLength(ops: ArrayLike<number>): number {
  let n = 0;
  for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 1 || op === 4 || op === 5 || op === 7 || op === 8) n += ops[k] >>> 4; }
  return n;
}

/** Length of a record, as sequenced (from its CIGAR: M = X I S H; the sequence may be absent, '*', on secondary records). */
export function rawReadLength(r: Pick<RawRead, 'cigar' | 'seq'>): number {
  let n = 0;
  for (const [, len, op] of r.cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) if ('M=XISH'.includes(op)) n += +len;
  return n || r.seq.length;
}

/** the site a variant's support stands for, as alleleAt reads it */
const siteOf = (a: ArcSupport): VariantSite | null => (a.allele
  ? { pos: a.start, kind: a.allele.kind, ref: '', alt: a.allele.alt, length: a.allele.length, alt_count: 0, depth: 0, vaf: 0 }
  : null);

/** Whether a record carries the alternate allele of a variant: its CIGAR walked to the site. */
export function rawCarriesAllele(r: Pick<RawRead, 'start' | 'cigar' | 'seq' | 'qual' | 'mismatches'>, s: VariantSite, minBq = MIN_BQ): boolean {
  let rpos = r.start, qpos = 0;
  for (const [, n, op] of r.cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) {
    const len = +n;
    if (rpos > s.pos + 1) break;
    if (op === 'M' || op === '=' || op === 'X') {
      if (s.kind === 'snv' && s.pos >= rpos && s.pos < rpos + len) {
        const q = qpos + (s.pos - rpos);
        if (r.seq) return r.seq[q]?.toUpperCase() === s.alt && (r.qual && r.qual.length > q ? r.qual[q] : 30) >= minBq;
        const m = r.mismatches?.find(x => x[0] === s.pos);
        return !!m && m[1].toUpperCase() === s.alt && m[2] >= minBq;
      }
      rpos += len; qpos += len;
    } else if (op === 'I') {
      if (s.kind === 'ins' && rpos === s.pos && len === s.length) return true;
      qpos += len;
    } else if (op === 'D') {
      if (s.kind === 'del' && rpos === s.pos && len === s.length) return true;
      rpos += len;
    } else if (op === 'N') rpos += len;
    else if (op === 'S') qpos += len;
  }
  return false;
}

/** Whether a decoded read supports an arc or carries a variant. */
export function supportsAligned(r: AlignedRead, chrom: string, a: ArcSupport): boolean {
  const s = siteOf(a);
  if (!s) return supportsArc(arcReadFromAligned(r, chrom), chrom, a);
  if (a.strand && plusStrand(r.f) !== (a.strand === 'plus')) return false;
  return alleleAt(r, s, MIN_BQ) === 1;
}

/** Whether a record aligns a base at `pos` (an M, = or X operation over it): the CIGAR alone. */
export function rawAlignsBase(r: Pick<RawRead, 'start' | 'cigar'>, pos: number): boolean {
  let rpos = r.start;
  for (const [, n, op] of r.cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) {
    const len = +n;
    if (rpos > pos) return false;
    if (op === 'M' || op === '=' || op === 'X') { if (pos < rpos + len) return true; rpos += len; }
    else if (op === 'D' || op === 'N') rpos += len;
  }
  return false;
}

/**
 * Whether a record supports an arc or carries a variant. A record in its light form (scans decode no bases before
 * filtering) is decoded whole with `full` for an SNV, once its CIGAR puts a base at the site.
 */
export function supportsRaw(r: RawRead, chrom: string, a: ArcSupport, full?: () => RawRead): boolean {
  const s = siteOf(a);
  if (!s) return supportsArc(arcReadFromCigar(r), chrom, a);
  if (a.strand && plusStrand(r.flags) !== (a.strand === 'plus')) return false;
  if (s.kind === 'snv' && !r.seq && !r.mismatches && full) return rawAlignsBase(r, s.pos) && rawCarriesAllele(full(), s);
  return rawCarriesAllele(r, s);
}

/** A variant's reads are those carrying it: their mates are not added (an arc's are, to show its pairs whole). */
export const addsMates = (a: ArcSupport | undefined): boolean => !!a && a.kind !== 'variant';

/** The length floor of a request (undefined: none). */
export const minLengthOf = (opts?: Pick<ReadsOptions, 'minReadLength'>): number | undefined => (opts?.minReadLength && opts.minReadLength > 0 ? opts.minReadLength : undefined);
