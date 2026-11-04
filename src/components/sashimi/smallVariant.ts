/**
 * A small variant typed as genomic HGVS (`chr2:g.60461323del`, `2:g.123A>G`, `NC_000002.12:g.100_102delinsTT`): a
 * substitution, or a deletion, duplication, insertion, deletion-insertion or inversion under SMALL_MAX bp. The SV
 * browser marks it on the reference instead of rebuilding a chromosome; a larger one is read as a rearrangement.
 * A transcript description (`NM_022893.4:c.1592del`) is first moved to the genome (Ensembl Variant Recoder).
 */
import { chromOfAccession } from './knownVariants';

/** Under this size (bp deleted, duplicated, inverted or inserted), a variant is marked, not rebuilt. */
export const SMALL_MAX = 50;

export interface SmallVariant {
  chrom: string;
  /** 0-based half-open stretch it touches (an insertion: the two bases around it) */
  start: number; end: number;
  kind: 'snv' | 'del' | 'dup' | 'ins' | 'delins' | 'inv';
  /** bp deleted, duplicated, inverted or inserted (the larger for a deletion-insertion; 1 for a substitution) */
  size: number;
  /** the description as typed, with chrN for the reference */
  hgvs: string;
}

const REF_RE = /^(?:(NC_\d+)(?:\.\d+)?|chr([0-9XYM]{1,2})|([0-9XYM]{1,2}))(?:\([^)]*\))?:g\.(.+)$/i;

/** An inserted sequence's length: a sequence, a count (`ins12`, `insN[12]`); null for another stretch of the genome. */
function insLength(s: string): number | null {
  if (/^[ACGTN]+$/i.test(s)) return s.length;
  const n = s.match(/^(?:N\[)?(\d+)\]?$/i);
  return n ? Number(n[1]) : null;
}

/** The small variant a genomic HGVS describes, or null (not genomic HGVS, or SMALL_MAX bp or more: a rearrangement). */
export function parseSmallGenomic(raw: string): SmallVariant | null {
  const m = raw.trim().replace(/\s+/g, '').match(REF_RE);
  if (!m) return null;
  const c = m[1] ? chromOfAccession(m[1].toUpperCase()) : (m[2] ?? m[3]).toUpperCase();
  if (!c) return null;
  const chrom = `chr${c === 'MT' ? 'M' : c}`, body = m[4].replace(/,/g, '');
  const hgvs = `${chrom}:g.${body}`;
  let v: Omit<SmallVariant, 'chrom' | 'hgvs'> | null = null;
  let x: RegExpMatchArray | null;
  if ((x = body.match(/^(\d+)[ACGTN]>[ACGTN]$/i))) v = { start: +x[1] - 1, end: +x[1], kind: 'snv', size: 1 };
  else if ((x = body.match(/^(\d+)(?:_(\d+))?delins(.+)$/i))) {
    const a = +x[1], b = x[2] ? +x[2] : a, n = insLength(x[3]);
    if (n == null) return null;
    v = { start: a - 1, end: b, kind: 'delins', size: Math.max(b - a + 1, n) };
  } else if ((x = body.match(/^(\d+)(?:_(\d+))?del(?:[ACGTN]*|\d*)$/i))) { const a = +x[1], b = x[2] ? +x[2] : a; v = { start: a - 1, end: b, kind: 'del', size: b - a + 1 }; }
  else if ((x = body.match(/^(\d+)(?:_(\d+))?dup(?:[ACGTN]*|\d*)$/i))) { const a = +x[1], b = x[2] ? +x[2] : a; v = { start: a - 1, end: b, kind: 'dup', size: b - a + 1 }; }
  else if ((x = body.match(/^(\d+)_(\d+)ins(.+)$/i))) {
    const n = insLength(x[3]);
    if (n == null) return null;
    v = { start: +x[1] - 1, end: +x[2], kind: 'ins', size: n };
  } else if ((x = body.match(/^(\d+)_(\d+)inv$/i))) { const a = +x[1], b = +x[2]; v = { start: a - 1, end: b, kind: 'inv', size: b - a + 1 }; }
  if (!v || v.end <= v.start || v.size >= SMALL_MAX) return null;
  return { chrom, hgvs, ...v };
}

/** A transcript description to move to the genome: `NM_022893.4:c.1592del`, `NM_…(BCL11A):c.…`, `ENST…:c.…`, `NR_…:n.…`. */
export const TRANSCRIPT_HGVS_RE = /^(?:N[MR]_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?|LRG_\d+t\d+)(?:\([^)]*\))?:[cn]\.\S+$/i;

export const SMALL_KIND_WORD: Record<SmallVariant['kind'], string> = { snv: 'substitution', del: 'deletion', dup: 'duplication', ins: 'insertion', delins: 'deletion-insertion', inv: 'inversion' };
