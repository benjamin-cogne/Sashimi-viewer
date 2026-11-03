/**
 * ClinVar variants classified pathogenic or likely pathogenic, from the UCSC Genome Browser's ClinVar tracks (bigBed,
 * autoSql `clinVarBed` of UCSC's kent source, src/hg/utils/otto/clinvar/clinvar.as): `clinvarCnv` holds the variants of
 * 50 bp or more, `clinvarMain` the shorter ones (hg19 and hg38). `_clinSignCode` is PG (pathogenic) or LP (likely
 * pathogenic, also for "Pathogenic/Likely pathogenic") for the ones kept; conflicting interpretations are CF, left out.
 */
export type ClinvarKind = 'loss' | 'gain' | 'inv' | 'complex' | 'snv' | 'indel' | 'other';
export interface ClinvarVariant {
  /** the ClinVar VCV accession (or the item's name) */
  id: string;
  /** the ClinVar variation id (for the link) */
  variantId: number | null;
  /** 0-based half-open */
  start: number; end: number;
  /** PG or LP, and ClinVar's wording ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") */
  code: 'PG' | 'LP'; sign: string;
  kind: ClinvarKind; type: string;
  genes: string[];
  /** the conditions, as ClinVar names them */
  conditions: string;
  /** review stars (0 to 4), and the review status in words */
  stars: number; review: string;
  /** allele origin: germline, somatic, both, unknown */
  origin: string;
  /** a nucleotide HGVS, and a protein one when ClinVar gives it */
  hgvs: string | null; protein: string | null;
  len: number;
}

const KIND: Record<string, ClinvarKind> = { LOSS: 'loss', GAIN: 'gain', INV: 'inv', STRUCT: 'complex', SUBST: 'snv', INS: 'indel', SEQALT: 'other', SEQLEN: 'indel' };
const ORIGIN: Record<string, string> = { GERM: 'germline', SOM: 'somatic', GERMSOM: 'germline and somatic', UNK: 'unknown', NOVO: 'de novo' };
const strip = (s: unknown) => String(s ?? '').replace(/<[^>]*>/g, '').replace(/&nbsp;/g, ' ').replace(/\s+/g, ' ').trim();

/** One item of a UCSC ClinVar track, or null when it is not pathogenic / likely pathogenic (or unreadable). */
export function parseClinvarItem(it: Record<string, unknown>): ClinvarVariant | null {
  const code = String(it._clinSignCode ?? '').toUpperCase();
  if (code !== 'PG' && code !== 'LP') return null;
  const start = Number(it.chromStart), end = Number(it.chromEnd);
  if (!Number.isFinite(start) || !Number.isFinite(end)) return null;
  const typeCode = String(it._allTypeCode ?? '').toUpperCase();
  let kind = KIND[typeCode] ?? 'other';
  const type = strip(it.type);
  if (kind === 'other' && /deletion|loss/i.test(type)) kind = 'loss';
  else if (kind === 'other' && /duplication|gain/i.test(type)) kind = 'gain';
  // "GeneID|SYMBOL" (or the symbol alone), several separated by commas
  const genes = String(it.geneId ?? '').split(/[,;]\s*/).map(g => g.split('|').pop()!.trim()).filter(g => g && !/^\d+$/.test(g));
  let hgvs: string | null = null, protein: string | null = null;
  try { const t = JSON.parse(String(it._jsonHgvsTable || '[]')); if (Array.isArray(t) && Array.isArray(t[0])) { hgvs = t[0][0] ? String(t[0][0]) : null; protein = t[0][1] ? String(t[0][1]) : null; } } catch { /* none */ }
  const review = strip(it.reviewStatus).replace(/^.*?based on:\s*/i, '');
  const vid = Number(it._variantId);
  return {
    id: String(it.vcvId || it.name || `${start}-${end}`), variantId: Number.isFinite(vid) && vid > 0 ? vid : null,
    start, end, code, sign: strip(it.clinSign) || (code === 'PG' ? 'Pathogenic' : 'Likely pathogenic'),
    kind, type, genes, conditions: strip(it.phenotypeList), stars: Number(it._starCount ?? it.score) || 0, review,
    origin: ORIGIN[String(it._originCode ?? '').toUpperCase()] ?? strip(it.origin), hgvs, protein, len: Number(it._varLen) || end - start,
  };
}

/** ClinVar's page of a variant. */
export const clinvarUrl = (v: ClinvarVariant) => v.variantId ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.variantId}/` : `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(v.id)}`;
