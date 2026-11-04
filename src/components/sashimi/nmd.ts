/**
 * Whether a loss-of-function ClinVar variant is predicted to trigger nonsense-mediated mRNA decay (NMD), on one
 * transcript (the gene's MANE Select / principal one).
 *
 * The rule (Nagy & Maquat, Trends Biochem Sci 1998; the PVS1 decision tree, Abou Tayoun et al., Hum Mutat 2018): a
 * premature stop codon more than 50 nt upstream of the last exon-exon junction triggers NMD; one in the last exon, or in
 * the last 50 nt of the penultimate exon, escapes it. A single-exon transcript has no junction: no NMD.
 *
 * - nonsense: the stop is the variant's codon;
 * - frameshift: the new stop is `fsTerN` / `fs*N` codons on (N counts the first changed residue as 1), when the protein
 *   HGVS gives it; else only a variant already in the escape zone is known to escape (the new stop lies further on);
 * - canonical splice site: the adjacent exon is assumed skipped. In frame (its length a multiple of 3): the protein loses
 *   those residues, no premature stop, escape. Out of frame: the new stop is assumed in the exon after the skipped one,
 *   so it escapes only when that exon is the last;
 * - start lost: no premature stop, so no NMD (the protein may start at a downstream methionine).
 *
 * Without the sequence, the out-of-frame stop of a skipped exon and that of a frameshift without `TerN` are positions
 * assumed, not read: the call says so (`assumed`). Calls are on the genomic position, so they hold for the transcript
 * given, whatever transcript ClinVar's HGVS is written on.
 */
import type { SvTranscript } from './svGenes';
import type { ClinvarVariant } from './clinvar';

export interface NmdCall {
  /** 'nmd': predicted to trigger NMD; 'escape': predicted to escape it (or no premature stop) */
  call: 'nmd' | 'escape';
  /** a position or an outcome assumed rather than computed (no sequence, no transcript, no stop position given) */
  assumed: boolean;
  why: string;
}

/** The distance (nt) upstream of the last exon-exon junction within which a premature stop escapes NMD. */
export const NMD_BOUNDARY_NT = 50;

/** Exons in transcription order, with their start in transcript (mRNA) coordinates. */
function txLayout(tx: SvTranscript) {
  const exons = [...tx.exons].sort((a, b) => (tx.strand > 0 ? a.start - b.start : b.start - a.start));
  let at = 0;
  const out = exons.map(e => { const o = { ...e, m0: at, len: e.end - e.start }; at += e.end - e.start; return o; });
  return { exons: out, length: at };
}

/** mRNA coordinate of a genomic base (0-based), or null when it is not in an exon. */
function toMrna(lay: ReturnType<typeof txLayout>, strand: 1 | -1, g: number): number | null {
  for (const e of lay.exons) if (g >= e.start && g < e.end) return e.m0 + (strand > 0 ? g - e.start : e.end - 1 - g);
  return null;
}

/** `fsTer34`, `fs*34`, `fsX34`: the new stop's place, counting the first changed residue as 1. */
export function fsStop(protein: string | null): number | null {
  const m = /fs(?:Ter|\*|X)(\d+)/.exec(protein ?? '');
  return m ? Number(m[1]) : null;
}

export function predictNmd(v: ClinvarVariant, tx: SvTranscript | null | undefined): NmdCall {
  const c = (v.conseq || '').toLowerCase(), p = v.protein ?? '';
  if (/initiator codon|start lost/.test(c) || /^p\.\(?Met1[^0-9]/.test(p)) return { call: 'escape', assumed: false, why: 'start lost: no premature stop, no NMD (translation may start at a downstream methionine)' };
  if (!tx || tx.cdsStart == null || tx.cdsEnd == null) return { call: 'nmd', assumed: true, why: 'no coding transcript to place it on: NMD assumed' };
  const lay = txLayout(tx), last = lay.exons[lay.exons.length - 1];
  const id = tx.id;
  if (!last) return { call: 'nmd', assumed: true, why: `no exon structure for ${id}: NMD assumed` };
  if (lay.exons.length === 1) return { call: 'escape', assumed: false, why: `single-exon transcript (${id}): no exon junction, no NMD` };
  const lastEJ = last.m0, zone = lastEJ - NMD_BOUNDARY_NT;
  const cds0 = toMrna(lay, tx.strand, tx.strand > 0 ? tx.cdsStart : tx.cdsEnd - 1);
  const where = (m: number) => (m >= lastEJ ? 'in the last exon' : `${lastEJ - m} nt before the last exon junction`);
  // the variant's first base in transcription order
  const g = tx.strand > 0 ? v.start : v.end - 1;
  const splice = /splice (donor|acceptor)/.exec(c);
  if (splice) {
    // the exon whose boundary the site belongs to: a donor follows an exon's 3' end, an acceptor precedes its 5' start
    const donor = splice[1] === 'donor';
    let best: (typeof lay.exons)[number] | null = null, bestD = 11;
    for (const e of lay.exons) {
      const b = donor ? (tx.strand > 0 ? e.end : e.start - 1) : (tx.strand > 0 ? e.start - 1 : e.end);
      const d = Math.min(Math.abs(v.start - b), Math.abs(v.end - 1 - b));
      if (d < bestD) { bestD = d; best = e; }
    }
    if (!best) return { call: 'nmd', assumed: true, why: `splice site not next to an exon of ${id}: NMD assumed` };
    const rank = lay.exons.indexOf(best) + 1, n = lay.exons.length;
    if ((donor && rank === n) || (!donor && rank === 1)) return { call: 'nmd', assumed: true, why: `splice site of exon ${rank}/${n} with no exon to skip: NMD assumed` };
    const hasStart = cds0 != null && cds0 >= best.m0 && cds0 < best.m0 + best.len;
    const stopM = toMrna(lay, tx.strand, tx.strand > 0 ? tx.cdsEnd - 1 : tx.cdsStart);
    const hasStop = stopM != null && stopM >= best.m0 && stopM < best.m0 + best.len;
    if (hasStart || hasStop) return { call: 'nmd', assumed: true, why: `skipping exon ${rank}/${n} would remove the ${hasStart ? 'start' : 'stop'} codon: outcome not predicted, NMD assumed` };
    const cdsEndM = stopM ?? 0;
    if (cds0 != null && (best.m0 + best.len <= cds0 || best.m0 > cdsEndM)) return { call: 'nmd', assumed: true, why: `exon ${rank}/${n} is untranslated: outcome not predicted, NMD assumed` };
    if (best.len % 3 === 0) return { call: 'escape', assumed: true, why: `exon ${rank}/${n} skipping assumed: in frame (${best.len} nt, ${best.len / 3} residues lost), no premature stop` };
    if (rank + 1 === n) return { call: 'escape', assumed: true, why: `exon ${rank}/${n} skipping assumed: out of frame, the new stop falls in the last exon` };
    return { call: 'nmd', assumed: true, why: `exon ${rank}/${n} skipping assumed: out of frame, a new stop assumed before the last exon junction` };
  }
  const m = toMrna(lay, tx.strand, g);
  if (m == null) return { call: 'nmd', assumed: true, why: `not in an exon of ${id}: NMD assumed` };
  if (cds0 == null || m < cds0) return { call: 'nmd', assumed: true, why: `before the coding sequence of ${id}: NMD assumed` };
  const codon = cds0 + Math.floor((m - cds0) / 3) * 3;
  if (/frameshift/.test(c) || /fs/.test(p)) {
    const n = fsStop(p);
    if (n != null) {
      const stop = codon + (n - 1) * 3;
      if (stop >= lay.length) return { call: 'escape', assumed: false, why: `frameshift reading past the transcript's end (stop ${n} codons on): no NMD` };
      return stop >= zone ? { call: 'escape', assumed: false, why: `frameshift, new stop ${n} codons on, ${where(stop)}: escapes NMD` }
        : { call: 'nmd', assumed: false, why: `frameshift, new stop ${n} codons on, ${where(stop)}: NMD` };
    }
    return codon >= zone ? { call: 'escape', assumed: false, why: `frameshift ${where(codon)}: the new stop lies further on, escapes NMD` }
      : { call: 'nmd', assumed: true, why: `frameshift ${where(codon)}, new stop not given: NMD assumed` };
  }
  // nonsense (or a truncating change the protein HGVS shows as a stop)
  return codon >= zone ? { call: 'escape', assumed: false, why: `stop ${where(codon)}: escapes NMD` } : { call: 'nmd', assumed: false, why: `stop ${where(codon)}: NMD` };
}
