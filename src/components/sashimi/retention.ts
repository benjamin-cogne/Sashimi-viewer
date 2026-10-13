/**
 * A retained intron on the transcript: its HGVS description, and what the retention does to the coding sequence.
 *
 * - HGVS (RNA level): the intron's bases inserted between the two exon bases around it,
 *   r.<last exon base>_<first exon base>ins<last exon base>+1_<first exon base>-1, numbered as the c. (or n.) positions
 *   of the transcript (https://hgvs-nomenclature.org, RNA insertions: an intron retained is described as the insertion
 *   of the intronic sequence).
 * - Consequence, for an intron inside the CDS: in frame when its length is a multiple of 3, else a frameshift; the
 *   first stop codon read in frame from the exon base before the intron, through the intron and into the next exon. A
 *   stop with exon junctions still downstream (not the last intron) makes the transcript a likely target of
 *   nonsense-mediated decay (the 50–55 nt rule); in the last intron, it does not.
 */
import { cdnaPosition, translateCodon, type TxModel } from './geometry';

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const revComp = (s: string) => { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] ?? 'N'; return o; };

/** Bases of the exons read on each side of the intron with it, for the codon split by the intron. */
export const RETENTION_FLANK = 3;

export interface RetentionHgvs {
  /** e.g. NM_000492.4:r.1584_1585ins1584+1_1585-1 */
  r: string;
  /** the c. (n.) numbering of the exon bases around the intron, transcript order */
  before: string; after: string;
  /** where the intron sits in the transcript */
  region: 'cds' | 'utr5' | 'utr3' | 'noncoding' | 'mixed';
  /** introns of the transcript in transcription order: this one's index (0 = first) and their number */
  index: number; count: number;
}

/** The HGVS r. description of retaining the intron [start, end) (0-based half-open) of `tx`, or null when it is not one of its introns. */
export function retentionHgvs(start: number, end: number, tx: TxModel): RetentionHgvs | null {
  const ex = tx.exons;
  const k = ex.findIndex((e, i) => i + 1 < ex.length && e.end === start && ex[i + 1].start === end);
  if (k < 0) return null;
  const plus = tx.strand > 0;
  // the exon base before the intron and the one after it, in transcription order
  const a = cdnaPosition(plus ? start - 1 : end, tx), b = cdnaPosition(plus ? end : start - 1, tx);
  const region = a.kind === 'noncoding' ? 'noncoding' : a.kind === b.kind || (a.kind === 'cds' && b.kind === 'cds') ? (a.kind as RetentionHgvs['region']) : 'mixed';
  const n = ex.length - 1;
  return {
    // r. on coding and non-coding transcripts alike (an n. transcript numbers its r. positions as n.)
    r: `${tx.transcriptId}:r.${a.bare}_${b.bare}ins${a.bare}+1_${b.bare}-1`,
    before: a.label, after: b.label, region,
    index: plus ? k : n - 1 - k, count: n,
  };
}

export interface RetentionEffect {
  frame: 'in frame' | 'frameshift';
  /** the first stop codon read in frame: codons from the first one holding an intron base, and whether it lies in the intron */
  stop: { codon: number; inIntron: boolean } | null;
  text: string;
}

/**
 * What retaining a coding intron does to the protein. `seq` is the genomic sequence (+ strand) of
 * [start - RETENTION_FLANK, end + RETENTION_FLANK); `h` its description (the intron's place and the c. numbering).
 */
export function retentionEffect(start: number, end: number, seq: string, h: RetentionHgvs, tx: TxModel): RetentionEffect | null {
  if (h.region !== 'cds') return null;
  const L = end - start;
  if (seq.length !== L + 2 * RETENTION_FLANK) return null;
  const t = (tx.strand > 0 ? seq : revComp(seq)).toUpperCase();
  // c. number of the exon base before the intron: that many bases of the CDS read, so (n mod 3) of the codon the intron splits
  const before = parseInt(h.before.replace(/^c\./, ''), 10);
  if (!Number.isFinite(before)) return null;
  const inCodon = before % 3;
  // codons from the first one holding an intron base: it starts `inCodon` bases before the intron
  let stop: RetentionEffect['stop'] = null;
  for (let i = RETENTION_FLANK - inCodon, c = 1; i + 3 <= t.length; i += 3, c++) {
    if (translateCodon(t.slice(i, i + 3)) === '*') { stop = { codon: c, inIntron: i < RETENTION_FLANK + L }; break; }
  }
  const frame = L % 3 === 0 ? 'in frame' : 'frameshift';
  const last = h.index === h.count - 1;
  const nmd = last ? 'in the last intron: no exon junction downstream, nonsense-mediated decay not expected' : 'exon junctions downstream: nonsense-mediated decay likely (the 50–55 nt rule)';
  let text: string;
  if (stop && stop.inIntron) text = `${frame} (${L.toLocaleString('en-US')} nt); a stop codon ${stop.codon} codon${stop.codon > 1 ? 's' : ''} into the retained intron: premature, ${nmd}`;
  else if (frame === 'in frame') text = `in frame: ${(L / 3).toLocaleString('en-US')} amino acids inserted, no stop codon in the intron`;
  else text = `frameshift (${L.toLocaleString('en-US')} nt, not a multiple of 3), no stop codon in the intron: the shifted frame runs on into the next exon`;
  return { frame, stop, text };
}
