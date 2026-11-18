/**
 * Array-CGH: the exons of a gene's principal transcript that a call cuts, and what that does to the reading frame.
 * The call's ends (its first and last probes) are taken as the breakpoints: the true ones lie somewhere between probes,
 * so this is a theoretical reading, never a call.
 */
import { ranksText, type SvTranscript } from './svGenes';

export type ExonEffect = 'in-frame' | 'out-of-frame' | 'utr' | 'start-lost' | 'stop-lost' | 'start-stop-lost'
  | 'start-dup' | 'stop-dup' | 'undetermined' | 'truncated' | 'partial-copy' | 'non-coding' | 'intronic' | 'none';
export interface CnvExons {
  tx: string; strand: 1 | -1;
  /** exon count of the transcript */
  n: number;
  /** ranks wholly inside the call, ranks only partly inside */
  inside: number[]; partial: number[];
  /** coding bases of the exons wholly inside (null: non-coding transcript) */
  coding: number | null;
  where: 'intragenic' | '5p' | '3p';
  /** both ends in one intron: the exons around it */
  intron: [number, number] | null;
  effect: ExonEffect;
  /** short text for the gene table, and whether it reads as harmful */
  text: string; severe: boolean;
}

const overlap = (a0: number, a1: number, b0: number, b1: number) => Math.max(0, Math.min(a1, b1) - Math.max(a0, b0));
const codingOf = (tx: SvTranscript, e: { start: number; end: number }) =>
  tx.cdsStart == null || tx.cdsEnd == null ? 0 : overlap(e.start, e.end, tx.cdsStart, tx.cdsEnd);

/** What a call [s0, e0) (0-based, half-open) does to the exons of a transcript it cuts; null when it holds the whole transcript or misses it. */
export function exonsInCall(tx: SvTranscript, s0: number, e0: number, kind: 'loss' | 'gain'): CnvExons | null {
  if (e0 <= tx.start || s0 >= tx.end || (s0 <= tx.start && e0 >= tx.end)) return null;
  const n = tx.exons.length, plus = tx.strand > 0;
  const inside = tx.exons.filter(e => e.start >= s0 && e.end <= e0).map(e => e.rank).sort((a, b) => a - b);
  const partial = tx.exons.filter(e => overlap(e.start, e.end, s0, e0) > 0 && !(e.start >= s0 && e.end <= e0)).map(e => e.rank).sort((a, b) => a - b);
  const coding = tx.cdsStart == null || tx.cdsEnd == null || tx.cdsEnd <= tx.cdsStart ? null
    : tx.exons.filter(e => inside.includes(e.rank)).reduce((k, e) => k + codingOf(tx, e), 0);
  const leftIn = s0 > tx.start, rightIn = e0 < tx.end;
  // the side of the gene the call reaches past: its left end is the 5′ end on the + strand
  const where: CnvExons['where'] = leftIn && rightIn ? 'intragenic' : (!leftIn) === plus ? '5p' : '3p';
  const loss = kind === 'loss';
  let intron: CnvExons['intron'] = null;
  if (where === 'intragenic' && !inside.length && !partial.length) {
    const before = tx.exons.filter(e => e.end <= s0), after = tx.exons.filter(e => e.start >= e0);
    const a = before.sort((x, y) => y.end - x.end)[0], b = after.sort((x, y) => x.start - y.start)[0];
    if (a && b) intron = [Math.min(a.rank, b.rank), Math.max(a.rank, b.rank)];
  }
  const all = [...inside, ...partial].sort((a, b) => a - b);
  const run = inside.length ? `${ranksText(inside)} of ${n}` : '';
  const cut = partial.length ? `breakpoint inside exon ${partial.join(', ')}` : '';
  let effect: ExonEffect, label: string;
  if (!all.length) { effect = intron ? 'intronic' : 'none'; label = ''; }
  else if (where !== 'intragenic') {
    effect = coding == null ? 'non-coding' : loss ? 'truncated' : 'partial-copy';
    label = loss ? 'gene truncated' : 'partial copy, the gene itself may be intact (if in tandem or outside the gene)';
  } else if (coding == null) { effect = 'non-coding'; label = 'non-coding transcript'; }
  else if (partial.length) { effect = 'undetermined'; label = 'frame undetermined'; }
  else {
    const codingRanks = tx.exons.filter(e => codingOf(tx, e) > 0).map(e => e.rank);
    const first = Math.min(...codingRanks), last = Math.max(...codingRanks);
    const hasStart = inside.includes(first), hasStop = inside.includes(last);
    const tandem = loss ? '' : 'if in tandem: ';
    if (coding === 0) { effect = 'utr'; label = 'non-coding exons (UTR)'; }
    else if (hasStart || hasStop) {
      effect = loss ? (hasStart && hasStop ? 'start-stop-lost' : hasStart ? 'start-lost' : 'stop-lost') : hasStart ? 'start-dup' : 'stop-dup';
      label = loss ? (hasStart && hasStop ? 'start and stop codons lost' : hasStart ? 'start codon lost' : 'stop codon lost')
        : `${hasStart ? 'start' : 'stop'} codon exon duplicated, frame undetermined`;
    } else { effect = coding % 3 === 0 ? 'in-frame' : 'out-of-frame'; label = `${tandem}${coding % 3 === 0 ? 'in frame' : 'out of frame'}`; }
  }
  const side = where === '5p' ? `5′ part (${ranksText(all)} of ${n}, with the promoter)` : where === '3p' ? `3′ part (${ranksText(all)} of ${n})` : run;
  const text = [side, cut, label].filter(Boolean).join(' · ');
  const severe = ['out-of-frame', 'start-lost', 'stop-lost', 'start-stop-lost', 'truncated'].includes(effect);
  return { tx: tx.id, strand: tx.strand > 0 ? 1 : -1, n, inside, partial, coding, where, intron, effect, text, severe };
}
