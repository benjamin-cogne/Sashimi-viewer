/**
 * Where a call lies on its chromosome: terminal (pter / qter), a whole arm, a whole chromosome, and, per patient, the
 * pairs of terminal losses and gains that suggest a rearrangement behind them. Only hints: a "possible …" to check by
 * karyotype / FISH, never a class and never part of the priority.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { chromLength } from '../../standalone/assembly';
import type { CghCnv } from './cghIntervals';
import type { CghState } from './cghStore';

/** a call this close to a chromosome end is terminal: arrays carry no probes in the telomeric repeats */
export const TERMINAL_BP = 500_000;
/** the share of an arm a call must cover to be a whole arm */
export const ARM_COVER = 0.9;
const ACROCENTRIC = new Set(['13', '14', '15', '21', '22']);

export interface CnvTopo {
  /** the chromosome ends the call lies on: p = pter, q = qter */
  ends: ('p' | 'q')[];
  /** a whole arm, or the whole chromosome */
  span?: 'p' | 'q' | 'chr';
}
export const topoLabel = (t: CnvTopo | undefined) => (t?.span === 'chr' ? 'chr' : t?.span ? 'arm' : t?.ends.length ? 'ter' : '');
export const topoWords = (t: CnvTopo | undefined, chrom: string) => {
  const n = chrom.replace(/^chr/, '');
  const out: string[] = [];
  if (t?.span === 'chr') out.push(`whole chromosome ${n}`);
  else if (t?.span) out.push(`whole ${t.span} arm of ${n} (≥ ${Math.round(ARM_COVER * 100)} % of it)`);
  for (const e of t?.ends ?? []) out.push(`terminal ${e}ter (within ${TERMINAL_BP / 1000} kb of the chromosome end)`);
  return out;
};

/** `bands` are the cytobands (0-based), null when they could not be read: then the arms are unknown, the ends are not */
export function topoOf(c: Pick<CghCnv, 'chrom' | 'start' | 'end'>, build: GenomeBuild, bands: { start: number; end: number; stain: string }[] | null): CnvTopo | undefined {
  const len = chromLength(build, c.chrom);
  if (!len) return undefined;
  const s0 = c.start - 1, e0 = c.end;
  const ends: CnvTopo['ends'] = [];
  if (s0 <= TERMINAL_BP) ends.push('p');
  if (e0 >= len - TERMINAL_BP) ends.push('q');
  let span: CnvTopo['span'];
  const acen = bands?.filter(b => b.stain === 'acen');
  if (acen?.length) {
    const pEnd = Math.min(...acen.map(b => b.start)), qStart = Math.max(...acen.map(b => b.end));
    const cover = (a: number, z: number) => Math.max(0, Math.min(e0, z) - Math.max(s0, a)) / Math.max(1, z - a);
    const acro = ACROCENTRIC.has(c.chrom.replace(/^chr/, ''));
    const p = !acro && cover(0, pEnd) >= ARM_COVER, q = cover(qStart, len) >= ARM_COVER;
    // an acrocentric's p arm is not read: its q arm alone is the chromosome
    span = (p || acro) && q ? 'chr' : p ? 'p' : q ? 'q' : undefined;
  }
  return ends.length || span ? { ends, span } : undefined;
}

export const HINT_TRANSLOCATION = 'terminal loss + terminal gain in this patient: a derivative of an unbalanced translocation is possible (parental karyotype / FISH)';
export const HINT_PERICENTRIC = 'terminal loss and terminal gain on opposite arms of the same chromosome: a recombinant of a pericentric inversion is possible (parental karyotype / FISH)';
export const HINT_INVDUPDEL = 'terminal loss with an adjacent gain on the same end: an inverted duplication-deletion is possible (parental karyotype / FISH)';

/** per call id, the hints of its patient's terminal losses and gains (among the analysed calls the annotation has read) */
export function rearrangementHints(s: CghState): Map<string, string[]> {
  const out = new Map<string, string[]>();
  const add = (c: CghCnv, h: string) => { const l = out.get(c.id) ?? []; if (!l.includes(h)) l.push(h); out.set(c.id, l); };
  for (const fl of s.files) for (const sm of fl.samples) {
    const read = sm.cnvs.filter(c => s.ann[c.id]?.state === 'done');
    // a whole chromosome is an aneuploidy, not a terminal event
    const ter = (c: CghCnv) => { const t = s.ann[c.id]?.topo; return t && t.span !== 'chr' ? t.ends : []; };
    const losses = read.filter(c => c.kind === 'loss' && ter(c).length), gains = read.filter(c => c.kind === 'gain');
    for (const l of losses) for (const g of gains) {
      const le = ter(l), ge = ter(g);
      if (l.chrom !== g.chrom) { if (ge.length) { add(l, HINT_TRANSLOCATION); add(g, HINT_TRANSLOCATION); } continue; }
      if (ge.length && le.some(e => !ge.includes(e)) && ge.some(e => !le.includes(e))) { add(l, HINT_PERICENTRIC); add(g, HINT_PERICENTRIC); }
      // the gain starts where the loss ends (or the reverse) within the terminal distance, on the loss's end
      const gap = le.includes('p') ? g.start - l.end : l.start - g.end;
      if (le.length === 1 && gap >= -1 && gap <= TERMINAL_BP) { add(l, HINT_INVDUPDEL); add(g, HINT_INVDUPDEL); }
    }
  }
  return out;
}
