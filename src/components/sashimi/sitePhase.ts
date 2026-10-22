/**
 * The haplotype of a variant site in a file phased by a tool (HP / PS haplotags), from its reads by haplotag
 * (VariantSite.hap, counted by the variants scan): within the site's phase set, the share of HP 1 reads and of HP 2 reads
 * that carry the allele.
 *
 * - on H1 (or H2): at least PHASE_HI of that haplotype's reads carry it, at most PHASE_LO of the other's;
 * - on both: at least PHASE_HI of each (homozygous);
 * - on neither: at most PHASE_LO of each; the allele sits in reads of no haplotype of the set (untagged reads, another
 *   phase set's), or it is too rare;
 * - check: anything else. Most often present in part of one haplotype's reads (a mosaic, post-zygotic change; an
 *   artefact; an indel in a homopolymer that noisy long reads only partly carry), or on both haplotypes without being
 *   homozygous (a mis-phased site, a collapsed duplication or paralog, a third haplotype);
 * - not judged: fewer than PHASE_MIN_DEPTH reads on either haplotype.
 *
 * PHASE_HI is 0.7 rather than the 0.8 of the haplotype consensus check (haplotypes.ts), so the partial indels of long reads
 * in homopolymers still land on their haplotype; the fractions are always shown.
 *
 * The two haplotypes of different phase sets are not linked: H1 of one set is not known to be H1 of the next one.
 */
import type { VariantSite } from './types';
import { fisher } from './siteQuality';

export const PHASE_MIN_DEPTH = 3, PHASE_HI = 0.7, PHASE_LO = 0.2;

export type PhaseCall = 'h1' | 'h2' | 'both' | 'neither' | 'check' | 'thin';
export interface SitePhase {
  call: PhaseCall;
  /** share of HP 1 / HP 2 reads with the allele (null without reads) */
  f: [number | null, number | null];
  /** two-sided Fisher exact test, alternate vs other reads × HP 1 vs HP 2 */
  p: number;
  /** the site's phase set (PS), null when the reads carry none */
  ps: number | null;
  /** in words, for the hover and the panel */
  note: string;
}

const pc = (f: number | null) => (f == null ? '–' : `${Math.round(f * 100)} %`);

export function sitePhase(s: VariantSite): SitePhase | null {
  const h = s.hap;
  if (!h) return null;
  const [a1, a2] = h.alt, [d1, d2] = h.depth;
  const f: [number | null, number | null] = [d1 ? a1 / d1 : null, d2 ? a2 / d2 : null];
  const p = fisher(a1, d1 - a1, a2, d2 - a2, 'two');
  const shares = `${pc(f[0])} of HP 1 reads (${a1}/${d1}), ${pc(f[1])} of HP 2 reads (${a2}/${d2})`;
  if (d1 < PHASE_MIN_DEPTH || d2 < PHASE_MIN_DEPTH)
    return { call: 'thin', f, p, ps: h.ps, note: `not judged: fewer than ${PHASE_MIN_DEPTH} reads on ${d1 < PHASE_MIN_DEPTH && d2 < PHASE_MIN_DEPTH ? 'either haplotype' : `HP ${d1 < PHASE_MIN_DEPTH ? 1 : 2}`} · ${shares}` };
  const [x1, x2] = [f[0]!, f[1]!];
  if (x1 >= PHASE_HI && x2 >= PHASE_HI) return { call: 'both', f, p, ps: h.ps, note: `on both haplotypes (homozygous) · ${shares}` };
  if (x1 >= PHASE_HI && x2 <= PHASE_LO) return { call: 'h1', f, p, ps: h.ps, note: `on H1 · ${shares}` };
  if (x2 >= PHASE_HI && x1 <= PHASE_LO) return { call: 'h2', f, p, ps: h.ps, note: `on H2 · ${shares}` };
  if (x1 <= PHASE_LO && x2 <= PHASE_LO) {
    const where = h.otherAlt ? `${h.otherAlt} of the allele's reads are tagged in another phase set` : 'the allele is mostly in untagged reads';
    return { call: 'neither', f, p, ps: h.ps, note: `on neither haplotype of the phase set: ${where} · ${shares}` };
  }
  const why = x1 > PHASE_LO && x2 > PHASE_LO
    ? 'on both haplotypes without being homozygous: a mis-phased site, a collapsed duplication or paralog, a third haplotype?'
    : `in part of the HP ${x1 > x2 ? 1 : 2} reads only: a mosaic (post-zygotic) change, an artefact, or an indel the reads carry in part (homopolymer)?`;
  return { call: 'check', f, p, ps: h.ps, note: `check: ${why} · ${shares}` };
}
