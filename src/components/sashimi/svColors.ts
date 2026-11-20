/**
 * Colours shared by the SV cartoon and the rearranged genome browser: one per segment (the flanks in their
 * chromosome's colour, inner segments from a palette), one per junction kind, one per gene impact.
 */
import type { SvReconstruction, SvJunction } from './svModel';
import type { ImpactKind } from './svGenes';

export const FLANK_COLORS = ['#64748b', '#7c3aed', '#0891b2'];
export const SEG_COLORS = ['#2a78d6', '#eda100', '#1baf7a', '#e87ba4', '#008300', '#eb6834', '#a16207', '#0ea5e9'];
export const KIND_COLOR: Record<SvJunction['kind'], string> = { deletion: '#b91c1c', duplication: '#15803d', inversion: '#2563eb', translocation: '#d97706' };

/** The colour of every segment of a reconstruction, by segment id. */
export function segmentColors(rec: SvReconstruction): Map<number, string> {
  const out = new Map<number, string>();
  let ci = 0;
  for (const s of rec.segs) out.set(s.id, s.first || s.last ? FLANK_COLORS[rec.chrs.indexOf(s.chrom) % FLANK_COLORS.length] : SEG_COLORS[ci++ % SEG_COLORS.length]);
  return out;
}

export const IMPACT: Record<ImpactKind, { color: string; label: string; short: string }> = {
  fusion: { color: '#7c3aed', label: 'fusion', short: 'fusion' },
  disrupted: { color: '#b91c1c', label: 'disrupted (cut, no intact copy)', short: 'disrupted' },
  deleted: { color: '#dc2626', label: 'whole gene deleted', short: 'deleted' },
  'intragenic-del': { color: '#ea580c', label: 'intragenic deletion', short: 'exons deleted' },
  'intragenic-dup': { color: '#15803d', label: 'intragenic duplication', short: 'exons duplicated' },
  'partial-copy': { color: '#65a30d', label: 'extra partial copy', short: 'partial copy' },
  duplicated: { color: '#16a34a', label: 'whole gene duplicated', short: 'duplicated' },
  domain: { color: '#b45309', label: 'chromatin domain changed', short: 'domain changed' },
  moved: { color: '#2563eb', label: 'intact, inverted or moved', short: 'inverted / moved' },
  nearby: { color: '#a16207', label: 'intact, near a breakpoint', short: 'near a breakpoint' },
  unaffected: { color: '#8796a8', label: 'unchanged', short: 'unchanged' },
};
