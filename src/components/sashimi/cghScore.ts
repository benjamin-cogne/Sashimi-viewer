/**
 * A priority for the review of an array-CGH call (HIGH, MODERATE, LOW, none): which calls to look at first, from what
 * its genes are known for. Not a classification: the class stays the reader's.
 *
 * Each rule says what it looks for; a call takes the highest level of the rules it meets, and lists them. The rules
 * follow the evidence the ACMG / ClinGen technical standards for CNVs weigh (Riggs ER et al., Genet Med 2020;22:245-57):
 * established dosage-sensitive genes (section 2: ClinGen haploinsufficiency / triplosensitivity scores, 3 sufficient,
 * 2 some evidence), predictors of haploinsufficiency (2H: gnomAD LoF constraint), the number of protein-coding genes
 * (section 3: losses 25-34 genes 0.45 point, ≥ 35 0.90; gains 35-49 0.45, ≥ 50 0.90), and a breakpoint inside a gene
 * (2I-2L for gains: a partial duplication may disrupt it). Genes of PanelApp rated green (diagnostic grade) with a
 * mode of inheritance that a loss or both copies lost can explain, and the recurrent CNVs the AChro-Puce network
 * classifies (class 5 / 4 pathogenic or likely pathogenic, PIEV susceptibility), are added as the laboratory uses them.
 *
 * The rules are the array below: change a level, a threshold or a test there.
 */
import type { MoiClass } from '../../standalone/panelapp';
import type { AchroClass } from './achropuce';
import type { CghSex } from './cghIntervals';

export type Level = 'high' | 'moderate' | 'low' | 'none';
export const LEVEL_RANK: Record<Level, number> = { high: 3, moderate: 2, low: 1, none: 0 };
export const LEVEL_STYLE: Record<Level, { label: string; color: string; bg: string }> = {
  high: { label: 'HIGH', color: '#b91c1c', bg: '#fee2e2' },
  moderate: { label: 'MODERATE', color: '#c2410c', bg: '#ffedd5' },
  low: { label: 'LOW', color: '#475569', bg: '#eceff4' },
  none: { label: 'none', color: '#8796a8', bg: 'transparent' },
};

/** a gene of the call: in it whole, or cut by a breakpoint; whether an exon of one of its models is in the call */
export interface CnvGene { symbol: string; whole: boolean; exons: boolean; start: number; end: number; tx: string }
/** PanelApp's entries of a gene: its best rating (0: in no panel), the modes of its green and amber entries */
export interface PanelView { best: 0 | 1 | 2 | 3; green: MoiClass[]; amber: MoiClass[]; diseases: string[] }
export interface GeneView {
  g: CnvGene;
  hi: number | null; ts: number | null;
  loeuf: number | null; pli: number | null; version: string | null;
  /** null while PanelApp is asked, or when it failed (`paState` says which) */
  pa: PanelView | null;
  paState?: 'asked' | 'failed' | 'not asked' | 'done';
}
export interface AchroHit { id: string; name: string; type: 'del' | 'dup'; cls: AchroClass[]; how: string }
export interface ScoreInput {
  kind: 'loss' | 'gain'; copies: number | null; chrom: string;
  /** the patient's sex, null when not given */
  sex: CghSex | null;
  /** the call's size, kb (the file's Size(kb)) */
  sizeKb: number;
  genes: GeneView[];
  /** whether the genes are RefSeq's (counted) or only the file's names (not counted) */
  refseq: boolean;
  achro: AchroHit[];
}
export interface Reason { id: string; level: Level; text: string }
export interface Score { level: Level; reasons: Reason[] }

/** a call larger than this (kb) is high priority whatever its genes (the AChro-Puce network's rule) */
export const BIG_KB = 1000;
const MONO: MoiClass[] = ['AD', 'AD+AR', 'XLD'];
const BI: MoiClass[] = ['AR', 'AD+AR'];
const any = (xs: MoiClass[] | undefined, of: MoiClass[]) => !!xs?.some(x => of.includes(x));
/** LoF-intolerant: LOEUF below the threshold of its gnomAD version (v2.1.1 0.35, v4 0.45), or pLI ≥ 0.9 */
const intolerant = (v: GeneView) => (v.loeuf != null && v.loeuf < (v.version === 'v4' ? 0.45 : 0.35)) || (v.pli != null && v.pli >= 0.9);

interface Rule {
  id: string; level: Level; kind: 'loss' | 'gain';
  /** what the rule looks for, as listed in the page */
  label: string;
  /** the genes that meet it */
  gene?: (v: GeneView, i: ScoreInput) => boolean;
  /** or the call as a whole: a text when it meets it */
  call?: (i: ScoreInput) => string | null;
}

export const SCORE_RULES: Rule[] = [
  // ---- losses ----
  { id: 'L1', level: 'high', kind: 'loss', label: 'exons of a gene with a ClinGen haploinsufficiency score of 3 (sufficient evidence)', gene: v => v.g.exons && v.hi === 3 },
  { id: 'L2', level: 'high', kind: 'loss', label: 'exons of a PanelApp green gene with a monoallelic mode (AD, AD/AR, X-linked dominant)', gene: v => v.g.exons && any(v.pa?.green, MONO) },
  { id: 'L3', level: 'high', kind: 'loss', label: 'both copies lost (x0) of exons of a PanelApp green gene with a biallelic mode (AR)', gene: (v, i) => i.copies === 0 && v.g.exons && any(v.pa?.green, BI) },
  { id: 'L3X', level: 'high', kind: 'loss', label: 'exons of a PanelApp green X-linked recessive gene lost in a male (hemizygous), or both copies lost (x0)', gene: (v, i) => i.chrom === 'chrX' && (i.sex === 'male' || i.copies === 0) && v.g.exons && any(v.pa?.green, ['XLR']) },
  { id: 'L4', level: 'high', kind: 'loss', label: 'an AChro-Puce recurrent deletion of class 5 or 4', call: i => hits(i, 'del', ['5', '4']) },
  { id: 'L5', level: 'high', kind: 'loss', label: `35 protein-coding genes or more (ACMG/ClinGen section 3C: 0.90 point), and / or more than ${BIG_KB / 1000} Mb (AChro-Puce)`, call: i => big(i, 35) },
  { id: 'L6', level: 'moderate', kind: 'loss', label: 'exons of a gene with a ClinGen haploinsufficiency score of 2 (some evidence)', gene: v => v.g.exons && v.hi === 2 },
  { id: 'L7', level: 'moderate', kind: 'loss', label: 'exons of a LoF-intolerant gene (gnomAD LOEUF < 0.35 in v2.1.1, < 0.45 in v4, or pLI ≥ 0.9)', gene: v => v.g.exons && intolerant(v) },
  { id: 'L8', level: 'moderate', kind: 'loss', label: 'one copy of a PanelApp green gene with a biallelic mode (AR) lost: look for a variant on the other copy', gene: (v, i) => i.copies !== 0 && v.g.exons && any(v.pa?.green, BI) && !any(v.pa?.green, MONO) },
  { id: 'L9', level: 'moderate', kind: 'loss', label: 'one copy of a PanelApp green X-linked recessive gene lost in a female (carrier), or the sex not given', gene: (v, i) => i.chrom === 'chrX' && i.sex !== 'male' && i.copies !== 0 && v.g.exons && any(v.pa?.green, ['XLR']) },
  { id: 'L10', level: 'moderate', kind: 'loss', label: 'a PanelApp green gene whose mode is not given, or an amber gene with a monoallelic mode', gene: v => v.g.exons && (any(v.pa?.green, ['?']) || any(v.pa?.amber, MONO)) },
  { id: 'L11', level: 'moderate', kind: 'loss', label: 'an AChro-Puce recurrent deletion classed PIEV (susceptibility)', call: i => hits(i, 'del', ['PIEV']) },
  { id: 'L12', level: 'moderate', kind: 'loss', label: '25 to 34 protein-coding genes (section 3B: 0.45 point)', call: i => (i.refseq && i.genes.length >= 25 && i.genes.length < 35 ? `${i.genes.length} protein-coding genes` : null) },
  { id: 'L13', level: 'low', kind: 'loss', label: 'exons of protein-coding genes lost, none of the above', gene: v => v.g.exons },
  { id: 'L14', level: 'low', kind: 'loss', label: 'genes cut without an exon lost (introns, ends of genes)', gene: v => !v.g.exons },
  { id: 'L15', level: 'low', kind: 'loss', label: 'an AChro-Puce recurrent deletion of class 3 (uncertain significance)', call: i => hits(i, 'del', ['3']) },
  { id: 'L16', level: 'low', kind: 'loss', label: 'an AChro-Puce recurrent deletion of class 2 or 1 (likely benign, benign)', call: i => hits(i, 'del', ['2', '1']) },
  // ---- gains ----
  { id: 'G1', level: 'high', kind: 'gain', label: 'a whole gene with a ClinGen triplosensitivity score of 3 (sufficient evidence)', gene: v => v.g.whole && v.ts === 3 },
  { id: 'G2', level: 'high', kind: 'gain', label: 'an AChro-Puce recurrent duplication of class 5 or 4', call: i => hits(i, 'dup', ['5', '4']) },
  { id: 'G3', level: 'high', kind: 'gain', label: `50 protein-coding genes or more (section 3C for gains: 0.90 point), and / or more than ${BIG_KB / 1000} Mb (AChro-Puce)`, call: i => big(i, 50) },
  { id: 'G4', level: 'moderate', kind: 'gain', label: 'a whole gene with a ClinGen triplosensitivity score of 2 (some evidence)', gene: v => v.g.whole && v.ts === 2 },
  { id: 'G5', level: 'moderate', kind: 'gain', label: 'a breakpoint inside a gene with a haploinsufficiency score of 3 or a PanelApp green monoallelic gene: a partial duplication may disrupt it (2I-2L)', gene: v => !v.g.whole && v.g.exons && (v.hi === 3 || any(v.pa?.green, MONO)) },
  { id: 'G6', level: 'moderate', kind: 'gain', label: 'an AChro-Puce recurrent duplication classed PIEV (susceptibility)', call: i => hits(i, 'dup', ['PIEV']) },
  { id: 'G7', level: 'moderate', kind: 'gain', label: '35 to 49 protein-coding genes (section 3B for gains: 0.45 point)', call: i => (i.refseq && i.genes.length >= 35 && i.genes.length < 50 ? `${i.genes.length} protein-coding genes` : null) },
  { id: 'G8', level: 'low', kind: 'gain', label: 'protein-coding genes gained, none of the above', gene: () => true },
  { id: 'G9', level: 'low', kind: 'gain', label: 'an AChro-Puce recurrent duplication of class 3 (uncertain significance)', call: i => hits(i, 'dup', ['3']) },
  { id: 'G10', level: 'low', kind: 'gain', label: 'an AChro-Puce recurrent duplication of class 2 or 1 (likely benign, benign)', call: i => hits(i, 'dup', ['2', '1']) },
];

/** a large call: more than BIG_KB, and / or this many protein-coding genes (counted only when they are RefSeq's) */
function big(i: ScoreInput, genes: number): string | null {
  const parts = [...(i.refseq && i.genes.length >= genes ? [`${i.genes.length} protein-coding genes`] : []), ...(i.sizeKb > BIG_KB ? [`${(i.sizeKb / 1000).toFixed(3)} Mb`] : [])];
  return parts.length ? parts.join(', ') : null;
}
function hits(i: ScoreInput, type: 'del' | 'dup', cls: AchroClass[]): string | null {
  const h = i.achro.filter(a => a.type === type && a.cls.some(c => cls.includes(c)));
  return h.length ? h.map(a => `${a.name} (${a.cls.join(' / ')}; ${a.how})`).join('; ') : null;
}

/** The call's priority: the highest level among the rules it meets, each with the genes or the region that meet it. */
export function scoreCnv(i: ScoreInput): Score {
  const reasons: Reason[] = [];
  for (const r of SCORE_RULES) {
    if (r.kind !== i.kind) continue;
    if (r.call) { const t = r.call(i); if (t) reasons.push({ id: r.id, level: r.level, text: `${r.label}: ${t}` }); continue; }
    const gs = i.genes.filter(v => r.gene!(v, i)).map(v => v.g.symbol);
    if (gs.length) reasons.push({ id: r.id, level: r.level, text: `${r.label}: ${gs.slice(0, 8).join(', ')}${gs.length > 8 ? ` +${gs.length - 8}` : ''}` });
  }
  // a call meeting no HIGH or MODERATE rule is LOW, whatever it holds
  if (!reasons.length) reasons.push({ id: i.kind === 'loss' ? 'L0' : 'G0', level: 'low', text: i.genes.length ? 'no HIGH or MODERATE rule met' : 'no protein-coding gene in the call' });
  // a low reason only says what is left: dropped when a higher one is there
  const top = reasons.reduce<Level>((m, r) => (LEVEL_RANK[r.level] > LEVEL_RANK[m] ? r.level : m), 'low');
  return { level: top, reasons: top === 'low' ? reasons : reasons.filter(r => r.level !== 'low') };
}
