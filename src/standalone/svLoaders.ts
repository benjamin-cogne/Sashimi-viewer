/**
 * What the SV browser reads, for one genome build: the RefSeq models of a window, ClinGen dosage and gnomAD
 * constraint, the cytogenetic bands, the DGV Gold Standard, and where a gene lies. Straight from the UCSC and gnomAD
 * APIs, so the SV browser can switch build whatever the session's reference.
 */
import type { GenomeBuild } from './ensembl';
import { chromLength } from './assembly';
import { getRegionTranscripts, getDosageScores, getCytobands, getDgvGold, getTranscript } from './ucsc';
import { getConstraint } from './gnomad';
import { parseDgvItem, type PopulationSv } from '../components/sashimi/svPopulation';
import type { SvGeneLoaders } from '../components/sashimi/SvGenomeBrowser';
import type { Constraint, DosageScore } from '../components/sashimi/svGenes';

/** ClinGen dosage of a window and gnomAD constraint of genes; each part fails on its own (its map left empty). */
export async function geneAnnotations(build: GenomeBuild, chrom: string, start: number, end: number, symbols: string[]): Promise<{ dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[] }> {
  const notes: string[] = [];
  const [dosage, constraint] = await Promise.all([
    getDosageScores(build, chrom, start, end).then(d => { notes.push('ClinGen dosage: UCSC Genome Browser'); return d; }).catch(e => { notes.push(`ClinGen dosage unavailable (${e?.message || e})`); return new Map<string, DosageScore>(); }),
    getConstraint(build, symbols, { chrom, start, end }).then(c => { notes.push(c.source); return c.values as Map<string, Constraint>; }).catch(e => { notes.push(`gnomAD constraint unavailable (${e?.message || e})`); return new Map<string, Constraint>(); }),
  ]);
  return { dosage, constraint, notes };
}

export function svLoadersFor(build: GenomeBuild): Omit<SvGeneLoaders, 'fallback'> {
  return {
    getTranscripts: (c, a, b) => getRegionTranscripts(build, c, a, b),
    getAnnotations: (c, a, b, sy) => geneAnnotations(build, c, a, b, sy),
    getCytobands: c => getCytobands(build, c, chromLength(build, c) ?? 300_000_000),
    getPopulationSvs: async (c, a, b) => (await getDgvGold(build, c, a, b)).map(parseDgvItem).filter((v): v is PopulationSv => !!v),
    locateGene: async symbol => {
      const t = await getTranscript(build, symbol);
      // 1-based inclusive → 0-based half-open
      return { chrom: t.chrom.startsWith('chr') ? t.chrom : `chr${t.chrom}`, start: t.start - 1, end: t.end, name: t.gene_name };
    },
  };
}
