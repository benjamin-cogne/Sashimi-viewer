/**
 * Data source used by the Sashimi viewer. The web application implements it with the
 * FastAPI backend (see apiDataSource.ts); the standalone HTML viewer implements it in the
 * browser on top of local BAM/CRAM files and the UCSC / Ensembl REST APIs (src/standalone).
 */
import type { TranscriptData, SampleCoverage, BoundaryHint, ReadsResponse, AllTranscripts, GeneModel, ExonUsageResponse, ProteinDomain, CommonSnp, GtexTissue, GtexProfile, RegionHint, ProteinModelRef, KnownVariant, LibraryEvidence } from './types';

export interface SampleRef { id: number; name: string }

/** Noise handling of a reads request. */
export interface ReadsOptions {
  /** for long reads (median aligned length above 1 kb): indels shorter than this are neither called nor collapsed on */
  longReadMinIndel?: number;
  /** for long reads: floor of the alternate-allele fraction a site needs (their error rate makes the short-read threshold too low) */
  longReadMinVaf?: number;
}
/** Budget of a coverage request. */
export interface CoverageOptions {
  /** the part of the window that must be read whatever its depth (the visible view); the rest is margin the source may shrink */
  core?: { start: number; end: number };
  /** reads decoded at most; a deeper window is sampled systematically (1 read in 2, 4, 8…) and its counts scaled back */
  maxReads?: number;
  /** also collect the structural evidence of a genomic library (deletions, split reads, soft clips, discordant pairs) */
  structural?: boolean;
}

export interface SashimiDataSource {
  /** Displayed model of a gene (MANE Select / RefSeq Select). `hint` is the gene span when already known (skips the symbol lookup). */
  getTranscript(geneName: string, geneId?: string, hint?: RegionHint): Promise<TranscriptData>;
  getAllTranscripts(geneName: string, geneId?: string, hint?: RegionHint): Promise<AllTranscripts>;
  /** Coverage runs and junctions for a 0-based half-open window; `boundaries` asks for unspliced reads through those exon–intron boundaries too (junction ends are always counted). */
  getCoverage(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint, opts?: CoverageOptions): Promise<SampleCoverage>;
  getReads(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, maxReads: number,
    mode: 'reads' | 'collapsed', minSupport: number, minVaf: number, opts?: ReadsOptions): Promise<ReadsResponse>;
  /** All samples that can be added as tracks. */
  getRunSamples(runId: number): Promise<SampleRef[]>;
  /** One comparison sample to load next to the primary one; may reject when none exists. */
  getRandomSample(runId: number, excludeId: number): Promise<SampleRef>;
  /** Genes overlapping a 1-based inclusive region with their canonical transcript, `exclude` (name or id) left out. */
  getRegionGenes(chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]>;
  /**
   * Median depth and read count of each exon (0-based half-open) for every sample of the run,
   * the cohort that normalises exon usage. `strand` is the gene strand (sense reads only when the
   * library is detected as stranded).
   */
  getExonUsage(runId: number, chrom: string, strand: number, exons: [number, number][], uniqueOnly: boolean): Promise<ExonUsageResponse>;
  /** Reference sequence (+ strand, upper case) of a 0-based half-open window; null when no reference is available. */
  getReference(chrom: string, start: number, end: number): Promise<string | null>;
  /** Protein domains (UniProt / Pfam) of a coding model, amino-acid coordinates. */
  getProteinDomains(model: ProteinModelRef): Promise<ProteinDomain[]>;
  /** Library type of a sample from what the source knows before any read is decoded (the aligner named in the header); absent when it cannot tell. */
  getLibraryType?(sampleId: number): Promise<LibraryEvidence>;
  /** Variants previously identified in a sample (clinical indication, diagnostic, chromosome map); absent when the source has no such record. */
  getKnownVariants?(sampleId: number): Promise<KnownVariant[]>;
  /** Common variants (dbSNP 155 common via UCSC, Ensembl fallback) overlapping a 0-based half-open window. */
  getCommonSnps(chrom: string, start: number, end: number): Promise<CommonSnp[]>;
  /** GTEx tissues available for tissue tracks (rejects when GTEx cannot be used, e.g. GRCh37). */
  getGtexTissues(): Promise<GtexTissue[]>;
  /** Median junction / exon profile of the gene in one tissue. */
  getGtexProfile(geneName: string, geneId: string | undefined, tissue: GtexTissue): Promise<GtexProfile>;
}
