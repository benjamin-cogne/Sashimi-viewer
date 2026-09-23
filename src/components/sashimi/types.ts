/**
 * Shared Sashimi data types (transcripts, coverage, junctions, reads, variant sites, consensus groups).
 * Used by the web application, its API client and the standalone browser viewer.
 */
/** Displayed transcript model (RefSeq NM_/NR_ from UCSC, or Ensembl ENST as fallback). Coordinates are 1-based inclusive. */
export interface TranscriptData {
  gene_name: string; transcript_id: string; translation_id?: string | null; is_mane_select?: boolean;
  /** how the displayed model was chosen: MANE Select, RefSeq Select / Ensembl canonical, or the longest CDS / transcript */
  model_kind?: 'mane' | 'canonical' | 'longest' | 'chosen'; biotype?: string; source?: 'refseq' | 'ensembl'; chrom: string; strand: number;
  start: number; end: number; exons: { start: number; end: number; rank: number }[];
  /** Genomic CDS bounds (1-based inclusive); null/undefined for non-coding transcripts. */
  cds_start?: number | null; cds_end?: number | null;
}
/** One transcript model of a gene (RefSeq or Ensembl), 1-based inclusive coordinates. */
export interface TranscriptModel {
  id: string; name: string; source: 'refseq' | 'ensembl'; biotype: string;
  start: number; end: number; strand: number; exons: { start: number; end: number }[];
  cds_start?: number | null; cds_end?: number | null; is_mane: boolean; is_canonical: boolean;
}
export interface AllTranscripts { gene_name: string; chrom: string; strand: number; source: 'refseq' | 'ensembl'; transcripts: TranscriptModel[] }

/** Run-length encoded coverage: 0-based half-open [start, end) at constant depth. */
export interface CoverageRun { start: number; end: number; depth: number; }
/** Splice junction = intron interval, 0-based half-open [start, end). */
export interface JunctionArc { start: number; end: number; count: number; }
/**
 * A structural arc: its breakpoints (0-based half-open) and supporting reads. Arcs of one event whose breakpoints lie
 * within a few tens of bases (svmerge.ts) are merged: the event then says what it was made of.
 */
export interface SvArc extends JunctionArc {
  /** evidence units by source: cigar (CIGAR D), split (SA chain), clip (clip cluster placed by realignment), rescued, pair (discordant), '+/-' and '-/+' (inversion junctions) */
  sources?: Record<string, number>;
  /** breakpoints of the merged arcs: first and last start, first and last end */
  spread?: [number, number, number, number];
  /** arcs merged into this one */
  merged?: number;
}
/**
 * Reads that continue through an exon–intron boundary unspliced (intron retention / pre-mRNA):
 * one aligned block covering at least 6 bases on the exon side and 10 on the intron side.
 * Keyed by the boundary position: an intron start is the first intronic base (= exon end,
 * 0-based half-open), an intron end the first exonic base after the intron.
 */
export interface BoundarySpanning { intronStart: Record<number, number>; intronEnd: Record<number, number> }
/** Boundaries a coverage request wants spanning counts for (the displayed model's exons); junction ends are always included. */
export interface BoundaryHint { intronStarts: number[]; intronEnds: number[] }
/** What a library is: RNA-seq (spliced reads, junction arcs) or genomic DNA (exome, genome, long reads). */
export type LibraryType = 'rna' | 'dna' | 'unknown';
/** How a sample's library type was decided. */
export interface LibraryEvidence { type: LibraryType; source: 'header' | 'reads' | 'user' | 'none'; note: string }

/** A soft-clip cluster: reads clipped on the same side at the same position (a breakpoint candidate). */
export interface ClipCluster { pos: number; side: 'left' | 'right'; count: number; /** hard-clipped records without SA tag among the count (no sequence of their own) */ hard?: number }
/** A breakpoint used as a target for rescuing clipped reads: the ends of an arc and its kind. */
export interface Breakpoint { start: number; end: number; kind: 'deletion' | 'split' | 'duplication' | 'inversion' }
/**
 * Clipped reads whose clipped bases match the reference at the other end of a known breakpoint (second-pass style).
 * `own` breakpoints are this sample's arcs (the reads are added to the arc); borrowed ones come from another sample
 * shown next to it and are reported without an arc of their own.
 */
export interface RescuedClips extends Breakpoint { count: number; /** hard-clipped records without SA tag at the breakpoint end, attached by position only */ hard: number; own: boolean }
/** A clip cluster whose clipped consensus was placed on the reference of the window by realignment: it counts in the arc named here. */
export interface RealignedClip {
  pos: number; side: 'left' | 'right'; count: number; hard: number;
  /** where the clipped sequence was placed (0-based start) and on which strand, and the length matched */
  target: number; strand: '+' | '-'; matched: number;
  /** the arc it joined */
  arc: { start: number; end: number; kind: 'deletion' | 'split' | 'duplication' | 'inversion' };
}
/** Mates or split alignments on another chromosome, by position in the window (split: the breakpoint, rounded to 5 bp; pair: the 500 bp bin of the read start) and target chromosome. */
export interface ElsewhereLink { kind: 'split' | 'pair'; pos: number; chrom: string; count: number }
/**
 * Structural evidence of a genomic window (DNA libraries), 0-based half-open positions, counts scaled like
 * the junctions when the window was sampled.
 */
export interface StructuralEvidence {
  /** deletions of at least 50 bp inside reads (CIGAR D), by span */
  deletions: SvArc[];
  /** split reads whose next part continues further along the chromosome on the same strand (deletion-type), from the chain of alignments of each read (primary + SA parts ordered along the read), breakpoints rounded to 5 bp */
  splits: SvArc[];
  /** split reads whose next part goes back (tandem duplication-type) */
  duplications: SvArc[];
  /** split reads whose next part is on the other strand (inversion breakpoints) */
  inversions: SvArc[];
  /** unaligned stretch of the read between two adjacent parts on the reference: an insertion of about `len` bases */
  insertions: { pos: number; len: number; count: number }[];
  /** discordant pairs on the same chromosome (insert size far above the median, or mates on the same strand), both ends binned to 500 bp */
  discordant: JunctionArc[];
  /** split alignments and mates on other chromosomes */
  elsewhere: ElsewhereLink[];
  /** soft-clip clusters of at least 3 reads clipped by 20 bases or more, whose clipped sequence could not be placed in the window */
  clips: ClipCluster[];
  /** clip clusters placed by realignment of their clipped consensus: their reads count in the split-read arcs */
  realigned?: RealignedClip[];
  /** clipped reads rescued at a known breakpoint: their clipped bases match the reference at the other end of an arc */
  rescued?: RescuedClips[];
  /** median insert size of the proper pairs of the window (paired libraries) */
  insertMedian: number | null;
  /** reads the evidence was gathered from (after sampling, for a source that samples) */
  reads: number;
}

export interface SampleCoverage {
  sample_id: number; sample_name: string;
  coverage: CoverageRun[]; junctions: JunctionArc[];
  /** unspliced reads through the exon–intron boundaries (absent from sources that do not compute it) */
  spanning?: BoundarySpanning;
  /** the window actually read (0-based half-open) when the source shrank the requested margins to stay within its read budget */
  window?: { start: number; end: number };
  /** set when only every `rate`-th read was decoded: depths, junction and boundary counts are scaled back by `rate` (estimates) */
  sampled?: { rate: number; total: number; decoded: number };
  /**
   * Reads counted and the fraction of them carrying a splice gap (CIGAR N): library-type evidence. A source that keeps its
   * counts between requests gives it over everything it has counted for the sample on this chromosome, not just the window.
   */
  spliced?: { reads: number; fraction: number };
  /** structural evidence, when the caller asked for it (DNA samples) */
  structural?: StructuralEvidence;
  error?: string;
}

/** One alignment for the reads track (compact keys, all coordinates 0-based half-open).
 *  n name · s/e reference span · r reverse (1/0) · q MAPQ · f SAM flag · nh NH tag ·
 *  b aligned blocks · d deletions · i insertions [pos, len] · m mismatches [pos, base, qual] ·
 *  c soft clips [left, right]. */
export interface AlignedRead {
  n: string; s: number; e: number; r: 0 | 1; q: number; f: number; nh: number | null;
  b: [number, number][]; d: [number, number][]; i: [number, number][]; m: [number, string, number][]; c: [number, number];
  /** mate: 0-based start (`mp`), chromosome when not the read's own (`mc`), template length as the aligner set it (`tl`); absent when unpaired or unknown */
  mp?: number; mc?: string; tl?: number;
  /** soft-clipped bases at the left and right ends of the alignment (read orientation as stored, i.e. the reference strand); absent when the sequence was not available */
  cs?: [string, string];
  /** hard-clipped lengths at the left and right ends (bases the record does not carry: they sit in the read's primary record) */
  h?: [number, number];
  /** inserted bases, one string per entry of `i`; absent when the sequence was not available */
  is?: string[];
  /** SA tag of a split read: "rname,pos,strand,CIGAR,mapQ,NM;" per other part */
  sa?: string;
  /**
   * haplotag written by a phasing tool (WhatsHap or LongPhase haplotag, PacBio HiPhase, DRAGEN): HP the haplotype (1, 2;
   * DRAGEN's "copy label" may go higher), PS the phase set it belongs to, PC the Phred-scaled confidence; absent when untagged
   */
  hp?: number; ps?: number; pc?: number;
  /** CpG calls of a long read with base-modification tags, when asked for: position of the CpG's C, P(5mC) × 255, flat */
  me?: number[];
}
/** A variable site called from the reads: SNV, insertion or deletion above the support and fraction thresholds. */
export interface VariantSite {
  pos: number; kind: 'snv' | 'ins' | 'del'; ref: string; alt: string; length: number;
  alt_count: number; depth: number; vaf: number;
}
/** One consensus group: a local haplotype × splice pattern with its supporting read count. */
export interface ReadGroup {
  id: string; kind: 'consensus' | 'ambiguous' | 'minor'; n: number; frac: number;
  chain: [number, number][]; alleles: string[]; blocks: [number, number][]; dense: [number, number][];
  absorbed: number; compatible?: string[]; patterns?: number;
}
/** One haplotype block of the read-based phasing: the heterozygous sites the reads tie together, and the two haplotypes over them. */
export interface PhaseBlock {
  id: string;
  /** genomic span of the block's sites, 0-based half-open */
  start: number; end: number;
  /** indices into PhaseResult.sites, in position order */
  sites: number[];
  /** allele of each site on haplotype 1 / 2 */
  h1: ('ref' | 'alt')[]; h2: ('ref' | 'alt')[];
  /** fragments (read + mate) assigned to haplotype 1 / 2 */
  support: [number, number];
  /** fragments matching both haplotypes equally */
  ambiguous: number;
  /** fragments assigned to a haplotype but disagreeing with it at one site or more (errors, mosaic alleles, a third haplotype) */
  conflicting: number;
  /** fragments linking two sites of the block: seeing the same phase (ref–ref or alt–alt) or the opposite */
  links: { a: number; b: number; same: number; diff: number }[];
  /** why the previous block ended before this one */
  breakBefore?: 'no link' | 'conflict';
}
export interface UnphasedSite { site: number; reason: 'low' | 'unlinked' | 'conflict' }
export interface PhaseResult {
  sites: VariantSite[];
  /** indices of the homozygous sites (on both haplotypes) */
  hom: number[];
  blocks: PhaseBlock[];
  unphased: UnphasedSite[];
  /** fragments (reads with their mates) examined */
  fragments: number;
  /** heterozygous sites considered */
  het: number;
}
/** One haplotype's consensus over a phase set: where at least HAP_MIN_DEPTH of its reads cover, and what most of them carry. */
export interface HaplotypeConsensus {
  /** 1, 2 (HP tag, or the block's haplotype); higher for DRAGEN copy labels */
  hap: number;
  /** reads (tags) or fragments (in-page phasing) of this haplotype */
  reads: number;
  /** stretches the haplotype's reads cover, deletions included, 0-based half-open */
  covered: [number, number][];
  /** variants carried by at least half of the haplotype's reads over them */
  sites: VariantSite[];
  /** median PC (assignment confidence) of its tagged reads, when the file gives one */
  pc?: number;
}
/** A phase set: haplotypes linked within it, not across (H1 here is unrelated to H1 of the next set). */
export interface HaplotypeSet { id: string; ps: number | null; start: number; end: number; haps: HaplotypeConsensus[] }
/** The haplotypes of a reads window, from the file's haplotags or from the in-page read-based phasing. */
export interface HaplotypeView {
  source: 'tags' | 'reads';
  sets: HaplotypeSet[];
  /** reads (tags) or fragments (reads) given a haplotype, and not */
  assigned: number; unassigned: number;
  /** heterozygous sites of the window that the two haplotypes of their set do not split (both, or neither, carry the alt) */
  notSplit: { pos: number; kind: VariantSite['kind']; alt: string; fractions: number[]; set: string }[];
  /** heterozygous sites checked, and the reads whose allele at one of them contradicts their haplotype */
  checked: number; conflicting: number;
}
export interface ReadsResponse {
  sample_id: number; sample_name: string;
  reads: AlignedRead[]; total: number; shown: number;
  sites: VariantSite[]; groups: ReadGroup[];
  /** read-based phasing of the window (collapsed mode with two haplotypes) */
  phase?: PhaseResult;
  /** the two haplotypes of the window as consensus rows (collapsed mode with two haplotypes) */
  haplotypes?: HaplotypeView;
  /** reads of the window carrying a haplotag (HP), and the phase sets (PS) among them: the file is phased */
  haplotags?: { tagged: number; sets: number };
  /** the reads are long (median aligned length above 1 kb): noise filters apply */
  long_reads?: boolean;
  /** Reference sequence covering the window (plus margin), or null when no source is available. */
  reference: { start: number; seq: string } | null;
  reference_source: 'fasta' | 'ensembl' | 'browser' | null;
}

/** A neighbouring gene with its canonical transcript model (1-based inclusive coordinates). */
export interface GeneModel {
  gene_id: string; gene_name: string; biotype: string; strand: number; start: number; end: number;
  transcript_id: string; is_canonical: boolean; exons: { start: number; end: number }[];
  cds_start?: number | null; cds_end?: number | null;
}

/** Depth statistics of one exon in one sample (aligned bases only; reads with a block on the exon). */
export interface ExonDepth { median: number; mean: number; reads: number }

export type Strandness = 'firststrand' | 'secondstrand' | 'unstranded' | 'unknown';

export interface SampleExonDepths {
  sample_id: number;
  sample_name: string;
  strandness: Strandness;
  strand_fraction: number | null;
  /** Same order as the request; empty when `error` is set. */
  exons: ExonDepth[];
  error?: string;
}

export interface ExonUsageResponse {
  run_id: number;
  chrom: string;
  /** 0-based half-open exon intervals as evaluated */
  exons: [number, number][];
  samples: SampleExonDepths[];
}

/** Protein feature (UniProt / Pfam domain via UCSC, or Ensembl protein_feature), 1-based inclusive amino-acid coordinates. */
export interface ProteinDomain { type: string; id: string; description: string; interpro?: string | null; start: number; end: number }

/** A gene span (1-based inclusive) already known to the caller, e.g. from the outlier tables: saves the symbol lookup. */
export interface RegionHint { chrom: string; start: number; end: number }

/** The coding model whose protein domains are wanted: exons and CDS bounds 0-based half-open. */
export interface ProteinModelRef {
  transcriptId: string; translationId?: string | null; chrom: string; strand: number;
  exons: { start: number; end: number }[]; cdsStart: number; cdsEnd: number;
}

/** A variant previously identified in a sample (clinical indication, diagnostic conclusion or chromosome map). */
export type KnownVariantKind = 'snv' | 'indel' | 'del' | 'dup' | 'inv' | 'ins' | 'cnv' | 'bnd' | 'other';
export interface KnownVariant {
  id: string; kind: KnownVariantKind;
  /** `chr`-prefixed; empty when the notation carried no chromosome (a bare `g.` with the gene name only) */
  chrom: string;
  /** 0-based half-open; an SNV spans one base */
  start: number; end: number;
  /** short text drawn next to the marker (c. notation, else the genomic notation) */
  label: string;
  gene?: string; cdna?: string; protein?: string;
  /** the notation as typed */
  text: string;
  source: 'indication' | 'diagnostic' | 'chromosome_map';
  /** build inferred from the notation (NC_ accession version, ISCN bracket); null when unknown */
  build: 'GRCh38' | 'GRCh37' | null;
}

/** A common variant of the SNP track (0-based half-open), with one allele frequency per frequency project. */
export interface CommonSnp {
  id: string; start: number; end: number; ref: string; alts: string[];
  cls: 'snv' | 'ins' | 'del' | 'delins' | 'mnv' | 'other';
  /** highest minor-allele frequency over the projects */
  maxAf: number;
  afs: { source: string; af: number }[];
  source: 'dbSNP155' | 'ensembl';
  impact?: string;
}

/** A GTEx tissue (tissueSiteDetail) with its portal colour and RNA-seq sample count. */
export interface GtexTissue { id: string; name: string; site: string; color: string; samples: number }
/** Median junction and exon read counts of a gene in one tissue, with the reads-per-base profile derived from the exons. */
export interface GtexProfile {
  dataset: 'gtex_v10' | 'gtex_v8';
  gencodeId: string;
  tissue: GtexTissue;
  /** 0-based half-open introns, count = median read count over the tissue's samples */
  junctions: JunctionArc[];
  /** collapsed gene-model exons, 0-based half-open, median read count */
  exons: { start: number; end: number; median: number }[];
  /** median reads per base over the exons (0 in introns) */
  coverage: CoverageRun[];
  unit: string;
  /** set when one of the two expression calls failed (the profile is partial) */
  warning?: string;
  /** median gene expression in the tissue (TPM), null when unavailable */
  tpm: number | null;
}
