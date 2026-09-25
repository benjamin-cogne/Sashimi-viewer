/**
 * Alignment file kinds beyond BAM and CRAM.
 *
 * A module under src/plugins/<name>/index.ts registers a kind at load time (src/plugins/index.ts loads every one it
 * finds, in the page and in the variant worker alike). The shell then accepts the kind's extensions, pairs the files
 * with their index (or takes a file that indexes itself alone), and LocalDataSource hands every request for such a
 * sample to the kind's provider instead of reading it as a BAM or CRAM. With no plugin present nothing changes.
 *
 * A provider answers in the viewer's own terms (SampleCoverage, ReadsResponse, VariantScan, MethylWindow), so a
 * registered kind gets the tracks, collapse, phasing, variants and methylation overlays of a BAM. readsWindow.ts has
 * the shared tail of a reads answer, methylation.ts `addCall` the counting of calls held in another form.
 */
import type { AlignedRead, BoundaryHint, Breakpoint, LibraryEvidence, ReadsResponse, RescuedClips, SampleCoverage } from '../components/sashimi/types';
import type { CoverageOptions, ReadsOptions, VariantScan, VariantScanOptions } from '../components/sashimi/datasource';
import type { LocalSample, ReferenceChoice } from './localSource';
import type { RawRead } from './alignments';
import type { MethylWindow } from './methylation';

/** 'bam' and 'cram' are read by LocalDataSource itself; any other kind names a registered FileKind. */
export type SampleKind = 'bam' | 'cram' | (string & {});

/** What LocalDataSource lends a provider. */
export interface ProviderHost {
  /** reference bases of [start, end), 0-based, upper case; null when no source has them */
  getReferenceSeq(chrom: string, start: number, end: number): Promise<string | null>;
  /** which source answered the last getReferenceSeq */
  readonly referenceSource: 'fasta' | 'ensembl' | null;
  readonly reference: ReferenceChoice;
}

/**
 * Reads one sample of a registered kind. Every method gets the sample as it is now (a rename replaces it). Only
 * coverage and reads are required; without the others the overlay says the file cannot give it, and exon usage
 * leaves the sample out.
 */
export interface SampleProvider {
  getCoverage(s: LocalSample, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint, opts?: CoverageOptions): Promise<SampleCoverage>;
  getReads(s: LocalSample, chrom: string, start: number, end: number, uniqueOnly: boolean, maxReads: number,
    mode: 'reads' | 'collapsed', minSupport: number, minVaf: number, opts?: ReadsOptions): Promise<ReadsResponse>;
  /** every variant site of the window (LocalDataSource.countVariants: in the worker when there is one) */
  countVariants?(s: LocalSample, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan>;
  /** CpG methylation of the window (LocalDataSource.countMethylation) */
  countMethylation?(s: LocalSample, chrom: string, start: number, end: number, opts?: { signal?: AbortSignal; onProgress?: (fraction: number) => void }): Promise<MethylWindow>;
  getLibraryType?(s: LocalSample): Promise<LibraryEvidence>;
  /** reads over one exon for exon usage and strandness, one in `rate` kept (as LocalDataSource scans them) */
  exonReads?(s: LocalSample, chrom: string, start: number, end: number, uniqueOnly: boolean, cap: number): Promise<{ total: number; rate: number; kept: RawRead[] }>;
  rescueClips?(s: LocalSample, chrom: string, breakpoints: Breakpoint[], uniqueOnly: boolean): Promise<RescuedClips[]>;
  getPrimaryRecord?(s: LocalSample, chrom: string, start: number, name: string): Promise<{ seq: string; flags: number; cigar: string } | null>;
  setReference?(reference: ReferenceChoice): void;
  /** drop what was kept for an overlay switched off, or decoded records */
  release?(what: 'methylation' | 'variants' | 'records'): void;
  /** the sample was removed or its files replaced */
  dispose?(): void;
}

export interface FileKind {
  /** the sample kind it registers ('bam' and 'cram' are taken) */
  id: string;
  /** shown in the sample list, e.g. "BAM" */
  label: string;
  /** file name endings, lower case, with the dot */
  extensions: string[];
  /** endings of the index file, appended to the data file's name; absent when the file indexes itself */
  index?: string[];
  /** a provider for one sample; nothing may be read before its first request */
  open(sample: LocalSample, host: ProviderHost): SampleProvider;
}

const kinds = new Map<string, FileKind>();

export function registerFileKind(k: FileKind): void {
  if (k.id === 'bam' || k.id === 'cram') throw new Error(`file kind "${k.id}" is built in`);
  kinds.set(k.id, k);
}
export const fileKind = (id: string): FileKind | undefined => kinds.get(id);
/** the registered kind a file name ends with */
export function fileKindOf(name: string): FileKind | undefined {
  const n = name.toLowerCase();
  for (const k of kinds.values()) if (k.extensions.some(x => n.endsWith(x))) return k;
  return undefined;
}
/** every ending a registered kind reads, data and index files */
export const kindExtensions = (): string[] => [...kinds.values()].flatMap(k => [...k.extensions, ...(k.index ?? [])]);
/** label of a sample kind: the registered one, else the kind in capitals */
export const kindLabel = (kind: string): string => kinds.get(kind)?.label ?? kind.toUpperCase();
