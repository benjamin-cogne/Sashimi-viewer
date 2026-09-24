/**
 * SashimiDataSource over local BAM/CRAM files, decoded in the browser with the GMOD
 * libraries. Nothing from the alignment files leaves the machine; only gene lookups and
 * reference-sequence requests go to the UCSC API (Ensembl REST as fallback), unless a local
 * FASTA is given.
 */
import { SharedBudget } from '@gmod/shared-read-cache';
import { BamFile } from '@gmod/bam';
import { IndexedCramFile, CraiIndex } from '@gmod/cram';
import { IndexedFasta, BgzipIndexedFasta } from '@gmod/indexedfasta';
import { BlobFile } from 'generic-filehandle2';
import { unzip } from '@gmod/bgzf-filehandle';
import type { AlignedRead, AllTranscripts, BoundarySpanning, ExonUsageResponse, GeneModel, GtexProfile, GtexTissue, KnownVariant, LibraryEvidence, ProteinDomain, ProteinModelRef, ReadsResponse, RegionHint, SampleCoverage, BoundaryHint, SampleExonDepths, TranscriptData, VariantSite } from '../components/sashimi/types';
import type { CoverageOptions, ReadsOptions, SashimiDataSource, SampleRef, VariantScan, VariantScanOptions } from '../components/sashimi/datasource';
import type { Breakpoint, RescuedClips } from '../components/sashimi/types';
import { RESCUE_MIN_CLIP, RESCUE_TOLERANCE_BP, SV_MIN_DELETION, clipEnds, cramCigar, cramMismatches, detectStrandness, encodeRead, exonDepth, hasRealignableClips, hasRescuableClips, keepFlags, rescueClipEnds, strandKeeper, structuralEvidence, uniqueFrom, type RawRead, type StrandnessCall } from './alignments';
import { CoverageState, Layer, packCigar, readSlice, type CoverageSlice } from './coverage';
import { AlleleLayer, AlleleState, refWindow, sitesFromCounts, type RefWindow } from './alleles';
import { MethylCounts, MethylState, countRead, cpgSites, methylWindow, newModScratch, visitReadCalls, type MethylWindow } from './methylation';
import { callSites, collapseReads } from './collapse';
import { phaseReads } from './phasing';
import { haplotagCounts, windowHaplotypes } from './haplotypes';
import type { GenomeBuild } from './ensembl';
import { getAllTranscripts, getProteinDomains, getReference, getRegionGenes, getTranscript } from './ucsc';
import { getCommonSnps } from './snps';
import { getGtexProfile, getGtexTissues } from './gtex';

export interface LocalSample { id: number; name: string; kind: 'bam' | 'cram'; file: File; index: File; /** paths relative to the run folder, when the files came from one */ path?: string; indexPath?: string; /** a sample of an exported page: no file, its regions are embedded in the page */ embedded?: boolean; /** RNA-seq or genomic DNA, and how that was decided */ lib?: LibraryEvidence; /** the file is being opened (header and index read) */ pending?: boolean }
export interface ReferenceChoice { build: GenomeBuild; fasta?: { fa: File; fai: File; gzi?: File } }

type Opened =
  | { kind: 'bam'; bam: BamFile; refNames: string[]; header: string }
  | { kind: 'cram'; cram: IndexedCramFile; refNames: string[]; header: string };

/**
 * Library type from the SAM header: the aligner named in the @PG lines (program name and command line).
 * Spliced aligners mean RNA-seq, genome aligners mean DNA; DRAGEN does both and says which on its command line.
 */
export function classifyHeader(header: string): LibraryEvidence {
  const pg = header.split('\n').filter(l => l.startsWith('@PG'));
  const text = pg.join(' ').toLowerCase();
  const name = (re: RegExp) => re.test(text);
  if (!pg.length) return { type: 'unknown', source: 'none', note: 'no @PG line in the header' };
  if (name(/\bstar\b|starsolo|hisat2|hisat|tophat|\bsubjunc\b|gsnap|olego|mapsplice|crac\b/) || name(/minimap2[^\n]*(-ax? ?splice|-x ?splice)/) || name(/dragen[^\n]*(--enable-rna(?:\s+|=)true|rna)/))
    return { type: 'rna', source: 'header', note: `spliced aligner in the header: ${pgName(pg)}` };
  if (name(/\bbwa\b|bwa-mem2|bwa mem|bowtie2|bowtie|minimap2|isaac|novoalign|ngmlr|winnowmap|pbmm2|dragen/))
    return { type: 'dna', source: 'header', note: `genome aligner in the header: ${pgName(pg)}` };
  return { type: 'unknown', source: 'none', note: `aligner not recognised: ${pgName(pg)}` };
}
const pgName = (pg: string[]): string => {
  const names = pg.map(l => /\tPN:([^\t]+)/.exec(l)?.[1] ?? /\tID:([^\t]+)/.exec(l)?.[1] ?? '').filter(Boolean);
  return [...new Set(names)].slice(0, 3).join(', ') || 'unnamed program';
};

const MAX_REGION_BP = 5_000_000;
const MAX_READS_REGION_BP = 250_000;

// ---------------- Deep regions ----------------
// A very deep library (targeted RNA-seq, a highly expressed gene) can hold millions of records over one gene,
// and decoding them all froze the page. Requests are therefore budgeted: the window is sized from the index
// before anything is decoded, records are scanned tile by tile with only the fields the filter needs, and past
// a cap every k-th read is kept (k = 2, 4, 8…) with the counts scaled back by k.
/** Reads decoded per exon for the exon-usage statistics (fractions and medians only need a sample). */
const EXON_USAGE_MAX_READS = 100_000;
/** widest window whose reference is fetched from the web APIs to place clipped sequences (a local FASTA has no limit) */
const REALIGN_MAX_BP = 500_000;
/** bases read on each side of a breakpoint end when rescuing another sample's clipped reads (longer than a short read, shorter than most long-read clips matter) */
const RESCUE_SPAN_BP = 400;
/**
 * A window is scanned in tiles so that one tile's decoded records can be released before the next is read.
 * The tile is sized from the index (see tileSize) so that about TILE_RECORDS records are decoded at a time
 * whatever the depth: the read cap bounds what a scan *keeps*, never what the library decodes, so a fixed
 * 250 kb tile over a capture panel at a few thousand × materialises millions of records at once.
 * A tile is also decoded in one go on the page's thread, so its size is the longest the page freezes while
 * a window is read: 20 000 records is ~60 ms of decoding (50 000 left 100–380 ms frames when zooming out
 * on 120× genome-like data) at the same throughput (1.9 s against 2.0 s for 1.4 Mb). Deep data sit at the
 * MIN_TILE_BP floor anyway, where this changes nothing.
 */
const TILE_RECORDS = 20_000;
/**
 * The floor is the BAI linear-index interval, 16 kb. A query cannot start reading later than the first
 * record that may overlap its 16 kb bin, so a narrower tile re-inflates (and copies out of the inflater)
 * up to 16 kb of data it has no use for, once per tile. Measured on a 1.5 M-read RNA gene at ~100 000×:
 * 2 kb tiles took 19.6 s against 3.2 s for 16 kb ones — the extra time all in BGZF chunk inflation
 * and buffer copies, not in the records.
 */
const MAX_TILE_BP = 250_000, MIN_TILE_BP = 16_384;
/** Compressed bytes per record assumed before a file has been scanned once (a scan then calibrates it). */
const BYTES_PER_READ: Record<'bam' | 'cram', number> = { bam: 60, cram: 30 };
/** Margins are halved while the window looks too deep; below this they are dropped altogether. */
const MIN_MARGIN_BP = 2_000;
/**
 * Reads the margins of a coverage request may cost, by the index's estimate, unless the caller says
 * otherwise. The view itself is always read in full and exactly; margins only make panning free, so
 * past this they shrink (about 5–10 s of decoding on one thread, measured at 200–400 k reads/s).
 */
const DEFAULT_MARGIN_READS = 2_000_000;
/** Longest stretch of one chromosome whose counts a sample keeps; a request further away starts over. */
const COVERAGE_MAX_SPAN = 8_000_000;
/** Counts kept across all samples; the least recently used samples' are dropped past it. */
const COVERAGE_CACHE_BYTES = 256 * 1024 * 1024;
/** A request this close to what a sample has counted extends it (the gap is decoded) rather than starting over. */
const COVERAGE_MIN_GAP = 100_000;
/** Longest stretch of one chromosome whose allele counts a sample keeps (they weigh ~9× the coverage counts per base). */
const ALLELE_MAX_SPAN = 4_000_000;
/** Widest window whose methylation is asked at once, longest stretch a sample keeps, and the counts kept across samples. */
export const METHYL_MAX_BP = 2_000_000;
const METHYL_MAX_SPAN = 6_000_000, METHYL_CACHE_BYTES = 256 * 1024 * 1024;
/** Allele counts kept across all samples; the least recently used samples' are dropped past it. */
const ALLELE_CACHE_BYTES = 384 * 1024 * 1024;
/** Partial coverage is handed to the caller at most this often while a deep view fills. */
const PROGRESS_MS = 200;
/** Records looked at for an NH tag before a file without one is taken as having none (uniqueness then from MAPQ). */
const NH_PROBE = 2_000;
/**
 * Decoded records the alignment libraries keep, all files of this data source together (the page and the variant
 * worker each have one). The libraries default to 1 GB per file, sized for a genome browser that decodes the same
 * chunks again at every pan; here coverage, variants and methylation are counted once into states of their own, and
 * only the reads track reads a window again, so the cache mostly held chunks nothing would reuse: 210 MB on a
 * 1 000× gene before any read was shown, times two (page and worker), times every sample. One budget across files
 * (SharedBudget: the least recently used chunk of any file goes first), and a chunk not looked at for
 * RECORD_CACHE_IDLE_MS dropped, so a parked tab gives it back.
 */
const RECORD_CACHE_BYTES = 128 * 1024 * 1024, RECORD_CACHE_IDLE_MS = 45_000;

/**
 * Uniform access to one decoded record of either library. Nothing heavy (name, sequence, qualities) is touched
 * unless `raw` is asked for a full record: BAM records decode their fields on demand, so a scan that reads only
 * flags and MAPQ stays cheap even over millions of records.
 */
interface RecordView<R> {
  start(r: R): number;
  flags(r: R): number;
  mapq(r: R): number;
  nh(r: R): number | null;
  /** the CIGAR packed as BAM stores it (length << 4 | op) */
  ops(r: R): ArrayLike<number>;
  /** reference id of the mate, −1 when none */
  mateRef(r: R): number;
  tlen(r: R): number;
  sa(r: R): unknown;
  /** the read's bases as character codes into `buf.a` (grown when too short); returns their number, 0 without a sequence */
  seqCodes(r: R, buf: { a: Uint8Array }): number;
  /** base qualities, null when the record has none */
  quals(r: R): ArrayLike<number> | null;
  /** the base-modification tags (MM / ML, or the legacy Mm / Ml), null when the record has none */
  mods(r: R): { mm: string; ml: ArrayLike<number> | null; mn: number | null } | null;
  /** `light` leaves out name, sequence and qualities (coverage only needs the alignment blocks); `structural` adds the pair and SA fields. */
  raw(r: R, light: boolean, structural: boolean, refNames: string[]): RawRead;
}
/** The MM / ML tags of a record (the legacy Mm / Ml names too), with MN (the SEQ length they describe) when present. */
const modTags = (get: (tag: string) => unknown): { mm: string; ml: ArrayLike<number> | null; mn: number | null } | null => {
  const mm = get('MM') ?? get('Mm');
  if (typeof mm !== 'string' || !mm) return null;
  const ml = get('ML') ?? get('Ml');
  return { mm, ml: ml != null && typeof ml === 'object' && 'length' in (ml as object) ? ml as ArrayLike<number> : null, mn: tagNumber(get('MN')) };
};
/** The haplotag of a record of a phased file: HP (haplotype), PS (phase set), PC (confidence); nothing when untagged. */
const haplotag = (get: (tag: string) => unknown): { hp?: number; ps?: number; pc?: number } => {
  const hp = tagNumber(get('HP'));
  if (hp == null || hp <= 0) return {};
  const ps = tagNumber(get('PS')), pc = tagNumber(get('PC'));
  return { hp, ...(ps != null ? { ps } : {}), ...(pc != null ? { pc } : {}) };
};
const mateFields = (chromOf: (id: number) => string, mateId: number, matePos: number, tlen: number, sa: unknown) => ({
  tlen, mateChrom: mateId >= 0 ? chromOf(mateId) : '', matePos: mateId >= 0 ? matePos : undefined, sa: typeof sa === 'string' ? sa : null,
});
/** A CIGAR with a soft clip of RESCUE_MIN_CLIP bases or more: the light structural scan decodes this record's sequence so the clip can be placed by realignment or rescued at a known breakpoint. */
const CLIP_RE = new RegExp(`(?:^|[A-Z=])(\\d+)S`, 'g');
const bigClip = (cigar: string) => { CLIP_RE.lastIndex = 0; let m: RegExpExecArray | null; while ((m = CLIP_RE.exec(cigar))) if (parseInt(m[1]) >= RESCUE_MIN_CLIP) return true; return false; };
/** BAM 4-bit base codes (SAM spec: =ACMGRSVTWYHKDBN) as character codes. */
const NIBBLE_CODES = Uint8Array.from('=ACMGRSVTWYHKDBN', c => c.charCodeAt(0));
const BAM_VIEW: RecordView<any> = {
  start: r => r.start, flags: r => r.flags, mapq: r => r.mq ?? 255, nh: r => tagNumber(r.getTag('NH')),
  ops: r => r.NUMERIC_CIGAR, mateRef: r => r.next_refid, tlen: r => r.template_length, sa: r => r.getTag('SA'),
  seqCodes: (r, buf) => {
    const n = r.seq_length ?? 0;
    if (!n) return 0;
    if (buf.a.length < n) buf.a = new Uint8Array(n * 2);
    const packed: Uint8Array = r.NUMERIC_SEQ, a = buf.a;
    for (let i = 0; i < n; i++) { const b = packed[i >> 1]; a[i] = NIBBLE_CODES[i & 1 ? b & 15 : b >> 4]; }
    return n;
  },
  quals: r => r.qual ?? null,
  mods: r => modTags(t => r.getTag(t)),
  raw: (r, light, structural, refNames) => {
    const sa = structural ? r.getTag('SA') : undefined;
    const withSeq = !light || (structural && typeof sa !== 'string' && bigClip(r.CIGAR));
    // the name ties the parts of a split read, and the two mates of a pair, together: kept in the light structural scan
    return { name: light && !structural ? '' : r.name, start: r.start, cigar: r.CIGAR, seq: withSeq ? r.seq : '', qual: light ? null : r.qual, flags: r.flags, mapq: r.mq ?? 255, nh: tagNumber(r.getTag('NH')),
      ...(light ? {} : { ...haplotag(tag => r.getTag(tag)), mods: modTags(t => r.getTag(t)) }),
      ...(structural ? mateFields(id => refNames[id] ?? '', r.next_refid, r.next_pos, r.template_length, sa) : {}) };
  },
};
/**
 * A CRAM record's bases. @gmod/cram stores only the differences from the reference and rebuilds the read in
 * getReadBases(); its `readBases` property stays undefined until that has been called.
 */
const cramBases = (r: any): string => (typeof r.getReadBases === 'function' ? r.getReadBases() : r.readBases) ?? '';
const CRAM_VIEW: RecordView<any> = {
  start: r => r.start, flags: r => r.flags, mapq: r => r.mappingQuality ?? 255, nh: r => tagNumber(r.getTag('NH')),
  ops: r => packCigar(cramCigar(r.readFeatures, r.readLength, r.lengthOnRef ?? 0)),
  mateRef: r => r.nextSequenceId ?? -1, tlen: r => r.templateLength ?? r.templateSize ?? 0, sa: r => r.getTag('SA'),
  seqCodes: (r, buf) => {
    const seq = cramBases(r);
    if (buf.a.length < seq.length) buf.a = new Uint8Array(seq.length * 2);
    for (let i = 0; i < seq.length; i++) buf.a[i] = seq.charCodeAt(i);
    return seq.length;
  },
  quals: r => r.qualityScores ?? null,
  mods: r => modTags(t => r.getTag(t)),
  raw: (r, light, structural, refNames) => {
    const feats = r.readFeatures as any;
    const qual = r.qualityScores ?? null;
    const cigar = cramCigar(feats, r.readLength, r.lengthOnRef ?? 0);
    const sa = structural ? r.getTag('SA') : undefined;
    const withSeq = !light || (structural && typeof sa !== 'string' && bigClip(cigar));
    return { name: light && !structural ? '' : (r.readName ?? ''), start: r.start, cigar, seq: withSeq ? cramBases(r) : '', qual: light ? null : qual, flags: r.flags, mapq: r.mappingQuality ?? 255, nh: tagNumber(r.getTag('NH')),
      mismatches: light ? undefined : cramMismatches(feats, qual),
      ...(light ? {} : { ...haplotag(tag => r.getTag(tag)), mods: modTags(t => r.getTag(t)) }),
      ...(structural ? mateFields(id => refNames[id] ?? '', r.nextSequenceId ?? -1, r.nextStart ?? 0, r.templateLength ?? r.templateSize ?? 0, sa) : {}) };
  },
};

/** Result of a budgeted scan: `kept` holds every `rate`-th of the `total` reads that passed the filters. */
interface Scan { total: number; rate: number; kept: RawRead[] }
/** What a scan needs beyond its window: the structural fields of each record, and the caller's abort signal. */
interface ScanOptions { structural?: boolean; signal?: AbortSignal }

/** Rejects with the standard AbortError once the caller has given up on the request. */
function throwIfAborted(signal?: AbortSignal): void {
  if (signal?.aborted) throw new DOMException('Request aborted', 'AbortError');
}
/**
 * Hands the thread back to the browser between tiles, so a deep window keeps painting and the buttons
 * stay clickable while it is read. scheduler.yield() is Chromium-only (Chrome / Edge 129+), hence the fallback.
 */
const yieldToUi = (): Promise<void> => {
  const s = (globalThis as { scheduler?: { yield?: () => Promise<void> } }).scheduler;
  return s?.yield ? s.yield() : new Promise<void>(resolve => { setTimeout(resolve, 0); });
};

function resolveName(names: string[], chrom: string): string | null {
  if (names.includes(chrom)) return chrom;
  const alt = chrom.startsWith('chr') ? chrom.slice(3) : `chr${chrom}`;
  return names.includes(alt) ? alt : null;
}

function tagNumber(v: unknown): number | null {
  return typeof v === 'number' && Number.isFinite(v) ? v : null;
}
/** Long reads (ONT, PacBio): the median aligned length of the window's reads is above 1 kb. */
/** Adds the pair fields of a paired, mate-mapped record to its encoded read: mate start, template length, mate chromosome when it differs. */
/**
 * The CpG calls of each read (`me`: position, P(5mC) × 255, flat), for the methylation colours of the reads track; the
 * same rules as the counting (reference CpGs only, both strands on the C of the + strand, MN / hard-clip check).
 */
function readMethylation(reads: AlignedRead[], raw: RawRead[], cpgs: Int32Array): void {
  const sc = newModScratch();
  let codes = new Uint8Array(4096);
  raw.forEach((r, i) => {
    const m = r.mods, n = r.seq.length;
    if (!m || !n || !keepFlags(r.flags)) return;
    const ops = packCigar(r.cigar);
    if (m.mn != null ? m.mn !== n : ops.some(v => (v & 15) === 5)) return;
    if (codes.length < n) codes = new Uint8Array(n * 2);
    for (let k = 0; k < n; k++) codes[k] = r.seq.charCodeAt(k);
    const me: number[] = [];
    visitReadCalls(r.start, ops, codes, n, (r.flags & 16) !== 0, m.mm, m.ml, cpgs, sc, (cpg, mod, conf) => { me.push(cpg, Math.round((mod ? conf : 1 - conf) * 255)); });
    if (me.length) reads[i].me = me;
  });
}

function withMate(a: AlignedRead, r: RawRead, own: string): AlignedRead {
  if (!(r.flags & 1) || r.flags & 8 || r.matePos == null || r.matePos < 0) return a;
  a.mp = r.matePos; a.tl = r.tlen ?? 0;
  if (r.mateChrom && r.mateChrom !== own) a.mc = r.mateChrom;
  return a;
}

export function isLongRead(reads: { s: number; e: number }[]): boolean {
  if (!reads.length) return false;
  const lens = reads.map(r => r.e - r.s).sort((a, b) => a - b);
  return lens[lens.length >> 1] > 1000;
}

/** Runs full variant scans away from the page's thread (see variantClient.ts); the data source hands them over. */
export interface VariantScanner {
  scan(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan>;
  /** the sample's files changed or it was removed: what was counted for it no longer holds */
  forget(sampleId: number): void;
  /** the reference changed: mismatches are read against it */
  reference(reference: ReferenceChoice): void;
  /** CpG methylation counted in the worker (absent from scanners that do not do it) */
  methyl?(sampleId: number, chrom: string, start: number, end: number, opts?: { signal?: AbortSignal; onProgress?: (fraction: number) => void }): Promise<MethylWindow>;
  /** drop the counts kept for these */
  release?(what: 'methylation' | 'variants'): void;
}

export class LocalDataSource implements SashimiDataSource {
  private samples = new Map<number, LocalSample>();
  private opened = new Map<number, Promise<Opened>>();
  private fasta: Promise<{ fa: IndexedFasta | BgzipIndexedFasta; names: string[] }> | null = null;
  reference: ReferenceChoice;
  /** Set by the shell to know which reference source actually answered. */
  lastReferenceSource: 'fasta' | 'ensembl' | null = null;

  constructor(reference: ReferenceChoice) {
    this.reference = reference;
  }

  setReference(reference: ReferenceChoice) {
    this.reference = reference;
    this.fasta = null;
    // CRAM decoding depends on the reference: reopen files, and forget what was counted from them
    for (const [id, s] of this.samples) if (s.kind === 'cram') { this.opened.delete(id); this.coverage.delete(id); }
    // mismatches are read against the reference: every sample's allele counts depend on it
    this.alleleStates.clear();
    this.variantScanner?.reference(reference);
  }

  /** Variants handed over by the page URL (deep link), drawn on every sample. */
  knownVariants: KnownVariant[] = [];
  async getKnownVariants(_sampleId: number): Promise<KnownVariant[]> { return this.knownVariants; }

  addSample(s: LocalSample) { this.samples.set(s.id, s); this.coverage.delete(s.id); this.alleleStates.delete(s.id); this.methylStates.delete(s.id); this.nhMode.delete(s.id); this.variantScanner?.forget(s.id); }
  renameSample(id: number, name: string) { const s = this.samples.get(id); if (s) this.samples.set(id, { ...s, name }); }
  removeSample(id: number) { this.samples.delete(id); this.opened.delete(id); this.headers.delete(id); this.coverage.delete(id); this.alleleStates.delete(id); this.methylStates.delete(id); this.nhMode.delete(id); this.variantScanner?.forget(id); }
  /** one byte budget for the decoded records of every alignment file of this source (see RECORD_CACHE_BYTES) */
  private recordBudget = new SharedBudget(RECORD_CACHE_BYTES);
  /**
   * Drops what was kept for an overlay the viewer switched off (it is read or counted again when it comes back): the
   * methylation or allele counts (here and in the worker), or the decoded records of every file (`records`: the
   * reads track was closed, or the coverage they were decoded for is counted).
   */
  release(what: 'methylation' | 'variants' | 'records') {
    if (what === 'records') { this.clearRecordCaches(); return; }
    if (what === 'methylation') this.methylStates.clear(); else this.alleleStates.clear();
    this.variantScanner?.release?.(what);
  }
  /** Empties the alignment libraries' record caches (every opened file). */
  clearRecordCaches(): void {
    for (const o of this.opened.values()) o.then(f => (f.kind === 'bam' ? f.bam.clearFeatureCache() : f.cram.clearFeatureCache())).catch(() => { /* a file that failed to open holds nothing */ });
  }
  /** The files behind a sample, for a variant scanner that reads them elsewhere (a worker). */
  sampleFiles(id: number): LocalSample | undefined { return this.samples.get(id); }
  list(): SampleRef[] { return [...this.samples.values()].map(s => ({ id: s.id, name: s.name })); }

  // ---- reference ----
  private async openFasta() {
    if (!this.fasta) {
      const f = this.reference.fasta!;
      this.fasta = (async () => {
        const fa = f.gzi
          ? new BgzipIndexedFasta({ fasta: new BlobFile(f.fa), fai: new BlobFile(f.fai), gzi: new BlobFile(f.gzi) })
          : new IndexedFasta({ fasta: new BlobFile(f.fa), fai: new BlobFile(f.fai) });
        const names = await fa.getSequenceNames();
        return { fa, names };
      })();
    }
    return this.fasta;
  }

  async getReferenceSeq(chrom: string, start: number, end: number): Promise<string | null> {
    if (end <= start) return '';
    if (this.reference.fasta) {
      try {
        const { fa, names } = await this.openFasta();
        const name = resolveName(names, chrom);
        if (name) {
          const s = await fa.getSequence(name, Math.max(0, start), end);
          if (s != null) { this.lastReferenceSource = 'fasta'; return s.toUpperCase(); }
        }
      } catch (e) { console.warn('FASTA read failed, falling back to the UCSC / Ensembl APIs:', e); }
    }
    const s = await getReference(this.reference.build, chrom, Math.max(0, start), end);
    this.lastReferenceSource = s == null ? null : 'ensembl';
    return s;
  }

  // ---- files ----
  /** SAM header text of each sample, read without its index (memoised). */
  private headers = new Map<number, Promise<string>>();

  /**
   * The SAM header of a BAM from the start of the file alone: the header sits before the first record, so a few bgzf
   * blocks give it. The library's own header read first parses the whole index to size the header (5–10 MB for a
   * genome), which is what made adding a file on a slow disk or share wait; the index is parsed later, at the first
   * region request. A CRAM's header is in its first container, read the same way by the CRAM library.
   */
  private headerText(id: number): Promise<string> {
    const s = this.samples.get(id);
    if (!s) return Promise.reject(new Error('Sample not found'));
    if (!this.headers.has(id)) {
      this.headers.set(id, (async () => {
        if (s.kind === 'cram') return (await this.open(id)).header;
        // read, inflate, and grow the read until the header text is complete (64 kB covers most files; a header of
        // thousands of contigs and hundreds of @PG lines takes a few hundred kB)
        for (let len = 64 * 1024; ; len *= 4) {
          const bytes = new Uint8Array(await s.file.slice(0, Math.min(len, s.file.size)).arrayBuffer());
          const data = await unzip(bytes);   // whole blocks only; a block cut by the read end is left for the next round
          if (data.length >= 8) {
            if (!(data[0] === 66 && data[1] === 65 && data[2] === 77 && data[3] === 1)) throw new Error(`${s.file.name} is not a BAM file`);
            const lText = new DataView(data.buffer, data.byteOffset, data.byteLength).getInt32(4, true);
            if (lText < 0) throw new Error(`${s.file.name}: invalid BAM header`);
            if (8 + lText <= data.length) return new TextDecoder().decode(data.subarray(8, 8 + lText)).replace(/\0+$/, '');
          }
          if (len >= s.file.size || len >= 64 * 1024 * 1024) throw new Error(`${s.file.name}: BAM header not found`);
        }
      })().catch(e => { this.headers.delete(id); throw e; }));
    }
    return this.headers.get(id)!;
  }

  private open(id: number): Promise<Opened> {
    const s = this.samples.get(id);
    if (!s) return Promise.reject(new Error('Sample not found'));
    if (!this.opened.has(id)) {
      this.opened.set(id, (async (): Promise<Opened> => {
        if (s.kind === 'bam') {
          const bam = new BamFile({ bamFilehandle: new BlobFile(s.file), baiFilehandle: new BlobFile(s.index), maxCacheBytes: RECORD_CACHE_BYTES, cacheIdleTimeoutMs: RECORD_CACHE_IDLE_MS, cacheBudget: this.recordBudget });
          await bam.getHeader();
          const refNames = (bam.indexToChr || []).map(r => r.refName);
          const header = (await bam.getHeaderText().catch(() => '')) ?? '';
          return { kind: 'bam', bam, refNames, header };
        }
        const cram = new IndexedCramFile({
          cramFilehandle: new BlobFile(s.file),
          index: new CraiIndex({ filehandle: new BlobFile(s.index) }),
          checkSequenceMD5: false,
          useSliceWorkerPool: false,
          maxCacheBytes: RECORD_CACHE_BYTES, cacheIdleTimeoutMs: RECORD_CACHE_IDLE_MS, cacheBudget: this.recordBudget,
          fetchReferenceSequence: async (seqId: number, start: number, end: number, refName?: string) => {
            const info = await cram.cram.getReferenceInfo();
            const name = refName ?? info[seqId]?.name ?? '';
            const seq = await this.getReferenceSeq(name, start, end);
            if (seq == null || seq.length !== end - start) throw new Error(`No reference sequence for ${name}:${start}-${end} (CRAM decoding needs it)`);
            return seq;
          },
        });
        const info = await cram.cram.getReferenceInfo();
        const header = (await cram.cram.getHeaderText().catch(() => '')) ?? '';
        return { kind: 'cram', cram, refNames: info.map(r => r.name), header };
      })().catch(e => { this.opened.delete(id); throw e; }));
    }
    return this.opened.get(id)!;
  }

  /** Compressed bytes per record of a file, from its last scan (a size estimate before anything is decoded). */
  private bytesPerRead = new Map<number, number>();

  private async locate(id: number, chrom: string): Promise<{ o: Opened; name: string; seqId: number } | null> {
    const o = await this.open(id);
    const name = resolveName(o.refNames, chrom);
    return name ? { o, name, seqId: o.refNames.indexOf(name) } : null;
  }

  /**
   * Compressed bytes the index says [start, end) occupies: BAI/CSI for BAM, CRAI slices for CRAM.
   * `estimatedBytesForRegions` is what the sizing wants — summing every chunk `blocksForRange` returns
   * charges a narrow window for every chunk of every overlapping bin, which on a long-read file is
   * documented as up to 5.6× the bytes the query really reads, and that is exactly the window a reader
   * spends their time in. The sum stays as the fallback for an index the library did not expose.
   */
  private async indexBytes(loc: { o: Opened; name: string; seqId: number }, start: number, end: number, signal?: AbortSignal): Promise<number> {
    throwIfAborted(signal);
    if (loc.o.kind === 'bam') {
      const index = loc.o.bam.index;
      if (index) return index.estimatedBytesForRegions([{ refId: loc.seqId, start, end }], { signal });
      const chunks = await loc.o.bam.blocksForRange(loc.name, start, end, { signal });
      return chunks.reduce((a, c) => a + c.fetchedSize(), 0);
    }
    const slices = await loc.o.cram.index.getEntriesForRange(loc.seqId, start, end);
    return slices.reduce((a, sl) => a + sl.sliceBytes, 0);
  }

  /**
   * Tile width for a scan of [start, end): the window narrowed so that about TILE_RECORDS records are
   * decoded (and released) at a time, from the compressed bytes the index reports for the whole window.
   *
   * Narrow tiles cost twice. Each is read with every record overlapping it, so a read on a boundary is
   * decoded once per tile it spans (spliced reads with long N gaps span several); and each starts at the
   * index's 16 kb granularity, so the bytes before it are inflated again (see MIN_TILE_BP). The floor
   * keeps both small; a long-read library never gets near it, the width falling only when the window
   * holds many records.
   */
  private tileSize(id: number, start: number, end: number, bytes: number): number {
    const s = this.samples.get(id);
    const bpr = this.bytesPerRead.get(id) ?? BYTES_PER_READ[s?.kind ?? 'bam'];
    const records = bytes / Math.max(1, bpr);
    if (!(records > TILE_RECORDS)) return MAX_TILE_BP;
    const span = Math.max(1, end - start);
    return Math.max(MIN_TILE_BP, Math.min(MAX_TILE_BP, Math.floor((span * TILE_RECORDS) / records)));
  }

  /**
   * Reads of [start, end) passing the flag (and uniqueness) filters, decoded tile by tile, at most about `cap`
   * of them: when the kept reads outgrow the cap they are thinned to every other one and the rate doubles, so
   * `kept` is always the reads whose rank (among the passing reads, in file order) is a multiple of `rate` — a
   * systematic sample that is exact (rate 1) whenever the region holds no more than `cap` reads. A read spanning
   * two tiles is counted in the tile holding its start (or the first tile when it starts before the window).
   * The tile width comes from the index, and the thread is handed back between tiles, so neither the memory
   * nor the pauses grow with the depth of the library. `opts.signal` drops the decoding and the fetch in flight.
   */
  private async scan(id: number, chrom: string, start: number, end: number, uniqueOnly: boolean, cap: number, light: boolean, opts: ScanOptions = {}): Promise<Scan> {
    const { structural = false, signal } = opts;
    const loc = await this.locate(id, chrom);
    if (!loc) return { total: 0, rate: 1, kept: [] };
    const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
    let kept: RawRead[] = [], total = 0, rate = 1, seen = 0;
    const bytes = await this.indexBytes(loc, start, end, signal);
    const tileBp = this.tileSize(id, start, end, bytes);
    for (let ts = start; ts < end; ts += tileBp) {
      throwIfAborted(signal);
      const te = Math.min(end, ts + tileBp);
      const recs: any[] = loc.o.kind === 'bam'
        ? await loc.o.bam.getRecordsForRange(loc.name, ts, te, { signal })
        : await loc.o.cram.getRecordsForRange(loc.seqId, ts, te, { signal });
      for (const r of recs) {
        if (ts > start && view.start(r) < ts) continue; // already counted in the previous tile
        seen++;
        if (!keepFlags(view.flags(r))) continue;
        if (uniqueOnly && !uniqueFrom(view.nh(r), view.mapq(r))) continue;
        if (total % rate === 0) {
          kept.push(view.raw(r, light, structural, loc.o.refNames));
          if (kept.length > cap) { kept = kept.filter((_, i) => i % 2 === 0); rate *= 2; }
        }
        total++;
      }
      if (te < end) await yieldToUi();
    }
    if (seen >= 1000 && bytes > 0) this.bytesPerRead.set(id, bytes / seen);
    return { total, rate, kept };
  }

  /**
   * The window to read for a request: the core is always read; the margins around it are halved while the
   * index suggests they hold more reads than the budget, and dropped once they get small. Only the margins
   * count against it: around a capture target at thousands of ×, margins that hold next to nothing stay.
   */
  private async budgetWindow(id: number, chrom: string, start: number, end: number, core: { start: number; end: number }, cap: number, signal?: AbortSignal): Promise<{ start: number; end: number }> {
    const loc = await this.locate(id, chrom);
    const s = this.samples.get(id);
    if (!loc || !s) return { start, end };
    const bpr = this.bytesPerRead.get(id) ?? BYTES_PER_READ[s.kind];
    const coreStart = Math.max(start, Math.min(end, core.start)), coreEnd = Math.max(coreStart, Math.min(end, core.end));
    let left = coreStart - start, right = end - coreEnd;
    let win = { start, end };
    const coreBytes = await this.indexBytes(loc, coreStart, coreEnd, signal);
    while ((left > 0 || right > 0) && ((await this.indexBytes(loc, win.start, win.end, signal)) - coreBytes) / bpr > cap) {
      left = left >= 2 * MIN_MARGIN_BP ? Math.floor(left / 2) : 0;
      right = right >= 2 * MIN_MARGIN_BP ? Math.floor(right / 2) : 0;
      win = { start: coreStart - left, end: coreEnd + right };
    }
    return win;
  }

  // ---- SashimiDataSource ----
  /** Clipped reads of this sample rescued at breakpoints seen in other samples: only the reads around the breakpoint ends are read. */
  async rescueClips(sampleId: number, chrom: string, breakpoints: Breakpoint[], uniqueOnly: boolean): Promise<RescuedClips[]> {
    if (!breakpoints.length) return [];
    const loc = await this.locate(sampleId, chrom);
    if (!loc) return [];
    const ends = [...new Set(breakpoints.flatMap(b => [b.start, b.end]))].sort((a, b) => a - b);
    const lo = Math.max(0, ends[0] - RESCUE_SPAN_BP), hi = ends[ends.length - 1] + RESCUE_SPAN_BP;
    if (!this.reference.fasta && hi - lo > REALIGN_MAX_BP) return [];
    // reads clipped at an end: those overlapping a short stretch around each end (merged when close)
    const ranges: { start: number; end: number }[] = [];
    for (const e of ends) { const a = Math.max(0, e - RESCUE_SPAN_BP), b = e + RESCUE_SPAN_BP; const last = ranges[ranges.length - 1]; if (last && a <= last.end) last.end = b; else ranges.push({ start: a, end: b }); }
    // full records (names included) over short stretches: a read overlapping two stretches is kept once
    const seen = new Set<string>(); const reads: RawRead[] = [];
    for (const r of ranges) {
      const { kept } = await this.scan(sampleId, chrom, r.start, r.end, uniqueOnly, Number.MAX_SAFE_INTEGER, false, { structural: true });
      for (const x of kept) { const k = `${x.name}:${x.start}:${x.flags}:${x.cigar}`; if (!seen.has(k)) { seen.add(k); reads.push(x); } }
    }
    const clipped = clipEnds(reads).filter(e => ends.some(p => Math.abs(p - e.pos) <= RESCUE_TOLERANCE_BP));
    if (!clipped.length) return [];
    const seq = await this.getReferenceSeq(chrom, lo, hi);
    const found = rescueClipEnds(clipped, breakpoints, seq ? { start: lo, seq } : null);
    return [...found].map(([k, n]) => { const [kind, span] = k.split(':'); const [a, b] = span.split('-').map(Number); return { start: a, end: b, kind: kind as Breakpoint['kind'], count: n.count, hard: n.hard, own: false }; });
  }
  async getPrimaryRecord(sampleId: number, chrom: string, start: number, name: string): Promise<{ seq: string; flags: number; cigar: string } | null> {
    const { kept } = await this.scan(sampleId, chrom, start, start + 1, false, Number.MAX_SAFE_INTEGER, false, { structural: true });
    const r = kept.find(x => x.start === start && x.name === name && !(x.flags & 2048) && x.seq);
    return r ? { seq: r.seq, flags: r.flags, cigar: r.cigar } : null;
  }
  // ---------------- Variant sites: allele counts, streamed and kept (see alleles.ts) ----------------

  /**
   * Where full variant scans run. The shell points this at a Web Worker that hosts its own LocalDataSource (see
   * variantWorker.ts): the scan then never holds the page's thread, however deep the window. Unset (and inside
   * the worker), scans run here, tile by tile.
   */
  variantScanner?: VariantScanner;
  /** What each sample has counted for its variants, one chromosome stretch per sample. */
  private alleleStates = new Map<number, AlleleState>();
  private methylStates = new Map<number, MethylState>();
  private methylLocks = new Map<number, Promise<unknown>>();
  private alleleLocks = new Map<number, Promise<unknown>>();

  /**
   * Every variant site above the thresholds in [start, end), from every read of the window: exact counts, no
   * sampling. Handed to the variant scanner when one is set, counted here otherwise.
   */
  getVariantSites(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan> {
    const s = this.samples.get(sampleId);
    if (!s) return Promise.reject(new Error('Sample not found'));
    if (this.variantScanner && !s.embedded) return this.variantScanner.scan(sampleId, chrom, start, end, uniqueOnly, minVaf, opts);
    return this.countVariants(sampleId, chrom, start, end, uniqueOnly, minVaf, opts);
  }

  /** The sample's allele state for a request on `chrom` around [start, end): kept when the request is near it, started over otherwise. */
  private alleleState(id: number, chrom: string, start: number, end: number): AlleleState {
    let st = this.alleleStates.get(id);
    const gap = Math.max(end - start, COVERAGE_MIN_GAP);
    if (!st || st.chrom !== chrom || (!st.empty && (start > st.pe + gap || end < st.ps - gap))
      || (!st.empty && Math.max(st.pe, end) - Math.min(st.ps, start) > ALLELE_MAX_SPAN)) {
      st = new AlleleState(chrom);
      this.alleleStates.set(id, st);
    }
    st.lastUsed = Date.now();
    let total = 0;
    for (const x of this.alleleStates.values()) total += x.bytes;
    for (const [k, x] of [...this.alleleStates.entries()].filter(([k]) => k !== id).sort((p, q) => p[1].lastUsed - q[1].lastUsed)) {
      if (total <= ALLELE_CACHE_BYTES) break;
      total -= x.bytes; this.alleleStates.delete(k);
    }
    return st;
  }

  /** CpG methylation of a window: counted in the worker when there is one (the reads carry their whole sequence), else here. */
  getMethylation(sampleId: number, chrom: string, start: number, end: number, opts?: { signal?: AbortSignal; onProgress?: (fraction: number) => void }): Promise<MethylWindow> {
    const s = this.samples.get(sampleId);
    if (!s) return Promise.reject(new Error('Sample not found'));
    if (s.embedded) return Promise.reject(new Error('methylation is not part of exported pages: add the alignment file'));
    if (end - start > METHYL_MAX_BP) return Promise.reject(new Error(`Region too large for methylation (${(end - start).toLocaleString()} bp; up to ${METHYL_MAX_BP.toLocaleString()})`));
    if (this.variantScanner?.methyl) return this.variantScanner.methyl(sampleId, chrom, start, end, opts);
    return this.countMethylation(sampleId, chrom, start, end, opts);
  }

  /** The sample's methylation state for a request around [start, end): kept when near, started over otherwise; old samples' states dropped past the cache. */
  private methylState(id: number, chrom: string, start: number, end: number): MethylState {
    let st = this.methylStates.get(id);
    const gap = Math.max(end - start, COVERAGE_MIN_GAP);
    if (!st || st.chrom !== chrom || (!st.empty && (start > st.pe + gap || end < st.ps - gap))
      || (!st.empty && Math.max(st.pe, end) - Math.min(st.ps, start) > METHYL_MAX_SPAN)) {
      st = new MethylState(chrom);
      this.methylStates.set(id, st);
    }
    st.lastUsed = Date.now();
    let total = 0;
    for (const x of this.methylStates.values()) total += x.bytes;
    for (const [k, x] of [...this.methylStates.entries()].filter(([k]) => k !== id).sort((p, q) => p[1].lastUsed - q[1].lastUsed)) {
      if (total <= METHYL_CACHE_BYTES) break;
      total -= x.bytes; this.methylStates.delete(k);
    }
    return st;
  }

  /**
   * The methylation counting itself (methylation.ts): grows the sample's state over [start, end), each tile's reads
   * counted against the CpGs of the reference under them, then reads the window back with the filter threshold.
   */
  countMethylation(sampleId: number, chrom: string, start: number, end: number, opts?: { signal?: AbortSignal; onProgress?: (fraction: number) => void }): Promise<MethylWindow> {
    const signal = opts?.signal;
    return this.locked(this.methylLocks, sampleId, async () => {
      throwIfAborted(signal);
      const loc = await this.locate(sampleId, chrom);
      if (!loc) return methylWindow([], start, end, null, start);
      const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
      const st = this.methylState(sampleId, loc.name, start, end);
      const scratch = { a: new Uint8Array(4096) }, mods = newModScratch();
      let cpgs: Int32Array = new Int32Array(0);
      // the CpGs of the reference under the reads of the tile about to be counted, before the (synchronous) counting
      const prepare = async (recs: any[]) => {
        let lo = Infinity, hi = -Infinity;
        for (const r of recs) {
          const rs = view.start(r);
          if (rs < lo) lo = rs;
          const ops = view.ops(r);
          let e = rs;
          for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) e += ops[k] >>> 4; }
          if (e > hi) hi = e;
        }
        const seq = hi > lo ? await this.getReferenceSeq(chrom, lo, hi + 1) : null;
        cpgs = seq ? cpgSites(lo, seq) : new Int32Array(0);
      };
      const count = (layer: MethylCounts, r: any) => {
        if (!keepFlags(view.flags(r))) return;
        const m = view.mods(r);
        if (!m) return;
        const ops = view.ops(r);
        const n = view.seqCodes(r, scratch);
        if (!n) return;
        // the tags describe the read as sequenced: a hard-clipped record matches them only when MN says so
        if (m.mn != null ? m.mn !== n : Array.from(ops).some(v => (v & 15) === 5)) return;
        const hp = tagNumber(r.getTag('HP')) ?? 0;
        countRead(layer, view.start(r), ops, scratch.a, n, (view.flags(r) & 16) !== 0, m.mm, m.ml, cpgs, hp, mods);
      };
      const span = Math.max(1, end - start);
      const progress = () => opts?.onProgress?.(Math.min(1, Math.max(0, Math.min(end, st.pe) - Math.max(start, st.ps)) / span));
      await this.fillStretch(sampleId, st, loc, start, end, () => new MethylCounts(), count, signal, progress, prepare);
      progress();
      const ref = await this.getReferenceSeq(chrom, start, end + 1);
      throwIfAborted(signal);
      return methylWindow([st.owned, st.spill], start, end, ref, start);
    });
  }

  /** The allele counting itself: grows the sample's allele state over [start, end), then reads the sites back. */
  countVariants(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan> {
    const signal = opts?.signal;
    return this.locked(this.alleleLocks, sampleId, async () => {
      throwIfAborted(signal);
      const loc = await this.locate(sampleId, chrom);
      if (!loc) return { sites: [], total: 0, long_reads: false };
      const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
      const st = this.alleleState(sampleId, loc.name, start, end);
      const scratch = { a: new Uint8Array(512) };
      // the reference over the records of the tile about to be counted, fetched before the (synchronous) counting
      let ref: RefWindow | null = null;
      const prepare = async (recs: any[]) => {
        let lo = Infinity, hi = -Infinity;
        for (const r of recs) {
          const rs = view.start(r);
          if (rs < lo) lo = rs;
          const ops = view.ops(r);
          let e = rs;
          for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) e += ops[k] >>> 4; }
          if (e > hi) hi = e;
        }
        ref = hi > lo ? refWindow(lo, await this.getReferenceSeq(chrom, lo, hi)) : null;
      };
      const count = (layer: AlleleLayer, r: any) => {
        const flags = view.flags(r);
        if (!keepFlags(flags)) return;
        // every read is counted, unique or not ("unique only" is read back from the counts, not scanned again)
        const unique = this.isUnique(sampleId, view, r);
        const ops = view.ops(r);
        const n = view.seqCodes(r, scratch);
        const start0 = view.start(r);
        layer.add(start0, ops, scratch.a, n, view.quals(r), ref, unique, (flags & 16) !== 0, view.mapq(r));
        if (st.longReads == null) {
          let span = 0;
          for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) span += ops[k] >>> 4; }
          st.spans.push(span);
        }
      };
      const span = Math.max(1, end - start);
      const progress = () => {
        const done = Math.max(0, Math.min(end, st.pe) - Math.max(start, st.ps));
        opts?.onProgress?.(Math.min(1, done / span));
        if (st.longReads == null && st.spans.length) { const sp = [...st.spans].sort((x, y) => x - y); st.longReads = sp[sp.length >> 1] > 1000; st.spans = []; }
      };
      await this.fillStretch(sampleId, st, loc, start, end, () => new AlleleLayer(), count, signal, progress, prepare);
      progress();
      const longReads = st.longReads ?? false;
      const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
      const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
      const siteRef = refWindow(start, await this.getReferenceSeq(chrom, start, end));
      throwIfAborted(signal);
      const { sites, reads } = sitesFromCounts([st.owned, st.spill], uniqueOnly, start, end, siteRef, vaf, minIndel);
      return { sites, total: reads, long_reads: longReads };
    });
  }

  async getLibraryType(sampleId: number): Promise<LibraryEvidence> {
    if (!this.samples.has(sampleId)) return { type: 'unknown', source: 'none', note: 'sample not found' };
    return classifyHeader(await this.headerText(sampleId));
  }
  getTranscript(geneName: string, geneId?: string, hint?: RegionHint): Promise<TranscriptData> { return getTranscript(this.reference.build, geneName, geneId, hint); }
  getAllTranscripts(geneName: string, geneId?: string, hint?: RegionHint): Promise<AllTranscripts> { return getAllTranscripts(this.reference.build, geneName, geneId, hint); }
  async getRunSamples(): Promise<SampleRef[]> { return this.list(); }
  getRegionGenes(chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]> { return getRegionGenes(this.reference.build, chrom, start, end, exclude); }
  getProteinDomains(model: ProteinModelRef): Promise<ProteinDomain[]> { return getProteinDomains(this.reference.build, model); }
  getCommonSnps(chrom: string, start: number, end: number) { return getCommonSnps(this.reference.build, chrom, start, end); }
  getGtexTissues(): Promise<GtexTissue[]> {
    if (this.reference.build !== 'GRCh38') return Promise.reject(new Error('GTEx tissue tracks need GRCh38 (GTEx v8/v10 are hg38)'));
    return getGtexTissues();
  }
  getGtexProfile(geneName: string, geneId: string | undefined, tissue: GtexTissue): Promise<GtexProfile> { return getGtexProfile(geneName, geneId, tissue); }
  async getReference(chrom: string, start: number, end: number): Promise<string | null> { const s = await this.getReferenceSeq(chrom, start, end); return s == null ? null : s.toUpperCase(); }
  async getRandomSample(_runId: number, excludeId: number): Promise<SampleRef> {
    const other = this.list().find(s => s.id !== excludeId);
    if (!other) throw new Error('No other sample');
    return other;
  }

  private strandCalls = new Map<number, StrandnessCall>();

  /** Every loaded file is the "run": per-exon depth of each of them, strand-aware when the library is detected as stranded. */
  async getExonUsage(_runId: number, chrom: string, strand: number, exons: [number, number][], uniqueOnly: boolean): Promise<ExonUsageResponse> {
    const geneStrand = strand >= 0 ? 1 : -1;
    const samples: SampleExonDepths[] = [];
    for (const s of this.list()) {
      try {
        const perExon: Scan[] = [];
        for (const [a, b] of exons) perExon.push(await this.scan(s.id, chrom, a, b, uniqueOnly, EXON_USAGE_MAX_READS, true));
        let strandness = this.strandCalls.get(s.id), fraction: number | null = null;
        if (!strandness) {
          const call = detectStrandness(perExon.flatMap(x => x.kept), geneStrand);
          strandness = call.strandness; fraction = call.fraction;
          if (strandness === 'firststrand' || strandness === 'secondstrand') this.strandCalls.set(s.id, strandness);
        }
        const keep = strandKeeper(strandness, geneStrand);
        samples.push({
          sample_id: s.id, sample_name: s.name, strandness, strand_fraction: fraction,
          exons: perExon.map((sc, i) => {
            const d = exonDepth(sc.kept.filter(keep).map(r => encodeRead(r, null, 0)), exons[i][0], exons[i][1]);
            return sc.rate === 1 ? d : { median: d.median * sc.rate, mean: d.mean * sc.rate, reads: d.reads * sc.rate };
          }),
        });
      } catch (e: any) {
        samples.push({ sample_id: s.id, sample_name: s.name, strandness: 'unknown', strand_fraction: null, exons: [], error: e?.message || String(e) });
      }
    }
    return { run_id: 0, chrom, exons, samples };
  }

  // ---------------- Coverage: exact, streamed, kept ----------------
  // Coverage requests count every read of their window straight from its packed CIGAR into a per-sample
  // CoverageState (see coverage.ts): no read object, no sampling, and what was counted once is not decoded
  // again — panning back, a later request inside the same stretch, another gene model's boundaries, or
  // "unique only" toggled are answered from the counts.

  /** What each sample has counted, one chromosome stretch per sample. */
  private coverage = new Map<number, CoverageState>();
  /** Requests of one sample run one after the other: they grow the same state. */
  private coverageLocks = new Map<number, Promise<unknown>>();
  /** Whether a file's records carry NH (uniqueness from it) or not (from MAPQ), or how many were probed so far. */
  private nhMode = new Map<number, 'tag' | 'mapq' | number>();

  private locked<T>(locks: Map<number, Promise<unknown>>, id: number, job: () => Promise<T>): Promise<T> {
    const run = (locks.get(id) ?? Promise.resolve()).catch(() => undefined).then(job);
    locks.set(id, run.catch(() => undefined));
    return run;
  }

  /** Uniquely mapped, as uniqueFrom decides; NH is not looked for in files that do not carry it (a tag lookup per record costs). */
  private isUnique(id: number, view: RecordView<any>, r: any): boolean {
    const mode = this.nhMode.get(id);
    if (mode === 'mapq') return view.mapq(r) >= 30;
    const nh = view.nh(r);
    if (mode !== 'tag') {
      if (nh != null) this.nhMode.set(id, 'tag');
      else { const seen = (typeof mode === 'number' ? mode : 0) + 1; this.nhMode.set(id, seen >= NH_PROBE ? 'mapq' : seen); }
    }
    return uniqueFrom(nh, view.mapq(r));
  }

  /**
   * Counts one record into a layer. With `structural`, the records the structural evidence reads
   * (SA tag, a clip of RESCUE_MIN_CLIP bases or more, a deletion of SV_MIN_DELETION or more, a mate
   * elsewhere, on the same strand or far away) are kept whole; the others would add nothing to it
   * but their insert size, which the layer histograms instead.
   */
  private countRecord(id: number, layer: Layer, view: RecordView<any>, r: any, seqId: number, refNames: string[], structural: boolean): void {
    const flags = view.flags(r);
    if (!keepFlags(flags)) return;
    const unique = this.isUnique(id, view, r);
    const ops = view.ops(r);
    const start = view.start(r);
    const insert = (flags & 1) && (flags & 2) ? Math.abs(view.tlen(r)) : 0;
    layer.add(start, ops, unique, insert);
    if (!structural) return;
    let keep = false, refLen = 0;
    for (let k = 0; k < ops.length; k++) {
      const len = ops[k] >>> 4, op = ops[k] & 15;
      if ((op === 4 || op === 5) && len >= RESCUE_MIN_CLIP) keep = true;
      else if (op === 2 && len >= SV_MIN_DELETION) keep = true;
      if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) refLen += len;
    }
    if (!keep && (flags & 1) && !(flags & 8)) {
      const mate = view.mateRef(r);
      keep = mate >= 0 && (mate !== seqId || ((flags & 16) !== 0) === ((flags & 32) !== 0) || Math.abs(view.tlen(r)) > 1000);
    }
    if (!keep) keep = typeof view.sa(r) === 'string';
    if (keep) layer.sv.push({ r: view.raw(r, true, true, refNames), unique, end: start + refLen });
  }

  /**
   * Grows a sample's counted stretch to cover [from, to), tile by tile, each tile's records counted in
   * one synchronous step after they arrive: an abort between tiles leaves the state whole. Growing right
   * adds the reads starting in the new tiles; growing left also rebuilds the spill (see CoverageState).
   */
  private fill(id: number, st: CoverageState, loc: { o: Opened; name: string; seqId: number }, from: number, to: number, signal?: AbortSignal, onTile?: () => void | Promise<void>): Promise<void> {
    const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
    return this.fillStretch(id, st, loc, from, to, () => new Layer(),
      (layer, r) => this.countRecord(id, layer, view, r, loc.seqId, loc.o.refNames, st.structural), signal, onTile);
  }

  /**
   * Grows a counted stretch (coverage, or allele counts) to cover [from, to), tile by tile. Each tile's records are
   * counted in one synchronous step after they arrive (and after `prepare`, when given, has fetched what the counting
   * needs), so an abort between tiles leaves the state whole. Growing right adds the reads starting in the new tiles;
   * growing left also rebuilds the spill (see CoverageState).
   */
  private async fillStretch<L>(id: number, st: { owned: L; spill: L; ps: number; pe: number; readonly empty: boolean }, loc: { o: Opened; name: string; seqId: number },
    from: number, to: number, newLayer: () => L, count: (layer: L, r: any) => void, signal?: AbortSignal, onTile?: () => void | Promise<void>,
    prepare?: (recs: any[]) => Promise<void>): Promise<void> {
    const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
    const records = async (a: number, b: number): Promise<any[]> => {
      const recs: any[] = loc.o.kind === 'bam'
        ? await loc.o.bam.getRecordsForRange(loc.name, a, b, { signal })
        : await loc.o.cram.getRecordsForRange(loc.seqId, a, b, { signal });
      if (prepare) { await prepare(recs); throwIfAborted(signal); }
      return recs;
    };
    let seen = 0;
    const calibrate = (bytes: number) => { if (seen >= 1000 && bytes > 0) this.bytesPerRead.set(id, bytes / seen); };
    if (st.empty) { st.ps = from; st.pe = from; }
    if (to > st.pe) {
      const a0 = st.pe, bytes = await this.indexBytes(loc, a0, to, signal);
      const tile = this.tileSize(id, a0, to, bytes);
      seen = 0;
      while (st.pe < to) {
        throwIfAborted(signal);
        const a = st.pe, b = Math.min(to, a + tile);
        const recs = await records(a, b);
        const first = st.pe === st.ps;   // nothing owned yet: the reads reaching in from the left are the spill
        for (const r of recs) {
          const rs = view.start(r);
          if (rs < a) { if (first) count(st.spill, r); continue; }
          seen++;
          count(st.owned, r);
        }
        st.pe = b;
        if (b < to) { await onTile?.(); await yieldToUi(); }
      }
      calibrate(bytes);
    }
    if (from < st.ps) {
      const b0 = st.ps, bytes = await this.indexBytes(loc, from, b0, signal);
      const tile = this.tileSize(id, from, b0, bytes);
      seen = 0;
      while (st.ps > from) {
        throwIfAborted(signal);
        const b = st.ps, a = Math.max(from, b - tile);
        const recs = await records(a, b);
        const spill = newLayer();
        for (const r of recs) {
          const rs = view.start(r);
          if (rs >= a) seen++;
          count(rs < a ? spill : st.owned, r);
        }
        st.spill = spill;
        st.ps = a;
        if (a > from) { await onTile?.(); await yieldToUi(); }
      }
      calibrate(bytes);
    }
  }

  /** The sample's state for a request on `chrom` around [start, end): kept when the request is near it, started over otherwise. */
  private coverageState(id: number, chrom: string, start: number, end: number, structural: boolean): CoverageState {
    let st = this.coverage.get(id);
    const gap = Math.max(end - start, COVERAGE_MIN_GAP);
    if (!st || st.chrom !== chrom || (structural && !st.structural)
      || (!st.empty && (start > st.pe + gap || end < st.ps - gap))
      || (!st.empty && Math.max(st.pe, end) - Math.min(st.ps, start) > COVERAGE_MAX_SPAN)) {
      st = new CoverageState(chrom, structural);
      this.coverage.set(id, st);
    }
    st.lastUsed = Date.now();
    return st;
  }

  /** Drops the least recently used samples' counts while all of them together are over COVERAGE_CACHE_BYTES. */
  private trimCoverage(keep: number): void {
    let total = 0;
    for (const st of this.coverage.values()) total += st.bytes;
    const others = [...this.coverage.entries()].filter(([id]) => id !== keep).sort((a, b) => a[1].lastUsed - b[1].lastUsed);
    for (const [id, st] of others) { if (total <= COVERAGE_CACHE_BYTES) break; total -= st.bytes; this.coverage.delete(id); }
  }

  /**
   * Exact coverage, junctions and boundary-spanning counts of [start, end). The core (the view) is
   * counted first, in full whatever its depth, and handed to `opts.onProgress` as it fills; the
   * margins follow within `opts.maxReads` (an index estimate of what they cost) and the returned
   * `window` says how far they went. Everything counted stays with the sample for later requests.
   */
  async getCoverage(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint, opts?: CoverageOptions): Promise<SampleCoverage> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    if (end - start > MAX_REGION_BP) throw new Error(`Region too large (${(end - start).toLocaleString()} bp); maximum is ${MAX_REGION_BP.toLocaleString()} bp`);
    const signal = opts?.signal;
    const loc = await this.locate(sampleId, chrom);
    const result = (w: { start: number; end: number }, sl: CoverageSlice): SampleCoverage => ({
      sample_id: sampleId, sample_name: s.name, coverage: sl.coverage, junctions: sl.junctions, spanning: sl.spanning, window: w, spliced: sl.spliced,
    });
    if (!loc) return result({ start, end }, await readSlice([], uniqueOnly, start, end, boundaries));
    return this.locked(this.coverageLocks, sampleId, async () => {
      throwIfAborted(signal);
      const core = { start: Math.max(start, Math.min(end, opts?.core?.start ?? start)), end: Math.min(end, Math.max(start, opts?.core?.end ?? end)) };
      if (core.end <= core.start) { core.start = start; core.end = end; }
      const structural = !!opts?.structural;
      const st = this.coverageState(sampleId, loc.name, start, end, structural);
      const exact = () => ({ start: Math.max(start, st.ps), end: Math.min(end, st.pe) });
      // a partial answer costs a slice of everything counted so far, which on a wide, deep window is ~10⁶ runs: the next
      // one waits at least PROGRESS_MS and 4× what the last one took, so the partials never take over the reading
      let next = 0;
      const progress = async (force: boolean) => {
        if (!opts?.onProgress) return;
        const t0 = performance.now();
        if (!force && t0 < next) return;
        const w = exact();
        if (w.end <= w.start) return;
        const sl = await st.slice(uniqueOnly, w.start, w.end, boundaries, yieldToUi);
        throwIfAborted(signal);
        const t1 = performance.now();
        next = t1 + Math.max(PROGRESS_MS, 4 * (t1 - t0));
        opts.onProgress(result(w, sl));
      };
      // 1. the view, whatever its depth
      await this.fill(sampleId, st, loc, core.start, core.end, signal, () => progress(false));
      // 2. the margins, within their budget; what was already counted costs nothing and is not trimmed
      if (!st.covers(start, end)) {
        const cap = Math.max(1000, opts?.maxReads ?? DEFAULT_MARGIN_READS);
        const win = await this.budgetWindow(sampleId, chrom, start, end, core, cap, signal);
        if (!st.covers(win.start, win.end)) {
          await progress(true);
          await this.fill(sampleId, st, loc, win.start, win.end, signal, () => progress(false));
        }
      }
      this.trimCoverage(sampleId);
      const w = exact();
      const sl = await st.slice(uniqueOnly, w.start, w.end, boundaries, yieldToUi);
      throwIfAborted(signal);
      const out = result(w, sl);
      if (structural) {
        // the reference of the window lets clip clusters be placed by realignment and clipped reads be rescued at the
        // breakpoints seen: fetched when the window has something to place or rescue, always with a FASTA, up to
        // REALIGN_MAX_BP through the web APIs
        let ev = structuralEvidence(sl.sv, loc.name, w.start, w.end, 1, null, sl.insertMedian);
        if (this.reference.fasta || w.end - w.start <= REALIGN_MAX_BP) {
          const hasArcs = ev.splits.length + ev.deletions.length + (ev.duplications?.length ?? 0) + (ev.inversions?.length ?? 0) > 0;
          if (hasRealignableClips(sl.sv) || (hasArcs && hasRescuableClips(sl.sv))) {
            try {
              const seq = await this.getReferenceSeq(chrom, w.start, w.end);
              if (seq) ev = structuralEvidence(sl.sv, loc.name, w.start, w.end, 1, { start: w.start, seq }, sl.insertMedian);
            } catch (e) { console.warn('reference for clip realignment not available:', e); }
            // outside the catch: a reference that cannot be fetched is not fatal, a caller that gave up is
            throwIfAborted(signal);
          }
        }
        // only the records carrying evidence went in; the reads the evidence was gathered from are all those counted
        out.structural = { ...ev, reads: sl.spliced.reads };
      }
      return out;
    });
  }

  async getReads(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, maxReads: number,
    mode: 'reads' | 'collapsed', minSupport: number, minVaf: number, opts?: ReadsOptions): Promise<ReadsResponse> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    if (end - start > MAX_READS_REGION_BP) throw new Error(`Region too large for reads (${(end - start).toLocaleString()} bp)`);
    const collapsed = mode === 'collapsed';
    const cap = collapsed ? Math.max(40000, maxReads) : Math.max(100, maxReads);
    // filtered and sampled before names, sequences and qualities are decoded: only the reads shown pay for them; the
    // pair fields (mate position, template length) come along so that mates can be drawn linked
    const { total, kept: raw } = await this.scan(sampleId, chrom, start, end, uniqueOnly, cap, false, { structural: true, signal: opts?.signal });
    const refStart = Math.max(0, start - 500);
    const ref = await this.getReferenceSeq(chrom, refStart, end + 500);
    throwIfAborted(opts?.signal);
    const own = (await this.locate(sampleId, chrom))?.name ?? chrom;
    const reads: AlignedRead[] = raw.map(r => withMate(encodeRead(r, ref, refStart), r, own));
    if (opts?.methylation && !collapsed && ref) readMethylation(reads, raw, cpgSites(refStart, ref));
    const longReads = isLongRead(reads);
    const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
    const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
    const base = { sample_id: sampleId, sample_name: s.name, total, shown: reads.length, long_reads: longReads,
      reference: ref != null ? { start: refStart, seq: ref } : null, reference_source: ref != null ? this.lastReferenceSource : null, haplotags: haplotagCounts(reads) };
    if (collapsed) {
      if (opts?.haplotypes !== 'any') {
        const phaseOf = () => phaseReads(reads, start, end, ref, refStart, 3, vaf, 20, minIndel);
        const { phase, haplotypes } = windowHaplotypes(reads, start, end, ref, refStart, vaf, minIndel, opts?.phaseSource ?? 'auto', phaseOf);
        return { ...base, reads: [], sites: phase?.sites ?? callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel), groups: [], phase, haplotypes };
      }
      const summary = collapseReads(reads, start, end, ref, refStart, 3, vaf, 20, Math.max(1, minSupport), minIndel);
      return { ...base, reads: [], sites: summary.sites, groups: summary.groups };
    }
    return { ...base, reads, sites: callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel), groups: [] };
  }
}
