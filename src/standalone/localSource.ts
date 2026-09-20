/**
 * SashimiDataSource over local BAM/CRAM files, decoded in the browser with the GMOD
 * libraries. Nothing from the alignment files leaves the machine; only gene lookups and
 * reference-sequence requests go to the UCSC API (Ensembl REST as fallback), unless a local
 * FASTA is given.
 */
import { BamFile } from '@gmod/bam';
import { IndexedCramFile, CraiIndex } from '@gmod/cram';
import { IndexedFasta, BgzipIndexedFasta } from '@gmod/indexedfasta';
import { BlobFile } from 'generic-filehandle2';
import { unzip } from '@gmod/bgzf-filehandle';
import type { AlignedRead, AllTranscripts, BoundarySpanning, ExonUsageResponse, GeneModel, GtexProfile, GtexTissue, KnownVariant, LibraryEvidence, ProteinDomain, ProteinModelRef, ReadsResponse, RegionHint, SampleCoverage, BoundaryHint, SampleExonDepths, TranscriptData, VariantSite } from '../components/sashimi/types';
import type { CoverageOptions, ReadsOptions, SashimiDataSource, SampleRef, VariantScan, VariantScanOptions } from '../components/sashimi/datasource';
import { SV_MIN_CLIP, boundarySpanning, coverageRuns, cramCigar, cramMismatches, detectStrandness, encodeRead, exonDepth, hasRealignableClips, junctionCounts, keepFlags, strandKeeper, structuralEvidence, uniqueFrom, type RawRead, type StrandnessCall } from './alignments';
import { callSites, collapseReads } from './collapse';
import { phaseReads } from './phasing';
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
/** Tile of the full variant scan: the reads of one tile are decoded, called and dropped before the next, so memory stays bounded on any window. */
const VARIANT_TILE_BP = 100_000;

// ---------------- Deep regions ----------------
// A very deep library (targeted RNA-seq, a highly expressed gene) can hold millions of records over one gene,
// and decoding them all froze the page. Requests are therefore budgeted: the window is sized from the index
// before anything is decoded, records are scanned tile by tile with only the fields the filter needs, and past
// a cap every k-th read is kept (k = 2, 4, 8…) with the counts scaled back by k.
/** Reads decoded per coverage request unless the caller says otherwise. */
const DEFAULT_MAX_READS = 250_000;
/** Reads decoded per exon for the exon-usage statistics (fractions and medians only need a sample). */
const EXON_USAGE_MAX_READS = 100_000;
/** A window is scanned in tiles so that one tile's decoded records can be released before the next is read. */
/** widest window whose reference is fetched from the web APIs to place clipped sequences (a local FASTA has no limit) */
const REALIGN_MAX_BP = 500_000;
const TILE_BP = 250_000;
/** Compressed bytes per record assumed before a file has been scanned once (a scan then calibrates it). */
const BYTES_PER_READ: Record<'bam' | 'cram', number> = { bam: 60, cram: 30 };
/** Margins are halved while the window looks too deep; below this they are dropped altogether. */
const MIN_MARGIN_BP = 2_000;
/** Decoded BAM chunks the library keeps in memory (its default is 1 GB, too much for a browser tab). */
const BAM_CACHE_BYTES = 256 * 1024 * 1024;

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
  /** `light` leaves out name, sequence and qualities (coverage only needs the alignment blocks); `structural` adds the pair and SA fields. */
  raw(r: R, light: boolean, structural: boolean, refNames: string[]): RawRead;
}
const mateFields = (chromOf: (id: number) => string, mateId: number, matePos: number, tlen: number, sa: unknown) => ({
  tlen, mateChrom: mateId >= 0 ? chromOf(mateId) : '', matePos: mateId >= 0 ? matePos : undefined, sa: typeof sa === 'string' ? sa : null,
});
/** A CIGAR with a soft clip of SV_MIN_CLIP bases or more: the light structural scan decodes this record's sequence so the clip can be placed by realignment. */
const CLIP_RE = new RegExp(`(?:^|[A-Z=])(\\d+)S`, 'g');
const bigClip = (cigar: string) => { CLIP_RE.lastIndex = 0; let m: RegExpExecArray | null; while ((m = CLIP_RE.exec(cigar))) if (parseInt(m[1]) >= SV_MIN_CLIP) return true; return false; };
const BAM_VIEW: RecordView<any> = {
  start: r => r.start, flags: r => r.flags, mapq: r => r.mq ?? 255, nh: r => tagNumber(r.getTag('NH')),
  raw: (r, light, structural, refNames) => {
    const sa = structural ? r.getTag('SA') : undefined;
    const withSeq = !light || (structural && typeof sa !== 'string' && bigClip(r.CIGAR));
    // the name is what ties the parts of a split read together: kept for reads with an SA tag even in the light scan
    return { name: light && typeof sa !== 'string' ? '' : r.name, start: r.start, cigar: r.CIGAR, seq: withSeq ? r.seq : '', qual: light ? null : r.qual, flags: r.flags, mapq: r.mq ?? 255, nh: tagNumber(r.getTag('NH')),
      ...(structural ? mateFields(id => refNames[id] ?? '', r.next_refid, r.next_pos, r.template_length, sa) : {}) };
  },
};
const CRAM_VIEW: RecordView<any> = {
  start: r => r.start, flags: r => r.flags, mapq: r => r.mappingQuality ?? 255, nh: r => tagNumber(r.getTag('NH')),
  raw: (r, light, structural, refNames) => {
    const feats = r.readFeatures as any;
    const qual = r.qualityScores ?? null;
    const cigar = cramCigar(feats, r.readLength, r.lengthOnRef ?? 0);
    const sa = structural ? r.getTag('SA') : undefined;
    const withSeq = !light || (structural && typeof sa !== 'string' && bigClip(cigar));
    return { name: light && typeof sa !== 'string' ? '' : (r.readName ?? ''), start: r.start, cigar, seq: withSeq ? (r.readBases ?? '') : '', qual: light ? null : qual, flags: r.flags, mapq: r.mappingQuality ?? 255, nh: tagNumber(r.getTag('NH')),
      mismatches: light ? undefined : cramMismatches(feats, qual),
      ...(structural ? mateFields(id => refNames[id] ?? '', r.nextSequenceId ?? -1, (r.nextStart ?? 0) - 1, r.templateLength ?? r.templateSize ?? 0, sa) : {}) };
  },
};

/** Result of a budgeted scan: `kept` holds every `rate`-th of the `total` reads that passed the filters. */
interface Scan { total: number; rate: number; kept: RawRead[] }
const scaleRuns = <T extends { depth: number }>(runs: T[], k: number): T[] => (k === 1 ? runs : runs.map(r => ({ ...r, depth: r.depth * k })));
const scaleCounts = <T extends { count: number }>(xs: T[], k: number): T[] => (k === 1 ? xs : xs.map(x => ({ ...x, count: x.count * k })));
const scaleRecord = (o: Record<number, number>, k: number): Record<number, number> => (k === 1 ? o : Object.fromEntries(Object.entries(o).map(([p, n]) => [p, n * k])));
const scaleSpanning = (sp: BoundarySpanning, k: number): BoundarySpanning => ({ intronStart: scaleRecord(sp.intronStart, k), intronEnd: scaleRecord(sp.intronEnd, k) });

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
    // CRAM decoding depends on the reference: reopen files
    for (const [id, s] of this.samples) if (s.kind === 'cram') this.opened.delete(id);
  }

  /** Variants handed over by the page URL (deep link), drawn on every sample. */
  knownVariants: KnownVariant[] = [];
  async getKnownVariants(_sampleId: number): Promise<KnownVariant[]> { return this.knownVariants; }

  addSample(s: LocalSample) { this.samples.set(s.id, s); }
  renameSample(id: number, name: string) { const s = this.samples.get(id); if (s) this.samples.set(id, { ...s, name }); }
  removeSample(id: number) { this.samples.delete(id); this.opened.delete(id); this.headers.delete(id); }
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
          const bam = new BamFile({ bamFilehandle: new BlobFile(s.file), baiFilehandle: new BlobFile(s.index), maxCacheBytes: BAM_CACHE_BYTES });
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

  /** Compressed bytes the index says [start, end) occupies: BAI chunks for BAM, CRAI slices for CRAM. */
  private async indexBytes(loc: { o: Opened; name: string; seqId: number }, start: number, end: number): Promise<number> {
    if (loc.o.kind === 'bam') {
      const chunks = await loc.o.bam.blocksForRange(loc.name, start, end);
      return chunks.reduce((a, c) => a + c.fetchedSize(), 0);
    }
    const slices = await loc.o.cram.index.getEntriesForRange(loc.seqId, start, end);
    return slices.reduce((a, sl) => a + sl.sliceBytes, 0);
  }

  /**
   * Reads of [start, end) passing the flag (and uniqueness) filters, decoded tile by tile, at most about `cap`
   * of them: when the kept reads outgrow the cap they are thinned to every other one and the rate doubles, so
   * `kept` is always the reads whose rank (among the passing reads, in file order) is a multiple of `rate` — a
   * systematic sample that is exact (rate 1) whenever the region holds no more than `cap` reads. A read spanning
   * two tiles is counted in the tile holding its start (or the first tile when it starts before the window).
   */
  private async scan(id: number, chrom: string, start: number, end: number, uniqueOnly: boolean, cap: number, light: boolean, structural = false): Promise<Scan> {
    const loc = await this.locate(id, chrom);
    if (!loc) return { total: 0, rate: 1, kept: [] };
    const view: RecordView<any> = loc.o.kind === 'bam' ? BAM_VIEW : CRAM_VIEW;
    let kept: RawRead[] = [], total = 0, rate = 1, seen = 0;
    const bytes = await this.indexBytes(loc, start, end);
    for (let ts = start; ts < end; ts += TILE_BP) {
      const te = Math.min(end, ts + TILE_BP);
      const recs: any[] = loc.o.kind === 'bam' ? await loc.o.bam.getRecordsForRange(loc.name, ts, te) : await loc.o.cram.getRecordsForRange(loc.seqId, ts, te);
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
    }
    if (seen >= 1000 && bytes > 0) this.bytesPerRead.set(id, bytes / seen);
    return { total, rate, kept };
  }

  /**
   * The window to read for a request: the core is always read; the margins around it are halved while the
   * index suggests more reads than the budget, and dropped once they get small.
   */
  private async budgetWindow(id: number, chrom: string, start: number, end: number, core: { start: number; end: number }, cap: number): Promise<{ start: number; end: number }> {
    const loc = await this.locate(id, chrom);
    const s = this.samples.get(id);
    if (!loc || !s) return { start, end };
    const bpr = this.bytesPerRead.get(id) ?? BYTES_PER_READ[s.kind];
    const coreStart = Math.max(start, Math.min(end, core.start)), coreEnd = Math.max(coreStart, Math.min(end, core.end));
    let left = coreStart - start, right = end - coreEnd;
    let win = { start, end };
    while ((left > 0 || right > 0) && (await this.indexBytes(loc, win.start, win.end)) / bpr > cap) {
      left = left >= 2 * MIN_MARGIN_BP ? Math.floor(left / 2) : 0;
      right = right >= 2 * MIN_MARGIN_BP ? Math.floor(right / 2) : 0;
      win = { start: coreStart - left, end: coreEnd + right };
    }
    return win;
  }

  // ---- SashimiDataSource ----
  async getPrimaryRecord(sampleId: number, chrom: string, start: number, name: string): Promise<{ seq: string; flags: number; cigar: string } | null> {
    const { kept } = await this.scan(sampleId, chrom, start, start + 1, false, Number.MAX_SAFE_INTEGER, false, true);
    const r = kept.find(x => x.start === start && x.name === name && !(x.flags & 2048) && x.seq);
    return r ? { seq: r.seq, flags: r.flags, cigar: r.cigar } : null;
  }
  /** Every site above the thresholds from every read of the window, one tile at a time (no read cap, any window width). */
  async getVariantSites(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    const sites: VariantSite[] = [];
    let total = 0, longReads: boolean | null = null;
    for (let ts = start; ts < end; ts += VARIANT_TILE_BP) {
      if (opts?.signal?.aborted) throw new DOMException('Variant scan cancelled', 'AbortError');
      const te = Math.min(end, ts + VARIANT_TILE_BP);
      // every read overlapping the tile (those starting before it too: they cover its first positions), decoded in full
      const { kept } = await this.scan(sampleId, chrom, ts, te, uniqueOnly, Number.MAX_SAFE_INTEGER, false);
      for (const r of kept) if (ts === start || r.start >= ts) total++;   // a read spanning two tiles is counted once
      if (kept.length) {
        const refStart = Math.max(0, ts - 500);
        const ref = await this.getReferenceSeq(chrom, refStart, te + 500);
        const reads = kept.map(r => encodeRead(r, ref, refStart));
        if (longReads == null) longReads = isLongRead(reads);   // decided on the first tile holding reads, kept for the whole window
        const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
        const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
        sites.push(...callSites(reads, ts, te, ref, refStart, 3, vaf, 20, minIndel));
      }
      opts?.onProgress?.((te - start) / (end - start));
    }
    return { sites, total, long_reads: longReads ?? false };
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

  async getCoverage(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint, opts?: CoverageOptions): Promise<SampleCoverage> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    if (end - start > MAX_REGION_BP) throw new Error(`Region too large (${(end - start).toLocaleString()} bp); maximum is ${MAX_REGION_BP.toLocaleString()} bp`);
    const cap = Math.max(1000, opts?.maxReads ?? DEFAULT_MAX_READS);
    const win = await this.budgetWindow(sampleId, chrom, start, end, opts?.core ?? { start, end }, cap);
    const { total, rate, kept } = await this.scan(sampleId, chrom, win.start, win.end, uniqueOnly, cap, true, !!opts?.structural);
    const reads = kept.map(r => encodeRead(r, null, 0));
    const loc = await this.locate(sampleId, chrom);
    // the reference of the window lets clip clusters be placed by realignment: always with a FASTA, up to REALIGN_MAX_BP through the web APIs
    let ref: { start: number; seq: string } | null = null;
    if (opts?.structural && hasRealignableClips(kept) && (this.reference.fasta || win.end - win.start <= REALIGN_MAX_BP)) {
      try { const seq = await this.getReferenceSeq(chrom, win.start, win.end); if (seq) ref = { start: win.start, seq }; } catch (e) { console.warn('reference for clip realignment not available:', e); }
    }
    const structural = opts?.structural ? structuralEvidence(kept, loc?.name ?? chrom, win.start, win.end, rate, ref) : undefined;
    const splicedReads = kept.reduce((n, r) => n + (/\d+N/.test(r.cigar) ? 1 : 0), 0);
    const junctions = junctionCounts(reads, win.start, win.end);
    // unspliced reads through every splice site seen in the reads, plus the boundaries the caller asked for (annotated exons)
    const spanning = boundarySpanning(reads,
      [...junctions.map(j => j.start), ...(boundaries?.intronStarts ?? [])].filter(p => p >= win.start && p < win.end),
      [...junctions.map(j => j.end), ...(boundaries?.intronEnds ?? [])].filter(p => p > win.start && p <= win.end));
    return {
      sample_id: sampleId, sample_name: s.name,
      coverage: scaleRuns(coverageRuns(reads, win.start, win.end), rate), junctions: scaleCounts(junctions, rate), spanning: scaleSpanning(spanning, rate),
      window: win, sampled: rate > 1 ? { rate, total, decoded: kept.length } : undefined,
      spliced: { reads: kept.length, fraction: kept.length ? splicedReads / kept.length : 0 },
      structural,
    };
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
    const { total, kept: raw } = await this.scan(sampleId, chrom, start, end, uniqueOnly, cap, false, true);
    const refStart = Math.max(0, start - 500);
    const ref = await this.getReferenceSeq(chrom, refStart, end + 500);
    const own = (await this.locate(sampleId, chrom))?.name ?? chrom;
    const reads: AlignedRead[] = raw.map(r => withMate(encodeRead(r, ref, refStart), r, own));
    const longReads = isLongRead(reads);
    const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
    const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
    const base = { sample_id: sampleId, sample_name: s.name, total, shown: reads.length, long_reads: longReads,
      reference: ref != null ? { start: refStart, seq: ref } : null, reference_source: ref != null ? this.lastReferenceSource : null };
    if (collapsed) {
      if (opts?.haplotypes !== 'any') {
        const phase = phaseReads(reads, start, end, ref, refStart, 3, vaf, 20, minIndel);
        return { ...base, reads: [], sites: phase.sites, groups: [], phase };
      }
      const summary = collapseReads(reads, start, end, ref, refStart, 3, vaf, 20, Math.max(1, minSupport), minIndel);
      return { ...base, reads: [], sites: summary.sites, groups: summary.groups };
    }
    return { ...base, reads, sites: callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel), groups: [] };
  }
}
