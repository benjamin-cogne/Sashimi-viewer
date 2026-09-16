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
import type { AlignedRead, AllTranscripts, ExonUsageResponse, GeneModel, GtexProfile, GtexTissue, KnownVariant, ProteinDomain, ProteinModelRef, ReadsResponse, RegionHint, SampleCoverage, BoundaryHint, SampleExonDepths, TranscriptData } from '../components/sashimi/types';
import type { SashimiDataSource, SampleRef } from '../components/sashimi/datasource';
import { boundarySpanning, coverageRuns, cramCigar, cramMismatches, detectStrandness, encodeRead, exonDepth, isUnique, junctionCounts, keepRead, strandKeeper, type RawRead, type StrandnessCall } from './alignments';
import { callSites, collapseReads } from './collapse';
import type { GenomeBuild } from './ensembl';
import { getAllTranscripts, getProteinDomains, getReference, getRegionGenes, getTranscript } from './ucsc';
import { getCommonSnps } from './snps';
import { getGtexProfile, getGtexTissues } from './gtex';

export interface LocalSample { id: number; name: string; kind: 'bam' | 'cram'; file: File; index: File; /** paths relative to the run folder, when the files came from one */ path?: string; indexPath?: string }
export interface ReferenceChoice { build: GenomeBuild; fasta?: { fa: File; fai: File; gzi?: File } }

type Opened =
  | { kind: 'bam'; bam: BamFile; refNames: string[] }
  | { kind: 'cram'; cram: IndexedCramFile; refNames: string[] };

const MAX_REGION_BP = 5_000_000;
const MAX_READS_REGION_BP = 250_000;

function resolveName(names: string[], chrom: string): string | null {
  if (names.includes(chrom)) return chrom;
  const alt = chrom.startsWith('chr') ? chrom.slice(3) : `chr${chrom}`;
  return names.includes(alt) ? alt : null;
}

function tagNumber(v: unknown): number | null {
  return typeof v === 'number' && Number.isFinite(v) ? v : null;
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
  removeSample(id: number) { this.samples.delete(id); this.opened.delete(id); }
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
  private open(id: number): Promise<Opened> {
    const s = this.samples.get(id);
    if (!s) return Promise.reject(new Error('Sample not found'));
    if (!this.opened.has(id)) {
      this.opened.set(id, (async (): Promise<Opened> => {
        if (s.kind === 'bam') {
          const bam = new BamFile({ bamFilehandle: new BlobFile(s.file), baiFilehandle: new BlobFile(s.index) });
          await bam.getHeader();
          const refNames = (bam.indexToChr || []).map(r => r.refName);
          return { kind: 'bam', bam, refNames };
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
        return { kind: 'cram', cram, refNames: info.map(r => r.name) };
      })().catch(e => { this.opened.delete(id); throw e; }));
    }
    return this.opened.get(id)!;
  }

  private async records(id: number, chrom: string, start: number, end: number): Promise<RawRead[]> {
    const o = await this.open(id);
    const name = resolveName(o.refNames, chrom);
    if (!name) return [];
    const out: RawRead[] = [];
    if (o.kind === 'bam') {
      const recs = await o.bam.getRecordsForRange(name, start, end);
      for (const r of recs) {
        out.push({ name: r.name, start: r.start, cigar: r.CIGAR, seq: r.seq, qual: r.qual, flags: r.flags, mapq: r.mq ?? 255, nh: tagNumber(r.getTag('NH')) });
      }
    } else {
      const seqId = o.refNames.indexOf(name);
      const recs = await o.cram.getRecordsForRange(seqId, start, end);
      for (const r of recs) {
        const feats = r.readFeatures as any;
        const qual = r.qualityScores ?? null;
        const cigar = cramCigar(feats, r.readLength, r.lengthOnRef ?? 0);
        out.push({ name: r.readName ?? '', start: r.start, cigar, seq: r.readBases ?? '', qual, flags: r.flags, mapq: r.mappingQuality ?? 255, nh: tagNumber(r.getTag('NH')),
          mismatches: cramMismatches(feats, qual) });
      }
    }
    return out;
  }

  private filtered(raw: RawRead[], uniqueOnly: boolean): RawRead[] {
    return raw.filter(r => keepRead(r) && (!uniqueOnly || isUnique(r)));
  }

  // ---- SashimiDataSource ----
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
        const perExon = await Promise.all(exons.map(([a, b]) => this.records(s.id, chrom, a, b).then(raw => this.filtered(raw, uniqueOnly))));
        let strandness = this.strandCalls.get(s.id), fraction: number | null = null;
        if (!strandness) {
          const call = detectStrandness(perExon.flat(), geneStrand);
          strandness = call.strandness; fraction = call.fraction;
          if (strandness === 'firststrand' || strandness === 'secondstrand') this.strandCalls.set(s.id, strandness);
        }
        const keep = strandKeeper(strandness, geneStrand);
        samples.push({
          sample_id: s.id, sample_name: s.name, strandness, strand_fraction: fraction,
          exons: perExon.map((raw, i) => exonDepth(raw.filter(keep).map(r => encodeRead(r, null, 0)), exons[i][0], exons[i][1])),
        });
      } catch (e: any) {
        samples.push({ sample_id: s.id, sample_name: s.name, strandness: 'unknown', strand_fraction: null, exons: [], error: e?.message || String(e) });
      }
    }
    return { run_id: 0, chrom, exons, samples };
  }

  async getCoverage(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint): Promise<SampleCoverage> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    if (end - start > MAX_REGION_BP) throw new Error(`Region too large (${(end - start).toLocaleString()} bp); maximum is ${MAX_REGION_BP.toLocaleString()} bp`);
    const raw = this.filtered(await this.records(sampleId, chrom, start, end), uniqueOnly);
    const reads = raw.map(r => encodeRead(r, null, 0));
    const junctions = junctionCounts(reads, start, end);
    // unspliced reads through every splice site seen in the reads, plus the boundaries the caller asked for (annotated exons)
    const spanning = boundarySpanning(reads,
      [...junctions.map(j => j.start), ...(boundaries?.intronStarts ?? [])].filter(p => p >= start && p < end),
      [...junctions.map(j => j.end), ...(boundaries?.intronEnds ?? [])].filter(p => p > start && p <= end));
    return { sample_id: sampleId, sample_name: s.name, coverage: coverageRuns(reads, start, end), junctions, spanning };
  }

  async getReads(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, maxReads: number,
    mode: 'reads' | 'collapsed', minSupport: number, minVaf: number): Promise<ReadsResponse> {
    const s = this.samples.get(sampleId);
    if (!s) throw new Error('Sample not found');
    if (end - start > MAX_READS_REGION_BP) throw new Error(`Region too large for reads (${(end - start).toLocaleString()} bp)`);
    const collapsed = mode === 'collapsed';
    let raw = this.filtered(await this.records(sampleId, chrom, start, end), uniqueOnly);
    const total = raw.length;
    const cap = collapsed ? 40000 : Math.max(100, Math.min(maxReads, 10000));
    if (total > cap) { const step = total / cap; raw = Array.from({ length: cap }, (_, i) => raw[Math.floor(i * step)]); }
    const refStart = Math.max(0, start - 500);
    const ref = await this.getReferenceSeq(chrom, refStart, end + 500);
    const reads: AlignedRead[] = raw.map(r => encodeRead(r, ref, refStart));
    const base = { sample_id: sampleId, sample_name: s.name, total, shown: reads.length,
      reference: ref != null ? { start: refStart, seq: ref } : null, reference_source: ref != null ? this.lastReferenceSource : null };
    if (collapsed) {
      const summary = collapseReads(reads, start, end, ref, refStart, 3, minVaf, 20, Math.max(1, minSupport));
      return { ...base, reads: [], sites: summary.sites, groups: summary.groups };
    }
    return { ...base, reads, sites: callSites(reads, start, end, ref, refStart, 3, minVaf, 20), groups: [] };
  }
}
