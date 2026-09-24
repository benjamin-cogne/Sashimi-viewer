/**
 * Exported viewer: one HTML file carrying the data of every registered view, for a reader who does
 * not have the alignment files.
 *
 * "Export HTML" writes a copy of this very page with a JSON payload in a `<script type="application/json">`
 * tag: per view, the gene models (reference transcript, every transcript, neighbouring genes) and, per
 * sample, the coverage, junctions and boundary-spanning counts of the window already fetched around the
 * view (view plus margins). When the page starts with such a payload, an EmbeddedDataSource serves those
 * regions from memory: the reader switches views, zooms and pans inside the exported windows; anything
 * else (another gene, the reads track, exon-usage statistics) needs the original files, which the reader
 * can still add. Gene lookups, common SNPs and GTEx go to the network as usual when it is available.
 */
import type { AlignedRead, Breakpoint, AllTranscripts, BoundaryHint, BoundarySpanning, ExonUsageResponse, GeneModel, KnownVariant, LibraryEvidence, ReadsResponse, RegionHint, SampleCoverage, StructuralEvidence, TranscriptData } from '../components/sashimi/types';
import type { CoverageOptions, ReadsOptions, SampleRef, VariantScan, VariantScanOptions } from '../components/sashimi/datasource';
import { LocalDataSource, isLongRead, type LocalSample, type ReferenceChoice } from './localSource';
import { callSites, collapseReads } from './collapse';
import { arcReadFromAligned, supportsArc } from './arcSupport';
import { phaseReads } from './phasing';
import { haplotagCounts, windowHaplotypes } from './haplotypes';
import { encodeCoverage as encodeCoverageColumns, decodeCoverage as decodeCoverageColumns, encodeReads as encodeReadsColumns, decodeReads as decodeReadsColumns, toBase64, fromBase64, type ReadsPayload } from './columnar';
import type { SessionFile } from './session';
import type { GenomeBuild } from './ensembl';

export const EMBEDDED_APP = 'sashimi-viewer-export';
export const EMBEDDED_VERSION = 2;
export const EMBEDDED_TAG_ID = 'sashimi-embedded';

/** Coverage of one sample over one window, compact: run lengths and depths instead of one object per run. */
export interface EncodedCoverage {
  /** 0-based start of the first run */
  start: number;
  len: number[];
  depth: number[];
  /** [start, end, count] per junction, 0-based half-open */
  junctions: [number, number, number][];
  spanning?: BoundarySpanning;
  window: { start: number; end: number };
  sampled?: { rate: number; total: number; decoded: number };
  spliced?: { reads: number; fraction: number };
  structural?: StructuralEvidence;
  error?: string;
}
/**
 * Coverage of one sample in the columnar form (version 2 exports): the runs and junctions as one compressed stream in
 * base64 (see docs/embedded-format.md); the small facts stay JSON.
 */
export interface EncodedCoverageV2 {
  bin: string;
  window: { start: number; end: number };
  spanning?: BoundarySpanning;
  sampled?: { rate: number; total: number; decoded: number };
  spliced?: { reads: number; fraction: number };
  structural?: StructuralEvidence;
  error?: string;
}
/** Reads of one sample in the columnar form: the stream carries the window, the reads with their pairs and the reference bases. */
export interface EncodedReadsV2 { bin: string; window: { start: number; end: number } }
const isBin = (x: unknown): x is { bin: string } => !!x && typeof (x as any).bin === 'string';

export async function encodeCoverageV2(c: SampleCoverage, window: { start: number; end: number }): Promise<EncodedCoverageV2> {
  return { bin: toBase64(await encodeCoverageColumns({ runs: c.coverage, junctions: c.junctions })), window: c.window ?? window, spanning: c.spanning, sampled: c.sampled, spliced: c.spliced, structural: c.structural, error: c.error };
}
export async function encodeReadsV2(p: ReadsPayload): Promise<EncodedReadsV2> {
  return { bin: toBase64(await encodeReadsColumns(p)), window: p.window };
}

/** Reads of one sample over one window (names replaced by numbers; mismatches and reference bases kept). */
export interface EncodedReads {
  window: { start: number; end: number };
  /** reads passing the filters in the window, before the cap */
  total: number;
  reads: AlignedRead[];
  reference: { start: number; seq: string } | null;
  reference_source: ReadsResponse['reference_source'];
}
export interface EmbeddedView {
  label: string;
  gene: { name: string; id?: string; chrom: string; start: number; end: number };
  transcript: TranscriptData;
  allTranscripts?: AllTranscripts;
  /** neighbouring genes of the window (the queried gene left out) */
  regionGenes?: GeneModel[];
  /** the window fetched (0-based half-open), the view plus margins */
  window: { chrom: string; start: number; end: number };
  uniqueOnly: boolean;
  /** by sample id; the JSON form of version 1 exports or the columnar form of version 2 */
  coverage: Record<string, EncodedCoverage | EncodedCoverageV2>;
  /** reads by sample id, for views exported with their reads track on */
  reads?: Record<string, EncodedReads | EncodedReadsV2>;
}
export interface EmbeddedExport {
  app: typeof EMBEDDED_APP;
  version: number;
  saved: string;
  build: GenomeBuild;
  samples: { id: number; name: string; kind: 'bam' | 'cram'; file: string; index: string; library?: LibraryEvidence }[];
  /** the session (views, options, groups) to restore; its sample names match `samples` */
  session: SessionFile;
  views: EmbeddedView[];
  knownVariants: KnownVariant[];
}

export function encodeCoverage(c: SampleCoverage, window: { start: number; end: number }): EncodedCoverage {
  const runs = c.coverage;
  return {
    start: runs.length ? runs[0].start : window.start,
    len: runs.map(r => r.end - r.start), depth: runs.map(r => r.depth),
    junctions: c.junctions.map(j => [j.start, j.end, j.count]),
    spanning: c.spanning, window: c.window ?? window, sampled: c.sampled, spliced: c.spliced, structural: c.structural, error: c.error,
  };
}
export function decodeCoverage(e: EncodedCoverage, sampleId: number, sampleName: string): SampleCoverage {
  const coverage: SampleCoverage['coverage'] = [];
  let pos = e.start;
  for (let i = 0; i < e.len.length; i++) { coverage.push({ start: pos, end: pos + e.len[i], depth: e.depth[i] }); pos += e.len[i]; }
  return {
    sample_id: sampleId, sample_name: sampleName, coverage,
    junctions: e.junctions.map(([start, end, count]) => ({ start, end, count })),
    spanning: e.spanning, window: e.window, sampled: e.sampled, spliced: e.spliced, structural: e.structural, error: e.error,
  };
}

/** The payload embedded in this page, if it is an exported viewer. */
export function readEmbedded(): EmbeddedExport | null {
  const tag = document.getElementById(EMBEDDED_TAG_ID);
  if (!tag?.textContent) return null;
  try {
    const p = JSON.parse(tag.textContent);
    if (!p || p.app !== EMBEDDED_APP || !Array.isArray(p.views) || !Array.isArray(p.samples)) return null;
    if (typeof p.version === 'number' && p.version > EMBEDDED_VERSION) console.warn(`[sashimi] exported page version ${p.version} is newer than this viewer (${EMBEDDED_VERSION}): reading what it can`);
    return p as EmbeddedExport;
  } catch { return null; }
}

/** True when this page is the unbuilt development entry, whose scripts live in other files: it cannot be exported as one file. */
export const pageIsUnbuilt = (): boolean => !!document.querySelector('script[type="module"][src]');

/**
 * The HTML of this page with the payload embedded. Served over http(s) the pristine source is fetched;
 * opened from disk it is rebuilt from the live document (the rendered application removed, the inline
 * scripts and styles kept, which is everything the single-file build needs).
 */
export async function buildExportHtml(payload: EmbeddedExport): Promise<string> {
  let html: string | null = null;
  if (/^https?:$/.test(location.protocol)) {
    try {
      const r = await fetch(location.href, { cache: 'no-store' });
      if (r.ok) html = await r.text();
    } catch { /* rebuilt from the DOM below */ }
  }
  if (html == null || /<script[^>]*\ssrc=/.test(html)) {
    const doc = document.documentElement.cloneNode(true) as HTMLElement;
    doc.querySelector('#root')?.replaceChildren();
    doc.querySelector(`#${EMBEDDED_TAG_ID}`)?.remove();
    const fb = doc.querySelector('#sashimi-fallback') as HTMLElement | null;
    if (fb) { fb.style.display = 'none'; while (fb.children.length > 1) fb.lastElementChild!.remove(); }
    html = '<!doctype html>\n' + doc.outerHTML;
  } else {
    html = html.replace(new RegExp(`<script id="${EMBEDDED_TAG_ID}"[^>]*>[\\s\\S]*?</script>\\s*`), '');
  }
  const json = JSON.stringify(payload).replace(/</g, '\\u003c');
  const tag = `<script id="${EMBEDDED_TAG_ID}" type="application/json">${json}</script>`;
  const i = html.lastIndexOf('</body>');
  return i >= 0 ? `${html.slice(0, i)}${tag}\n${html.slice(i)}` : html + tag;
}

/** Placeholder sample records of an export: the page shows them as chips, the data source serves them from the payload. */
export function embeddedSamples(p: EmbeddedExport): LocalSample[] {
  return p.samples.map(s => ({ id: s.id, name: s.name, kind: s.kind, file: new File([], s.file), index: new File([], s.index), embedded: true, lib: s.library }));
}

/**
 * Data source of an exported page: embedded regions from memory, everything else through the local source
 * (files the reader adds, network annotations).
 */
export class EmbeddedDataSource extends LocalDataSource {
  private payload: EmbeddedExport;
  private names = new Map<number, string>();
  constructor(payload: EmbeddedExport, reference: ReferenceChoice) {
    super(reference);
    this.payload = payload;
    for (const s of payload.samples) this.names.set(s.id, s.name);
    this.knownVariants = payload.knownVariants ?? [];
  }
  isEmbedded(id: number) { return this.names.has(id); }

  // ---- lazy decoding: a compressed block is decoded the first time its view and sample are asked for, then kept ----
  private covCache = new Map<string, Promise<{ coverage: SampleCoverage['coverage']; junctions: SampleCoverage['junctions'] }>>();
  private readsCache = new Map<string, Promise<ReadsPayload>>();
  private decodedCoverage(vi: number, sid: number): Promise<{ coverage: SampleCoverage['coverage']; junctions: SampleCoverage['junctions'] }> {
    const key = `${vi}|${sid}`;
    let p = this.covCache.get(key);
    if (!p) {
      const c = this.payload.views[vi].coverage[String(sid)];
      p = isBin(c) ? decodeCoverageColumns(fromBase64(c.bin)).then(d => ({ coverage: d.runs, junctions: d.junctions }))
        : Promise.resolve((() => { const d = decodeCoverage(c as EncodedCoverage, sid, ''); return { coverage: d.coverage, junctions: d.junctions }; })());
      p.catch(() => this.covCache.delete(key));
      this.covCache.set(key, p);
    }
    return p;
  }
  private decodedReads(vi: number, sid: number): Promise<ReadsPayload> {
    const key = `${vi}|${sid}`;
    let p = this.readsCache.get(key);
    if (!p) {
      const r = this.payload.views[vi].reads![String(sid)];
      p = isBin(r) ? decodeReadsColumns(fromBase64(r.bin)) : Promise.resolve(r as EncodedReads);
      p.catch(() => this.readsCache.delete(key));
      this.readsCache.set(key, p);
    }
    return p;
  }
  /** Decodes every block of one view (its samples' coverage and reads), so that showing it costs nothing. */
  prefetch(vi: number): Promise<void> {
    const v = this.payload.views[vi];
    if (!v) return Promise.resolve();
    const jobs: Promise<unknown>[] = [];
    for (const sid of Object.keys(v.coverage)) jobs.push(this.decodedCoverage(vi, Number(sid)).catch(() => undefined));
    for (const sid of Object.keys(v.reads ?? {})) jobs.push(this.decodedReads(vi, Number(sid)).catch(() => undefined));
    return Promise.all(jobs).then(() => undefined);
  }
  /** The active view first, then the others one by one in idle moments, so the page stays fluid. */
  prefetchAll(active: number) {
    const rest = this.payload.views.map((_, i) => i).filter(i => i !== active);
    const idle = (f: () => void) => (typeof (globalThis as any).requestIdleCallback === 'function' ? (globalThis as any).requestIdleCallback(f, { timeout: 2000 }) : setTimeout(f, 200));
    const next = () => { const i = rest.shift(); if (i == null) return; void this.prefetch(i).then(() => idle(next)); };
    void this.prefetch(active).then(() => idle(next));
  }
  /** An exported sample has no file to scan: its variant sites come from the exported reads of the window, when they were exported. */
  override async getVariantSites(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan> {
    if (!this.names.has(sampleId)) return super.getVariantSites(sampleId, chrom, start, end, uniqueOnly, minVaf, opts);
    const r = await this.getReads(sampleId, chrom, start, end, uniqueOnly, Number.MAX_SAFE_INTEGER, 'reads', 1, minVaf, opts);
    opts?.onProgress?.(1);
    return { sites: r.sites, total: r.total, long_reads: !!r.long_reads };
  }
  override async rescueClips(sampleId: number, chrom: string, breakpoints: Breakpoint[], uniqueOnly: boolean) {
    if (this.names.has(sampleId)) return [];   // an embedded sample has no reads to go back to
    return super.rescueClips(sampleId, chrom, breakpoints, uniqueOnly);
  }
  override async getPrimaryRecord(sampleId: number, chrom: string, start: number, name: string) {
    if (this.names.has(sampleId)) return null;   // embedded reads carry no names and no file to go back to
    return super.getPrimaryRecord(sampleId, chrom, start, name);
  }
  override async getLibraryType(id: number): Promise<LibraryEvidence> {
    if (!this.names.has(id)) return super.getLibraryType(id);
    return this.payload.samples.find(s => s.id === id)?.library ?? { type: 'unknown', source: 'none', note: 'not recorded in the exported file' };
  }
  override renameSample(id: number, name: string) { if (this.names.has(id)) this.names.set(id, name); else super.renameSample(id, name); }
  override removeSample(id: number) { this.names.delete(id); super.removeSample(id); }
  override list(): SampleRef[] {
    return [...[...this.names].map(([id, name]) => ({ id, name })), ...super.list()];
  }
  private viewsOfGene(geneName: string, geneId?: string) {
    const n = geneName.toUpperCase();
    return this.payload.views.filter(v => v.gene.name.toUpperCase() === n || (geneId && v.gene.id === geneId) || v.transcript.gene_name.toUpperCase() === n);
  }
  override async getTranscript(geneName: string, geneId?: string, hint?: RegionHint): Promise<TranscriptData> {
    const v = this.viewsOfGene(geneName, geneId)[0];
    if (v) return structuredClone(v.transcript);
    return super.getTranscript(geneName, geneId, hint);
  }
  override async getAllTranscripts(geneName: string, geneId?: string, hint?: RegionHint): Promise<AllTranscripts> {
    const v = this.viewsOfGene(geneName, geneId).find(x => x.allTranscripts);
    if (v?.allTranscripts) return structuredClone(v.allTranscripts);
    return super.getAllTranscripts(geneName, geneId, hint);
  }
  override async getRegionGenes(chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]> {
    // a window that holds the whole request answers from memory (1-based inclusive request, 0-based window)
    const v = this.payload.views.find(x => x.regionGenes && sameChrom(x.window.chrom, chrom) && x.window.start <= start - 1 && x.window.end >= end);
    if (v?.regionGenes) return structuredClone(v.regionGenes.filter(g => !exclude || (g.gene_name !== exclude && g.gene_id !== exclude)));
    try { return await super.getRegionGenes(chrom, start, end, exclude); } catch { return []; }
  }
  override async getCoverage(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, boundaries?: BoundaryHint, opts?: CoverageOptions): Promise<SampleCoverage> {
    if (!this.names.has(sampleId)) return super.getCoverage(sampleId, chrom, start, end, uniqueOnly, boundaries, opts);
    const name = this.names.get(sampleId)!;
    const core = opts?.core ?? { start, end };
    // the exported window overlapping the view the most (several views may hold the same sample)
    let bestOv = 0, bestVi = -1;
    this.payload.views.forEach((v, vi) => {
      const c = v.coverage[String(sampleId)];
      if (!c || !sameChrom(v.window.chrom, chrom)) return;
      const ov = Math.min(c.window.end, core.end) - Math.max(c.window.start, core.start);
      if (ov > bestOv) { bestOv = ov; bestVi = vi; }
    });
    if (bestVi < 0) return { sample_id: sampleId, sample_name: name, coverage: [], junctions: [], window: { start, end }, error: `not in this exported file (${chrom}:${(start + 1).toLocaleString('en-US')}-${end.toLocaleString('en-US')}); add the alignment files to see it` };
    const b = this.payload.views[bestVi].coverage[String(sampleId)];
    const { coverage, junctions } = await this.decodedCoverage(bestVi, sampleId);
    return { sample_id: sampleId, sample_name: name, coverage, junctions, spanning: b.spanning, window: b.window, sampled: b.sampled, spliced: b.spliced, structural: b.structural, error: b.error };
  }
  override async getReads(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, maxReads: number,
    mode: 'reads' | 'collapsed', minSupport: number, minVaf: number, opts?: ReadsOptions): Promise<ReadsResponse> {
    if (!this.names.has(sampleId)) return super.getReads(sampleId, chrom, start, end, uniqueOnly, maxReads, mode, minSupport, minVaf, opts);
    const name = this.names.get(sampleId)!;
    // the exported reads window overlapping the request the most
    let bestOv = 0, bestVi = -1;
    this.payload.views.forEach((v, vi) => {
      const r = v.reads?.[String(sampleId)];
      if (!r || !sameChrom(v.window.chrom, chrom)) return;
      const ov = Math.min(r.window.end, end) - Math.max(r.window.start, start);
      if (ov > bestOv) { bestOv = ov; bestVi = vi; }
    });
    if (bestVi < 0) throw new Error(`The reads of this window are not part of this exported file (${chrom}:${(start + 1).toLocaleString('en-US')}-${end.toLocaleString('en-US')}); add the alignment files to see them`);
    const best = await this.decodedReads(bestVi, sampleId);
    const collapsed = mode === 'collapsed';
    let reads = best.reads.filter(r => r.e > start && r.s < end);
    const support = !collapsed ? opts?.support : undefined;
    let total = reads.length, mates = 0;
    if (support) {
      // the supporting reads (every k-th past maxReads), then their mates in the window
      let hits = reads.filter(r => supportsArc(arcReadFromAligned(r, chrom), chrom, support));
      total = hits.length;
      const cap = Math.max(1, maxReads);
      if (hits.length > cap) { const step = hits.length / cap; hits = Array.from({ length: cap }, (_, i) => hits[Math.floor(i * step)]); }
      const kept = new Set(hits), names = new Set(hits.filter(r => r.mp != null).map(r => r.n));
      const extra = reads.filter(r => !kept.has(r) && names.has(r.n));
      mates = extra.length;
      reads = [...hits, ...extra].sort((a, b) => a.s - b.s);
    } else {
      const cap = collapsed ? 40000 : Math.max(100, maxReads);
      if (reads.length > cap) { const step = reads.length / cap; reads = Array.from({ length: cap }, (_, i) => reads[Math.floor(i * step)]); }
    }
    const ref = best.reference?.seq ?? null, refStart = best.reference?.start ?? 0;
    const longReads = isLongRead(reads);
    const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
    const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
    const base = { sample_id: sampleId, sample_name: name, total, shown: reads.length - mates, ...(support ? { supporting: { mates } } : {}), long_reads: longReads, reference: best.reference, reference_source: best.reference_source, haplotags: haplotagCounts(reads) };
    if (collapsed) {
      if (opts?.haplotypes !== 'any') {
        // the window's sites, called once: the phasing, the haplotypes' checks and the answer share them
        const sites = callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel);
        const phaseOf = () => phaseReads(reads, start, end, ref, refStart, 3, vaf, 20, minIndel, sites);
        const { phase, haplotypes } = windowHaplotypes(reads, start, end, ref, refStart, vaf, minIndel, opts?.phaseSource ?? 'auto', phaseOf, sites, longReads);
        return { ...base, reads: [], sites, groups: [], phase, haplotypes };
      }
      const summary = collapseReads(reads, start, end, ref, refStart, 3, vaf, 20, Math.max(1, minSupport), minIndel, longReads);
      return { ...base, reads: [], sites: summary.sites, groups: summary.groups };
    }
    return { ...base, reads, sites: callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel), groups: [] };
  }
  override async getExonUsage(runId: number, chrom: string, strand: number, exons: [number, number][], uniqueOnly: boolean): Promise<ExonUsageResponse> {
    const local = super.list().length ? await super.getExonUsage(runId, chrom, strand, exons, uniqueOnly) : { run_id: runId, chrom, exons, samples: [] };
    const embedded = [...this.names].map(([id, name]) => ({ sample_id: id, sample_name: name, strandness: 'unknown' as const, strand_fraction: null, exons: [], error: 'exon depths are not part of this exported file' }));
    return { ...local, samples: [...embedded, ...local.samples] };
  }
}

const sameChrom = (a: string, b: string) => a.replace(/^chr/i, '').toUpperCase() === b.replace(/^chr/i, '').toUpperCase();
