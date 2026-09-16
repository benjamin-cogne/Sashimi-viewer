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
import type { AllTranscripts, BoundaryHint, BoundarySpanning, ExonUsageResponse, GeneModel, KnownVariant, ReadsResponse, RegionHint, SampleCoverage, TranscriptData } from '../components/sashimi/types';
import type { CoverageOptions, SampleRef } from '../components/sashimi/datasource';
import { LocalDataSource, type LocalSample, type ReferenceChoice } from './localSource';
import type { SessionFile } from './session';
import type { GenomeBuild } from './ensembl';

export const EMBEDDED_APP = 'sashimi-viewer-export';
export const EMBEDDED_VERSION = 1;
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
  error?: string;
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
  /** by sample id */
  coverage: Record<string, EncodedCoverage>;
}
export interface EmbeddedExport {
  app: typeof EMBEDDED_APP;
  version: number;
  saved: string;
  build: GenomeBuild;
  samples: { id: number; name: string; kind: 'bam' | 'cram'; file: string; index: string }[];
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
    spanning: c.spanning, window: c.window ?? window, sampled: c.sampled, error: c.error,
  };
}
export function decodeCoverage(e: EncodedCoverage, sampleId: number, sampleName: string): SampleCoverage {
  const coverage: SampleCoverage['coverage'] = [];
  let pos = e.start;
  for (let i = 0; i < e.len.length; i++) { coverage.push({ start: pos, end: pos + e.len[i], depth: e.depth[i] }); pos += e.len[i]; }
  return {
    sample_id: sampleId, sample_name: sampleName, coverage,
    junctions: e.junctions.map(([start, end, count]) => ({ start, end, count })),
    spanning: e.spanning, window: e.window, sampled: e.sampled, error: e.error,
  };
}

/** The payload embedded in this page, if it is an exported viewer. */
export function readEmbedded(): EmbeddedExport | null {
  const tag = document.getElementById(EMBEDDED_TAG_ID);
  if (!tag?.textContent) return null;
  try {
    const p = JSON.parse(tag.textContent);
    if (!p || p.app !== EMBEDDED_APP || !Array.isArray(p.views) || !Array.isArray(p.samples)) return null;
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
  return p.samples.map(s => ({ id: s.id, name: s.name, kind: s.kind, file: new File([], s.file), index: new File([], s.index), embedded: true }));
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
    let best: EncodedCoverage | null = null, bestOv = 0;
    for (const v of this.payload.views) {
      const c = v.coverage[String(sampleId)];
      if (!c || !sameChrom(v.window.chrom, chrom)) continue;
      const ov = Math.min(c.window.end, core.end) - Math.max(c.window.start, core.start);
      if (ov > bestOv) { best = c; bestOv = ov; }
    }
    if (!best) return { sample_id: sampleId, sample_name: name, coverage: [], junctions: [], window: { start, end }, error: `not in this exported file (${chrom}:${(start + 1).toLocaleString()}-${end.toLocaleString()}); add the alignment files to see it` };
    return decodeCoverage(best, sampleId, name);
  }
  override async getReads(sampleId: number, ...rest: [string, number, number, boolean, number, 'reads' | 'collapsed', number, number]): Promise<ReadsResponse> {
    if (this.names.has(sampleId)) throw new Error('The reads are not part of this exported file; add the alignment files to see them');
    return super.getReads(sampleId, ...rest);
  }
  override async getExonUsage(runId: number, chrom: string, strand: number, exons: [number, number][], uniqueOnly: boolean): Promise<ExonUsageResponse> {
    const local = super.list().length ? await super.getExonUsage(runId, chrom, strand, exons, uniqueOnly) : { run_id: runId, chrom, exons, samples: [] };
    const embedded = [...this.names].map(([id, name]) => ({ sample_id: id, sample_name: name, strandness: 'unknown' as const, strand_fraction: null, exons: [], error: 'exon depths are not part of this exported file' }));
    return { ...local, samples: [...embedded, ...local.samples] };
  }
}

const sameChrom = (a: string, b: string) => a.replace(/^chr/i, '').toUpperCase() === b.replace(/^chr/i, '').toUpperCase();
