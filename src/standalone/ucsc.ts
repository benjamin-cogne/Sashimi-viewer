/**
 * UCSC Genome Browser REST client for the Sashimi viewer (browser-side, like the dbSNP track).
 *
 * Gene models come from the NCBI RefSeq tracks (`ncbiRefSeqCurated`: NM_/NR_ accessions; the
 * `ncbiRefSeqSelect` subset marks RefSeq Select, which on GRCh38 is the MANE Select set for coding
 * genes; the `mane` track, when it answers, flags MANE Select / Plus Clinical exactly). Reference
 * bases come from `getData/sequence`, protein domains from the UniProt / Pfam tracks, gene symbols
 * are resolved with `/search`. Identifiers are therefore RefSeq (NM_000000.0) throughout.
 *
 * The Ensembl REST client (ensembl.ts) stays as the fallback: symbol resolution when `/search`
 * has no hit (aliases, previous symbols), and every call when the UCSC API is unreachable.
 *
 * Coordinates: the UCSC API is 0-based half-open; the viewer types are 1-based inclusive, as the
 * Ensembl payloads were, so the conversion happens here.
 */
import type { TranscriptData, AllTranscripts, TranscriptModel, GeneModel, ProteinDomain, ProteinModelRef, RegionHint } from '../components/sashimi/types';
import * as ensembl from './ensembl';
import type { GenomeBuild } from './ensembl';

export const UCSC_API = 'https://api.genome.ucsc.edu';
const GENOME: Record<GenomeBuild, string> = { GRCh38: 'hg38', GRCh37: 'hg19' };

/** Set when the UCSC API refused or failed a request at the network level; later calls go to Ensembl. */
let ucscBroken: string | null = null;
export function annotationSourceLabel(): string { return ucscBroken ? 'Ensembl (UCSC API unavailable)' : 'RefSeq (UCSC)'; }

const ucscChrom = (chrom: string): string => {
  const bare = chrom.startsWith('chr') ? chrom.slice(3) : chrom;
  return `chr${bare === 'MT' ? 'M' : bare}`;
};
const viewerChrom = (chrom: string): string => (chrom.startsWith('chr') ? chrom : `chr${chrom}`);

class UcscError extends Error { constructor(msg: string, public network: boolean) { super(msg); } }

async function getJson(url: string): Promise<any> {
  let r: Response;
  try { r = await fetch(url, { headers: { Accept: 'application/json' } }); }
  catch (e: any) { throw new UcscError(`UCSC API unreachable from the browser: ${e?.message || e}`, true); }
  if (!r.ok) {
    let detail = '';
    try { const j = await r.json(); detail = j?.error ? `: ${String(j.error).slice(0, 200)}` : ''; } catch { /* no body */ }
    throw new UcscError(`UCSC API ${r.status}${detail}`, r.status >= 500);
  }
  return r.json();
}

/** Items of a `getData/track` answer: `{track: [...]}`, `{track: {chr1: [...]}}` or a bare array. */
function trackItems(data: any, track: string): any[] {
  const v = data?.[track] ?? data;
  if (Array.isArray(v)) return v;
  if (v && typeof v === 'object') {
    const arrays = Object.values(v).filter(Array.isArray) as any[][];
    if (arrays.length) return arrays.flat();
  }
  if (data?.error) throw new UcscError(`UCSC API: ${String(data.error).slice(0, 200)}`, false);
  return [];
}

// ======================== genePred models (RefSeq tracks) ========================

/** One RefSeq transcript as the genePred tracks describe it; 0-based half-open. */
interface GenePred {
  name: string; symbol: string; chrom: string; strand: 1 | -1;
  txStart: number; txEnd: number; cdsStart: number; cdsEnd: number;
  exons: { start: number; end: number }[];
}

const intList = (v: unknown): number[] =>
  Array.isArray(v) ? v.map(Number) : typeof v === 'string' ? v.split(',').filter(s => s !== '').map(Number) : [];

function parseGenePred(it: any): GenePred | null {
  const name = String(it?.name ?? '');
  const starts = intList(it?.exonStarts), ends = intList(it?.exonEnds);
  if (!name || !starts.length || starts.length !== ends.length) return null;
  const exons = starts.map((s, i) => ({ start: s, end: ends[i] })).filter(e => Number.isFinite(e.start) && e.end > e.start).sort((a, b) => a.start - b.start);
  if (!exons.length) return null;
  return {
    name, symbol: String(it.name2 ?? it.geneName ?? ''), chrom: viewerChrom(String(it.chrom ?? '')), strand: String(it.strand) === '-' ? -1 : 1,
    txStart: Number(it.txStart ?? exons[0].start), txEnd: Number(it.txEnd ?? exons[exons.length - 1].end),
    cdsStart: Number(it.cdsStart ?? 0), cdsEnd: Number(it.cdsEnd ?? 0), exons,
  };
}

const accession = (name: string) => name.split('.')[0];
const isCoding = (g: GenePred) => g.cdsEnd > g.cdsStart;
const cdsLength = (g: GenePred) => isCoding(g) ? g.exons.reduce((a, e) => a + Math.max(0, Math.min(e.end, g.cdsEnd) - Math.max(e.start, g.cdsStart)), 0) : 0;
const txLength = (g: GenePred) => g.exons.reduce((a, e) => a + e.end - e.start, 0);

const trackCache = new Map<string, Promise<any[]>>();
/** Raw items of a track over a 0-based half-open window (cached per exact window). */
function trackWindow(build: GenomeBuild, track: string, chrom: string, start: number, end: number): Promise<any[]> {
  const key = `${build}|${track}|${chrom}|${start}-${end}`;
  if (!trackCache.has(key)) {
    const url = `${UCSC_API}/getData/track?genome=${GENOME[build]};track=${track};chrom=${encodeURIComponent(ucscChrom(chrom))};start=${Math.max(0, start)};end=${end};maxItemsOutput=100000`;
    trackCache.set(key, getJson(url).then(d => trackItems(d, track)).catch(e => { trackCache.delete(key); throw e; }));
  }
  return trackCache.get(key)!;
}

async function genePreds(build: GenomeBuild, track: string, chrom: string, start: number, end: number): Promise<GenePred[]> {
  const items = await trackWindow(build, track, chrom, start, end);
  return items.map(parseGenePred).filter((g): g is GenePred => !!g);
}

/** Accessions (without version) of the RefSeq Select models in a window; empty when the track does not answer. */
async function selectAccessions(build: GenomeBuild, chrom: string, start: number, end: number): Promise<Set<string>> {
  try { return new Set((await genePreds(build, 'ncbiRefSeqSelect', chrom, start, end)).map(g => accession(g.name))); }
  catch { return new Set(); }
}

/**
 * MANE status by RefSeq accession from the `mane` track (GRCh38 only): 'select' or 'plus'. The track's
 * field names are not relied upon: any NM_/NR_ accession found in an item is taken, and an item whose
 * text mentions "Plus Clinical" is a MANE Plus Clinical model.
 */
async function maneStatus(build: GenomeBuild, chrom: string, start: number, end: number): Promise<Map<string, 'select' | 'plus'>> {
  const out = new Map<string, 'select' | 'plus'>();
  if (build !== 'GRCh38') return out;
  let items: any[];
  try { items = await trackWindow(build, 'mane', chrom, start, end); } catch { return out; }
  for (const it of items) {
    if (!it || typeof it !== 'object') continue;
    const values = Object.values(it).filter((v): v is string => typeof v === 'string');
    const acc = values.map(v => v.match(/^N[MR]_\d+/)?.[0]).find(Boolean);
    if (!acc) continue;
    const plus = values.some(v => /plus\s*clinical/i.test(v));
    if (!out.has(acc) || !plus) out.set(acc, plus ? 'plus' : 'select');
  }
  return out;
}

// ======================== Locating a gene ========================

interface Locus { chrom: string; start: number; end: number; symbol: string }   // 1-based inclusive

/** Objects of a `/search` answer carrying a `position` ("chr17:43044295-43125483") and a name. */
function searchHits(node: any, out: { name: string; chrom: string; start: number; end: number }[] = []): typeof out {
  if (Array.isArray(node)) { node.forEach(n => searchHits(n, out)); return out; }
  if (!node || typeof node !== 'object') return out;
  const pos = typeof node.position === 'string' ? node.position.match(/^(chr[\w.]+):([\d,]+)-([\d,]+)$/) : null;
  if (pos) {
    const name = String(node.posName ?? node.name ?? node.value ?? node.hgFindMatches ?? '');
    out.push({ name, chrom: pos[1], start: Number(pos[2].replace(/,/g, '')), end: Number(pos[3].replace(/,/g, '')) });
  }
  for (const v of Object.values(node)) if (v && typeof v === 'object') searchHits(v, out);
  return out;
}

const locusCache = new Map<string, Promise<Locus>>();
/** Gene span by symbol (or Ensembl id): UCSC `/search`, then Ensembl (aliases, previous symbols, ENSG ids). */
function locateGene(build: GenomeBuild, geneName: string, geneId?: string): Promise<Locus> {
  const key = `${build}|${geneId?.startsWith('ENSG') ? geneId : geneName}`;
  if (!locusCache.has(key)) {
    const run = async (): Promise<Locus> => {
      const q = geneName || geneId || '';
      if (!ucscBroken && q && !q.toUpperCase().startsWith('ENSG')) {
        try {
          const d = await getJson(`${UCSC_API}/search?genome=${GENOME[build]};search=${encodeURIComponent(q)}`);
          const hits = searchHits(d).filter(h => h.name.toUpperCase() === q.toUpperCase() && h.end > h.start);
          if (hits.length) {
            // one span covering every exact hit on the most frequent chromosome (transcripts of the same gene)
            const byChrom = new Map<string, typeof hits>();
            for (const h of hits) byChrom.set(h.chrom, [...(byChrom.get(h.chrom) || []), h]);
            const best = [...byChrom.values()].sort((a, b) => b.length - a.length)[0];
            return { chrom: viewerChrom(best[0].chrom), start: Math.min(...best.map(h => h.start)), end: Math.max(...best.map(h => h.end)), symbol: q };
          }
        } catch (e: any) { if (e instanceof UcscError && e.network) ucscBroken = e.message; }
      }
      const g = await ensembl.locateGene(build, geneName, geneId);
      return { chrom: g.chrom, start: g.start, end: g.end, symbol: g.symbol || geneName };
    };
    locusCache.set(key, run().catch(e => { locusCache.delete(key); throw e; }));
  }
  return locusCache.get(key)!;
}

// ======================== Transcript models of a gene ========================

interface GeneModels { symbol: string; chrom: string; strand: 1 | -1; start: number; end: number; models: GenePred[]; select: Set<string>; mane: Map<string, 'select' | 'plus'>; /** GRCh38: RefSeq Select is the MANE Select set for coding genes */ build38: boolean }

const geneCache = new Map<string, Promise<GeneModels>>();

/**
 * Every RefSeq model of a gene. With a region hint (the gene span of the outlier tables) no symbol
 * lookup is needed: the curated models of that region whose symbol matches are taken; when none
 * matches (symbol changed between annotation releases) the gene with the largest overlap is used.
 */
function geneModels(build: GenomeBuild, geneName: string, geneId?: string, hint?: RegionHint): Promise<GeneModels> {
  const key = `${build}|${(geneId?.startsWith('ENSG') ? geneId : geneName).toUpperCase()}|${hint ? `${hint.chrom}:${hint.start}-${hint.end}` : ''}`;
  if (!geneCache.has(key)) {
    const run = async (): Promise<GeneModels> => {
      let locus: Locus;
      if (hint && hint.chrom && hint.end > hint.start) locus = { chrom: viewerChrom(hint.chrom), start: hint.start, end: hint.end, symbol: geneName };
      else locus = await locateGene(build, geneName, geneId);
      const wanted = (locus.symbol || geneName).toUpperCase();
      // margin: a model may start a little before the span the tables carry
      const wStart = Math.max(0, locus.start - 1 - 5000), wEnd = locus.end + 5000;
      let all = await genePreds(build, 'ncbiRefSeqCurated', locus.chrom, wStart, wEnd);
      let mine = all.filter(g => g.symbol.toUpperCase() === wanted);
      if (!mine.length && !hint && geneName && locus.symbol.toUpperCase() !== geneName.toUpperCase()) mine = all.filter(g => g.symbol.toUpperCase() === geneName.toUpperCase());
      if (!mine.length) {
        // no curated model: predicted models (XM_/XR_) of the full RefSeq track
        try { all = await genePreds(build, 'ncbiRefSeq', locus.chrom, wStart, wEnd); mine = all.filter(g => g.symbol.toUpperCase() === wanted); } catch { /* keep curated */ }
      }
      if (!mine.length && hint) {
        // symbol mismatch: the gene overlapping the hinted span the most
        const ov = (g: GenePred) => Math.max(0, Math.min(g.txEnd, locus.end) - Math.max(g.txStart, locus.start - 1));
        const byGene = new Map<string, number>();
        for (const g of all) if (g.symbol) byGene.set(g.symbol, (byGene.get(g.symbol) || 0) + ov(g));
        const best = [...byGene.entries()].filter(([, o]) => o > 0).sort((a, b) => b[1] - a[1])[0];
        if (best) mine = all.filter(g => g.symbol === best[0]);
      }
      if (!mine.length) throw new Error(`no RefSeq model for ${geneName || geneId} in ${locus.chrom}:${locus.start.toLocaleString('en-US')}-${locus.end.toLocaleString('en-US')} (UCSC ncbiRefSeq)`);
      const strand = mine.filter(g => g.strand === -1).length > mine.length / 2 ? -1 : 1;
      const start = Math.min(...mine.map(g => g.txStart)), end = Math.max(...mine.map(g => g.txEnd));
      const [select, mane] = await Promise.all([selectAccessions(build, locus.chrom, start, end), maneStatus(build, locus.chrom, start, end)]);
      return { symbol: mine[0].symbol || geneName, chrom: locus.chrom, strand, start, end, models: mine, select, mane, build38: build === 'GRCh38' };
    };
    geneCache.set(key, run().catch(e => { geneCache.delete(key); throw e; }));
  }
  return geneCache.get(key)!;
}

const isMane = (g: GeneModels, m: GenePred) => g.mane.get(accession(m.name)) === 'select';
const isSelect = (g: GeneModels, m: GenePred) => g.select.has(accession(m.name));
const biotypeOf = (m: GenePred) => isCoding(m) ? 'protein_coding' : m.name.startsWith('NR_') || m.name.startsWith('XR_') ? 'ncRNA' : 'non_coding';

/** MANE Select, else RefSeq Select, else the coding model with the longest CDS, else the longest model. */
function pickModel(g: GeneModels): { model: GenePred; kind: 'mane' | 'canonical' | 'longest' } {
  const mane = g.models.find(m => isMane(g, m));
  if (mane) return { model: mane, kind: 'mane' };
  const sel = g.models.filter(m => isSelect(g, m)).sort((a, b) => cdsLength(b) - cdsLength(a))[0];
  if (sel) return { model: sel, kind: g.build38 ? 'mane' : 'canonical' };
  const coding = g.models.filter(isCoding);
  if (coding.length) return { model: coding.reduce((b, m) => (cdsLength(m) > cdsLength(b) ? m : b)), kind: 'longest' };
  return { model: g.models.reduce((b, m) => (txLength(m) > txLength(b) ? m : b)), kind: 'longest' };
}

function toTranscriptData(g: GeneModels): TranscriptData {
  const { model: m, kind } = pickModel(g);
  const exons = m.exons.map((e, i) => ({ start: e.start + 1, end: e.end, rank: m.strand > 0 ? i + 1 : m.exons.length - i }));
  return {
    gene_name: g.symbol, transcript_id: m.name, translation_id: null, is_mane_select: kind === 'mane', model_kind: kind, biotype: biotypeOf(m),
    // the model's own span, not the gene's: another isoform starting upstream must not draw as an intron before exon 1
    source: 'refseq', chrom: g.chrom, strand: m.strand, start: m.txStart + 1, end: m.txEnd, exons,
    cds_start: isCoding(m) ? m.cdsStart + 1 : null, cds_end: isCoding(m) ? m.cdsEnd : null,
  };
}

function toModels(g: GeneModels): TranscriptModel[] {
  const rank = (m: TranscriptModel) => (m.is_mane ? 0 : m.is_canonical ? 1 : 2) * 10 + (m.biotype === 'protein_coding' ? 0 : 1);
  return g.models.map((m): TranscriptModel => ({
    id: m.name, name: m.name, source: 'refseq', biotype: biotypeOf(m), start: m.txStart + 1, end: m.txEnd, strand: m.strand,
    exons: m.exons.map(e => ({ start: e.start + 1, end: e.end })),
    cds_start: isCoding(m) ? m.cdsStart + 1 : null, cds_end: isCoding(m) ? m.cdsEnd : null,
    is_mane: isMane(g, m) || (g.build38 && isSelect(g, m)), is_canonical: isSelect(g, m),
  })).sort((a, b) => rank(a) - rank(b) || a.id.localeCompare(b.id, undefined, { numeric: true }));
}

/** Falls back on Ensembl when the UCSC API cannot be reached (identifiers are then ENST). */
async function withFallback<T>(ucsc: () => Promise<T>, fallback: () => Promise<T>): Promise<T> {
  if (ucscBroken) return fallback();
  try { return await ucsc(); }
  catch (e: any) {
    if (e instanceof UcscError && e.network) { ucscBroken = e.message; console.warn(`[sashimi] UCSC API unavailable (${e.message}); using Ensembl REST`); return fallback(); }
    throw e;
  }
}

export function getTranscript(build: GenomeBuild, geneName: string, geneId?: string, hint?: RegionHint): Promise<TranscriptData> {
  return withFallback(async () => toTranscriptData(await geneModels(build, geneName, geneId, hint)), () => ensembl.getTranscript(build, geneName, geneId));
}

export function getAllTranscripts(build: GenomeBuild, geneName: string, geneId?: string, hint?: RegionHint): Promise<AllTranscripts> {
  return withFallback(async () => {
    const g = await geneModels(build, geneName, geneId, hint);
    return { gene_name: g.symbol, chrom: g.chrom, strand: g.strand, source: 'refseq', transcripts: toModels(g) };
  }, () => ensembl.getAllTranscripts(build, geneName, geneId));
}

// ======================== Genes around a region (neighbours), 500 kb chunks cached ========================

const REGION_CHUNK = 500_000;
const regionCache = new Map<string, Promise<GeneModel[]>>();

function regionChunk(build: GenomeBuild, chrom: string, chunk: number): Promise<GeneModel[]> {
  const key = `${build}|${chrom}|${chunk}`;
  if (!regionCache.has(key)) {
    regionCache.set(key, (async () => {
      const start = chunk * REGION_CHUNK, end = (chunk + 1) * REGION_CHUNK;
      const [models, select] = await Promise.all([genePreds(build, 'ncbiRefSeqCurated', chrom, start, end), selectAccessions(build, chrom, start, end)]);
      const byGene = new Map<string, GenePred[]>();
      for (const m of models) if (m.symbol) byGene.set(m.symbol, [...(byGene.get(m.symbol) || []), m]);
      const out: GeneModel[] = [];
      for (const [symbol, ms] of byGene) {
        const sel = ms.filter(m => select.has(accession(m.name))).sort((a, b) => cdsLength(b) - cdsLength(a))[0];
        const coding = ms.filter(isCoding);
        const m = sel ?? (coding.length ? coding.reduce((b, x) => (cdsLength(x) > cdsLength(b) ? x : b)) : ms.reduce((b, x) => (txLength(x) > txLength(b) ? x : b)));
        out.push({
          gene_id: symbol, gene_name: symbol, biotype: biotypeOf(m), strand: m.strand,
          start: Math.min(...ms.map(x => x.txStart)) + 1, end: Math.max(...ms.map(x => x.txEnd)),
          transcript_id: m.name, is_canonical: !!sel, exons: m.exons.map(e => ({ start: e.start + 1, end: e.end })),
          cds_start: isCoding(m) ? m.cdsStart + 1 : null, cds_end: isCoding(m) ? m.cdsEnd : null,
        });
      }
      return out;
    })().catch(err => { regionCache.delete(key); throw err; }));
  }
  return regionCache.get(key)!;
}

/** Genes overlapping [start, end] (1-based inclusive) with their RefSeq Select model, `exclude` left out. */
export function getRegionGenes(build: GenomeBuild, chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]> {
  return withFallback(async () => {
    if (end < start) return [];
    const seen = new Map<string, GeneModel>();
    const chunks: Promise<GeneModel[]>[] = [];
    for (let chunk = Math.floor((start - 1) / REGION_CHUNK); chunk <= Math.floor((end - 1) / REGION_CHUNK); chunk++) chunks.push(regionChunk(build, chrom, chunk));
    for (const list of await Promise.all(chunks)) {
      for (const g of list) {
        if (g.end < start || g.start > end) continue;
        if (exclude && g.gene_name.toUpperCase() === exclude.toUpperCase()) continue;
        if (!seen.has(g.gene_id)) seen.set(g.gene_id, g);
      }
    }
    return [...seen.values()].sort((a, b) => a.start - b.start);
  }, () => ensembl.getRegionGenes(build, chrom, start, end, exclude));
}

// ======================== Reference sequence (128 kb chunks, fetched in parallel, cached) ========================

const CHUNK = 1 << 17;
const MAX_PARALLEL = 6;
const chunkCache = new Map<string, Promise<string | null>>();
let running = 0;
const waiting: (() => void)[] = [];
const acquire = () => new Promise<void>(res => { if (running < MAX_PARALLEL) { running++; res(); } else waiting.push(() => { running++; res(); }); });
const release = () => { running--; const next = waiting.shift(); if (next) next(); };

function fetchChunk(build: GenomeBuild, chrom: string, idx: number): Promise<string | null> {
  const key = `${build}|${chrom}|${idx}`;
  if (!chunkCache.has(key)) {
    chunkCache.set(key, (async () => {
      await acquire();
      try {
        const url = `${UCSC_API}/getData/sequence?genome=${GENOME[build]};chrom=${encodeURIComponent(ucscChrom(chrom))};start=${idx * CHUNK};end=${(idx + 1) * CHUNK}`;
        const d = await getJson(url);
        const dna = typeof d?.dna === 'string' ? d.dna : null;
        if (dna == null) throw new UcscError(d?.error ? `UCSC API: ${String(d.error).slice(0, 200)}` : 'UCSC API: no sequence in the answer', false);
        return dna.toUpperCase();
      } finally { release(); }
    })().catch(e => { chunkCache.delete(key); throw e; }));
  }
  return chunkCache.get(key)!;
}

/** Reference bases for [start, end) (0-based half-open, upper case); Ensembl when UCSC fails; null when neither answers. */
export async function getReference(build: GenomeBuild, chrom: string, start: number, end: number): Promise<string | null> {
  if (end <= start) return '';
  start = Math.max(0, start);
  if (!ucscBroken) {
    try {
      const idx: number[] = [];
      for (let i = Math.floor(start / CHUNK); i <= Math.floor((end - 1) / CHUNK); i++) idx.push(i);
      const parts = await Promise.all(idx.map(i => fetchChunk(build, chrom, i)));
      if (parts.every((p): p is string => p != null)) {
        const offset = start - idx[0] * CHUNK;
        const seq = parts.join('').substring(offset, offset + (end - start));
        if (seq.length === end - start || parts[parts.length - 1]!.length < CHUNK) return seq;   // short last chunk = chromosome end
      }
    } catch (e: any) {
      if (e instanceof UcscError && e.network) ucscBroken = e.message;
      console.warn(`[sashimi] UCSC sequence unavailable (${e?.message || e}); trying Ensembl`);
    }
  }
  try { return await ensembl.getReference(build, chrom, start, end); } catch { return null; }
}

// ======================== Protein domains (UniProt / Pfam tracks mapped onto the CDS) ========================

/**
 * Genomic → amino-acid mapping of one transcript: coding offset (0-based, in transcription order) of a
 * genomic position, positions inside introns snapping to the next coding base.
 */
function codingOffset(model: ProteinModelRef): { at: (pos: number) => number; length: number } {
  const cds = model.exons.map(e => ({ start: Math.max(e.start, model.cdsStart), end: Math.min(e.end, model.cdsEnd) })).filter(e => e.end > e.start).sort((a, b) => a.start - b.start);
  const length = cds.reduce((a, e) => a + e.end - e.start, 0);
  const forward = (pos: number) => {
    let off = 0;
    for (const e of cds) {
      if (pos >= e.end) { off += e.end - e.start; continue; }
      return off + Math.max(0, pos - e.start);
    }
    return length;
  };
  return { at: model.strand < 0 ? (pos: number) => Math.max(0, length - 1 - forward(pos)) : forward, length };
}

/** Domain items of a bed-like track over the CDS, as amino-acid intervals on the model. */
async function domainTrack(build: GenomeBuild, track: string, type: string, model: ProteinModelRef): Promise<ProteinDomain[]> {
  const items = await trackWindow(build, track, model.chrom, model.cdsStart, model.cdsEnd);
  const { at, length } = codingOffset(model);
  const aaLen = Math.max(1, Math.floor(length / 3));
  const out: ProteinDomain[] = [];
  for (const it of items) {
    const s = Number(it?.chromStart), e = Number(it?.chromEnd);
    if (!Number.isFinite(s) || !Number.isFinite(e) || e <= s) continue;
    if (e <= model.cdsStart || s >= model.cdsEnd) continue;
    const strandOk = !it.strand || String(it.strand) === (model.strand < 0 ? '-' : '+');
    if (!strandOk) continue;
    const lo = Math.max(s, model.cdsStart), hi = Math.min(e, model.cdsEnd) - 1;
    const a = at(lo), b = at(hi);
    const aaStart = Math.min(aaLen, Math.floor(Math.min(a, b) / 3) + 1), aaEnd = Math.min(aaLen, Math.floor(Math.max(a, b) / 3) + 1);
    if (aaEnd < aaStart) continue;
    const name = String(it.name ?? it.description ?? '').trim();
    if (!name) continue;
    out.push({ type, id: name, description: String(it.description ?? it.longName ?? it.uniprotName ?? name), interpro: null, start: aaStart, end: aaEnd });
  }
  return out;
}

/** Protein domains of a coding model: UniProt domains first, Pfam when UniProt has none; empty when the tracks fail. */
export async function getProteinDomains(build: GenomeBuild, model: ProteinModelRef): Promise<ProteinDomain[]> {
  if (model.cdsEnd <= model.cdsStart) return [];
  if (!ucscBroken) {
    try {
      const unip = await domainTrack(build, 'unipDomain', 'UniProt', model);
      if (unip.length) return unip;
      return await domainTrack(build, 'ucscGenePfam', 'Pfam', model);
    } catch (e: any) {
      if (e instanceof UcscError && e.network) ucscBroken = e.message;
      console.warn(`[sashimi] UCSC domain tracks unavailable (${e?.message || e})`);
    }
  }
  if (model.translationId) return ensembl.getProteinDomains(build, model.translationId);
  return [];
}

/** Pure helpers exposed for tests (not part of the viewer API). */
export const _internal = { parseGenePred, codingOffset, searchHits, trackItems };
