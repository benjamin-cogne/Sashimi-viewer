/**
 * Ensembl REST client (transcripts, all transcript models, reference sequence, protein features).
 * Mirrors backend/transcript_lookup.py and backend/reference.py. Since the viewer moved to the
 * UCSC API (ucsc.ts, RefSeq identifiers) this is the fallback: gene symbols the UCSC search does
 * not know (aliases, previous symbols, ENSG ids) and every call when the UCSC API is unreachable.
 * Ensembl REST answers cross-origin requests, so this runs straight from the browser.
 */
import type { TranscriptData, AllTranscripts, TranscriptModel, GeneModel, ProteinDomain } from '../components/sashimi/types';

export type GenomeBuild = 'GRCh38' | 'GRCh37';
const HOSTS: Record<GenomeBuild, string> = { GRCh38: 'https://rest.ensembl.org', GRCh37: 'https://grch37.rest.ensembl.org' };

async function getJson(url: string): Promise<any> {
  const resp = await fetch(url, { headers: { Accept: 'application/json' } });
  if (!resp.ok) throw new Error(`Ensembl ${resp.status} for ${url.replace(/\?.*/, '')}`);
  return resp.json();
}

const cdsLength = (tx: any): number => {
  const tr = tx.Translation || {};
  if (tr.start == null || tr.end == null) return 0;
  const lo = Math.min(tr.start, tr.end), hi = Math.max(tr.start, tr.end);
  return (tx.Exon || []).reduce((a: number, e: any) => a + Math.max(0, Math.min(e.end, hi) - Math.max(e.start, lo) + 1), 0);
};
const txLength = (tx: any): number => (tx.Exon || []).reduce((a: number, e: any) => a + e.end - e.start + 1, 0);

/** MANE Select, else Ensembl canonical, else the protein-coding transcript with the longest CDS, else the longest. */
function pickTranscript(transcripts: any[]): any {
  for (const key of ['is_mane_select', 'is_canonical']) { const t = transcripts.find(x => x[key] && (x.Exon || []).length); if (t) return t; }
  const coding = transcripts.filter(x => x.biotype === 'protein_coding' && (x.Exon || []).length);
  if (coding.length) return coding.reduce((b, x) => (cdsLength(x) > cdsLength(b) ? x : b));
  const withExons = transcripts.filter(x => (x.Exon || []).length);
  return withExons.length ? withExons.reduce((b, x) => (txLength(x) > txLength(b) ? x : b)) : transcripts[0];
}

function buildTranscript(data: any, label: string): TranscriptData {
  const transcripts: any[] = data.Transcript || [];
  if (!transcripts.length) throw new Error(`No transcripts for ${label}`);
  const tx = pickTranscript(transcripts);
  const strand = data.strand ?? 1;
  // Exon numbering follows transcription: on the minus strand exon 1 is the rightmost one.
  const sortedExons = [...(tx.Exon || [])].sort((a, b) => a.start - b.start);
  const exons = sortedExons.map((e, i) => ({ start: e.start, end: e.end, rank: strand >= 0 ? i + 1 : sortedExons.length - i }));
  const tr = tx.Translation || {};
  let cs = tr.start ?? null, ce = tr.end ?? null;
  if (cs != null && ce != null && cs > ce) [cs, ce] = [ce, cs];
  const seq = String(data.seq_region_name || '');
  return {
    gene_name: data.display_name || label, transcript_id: tx.id || '', translation_id: tr.id ?? null, is_mane_select: !!tx.is_mane_select,
    model_kind: tx.is_mane_select ? 'mane' : tx.is_canonical ? 'canonical' : 'longest', biotype: tx.biotype || '', source: 'ensembl',
    // the transcript's own span (data.start/end is the gene's, which may extend beyond this model)
    chrom: seq.startsWith('chr') ? seq : `chr${seq}`, strand, start: tx.start ?? exons[0]?.start ?? data.start ?? 0, end: tx.end ?? exons[exons.length - 1]?.end ?? data.end ?? 0,
    exons, cds_start: cs, cds_end: ce,
  };
}

const geneCache = new Map<string, Promise<any>>();
/** /lookup with expand=1: null on 404 (unknown id or symbol), throws on API errors. */
async function lookupRaw(host: string, path: string): Promise<any | null> {
  const resp = await fetch(`${host}${path}?expand=1;content-type=application/json`, { headers: { Accept: 'application/json' } });
  if (resp.status === 404 || resp.status === 400) return null;
  if (!resp.ok) throw new Error(`Ensembl REST error ${resp.status}`);
  return resp.json();
}
/**
 * Gene record by Ensembl id, then symbol, then symbol aliases / previous symbols (xrefs/symbol),
 * so renamed genes and old symbols in the outlier tables still get a model.
 */
async function lookupGene(build: GenomeBuild, geneName: string, geneId?: string): Promise<any> {
  const host = HOSTS[build];
  const key = `${build}|${geneId?.startsWith('ENSG') ? geneId : geneName}`;
  if (!geneCache.has(key)) {
    const run = async () => {
      const tried: string[] = [];
      if (geneId?.startsWith('ENSG')) { const d = await lookupRaw(host, `/lookup/id/${geneId.split('.')[0]}`); tried.push(geneId); if (d) return d; }
      if (geneName.startsWith('ENSG')) { const d = await lookupRaw(host, `/lookup/id/${geneName.split('.')[0]}`); tried.push(geneName); if (d) return d; }
      else if (geneName) {
        const d = await lookupRaw(host, `/lookup/symbol/homo_sapiens/${encodeURIComponent(geneName)}`); tried.push(`symbol ${geneName}`); if (d) return d;
        const xr = await getJson(`${host}/xrefs/symbol/homo_sapiens/${encodeURIComponent(geneName)}?object_type=gene;content-type=application/json`).catch(() => []);
        const ids: string[] = [...new Set((Array.isArray(xr) ? xr : []).map((x: any) => String(x.id || '')).filter(id => id.startsWith('ENSG')))].slice(0, 5);
        const cands = (await Promise.all(ids.map(id => lookupRaw(host, `/lookup/id/${id}`).catch(() => null)))).filter(Boolean);
        tried.push(`aliases of ${geneName} (${cands.length} gene${cands.length === 1 ? '' : 's'})`);
        if (cands.length) {
          cands.sort((a: any, b: any) => Number(a.biotype !== 'protein_coding') - Number(b.biotype !== 'protein_coding') || Number((a.display_name || '').toUpperCase() !== geneName.toUpperCase()) - Number((b.display_name || '').toUpperCase() !== geneName.toUpperCase()));
          return cands[0];
        }
      }
      throw new Error(`not found in Ensembl (tried ${tried.join('; ')})`);
    };
    geneCache.set(key, run().catch(e => { geneCache.delete(key); throw e; }));
  }
  return geneCache.get(key)!;
}

/** Gene span (1-based inclusive) and current symbol, by Ensembl id, symbol or alias. */
export async function locateGene(build: GenomeBuild, geneName: string, geneId?: string): Promise<{ chrom: string; start: number; end: number; strand: number; symbol: string }> {
  let data: any = null;
  if (geneId?.startsWith('ENSG')) { try { data = await lookupGene(build, geneName, geneId); } catch { data = null; } }
  if (!data) data = await lookupGene(build, geneName);
  const seq = String(data.seq_region_name || '');
  return { chrom: seq.startsWith('chr') ? seq : `chr${seq}`, start: data.start ?? 0, end: data.end ?? 0, strand: data.strand ?? 1, symbol: data.display_name || geneName };
}

export async function getTranscript(build: GenomeBuild, geneName: string, geneId?: string): Promise<TranscriptData> {
  let data: any = null;
  if (geneId?.startsWith('ENSG')) { try { data = await lookupGene(build, geneName, geneId); } catch { data = null; } }
  if (!data) data = await lookupGene(build, geneName);
  return buildTranscript(data, geneName);
}

function ensemblModels(data: any): TranscriptModel[] {
  return (data.Transcript || []).map((tx: any) => {
    const exons = [...(tx.Exon || [])].sort((a, b) => a.start - b.start).map(e => ({ start: e.start, end: e.end }));
    const tr = tx.Translation || {};
    return {
      id: tx.id || '', name: tx.display_name || tx.id || '', source: 'ensembl' as const, biotype: tx.biotype || '',
      start: tx.start ?? exons[0]?.start ?? 0, end: tx.end ?? exons[exons.length - 1]?.end ?? 0, strand: tx.strand ?? data.strand ?? 1,
      exons, cds_start: tr.start ?? null, cds_end: tr.end ?? null, is_mane: !!tx.is_mane_select, is_canonical: !!tx.is_canonical,
    };
  }).filter((m: TranscriptModel) => m.exons.length > 0);
}

async function refseqModels(build: GenomeBuild, chrom: string, start: number, end: number, strand: number): Promise<TranscriptModel[]> {
  const host = HOSTS[build];
  const bare = chrom.startsWith('chr') ? chrom.slice(3) : chrom;
  const region = `${bare}:${start}-${end}`;
  const base = 'db_type=otherfeatures;logic_name=refseq_import;content-type=application/json';
  let txs: any[];
  try { txs = await getJson(`${host}/overlap/region/human/${region}?feature=transcript;${base}`); } catch { return []; }
  if (!Array.isArray(txs) || !txs.length) return [];
  const [exons, cdss] = await Promise.all([
    getJson(`${host}/overlap/region/human/${region}?feature=exon;${base}`).catch(() => []),
    getJson(`${host}/overlap/region/human/${region}?feature=cds;${base}`).catch(() => []),
  ]);
  const byParent = new Map<string, { start: number; end: number }[]>();
  for (const e of Array.isArray(exons) ? exons : []) { const p = String(e.Parent || ''); byParent.set(p, [...(byParent.get(p) || []), { start: e.start, end: e.end }]); }
  const cdsBy = new Map<string, [number, number]>();
  for (const c of Array.isArray(cdss) ? cdss : []) { const p = String(c.Parent || ''); const cur = cdsBy.get(p); cdsBy.set(p, cur ? [Math.min(cur[0], c.start), Math.max(cur[1], c.end)] : [c.start, c.end]); }
  const out: TranscriptModel[] = [];
  for (const tx of txs) {
    const tid = String(tx.id || '');
    if (!tid || Number(tx.strand) !== Number(strand)) continue;
    if ((tx.end ?? 0) < start || (tx.start ?? 0) > end) continue;
    if ((tx.start ?? start) < start - 10000 || (tx.end ?? end) > end + 10000) continue;
    const ex = (byParent.get(tid) || []).sort((a, b) => a.start - b.start);
    if (!ex.length) continue;
    const id = tx.version ? `${tid}.${tx.version}` : tid;
    const cds = cdsBy.get(tid);
    out.push({ id, name: tx.external_name || id, source: 'refseq', biotype: tx.biotype || '', start: tx.start ?? ex[0].start, end: tx.end ?? ex[ex.length - 1].end,
      strand: Number(tx.strand ?? strand), exons: ex, cds_start: cds ? cds[0] : null, cds_end: cds ? cds[1] : null, is_mane: false, is_canonical: false });
  }
  return out;
}

const structureKey = (m: TranscriptModel) => m.exons.map(e => `${e.start}-${e.end}`).join(',');

export async function getAllTranscripts(build: GenomeBuild, geneName: string, geneId?: string): Promise<AllTranscripts> {
  let data: any = null;
  if (geneId?.startsWith('ENSG')) { try { data = await lookupGene(build, geneName, geneId); } catch { data = null; } }
  if (!data) data = await lookupGene(build, geneName);
  const ensembl = ensemblModels(data);
  if (!ensembl.length) throw new Error(`No transcripts for ${geneName}`);
  const seq = String(data.seq_region_name || '');
  const chrom = seq.startsWith('chr') ? seq : `chr${seq}`;
  const strand = data.strand ?? 1;
  // Ensembl models only (ENST): identifiers stay Ensembl throughout the viewer
  const models = ensembl;
  const order: Record<string, number> = { protein_coding: 0, mRNA: 0 };
  models.sort((a, b) => Number(!a.is_mane) - Number(!b.is_mane) || (order[a.biotype] ?? 1) - (order[b.biotype] ?? 1) || Number(!a.is_canonical) - Number(!b.is_canonical) || a.id.localeCompare(b.id));
  return { gene_name: data.display_name || geneName, chrom, strand, source: 'ensembl', transcripts: models };
}

// ======================== Reference sequence (64 kb chunks, cached) ========================

const CHUNK = 1 << 16;
const chunkCache = new Map<string, Promise<string | null>>();

async function fetchChunk(build: GenomeBuild, chrom: string, idx: number): Promise<string | null> {
  const key = `${build}|${chrom}|${idx}`;
  if (!chunkCache.has(key)) {
    const bare = chrom.startsWith('chr') ? chrom.slice(3) : chrom;
    const url = `${HOSTS[build]}/sequence/region/human/${bare}:${idx * CHUNK + 1}..${(idx + 1) * CHUNK}:1?content-type=text/plain`;
    chunkCache.set(key, fetch(url, { headers: { Accept: 'text/plain' } }).then(async r => (r.ok ? (await r.text()).trim().toUpperCase() : null)).catch(() => null));
  }
  return chunkCache.get(key)!;
}

/** Reference bases for [start, end) from Ensembl, or null when unavailable. */
export async function getReference(build: GenomeBuild, chrom: string, start: number, end: number): Promise<string | null> {
  if (end <= start) return '';
  start = Math.max(0, start);
  const parts: string[] = [];
  for (let idx = Math.floor(start / CHUNK); idx <= Math.floor((end - 1) / CHUNK); idx++) {
    const c = await fetchChunk(build, chrom, idx);
    if (c == null) return null;
    parts.push(c);
  }
  const joined = parts.join('');
  const offset = start - Math.floor(start / CHUNK) * CHUNK;
  return joined.substring(offset, offset + (end - start));
}

// ======================== Genes around a region (neighbours), 500 kb chunks cached ========================

const REGION_CHUNK = 500_000;
const regionCache = new Map<string, Promise<GeneModel[]>>();

async function regionChunk(build: GenomeBuild, chrom: string, chunk: number): Promise<GeneModel[]> {
  const key = `${build}|${chrom}|${chunk}`;
  if (!regionCache.has(key)) {
    regionCache.set(key, (async () => {
      const host = HOSTS[build];
      const bare = chrom.startsWith('chr') ? chrom.slice(3) : chrom;
      const region = `${bare}:${chunk * REGION_CHUNK + 1}-${(chunk + 1) * REGION_CHUNK}`;
      const get = (feature: string) => getJson(`${host}/overlap/region/human/${region}?feature=${feature};content-type=application/json`).catch(() => []);
      const [genes, txs, exons, cdss] = await Promise.all([get('gene'), get('transcript'), get('exon'), get('cds')]);
      const exonsBy = new Map<string, { start: number; end: number }[]>();
      for (const x of Array.isArray(exons) ? exons : []) { const p = String(x.Parent || ''); exonsBy.set(p, [...(exonsBy.get(p) || []), { start: x.start, end: x.end }]); }
      const cdsBy = new Map<string, [number, number]>();
      for (const c of Array.isArray(cdss) ? cdss : []) { const p = String(c.Parent || ''); const cur = cdsBy.get(p); cdsBy.set(p, cur ? [Math.min(cur[0], c.start), Math.max(cur[1], c.end)] : [c.start, c.end]); }
      const txBy = new Map<string, any[]>();
      for (const t of Array.isArray(txs) ? txs : []) { const p = String(t.Parent || ''); txBy.set(p, [...(txBy.get(p) || []), t]); }
      const out: GeneModel[] = [];
      for (const g of Array.isArray(genes) ? genes : []) {
        const gid = String(g.id || '');
        const cands = txBy.get(gid) || [];
        if (!cands.length) continue;
        let canon = cands.filter((t: any) => t.is_canonical);
        if (!canon.length) canon = [...(cands.filter((t: any) => t.biotype === 'protein_coding').length ? cands.filter((t: any) => t.biotype === 'protein_coding') : cands)].sort((a: any, b: any) => (b.end - b.start) - (a.end - a.start));
        const t = canon[0];
        const tid = String(t.id || '');
        const ex = (exonsBy.get(tid) || []).sort((a, b) => a.start - b.start);
        if (!ex.length) continue;
        const cds = cdsBy.get(tid);
        out.push({ gene_id: gid, gene_name: g.external_name || gid, biotype: g.biotype || '', strand: Number(g.strand ?? 1), start: g.start ?? ex[0].start, end: g.end ?? ex[ex.length - 1].end,
          transcript_id: tid, is_canonical: !!t.is_canonical, exons: ex, cds_start: cds ? cds[0] : null, cds_end: cds ? cds[1] : null });
      }
      return out;
    })().catch(err => { regionCache.delete(key); throw err; }));
  }
  return regionCache.get(key)!;
}

/** Genes overlapping [start, end] (1-based inclusive) with their canonical transcript, `exclude` left out. */
export async function getRegionGenes(build: GenomeBuild, chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]> {
  if (end < start) return [];
  const seen = new Map<string, GeneModel>();
  for (let chunk = Math.floor((start - 1) / REGION_CHUNK); chunk <= Math.floor((end - 1) / REGION_CHUNK); chunk++) {
    for (const g of await regionChunk(build, chrom, chunk)) {
      if (g.end < start || g.start > end) continue;
      if (exclude && [g.gene_name.toUpperCase(), g.gene_id.toUpperCase().split('.')[0]].includes(exclude.toUpperCase())) continue;
      if (!seen.has(g.gene_id)) seen.set(g.gene_id, g);
    }
  }
  return [...seen.values()].sort((a, b) => a.start - b.start);
}

// ======================== Protein features of a translation ========================

export async function getProteinDomains(build: GenomeBuild, translationId: string): Promise<ProteinDomain[]> {
  const r = await fetch(`${HOSTS[build]}/overlap/translation/${encodeURIComponent(translationId)}?feature=protein_feature`, { headers: { Accept: 'application/json' } });
  if (!r.ok) throw new Error(`Ensembl protein features: HTTP ${r.status}`);
  const data = await r.json();
  const skip = new Set(['low_complexity', 'seg', 'ncoils', 'sigp', 'tmhmm', 'Seg']);
  return (Array.isArray(data) ? data : [])
    .filter((f: any) => !skip.has(String(f.type || '')) && f.start != null && f.end != null)
    .map((f: any) => ({ type: String(f.type || ''), id: String(f.id || f.hseqname || ''), description: String(f.description || ''), interpro: f.interpro ?? null, start: Number(f.start), end: Number(f.end) }));
}
