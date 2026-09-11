/**
 * GTEx Portal API v2 client for the Sashimi "GTEx tissue" tracks (browser-side, like the Ensembl
 * and dbSNP lookups). For a gene and a tissue it returns the median junction read counts (drawn as
 * junction arcs) and the median exon read counts of the collapsed gene model, turned into a
 * reads-per-base profile that stands in for coverage. GTEx v10 (GENCODE v39) is tried first, v8
 * (GENCODE v26) when v10 has no data for the gene. hg38 only.
 */
import type { CoverageRun, JunctionArc, GtexTissue, GtexProfile } from '../components/sashimi/types';

const API = 'https://gtexportal.org/api/v2';
export type GtexDataset = 'gtex_v10' | 'gtex_v8';
const GENCODE: Record<GtexDataset, string> = { gtex_v10: 'v39', gtex_v8: 'v26' };
const DATASETS: GtexDataset[] = ['gtex_v10', 'gtex_v8'];

async function getJson(url: string): Promise<any> {
  let r: Response;
  try { r = await fetch(url, { headers: { Accept: 'application/json' } }); }
  catch (e: any) { throw new Error(`GTEx API unreachable from the browser (network or CORS): ${e?.message || e}`); }
  if (!r.ok) {
    let detail = '';
    try { const j = await r.json(); detail = j?.detail ? ` · ${typeof j.detail === 'string' ? j.detail : JSON.stringify(j.detail).slice(0, 160)}` : ''; } catch { /* no body */ }
    throw new Error(`GTEx API ${r.status} on ${url.replace(API, '').replace(/\?.*/, '')}${detail}`);
  }
  return r.json();
}
const rows = (d: any): any[] => (Array.isArray(d?.data) ? d.data : Array.isArray(d) ? d : []);
/** Every page of a paginated GTEx answer (paging_info.numberOfPages). */
async function getAllRows(url: string, perPage = 500): Promise<any[]> {
  const out: any[] = [];
  for (let page = 0; page < 40; page++) {
    const d = await getJson(`${url}&page=${page}&itemsPerPage=${perPage}`);
    out.push(...rows(d));
    const n = Number(d?.paging_info?.numberOfPages ?? 1);
    if (!(page + 1 < n)) break;
  }
  return out;
}
const num = (v: unknown): number => { const n = typeof v === 'string' ? parseFloat(v) : Number(v); return Number.isFinite(n) ? n : NaN; };
/** median-like value of an expression row whatever the field name */
const medianOf = (r: any): number => num(r?.median ?? r?.medianExpression ?? r?.value ?? r?.data);

// ---- tissues (54 in v8/v10), cached ----
let tissuesCache: Promise<GtexTissue[]> | null = null;
export function getGtexTissues(): Promise<GtexTissue[]> {
  if (!tissuesCache) {
    tissuesCache = (async () => {
      const d = await getAllRows(`${API}/dataset/tissueSiteDetail?datasetId=gtex_v10`).catch(() => getAllRows(`${API}/dataset/tissueSiteDetail?datasetId=gtex_v8`));
      const list: GtexTissue[] = d.map((t: any) => ({
        id: String(t.tissueSiteDetailId || ''), name: String(t.tissueSiteDetail || t.tissueSiteDetailId || ''), site: String(t.tissueSite || ''),
        color: `#${String(t.colorHex || '888888').replace('#', '')}`, samples: Number(t.rnaSeqSampleSummary?.totalCount ?? t.rnaSeqAndGenotypeSampleCount ?? t.rnaSeqSampleCount ?? 0),
      })).filter((t: GtexTissue) => t.id);
      list.sort((a, b) => a.site.localeCompare(b.site) || a.name.localeCompare(b.name));
      if (!list.length) throw new Error('GTEx API returned no tissue');
      return list;
    })().catch(e => { tissuesCache = null; throw e; });
  }
  return tissuesCache;
}

// ---- gene → versioned gencode id per dataset, cached ----
let lastResolveError = '';
const geneCache = new Map<string, Promise<{ gencodeId: string; dataset: GtexDataset } | null>>();
async function resolveGene(geneName: string, geneId: string | undefined, dataset: GtexDataset): Promise<{ gencodeId: string; dataset: GtexDataset } | null> {
  const key = `${dataset}|${geneId || geneName}`;
  if (!geneCache.has(key)) {
    geneCache.set(key, (async () => {
      const queries = [geneId?.split('.')[0], geneName].filter((x): x is string => !!x);
      for (const q of queries) {
        const d = await getJson(`${API}/reference/gene?geneId=${encodeURIComponent(q)}&gencodeVersion=${GENCODE[dataset]}&genomeBuild=GRCh38%2Fhg38&page=0&itemsPerPage=50`).catch(e => { lastResolveError = e?.message || String(e); return null; });
        const hits = rows(d).filter((g: any) => g.gencodeId);
        const best = hits.find((g: any) => geneId && String(g.gencodeId).split('.')[0] === geneId.split('.')[0])
          || hits.find((g: any) => String(g.geneSymbol || '').toUpperCase() === geneName.toUpperCase()) || hits[0];
        if (best) return { gencodeId: String(best.gencodeId), dataset };
      }
      return null;
    })().catch(e => { geneCache.delete(key); throw e; }));
  }
  return geneCache.get(key)!;
}

// ---- collapsed exons of the gene model (1-based) per dataset, cached ----
interface RefExon { id: string; start: number; end: number }
interface RefExons {
  /** exons in transcription order (5′ → 3′) */
  ordered: RefExon[];
  /** 'collapsed' = GTEx collapsed gene model (ids ENSG…_N, what the exon expression rows use);
   *  'merged' = union of all transcript exons merged into non-overlapping intervals (ids assigned by rank) */
  source: 'collapsed' | 'merged';
  strand: number;
}
const exonCache = new Map<string, Promise<RefExons>>();

const strandOf = (v: unknown): number => (v === '-' || v === -1 || v === '-1' ? -1 : 1);

/**
 * Exons of the GTEx collapsed gene model, which the exon expression rows (ENSG…_N) refer to.
 * First choice: the collapsedGeneModelExon endpoint. Otherwise the transcript exons of the gene
 * (reference/exon lists every exon of every transcript, Ensembl ENSE ids) merged into non-overlapping
 * intervals, the way GTEx's collapse_annotation builds its model; those are matched by rank.
 */
function geneExons(gencodeId: string, dataset: GtexDataset): Promise<RefExons> {
  const key = `${dataset}|${gencodeId}`;
  if (!exonCache.has(key)) {
    exonCache.set(key, (async (): Promise<RefExons> => {
      const gc = encodeURIComponent(gencodeId);
      const toExon = (e: any): RefExon => ({ id: String(e.exonId ?? e.id ?? ''), start: num(e.start ?? e.exonStart), end: num(e.end ?? e.exonEnd) });
      const valid = (e: RefExon) => Number.isFinite(e.start) && Number.isFinite(e.end) && e.end >= e.start;
      const order = (list: RefExon[], strand: number) => [...list].sort((a, b) => (strand < 0 ? b.start - a.start : a.start - b.start));
      for (const q of [`datasetId=${dataset}`, `gencodeVersion=${GENCODE[dataset]}&genomeBuild=GRCh38%2Fhg38`]) {
        const rows = await getAllRows(`${API}/reference/collapsedGeneModelExon?gencodeId=${gc}&${q}`).catch(() => [] as any[]);
        const list = rows.map(toExon).filter(valid);
        if (list.length && list.some(e => /_\d+$/.test(e.id))) return { ordered: order(list, strandOf(rows[0]?.strand)), source: 'collapsed', strand: strandOf(rows[0]?.strand) };
      }
      const q2 = `gencodeVersion=${GENCODE[dataset]}&genomeBuild=GRCh38%2Fhg38`;
      const [rows, trRows] = await Promise.all([
        getAllRows(`${API}/reference/exon?gencodeId=${gc}&${q2}`),
        getAllRows(`${API}/reference/transcript?gencodeId=${gc}&${q2}`).catch(() => [] as any[]),
      ]);
      const strand = strandOf(rows[0]?.strand);
      // GTEx's collapse leaves retained-intron transcripts out (their exons would bridge neighbouring model exons)
      const trKey = (r: any) => Object.keys(r || {}).find(k => /^transcript(_?id)?$/i.test(k) || /transcriptId/i.test(k));
      const excluded = new Set(trRows.filter(r => JSON.stringify(r).toLowerCase().includes('retained_intron')).map(r => String(r[trKey(r) || 'transcriptId'] ?? '')).filter(Boolean));
      const kept = rows.filter(r => { const k = trKey(r); return !k || !excluded.has(String(r[k])); });
      const all = (kept.length ? kept : rows).map(toExon).filter(valid).sort((a, b) => a.start - b.start);
      const merged: RefExon[] = [];
      for (const e of all) {
        const last = merged[merged.length - 1];
        if (last && e.start <= last.end) last.end = Math.max(last.end, e.end);
        else merged.push({ id: '', start: e.start, end: e.end });
      }
      // ids deliberately unlike ENSG…_N so that only the rank logic (numbering base inferred from the rows) applies
      const ordered = order(merged, strand).map((e, i) => ({ ...e, id: `rank:${i}` }));
      return { ordered, source: 'merged', strand };
    })().catch(e => { exonCache.delete(key); throw e; }));
  }
  return exonCache.get(key)!;
}

/** "chr7_100401_104000" (1-based intron start and end, inclusive) → 0-based half-open */
function parseJunction(id: string): { start: number; end: number } | null {
  // "chr7_100401_104000", "chr7:100401-104000", "7_100401_104000" … (1-based intron start/end, inclusive)
  const m = String(id).match(/^(?:chr)?[0-9XYMxym]+[_:](\d+)[_\-:](\d+)(?:[_:][+\-.])?$/);
  return m ? { start: Number(m[1]) - 1, end: Number(m[2]) } : null;
}

/** Intron of a junction row: junctionId in any known spelling, else explicit start/end fields, else any string field that looks like a locus. */
function junctionOf(r: any): { start: number; end: number } | null {
  for (const k of ['junctionId', 'id', 'junction', 'name']) { const p = r?.[k] != null ? parseJunction(String(r[k])) : null; if (p) return p; }
  const st = num(r?.junctionStart ?? r?.start ?? r?.intronStart), en = num(r?.junctionEnd ?? r?.end ?? r?.intronEnd);
  if (Number.isFinite(st) && Number.isFinite(en) && en > st) return { start: st - 1, end: en };
  for (const v of Object.values(r || {})) { if (typeof v === 'string') { const p = parseJunction(v); if (p) return p; } }
  return null;
}

/** Median of an expression row: a key containing "median", else value/count/expression, else the only numeric field that is not a coordinate. */
function medianAny(r: any): number {
  const m = medianOf(r);
  if (Number.isFinite(m)) return m;
  const keys = Object.keys(r || {});
  const k1 = keys.find(k => /median/i.test(k) && Number.isFinite(num(r[k])));
  if (k1) return num(r[k1]);
  const k2 = keys.find(k => /^(value|count|expression|tpm|reads?)$/i.test(k) && Number.isFinite(num(r[k])));
  if (k2) return num(r[k2]);
  const numeric = keys.filter(k => !/start|end|pos|number|page|length|strand|chrom/i.test(k) && typeof r[k] === 'number');
  return numeric.length === 1 ? num(r[numeric[0]]) : NaN;
}

const profileCache = new Map<string, Promise<GtexProfile>>();

/** Median junction and exon read counts of a tissue for the gene, as arcs and a reads-per-base coverage profile. */
export function getGtexProfile(geneName: string, geneId: string | undefined, tissue: GtexTissue): Promise<GtexProfile> {
  const key = `${geneId || geneName}|${tissue.id}`;
  if (!profileCache.has(key)) {
    profileCache.set(key, (async () => {
      let lastError = '';
      let resolved = false;
      for (const dataset of DATASETS) {
        const gene = await resolveGene(geneName, geneId, dataset);
        if (!gene) continue;
        resolved = true;
        const q = `gencodeId=${encodeURIComponent(gene.gencodeId)}&datasetId=${dataset}&tissueSiteDetailId=${encodeURIComponent(tissue.id)}`;
        const [jrows, erows, refExons, grows] = await Promise.all([
          getAllRows(`${API}/expression/medianJunctionExpression?${q}`).catch(e => { lastError = e.message; return [] as any[]; }),
          getAllRows(`${API}/expression/medianExonExpression?${q}`).catch(e => { lastError = e.message; return [] as any[]; }),
          geneExons(gene.gencodeId, dataset).catch(e => { lastError = e.message; return { ordered: [] as RefExon[], source: 'merged' as const, strand: 1 }; }),
          getAllRows(`${API}/expression/medianGeneExpression?${q}`).catch(() => [] as any[]),
        ]);
        const tpmRow = grows.find((r: any) => Number.isFinite(medianAny(r)));
        const tpm = tpmRow ? medianAny(tpmRow) : null;
        if (!jrows.length && !erows.length) continue;   // no data in this release: try the next
        const junctions: JunctionArc[] = jrows.map((r: any) => { const p = junctionOf(r); const m = medianAny(r); return p && m > 0 ? { ...p, count: m } : null; }).filter((x): x is JunctionArc => !!x).sort((a, b) => a.start - b.start);
        // Exon rows (ENSG…_N) join the collapsed model by id (gene versions ignored) or by suffix; with a merged
        // model the suffix is a rank in transcription order (numbering base taken from the smallest suffix seen).
        const exons = refExons.ordered;
        const norm = (id: string) => id.replace(/\.\d+(?=_|$)/g, '').toLowerCase();
        const byExon = new Map(exons.map(e => [norm(e.id), e]));
        const bySuffix = new Map(exons.map(e => [String(e.id).split('_').pop() || '', e]));
        const suffixOf = (r: any): number => { const m = String(r?.exonId ?? r?.id ?? '').match(/_(\d+)$/); return m ? Number(m[1]) : num(r?.exonNumber); };
        const suffixes = erows.map(suffixOf).filter(Number.isFinite);
        const base = suffixes.length ? Math.min(...suffixes) : 0;
        const locate = (r: any, i: number): { start: number; end: number } | null => {
          const id = String(r?.exonId ?? r?.id ?? '');
          let e = byExon.get(norm(id));
          if (!e && refExons.source === 'collapsed') e = bySuffix.get(id.split('_').pop() || '') || (r?.exonNumber != null ? bySuffix.get(String(r.exonNumber)) : undefined);
          if (!e && refExons.source === 'merged') {
            const k = suffixOf(r);
            e = Number.isFinite(k) && exons.length === erows.length ? exons[k - base] : Number.isFinite(k) && Math.abs(exons.length - erows.length) <= 2 ? exons[k - base] : erows.length === exons.length ? exons[i] : undefined;
          }
          if (e) return { start: e.start - 1, end: e.end };
          const st = num(r?.start ?? r?.exonStart), en = num(r?.end ?? r?.exonEnd);
          return Number.isFinite(st) && Number.isFinite(en) && en > st ? { start: st - 1, end: en } : null;
        };
        const ex = erows.map((r: any, i: number) => { const e = locate(r, i); const m = medianAny(r); return e && Number.isFinite(m) ? { start: e.start, end: e.end, median: m } : null; })
          .filter((x): x is { start: number; end: number; median: number } => !!x && x.end > x.start).sort((a, b) => a.start - b.start);
        if (!junctions.length && !ex.length) {
          // the API answered but nothing could be interpreted: say what came back so the parser can be adjusted
          const sample = (r: any) => (r ? Object.keys(r).slice(0, 8).join(',') : 'none');
          console.warn('[sashimi] GTEx rows not understood', { junction: jrows[0], exon: erows[0], exonRef: exons[0] });
          throw new Error(`GTEx ${dataset} answered (${jrows.length} junction rows, ${erows.length} exon rows, ${exons.length} model exons) but nothing could be parsed · junction keys: ${sample(jrows[0])} · exon keys: ${sample(erows[0])}`);
        }
        // coverage profile: median reads of the exon spread over its length (reads per base); introns at 0
        const coverage: CoverageRun[] = [];
        let cursor = ex.length ? ex[0].start : 0;
        for (const e of ex) {
          if (e.start > cursor) coverage.push({ start: cursor, end: e.start, depth: 0 });
          const s = Math.max(e.start, cursor);
          if (e.end > s) coverage.push({ start: s, end: e.end, depth: e.median / (e.end - e.start) });
          cursor = Math.max(cursor, e.end);
        }
        const brief = (r: any) => (r ? JSON.stringify(r).slice(0, 140) : 'none');
        if (erows.length && !ex.length) console.warn('[sashimi] GTEx exon rows could not be placed', { exonRow: erows[0], referenceExon: exons[0], exonRows: erows.length, referenceExons: exons.length, source: refExons.source });
        const warning = lastError ? `partial: ${lastError}`
          : !junctions.length ? `no junction data (${jrows.length} rows, first: ${brief(jrows[0])})`
          : !ex.length ? `no exon data: ${erows.length} expression rows (first: ${brief(erows[0])}) vs ${exons.length} ${refExons.source} model exons (first: ${brief(exons[0])})`
          : undefined;
        return { dataset, gencodeId: gene.gencodeId, tissue, junctions, exons: ex, coverage, unit: String(jrows[0]?.unit || erows[0]?.unit || 'median read count'), warning, tpm };
      }
      throw new Error(lastError ? `${lastError} (GTEx v10 and v8)` : resolved ? `no GTEx expression data for ${geneName} in ${tissue.name} (v10 and v8)` : `${geneName} not found in GTEx (v39 / v26 gene models)${lastResolveError ? ` · ${lastResolveError}` : ''}`);
    })().catch(e => { profileCache.delete(key); throw e; }));
  }
  return profileCache.get(key)!;
}

/** Default favourite tissues for the picker (most used in diagnostic RNA-seq). */
export const GTEX_DEFAULT_FAVOURITES = ['Whole_Blood', 'Cells_Cultured_fibroblasts', 'Muscle_Skeletal', 'Brain_Cortex'];
