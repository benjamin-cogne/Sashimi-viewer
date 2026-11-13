/**
 * What an array-CGH call holds, by the SV browser's own routes, on the file's build (GRCh37): the RefSeq NM_ genes of
 * its interval (UCSC ncbiRefSeqCurated; the file's Gene Name only when they cannot be read, and said so), their ClinGen
 * dosage scores and gnomAD constraint (UCSC, gnomAD), their PanelApp entries (Genomics England, every panel: rating and
 * mode of inheritance), and the AChro-Puce recurrent CNVs it matches. Calls are asked two at a time, PanelApp 25 genes
 * a request, one request at a time (a batch not answered asked again twice); each answer is kept for the tab's life.
 */
import { PANEL_SOURCES, classifyMoi, summarizeModes, type GeneModes, type GeneMode } from '../../standalone/panelapp';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { SvGeneLoaders } from './SvGenomeBrowser';
import type { SvTranscript, DosageScore, Constraint } from './svGenes';
import { ACHRO_REGIONS, type AchroClass } from './achropuce';
import { RO_MATCH, type PopulationSv } from './svPopulation';
import type { CghCnv } from './cghIntervals';
import { scoreCnv, type RegionHit, type CnvGene, type AchroHit, type GeneView, type PanelView, type Score } from './cghScore';
import { cgh, type CghState } from './cghStore';
import { topoOf, type CnvTopo } from './cghTopo';

type Loaders = Omit<SvGeneLoaders, 'fallback'>;
export interface CnvAnn {
  state: 'loading' | 'done';
  /** RefSeq's genes, or the file's names when RefSeq could not be read */
  source: 'refseq' | 'file';
  genes: CnvGene[];
  dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>;
  achro: AchroHit[];
  /** the ClinGen dosage regions (ISCA-…) the call overlaps, with the share of each region it covers */
  regions: RegionHit[];
  /** the sources that did not answer */
  missing: string[];
  /** the DGV Gold CNVs of the same type sharing DGV_COVER of the call or more, and those lying inside it whatever their share; not in the priority */
  dgv: DgvHit[];
  dgvError?: string;
  /** the call's band from the cytogenetic bands (UCSC), when the file gives none (a BED) */
  band?: string;
  /** terminal / whole arm / whole chromosome, from the chromosome length and the cytobands */
  topo?: CnvTopo;
}
/** the band of a call, as an interval table writes it: p21.1, or p21.2 - p21.1 */
export function bandOf(bands: { start: number; end: number; name: string }[], s0: number, e0: number): string | undefined {
  const at = (pos: number) => bands.find(b => b.start <= pos && pos < b.end)?.name;
  const a = at(s0), z = at(Math.max(s0, e0 - 1));
  return !a || !z ? undefined : a === z ? a : `${a} - ${z}`;
}
/** a call's band: the file's, else the one read from the bands */
export const bandOfCall = (c: CghCnv, ann?: CnvAnn) => c.cytoband || ann?.band || '';
/** how a DGV Gold CNV lies against the call: the call inside it, the call holding it, or the two overlapping on one side */
export type DgvRelation = 'within' | 'contains' | 'overlaps';
export interface DgvHit { id: string; type: string; rel: DgvRelation; cover: number; freq: string | null; samples: number | null; studies: number | null }
/** the share of the call a DGV Gold CNV must share with it (a call wholly inside one: 100 %) */
export const DGV_COVER = 0.8;
/**
 * The DGV Gold CNVs of the call's type (loss, gain, or both) sharing DGV_COVER of the call or more: the common stretch over
 * the call's length, with the record's inner or outer extent, the one sharing more; and how they lie against each other.
 * A record lying inside the call (one of its extents) is kept whatever its share, as 'contains'.
 */
export function dgvHits(kind: 'loss' | 'gain', s0: number, e0: number, svs: PopulationSv[]): DgvHit[] {
  const len = Math.max(1, e0 - s0);
  const out: DgvHit[] = [];
  for (const v of svs) {
    if (!(v.type === kind || v.type === 'both')) continue;
    const ext = [v.inner, v.outer].map(([a, b]) => ({ a, b, cover: Math.max(0, Math.min(e0, b) - Math.max(s0, a)) / len }));
    let best = ext.slice().sort((x, y) => y.cover - x.cover || (x.b - x.a) - (y.b - y.a))[0];
    if (best.cover < DGV_COVER) {
      // under the share: kept only when it lies inside the call (its larger extent that does)
      const inside = ext.filter(x => x.a >= s0 && x.b <= e0 && x.b > x.a).sort((x, y) => y.cover - x.cover)[0];
      if (!inside) continue;
      best = inside;
    }
    const rel: DgvRelation = best.a <= s0 && best.b >= e0 ? 'within' : s0 <= best.a && e0 >= best.b ? 'contains' : 'overlaps';
    out.push({ id: v.id, type: v.type, rel, cover: best.cover, freq: v.freqText, samples: v.samples, studies: v.studies });
  }
  return out.sort((a, b) => b.cover - a.cover);
}
/** every gene of a call is asked of PanelApp (one request per gene, cached); a bound against a whole-arm call only */
export const PANELAPP_MAX_GENES = 1000;
export const PANEL_SOURCE = PANEL_SOURCES[0];

/** the genes of a window from its RefSeq models: each gene's principal model (MANE / RefSeq Select, else the longest) */
export function genesOf(txs: SvTranscript[], s0: number, e0: number): CnvGene[] {
  const by = new Map<string, SvTranscript[]>();
  for (const t of txs) { if (t.end <= s0 || t.start >= e0) continue; const l = by.get(t.gene) ?? []; l.push(t); by.set(t.gene, l); }
  return [...by.entries()].map(([symbol, list]) => {
    const p = list.slice().sort((a, b) => Number(!a.status) - Number(!b.status) || (b.end - b.start) - (a.end - a.start))[0];
    return {
      symbol, tx: p.id, start: Math.min(...list.map(t => t.start)), end: Math.max(...list.map(t => t.end)),
      whole: p.start >= s0 && p.end <= e0,
      exons: list.some(t => t.exons.some(x => x.end > s0 && x.start < e0)),
    };
  }).sort((a, b) => a.start - b.start);
}

/** the AChro-Puce regions a call matches (its type), as the SV browser matches them: the sheet's critical genes, else its coordinates */
export function achroHits(c: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>, genes: CnvGene[], build: GenomeBuild): AchroHit[] {
  const t = c.kind === 'loss' ? 'del' : 'dup';
  const out: AchroHit[] = [];
  for (const r of ACHRO_REGIONS) {
    const rows = r.rows.filter(x => x.type === t);
    if (!rows.length) continue;
    const cls = [...new Set(rows.flatMap(x => [x.cls, ...(x.also ? [x.also.cls] : [])]))] as AchroClass[];
    const rule = r.typed?.[t];
    let how: string | null = null;
    if (rule?.genes?.length) {
      const inCall = rule.genes.filter(g => genes.some(x => x.symbol === g));
      if (rule.match === 'overlaps' && inCall.length) how = `takes ${inCall.join(', ')}`;
      if (rule.match === 'contains' && rule.genes.every(g => genes.some(x => x.symbol === g && x.whole))) how = `spans ${rule.genes.join(', ')}`;
    } else if (!r.coords?.some(x => x.build === build)) {
      // a region placed by its genes only (2p16.3, NRXN1): the call takes exons of one of them
      const hit = r.genes.filter(g => genes.some(x => x.symbol === g && x.exons));
      if (hit.length) how = `takes exons of ${hit.join(', ')}`;
    } else {
      const co = r.coords?.find(x => x.build === build && x.chrom === c.chrom);
      if (co) {
        const ov = Math.max(0, Math.min(c.end, co.end) - Math.max(c.start, co.start) + 1), len = co.end - co.start + 1;
        if (rule?.match === 'overlaps') how = ov > 0 ? `overlaps ${Math.round((ov / len) * 100)} % of it` : null;
        else if (rule?.match === 'contains') how = ov >= len ? 'spans its smallest region of overlap' : null;
        else { const ro = ov / Math.max(c.end - c.start + 1, len); how = ro >= RO_MATCH ? `${Math.round(ro * 100)} % reciprocal overlap` : null; }
      }
    }
    if (how) out.push({ id: r.id, name: r.name, type: t, cls, how });
  }
  return out;
}

// ---- the queue of calls ----
let loadersFor: ((b: GenomeBuild) => Loaders) | null = null;
const queue: { c: CghCnv; build: GenomeBuild }[] = [];
let running = 0;
/** annotate these calls (those not yet asked) */
export function annotate(items: { c: CghCnv; build: GenomeBuild }[], lf: (b: GenomeBuild) => Loaders) {
  loadersFor = lf;
  const known = cgh.get().ann;
  for (const it of items) if (!known[it.c.id] && !queue.some(q => q.c.id === it.c.id)) queue.push(it);
  pump();
}
function pump() {
  while (running < 2 && queue.length) {
    const it = queue.shift()!;
    running++;
    void run(it).finally(() => { running--; pump(); });
  }
}
const msg = (e: unknown) => String((e as Error)?.message ?? e);
async function run({ c, build }: { c: CghCnv; build: GenomeBuild }) {
  cgh.setAnn(c.id, { state: 'loading', source: 'refseq', genes: [], dosage: new Map(), constraint: new Map(), achro: [], regions: [], missing: [], dgv: [] });
  const ld = loadersFor!(build), s0 = c.start - 1, e0 = c.end, missing: string[] = [];
  let genes: CnvGene[] = [], source: CnvAnn['source'] = 'refseq';
  try {
    if (!ld.getTranscripts) throw new Error('no gene source for this build');
    genes = genesOf(await ld.getTranscripts(c.chrom, s0, e0), s0, e0);
  } catch (e) {
    source = 'file';
    genes = c.genes.map(symbol => ({ symbol, tx: '', start: s0, end: e0, whole: false, exons: true }));
    missing.push(`RefSeq genes (${msg(e)}): the file's Gene Name instead`);
  }
  let dosage = new Map<string, DosageScore>(), constraint = new Map<string, Constraint>();
  if (ld.getAnnotations && genes.length) {
    try {
      const a = await ld.getAnnotations(c.chrom, s0, e0, genes.map(g => g.symbol));
      dosage = a.dosage; constraint = a.constraint;
      missing.push(...a.notes.filter(n => /unavailable/i.test(n)));
    } catch (e) { missing.push(`ClinGen dosage and gnomAD constraint (${msg(e)})`); }
  }
  // ClinGen dosage regions: those the call overlaps, with the share of the region it covers
  let regions: RegionHit[] = [];
  if (ld.getDosageRegions) {
    try {
      regions = (await ld.getDosageRegions(c.chrom, s0, e0)).flatMap(r => {
        const ov = Math.min(e0, r.end) - Math.max(s0, r.start);
        return ov > 0 && r.end > r.start ? [{ id: r.id, label: r.label, hi: r.hi, ts: r.ts, covered: Math.min(1, ov / (r.end - r.start)) }] : [];
      });
    } catch (e) { missing.push(`ClinGen dosage regions (${msg(e)})`); }
  }
  // DGV Gold: the records sharing 80 % of the call or more, the call inside one counting 100 %
  let dgv: DgvHit[] = [], dgvError: string | undefined;
  if (ld.getPopulationSvs) {
    try { dgv = dgvHits(c.kind, s0, e0, await ld.getPopulationSvs(c.chrom, s0, e0)); } catch (e) { dgvError = msg(e); }
  }
  // a BED gives no band: read it from the cytogenetic bands (as the browser's ideogram)
  let band: string | undefined;
  let bands: Awaited<ReturnType<NonNullable<Loaders['getCytobands']>>> | null = null;
  if (ld.getCytobands) { try { bands = await ld.getCytobands(c.chrom); } catch { /* the band stays unknown */ } }
  if (!c.cytoband && bands) band = bandOf(bands, s0, e0);
  cgh.setAnn(c.id, { state: 'done', source, genes, dosage, constraint, achro: achroHits(c, genes, build), regions, missing, dgv, dgvError, band, topo: topoOf(c, build, bands) });
  wantModes(genes.slice(0, PANELAPP_MAX_GENES).map(g => g.symbol));
}

// ---- PanelApp, by batches of genes ----
const modeQueue: string[] = [];
let modeRunning = 0;
let pendingModes: Record<string, GeneModes | 'loading' | 'error'> = {};
let flushTimer: ReturnType<typeof setTimeout> | null = null;
/** the answers are written to the store a few at a time (one redraw for many genes) */
const putMode = (sym: string, m: GeneModes | 'loading' | 'error') => {
  pendingModes[sym] = m;
  if (!flushTimer) flushTimer = setTimeout(() => { flushTimer = null; const p = pendingModes; pendingModes = {}; cgh.setModes(p); }, 150);
};
function wantModes(symbols: string[]) {
  const known = cgh.get().modes;
  for (const s of symbols) if (!(s in known) && !(s in pendingModes) && !modeQueue.includes(s)) { modeQueue.push(s); putMode(s, 'loading'); }
  pumpModes();
}
/**
 * PanelApp is asked for 25 genes in one request (`/genes/?entity_name=A,B,…`, every page), one request at a time with
 * 300 ms between two: it refuses bursts (HTTP 429, which the browser sees as a network error). A batch not answered is
 * asked again twice, after 5 s then 15 s, its genes then marked not answered. A gene in no panel has no entry.
 */
const MODE_BATCH = 25, MODE_GAP_MS = 300, MODE_RETRY_MS = [5000, 15000];
let modeLast = 0, modeTimer: ReturnType<typeof setTimeout> | null = null;
const conf = (v: unknown): 1 | 2 | 3 | null => { const n = Number(v); return n >= 3 ? 3 : n === 2 ? 2 : n === 1 ? 1 : null; };
async function batchModes(symbols: string[]): Promise<Map<string, GeneModes>> {
  const by = new Map<string, Map<string, GeneMode>>(symbols.map(x => [x.toUpperCase(), new Map()]));
  let url: string | null = `${PANEL_SOURCE.base}/genes/?entity_name=${symbols.map(encodeURIComponent).join(',')}`;
  for (let guard = 0; url && guard < 50; guard++) {
    const r: Response = await fetch(url, { headers: { Accept: 'application/json' } });
    if (!r.ok) throw new Error(`${PANEL_SOURCE.label}: HTTP ${r.status}`);
    const j: any = await r.json();
    for (const g of j?.results ?? []) {
      const sym = String(g?.entity_name ?? '').toUpperCase(), c = conf(g?.confidence_level), pid = String(g?.panel?.id ?? '');
      const m = by.get(sym);
      if (!m || !c || !pid || m.has(pid) || (g?.entity_type && g.entity_type !== 'gene')) continue;
      m.set(pid, { panel: String(g?.panel?.name ?? `panel ${pid}`), panelId: pid, confidence: c, moi: String(g?.mode_of_inheritance ?? ''), cls: classifyMoi(g?.mode_of_inheritance),
        phenotypes: Array.isArray(g?.phenotypes) ? g.phenotypes.map(String) : [] });
    }
    // the next page, on the same (https) address
    url = typeof j?.next === 'string' ? j.next.replace(/^http:\/\//, 'https://') : null;
  }
  return new Map(symbols.map(x => [x, summarizeModes(x, PANEL_SOURCE.label, [...by.get(x.toUpperCase())!.values()])]));
}
function pumpModes() {
  if (modeRunning || !modeQueue.length) return;
  const wait = modeLast + MODE_GAP_MS - Date.now();
  if (wait > 0) { if (!modeTimer) modeTimer = setTimeout(() => { modeTimer = null; pumpModes(); }, wait); return; }
  const batch = modeQueue.splice(0, MODE_BATCH);
  modeRunning = 1; modeLast = Date.now();
  const ask = (n: number) => batchModes(batch).then(ms => { for (const [k, m] of ms) putMode(k, m); }, () => {
    if (n >= MODE_RETRY_MS.length) { for (const k of batch) putMode(k, 'error'); return; }
    return new Promise<void>(res => setTimeout(() => res(ask(n + 1)), MODE_RETRY_MS[n]));
  });
  void ask(0).finally(() => { modeRunning = 0; modeLast = Date.now(); pumpModes(); });
}

/** a gene's PanelApp entries as the score reads them; null while asked or when PanelApp failed */
export function panelView(m: GeneModes | 'loading' | 'error' | undefined): PanelView | null {
  if (!m || m === 'loading' || m === 'error') return null;
  const best = m.entries.reduce<0 | 1 | 2 | 3>((b, e) => (e.confidence > b ? e.confidence : b), 0);
  return {
    best,
    green: [...new Set(m.entries.filter(e => e.confidence === 3).map(e => e.cls))],
    amber: [...new Set(m.entries.filter(e => e.confidence === 2).map(e => e.cls))],
    diseases: m.diseases.map(d => d.name),
  };
}

/** a call's genes with their scores, its priority, and what is still asked or missing */
export function viewOf(s: CghState, c: CghCnv): { ann: CnvAnn | undefined; genes: GeneView[]; score: Score | null; waiting: number; missing: string[] } {
  const ann = s.ann[c.id];
  if (!ann || ann.state === 'loading') return { ann, genes: [], score: null, waiting: 0, missing: [] };
  const genes: GeneView[] = ann.genes.map((g, i) => {
    const d = ann.dosage.get(g.symbol), k = ann.constraint.get(g.symbol);
    return { g, hi: d?.hi ?? null, ts: d?.ts ?? null, loeuf: k?.loeuf ?? null, pli: k?.pli ?? null, version: k?.version ?? null,
      pa: i < PANELAPP_MAX_GENES ? panelView(s.modes[g.symbol]) : null,
      paState: i >= PANELAPP_MAX_GENES ? 'not asked' : s.modes[g.symbol] === 'error' ? 'failed' : !s.modes[g.symbol] || s.modes[g.symbol] === 'loading' ? 'asked' : 'done' };
  });
  const asked = ann.genes.slice(0, PANELAPP_MAX_GENES).map(g => s.modes[g.symbol]);
  const waiting = asked.filter(m => m === 'loading' || m === undefined).length;
  const failed = asked.filter(m => m === 'error').length;
  const missing = [...ann.missing, ...(failed ? [`PanelApp (${failed} gene${failed > 1 ? 's' : ''} not answered)`] : []), ...(ann.genes.length > PANELAPP_MAX_GENES ? [`PanelApp asked for ${PANELAPP_MAX_GENES} of the ${ann.genes.length} genes`] : [])];
  const score = scoreCnv({ kind: c.kind, copies: c.copies, chrom: c.chrom, sex: s.sex[c.sampleKey] ?? null, sizeKb: c.sizeKb, genes, refseq: ann.source === 'refseq', achro: ann.achro, regions: ann.regions });
  return { ann, genes, score, waiting, missing };
}
