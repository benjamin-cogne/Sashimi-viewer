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
import { exonsInCall } from './cghExons';
import { scoreCnv, type RegionHit, type CnvGene, type AchroHit, type GeneView, type PanelView, type Score } from './cghScore';
import { cgh, type CghState } from './cghStore';
import { topoOf, type CnvTopo } from './cghTopo';

type Loaders = Omit<SvGeneLoaders, 'fallback'>;
/** the parts of a call's annotation read (and retried) on their own */
export type CghSource = 'genes' | 'dosage' | 'regions' | 'dgv' | 'bands';
export const CGH_SOURCES: CghSource[] = ['genes', 'dosage', 'regions', 'dgv', 'bands'];
/** the sources whose failure leaves the priority incomplete (DGV Gold and the bands are not in it) */
const PRIORITY_SOURCES: CghSource[] = ['genes', 'dosage', 'regions'];
export interface CnvAnn {
  state: 'loading' | 'done';
  /** RefSeq's genes, or the file's names when RefSeq could not be read */
  source: 'refseq' | 'file';
  genes: CnvGene[];
  dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>;
  achro: AchroHit[];
  /** the ClinGen dosage regions (ISCA-…) the call overlaps, with the share of each region it covers */
  regions: RegionHit[];
  /** the sources that did not answer and weigh on the priority (genes, dosage, regions), readable; kept from `failed` */
  missing: string[];
  /** each source that did not answer, its readable reason: what a manual retry asks again (PanelApp: per gene, modes[sym] === 'error') */
  failed: Partial<Record<CghSource, string>>;
  /** the sources being asked again */
  asking: CghSource[];
  /** the DGV Gold CNVs of the same type sharing DGV_COVER of the call or more, and those lying inside it whatever their share; not in the priority */
  dgv: DgvHit[];
  dgvError?: string;
  /** the call's band from the cytogenetic bands (UCSC), when the file gives none (a BED) */
  band?: string;
  /** terminal / whole arm / whole chromosome, from the chromosome length and the cytobands */
  topo?: CnvTopo;
}
/** the band of a call, as an interval table writes it: p21.1, or p21.2 - p21.1 */
function bandOf(bands: { start: number; end: number; name: string }[], s0: number, e0: number): string | undefined {
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
function dgvHits(kind: 'loss' | 'gain', s0: number, e0: number, svs: PopulationSv[]): DgvHit[] {
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
function genesOf(txs: SvTranscript[], s0: number, e0: number, kind?: 'loss' | 'gain'): CnvGene[] {
  const by = new Map<string, SvTranscript[]>();
  for (const t of txs) { if (t.end <= s0 || t.start >= e0) continue; const l = by.get(t.gene) ?? []; l.push(t); by.set(t.gene, l); }
  return [...by.entries()].map(([symbol, list]) => {
    const p = list.slice().sort((a, b) => Number(!a.status) - Number(!b.status) || (b.end - b.start) - (a.end - a.start))[0];
    return {
      symbol, tx: p.id, start: Math.min(...list.map(t => t.start)), end: Math.max(...list.map(t => t.end)),
      whole: p.start >= s0 && p.end <= e0,
      exons: list.some(t => t.exons.some(x => x.end > s0 && x.start < e0)),
      ...(kind ? { cut: exonsInCall(p, s0, e0, kind) } : {}),
    };
  }).sort((a, b) => a.start - b.start);
}

/** the AChro-Puce regions a call matches (its type), as the SV browser matches them: the sheet's critical genes, else its coordinates */
function achroHits(c: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>, genes: CnvGene[], build: GenomeBuild): AchroHit[] {
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
/** a call to read: every source, or only those `only` names (a manual retry) */
type Job = { c: CghCnv; build: GenomeBuild; only?: CghSource[] };
const queue: Job[] = [];
let running = 0;
/** annotate these calls (those not yet asked) */
export function annotate(items: { c: CghCnv; build: GenomeBuild }[], lf: (b: GenomeBuild) => Loaders) {
  loadersFor = lf;
  const known = cgh.get().ann, put: Record<string, CnvAnn> = {};
  for (const it of items) if (!known[it.c.id] && !queue.some(q => q.c.id === it.c.id)) { queue.push(it); put[it.c.id] = placeholder(); }
  // a call queued shows as waiting (…), not as not annotated, and is never asked twice
  if (Object.keys(put).length) cgh.setAnns(put);
  pump();
}
/**
 * Ask again, for these calls, only the sources that failed (PanelApp: only the genes not answered), through the same
 * queue (two calls at a time) and PanelApp's batches: what answered is kept and not asked again.
 */
export function retry(items: { c: CghCnv; build: GenomeBuild }[], sources: (CghSource | 'panelapp')[], lf: (b: GenomeBuild) => Loaders) {
  loadersFor = lf;
  const st = cgh.get(), put: Record<string, CnvAnn> = {};
  if (sources.includes('panelapp')) {
    const syms = items.flatMap(({ c }) => (st.ann[c.id]?.genes ?? []).slice(0, PANELAPP_MAX_GENES).map(g => g.symbol)).filter(x => st.modes[x] === 'error');
    wantModes([...new Set(syms)], true);
  }
  for (const it of items) {
    const a = st.ann[it.c.id];
    if (!a || a.state !== 'done') continue;
    const only = (sources.filter(x => x !== 'panelapp') as CghSource[]).filter(x => a.failed[x] && !a.asking.includes(x));
    if (!only.length) continue;
    const q = queue.find(j => j.c.id === it.c.id);
    if (q) q.only = q.only ? [...new Set([...q.only, ...only])] : undefined; else queue.push({ ...it, only });
    put[it.c.id] = { ...a, asking: [...a.asking, ...only] };
  }
  if (Object.keys(put).length) cgh.setAnns(put);
  pump();
}
function pump() {
  while (running < 2 && queue.length) {
    const it = queue.shift()!;
    running++;
    void run(it).finally(() => { running--; pump(); });
  }
}
const placeholder = (): CnvAnn => ({ state: 'loading', source: 'refseq', genes: [], dosage: new Map(), constraint: new Map(), achro: [], regions: [], missing: [], failed: {}, asking: [], dgv: [] });
const msg = (e: unknown) => String((e as Error)?.message ?? e);
type Bands = Awaited<ReturnType<NonNullable<Loaders['getCytobands']>>>;
/** the readable `missing` from `failed`: the sources weighing on the priority */
const missingOf = (f: CnvAnn['failed']) => PRIORITY_SOURCES.flatMap(x => (f[x] ? [f[x]!] : []));
/** each source read into `a`, its failure written to a.failed[source] (cleared when it answers) */
const STEPS: Record<CghSource, (c: CghCnv, ld: Loaders, a: CnvAnn, out: { bands?: Bands }) => Promise<void>> = {
  async genes(c, ld, a) {
    const s0 = c.start - 1, e0 = c.end;
    try {
      if (!ld.getTranscripts) throw new Error('no gene source for this build');
      a.genes = genesOf(await ld.getTranscripts(c.chrom, s0, e0), s0, e0, c.kind); a.source = 'refseq'; delete a.failed.genes;
    } catch (e) {
      a.source = 'file';
      a.genes = c.genes.map(symbol => ({ symbol, tx: '', start: s0, end: e0, whole: false, exons: true }));
      a.failed.genes = `RefSeq genes (${msg(e)}): the file's Gene Name instead`;
    }
  },
  async dosage(c, ld, a) {
    delete a.failed.dosage;
    if (!ld.getAnnotations || !a.genes.length) { a.dosage = new Map(); a.constraint = new Map(); return; }
    try {
      const r = await ld.getAnnotations(c.chrom, c.start - 1, c.end, a.genes.map(g => g.symbol));
      a.dosage = r.dosage; a.constraint = r.constraint;
      const un = r.notes.filter(n => /unavailable/i.test(n));
      if (un.length) a.failed.dosage = un.join('; ');
    } catch (e) { a.failed.dosage = `ClinGen dosage and gnomAD constraint (${msg(e)})`; }
  },
  // ClinGen dosage regions: those the call overlaps, with the share of the region it covers
  async regions(c, ld, a) {
    if (!ld.getDosageRegions) return;
    const s0 = c.start - 1, e0 = c.end;
    try {
      a.regions = (await ld.getDosageRegions(c.chrom, s0, e0)).flatMap(r => {
        const ov = Math.min(e0, r.end) - Math.max(s0, r.start);
        return ov > 0 && r.end > r.start ? [{ id: r.id, label: r.label, hi: r.hi, ts: r.ts, covered: Math.min(1, ov / (r.end - r.start)) }] : [];
      });
      delete a.failed.regions;
    } catch (e) { a.failed.regions = `ClinGen dosage regions (${msg(e)})`; }
  },
  // DGV Gold: the records sharing 80 % of the call or more, the call inside one counting 100 %
  async dgv(c, ld, a) {
    if (!ld.getPopulationSvs) return;
    try { a.dgv = dgvHits(c.kind, c.start - 1, c.end, await ld.getPopulationSvs(c.chrom, c.start - 1, c.end)); a.dgvError = undefined; delete a.failed.dgv; }
    catch (e) { a.dgvError = msg(e); a.failed.dgv = `DGV Gold (${a.dgvError})`; }
  },
  // a BED gives no band: read it from the cytogenetic bands (as the browser's ideogram); they also place the call (terminal, arm)
  async bands(c, ld, a, out) {
    if (!ld.getCytobands) return;
    try { const b = await ld.getCytobands(c.chrom); out.bands = b; if (!c.cytoband) a.band = bandOf(b, c.start - 1, c.end); delete a.failed.bands; }
    catch (e) { a.failed.bands = `cytogenetic bands (${msg(e)})`; }
  },
};
async function run({ c, build, only }: Job) {
  const prev = cgh.get().ann[c.id];
  // a retry starts from what answered (copied: the store's object is never changed in place); a first read from nothing
  const a: CnvAnn = only && prev?.state === 'done' ? { ...prev, failed: { ...prev.failed }, asking: [] } : { ...placeholder(), state: 'done' };
  if (!only) cgh.setAnn(c.id, placeholder());
  const steps = new Set<CghSource>(only ?? CGH_SOURCES), out: { bands?: Bands } = {}, done = new Set<CghSource>();
  let step: CghSource | null = null;
  try {
    const ld = loadersFor!(build);
    for (const x of CGH_SOURCES) {
      if (!steps.has(x)) continue;
      step = x;
      const before = x === 'genes' ? a.genes.map(g => g.symbol).join(',') + a.source : '';
      await STEPS[x](c, ld, a, out);
      done.add(x);
      // genes read again (after a failure): their dosage and constraint are asked again for the new list; their exons
      // come with them (genesOf), AChro-Puce by genes and PanelApp for the new genes follow below
      if (x === 'genes' && only && before !== a.genes.map(g => g.symbol).join(',') + a.source) steps.add('dosage');
    }
  } catch (e) {
    // an unexpected throw: what was read is kept, the step it stopped at and those not read are said, and can be retried
    for (const x of steps) if (!done.has(x)) a.failed[x] = `${x === step ? '' : 'not read: '}annotation stopped (${msg(e)})`;
  }
  try { a.achro = achroHits(c, a.genes, build); if (out.bands || !only) a.topo = topoOf(c, build, out.bands ?? null); } catch { /* left as it was */ }
  a.missing = missingOf(a.failed);
  // a retry's sources asked again since this job started stay marked as asked
  const now = cgh.get().ann[c.id];
  a.asking = only && now ? now.asking.filter(x => !steps.has(x)) : [];
  cgh.setAnn(c.id, a);
  wantModes(a.genes.slice(0, PANELAPP_MAX_GENES).map(g => g.symbol));
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
/** `again`: the genes PanelApp did not answer ('error') are asked again too (a manual retry) */
function wantModes(symbols: string[], again = false) {
  const known = cgh.get().modes;
  for (const s of symbols) if ((!(s in known) || (again && known[s] === 'error')) && !(s in pendingModes && pendingModes[s] !== 'error') && !modeQueue.includes(s)) { modeQueue.push(s); putMode(s, 'loading'); }
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
function panelView(m: GeneModes | 'loading' | 'error' | undefined): PanelView | null {
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
export function viewOf(s: CghState, c: CghCnv): { ann: CnvAnn | undefined; genes: GeneView[]; score: Score | null; waiting: number; missing: string[]; paFailed: number; /** priority sources asked again: the priority stays incomplete until they answer */ asking: number } {
  const ann = s.ann[c.id];
  if (!ann || ann.state === 'loading') return { ann, genes: [], score: null, waiting: 0, missing: [], paFailed: 0, asking: 0 };
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
  return { ann, genes, score, waiting, missing, paFailed: failed, asking: ann.asking.filter(x => PRIORITY_SOURCES.includes(x)).length };
}
