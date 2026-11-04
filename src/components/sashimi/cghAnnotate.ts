/**
 * What an array-CGH call holds, by the SV browser's own routes, on the file's build (GRCh37): the RefSeq NM_ genes of
 * its interval (UCSC ncbiRefSeqCurated; the file's Gene Name only when they cannot be read, and said so), their ClinGen
 * dosage scores and gnomAD constraint (UCSC, gnomAD), their PanelApp entries (Genomics England, every panel: rating and
 * mode of inheritance), and the AChro-Puce recurrent CNVs it matches. Calls are asked two at a time, PanelApp four genes
 * at a time; each answer is kept for the tab's life.
 */
import { geneModes, PANEL_SOURCES, type GeneModes } from '../../standalone/panelapp';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { SvGeneLoaders } from './SvGenomeBrowser';
import type { SvTranscript, DosageScore, Constraint } from './svGenes';
import { ACHRO_REGIONS, type AchroClass } from './achropuce';
import { RO_MATCH } from './svPopulation';
import type { CghCnv } from './cghIntervals';
import { scoreCnv, type CnvGene, type AchroHit, type GeneView, type PanelView, type Score } from './cghScore';
import { cgh, type CghState } from './cghStore';

type Loaders = Omit<SvGeneLoaders, 'fallback'>;
export interface CnvAnn {
  state: 'loading' | 'done';
  /** RefSeq's genes, or the file's names when RefSeq could not be read */
  source: 'refseq' | 'file';
  genes: CnvGene[];
  dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>;
  achro: AchroHit[];
  /** the sources that did not answer */
  missing: string[];
}
/** the genes PanelApp is asked for, at most, per call (the largest calls are high priority by their gene count anyway) */
export const PANELAPP_MAX_GENES = 150;
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
  cgh.setAnn(c.id, { state: 'loading', source: 'refseq', genes: [], dosage: new Map(), constraint: new Map(), achro: [], missing: [] });
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
  cgh.setAnn(c.id, { state: 'done', source, genes, dosage, constraint, achro: achroHits(c, genes, build), missing });
  wantModes(genes.slice(0, PANELAPP_MAX_GENES).map(g => g.symbol));
}

// ---- PanelApp, gene by gene ----
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
function pumpModes() {
  while (modeRunning < 4 && modeQueue.length) {
    const s = modeQueue.shift()!;
    modeRunning++;
    geneModes(PANEL_SOURCE, s).then(m => putMode(s, m), () => putMode(s, 'error')).finally(() => { modeRunning--; pumpModes(); });
  }
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
  const missing = [...ann.missing, ...(failed ? [`PanelApp (${failed} gene${failed > 1 ? 's' : ''} not answered)`] : []), ...(ann.genes.length > PANELAPP_MAX_GENES ? [`PanelApp asked for the first ${PANELAPP_MAX_GENES} genes only`] : [])];
  const score = scoreCnv({ kind: c.kind, copies: c.copies, chrom: c.chrom, genes, refseq: ann.source === 'refseq', achro: ann.achro });
  return { ann, genes, score, waiting, missing };
}
