/**
 * The SV browser: a genome browser over a whole chromosome, with the rearrangement when there is one. On top, the
 * chromosome's ideogram (its cytogenetic bands) with the part in view, a ruler, the genes, the population CNVs (DGV Gold
 * Standard), then the reference cut into the event's segments, a ribbon from every segment to wherever it went
 * (twisted when inverted, doubled when duplicated, none when lost), and the rebuilt chromosome with its genes in their
 * new coordinates. The rebuilt chromosome follows the reference: under the place in view, its own copy of it (it slides
 * where the two part ways, an inverted segment running the other way).
 *
 * Anywhere on the chromosome: data come by tiles of the span in view (browserTiles.ts) and stay in a bounded cache;
 * the level of detail follows the zoom (exons and strands up close, gene bars at a few Mb, a density profile beyond).
 * Only what is in view is drawn, a transcript as three paths. Panning and zooming update the drawing once per frame;
 * a fling goes on after a quick drag; the buttons and the ideogram fly there. Genes are coloured by their gnomAD LoF
 * constraint or by what the rearrangement does to them; a cut gene shows torn ends, a fusion its two partners joined
 * under one label. Chromatin domains from a BED the user loads are drawn over the reference and under the rebuilt
 * chromosome, new domains outlined.
 */
import React, { useEffect, useLayoutEffect, useMemo, useReducer, useRef, useState } from 'react';
import { reconstruct, type SvReconstruction } from './svModel';
import type { Cytoband } from './svNomenclature';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';
import { geneStory, dosageReading, toDer, LOEUF_INTOLERANT, LOEUF_CONSTRAINED, LOEUF_CLASS, loeufClass, PLI_INTOLERANT, NEAR_BREAKPOINT_DEFAULT,
  type SvGene, type SvTranscript, type GeneAnnotations, type GenePanel, type Domain, type GeneImpact, type GeneStory, type Chain, type Fragment, type DerLayout, type DerPiece, type DosageScore, type Constraint, type ImpactKind } from './svGenes';
import { IMPACT, KIND_COLOR } from './svColors';
import { cnvEvents, matchesOf, RO_MATCH, type PopulationSv, type PopMatch, type CnvEvent } from './svPopulation';
import { TileCache, tileSizeFor } from './browserTiles';
import { niceStep } from './geometry';
import { PANEL_SOURCES, listPanels, loadPanel, parsePanelFile, type PanelSummary, type LoadedPanel } from '../../standalone/panelapp';

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const MONO = 'ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
/**
 * The drawing's width in px (one unit per pixel, so a wider screen shows more of the genome at the same text size),
 * the plot's left and right edges, its width and centre. Set by the browser from its box at each render (one browser
 * drawn at a time).
 */
const X0 = 112;
let VW = 1200, X1 = VW - 14, PW = X1 - X0, XC = (X0 + X1) / 2;
const setWidth = (w: number) => { VW = w; X1 = w - 14; PW = X1 - X0; XC = (X0 + X1) / 2; };
const DGV_LANE_H = 9, DGV_MAX_LANES = 8;
const DGV_COLOR = { loss: '#dc2626', gain: '#2563eb', both: '#7c3aed', unknown: '#64748b' } as const;
const TX_H = 13, LABEL_H = 13, LANE_GAP = 5, MAX_LANES = 8, TAD_H = 58, BAND_H = 20, RIB_H = 118;
const IDEO_Y = 22, IDEO_H = 13, LENS_H = 22, RULER_H = 24;
/** the deepest zoom (bp across the plot), and the spans past which genes become bars, then a density profile */
const MIN_SPAN = 120, LOD_FULL = 3e6, LOD_BARS = 30e6;
/** population CNVs are drawn up to this span (beyond, a chromosome's thousands of records would be a smear) */
const DGV_MAX_SPAN = 30e6;
const C = { text: '#1f2937', muted: '#5b6b7e', faint: '#8796a8', grid: '#e3e8ee', line: '#d3dbe5', lost: '#dc2626', accent: '#0f766e', fusion5: '#7c3aed', fusion3: '#a855f7' };
const CONF = { 3: { label: 'Green', color: '#16a34a', bg: '#dcfce7' }, 2: { label: 'Amber', color: '#d97706', bg: '#fef3c7' }, 1: { label: 'Red', color: '#dc2626', bg: '#fee2e2' } } as const;
/** cytogenetic band stains (UCSC gieStain) */
const STAIN: Record<string, string> = { gneg: '#f8fafc', gpos25: '#d6dde6', gpos33: '#cbd3de', gpos50: '#a9b5c4', gpos66: '#8f9db0', gpos75: '#768699', gpos100: '#4b5869', gvar: '#c7d0f0', stalk: '#e2e8f0', acen: '#f08080' };
const fmt = (n: number) => Math.round(n).toLocaleString('en-US');
const widths = new Map<string, number>();
let measureCtx: CanvasRenderingContext2D | null | undefined;
/** The width of a gene name in bold italic 11 px (measured by the browser once per name; estimated when it cannot). */
function textWidth(s: string): number {
  const hit = widths.get(s); if (hit != null) return hit;
  if (measureCtx === undefined) { try { measureCtx = document.createElement('canvas').getContext('2d'); if (measureCtx) measureCtx.font = `italic 700 11px ${FONT}`; } catch { measureCtx = null; } }
  const w = measureCtx ? measureCtx.measureText(s).width : s.length * 7;
  if (widths.size > 20000) widths.clear();
  widths.set(s, w);
  return w;
}
const kb = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} Mb` : n >= 1000 ? `${(n / 1000).toFixed(n >= 1e4 ? 0 : 1)} kb` : `${Math.round(n)} bp`);
/** a ruler label: Mb with the decimals the step needs, bases once the step is under 10 kb */
const tickText = (t: number, step: number) => (step >= 1e4 ? `${(t / 1e6).toFixed(Math.max(0, Math.min(3, Math.ceil(-Math.log10(step / 1e6) - 1e-9))))} Mb` : fmt(t));
/** keeps a coordinate within what an SVG renderer draws exactly (a far end of a long arc or ribbon) */
const cx = (x: number) => Math.max(-4000, Math.min(VW + 4000, x));
const easeInOut = (t: number) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

export interface SvGeneLoaders {
  getTranscripts?: (chrom: string, start: number, end: number) => Promise<SvTranscript[]>;
  getAnnotations?: (chrom: string, start: number, end: number, symbols: string[]) => Promise<{ dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[] }>;
  /** the cytogenetic bands of a chromosome (the ideogram, the ISCN of the event) */
  getCytobands?: (chrom: string) => Promise<Cytoband[]>;
  /** CNV regions of the general population (DGV Gold Standard) over a 0-based half-open window */
  getPopulationSvs?: (chrom: string, start: number, end: number) => Promise<PopulationSv[]>;
  /** where a gene lies (0-based half-open): to go to a gene typed by name */
  locateGene?: (symbol: string) => Promise<{ chrom: string; start: number; end: number; name: string }>;
  /** the models the viewer already has (the gene in view): used when the source cannot list a window's transcripts */
  fallback: SvTranscript[];
}

interface Props {
  /** the rearrangement, or null for the reference alone (a region typed as coordinates) */
  rec: SvReconstruction | null; chrom: string;
  /** the event's window (or the region to open on, without rearrangement) */
  window: { start: number; end: number };
  segColors: Map<number, string>; flankCopies: number; zygosity: 'het' | 'hom'; build: string;
  loaders: SvGeneLoaders;
  /** the same derivative over a wider window; by default the junctions are walked again */
  rebuild?: (w: { start: number; end: number }) => SvReconstruction;
  /** a place to go to: each new `key` flies there */
  focus?: { start: number; end: number; key: number } | null;
}

// ======================== Genes of the window ========================

/** One gene per symbol, its principal model the MANE Select (else RefSeq Select, else the longest coding sequence). */
function groupGenes(txs: SvTranscript[]): SvGene[] {
  const by = new Map<string, SvTranscript[]>();
  for (const t of txs) { const l = by.get(t.gene); if (l) l.push(t); else by.set(t.gene, [t]); }
  const cds = (t: SvTranscript) => (t.cdsStart != null && t.cdsEnd != null ? t.exons.reduce((a, e) => a + Math.max(0, Math.min(e.end, t.cdsEnd!) - Math.max(e.start, t.cdsStart!)), 0) : 0);
  const rank = (t: SvTranscript) => (t.status === 'mane' ? 0 : t.status === 'select' ? 1 : t.status === 'mane-plus' ? 2 : 3);
  return [...by.entries()].map(([symbol, list]) => {
    const sorted = list.length > 1 ? [...list].sort((a, b) => rank(a) - rank(b) || cds(b) - cds(a) || (b.end - b.start) - (a.end - a.start)) : list;
    let start = Infinity, end = -Infinity;
    for (const t of list) { if (t.start < start) start = t.start; if (t.end > end) end = t.end; }
    return { symbol, strand: sorted[0].strand, start, end, transcripts: sorted, principal: sorted[0] };
  }).sort((a, b) => a.start - b.start);
}
/** What a 16 or 64 Mb tile keeps: one model per gene, spanning the gene, without its exons (drawn as a bar that far out). */
const reduceTile = (list: SvTranscript[]): SvTranscript[] => groupGenes(list).map(g => ({ ...g.principal, start: g.start, end: g.end, exons: [] }));

/** Domains of a BED file: chrom, start, end (0-based half-open) and an optional name; `track` / `browser` / `#` lines skipped. */
export function parseDomainsBed(text: string): Map<string, Domain[]> {
  const out = new Map<string, Domain[]>();
  for (const line of text.split(/\r?\n/)) {
    if (!line.trim() || /^(track|browser|#)/.test(line)) continue;
    const f = line.split(/\t|\s+/);
    const s = Number(f[1]), e = Number(f[2]);
    if (f.length < 3 || !Number.isFinite(s) || !Number.isFinite(e) || e <= s) continue;
    const c = f[0].startsWith('chr') ? f[0] : `chr${f[0]}`;
    out.set(c, [...(out.get(c) ?? []), { start: s, end: e, name: f[3] }]);
  }
  for (const list of out.values()) list.sort((a, b) => a.start - b.start);
  return out;
}

// ======================== Persistence (per browser, best effort) ========================

const store = {
  get<T>(k: string): T | null { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) as T : null; } catch { return null; } },
  set(k: string, v: unknown) { try { if (v == null) localStorage.removeItem(k); else localStorage.setItem(k, JSON.stringify(v)); } catch { /* storage full or blocked */ } },
};
const toPanel = (p: LoadedPanel): GenePanel => ({ source: p.source, id: p.id, name: p.name, version: p.version, genes: new Map(p.genes.map(g => [g.symbol.toUpperCase(), g])) });

// ======================== Drawing helpers ========================

/**
 * Items in lanes, each with its label; when that would take more than MAX_LANES lanes (a gene-dense region seen whole),
 * glyphs alone are packed (sharing a lane past the cap), and a label is drawn only where it fits in its lane, clear of
 * the other glyphs and labels: the items worth a label first, then the others (`labeled`; the rest show on hover).
 */
function pack<T extends { x0: number; x1: number; xg: number; h: number; worth: boolean; lw: number }>(items: T[], maxLanes = MAX_LANES): { items: (T & { lane: number; labeled: boolean })[]; laneY: number[]; height: number } {
  const lanes: number[] = [], laneH0: number[] = [];
  const sorted = items.map((it, i) => ({ it, i })).sort((a, b) => a.it.x0 - b.it.x0);
  const lane0: number[] = new Array(items.length);
  for (const { it, i } of sorted) {
    let l = lanes.findIndex(end => end + 8 < it.x0);
    if (l < 0) { l = lanes.length; lanes.push(-Infinity); laneH0.push(0); }
    lanes[l] = Math.max(lanes[l], it.x1); laneH0[l] = Math.max(laneH0[l], it.h); lane0[i] = l;
  }
  let lane = lane0, laneH = laneH0;
  const labeled = new Uint8Array(items.length);
  if (laneH0.length <= maxLanes) labeled.fill(1);
  else {
    const spans: [number, number][][] = []; laneH = []; lane = new Array(items.length);
    const fits = (l: number, a: number, b: number) => spans[l].every(([x, y]) => y + 8 < a || x - 8 > b);
    const put = (a: number, b: number): number => {
      for (let l = 0; l < spans.length; l++) if (fits(l, a, b)) return l;
      if (spans.length < maxLanes) { spans.push([]); laneH.push(0); return spans.length - 1; }
      return -1;
    };
    const order = items.map((it, i) => ({ it, i })).sort((p, q) => Number(q.it.worth) - Number(p.it.worth) || p.it.x0 - q.it.x0);
    for (const { it, i } of order) {
      let l = put(it.x0, Math.max(it.xg, it.x0 + it.lw));
      if (l >= 0) { labeled[i] = 1; spans[l].push([it.x0, Math.max(it.xg, it.x0 + it.lw)]); }
      else {
        l = put(it.x0, it.xg);
        if (l < 0) l = spans.reduce((best, sp, k) => (sp.length < spans[best].length ? k : best), 0);
        spans[l].push([it.x0, it.xg]);
      }
      lane[i] = l; laneH[l] = Math.max(laneH[l], it.h);
    }
  }
  const laneY: number[] = []; let y = 0;
  for (const h of laneH) { laneY.push(y); y += h + LANE_GAP; }
  return { items: items.map((it, i) => ({ ...it, lane: lane[i], labeled: !!labeled[i] })), laneY, height: Math.max(26, y) };
}

/**
 * A transcript (or its part [lo, hi]) as three paths: its line with strand chevrons, its UTRs, its coding exons. `X`
 * maps reference positions to the drawing; only the part in the plot gets chevrons.
 */
function glyphPaths(t: SvTranscript, X: (p: number) => number, yy: number, lo = t.start, hi = t.end, dir: 1 | -1 = t.strand): { line: string; utr: string; cds: string } {
  const a = cx(X(lo)), b = cx(X(hi)), l = Math.min(a, b), r = Math.max(a, b);
  let line = `M${l.toFixed(1)},${yy}H${r.toFixed(1)}`;
  const c0 = Math.max(l + 10, X0 - ((X0 - l) % 16)), c1 = Math.min(r - 6, X1 + 8);
  for (let x = c0; x < c1; x += 16) line += `M${(x - 2.5 * dir).toFixed(1)},${yy - 3}l${4 * dir},3l${-4 * dir},3`;
  let utr = '', cds = '';
  const box = (s: number, e: number, h: number) => { const xa = cx(Math.min(X(s), X(e))), w = Math.max(1.2, Math.abs(cx(X(e)) - cx(X(s)))); return `M${xa.toFixed(1)},${yy - h}h${w.toFixed(1)}v${2 * h}h${(-w).toFixed(1)}z`; };
  for (const e of t.exons) {
    const s = Math.max(e.start, lo), en = Math.min(e.end, hi);
    if (en <= s) continue;
    if (Math.max(X(s), X(en)) < X0 - 4 || Math.min(X(s), X(en)) > X1 + 4) continue;
    utr += box(s, en, 3.5);
    if (t.cdsStart != null && t.cdsEnd != null) { const cs = Math.max(s, t.cdsStart), ce = Math.min(en, t.cdsEnd); if (ce > cs) cds += box(cs, ce, 5.5); }
  }
  // a model without exons (a far-out tile): one bar for the whole span
  if (!t.exons.length) utr = `M${l.toFixed(1)},${yy - 2.5}h${Math.max(1.5, r - l).toFixed(1)}v5h${(-Math.max(1.5, r - l)).toFixed(1)}z`;
  return { line, utr, cds };
}

/** Breakpoint arcs and segments of the rebuilt chromosome: each piece's reference stretch, the first and last reaching the chromosome's ends. */
const pieceRange = (dp: DerPiece, L: number): [number, number] => [dp.seg.first ? 0 : dp.seg.start, dp.seg.last ? L : dp.seg.end];

/**
 * Where the rebuilt chromosome should sit under reference position `c`: its own position of `c` (the nearest retained
 * edge when `c` was lost); among several copies, the one nearest `current`, so the drawing does not jump.
 */
function derTarget(d: DerLayout, c: number, current: number | null, L: number): { piece: number; pos: number } {
  let best: { piece: number; pos: number; pen: number } | null = null;
  for (const dp of d.pieces) {
    const [lo, hi] = pieceRange(dp, L);
    const cc = Math.min(hi, Math.max(lo, c)), pen = Math.abs(c - cc), pos = toDer(dp, cc);
    if (!best || pen < best.pen - 1 || (Math.abs(pen - best.pen) <= 1 && current != null && Math.abs(pos - current) < Math.abs(best.pos - current))) best = { piece: dp.index, pos, pen };
  }
  return best ? { piece: best.piece, pos: best.pos } : { piece: -1, pos: c };
}

const EMPTY_STORY: GeneStory = { ders: [], genes: [], chains: [], fragments: [], domains: [], derDomains: [] };
type Vp = { s: number; e: number };
type TipAt = { x: number; y: number; gene: string; dgv?: PopulationSv };

// ======================== Component ========================

export default function SvGenomeBrowser({ rec: eventRec, chrom, window: eventWin, segColors, flankCopies, zygosity, build, loaders, rebuild, focus }: Props) {
  // the drawing as wide as its box, in px
  const [vw, setVw] = useState(1200);
  setWidth(vw);
  const widthBox = useRef<HTMLDivElement>(null);
  useLayoutEffect(() => {
    const el = widthBox.current; if (!el) return;
    const measure = () => setVw(Math.max(760, Math.round(el.clientWidth)));
    measure();
    const ro = new ResizeObserver(measure); ro.observe(el);
    return () => ro.disconnect();
  }, []);
  const L = useMemo(() => chromLength(build as GenomeBuild, chrom) ?? Math.max(eventWin.end * 1.5, 2e6), [build, chrom, eventWin.end]);
  const clampVp = (s: number, e: number): Vp => {
    const span = Math.min(L, Math.max(MIN_SPAN, e - s));
    const c = Math.min(L - span / 2, Math.max(span / 2, (s + e) / 2));
    return { s: c - span / 2, e: c + span / 2 };
  };
  const [mode, setMode] = useState<'mane' | 'all'>('mane');
  const [txs, setTxs] = useState<{ list: SvTranscript[]; loading: boolean; error?: string; source: string }>({ list: [], loading: !!eventRec, source: '' });
  const [ann, setAnn] = useState<{ dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[]; loading: boolean }>({ dosage: new Map(), constraint: new Map(), notes: [], loading: false });
  const [panelData, setPanelData] = useState<LoadedPanel | null>(() => store.get<LoadedPanel>('sashimi.svPanel'));
  const [domainsFile, setDomainsFile] = useState<{ name: string; text: string } | null>(() => store.get('sashimi.svDomains'));
  const [tip, setTip] = useState<TipAt | null>(null);
  /** a card kept open by a click (the hover card follows the pointer and goes with it) */
  const [pinned, setPinned] = useState<TipAt | null>(null);
  const tipEl = useRef<HTMLDivElement | null>(null);
  const [colorBy, setColorBy] = useState<'loeuf' | 'impact'>(() => store.get<'loeuf' | 'impact'>('sashimi.svColorBy') ?? 'loeuf');
  const [near, setNear] = useState<number>(() => store.get<number>('sashimi.svNear') ?? NEAR_BREAKPOINT_DEFAULT);
  /** the population CNVs (DGV Gold Standard): drawn under the reference genes, and compared with the event's losses and gains */
  const [showDgv, setShowDgv] = useState<boolean>(() => store.get<boolean>('sashimi.svDgv') ?? true);
  const [pop, setPop] = useState<{ list: PopulationSv[]; loading: boolean; error?: string }>({ list: [], loading: false });
  const [bands, setBands] = useState<Cytoband[]>([]);
  const svgRef = useRef<SVGSVGElement>(null);
  const boxRef = useRef<HTMLDivElement>(null);
  const guideRef = useRef<SVGGElement>(null);

  // ======================== the view: a stretch of the reference, moved once per frame ========================
  const initialVp = (): Vp => { const pad = Math.max(2000, (eventWin.end - eventWin.start) * (eventRec ? 0.12 : 0)); return clampVp(eventWin.start - pad, eventWin.end + pad); };
  const [vp, setVpState] = useState<Vp>(initialVp);
  const vpRef = useRef(vp);
  const raf = useRef(0);
  const anim = useRef<{ fly?: { from: Vp; to: Vp; t0: number; dur: number }; fling?: { v: number; t: number } }>({});
  /** the rebuilt chromosomes' offset from where they should sit, decaying to 0 (a slide instead of a jump) */
  const derState = useRef(new Map<string, { piece: number; target: number; offset: number }>());
  const lastMove = useRef(0);
  const [, settle] = useReducer((x: number) => x + 1, 0);
  const settleTimer = useRef(0);
  const frame = () => {
    raf.current = 0;
    const now = performance.now();
    let next = vpRef.current, more = false;
    const a = anim.current;
    if (a.fly) {
      const t = Math.min(1, (now - a.fly.t0) / a.fly.dur), k2 = easeInOut(t);
      const c0 = (a.fly.from.s + a.fly.from.e) / 2, c1 = (a.fly.to.s + a.fly.to.e) / 2;
      const l0 = Math.log(a.fly.from.e - a.fly.from.s), l1 = Math.log(a.fly.to.e - a.fly.to.s);
      const span = Math.exp(l0 + (l1 - l0) * k2), c = c0 + (c1 - c0) * k2;
      next = clampVp(c - span / 2, c + span / 2);
      if (t < 1) more = true; else a.fly = undefined;
    } else if (a.fling) {
      const dt = Math.min(48, now - a.fling.t); a.fling.t = now;
      const sh = a.fling.v * dt;
      next = clampVp(next.s + sh, next.e + sh);
      a.fling.v *= Math.pow(0.9, dt / 16);
      if (Math.abs(a.fling.v) * 16 > (next.e - next.s) / PW * 0.25 && next.s > 0 && next.e < L) more = true; else a.fling = undefined;
    }
    const bpPerPx = (next.e - next.s) / PW;
    for (const st of derState.current.values()) {
      if (Math.abs(st.offset) > bpPerPx * 0.5) { st.offset *= 0.8; more = true; } else st.offset = 0;
    }
    vpRef.current = next;
    lastMove.current = now;
    setVpState(next);
    if (more) raf.current = requestAnimationFrame(frame);
    window.clearTimeout(settleTimer.current);
    settleTimer.current = window.setTimeout(settle, 450);
  };
  const kick = () => { if (!raf.current) raf.current = requestAnimationFrame(frame); };
  /** a move by the user: animations stop, the view follows */
  const moveTo = (s: number, e: number) => { anim.current = {}; vpRef.current = clampVp(s, e); kick(); };
  const flyTo = (s: number, e: number, dur = 480) => { anim.current = { fly: { from: vpRef.current, to: clampVp(s, e), t0: performance.now(), dur } }; kick(); };
  useEffect(() => () => { cancelAnimationFrame(raf.current); window.clearTimeout(settleTimer.current); }, []);
  // a new chromosome or event: start over on it
  useEffect(() => { derState.current.clear(); const v = initialVp(); anim.current = {}; vpRef.current = v; setVpState(v); }, [chrom, build, eventWin.start, eventWin.end]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (focus) flyTo(focus.start, focus.end, 700); }, [focus?.key]); // eslint-disable-line react-hooks/exhaustive-deps
  const span = vp.e - vp.s, ppb = PW / span;
  const refX = (pos: number) => X0 + (pos - vp.s) * ppb;
  const lod: 'full' | 'bars' | 'density' = span <= LOD_FULL ? 'full' : span <= LOD_BARS ? 'bars' : 'density';
  const interacting = performance.now() - lastMove.current < 420;

  // ======================== the event: its genes, worked out at full detail ========================
  const bpKey = eventRec ? eventRec.junctions.flatMap(j => [j.a, j.b]).filter(e => e.chrom === chrom).map(e => e.pos).sort((a, b) => a - b).join(',') : '';
  const win = useMemo(() => {
    const bps = bpKey ? bpKey.split(',').map(Number) : [];
    if (!eventRec || !bps.length || near <= 0) return eventWin;
    return { start: Math.max(0, Math.min(eventWin.start, bps[0] - near)), end: Math.min(Math.max(L, eventWin.end), Math.max(eventWin.end, bps[bps.length - 1] + near)) };
  }, [eventWin.start, eventWin.end, bpKey, near, L, !!eventRec]); // eslint-disable-line react-hooks/exhaustive-deps
  const rec = useMemo(() => (!eventRec ? null : win === eventWin ? eventRec : rebuild ? rebuild(win) : reconstruct([{ chrom, start: win.start, end: win.end }], eventRec.junctions, eventRec.choice)), [eventRec, win]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    let alive = true;
    if (!eventRec) { setTxs({ list: [], loading: false, source: '' }); return; }
    if (!loaders.getTranscripts) { setTxs({ list: loaders.fallback, loading: false, source: 'the gene in view only (this source cannot list the window\'s genes)' }); return; }
    setTxs(t => ({ ...t, loading: true, error: undefined }));
    loaders.getTranscripts(chrom, win.start, win.end)
      .then(list => { if (alive) setTxs({ list, loading: false, source: 'RefSeq NM_ models (UCSC ncbiRefSeqCurated), MANE status from the UCSC mane track' }); })
      .catch(e => { if (alive) setTxs({ list: loaders.fallback, loading: false, error: `gene models of the window unavailable (${e?.message || e}): the gene in view only`, source: '' }); });
    return () => { alive = false; };
  }, [chrom, win.start, win.end, !!eventRec]); // eslint-disable-line react-hooks/exhaustive-deps
  const genes = useMemo(() => groupGenes(txs.list), [txs.list]);
  const symbols = useMemo(() => genes.map(g => g.symbol), [genes]);
  /** annotations merge in as they come (the event's genes, then those of wherever the view goes) */
  const mergeAnn = (r: { dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[] }) =>
    setAnn(a => ({ dosage: new Map([...a.dosage, ...r.dosage]), constraint: new Map([...a.constraint, ...r.constraint]), notes: a.notes.length ? a.notes : r.notes, loading: false }));
  const annTried = useRef(new Set<string>());
  useEffect(() => { annTried.current = new Set(); setAnn({ dosage: new Map(), constraint: new Map(), notes: [], loading: false }); }, [chrom, build]);
  useEffect(() => {
    if (!loaders.getAnnotations || !symbols.length) return;
    let alive = true;
    symbols.forEach(s => annTried.current.add(s));
    setAnn(a => ({ ...a, loading: true }));
    loaders.getAnnotations(chrom, win.start, win.end, symbols).then(r => { if (alive) mergeAnn(r); }).catch(e => { if (alive) setAnn(a => ({ ...a, notes: [`annotations unavailable (${e?.message || e})`], loading: false })); });
    return () => { alive = false; };
  }, [symbols.join(','), chrom]); // eslint-disable-line react-hooks/exhaustive-deps
  // the event's population CNVs, to compare with its losses and gains
  useEffect(() => {
    if (!eventRec || !showDgv || !loaders.getPopulationSvs) { setPop({ list: [], loading: false }); return; }
    let alive = true;
    setPop(p => ({ ...p, loading: true, error: undefined }));
    loaders.getPopulationSvs(chrom, win.start, win.end)
      .then(list => { if (alive) setPop({ list, loading: false }); })
      .catch(e => { if (alive) setPop({ list: [], loading: false, error: String(e?.message || e) }); });
    return () => { alive = false; };
  }, [showDgv, chrom, win.start, win.end, !!eventRec]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    let alive = true; setBands([]);
    loaders.getCytobands?.(chrom).then(b => { if (alive) setBands(b); }).catch(() => { /* a plain bar */ });
    return () => { alive = false; };
  }, [chrom, loaders]);

  const domains = useMemo(() => (domainsFile ? parseDomainsBed(domainsFile.text).get(chrom.startsWith('chr') ? chrom : `chr${chrom}`) ?? [] : []), [domainsFile, chrom]);
  const story = useMemo(() => (rec ? geneStory(rec, genes, win, { flankCopies, zygosity, domains, near }) : EMPTY_STORY), [rec, genes, win, flankCopies, zygosity, domains, near]);
  const panel = useMemo(() => (panelData ? toPanel(panelData) : null), [panelData]);
  const cnv = useMemo(() => (rec ? cnvEvents(rec, chrom).map(ev => ({ ev, matches: matchesOf(ev, pop.list) })) : []), [rec, chrom, pop.list]);
  /** population CNVs that are the same event as a loss or gain (same type, reciprocal overlap ≥ RO_MATCH) */
  const popSame = useMemo(() => new Set(cnv.flatMap(c => c.matches.filter(m => m.sameType && m.ro >= RO_MATCH).map(m => `${m.sv.id}:${m.sv.outer[0]}`))), [cnv]);
  const annotations: GeneAnnotations = useMemo(() => ({ dosage: ann.dosage, constraint: ann.constraint, panel }), [ann, panel]);
  const impactOf = useMemo(() => new Map(story.genes.map(g => [g.gene.symbol, g])), [story]);
  const geneColor = (symbol: string, gi: GeneImpact | undefined) => (colorBy === 'loeuf' ? LOEUF_CLASS[loeufClass(ann.constraint.get(symbol))].color : IMPACT[gi?.kind ?? 'unaffected'].color);
  const panelDot = (sym: string) => { const pg = panel?.genes.get(sym.toUpperCase()); return pg ? CONF[pg.confidence].color : null; };
  /** worth a label in a crowded view: changed by the event, or in the panel, dosage-sensitive or LoF-intolerant */
  const worthSym = (sym: string) => { const gi = impactOf.get(sym); return (!!gi && gi.kind !== 'unaffected' && gi.kind !== 'nearby') || notableSym(sym, annotations); };

  // ======================== the rebuilt chromosomes: under the reference, sliding where the two part ways ========================
  const ders = story.ders.filter(d => d.rearranged);
  const c = (vp.s + vp.e) / 2;
  const derCenter = new Map<DerLayout, number>();
  for (const d of ders) {
    const st = derState.current.get(d.name);
    const t = derTarget(d, c, st ? st.target + st.offset : null, L);
    if (!st) derState.current.set(d.name, { piece: t.piece, target: t.pos, offset: 0 });
    else {
      if (st.piece !== t.piece) { st.offset = st.target + st.offset - t.pos; st.piece = t.piece; kick(); }
      st.target = t.pos;
    }
    const s2 = derState.current.get(d.name)!;
    derCenter.set(d, s2.target + s2.offset);
  }
  const derX = (d: DerLayout, dpos: number) => XC + (dpos - (derCenter.get(d) ?? 0)) * ppb;
  /** the stretches of the reference the rebuilt chromosomes show (their genes are fetched and drawn too) */
  let needS = vp.s, needE = vp.e;
  for (const d of ders) {
    const dc = derCenter.get(d) ?? 0, dl = dc + (X0 - XC) / ppb, dr = dc + (X1 - XC) / ppb;
    for (const dp of d.pieces) {
      const [lo, hi] = pieceRange(dp, L);
      const a = Math.min(toDer(dp, lo), toDer(dp, hi)), b = Math.max(toDer(dp, lo), toDer(dp, hi));
      const oa = Math.max(a, dl), ob = Math.min(b, dr);
      if (ob <= oa) continue;
      const back = (dpos: number) => (dp.piece.reversed ? dp.seg.end - (dpos - dp.derStart) : dp.seg.start + (dpos - dp.derStart));
      needS = Math.min(needS, back(oa), back(ob)); needE = Math.max(needE, back(oa), back(ob));
    }
  }
  needS = Math.max(0, needS); needE = Math.min(L, needE);

  // ======================== wherever the view goes: genes and population CNVs by tiles ========================
  const [tilesVersion, bumpTiles] = useReducer((x: number) => x + 1, 0);
  const geneTiles = useMemo(() => (loaders.getTranscripts ? new TileCache<SvTranscript>((a, b, size) => loaders.getTranscripts!(chrom, a, b).then(list => (size >= 16e6 ? reduceTile(list) : list)), L, 36, bumpTiles) : null), [loaders, chrom, L]);
  const dgvTiles = useMemo(() => (loaders.getPopulationSvs ? new TileCache<PopulationSv>((a, b) => loaders.getPopulationSvs!(chrom, a, b), L, 24, bumpTiles) : null), [loaders, chrom, L]);
  const tileSize = tileSizeFor(span);
  const firstFetch = useRef(true);
  // the tiles are asked for once the view rests a moment: a fast zoom does not fetch every level it crosses
  useEffect(() => {
    const go = () => { geneTiles?.want(needS, needE, tileSize); if (showDgv && span <= DGV_MAX_SPAN) dgvTiles?.want(vp.s, vp.e, tileSize); };
    if (firstFetch.current) { firstFetch.current = false; go(); return; }
    const id = window.setTimeout(go, 140);
    return () => window.clearTimeout(id);
  }, [vp.s, vp.e, needS, needE, tileSize, geneTiles, dgvTiles, showDgv]); // eslint-disable-line react-hooks/exhaustive-deps
  const i0 = Math.floor(needS / tileSize), i1 = Math.floor(needE / tileSize);
  const view = useMemo(() => {
    const a = (i0 - 1) * tileSize, b = (i1 + 2) * tileSize;
    const got = geneTiles ? geneTiles.collect(a, b, tileSize) : { items: loaders.fallback, loading: false, error: null };
    const byId = new Map<string, SvTranscript>();
    for (const t of got.items) { const o = byId.get(t.id); if (!o || (!o.exons.length && t.exons.length)) byId.set(t.id, t); }
    for (const t of txs.list) byId.set(t.id, t);
    return { genes: groupGenes([...byId.values()]), loading: got.loading, error: got.error };
  }, [tilesVersion, tileSize, i0, i1, txs.list, geneTiles]); // eslint-disable-line react-hooks/exhaustive-deps
  const dgvView = useMemo(() => {
    if (!dgvTiles || !showDgv || span > DGV_MAX_SPAN) return { list: [] as PopulationSv[], loading: false, error: null as string | null };
    const got = dgvTiles.collect((Math.floor(vp.s / tileSize) - 1) * tileSize, (Math.floor(vp.e / tileSize) + 2) * tileSize, tileSize);
    const seen = new Set<string>(), list: PopulationSv[] = [];
    for (const v of got.items) { const k = `${v.id}:${v.outer[0]}`; if (!seen.has(k)) { seen.add(k); list.push(v); } }
    return { list, loading: got.loading, error: got.error };
  }, [tilesVersion, tileSize, Math.floor(vp.s / tileSize), Math.floor(vp.e / tileSize), dgvTiles, showDgv, span > DGV_MAX_SPAN]); // eslint-disable-line react-hooks/exhaustive-deps
  const geneBySym = useMemo(() => new Map(view.genes.map(g => [g.symbol, g])), [view.genes]);
  const visGenes = view.genes.filter(g => g.end >= vp.s - span * 0.02 && g.start <= vp.e + span * 0.02);
  // constraint and dosage of the genes in view, asked once each, when the view rests
  const visKey = span <= LOD_BARS ? visGenes.map(g => g.symbol).join(',') : '';
  useEffect(() => {
    if (!loaders.getAnnotations || !visKey) return;
    const todo = visKey.split(',').filter(s => !annTried.current.has(s));
    if (!todo.length) return;
    const id = window.setTimeout(() => {
      todo.forEach(s => annTried.current.add(s));
      loaders.getAnnotations!(chrom, Math.floor(vpRef.current.s), Math.ceil(vpRef.current.e), todo).then(mergeAnn).catch(() => { /* colours stay unknown */ });
    }, 500);
    return () => window.clearTimeout(id);
  }, [visKey]); // eslint-disable-line react-hooks/exhaustive-deps

  const zoomTo = (lo: number, hi: number) => { const pad = (hi - lo) * 0.15 + 200; flyTo(lo - pad, hi + pad); };
  const eventView = () => { const pad = Math.max(2000, (eventWin.end - eventWin.start) * 0.12); flyTo(eventWin.start - pad, eventWin.end + pad); };

  // ======================== layout ========================
  const shown = (g: SvGene) => (mode === 'mane' || lod !== 'full' ? [g.principal] : g.transcripts);
  const labelW = (s: string) => textWidth(s) + 10;
  const inView = (a: number, b: number) => Math.max(a, b) >= X0 - 2 && Math.min(a, b) <= X1 + 2;
  /** sticky heights: a track grows at once but shrinks only once the view rests (no jumping while panning) */
  const sticky = useRef(new Map<string, number>());
  const keep = (key: string, h: number) => { const prev = sticky.current.get(key) ?? 0; const v = interacting ? Math.max(h, prev) : h; sticky.current.set(key, v); return v; };

  // reference genes
  const refItems = (() => {
    if (lod === 'density') {
      const worth = visGenes.filter(g => worthSym(g.symbol));
      return pack(worth.map(g => { const a = Math.max(X0 + 2, refX(g.start)); const lw = labelW(g.symbol); return { g, x0: a, x1: a + lw, xg: a + 3, h: LABEL_H + 2, worth: true, lw }; }), 3);
    }
    return pack(visGenes.filter(g => inView(refX(g.start), refX(g.end))).map(g => {
      const a = Math.max(X0 + 2, refX(g.start)), b = refX(g.end), lw = labelW(g.symbol) + (panelDot(g.symbol) ? 10 : 0);
      return { g, x0: a, x1: Math.max(b, a + lw), xg: Math.max(b, a + 3), h: LABEL_H + (lod === 'full' ? shown(g).length * TX_H : 8), worth: worthSym(g.symbol), lw };
    }));
  })();
  /** a density profile of the genes, for a chromosome seen whole */
  const density = useMemo(() => {
    if (lod !== 'density') return null;
    const BIN = 3, n = Math.ceil(PW / BIN), counts = new Float32Array(n);
    for (const g of view.genes) {
      const a = Math.floor((refX(g.start) - X0) / BIN), b = Math.floor((refX(g.end) - X0) / BIN);
      for (let i = Math.max(0, a); i <= Math.min(n - 1, b); i++) counts[i]++;
    }
    let max = 1; for (const v of counts) max = Math.max(max, v);
    let d = `M${X0},0`;
    for (let i = 0; i < n; i++) d += `L${X0 + i * BIN},${(-counts[i] / max).toFixed(3)}L${X0 + (i + 1) * BIN},${(-counts[i] / max).toFixed(3)}`;
    return { d: `${d}L${X1},0Z`, max };
  }, [lod, view.genes, vp.s, vp.e]); // eslint-disable-line react-hooks/exhaustive-deps

  // rebuilt genes: one block per chain (a fusion's two partners in one block); beyond the event's window, the genes of the view carried along
  interface DerBlock { der: DerLayout; chains: Chain[]; genes: string[]; label: string; x0: number; x1: number; xg: number; h: number; fusion: GeneImpact['principal']['fusions'][number] | null; worth: boolean; lw: number; flank?: { g: SvGene; dp: DerPiece } }
  const derPacked = ders.map(d => {
    const list: DerBlock[] = [];
    if (lod !== 'density') {
      const keepC = (ch: Chain) => ch.parts.every(p => (mode === 'mane' || lod !== 'full' ? p.tx === impactOf.get(p.tx.gene)?.gene.principal : true));
      const copyNo = new Map<string, number>();
      for (const ch of story.chains) {
        if (story.ders[ch.der] !== d || !keepC(ch)) continue;
        const geneList = [...new Set(ch.parts.map(p => p.tx.gene))];
        const xs = ch.parts.flatMap(p => [derX(d, toDer(p.dp, p.lo)), derX(d, toDer(p.dp, p.hi))]);
        if (!inView(Math.min(...xs), Math.max(...xs))) continue;
        const x0 = Math.max(X0 + 2, Math.min(...xs)), x1 = Math.max(...xs);
        const first = ch.parts[0], last = ch.parts[ch.parts.length - 1];
        const gi = impactOf.get(geneList[0]);
        const fusion = geneList.length > 1 ? gi?.principal.fusions.find(f => f.five.gene === geneList[0] && f.three.gene === geneList[1]) ?? impactOf.get(geneList[1])?.principal.fusions.find(f => f.five.gene === geneList[0]) ?? null : null;
        const five = first.tx.strand > 0 ? !first.cutLo || first.edgeLo : !first.cutHi || first.edgeHi;
        const three = last.tx.strand > 0 ? !last.cutHi || last.edgeHi : !last.cutLo || last.edgeLo;
        let label = geneList.join('::');
        if (geneList.length === 1) {
          if (!five && !three) label += ' · inner part'; else if (!five) label += ' · 3′ part'; else if (!three) label += ' · 5′ part';
          else { const key = `${ch.der}|${geneList[0]}|${first.tx.id}`; const n = (copyNo.get(key) ?? 0) + 1; copyNo.set(key, n); if ((gi?.principal.intact ?? 0) > 1 && n > 1) label += ` · copy ${n}`; }
        }
        const same = list.find(b => b.genes.join() === geneList.join() && Math.abs(b.x0 - x0) < 2 && Math.abs(b.x1 - x1) < 2);
        if (same) { same.chains.push(ch); same.h += TX_H; continue; }
        const lw = labelW(label) + (fusion ? 80 : 0) + (geneList.length === 1 && panelDot(geneList[0]) ? 10 : 0);
        list.push({ der: d, chains: [ch], genes: geneList, label, x0, x1: Math.max(x1, x0 + lw), xg: Math.max(x1, x0 + 3), h: LABEL_H + (lod === 'full' ? TX_H : 8), fusion, worth: geneList.length > 1 || worthSym(geneList[0]), lw });
      }
      // genes of the view in the flanks (the chromosome on either side of the event's window), carried as they are
      for (const dp of d.pieces) {
        if (!dp.seg.first && !dp.seg.last) continue;
        const [lo, hi] = pieceRange(dp, L);
        for (const g of view.genes) {
          if (impactOf.has(g.symbol) || g.start < lo || g.end > hi || (g.end > win.start && g.start < win.end)) continue;
          const xa = derX(d, toDer(dp, g.start)), xb = derX(d, toDer(dp, g.end));
          if (!inView(xa, xb)) continue;
          const x0 = Math.max(X0 + 2, Math.min(xa, xb)), x1 = Math.max(xa, xb), lw = labelW(g.symbol) + (panelDot(g.symbol) ? 10 : 0);
          list.push({ der: d, chains: [], genes: [g.symbol], label: g.symbol, x0, x1: Math.max(x1, x0 + lw), xg: Math.max(x1, x0 + 3), h: LABEL_H + (lod === 'full' ? TX_H : 8), fusion: null, worth: worthSym(g.symbol), lw, flank: { g, dp } });
        }
      }
    }
    return pack(list);
  });

  // ---- vertical layout ----
  let y = IDEO_Y + IDEO_H + LENS_H;
  const yRuler = y; y += RULER_H;
  const yRefTad = domains.length ? y : -1; if (domains.length) y += TAD_H + 6;
  const yRefGenes = y; y += keep('ref', lod === 'density' ? refItems.height + 40 : refItems.height) + 4;
  const yDens = yRefGenes + refItems.height + 2;
  const dgvShown = showDgv && !!loaders.getPopulationSvs;
  const dgvItems: { v: PopulationSv; lane: number; a: number; b: number; ia: number; ib: number }[] = [];
  if (dgvShown) {
    const ends: number[] = [];
    for (const v of dgvView.list) {
      const a = refX(v.outer[0]), b = Math.max(refX(v.outer[1]), a + 2);
      if (b < X0 || a > X1) continue;
      let lane = ends.findIndex(e => e + 3 < a);
      if (lane < 0) { if (ends.length < DGV_MAX_LANES) { lane = ends.length; ends.push(b); } else { lane = ends.indexOf(Math.min(...ends)); ends[lane] = Math.max(ends[lane], b); } }
      else ends[lane] = b;
      dgvItems.push({ v, lane, a: cx(a), b: cx(b), ia: cx(refX(v.inner[0])), ib: cx(refX(v.inner[1])) });
    }
  }
  const yDgv = dgvShown ? y : -1; if (dgvShown) y += keep('dgv', Math.max(1, ...dgvItems.map(d => d.lane + 1)) * DGV_LANE_H) + 10;
  const yRefBand = rec ? y : -1; if (rec) y += BAND_H;
  const yRib = y; if (ders.length) y += RIB_H;
  const derY = ders.map((d, i) => {
    const band = y; y += BAND_H;
    const ruler = y; y += 24;
    const genesY = y; y += keep(`der${i}`, derPacked[i].height) + 4;
    const tad = domains.length ? y : -1; if (domains.length) y += TAD_H;
    y += 18;
    return { band, ruler, genes: genesY, tad };
  });
  const H = y + 6;

  // ---- cards: hover follows the pointer (moved directly, no redraw), a click keeps one open ----
  const boxPos = (e: React.MouseEvent) => { const r = boxRef.current!.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const placeTip = (x: number, yy: number) => {
    const el = tipEl.current, box = boxRef.current; if (!el || !box) return;
    const w = box.clientWidth, h = box.clientHeight;
    el.style.left = `${Math.max(0, Math.min(x + 14, w - 400))}px`;
    if (yy > h - 160) { el.style.top = ''; el.style.bottom = `${h - yy + 14}px`; } else { el.style.bottom = ''; el.style.top = `${yy + 14}px`; }
  };
  const lastTipPos = useRef({ x: 0, y: 0 });
  const hoverOn = (key: string, dgv?: PopulationSv) => (e: React.MouseEvent) => {
    if (dragRef.current?.moved) return;
    const p = boxPos(e); lastTipPos.current = p;
    if (tip?.gene !== key) setTip({ ...p, gene: key, dgv }); else placeTip(p.x, p.y);
  };
  const hoverOff = () => setTip(null);
  useLayoutEffect(() => { if (tip) placeTip(lastTipPos.current.x, lastTipPos.current.y); });
  const pinOn = (key: string, dgv?: PopulationSv) => (e: React.MouseEvent) => {
    e.stopPropagation();
    if (suppressClick.current) return;
    const p = boxPos(e);
    setPinned(cur => (cur?.gene === key ? null : { ...p, gene: key, dgv }));
  };
  // the card points at its gene where it was: it closes when the view moves
  useEffect(() => { if (pinned && !anim.current.fly) setPinned(null); }, [Math.round(vp.s / Math.max(1, span) * 50), Math.round(Math.log(span) * 20)]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- gestures: drag to pan (a fling goes on), Ctrl / ⌘ + wheel or pinch to zoom, two-finger swipe to pan, double-click to zoom in,
  //      drag along the ruler to zoom on a stretch, click or drag on the ideogram to go there ----
  const dragRef = useRef<{ kind: 'pan' | 'band' | 'ideo'; x0: number; vp0: Vp; moved: boolean; samples: { x: number; t: number }[]; id: number } | null>(null);
  /** the click that ends a drag opens no card */
  const suppressClick = useRef(false);
  const [band, setBand] = useState<{ a: number; b: number } | null>(null);
  const toSvgX = (clientX: number) => { const r = svgRef.current!.getBoundingClientRect(); return ((clientX - r.left) / r.width) * VW; };
  const toSvgY = (clientY: number) => { const r = svgRef.current!.getBoundingClientRect(); return ((clientY - r.top) / r.height) * H; };
  const ideoPos = (x: number) => Math.max(0, Math.min(L, ((x - X0) / PW) * L));
  const onPointerDown = (e: React.PointerEvent) => {
    if (e.button !== 0) return;
    anim.current = {}; suppressClick.current = false;
    const x = toSvgX(e.clientX), yy = toSvgY(e.clientY);
    const kind = yy < IDEO_Y + IDEO_H + 4 ? 'ideo' : yy >= yRuler && yy < yRuler + RULER_H ? 'band' : 'pan';
    dragRef.current = { kind, x0: x, vp0: vpRef.current, moved: false, samples: [{ x: e.clientX, t: performance.now() }], id: e.pointerId };
    if (kind === 'ideo') { const sp = vpRef.current.e - vpRef.current.s, p = ideoPos(x); flyTo(p - sp / 2, p + sp / 2, 320); }
  };
  /** the guide line under the pointer and its position, moved without redrawing (and kept right as the view moves) */
  const guideX = useRef<number | null>(null);
  const updateGuide = () => {
    const gEl = guideRef.current, x = guideX.current; if (!gEl) return;
    if (x == null || x < X0 || x > X1) { gEl.style.display = 'none'; return; }
    gEl.setAttribute('transform', `translate(${x.toFixed(1)},0)`); gEl.style.display = '';
    const t = gEl.querySelector('text'); if (t) t.textContent = `${chrom}:${fmt(vpRef.current.s + (x - X0) / (PW / (vpRef.current.e - vpRef.current.s)) + 1)}`;
  };
  useLayoutEffect(updateGuide);
  const onPointerMove = (e: React.PointerEvent) => {
    const x = toSvgX(e.clientX);
    guideX.current = x; updateGuide();
    const d = dragRef.current; if (!d || d.id !== e.pointerId) return;
    if (!d.moved && Math.abs(x - d.x0) < 4) return;
    if (!d.moved) { d.moved = true; svgRef.current?.setPointerCapture(e.pointerId); setTip(null); }
    d.samples.push({ x: e.clientX, t: performance.now() }); if (d.samples.length > 6) d.samples.shift();
    if (d.kind === 'pan') { const sp = d.vp0.e - d.vp0.s, sh = -(x - d.x0) / PW * sp; moveTo(d.vp0.s + sh, d.vp0.e + sh); }
    else if (d.kind === 'band') setBand({ a: Math.max(X0, Math.min(X1, d.x0)), b: Math.max(X0, Math.min(X1, x)) });
    else { const sp = vpRef.current.e - vpRef.current.s, p = ideoPos(x); moveTo(p - sp / 2, p + sp / 2); }
  };
  const onPointerUp = (e: React.PointerEvent) => {
    const d = dragRef.current; if (!d || d.id !== e.pointerId) return;
    dragRef.current = null; suppressClick.current = d.moved;
    if (d.kind === 'band' && band) {
      const a = Math.min(band.a, band.b), b = Math.max(band.a, band.b); setBand(null);
      if (b - a > 4) { const sp = vpRef.current.e - vpRef.current.s, s0 = vpRef.current.s; flyTo(s0 + (a - X0) / PW * sp, s0 + (b - X0) / PW * sp); }
      return;
    }
    if (d.kind === 'pan' && d.moved && d.samples.length > 1) {
      const p = d.samples[0], q = d.samples[d.samples.length - 1], dt = q.t - p.t;
      if (dt > 0 && performance.now() - q.t < 60) {
        const r = svgRef.current!.getBoundingClientRect(), pxPerMs = (q.x - p.x) / dt * (VW / r.width);
        const sp = vpRef.current.e - vpRef.current.s;
        if (Math.abs(pxPerMs) > 0.25) { anim.current = { fling: { v: -pxPerMs / PW * sp, t: performance.now() } }; kick(); }
      }
    }
  };
  useEffect(() => {
    const el = svgRef.current; if (!el) return;
    const onWheel = (e: WheelEvent) => {
      const unit = e.deltaMode === 1 ? 33 : e.deltaMode === 2 ? 400 : 1;
      const r = el.getBoundingClientRect(), x = ((e.clientX - r.left) / r.width) * VW;
      const v = vpRef.current, sp = v.e - v.s;
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        const f = Math.exp(e.deltaY * unit * 0.0022), at = v.s + (Math.max(X0, Math.min(X1, x)) - X0) / PW * sp;
        moveTo(at - (at - v.s) * f, at + (v.e - at) * f);
      } else if (Math.abs(e.deltaX) > Math.abs(e.deltaY) || e.shiftKey) {
        e.preventDefault();
        const dx = (e.shiftKey && !e.deltaX ? e.deltaY : e.deltaX) * unit, sh = dx / r.width * VW / PW * sp;
        moveTo(v.s + sh, v.e + sh);
      }
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, [txs.loading, L]); // eslint-disable-line react-hooks/exhaustive-deps
  const onDoubleClick = (e: React.MouseEvent) => {
    const x = toSvgX(e.clientX); if (x < X0 || x > X1) return;
    const v = vpRef.current, sp = v.e - v.s, at = v.s + (x - X0) / PW * sp, f = e.shiftKey ? 2 : 0.5;
    flyTo(at - (at - v.s) * f, at + (v.e - at) * f, 320);
  };
  const onKeyDown = (e: React.KeyboardEvent) => {
    const v = vpRef.current, sp = v.e - v.s;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') { e.preventDefault(); const sh = (e.key === 'ArrowLeft' ? -1 : 1) * sp * 0.3; flyTo(v.s + sh, v.e + sh, 260); }
    else if (e.key === '+' || e.key === '=') { e.preventDefault(); flyTo(v.s + sp * 0.25, v.e - sp * 0.25, 260); }
    else if (e.key === '-' || e.key === '_') { e.preventDefault(); flyTo(v.s - sp * 0.5, v.e + sp * 0.5, 260); }
    else if ((e.key === 'e' || e.key === 'E') && rec) { e.preventDefault(); eventView(); }
  };

  // ---- ruler ----
  const step = niceStep(span, 9);
  const ticks: number[] = [];
  for (let t = Math.ceil(vp.s / step) * step; t <= vp.e; t += step) ticks.push(t);
  const ix = (pos: number) => X0 + (pos / L) * PW;

  const torn = (x: number, yy: number, color: string, key: string) => <path key={key} d={`M${x},${yy - 8} l-3,3 l3,3 l-3,3 l3,3 l-3,3`} fill="none" stroke={color} strokeWidth={1.4} />;
  const glyph = (t: SvTranscript, X: (p: number) => number, yy: number, color: string, key: string, lo?: number, hi?: number, dir?: 1 | -1) => {
    const p = glyphPaths(t, X, yy, lo, hi, dir);
    return (
      <g key={key}>
        <path d={p.line} fill="none" stroke={color} strokeWidth={1.1} strokeOpacity={0.85} />
        {p.utr && <path d={p.utr} fill={color} fillOpacity={0.5} />}
        {p.cds && <path d={p.cds} fill={color} />}
      </g>
    );
  };

  // ---- the SVG ----
  const vl = Math.max(2, ix(vp.e) - ix(vp.s));
  const svg = (
    <svg ref={svgRef} viewBox={`0 0 ${VW} ${H}`} width="100%" style={{ aspectRatio: `${VW} / ${H}`, cursor: dragRef.current?.moved && dragRef.current.kind === 'pan' ? 'grabbing' : 'grab', fontFamily: FONT, touchAction: 'pan-y' }} className="block select-none"
      onPointerDown={onPointerDown} onPointerMove={onPointerMove} onPointerUp={onPointerUp} onPointerCancel={onPointerUp}
      onPointerLeave={() => { guideX.current = null; updateGuide(); }}
      onClick={() => { if (!suppressClick.current) setPinned(null); }} onDoubleClick={onDoubleClick}>
      <defs>
        <clipPath id="svb-clip"><rect x={X0 - 2} y={0} width={PW + 4} height={H} /></clipPath>
        <clipPath id="svb-ideo"><rect x={X0} y={IDEO_Y} width={PW} height={IDEO_H} rx={IDEO_H / 2} /></clipPath>
        <pattern id="svb-lost" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#fef2f2" /><line x1="0" y1="0" x2="0" y2="6" stroke={C.lost} strokeWidth="1.6" /></pattern>
        <linearGradient id="svb-bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#fbfcfd" /><stop offset="1" stopColor="#f4f7fa" /></linearGradient>
        <linearGradient id="svb-lens" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor={C.accent} stopOpacity={0.16} /><stop offset="1" stopColor={C.accent} stopOpacity={0.03} /></linearGradient>
        <linearGradient id="svb-dens" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#64748b" stopOpacity={0.55} /><stop offset="1" stopColor="#64748b" stopOpacity={0.12} /></linearGradient>
      </defs>
      <rect width={VW} height={H} fill="url(#svb-bg)" rx={12} />
      {/* ideogram: the whole chromosome, the part in view boxed, the event's breakpoints marked */}
      <text x={14} y={IDEO_Y + 10} fontSize={11} fontWeight={700} fill={C.text}>{chrom}</text>
      <g clipPath="url(#svb-ideo)">
        <rect x={X0} y={IDEO_Y} width={PW} height={IDEO_H} fill={STAIN.gneg} />
        {bands.map(b => b.stain === 'acen'
          ? <path key={b.name + b.start} d={b.name.startsWith('p') ? `M${ix(b.start)},${IDEO_Y} L${ix(b.end)},${IDEO_Y + IDEO_H / 2} L${ix(b.start)},${IDEO_Y + IDEO_H} Z` : `M${ix(b.end)},${IDEO_Y} L${ix(b.start)},${IDEO_Y + IDEO_H / 2} L${ix(b.end)},${IDEO_Y + IDEO_H} Z`} fill={STAIN.acen} />
          : <rect key={b.name + b.start} x={ix(b.start)} y={IDEO_Y} width={Math.max(0.5, ix(b.end) - ix(b.start))} height={IDEO_H} fill={STAIN[b.stain] ?? STAIN.gneg}><title>{`${chrom.replace(/^chr/, '')}${b.name} · ${fmt(b.start + 1)}-${fmt(b.end)}`}</title></rect>)}
      </g>
      <rect x={X0} y={IDEO_Y} width={PW} height={IDEO_H} rx={IDEO_H / 2} fill="none" stroke="#94a3b8" strokeWidth={0.8} />
      {bands.filter(b => ix(b.end) - ix(b.start) > textWidth(b.name) * 0.78 + 6 && b.stain !== 'acen').map(b => <text key={`bl${b.name}${b.start}`} x={(ix(b.start) + ix(b.end)) / 2} y={IDEO_Y - 3} fontSize={8} textAnchor="middle" fill={C.faint}>{b.name}</text>)}
      {rec && bpKey.split(',').filter(Boolean).map(p => <path key={`ibp${p}`} d={`M${ix(+p)},${IDEO_Y - 1} l-3,-5 h6 z`} fill={C.lost} />)}
      <path d={`M${ix(vp.s)},${IDEO_Y + IDEO_H} L${X0},${IDEO_Y + IDEO_H + LENS_H} H${X1} L${ix(vp.s) + vl},${IDEO_Y + IDEO_H} Z`} fill="url(#svb-lens)" />
      <rect x={ix(vp.s) - 1} y={IDEO_Y - 2} width={vl + 2} height={IDEO_H + 4} rx={3} fill={C.accent} fillOpacity={0.12} stroke={C.accent} strokeWidth={1.4} />
      {/* track names */}
      <g fontSize={10.5} fill={C.muted}>
        {domains.length > 0 && <text x={14} y={yRefTad + TAD_H - 6}>TADs</text>}
        <text x={14} y={yRefGenes + 12} fontWeight={600} fill={C.text}>genes</text>
        <text x={14} y={yRefGenes + 25} fontSize={9.5}>{lod === 'density' ? 'density' : lod === 'bars' || mode === 'mane' ? 'MANE Select' : 'all NM_'}</text>
        {view.loading && <text x={14} y={yRefGenes + 38} fontSize={9} fill={C.faint}>loading…</text>}
        {yDgv >= 0 && <text x={14} y={yDgv + 8} fontWeight={600} fill={C.text}>DGV Gold<title>DGV Gold Standard: copy-number variants of the general population (curated, found in two studies and two samples at least)</title></text>}
        {rec && <text x={14} y={yRefBand + 14} fontSize={9.5}>reference</text>}
        {ders.map((d, i) => <g key={d.name}><text x={14} y={derY[i].band + 14} fontWeight={700} fill={C.accent}>{d.name}</text><text x={14} y={derY[i].band + 26} fontSize={9.5}>rebuilt</text></g>)}
      </g>
      {/* ruler: drag along it to zoom on a stretch */}
      <rect x={X0} y={yRuler} width={PW} height={RULER_H - 4} fill="#eef2f6" rx={4} style={{ cursor: 'col-resize' }}><title>Drag along the ruler to zoom on a stretch</title></rect>
      <g clipPath="url(#svb-clip)">
        {ticks.map(t => <g key={`tk${t}`}><line x1={refX(t)} x2={refX(t)} y1={yRuler + RULER_H - 9} y2={yRuler + RULER_H - 4} stroke={C.faint} />{refX(t) > X0 + 28 && refX(t) < X1 - 28 && <text x={refX(t)} y={yRuler + 11} fontSize={9} textAnchor="middle" fill={C.muted} fontFamily={MONO}>{tickText(t, step)}</text>}</g>)}
        {/* grid lines through the tracks */}
        {ticks.map(t => <line key={`gl${t}`} x1={refX(t)} x2={refX(t)} y1={yRuler + RULER_H} y2={rec ? yRefBand : H - 8} stroke={C.grid} strokeWidth={0.6} />)}
        {/* reference TADs */}
        {domains.filter(d => d.end > vp.s && d.start < vp.e).map((d, i) => {
          const a = cx(refX(d.start)), b = cx(refX(d.end)), h = Math.min(TAD_H - 6, (b - a) / 2);
          return <path key={`rt${d.start}`} d={`M${a},${yRefTad + TAD_H} L${(a + b) / 2},${yRefTad + TAD_H - h} L${b},${yRefTad + TAD_H} Z`} fill={i % 2 ? '#dfe5ee' : '#e9eef4'} stroke="#a9b6c6" strokeWidth={0.8}><title>{`domain ${d.name ?? ''}: ${chrom}:${fmt(d.start + 1)}-${fmt(d.end)} (${kb(d.end - d.start)})`}</title></path>;
        })}
        {/* breakpoints through the reference */}
        {rec && rec.segs.filter(s => !s.first && inView(refX(s.start), refX(s.start))).map(s => <line key={`bpl${s.id}`} x1={refX(s.start)} x2={refX(s.start)} y1={yRuler + RULER_H} y2={yRefBand + BAND_H} stroke={C.lost} strokeWidth={1} strokeDasharray="3 3" opacity={0.55} />)}
        {/* reference genes */}
        {density && <path d={density.d} transform={`translate(0 ${yDens + 34}) scale(1 34)`} fill="url(#svb-dens)" stroke="none"><title>{`genes per ${kb(3 / ppb)}: up to ${density.max}`}</title></path>}
        {refItems.items.map(({ g, lane, labeled }) => {
          const gi = impactOf.get(g.symbol), col = geneColor(g.symbol, gi);
          const yy = yRefGenes + refItems.laneY[lane];
          const dot = panelDot(g.symbol), lx = Math.max(X0 + 2, refX(g.start));
          return (
            <g key={`rg${g.symbol}`} onMouseEnter={hoverOn(g.symbol)} onMouseMove={hoverOn(g.symbol)} onMouseLeave={hoverOff} style={{ cursor: 'pointer' }} onClick={pinOn(g.symbol)}>
              {dot && <circle cx={lx + 3.5} cy={yy + 6.5} r={3.5} fill={dot} />}
              {labeled && <text x={lx + (dot ? 10 : 0)} y={yy + 10} fontSize={11} fontWeight={700} fontStyle="italic" fill={col}>{g.symbol}</text>}
              {lod === 'density'
                ? <line x1={lx} x2={lx} y1={yy + 12} y2={yDens + 6} stroke={col} strokeWidth={0.8} strokeOpacity={0.6} />
                : lod === 'bars'
                  ? <rect x={cx(refX(g.start))} y={yy + LABEL_H + 3} width={Math.max(1.5, cx(refX(g.end)) - cx(refX(g.start)))} height={6} rx={1.5} fill={col} />
                  : shown(g).map((t, i) => glyph(t, refX, yy + LABEL_H + 6 + i * TX_H, col, `r${t.id}`))}
              {gi?.kind === 'deleted' && lod !== 'density' && <line x1={cx(refX(g.start))} x2={cx(refX(g.end))} y1={yy + LABEL_H + 6} y2={yy + LABEL_H + 6} stroke={C.lost} strokeWidth={2.2} opacity={0.6} />}
            </g>
          );
        })}
        {/* population CNVs */}
        {yDgv >= 0 && dgvItems.map(({ v, lane, a, b, ia, ib }) => {
          const yy = yDgv + lane * DGV_LANE_H, col = DGV_COLOR[v.type], same = popSame.has(`${v.id}:${v.outer[0]}`), key = `dgv:${v.id}:${v.outer[0]}`;
          return (
            <g key={key} opacity={same ? 1 : 0.5} style={{ cursor: 'pointer' }} onMouseEnter={hoverOn(key, v)} onMouseMove={hoverOn(key, v)} onMouseLeave={hoverOff} onClick={pinOn(key, v)}>
              <rect x={a} y={yy} width={b - a} height={DGV_LANE_H - 1} fill="transparent" />
              <path d={`M${a},${yy + 4}H${b}M${a},${yy + 2}v4M${b},${yy + 2}v4`} stroke={col} strokeWidth={1} fill="none" />
              <rect x={ia} y={yy + 1.5} width={Math.max(1.5, ib - ia)} height={5} fill={col} stroke={same ? '#111827' : 'none'} strokeWidth={0.6} />
            </g>
          );
        })}
        {yDgv >= 0 && !dgvItems.length && <text x={X0 + 6} y={yDgv + 8} fontSize={9.5} fill={C.faint}>{span > DGV_MAX_SPAN ? `zoom in below ${kb(DGV_MAX_SPAN)} to see the population CNVs` : dgvView.loading ? 'loading the DGV Gold Standard…' : dgvView.error ? `DGV Gold Standard unavailable: ${dgvView.error}` : 'no DGV Gold Standard CNV here'}</text>}
        {/* reference band: the event's segments */}
        {rec && rec.segs.map(s => {
          const [lo, hi] = [s.first ? 0 : s.start, s.last ? L : s.end];
          const a = cx(refX(lo)), b = cx(refX(hi)), col = segColors.get(s.id)!, lost = rec.lost.includes(s);
          if (b < X0 - 2 || a > X1 + 2) return null;
          const va = Math.max(a, X0), vb = Math.min(b, X1);
          return (
            <g key={`rb${s.id}`}>
              <rect x={a} y={yRefBand} width={Math.max(1, b - a)} height={BAND_H} fill={lost ? 'url(#svb-lost)' : col} fillOpacity={lost ? 1 : 0.85} stroke="#fff" strokeWidth={1} />
              {vb - va > 14 && <text x={(va + vb) / 2} y={yRefBand + 14} fontSize={11.5} fontWeight={700} textAnchor="middle" fill={lost ? C.lost : '#fff'}>{s.letter}</text>}
              {lost && vb - va > 30 && <text x={(va + vb) / 2} y={yRefBand + BAND_H + 13} fontSize={10} fontWeight={700} textAnchor="middle" fill={C.lost}>lost · {kb(s.end - s.start)}</text>}
            </g>
          );
        })}
        {/* ribbons */}
        {ders.map((d, di) => d.pieces.map(dp => {
          const [lo, hi] = pieceRange(dp, L), col = segColors.get(dp.seg.id)!;
          const ra = cx(refX(lo)), rb = cx(refX(hi));
          let da = cx(derX(d, toDer(dp, lo))), db = cx(derX(d, toDer(dp, hi)));
          if (dp.piece.reversed) [da, db] = [db, da];
          if (Math.max(ra, rb, da, db) < X0 - 2 || Math.min(ra, rb, da, db) > X1 + 2) return null;
          const y0 = yRefBand + BAND_H, y1 = derY[di].band, m = (y0 + y1) / 2;
          return <path key={`rib${di}-${dp.index}`} d={`M${ra},${y0} C${ra},${m} ${da},${m} ${da},${y1} L${db},${y1} C${db},${m} ${rb},${m} ${rb},${y0} Z`} fill={col} fillOpacity={0.16} stroke={col} strokeOpacity={0.45} strokeWidth={0.9}><title>{`${dp.seg.letter}${dp.piece.reversed ? '′ (inverted)' : ''}: ${chrom}:${fmt(dp.seg.start + 1)}-${fmt(dp.seg.end)} → ${d.name}${dp.piece.copy ? ` (copy ${dp.piece.copy + 1})` : ''}`}</title></path>;
        }))}
        {/* rebuilt chromosomes */}
        {ders.map((d, di) => {
          const Y = derY[di], packed = derPacked[di];
          return (
            <g key={`der${d.name}`}>
              {d.pieces.map(dp => {
                const [lo, hi] = pieceRange(dp, L);
                let a = cx(derX(d, toDer(dp, lo))), b = cx(derX(d, toDer(dp, hi)));
                if (a > b) [a, b] = [b, a];
                if (b < X0 - 2 || a > X1 + 2) return null;
                const col = segColors.get(dp.seg.id)!, va = Math.max(a, X0), vb = Math.min(b, X1);
                return (
                  <g key={`db${dp.index}`}>
                    <rect x={a} y={Y.band} width={Math.max(1, b - a)} height={BAND_H} fill={col} fillOpacity={0.85} stroke="#fff" strokeWidth={1} />
                    {vb - va > 18 && <text x={(va + vb) / 2} y={Y.band + 14} fontSize={11.5} fontWeight={700} textAnchor="middle" fill="#fff">{dp.piece.reversed ? `◂ ${dp.seg.letter}′` : dp.seg.letter}</text>}
                  </g>
                );
              })}
              {d.pieces.slice(1).map((dp, i) => {
                const x = derX(d, dp.derStart), via = dp.piece.via, alt = i % 2;
                if (x < X0 - 2 || x > X1 + 2) return null;
                if (!via || via === 'ref') return <line key={`js${dp.index}`} x1={x} x2={x} y1={Y.band} y2={Y.band + BAND_H} stroke="#fff" strokeWidth={1.5} />;
                const col = KIND_COLOR[via.kind];
                return (
                  <g key={`js${dp.index}`}>
                    <line x1={x} x2={x} y1={Y.band - 3} y2={Y.genes + packed.height} stroke={col} strokeWidth={1} strokeDasharray="2 3" opacity={0.6} />
                    <circle cx={x} cy={Y.band + BAND_H / 2} r={5} fill={col} stroke="#fff" strokeWidth={1.8} />
                    <text x={x} y={Y.ruler + (alt ? 18 : 9)} fontSize={9.5} fontWeight={700} textAnchor="middle" fill={col}>{via.id}{via.reads ? ` · ${fmt(via.reads)}` : ''}</text>
                  </g>
                );
              })}
              {packed.items.map((b, bi) => {
                const yy = Y.genes + packed.laneY[b.lane];
                const gi = impactOf.get(b.genes[0]);
                const isFusion = b.genes.length > 1;
                const head = isFusion && colorBy === 'impact' ? C.fusion5 : geneColor(b.genes[0], gi);
                const lx = Math.max(X0 + 2, b.x0);
                const dot = !isFusion ? panelDot(b.genes[0]) : null;
                const yl = yy + LABEL_H + 6;
                return (
                  <g key={`dg${b.flank ? `f${b.genes[0]}` : bi}`} onMouseEnter={hoverOn(b.genes[0])} onMouseMove={hoverOn(b.genes[0])} onMouseLeave={hoverOff} style={{ cursor: 'pointer' }} onClick={pinOn(b.genes[0])}>
                    {dot && <circle cx={lx + 3.5} cy={yy + 6.5} r={3.5} fill={dot} />}
                    {b.labeled && <text x={lx + (dot ? 10 : 0)} y={yy + 10} fontSize={11} fontWeight={700} fontStyle="italic" fill={head}>{b.label}</text>}
                    {b.labeled && isFusion && b.fusion && <g>
                      <rect x={lx + labelW(b.label) - 4} y={yy + 1} width={b.fusion.kind === 'in-frame' ? 52 : b.fusion.kind === 'promoter-swap' ? 74 : 66} height={12} rx={6} fill={b.fusion.kind === 'in-frame' ? '#7c3aed' : '#ede9fe'} />
                      <text x={lx + labelW(b.label) + 2} y={yy + 10} fontSize={9} fontWeight={700} fill={b.fusion.kind === 'in-frame' ? '#fff' : '#6d28d9'}>{b.fusion.kind === 'in-frame' ? 'in frame' : b.fusion.kind === 'promoter-swap' ? 'promoter swap' : b.fusion.kind === 'out-of-frame' ? 'out of frame' : b.fusion.kind === 'five-cds-whole' ? '3′ UTR join' : 'frame ?'}</text>
                    </g>}
                    {b.flank && (() => {
                      const { g, dp } = b.flank, X = (pos: number) => derX(d, toDer(dp, pos));
                      return lod === 'full' ? glyph(g.principal, X, yl, head, `fl${g.symbol}`, g.principal.start, g.principal.end, (dp.piece.reversed ? -g.strand : g.strand) as 1 | -1)
                        : <rect x={cx(Math.min(X(g.start), X(g.end)))} y={yy + LABEL_H + 3} width={Math.max(1.5, Math.abs(cx(X(g.end)) - cx(X(g.start))))} height={6} rx={1.5} fill={head} />;
                    })()}
                    {b.chains.map((ch, ci) => {
                      const ylc = yl + ci * TX_H;
                      if (lod !== 'full') {
                        const xs = ch.parts.flatMap(p => [derX(d, toDer(p.dp, p.lo)), derX(d, toDer(p.dp, p.hi))]);
                        return <rect key={`ch${ci}`} x={cx(Math.min(...xs))} y={yy + LABEL_H + 3} width={Math.max(1.5, cx(Math.max(...xs)) - cx(Math.min(...xs)))} height={6} rx={1.5} fill={head} />;
                      }
                      return (
                        <g key={`ch${ci}`}>
                          {ch.parts.map((p: Fragment, pi) => {
                            const col = isFusion ? (colorBy === 'loeuf' ? geneColor(p.tx.gene, impactOf.get(p.tx.gene)) : p.tx.gene === b.genes[0] ? C.fusion5 : C.fusion3) : head;
                            const X = (pos: number) => derX(d, toDer(p.dp, pos));
                            const els = [glyph(p.tx, X, ylc, col, `d${pi}`, p.lo, p.hi, p.strand)];
                            const leftX = Math.min(X(p.lo), X(p.hi)), rightX = Math.max(X(p.lo), X(p.hi));
                            const cutL = p.dp.piece.reversed ? p.cutHi && !p.edgeHi : p.cutLo && !p.edgeLo, cutR = p.dp.piece.reversed ? p.cutLo && !p.edgeLo : p.cutHi && !p.edgeHi;
                            const isFirst = (ch.parts[0].strand > 0 ? pi === 0 : pi === ch.parts.length - 1), isLast = (ch.parts[0].strand > 0 ? pi === ch.parts.length - 1 : pi === 0);
                            if (cutL && (ch.parts[0].strand > 0 ? isFirst : isLast)) els.push(torn(leftX, ylc, C.lost, `tl${pi}`));
                            if (cutR && (ch.parts[0].strand > 0 ? isLast : isFirst)) els.push(torn(rightX, ylc, C.lost, `tr${pi}`));
                            return <g key={`p${pi}`}>{els}</g>;
                          })}
                          {isFusion && (() => {
                            const xs = ch.parts.flatMap(p => [derX(d, toDer(p.dp, p.lo)), derX(d, toDer(p.dp, p.hi))]);
                            return <rect x={cx(Math.min(...xs)) - 3} y={ylc - 6.5} width={cx(Math.max(...xs)) - cx(Math.min(...xs)) + 6} height={13} rx={6.5} fill="none" stroke={C.fusion5} strokeWidth={1.1} strokeDasharray="4 2" />;
                          })()}
                        </g>
                      );
                    })}
                  </g>
                );
              })}
              {domains.length > 0 && story.derDomains.filter(x => story.ders[x.der] === d).map((x, i) => {
                const a = cx(derX(d, x.derStart)), b = cx(derX(d, x.derEnd)), h = Math.min(TAD_H - 8, (b - a) / 2), top = Y.tad + 2;
                if (b < X0 || a > X1) return null;
                return (
                  <g key={`dt${i}`}>
                    <path d={`M${a},${top} L${(a + b) / 2},${top + h} L${b},${top} Z`} fill={x.neo ? '#ede9fe' : i % 2 ? '#dfe5ee' : '#e9eef4'} stroke={x.neo ? '#7c3aed' : '#a9b6c6'} strokeWidth={x.neo ? 1.4 : 0.8} strokeDasharray={x.neo ? '4 2' : undefined}>
                      <title>{x.neo ? `new domain: ${x.parts.map(p => `${kb(p.hi - p.lo)} of domain ${p.domain + 1}`).join(' + ')}` : `domain ${x.parts[0].domain + 1}`}</title>
                    </path>
                    {x.neo && b - a > 60 && <text x={(a + b) / 2} y={top + Math.min(h, 22) - 4} fontSize={9.5} fontWeight={700} textAnchor="middle" fill="#6d28d9">new domain</text>}
                  </g>
                );
              })}
            </g>
          );
        })}
      </g>
      {domains.length > 0 && ders.map((d, i) => <text key={`tadl${i}`} x={14} y={derY[i].tad + 16} fontSize={10.5} fill={C.muted}>TADs</text>)}
      {/* the stretch being chosen on the ruler */}
      {band && <g pointerEvents="none">
        <rect x={Math.min(band.a, band.b)} y={yRuler} width={Math.abs(band.b - band.a)} height={H - yRuler - 6} fill={C.accent} fillOpacity={0.08} stroke={C.accent} strokeDasharray="4 3" />
        <text x={(band.a + band.b) / 2} y={yRuler - 4} fontSize={10} fontWeight={700} textAnchor="middle" fill={C.accent}>{kb(Math.abs(band.b - band.a) / ppb)}</text>
      </g>}
      {/* the pointer's position, moved directly (data-export="skip": not in the figure's export) */}
      <g ref={guideRef} style={{ display: 'none' }} pointerEvents="none" data-export="skip">
        <line x1={0} x2={0} y1={yRuler + RULER_H} y2={H - 8} stroke={C.accent} strokeWidth={0.8} strokeOpacity={0.5} />
        <text x={4} y={yRuler + RULER_H + 10} fontSize={9.5} fill={C.accent} fontFamily={MONO} />
      </g>
    </svg>
  );

  // ---- the card of the gene or population CNV under the pointer, or kept open ----
  const card = (at: TipAt, isPin: boolean) => {
    const place: React.CSSProperties = { left: Math.max(0, Math.min(at.x + 14, (boxRef.current?.clientWidth ?? 900) - 400)), ...(at.y > (boxRef.current?.clientHeight ?? 400) - 160 ? { bottom: (boxRef.current?.clientHeight ?? 400) - at.y + 14 } : { top: at.y + 14 }) };
    const ref = isPin ? undefined : (el: HTMLDivElement | null) => { tipEl.current = el; };
    if (at.dgv) return <div key={isPin ? 'pin' : 'tip'} ref={ref} className={isPin ? '' : 'absolute z-20 pointer-events-none'}><DgvCard v={at.dgv} chrom={chrom} cnv={cnv} pinned={isPin} place={isPin ? place : {}} onClose={() => setPinned(null)} onZoom={() => zoomTo(at.dgv!.outer[0], at.dgv!.outer[1])} /></div>;
    const gi = impactOf.get(at.gene), g = gi?.gene ?? geneBySym.get(at.gene);
    if (!g) return null;
    const r = gi ? dosageReading(gi, annotations) : null;
    const c2 = ann.constraint.get(g.symbol), d2 = ann.dosage.get(g.symbol), pg = panel?.genes.get(g.symbol.toUpperCase());
    const parts = [
      c2 ? `gnomAD ${c2.version}: LOEUF ${c2.loeuf != null ? c2.loeuf.toFixed(2) : 'n/a'} · pLI ${c2.pli != null ? c2.pli.toFixed(2) : 'n/a'}` : 'gnomAD: no constraint (yet)',
      d2 ? `ClinGen HI ${d2.hi ?? '—'} · TS ${d2.ts ?? '—'}` : '',
      pg ? `${panel!.name}: ${CONF[pg.confidence].label}` : '',
    ].filter(Boolean);
    const dot = { intolerant: '#fca5a5', constrained: '#fdba74', tolerant: '#cbd5e1', unknown: '#cbd5e1' }[loeufClass(c2)];
    const tx = gi?.principal.tx ?? g.principal;
    return (
      <div key={isPin ? 'pin' : 'tip'} ref={ref} className={`absolute z-20 rounded-lg bg-gray-900 text-white text-[11px] px-3 py-2 shadow-xl max-w-[390px] ${isPin ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={isPin ? place : undefined}
        onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
        <div className="font-semibold flex items-start gap-2"><span className="flex-1"><i>{g.symbol}</i>{gi ? ` · ${IMPACT[gi.kind].label}` : ''}</span>
          {isPin && <button onClick={() => setPinned(null)} className="-mr-1 -mt-0.5 px-1 text-gray-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
        <div className="text-gray-300">{tx.id}{tx.status === 'mane' ? ' (MANE Select)' : tx.status === 'select' ? ' (RefSeq Select)' : ''} · {chrom}:{fmt(g.start + 1)}-{fmt(g.end)} ({kb(g.end - g.start)}, {g.strand > 0 ? '+' : '−'} strand)</div>
        {gi && <div className="mt-0.5">{gi.principal.text}</div>}
        {r?.text && <div className="mt-0.5 text-amber-200">{r.text}</div>}
        <div className="mt-0.5 text-gray-300"><span style={{ color: dot }}>●</span> {parts.join(' · ')}</div>
        {isPin
          ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={() => zoomTo(g.start, g.end)} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on <i>{g.symbol}</i></button></div>
          : <div className="mt-0.5 text-gray-400">click to keep this card open</div>}
      </div>
    );
  };

  const btn = 'min-w-[28px] px-1.5 h-6 rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 font-semibold';
  return (
    <div className="mt-2 border-t border-gray-200 pt-3 outline-none" style={{ fontFamily: FONT }} tabIndex={0} onKeyDown={onKeyDown}>
      {/* header */}
      <div className="px-5 flex flex-wrap items-center gap-2 text-xs">
        <div className="text-sm font-bold mr-1">{rec ? 'Genes on the rebuilt chromosome' : 'Genes'}</div>
        <span className="inline-flex rounded-md border border-gray-300 overflow-hidden" role="group" aria-label="Transcripts">
          {(['mane', 'all'] as const).map(m => <button key={m} onClick={() => setMode(m)} className={`px-2 py-0.5 ${mode === m ? 'bg-accent-600 text-white' : 'bg-white text-gray-600 hover:bg-gray-50'}`} title={m === 'mane' ? 'One model per gene: MANE Select (else RefSeq Select, else the longest coding sequence)' : `Every RefSeq NM_ model of each gene (no XM_ prediction, no NR_ non-coding): the consequence can differ between isoforms. Drawn below ${kb(LOD_FULL)} in view.`}>{m === 'mane' ? 'MANE Select' : 'all NM_ transcripts'}</button>)}
        </span>
        <span className="inline-flex rounded-md border border-gray-300 overflow-hidden" role="group" aria-label="Gene colour">
          {(['loeuf', 'impact'] as const).map(m => <button key={m} onClick={() => { setColorBy(m); store.set('sashimi.svColorBy', m); }} className={`px-2 py-0.5 ${colorBy === m ? 'bg-accent-600 text-white' : 'bg-white text-gray-600 hover:bg-gray-50'}`} title={m === 'loeuf' ? `Genes coloured by their gnomAD LoF constraint (LOEUF): below ${LOEUF_INTOLERANT.v4} red (Mendelian threshold), ${LOEUF_INTOLERANT.v4}–${LOEUF_CONSTRAINED.v4} orange (gnomAD v4 general cut-off), above grey; v2.1.1 scores: below ${LOEUF_INTOLERANT['v2.1.1']} red` : 'Genes coloured by what the rearrangement does to them'}>{m === 'loeuf' ? 'colour: LOEUF' : 'colour: impact'}</button>)}
        </span>
        <PanelPicker panel={panelData} onPanel={p => { setPanelData(p); store.set('sashimi.svPanel', p); }} />
        {loaders.getPopulationSvs && (
          <button onClick={() => { setShowDgv(v => { store.set('sashimi.svDgv', !v); return !v; }); }} aria-pressed={showDgv}
            className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-md border ${showDgv ? 'border-accent-300 bg-accent-50 text-accent-800' : 'border-gray-300 bg-white text-gray-600'} hover:bg-gray-50`}
            title="DGV Gold Standard (UCSC track dgvGold): curated copy-number variants of the general population, each found in two studies and two samples at least; drawn under the genes and compared with the event's losses and gains">
            ▤ DGV Gold
          </button>
        )}
        {rec && <label className="inline-flex items-center gap-1 px-2 py-0.5 rounded-md border border-gray-300 bg-white hover:bg-gray-50 cursor-pointer" title="Chromatin domains (TADs) as a BED file (chrom, start, end, name; 0-based, the same genome build), e.g. from Hi-C of a relevant tissue. They show which genes change domain: an intact gene joined to another domain may be misregulated (position effect).">
          ◭ {domainsFile ? `TADs: ${domainsFile.name} (${domains.length})` : 'Add chromatin domains (BED)…'}
          <input type="file" accept=".bed,.txt,.tsv" className="hidden" onChange={async e => { const f = e.target.files?.[0]; if (f) { const v = { name: f.name, text: await f.text() }; setDomainsFile(v); store.set('sashimi.svDomains', v.text.length < 3e6 ? v : null); } e.target.value = ''; }} />
        </label>}
        {rec && domainsFile && <button onClick={() => { setDomainsFile(null); store.set('sashimi.svDomains', null); }} className="text-gray-400 hover:text-red-600" title="Remove the domains">×</button>}
        {rec && <label className="inline-flex items-center gap-1 text-gray-600" title={`Without chromatin domains, an intact gene within this distance of a breakpoint is reported as a possible position effect. TADs average about 1 Mb; position effects have been documented 1.3 Mb from a breakpoint (SOX9) and further.${domains.length ? ' With domains loaded, the domains are used instead.' : ''}`}>
          near a breakpoint ≤
          <select value={near} onChange={e => { const v = Number(e.target.value); setNear(v); store.set('sashimi.svNear', v); }} className="border rounded px-1 py-0.5">
            {[0, 250e3, 500e3, 1e6, 2e6, 3e6].map(v => <option key={v} value={v}>{v === 0 ? 'off' : v >= 1e6 ? `${v / 1e6} Mb` : `${v / 1e3} kb`}</option>)}
          </select>
        </label>}
        <span className="ml-auto flex items-center gap-1">
          <span className="font-mono text-[11px] text-gray-600 mr-1 tabular-nums" title="The stretch in view">{chrom}:{fmt(vp.s + 1)}-{fmt(vp.e)} <span className="text-gray-400">· {kb(span)}</span></span>
          <button className={btn} onClick={() => flyTo(vp.s - span * 0.5, vp.e + span * 0.5, 300)} title="Zoom out (−, Ctrl + wheel, Shift + double-click)">−</button>
          <button className={btn} onClick={() => flyTo(vp.s + span * 0.25, vp.e - span * 0.25, 300)} title="Zoom in (+, Ctrl + wheel, double-click)">+</button>
          {rec && <button className={btn} onClick={eventView} title="Fit the event (E)">event</button>}
          <button className={btn} onClick={() => flyTo(0, L, 700)} title="The whole chromosome">{chrom.replace(/^chr/, 'chr ')}</button>
        </span>
      </div>
      {rec && showDgv && loaders.getPopulationSvs && <DgvSummary cnv={cnv} loading={pop.loading} error={pop.error} />}
      <div className="px-3 mt-2"><div ref={widthBox} className="w-full" /></div>
      <div className="px-3 mt-0 relative" ref={boxRef}>
        {svg}
        {pinned && card(pinned, true)}
        {tip && tip.gene !== pinned?.gene && card(tip, false)}
      </div>
      {/* legend */}
      <div className="px-5 mt-1 flex flex-wrap items-center gap-x-3 gap-y-1 text-[10.5px] text-gray-500">
        {colorBy === 'impact'
          ? (['fusion', 'disrupted', 'deleted', 'intragenic-del', 'intragenic-dup', 'partial-copy', 'duplicated', 'domain', 'moved', 'nearby', 'unaffected'] as ImpactKind[]).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: IMPACT[k2].color }} />{IMPACT[k2].label}</span>)
          : (['intolerant', 'constrained', 'tolerant', 'unknown'] as const).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: LOEUF_CLASS[k2].color }} />{LOEUF_CLASS[k2].label('v4')}</span>)}
        {colorBy === 'loeuf' && <span>(gnomAD v4; v2.1.1 scores: &lt; {LOEUF_INTOLERANT['v2.1.1']} intolerant)</span>}
        {rec && <span className="inline-flex items-center gap-1"><svg width="10" height="14"><path d="M6,1 l-3,3 l3,3 l-3,3 l3,3" fill="none" stroke={C.lost} strokeWidth="1.4" /></svg>cut by a breakpoint</span>}
        {showDgv && loaders.getPopulationSvs && <span className="inline-flex items-center gap-1"><svg width="26" height="10"><line x1="1" x2="25" y1="5" y2="5" stroke={DGV_COLOR.loss} /><line x1="1" x2="1" y1="2" y2="8" stroke={DGV_COLOR.loss} /><line x1="25" x2="25" y1="2" y2="8" stroke={DGV_COLOR.loss} /><rect x="8" y="2.5" width="10" height="5" fill={DGV_COLOR.loss} /></svg>DGV Gold CNV: outer extent (line), inner (bar); <span style={{ color: DGV_COLOR.loss }}>loss</span> <span style={{ color: DGV_COLOR.gain }}>gain</span> <span style={{ color: DGV_COLOR.both }}>both</span>{rec ? '; solid and outlined when it matches a loss or gain of the event' : ''}</span>}
        {panel && <span className="inline-flex items-center gap-1">● panel rating: <span style={{ color: CONF[3].color }}>green</span> <span style={{ color: CONF[2].color }}>amber</span> <span style={{ color: CONF[1].color }}>red</span></span>}
        <span>· drag to pan (it glides on), Ctrl / ⌘ + scroll or pinch to zoom, two-finger swipe or Shift + scroll to pan, double-click to zoom in (Shift: out), drag along the ruler to zoom on a stretch, click the ideogram to go there; ← → + − {rec ? 'E ' : ''}on the keyboard; hover a gene for its card, click to keep it open</span>
      </div>
      <div className="px-5 mt-2 pb-1">
        {rec && !story.genes.some(g => g.kind !== 'unaffected') && !txs.loading && <div className="text-gray-400 text-xs py-1">No gene of the window is changed by the selected junctions.</div>}
        <div className="text-[10.5px] text-gray-500 mt-2 leading-4">
          Genes: RefSeq NM_ models (UCSC ncbiRefSeqCurated), MANE status from the UCSC mane track, loaded by tiles around the view{view.error ? <span className="text-amber-700"> · {view.error}</span> : null}{txs.error ? <span className="text-amber-700"> · {txs.error}</span> : null}. {ann.loading ? 'Annotations loading… ' : ann.notes.length ? `${ann.notes.join(' · ')}. ` : ''}
          {rec && <>Copies count the {zygosity === 'het' ? `${flankCopies - 1} intact homologue${flankCopies - 1 === 1 ? '' : 's'} and the rearranged one` : 'rearranged homologues'}; a fusion needs both partners on the same strand of the rebuilt chromosome, its frame is read from the coding phases of the two exons joined (breakpoints in introns); a breakpoint inside an exon leaves the frame undetermined. </>}
          LoF-intolerant: LOEUF &lt; {LOEUF_INTOLERANT.v4} (gnomAD v4, the browser's threshold for Mendelian cases; &lt; {LOEUF_INTOLERANT['v2.1.1']} for v2.1.1) or pLI ≥ {PLI_INTOLERANT}; constrained: LOEUF &lt; {LOEUF_CONSTRAINED.v4} (gnomAD v4.0's general cut-off). LOEUF is best read as a continuous value: hover a gene for it. ClinGen scores: 3 sufficient, 2 some, 1 little, 0 no evidence, 30 autosomal recessive, 40 dosage sensitivity unlikely.
          {rec && (domains.length ? ' Domains: a junction joining two parts of domains makes one new domain (boundaries between them lost); genes whose domain gains foreign sequence or loses part of its own are flagged, whatever the size.' : near > 0 ? ` No domains loaded: intact genes within ${near >= 1e6 ? `${near / 1e6} Mb` : `${near / 1e3} kb`} of a breakpoint are flagged as a possible position effect; those in the panel, with a ClinGen score of 2 or 3, or LoF-intolerant keep their label in the figure.` : '')}
          {' '}Far out, a gene is drawn from its principal model only (bars past {kb(LOD_FULL)} in view, a density profile past {kb(LOD_BARS)}, with the genes worth a look named).
        </div>
      </div>
    </div>
  );
}

// ======================== Genes worth a look ========================

/** In the chosen panel, a ClinGen haploinsufficiency or triplosensitivity score of 2 or 3, or LoF-intolerant in gnomAD. */
function notableSym(sym: string, a: GeneAnnotations): boolean {
  const d = a.dosage.get(sym), c = a.constraint.get(sym);
  return !!a.panel?.genes.has(sym.toUpperCase()) || [2, 3].includes(d?.hi ?? -1) || [2, 3].includes(d?.ts ?? -1)
    || (c?.loeuf != null && c.loeuf < (LOEUF_INTOLERANT[c.version] ?? LOEUF_INTOLERANT.v4)) || (c?.pli != null && c.pli >= PLI_INTOLERANT);
}

// ======================== PanelApp picker ========================

const TYPE_LABEL = { loss: 'loss', gain: 'gain', both: 'gain and loss', unknown: 'type not given' } as const;
const pct = (v: number) => `${Math.round(v * 100)} %`;
const evName = (ev: CnvEvent) => `${ev.kind} ${ev.letters} (${kb(ev.end - ev.start)})`;

/** The card of a population CNV: its extents, its frequency, its overlap with the event's losses and gains, and every field the source gives. */
function DgvCard({ v, chrom, cnv, pinned, place, onClose, onZoom }: { v: PopulationSv; chrom: string; cnv: { ev: CnvEvent; matches: PopMatch[] }[]; pinned: boolean; place: React.CSSProperties; onClose: () => void; onZoom: () => void }) {
  const span = (r: [number, number]) => `${chrom}:${fmt(r[0] + 1)}-${fmt(r[1])} (${kb(r[1] - r[0])})`;
  const sameInner = v.inner[0] === v.outer[0] && v.inner[1] === v.outer[1];
  const rows = cnv.map(c => ({ ev: c.ev, m: c.matches.find(m => m.sv.id === v.id && m.sv.outer[0] === v.outer[0]) })).filter(x => x.m);
  return (
    <div className={`absolute z-20 rounded-md bg-gray-900 text-white text-[11px] px-2.5 py-1.5 shadow-lg max-w-[400px] ${pinned ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={place}
      onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
      <div className="font-semibold flex items-start gap-2"><span className="flex-1"><span style={{ color: { loss: '#fca5a5', gain: '#93c5fd', both: '#c4b5fd', unknown: '#cbd5e1' }[v.type] }}>●</span> {v.id} · population {TYPE_LABEL[v.type]} (DGV Gold Standard)</span>
        {pinned && <button onClick={onClose} className="-mr-1 -mt-0.5 px-1 text-gray-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
      <div className="text-gray-300">{sameInner ? span(v.outer) : <>outer {span(v.outer)}<br />inner {span(v.inner)}</>}</div>
      <div className="mt-0.5">{[v.freqText != null ? `frequency ${v.freqText}` : 'frequency not given', v.samples != null ? `${fmt(v.samples)} samples` : '', v.studies != null ? `${fmt(v.studies)} studies` : ''].filter(Boolean).join(' · ')}</div>
      {rows.map(({ ev, m }) => (
        <div key={ev.start} className={`mt-0.5 ${m!.sameType && m!.ro >= RO_MATCH ? 'text-emerald-300' : 'text-gray-300'}`}>
          vs the {evName(ev)}: reciprocal overlap {pct(m!.ro)}{sameInner ? '' : ` (inner ${pct(m!.roInner)}, outer ${pct(m!.roOuter)})`}{m!.covers ? ', its inner extent holds the whole event' : ''}{m!.sameType ? '' : ', other type'}
        </div>
      ))}
      {pinned && v.fields.length > 0 && (
        <details className="mt-1 text-gray-400"><summary className="cursor-pointer">all fields ({v.fields.length})</summary>
          <div className="mt-0.5 max-h-40 overflow-auto font-mono text-[10px] leading-snug">{v.fields.map(([k2, val]) => <div key={k2}><span className="text-gray-500">{k2}</span> {val.length > 160 ? `${val.slice(0, 160)}…` : val}</div>)}</div>
        </details>
      )}
      {pinned
        ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={onZoom} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on it</button></div>
        : <div className="mt-0.5 text-gray-400">click to keep this card open</div>}
    </div>
  );
}

/**
 * One line per loss or gain of the event: the population CNVs of the same type with a reciprocal overlap of RO_MATCH or
 * more (the same event seen in the general population: evidence towards benign), else the closest one.
 */
function DgvSummary({ cnv, loading, error }: { cnv: { ev: CnvEvent; matches: PopMatch[] }[]; loading: boolean; error?: string }) {
  if (!cnv.length) return null;
  const tip = `DGV Gold Standard: curated CNV regions of the general population, each found in two studies and two samples at least. A loss or gain counts as seen when a population CNV of the same type overlaps it reciprocally by ${pct(RO_MATCH)} or more (inner or outer extent). Seen in the population is evidence towards benign; not seen is not evidence of pathogenicity.`;
  return (
    <div className="px-5 mt-1.5 flex flex-wrap gap-1.5 text-[11.5px]" title={tip}>
      <span className="text-gray-500 self-center">DGV Gold:</span>
      {loading ? <span className="text-gray-500">loading…</span> : error ? <span className="text-gray-500">unavailable ({error})</span> : cnv.map(({ ev, matches }) => {
        const same = matches.filter(m => m.sameType && m.ro >= RO_MATCH);
        const best = matches.find(m => m.sameType);
        const seen = same.length > 0, top = same.sort((a, b) => (b.sv.freq ?? 0) - (a.sv.freq ?? 0))[0];
        return (
          <span key={ev.start} className={`inline-flex items-center gap-1 rounded-full px-2.5 py-0.5 border ${seen ? 'bg-emerald-50 border-emerald-300 text-emerald-800' : 'bg-gray-50 border-gray-300 text-gray-700'}`}>
            <b>{evName(ev)}</b>
            {seen
              ? <>seen in the population: {same.length} CNV{same.length > 1 ? 's' : ''} ≥ {pct(RO_MATCH)} overlap{top ? `, e.g. ${top.sv.id}${top.sv.freqText != null ? ` (frequency ${top.sv.freqText})` : ''}` : ''}</>
              : <>no population {ev.kind} ≥ {pct(RO_MATCH)} overlap{best ? `; closest ${best.sv.id}, ${pct(best.ro)}` : ''}</>}
          </span>
        );
      })}
    </div>
  );
}


function PanelPicker({ panel, onPanel }: { panel: LoadedPanel | null; onPanel: (p: LoadedPanel | null) => void }) {
  const [open, setOpen] = useState(false);
  const [src, setSrc] = useState(PANEL_SOURCES[0].id);
  const [q, setQ] = useState('');
  const [list, setList] = useState<{ items: PanelSummary[]; loading: boolean; error?: string }>({ items: [], loading: false });
  const [busy, setBusy] = useState<string | null>(null);
  const source = PANEL_SOURCES.find(s => s.id === src)!;
  useEffect(() => {
    if (!open) return;
    let alive = true;
    setList({ items: [], loading: true });
    listPanels(source).then(items => { if (alive) setList({ items, loading: false }); }).catch(e => { if (alive) setList({ items: [], loading: false, error: `${e?.message || e}: ${source.label} could not be reached from this page (network, proxy). Load the panel from a file instead.` }); });
    return () => { alive = false; };
  }, [open, src]); // eslint-disable-line react-hooks/exhaustive-deps
  const words = q.toLowerCase().split(/\s+/).filter(Boolean);
  const hits = list.items.filter(p => words.every(w => p.name.toLowerCase().includes(w) || p.disorders.some(d => d.toLowerCase().includes(w)))).slice(0, 80);
  return (
    <span className="relative">
      <button onClick={() => setOpen(o => !o)} className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-md border ${panel ? 'border-emerald-300 bg-emerald-50 text-emerald-800' : 'border-gray-300 bg-white text-gray-700'} hover:bg-gray-50`}
        title="Rate the genes against a PanelApp gene panel (green / amber / red)">
        ▤ {panel ? `${panel.name}${panel.version ? ` v${panel.version}` : ''} · ${panel.genes.length} genes` : 'Gene panel (PanelApp)…'}
      </button>
      {panel && <button onClick={() => onPanel(null)} className="ml-1 text-gray-400 hover:text-red-600" title="Remove the panel">×</button>}
      {open && (
        <div className="absolute z-40 mt-1 w-[460px] rounded-lg border border-gray-300 bg-white shadow-2xl p-2 text-xs" onMouseDown={e => e.stopPropagation()}>
          <div className="flex gap-1.5">
            <select value={src} onChange={e => setSrc(e.target.value)} className="border rounded px-1 py-1">{PANEL_SOURCES.map(s => <option key={s.id} value={s.id}>{s.label}</option>)}</select>
            <input autoFocus value={q} onChange={e => setQ(e.target.value)} placeholder="search a panel or a disorder" className="flex-1 border rounded px-2 py-1 outline-none focus:border-accent-300" />
            <button onClick={() => setOpen(false)} className="px-1.5 text-gray-400 hover:text-gray-700">×</button>
          </div>
          <div className="mt-1.5 max-h-[260px] overflow-y-auto">
            {list.loading && <div className="text-gray-500 p-2">loading the panels of {source.label}…</div>}
            {list.error && <div className="text-amber-700 p-2">{list.error}</div>}
            {hits.map(p => (
              <button key={p.id} disabled={!!busy} onClick={async () => { setBusy(p.id); try { const lp = await loadPanel(source, p.id); onPanel(lp); setOpen(false); } catch (e: any) { setList(l => ({ ...l, error: `${p.name}: ${e?.message || e}` })); } finally { setBusy(null); } }}
                className="w-full text-left px-2 py-1 rounded hover:bg-accent-50 flex justify-between gap-2">
                <span className="truncate">{p.name}</span><span className="text-gray-400 whitespace-nowrap">{busy === p.id ? 'loading…' : `v${p.version} · ${p.genes} genes`}</span>
              </button>
            ))}
            {!list.loading && !list.error && !hits.length && <div className="text-gray-400 p-2">no panel matches</div>}
          </div>
          <label className="mt-1.5 block border-t border-gray-100 pt-1.5 text-gray-600 cursor-pointer hover:text-accent-700">
            or load a file: a panel downloaded from PanelApp (TSV) or a list of gene symbols
            <input type="file" accept=".tsv,.txt,.csv" className="hidden" onChange={async e => { const f = e.target.files?.[0]; if (f) { onPanel(parsePanelFile(f.name, await f.text())); setOpen(false); } e.target.value = ''; }} />
          </label>
        </div>
      )}
    </span>
  );
}
