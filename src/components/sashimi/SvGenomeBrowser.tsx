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
import { Icon, Pill, Segmented } from './controls';
import { ACHRO_REGIONS, ACHRO_SOURCE, ACHRO_GUIDELINES, ACHRO_COLOR, ACHRO_CLASS_LABEL, achroShort, type AchroRegion, type AchroRow } from './achropuce';
import { niceStep } from './geometry';
import { clinvarUrl, type ClinvarVariant } from './clinvar';
import { predictNmd, NMD_BOUNDARY_NT, type NmdCall } from './nmd';
import { PANEL_SOURCES, listPanels, loadPanel, parsePanelFile, geneModes, classifyMoi, cleanDisease, type PanelSummary, type LoadedPanel, type GeneModes, type MoiClass } from '../../standalone/panelapp';

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
/** ClinVar (pathogenic / likely pathogenic): SVs drawn up to this span in view (as DGV), SNVs and short indels up to the second */
const CLIN_SV_MAX_SPAN = 30e6, CLIN_SNV_MAX_SPAN = 10e6, CLIN_SV_LANE_H = 9, CLIN_SV_MAX_LANES = 60, CLIN_SNV_H = 48;
/** the ClinVar SV summary: a profile of the SVs over each point (losses up, gains down), and the bars matching the event under it */
const CLIN_AGG_H = 30, CLIN_AGG_BIN = 2;
/** the ClinVar SVs are summed up when more than this many lie over one point (their bars would need more lanes) */
const CLIN_SV_STACK_MAX = 8;
/** the ClinVar SNVs are summed up (a mirrored histogram) past this span; their tiles stay at most 4 Mb whatever the zoom, as
 * UCSC sends every class and the page keeps the pathogenic ones (about 1-2 kB a variant) */
const CLIN_SNV_SUM_SPAN = 1e6, CLIN_SNV_TILE_MAX = 4e6, CLIN_SNV_BIN = 4;
const CLIN_KIND_COLOR: Record<string, string> = { loss: '#dc2626', gain: '#2563eb', inv: '#7c3aed', complex: '#7c3aed', other: '#64748b', snv: '#b91c1c', indel: '#b91c1c' };
/** SNVs: pathogenic dark red, likely pathogenic orange */
const CLIN_SIGN_COLOR = { PG: '#b91c1c', LP: '#f97316' } as const;
/** SNVs by what they do: loss of function rises above the line, missense hangs below it, in frame and other sit on it */
const CLIN_EFFECT_COLOR = { lof: '#be123c', missense: '#6d28d9', inframe: '#64748b', other: '#94a3b8' } as const;
/** loss of function predicted to escape NMD (the loss of function predicted to trigger it keeps CLIN_EFFECT_COLOR.lof) */
const CLIN_LOF_ESC = '#ea580c';
/** the ClinVar SNV track summed per gene: its height, and the glyphs' sizes (radius / half-width by count) */
const CLIN_SNV_GENE_H = 70;
const snvR = (n: number) => Math.min(13, 6.6 + 1.5 * Math.sqrt(n));
/** the missense ▼: half-width by count, flat enough (height 1.3 × half-width) for its number to sit in its wide part */
const snvTri = (n: number) => Math.min(15, 10 + 1.3 * Math.sqrt(n));
const CLIN_EFFECT_WORD = { lof: 'loss of function', missense: 'missense', inframe: 'in frame', other: 'other' } as const;
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
  /** a stretch moved between builds (1-based inclusive), for coordinates published in the other build */
  liftOver?: (from: GenomeBuild, to: GenomeBuild, chrom: string, start: number, end: number) => Promise<{ chrom: string; start: number; end: number } | null>;
  /**
   * a transcript HGVS (NM_…:c.…) on the genome: its genomic HGVS on the NC_ accession, its SPDI or VCF form when given,
   * the service's warnings and the service that answered; `onStep` names each service as it is asked
   */
  recodeHgvs?: (hgvs: string, onStep?: (service: string) => void) => Promise<{ hgvsg: string; spdi: string | null; vcf: string | null; warnings: string[]; source: string }>;
  /** ClinVar variants classified pathogenic / likely pathogenic over a 0-based half-open window: SVs (≥ 50 bp) or SNVs and short indels */
  getClinvar?: (which: 'sv' | 'snv', chrom: string, start: number, end: number) => Promise<ClinvarVariant[]>;
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
  /** an interval to highlight (as a Ctrl / ⌘ + drag would): each new `key` sets it */
  highlight?: { start: number; end: number; key: number } | null;
  /** a small variant marked on the reference (no rebuilt chromosome): its stretch (0-based half-open), a label, its card's lines */
  variant?: { start: number; end: number; key: string; label: string; lines: string[] };
}

/** the small variant's mark */
const VARIANT_COLOR = '#0f172a';
/** Where a position lies in a transcript: exon n of N (UTR or coding), intron n, upstream or downstream. */
function whereInTx(tx: SvTranscript, pos: number): string {
  const ex = [...tx.exons].sort((a, b) => a.rank - b.rank), n = ex.length;
  const e = ex.find(x => pos >= x.start && pos < x.end);
  if (e) {
    const utr = tx.cdsStart == null || tx.cdsEnd == null ? 'non-coding' : pos < tx.cdsStart ? (tx.strand > 0 ? '5′ UTR' : '3′ UTR') : pos >= tx.cdsEnd ? (tx.strand > 0 ? '3′ UTR' : '5′ UTR') : '';
    return `exon ${e.rank} of ${n}${utr ? ` (${utr})` : ''}`;
  }
  for (let i = 0; i + 1 < n; i++) {
    const lo = Math.min(ex[i].end, ex[i + 1].end), hi = Math.max(ex[i].start, ex[i + 1].start);
    if (pos >= lo && pos < hi) return `intron ${ex[i].rank} of ${n - 1}`;
  }
  const before = tx.strand > 0 ? pos < tx.start : pos >= tx.end;
  return before ? 'upstream' : 'downstream';
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
function pack<T extends { x0: number; x1: number; xg: number; h: number; worth: boolean; lw: number }>(items: T[], maxLanes = MAX_LANES, prev?: Map<string, number>, key?: (it: T) => string, moving = false): { items: (T & { lane: number; labeled: boolean })[]; laneY: number[]; height: number } {
  // `prev`: the lanes of the last drawing. The items drawn then are placed first, each in its lane while it still fits
  // there, the newcomers in the gaps left: zooming and panning do not shuffle the labels (the map is then updated).
  // `moving` (a zoom or a pan under way): an item keeps its lane, its label hidden while it does not fit, as maps do;
  // everything is laid out again once the view rests.
  // the lanes of the last drawing hold only while the view moves; at rest everything is packed afresh (once)
  const pref = (it: T) => (moving && prev && key ? prev.get(key(it)) : undefined);
  const order = items.map((it, i) => ({ it, i, p: pref(it) }))
    .sort((a, b) => Number(a.p == null) - Number(b.p == null) || Number(b.it.worth) - Number(a.it.worth) || a.it.x0 - b.it.x0);
  const place = (cap: number, compact: boolean) => {
    const spans: [number, number][][] = [], laneH: number[] = [], lane = new Array<number>(items.length), labeled = new Uint8Array(items.length);
    const fits = (l: number, a: number, b: number) => spans[l].every(([x, y]) => y + 8 < a || x - 8 > b);
    const put = (a: number, b: number, p?: number): number => {
      // its lane of the last drawing, created when the items placed so far did not need it yet
      if (p != null && p < cap) { while (spans.length <= p) { spans.push([]); laneH.push(0); } if (fits(p, a, b)) return p; }
      for (let l = 0; l < spans.length; l++) if (fits(l, a, b)) return l;
      if (spans.length < cap) { spans.push([]); laneH.push(0); return spans.length - 1; }
      return -1;
    };
    for (const { it, i, p } of order) {
      const withLabel = Math.max(it.xg, it.x1, it.x0 + it.lw);
      if (moving && p != null && p < cap) {
        while (spans.length <= p) { spans.push([]); laneH.push(0); }
        // its lane whatever (glyphs may touch for a moment), its label only where it fits
        const full = fits(p, it.x0, withLabel);
        if (full) labeled[i] = 1;
        spans[p].push([it.x0, full ? withLabel : it.xg]); lane[i] = p; laneH[p] = Math.max(laneH[p], it.h); continue;
      }
      let l = put(it.x0, withLabel, p);
      if (l >= 0) { labeled[i] = 1; spans[l].push([it.x0, withLabel]); }
      else if (!compact) return null;
      else {
        l = put(it.x0, it.xg, p);
        if (l < 0) l = spans.reduce((best, sp, k) => (sp.length < spans[best].length ? k : best), 0);
        spans[l].push([it.x0, it.xg]);
      }
      lane[i] = l; laneH[l] = Math.max(laneH[l], it.h);
    }
    return { lane, laneH, labeled };
  };
  // every label in its own lane if they fit in maxLanes lanes; else the glyphs packed and the labels drawn where they fit
  const r = place(maxLanes, false) ?? place(maxLanes, true)!;
  if (prev && key) { prev.clear(); items.forEach((it, i) => prev.set(key(it), r.lane[i])); }
  const laneY: number[] = []; let y = 0;
  for (const h of r.laneH) { laneY.push(y); y += h + LANE_GAP; }
  return { items: items.map((it, i) => ({ ...it, lane: r.lane[i], labeled: !!r.labeled[i] })), laneY, height: Math.max(26, y) };
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



/** Modes of inheritance (PanelApp), as colours and words. */
const MOI: Record<MoiClass | 'none', { color: string; short: string; label: string }> = {
  AD: { color: '#dc2626', short: 'AD', label: 'monoallelic (AD)' },
  'AD+AR': { color: '#7c3aed', short: 'AD/AR', label: 'mono- or biallelic' },
  AR: { color: '#2563eb', short: 'AR', label: 'biallelic (AR)' },
  XLR: { color: '#db2777', short: 'XLR', label: 'X-linked, females biallelic' },
  XLD: { color: '#9d174d', short: 'XLD', label: 'X-linked, females may be affected' },
  MT: { color: '#059669', short: 'MT', label: 'mitochondrial' },
  '?': { color: '#94a3b8', short: '?', label: 'unknown / other' },
  none: { color: '#cbd5e1', short: '—', label: 'not in PanelApp' },
};


/** Which genes are drawn: every criterion set must hold (empty or 'any' = not filtering on it). */
interface GeneFilter {
  panel: 'any' | 'in' | 'green';
  moi: (MoiClass | 'none')[];
  /** LOEUF strictly below this, or null */
  loeuf: number | null;
  changed: 'any' | 'changed' | 'near';
}
const NO_FILTER: GeneFilter = { panel: 'any', moi: [], loeuf: null, changed: 'any' };
const filterCount = (f: GeneFilter, hasPanel: boolean, hasEvent: boolean) => (hasPanel && f.panel !== 'any' ? 1 : 0) + (f.moi.length ? 1 : 0) + (f.loeuf != null ? 1 : 0) + (hasEvent && f.changed !== 'any' ? 1 : 0);

const EMPTY_STORY: GeneStory = { ders: [], genes: [], chains: [], fragments: [], domains: [], derDomains: [] };
type Vp = { s: number; e: number };
type TipAt = { x: number; y: number; gene: string; dgv?: PopulationSv; achro?: AchroPlaced; achroType?: 'del' | 'dup'; clin?: ClinvarVariant[];
  /** the ClinVar SVs over one point of the summary: that point (0-based) */
  clinAt?: number;
  /** the ClinVar SNVs of one bar of the summary: its stretch (0-based half-open) */
  clinRange?: [number, number];
  /** the ClinVar SNVs of one gene (or of neighbouring genes merged): their names, the transcripts NMD is predicted on, their span */
  clinGenes?: { names: string[]; tx: string[]; s: number; e: number };
  /** the small variant searched */
  variant?: boolean };
/** AChro-Puce bars: deletions red, duplications blue */
const ACHRO_TYPE_COLOR = { del: '#dc2626', dup: '#2563eb' } as const;
/** an AChro-Puce region placed on the chromosome in view: 0-based half-open, and how it was placed */
type AchroSpan = { start: number; end: number; approx: boolean; how: string };
/** a region placed: the union of its types' stretches, and each type's own (a deletion and a duplication may lie apart) */
type AchroPlaced = AchroSpan & { r: AchroRegion; types: Partial<Record<'del' | 'dup', AchroSpan>> };
/** the chromosome a region's name starts with (1q21.1… → chr1) */
const achroChrom = (r: AchroRegion) => `chr${r.name.match(/^(\d+|X|Y)[pq]/)?.[1] ?? ''}`;

// ======================== Component ========================

export default function SvGenomeBrowser({ rec: eventRec, chrom, window: eventWin, segColors, flankCopies, zygosity, build, loaders, rebuild, focus, highlight, variant }: Props) {
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
  /** how far the view may go: the reference, and past its ends as far as a longer derivative reaches (set below) */
  const boundsRef = useRef({ lo: 0, hi: L });
  const clampVp = (s: number, e: number): Vp => {
    const { lo, hi } = boundsRef.current;
    const span = Math.min(hi - lo, Math.max(MIN_SPAN, e - s));
    const c = Math.min(hi - span / 2, Math.max(lo + span / 2, (s + e) / 2));
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
  /** the gene joined to its copies: set by a click on it, kept when its card is closed, cleared by a click on the background */
  const [linked, setLinked] = useState<string | null>(null);
  const tipEl = useRef<HTMLDivElement | null>(null);
  const [colorBy, setColorBy] = useState<'loeuf' | 'impact' | 'moi'>(() => store.get<'loeuf' | 'impact' | 'moi'>('sashimi.svColorBy') ?? 'loeuf');
  const [near, setNear] = useState<number>(() => store.get<number>('sashimi.svNear') ?? NEAR_BREAKPOINT_DEFAULT);
  /** the population CNVs (DGV Gold Standard): drawn under the reference genes, and compared with the event's losses and gains */
  const [showDgv, setShowDgv] = useState<boolean>(() => store.get<boolean>('sashimi.svDgv') ?? true);
  const [pop, setPop] = useState<{ list: PopulationSv[]; loading: boolean; error?: string }>({ list: [], loading: false });
  /** the AChro-Puce recurrent CNVs (their sheets' classes), drawn under the population CNVs */
  const [showAchro, setShowAchro] = useState<boolean>(() => store.get<boolean>('sashimi.svAchro') ?? true);
  const [showClinSv, setShowClinSv] = useState<boolean>(() => store.get<boolean>('sashimi.svClinvarSv') ?? true);
  const [showClinSnv, setShowClinSnv] = useState<boolean>(() => store.get<boolean>('sashimi.svClinvarSnv') ?? true);
  /** the ClinVar SVs one by one (expanded), or summed up as a profile (the default) */
  // bars unless more than CLIN_SV_STACK_MAX are stacked over a point (then summed up), or every SV as bars when asked
  const [clinSvForce, setClinSvForce] = useState<'auto' | 'all'>(() => (store.get<string>('sashimi.svClinvarSvMode') === 'all' ? 'all' : 'auto'));
  /** ClinVar SNVs summed per gene (by default) or each variant on its own */
  const [clinSnvMode, setClinSnvMode] = useState<'gene' | 'each'>(() => (store.get<string>('sashimi.svClinvarSnvMode') === 'each' ? 'each' : 'gene'));
  /** the genes drawn: in the panel, of some modes of inheritance, under a LOEUF, changed by the SV (kept) */
  const [geneFilter, setGeneFilterState] = useState<GeneFilter>(() => { const f = { ...NO_FILTER, ...(store.get<GeneFilter>('sashimi.svGeneFilter') ?? {}) }; return { ...f, moi: (f.moi ?? []).filter(m => m in MOI) }; });
  const setGeneFilter = (f: GeneFilter) => { setGeneFilterState(f); store.set('sashimi.svGeneFilter', f); };
  /** coordinates published in the other build, moved to this one (Ensembl), by region id; null when they could not be */
  const [achroLifted, setAchroLifted] = useState<Map<string, { start: number; end: number } | null>>(new Map());
  const [bands, setBands] = useState<Cytoband[]>([]);
  const svgRef = useRef<SVGSVGElement>(null);
  const boxRef = useRef<HTMLDivElement>(null);

  // ======================== the view: a stretch of the reference, moved once per frame ========================
  const initialVp = (): Vp => { const pad = Math.max(2000, (eventWin.end - eventWin.start) * (eventRec ? 0.12 : 0)); return clampVp(eventWin.start - pad, eventWin.end + pad); };
  const [vp, setVpState] = useState<Vp>(initialVp);
  const vpRef = useRef(vp);
  const raf = useRef(0);
  const anim = useRef<{ fly?: { from: Vp; to: Vp; t0: number; dur: number; a0: number | null; a1: number | null }; fling?: { v: number; t: number } }>({});
  /**
   * Where along the reference the rebuilt chromosomes are lined up (their offset is read there): set by a zoom or a
   * jump, kept by a pan, so that dragging moves the reference and the rebuilt chromosomes together, the figure whole.
   * null: the view's centre.
   */
  const anchor = useRef<number | null>(null);
  /** the rebuilt chromosomes' offset from where they should sit, decaying to 0 (a slide instead of a jump) */

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
      // the line-up glides with the flight, from where it was to the target's
      const a0 = a.fly.a0 ?? c0, a1 = a.fly.a1 ?? c1;
      anchor.current = a0 + (a1 - a0) * k2;
      if (t < 1) more = true; else a.fly = undefined;
    } else if (a.fling) {
      const dt = Math.min(48, now - a.fling.t); a.fling.t = now;
      const sh = a.fling.v * dt;
      next = clampVp(next.s + sh, next.e + sh);
      a.fling.v *= Math.pow(0.9, dt / 16);
      if (Math.abs(a.fling.v) * 16 > (next.e - next.s) / PW * 0.25 && next.s > boundsRef.current.lo && next.e < boundsRef.current.hi) more = true; else a.fling = undefined;
    }
    vpRef.current = next;
    lastMove.current = now;
    setVpState(next);
    if (more) raf.current = requestAnimationFrame(frame);
    window.clearTimeout(settleTimer.current);
    settleTimer.current = window.setTimeout(settle, 450);
  };
  const kick = () => { if (!raf.current) raf.current = requestAnimationFrame(frame); };
  /**
   * A move by the user: animations stop, the view follows. `at`: where the rebuilt chromosomes line up afterwards
   * (a zoom: the point it zooms about); left out, the line-up stays (a pan moves the figure whole).
   */
  const moveTo = (s: number, e: number, at?: number) => {
    anim.current = {};
    const v = vpRef.current;
    // a pan from the view's centre line-up: that line-up is fixed first (else it would follow the centre, cancelling the pan)
    anchor.current = at !== undefined ? at : anchor.current ?? (v.s + v.e) / 2;
    vpRef.current = clampVp(s, e); kick();
  };
  /** a flight there; `at`: the line-up at the end ('keep': as it is, a pan), the target's centre by default */
  const flyTo = (s: number, e: number, dur = 480, at: number | 'keep' = (s + e) / 2) => {
    const from = vpRef.current, a0 = anchor.current ?? (from.s + from.e) / 2;
    anim.current = { fly: { from, to: clampVp(s, e), t0: performance.now(), dur, a0, a1: at === 'keep' ? a0 : at } }; kick();
  };
  useEffect(() => () => { cancelAnimationFrame(raf.current); window.clearTimeout(settleTimer.current); }, []);
  // a new chromosome or event: start over on it
  useEffect(() => { const v = eventRec ? clampVp(eventRange().s, eventRange().e) : initialVp(); anim.current = {}; anchor.current = null; vpRef.current = v; setVpState(v); }, [chrom, build, eventWin.start, eventWin.end]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (focus) flyTo(focus.start, focus.end, 700); }, [focus?.key]); // eslint-disable-line react-hooks/exhaustive-deps
  const span = vp.e - vp.s;
  /** px per base, shared by the reference and the rebuilt chromosomes (set once their stretches are known, below) */
  let ppb = PW / span;
  /** bases per px as drawn, for the gestures */
  const bppRef = useRef(span / PW);
  bppRef.current = span / PW;
  const refX = (pos: number) => X0 + (pos - vp.s) * ppb;
  // the level of detail changes 15 % past its thresholds, so a zoom hovering at one does not flip it back and forth
  const lodRef = useRef<'full' | 'bars' | 'density'>(span <= LOD_FULL ? 'full' : span <= LOD_BARS ? 'bars' : 'density');
  {
    const up = (t: number) => span > t * 1.15, down = (t: number) => span < t / 1.15, cur = lodRef.current;
    if (cur === 'full' && up(LOD_FULL)) lodRef.current = up(LOD_BARS) ? 'density' : 'bars';
    else if (cur === 'bars') { if (down(LOD_FULL)) lodRef.current = 'full'; else if (up(LOD_BARS)) lodRef.current = 'density'; }
    else if (cur === 'density' && down(LOD_BARS)) lodRef.current = down(LOD_FULL) ? 'full' : 'bars';
  }
  const lod = lodRef.current;
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
  // ---- modes of inheritance (PanelApp, every panel): asked for the event's changed genes, the genes in view in the
  //      inheritance colouring, and the gene of an open card; four requests at a time, cached ----
  const moiSrc = PANEL_SOURCES.find(s => s.label === panelData?.source) ?? PANEL_SOURCES[0];
  const modes = useRef(new Map<string, GeneModes | 'loading' | 'error'>());
  const [modesVersion, bumpModes] = useReducer((x: number) => x + 1, 0);
  const modeQueue = useRef<string[]>([]), modeActive = useRef(0), modeBump = useRef(0);
  useEffect(() => { modes.current = new Map(); modeQueue.current = []; bumpModes(); }, [moiSrc.id]);
  const pumpModes = () => {
    while (modeActive.current < 4 && modeQueue.current.length) {
      const sym = modeQueue.current.shift()!, src = moiSrc;
      modeActive.current++;
      geneModes(src, sym).then(m => { if (src === moiSrc) modes.current.set(sym, m); }, () => { modes.current.set(sym, 'error'); })
        .finally(() => { modeActive.current--; pumpModes(); window.clearTimeout(modeBump.current); modeBump.current = window.setTimeout(bumpModes, 120); });
    }
  };
  const wantModes = (syms: string[]) => { let any = false; for (const s2 of syms) if (!modes.current.has(s2)) { modes.current.set(s2, 'loading'); modeQueue.current.push(s2); any = true; } if (any) pumpModes(); };
  /** a gene's mode: the loaded panel's own value first, else the commonest among its PanelApp entries (green first) */
  /** a gene's diseases in PanelApp: the loaded panel's first, then those of its panels, green ones first (names cleaned of identifiers) */
  const diseasesOf = (sym: string): { name: string; panels: number; green: number }[] => {
    const out = new Map<string, { name: string; panels: number; green: number }>();
    for (const ph of panel?.genes.get(sym.toUpperCase())?.phenotypes ?? []) { const n = cleanDisease(ph); if (n && !out.has(n.toLowerCase())) out.set(n.toLowerCase(), { name: n, panels: 0, green: 0 }); }
    const m = modes.current.get(sym);
    if (m && m !== 'loading' && m !== 'error') for (const d of m.diseases) if (!out.has(d.name.toLowerCase())) out.set(d.name.toLowerCase(), d);
    return [...out.values()];
  };
  const modeOf = (sym: string): { cls: MoiClass | 'none'; m: GeneModes | null; fromPanel: string | null } => {
    const pg = panel?.genes.get(sym.toUpperCase()), m = modes.current.get(sym);
    const mm = m && typeof m === 'object' ? m : null;
    if (pg?.moi) return { cls: classifyMoi(pg.moi), m: mm, fromPanel: panel!.name };
    if (!mm) return { cls: m === 'loading' || m === undefined ? '?' : 'none', m: null, fromPanel: null };
    return { cls: mm.entries.length ? mm.main : 'none', m: mm, fromPanel: null };
  };
  void modesVersion;
  /** the mode's badge after a gene name: its width (0 when none) and the badge itself */
  const badgeText = (sym: string) => { const c = modeOf(sym).cls; return c === 'none' || c === '?' ? null : MOI[c].short; };
  const badgeW = (sym: string) => { const t = badgeText(sym); return t ? t.length * 5.6 + 9 : 0; };
  const badge = (sym: string, x: number, y: number) => {
    const t = badgeText(sym); if (!t) return null;
    const m = modeOf(sym), col = MOI[m.cls].color, w = t.length * 5.6 + 6;
    return <g key="moi"><rect x={x} y={y - 8.5} width={w} height={11} rx={5.5} fill={col} fillOpacity={0.13} stroke={col} strokeOpacity={0.55} strokeWidth={0.7} />
      <text x={x + w / 2} y={y} fontSize={8.5} fontWeight={700} textAnchor="middle" fill={col} fontStyle="normal">{t}</text>
      <title>{`Mode of inheritance (PanelApp): ${MOI[m.cls].label}${m.fromPanel ? ` in ${m.fromPanel}` : m.m ? ` · ${m.m.parts.map(p => `${MOI[p.cls].short} in ${p.panels} panel${p.panels > 1 ? 's' : ''}`).join(', ')}` : ''}${(() => { const ds = diseasesOf(sym); return ds.length ? `\nDisease: ${ds.slice(0, 5).map(d => d.name).join('; ')}${ds.length > 5 ? `; +${ds.length - 5} more` : ''}` : ''; })()}`}</title></g>;
  };
  useEffect(() => { wantModes(story.genes.filter(g => g.kind !== 'unaffected').map(g => g.gene.symbol)); }, [story, moiSrc.id]); // eslint-disable-line react-hooks/exhaustive-deps
  /** does a gene pass the gene filter (a criterion whose data is still on its way does not hold yet) */
  const passes = (sym: string): boolean => {
    const f = geneFilter;
    if (panel && f.panel !== 'any') { const pg = panel.genes.get(sym.toUpperCase()); if (!pg || (f.panel === 'green' && pg.confidence !== 3)) return false; }
    if (f.loeuf != null) { const c = ann.constraint.get(sym); if (c?.loeuf == null || c.loeuf >= f.loeuf) return false; }
    if (rec && f.changed !== 'any') { const k = impactOf.get(sym)?.kind; if (!k || k === 'unaffected' || (f.changed === 'changed' && k === 'nearby')) return false; }
    if (f.moi.length) {
      const m = modes.current.get(sym), pg = panel?.genes.get(sym.toUpperCase());
      if (!pg?.moi && (m === undefined || m === 'loading')) return false;
      const cls = modeOf(sym).cls;
      if (!f.moi.includes(cls)) return false;
    }
    return true;
  };
  const filtering = filterCount(geneFilter, !!panel, !!rec) > 0;
  const geneColor = (symbol: string, gi: GeneImpact | undefined) => (colorBy === 'loeuf' ? LOEUF_CLASS[loeufClass(ann.constraint.get(symbol))].color : colorBy === 'moi' ? MOI[modeOf(symbol).cls].color : IMPACT[gi?.kind ?? 'unaffected'].color);
  const panelDot = (sym: string) => { const pg = panel?.genes.get(sym.toUpperCase()); return pg ? CONF[pg.confidence].color : null; };
  /** worth a label in a crowded view: changed by the event, or in the panel, dosage-sensitive or LoF-intolerant */
  const worthSym = (sym: string) => { const gi = impactOf.get(sym); return (!!gi && gi.kind !== 'unaffected' && gi.kind !== 'nearby') || notableSym(sym, annotations); };

  // ======================== the rebuilt chromosomes: on the reference's axis ========================
  // One scale for all, and one axis: outside the event each rebuilt chromosome lies exactly under the reference it
  // carries (the flank before the event at its own offset, the flank after at its own, the difference being the event's
  // gain or loss); across the event the offset goes from one to the other in proportion to where the view's middle is.
  // So zooming never moves anything, panning slides the rebuilt chromosome smoothly across the event, and its length
  // shows: a longer one runs past the reference's end (or starts before its start), a shorter one stops before.
  const ders = story.ders.filter(d => d.rearranged);
  /** each rebuilt chromosome's whole extent in its own positions (its ends are the reference's), and its length */
  const derExtent = (d: DerLayout) => { let lo = Infinity, hi = -Infinity; for (const dp of d.pieces) { const [a, b] = pieceRange(dp, L); for (const v of [toDer(dp, a), toDer(dp, b)]) { lo = Math.min(lo, v); hi = Math.max(hi, v); } } return { lo, hi, len: hi - lo }; };
  const bpList = bpKey ? bpKey.split(',').map(Number) : [];
  const evL = bpList.length ? bpList[0] : 0, evR = bpList.length ? bpList[bpList.length - 1] : 0;
  /** the offset (derivative position − reference position) of each flank: the pieces holding the chromosome's two ends */
  const flankOff = useMemo(() => new Map(ders.map(d => {
    const at = (pos: number) => { for (const dp of d.pieces) { const [lo, hi] = pieceRange(dp, L); if (pos >= lo && pos <= hi && !dp.piece.reversed) return toDer(dp, pos) - pos; } return null; };
    const ext = derExtent(d);
    return [d.name, { l: at(Math.max(0, evL - 1)) ?? ext.lo, r: at(Math.min(L, evR + 1)) ?? ext.hi - L }];
  })), [story, L, evL, evR]); // eslint-disable-line react-hooks/exhaustive-deps
  const offsetAt = (d: DerLayout, c: number) => {
    const f = flankOff.get(d.name); if (!f) return 0;
    if (c <= evL || evR <= evL) return f.l;
    if (c >= evR) return f.r;
    return f.l + (f.r - f.l) * (c - evL) / (evR - evL);
  };
  const c = anchor.current ?? (vp.s + vp.e) / 2;
  const derOff = new Map(ders.map(d => [d, offsetAt(d, c)]));
  const derX = (d: DerLayout, dpos: number) => X0 + (dpos - (derOff.get(d) ?? 0) - vp.s) * ppb;
  /** the reference drawn: the view (it may reach past the reference's ends, where a longer derivative goes on) */
  const ve = Math.min(L, vp.e);
  // how far past the reference's ends the view may go: as far as the longest derivative reaches
  const extra = Math.max(0, ...ders.map(d => derExtent(d).len - L));
  boundsRef.current = { lo: -extra, hi: L + extra };
  /**
   * The event in full: its stretch of the reference, and that of each rebuilt chromosome (longer when it gained), as
   * they stand when the view is centred on the event.
   */
  const eventRange = (): Vp => {
    const pad = Math.max(2000, (eventWin.end - eventWin.start) * 0.12);
    let lo = eventWin.start, hi = eventWin.end;
    for (let it = 0; it < 2; it++) {
      const cc = (lo + hi) / 2;
      for (const d of ders) { const f = flankOff.get(d.name); if (!f) continue; const off = offsetAt(d, cc); lo = Math.min(lo, evL + f.l - off); hi = Math.max(hi, evR + f.r - off); }
    }
    return { s: lo - pad, e: hi + pad };
  };
  /** the whole chromosome, with its derivatives' ends, from where the view would stand */
  const wholeView = (): Vp => {
    let lo = 0, hi = L;
    for (let it = 0; it < 2; it++) {
      const cc = (lo + hi) / 2; let a = 0, b = L;
      for (const d of ders) { const e = derExtent(d), off = offsetAt(d, cc); a = Math.min(a, e.lo - off); b = Math.max(b, e.hi - off); }
      lo = a; hi = b;
    }
    const pad = (hi - lo) * 0.01;
    return { s: lo - pad, e: hi + pad };
  };

  /** the stretches of the reference the rebuilt chromosomes show (their genes are fetched and drawn too) */
  let needS = vp.s, needE = ve;
  for (const d of ders) {
    const off = derOff.get(d) ?? 0, dl = vp.s + off, dr = vp.e + off;
    for (const dp of d.pieces) {
      const [lo, hi] = pieceRange(dp, L);
      const a = Math.min(toDer(dp, lo), toDer(dp, hi)), b = Math.max(toDer(dp, lo), toDer(dp, hi));
      const oa = Math.max(a, dl), ob = Math.min(b, dr);
      if (ob <= oa) continue;
      const back = (dpos: number) => (dp.piece.reversed ? dp.seg.end - (dpos - dp.derStart) : dp.seg.start + (dpos - dp.derStart));
      needS = Math.min(needS, back(oa), back(ob)); needE = Math.max(needE, back(oa), back(ob));
    }
  }
  needS = Math.max(0, needS); needE = Math.min(L, Math.max(needE, needS + 1));

  // ======================== wherever the view goes: genes and population CNVs by tiles ========================
  const [tilesVersion, bumpTiles] = useReducer((x: number) => x + 1, 0);
  const geneTiles = useMemo(() => (loaders.getTranscripts ? new TileCache<SvTranscript>((a, b, size) => loaders.getTranscripts!(chrom, a, b).then(list => (size >= 16e6 ? reduceTile(list) : list)), L, 36, bumpTiles) : null), [loaders, chrom, L]);
  const dgvTiles = useMemo(() => (loaders.getPopulationSvs ? new TileCache<PopulationSv>((a, b) => loaders.getPopulationSvs!(chrom, a, b), L, 24, bumpTiles) : null), [loaders, chrom, L]);
  const clinSvTiles = useMemo(() => (loaders.getClinvar ? new TileCache<ClinvarVariant>((a, b) => loaders.getClinvar!('sv', chrom, a, b), L, 24, bumpTiles) : null), [loaders, chrom, L]);
  const clinSnvTiles = useMemo(() => (loaders.getClinvar ? new TileCache<ClinvarVariant>((a, b) => loaders.getClinvar!('snv', chrom, a, b), L, 24, bumpTiles) : null), [loaders, chrom, L]);
  const tileSize = tileSizeFor(span);
  const firstFetch = useRef(true);
  // the tiles are asked for once the view rests a moment: a fast zoom does not fetch every level it crosses
  useEffect(() => {
    const go = () => {
      geneTiles?.want(needS, needE, tileSize); if (showDgv && span <= DGV_MAX_SPAN) dgvTiles?.want(vp.s, ve, tileSize);
      if (showClinSv && span <= CLIN_SV_MAX_SPAN) clinSvTiles?.want(vp.s, ve, tileSize);
      if (showClinSnv && span <= CLIN_SNV_MAX_SPAN) clinSnvTiles?.want(vp.s, ve, Math.min(tileSize, CLIN_SNV_TILE_MAX));
    };
    if (firstFetch.current) { firstFetch.current = false; go(); return; }
    const id = window.setTimeout(go, 140);
    return () => window.clearTimeout(id);
  }, [vp.s, vp.e, needS, needE, tileSize, geneTiles, dgvTiles, showDgv, clinSvTiles, clinSnvTiles, showClinSv, showClinSnv]); // eslint-disable-line react-hooks/exhaustive-deps
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
    const got = dgvTiles.collect((Math.floor(vp.s / tileSize) - 1) * tileSize, (Math.floor(ve / tileSize) + 2) * tileSize, tileSize);
    const seen = new Set<string>(), list: PopulationSv[] = [];
    for (const v of got.items) { const k = `${v.id}:${v.outer[0]}`; if (!seen.has(k)) { seen.add(k); list.push(v); } }
    return { list, loading: got.loading, error: got.error };
  }, [tilesVersion, tileSize, Math.floor(vp.s / tileSize), Math.floor(ve / tileSize), dgvTiles, showDgv, span > DGV_MAX_SPAN]); // eslint-disable-line react-hooks/exhaustive-deps
  /** the ClinVar variants of the view (tiles of several sizes may hold one twice) */
  const clinView = (tiles: TileCache<ClinvarVariant> | null, on: boolean, max: number, ts = tileSize, lo = vp.s, hi = ve) => {
    if (!tiles || !on || span > max) return { list: [] as ClinvarVariant[], loading: false, error: null as string | null };
    const got = tiles.collect((Math.floor(Math.min(lo, vp.s) / ts) - 1) * ts, (Math.floor(Math.max(hi, ve) / ts) + 2) * ts, ts);
    const seen = new Set<string>(), list: ClinvarVariant[] = [];
    for (const v of got.items) { const k = `${v.id}:${v.start}`; if (!seen.has(k)) { seen.add(k); list.push(v); } }
    return { list, loading: got.loading, error: got.error };
  };
  const clinSvView = useMemo(() => clinView(clinSvTiles, showClinSv, CLIN_SV_MAX_SPAN), [tilesVersion, tileSize, Math.floor(vp.s / tileSize), Math.floor(ve / tileSize), clinSvTiles, showClinSv, span > CLIN_SV_MAX_SPAN]); // eslint-disable-line react-hooks/exhaustive-deps
  const snvTs = Math.min(tileSize, CLIN_SNV_TILE_MAX);
  // summed per gene, a gene counts all its variants: the SNVs load over the whole of the genes in view (within bounds)
  const snvExt = (() => {
    if (clinSnvMode !== 'gene') return [vp.s, ve] as const;
    let lo = vp.s, hi = ve;
    for (const g of view.genes) if (g.end > vp.s && g.start < ve) { lo = Math.min(lo, g.start); hi = Math.max(hi, g.end); }
    const room = Math.max(0, CLIN_SNV_MAX_SPAN - span) / 2;
    return [Math.max(0, lo, vp.s - room), Math.min(L, hi, ve + room)] as const;
  })();
  useEffect(() => {
    if (clinSnvMode === 'gene' && showClinSnv && span <= CLIN_SNV_MAX_SPAN && (snvExt[0] < vp.s || snvExt[1] > ve)) clinSnvTiles?.want(snvExt[0], snvExt[1], snvTs);
  }, [Math.floor(snvExt[0] / snvTs), Math.floor(snvExt[1] / snvTs), snvTs, clinSnvMode, showClinSnv, clinSnvTiles, span > CLIN_SNV_MAX_SPAN]); // eslint-disable-line react-hooks/exhaustive-deps
  const clinSnvView = useMemo(() => clinView(clinSnvTiles, showClinSnv, CLIN_SNV_MAX_SPAN, snvTs, snvExt[0], snvExt[1]), [tilesVersion, snvTs, Math.floor(vp.s / snvTs), Math.floor(ve / snvTs), Math.floor(snvExt[0] / snvTs), Math.floor(snvExt[1] / snvTs), clinSnvTiles, showClinSnv, span > CLIN_SNV_MAX_SPAN]); // eslint-disable-line react-hooks/exhaustive-deps
  const geneBySym = useMemo(() => new Map(view.genes.map(g => [g.symbol, g])), [view.genes]);
  // a small variant searched: the ClinVar variants around it each at its place (not counted per gene), and the view on
  // its gene once the gene models are in (unless the view was moved meanwhile, or the gene is very long)
  useEffect(() => { if (variant) setClinSnvMode('each'); }, [variant?.key]); // eslint-disable-line react-hooks/exhaustive-deps
  const variantFit = useRef<{ key: string; s: number; done: boolean } | null>(null);
  useEffect(() => {
    if (!variant) return;
    if (variantFit.current?.key !== variant.key) variantFit.current = { key: variant.key, s: vpRef.current.s, done: false };
    const f = variantFit.current;
    if (f.done) return;
    if (Math.abs(vpRef.current.s - f.s) > 1) { f.done = true; return; }
    const at = (variant.start + variant.end) / 2, hits = view.genes.filter(g => g.start <= at && g.end >= at);
    if (!hits.length) { if (!view.loading) f.done = true; return; }
    f.done = true;
    const g = hits.find(x => x.principal.status === 'mane') ?? hits[0];
    if (g.end - g.start > 2.5e6) return;
    const pad = Math.max(2000, (g.end - g.start) * 0.08);
    flyTo(Math.min(g.start, variant.start) - pad, Math.max(g.end, variant.end) + pad, 600);
  }, [view.genes, view.loading, variant?.key]); // eslint-disable-line react-hooks/exhaustive-deps
  // the ClinVar SNVs' genes need their exons at any zoom (the zoomed-out gene tiles keep only each gene's span): tiles of
  // at most CLIN_SNV_TILE_MAX over the SNVs' stretch, the view's own genes when its tiles are that size already
  const snvGenesOwn = tileSize > snvTs && showClinSnv && !!loaders.getClinvar && span <= CLIN_SNV_MAX_SPAN;
  useEffect(() => { if (snvGenesOwn) geneTiles?.want(snvExt[0], snvExt[1], snvTs); },
    [snvGenesOwn, Math.floor(snvExt[0] / snvTs), Math.floor(snvExt[1] / snvTs), snvTs, geneTiles]); // eslint-disable-line react-hooks/exhaustive-deps
  const snvGenes = useMemo(() => {
    if (!snvGenesOwn || !geneTiles) return view.genes;
    const got = geneTiles.collect(Math.floor(snvExt[0] / snvTs) * snvTs, (Math.floor(snvExt[1] / snvTs) + 1) * snvTs, snvTs);
    const byId = new Map<string, SvTranscript>();
    for (const t of got.items) { const o = byId.get(t.id); if (!o || (!o.exons.length && t.exons.length)) byId.set(t.id, t); }
    return groupGenes([...byId.values()]);
  }, [snvGenesOwn, tilesVersion, snvTs, Math.floor(snvExt[0] / snvTs), Math.floor(snvExt[1] / snvTs), geneTiles, view.genes]); // eslint-disable-line react-hooks/exhaustive-deps
  const snvGeneBySym = useMemo(() => new Map(snvGenes.map(g => [g.symbol, g])), [snvGenes]);
  /**
   * Each ClinVar SNV's gene: the one whose principal (MANE Select) coding stretch holds it, ClinVar's gene first among
   * several; else ClinVar's gene; else the gene it lies in. Each loss of function's NMD call on that gene's principal transcript.
   */
  const clinSnvGene = useMemo(() => {
    const geneOf = new Map<ClinvarVariant, SvGene | null>(), nmd = new Map<ClinvarVariant, NmdCall>();
    for (const v of clinSnvView.list) {
      const over = snvGenes.filter(g => g.start <= v.end + 10 && g.end >= v.start - 10);
      const inCds = over.filter(g => { const t = g.principal; return t.cdsStart != null && t.cdsEnd != null && v.start < t.cdsEnd + 10 && v.end > t.cdsStart - 10; });
      const named = (xs: SvGene[]) => xs.find(g => v.genes.includes(g.symbol));
      const g = named(inCds) ?? inCds[0] ?? named(over) ?? v.genes.map(sym => snvGeneBySym.get(sym)).find((x): x is SvGene => !!x) ?? over[0] ?? null;
      geneOf.set(v, g);
      if (v.effect === 'lof') nmd.set(v, predictNmd(v, g?.principal));
    }
    return { geneOf, nmd };
  }, [clinSnvView.list, snvGenes, snvGeneBySym]);
  const nmdOf = (v: ClinvarVariant) => clinSnvGene.nmd.get(v);
  // AChro-Puce regions of this chromosome published only in the other build: moved to this one once
  useEffect(() => {
    setAchroLifted(new Map());
    if (!loaders.liftOver) return;
    let alive = true;
    const b = build as GenomeBuild, other: GenomeBuild = b === 'GRCh38' ? 'GRCh37' : 'GRCh38';
    for (const r of ACHRO_REGIONS) {
      if (achroChrom(r) !== chrom || !r.coords || r.coords.some(c => c.build === b)) continue;
      const c = r.coords.find(x => x.build === other)!;
      loaders.liftOver(other, b, c.chrom, c.start, c.end)
        .then(m => { if (alive) setAchroLifted(prev => new Map(prev).set(r.id, m && m.chrom === chrom ? { start: m.start - 1, end: m.end } : null)); })
        .catch(() => { if (alive) setAchroLifted(prev => new Map(prev).set(r.id, null)); });
    }
    return () => { alive = false; };
  }, [chrom, build, loaders]); // eslint-disable-line react-hooks/exhaustive-deps
  /**
   * The regions of this chromosome, placed: the sheet's coordinates in this build; else those of the other build moved
   * here; else the span of the genes the sheet names, among the gene models loaded (so a region shows once its genes
   * are near the view). Positions given to the Mb, moved, or taken from genes are approximate.
   */
  const achroPlaced: AchroPlaced[] = (() => {
    if (!showAchro) return [];
    const out: AchroPlaced[] = [];
    const genesSpan = (names: string[]) => {
      const found = names.map(g => geneBySym.get(g)).filter((g): g is SvGene => !!g);
      return found.length ? { start: Math.min(...found.map(g => g.start)), end: Math.max(...found.map(g => g.end)), names: found.map(g => g.symbol).join(', ') } : null;
    };
    for (const r of ACHRO_REGIONS) {
      if (achroChrom(r) !== chrom) continue;
      let base: AchroSpan | null = null;
      const own = r.coords?.find(c => c.build === build), lifted = achroLifted.get(r.id), src = r.coords?.[0];
      if (own) base = { start: own.start - 1, end: own.end, approx: !!own.approx, how: `${own.build} coordinates of ${own.from === 'table2' ? ACHRO_GUIDELINES.short : 'the sheet'}${own.note ? ` (${own.note})` : ''}` };
      else if (lifted) base = { start: lifted.start, end: lifted.end, approx: true, how: `${src!.build} coordinates of the sheet${src!.note ? ` (${src!.note})` : ''}, moved to ${build} by Ensembl` };
      else { const g = genesSpan(r.genes); if (g) base = { start: g.start, end: g.end, approx: true, how: `span of the genes the sheet names (${g.names})` }; }
      // a type the sheet places on genes lies on their span (never on the region's coordinates, which are another stretch)
      const types: AchroPlaced['types'] = {};
      for (const t of ['del', 'dup'] as const) {
        if (!r.rows.some(x => x.type === t)) continue;
        const ty = r.typed?.[t];
        if (ty?.genes) { const g = genesSpan(ty.genes); if (g) types[t] = { start: g.start, end: g.end, approx: false, how: `span of ${g.names} (${ty.note})` }; }
        else if (base) types[t] = ty ? { ...base, how: `${base.how}; ${ty.note}` } : base;
      }
      const spans = Object.values(types);
      if (!spans.length) continue;
      out.push({ r, types, start: Math.min(...spans.map(x => x.start)), end: Math.max(...spans.map(x => x.end)), approx: spans.every(x => x.approx), how: base?.how ?? spans[0].how });
    }
    return out;
  })();
  const visGenes = view.genes.filter(g => g.end >= vp.s - span * 0.02 && g.start <= ve + span * 0.02);
  // the modes of the genes in view (their badges; the named ones when a crowd), once the view rests
  useEffect(() => {
    const id = window.setTimeout(() => { const vis = visGenes.length <= 250 ? visGenes : visGenes.filter(g => worthSym(g.symbol)); wantModes(vis.slice(0, 250).map(g => g.symbol)); }, 400);
    return () => window.clearTimeout(id);
  }, [visGenes.length, Math.round(vp.s / Math.max(1, span) * 20), Math.round(Math.log(span) * 10)]); // eslint-disable-line react-hooks/exhaustive-deps
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
  const eventView = () => { const r = eventRange(); flyTo(r.s, r.e); };

  // ======================== layout ========================
  const shown = (g: SvGene) => (mode === 'mane' || lod !== 'full' ? [g.principal] : g.transcripts);
  const labelW = (s: string) => textWidth(s) + 10;
  const inView = (a: number, b: number) => Math.max(a, b) >= X0 - 2 && Math.min(a, b) <= X1 + 2;
  /** sticky heights: a track grows at once but shrinks only once the view rests (no jumping while panning) */
  const sticky = useRef(new Map<string, number>());
  const keep = (key: string, h: number) => { const prev = sticky.current.get(key) ?? 0; const v = interacting ? Math.max(h, prev) : h; sticky.current.set(key, v); return v; };

  // reference genes (their lanes kept from one drawing to the next)
  const refLanes = useRef(new Map<string, number>()), derLanes = useRef(new Map<string, Map<string, number>>());
  const refItems = (() => {
    if (lod === 'density') {
      const worth = visGenes.filter(g => worthSym(g.symbol) && (!filtering || passes(g.symbol)));
      return pack(worth.map(g => { const a = Math.max(X0 + 2, refX(g.start)); const lw = labelW(g.symbol) + badgeW(g.symbol); return { g, x0: a, x1: a + lw, xg: a + 3, h: LABEL_H + 2, worth: true, lw }; }), 3);
    }
    return pack(visGenes.filter(g => inView(refX(g.start), refX(g.end)) && (!filtering || passes(g.symbol))).map(g => {
      const a = Math.max(X0 + 2, refX(g.start)), b = refX(g.end), lw = labelW(g.symbol) + (panelDot(g.symbol) ? 10 : 0) + badgeW(g.symbol);
      return { g, x0: a, x1: Math.max(b, a + lw), xg: Math.max(b, a + 3), h: LABEL_H + (lod === 'full' ? shown(g).length * TX_H : 8), worth: worthSym(g.symbol), lw };
    }), MAX_LANES, refLanes.current, it => it.g.symbol, interacting);
  })();
  /** a density profile of the genes, for a chromosome seen whole */
  const density = useMemo(() => {
    if (lod !== 'density') return null;
    const BIN = 3, n = Math.ceil(PW / BIN), counts = new Float32Array(n);
    for (const g of view.genes) {
      if (filtering && !passes(g.symbol)) continue;
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
    // a chromosome seen whole: as on the reference, only the genes worth a label (changed by the event, in the panel,
    // dosage-sensitive, LoF-intolerant), each at its place on the rebuilt chromosome, over a density profile (below)
    const dens = lod === 'density';
    {
      const keepC = (ch: Chain) => ch.parts.every(p => (mode === 'mane' || lod !== 'full' ? p.tx === impactOf.get(p.tx.gene)?.gene.principal : true));
      const copyNo = new Map<string, number>();
      for (const ch of story.chains) {
        if (story.ders[ch.der] !== d || !keepC(ch)) continue;
        const geneList = [...new Set(ch.parts.map(p => p.tx.gene))];
        if (filtering && !geneList.some(passes)) continue;
        if (dens && geneList.length === 1 && !worthSym(geneList[0])) continue;
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
        const lw = labelW(label) + (fusion ? 80 : 0) + (geneList.length === 1 ? (panelDot(geneList[0]) ? 10 : 0) + badgeW(geneList[0]) : 0);
        list.push({ der: d, chains: [ch], genes: geneList, label, x0, x1: dens ? x0 + lw : Math.max(x1, x0 + lw), xg: Math.max(x1, x0 + 3), h: dens ? LABEL_H + 2 : LABEL_H + (lod === 'full' ? TX_H : 8), fusion, worth: geneList.length > 1 || worthSym(geneList[0]), lw });
      }
      // genes of the view in the flanks (the chromosome on either side of the event's window), carried as they are
      for (const dp of d.pieces) {
        if (!dp.seg.first && !dp.seg.last) continue;
        const [lo, hi] = pieceRange(dp, L);
        for (const g of view.genes) {
          if (impactOf.has(g.symbol) || g.start < lo || g.end > hi || (g.end > win.start && g.start < win.end) || (filtering && !passes(g.symbol)) || (dens && !worthSym(g.symbol))) continue;
          const xa = derX(d, toDer(dp, g.start)), xb = derX(d, toDer(dp, g.end));
          if (!inView(xa, xb)) continue;
          const x0 = Math.max(X0 + 2, Math.min(xa, xb)), x1 = Math.max(xa, xb), lw = labelW(g.symbol) + (panelDot(g.symbol) ? 10 : 0) + badgeW(g.symbol);
          list.push({ der: d, chains: [], genes: [g.symbol], label: g.symbol, x0, x1: dens ? x0 + lw : Math.max(x1, x0 + lw), xg: Math.max(x1, x0 + 3), h: dens ? LABEL_H + 2 : LABEL_H + (lod === 'full' ? TX_H : 8), fusion: null, worth: worthSym(g.symbol), lw, flank: { g, dp } });
        }
      }
    }
    if (!derLanes.current.has(d.name)) derLanes.current.set(d.name, new Map());
    if (dens) return pack(list, 3);
    return pack(list, MAX_LANES, derLanes.current.get(d.name), b => `${b.label}|${b.flank ? 'flank' : b.chains[0]?.parts[0]?.dp.index ?? ''}`, interacting);
  });
  /**
   * The density profile of each rebuilt chromosome, for a chromosome seen whole: every gene of the view counted where it
   * lies on the rebuilt chromosome (through each piece holding it), so a deleted gene drops out and those beyond a
   * junction count at their new place.
   */
  const derDensity = lod !== 'density' ? ders.map(() => null) : ders.map(d => {
    const BIN = 3, n = Math.ceil(PW / BIN), counts = new Float32Array(n);
    for (const g of view.genes) {
      if (filtering && !passes(g.symbol)) continue;
      for (const dp of d.pieces) {
        const [lo, hi] = pieceRange(dp, L);
        const a0 = Math.max(lo, g.start), b0 = Math.min(hi, g.end);
        if (b0 <= a0) continue;
        const xa = derX(d, toDer(dp, a0)), xb = derX(d, toDer(dp, b0));
        const ia = Math.floor((Math.min(xa, xb) - X0) / BIN), ib = Math.floor((Math.max(xa, xb) - X0) / BIN);
        for (let i = Math.max(0, ia); i <= Math.min(n - 1, ib); i++) counts[i]++;
      }
    }
    let max = 1; for (const v of counts) max = Math.max(max, v);
    let path = `M${X0},0`;
    for (let i = 0; i < n; i++) path += `L${X0 + i * BIN},${(-counts[i] / max).toFixed(3)}L${X0 + (i + 1) * BIN},${(-counts[i] / max).toFixed(3)}`;
    return { d: `${path}L${X1},0Z`, max };
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
  // ClinVar SVs (≥ 50 bp, pathogenic / likely pathogenic). Summed up by default: how many cover each point, losses above
  // a line and gains below it (inversions and complex on it), with the SVs matching a loss or gain of the event as bars
  // under the profile; expanded, every SV a bar in lanes
  const clinSvShown = showClinSv && !!loaders.getClinvar;
  const clinMatch = (v: ClinvarVariant) => cnv.some(c => (c.ev.kind === 'loss' ? v.kind === 'loss' : v.kind === 'gain') && clinOverlap(c.ev, v) >= RO_MATCH);
  const clinSvInView = clinSvShown ? clinSvView.list.filter(v => refX(v.end) >= X0 && refX(v.start) <= X1) : [];
  // the most SVs overlapping over one point of the view (at the drawing's resolution)
  const clinSvStack = (() => {
    const n = Math.ceil(PW / CLIN_AGG_BIN), d = new Int32Array(n + 1);
    for (const v of clinSvInView) { const i0 = Math.max(0, Math.floor((refX(v.start) - X0) / CLIN_AGG_BIN)), i1 = Math.min(n - 1, Math.floor((refX(v.end) - X0) / CLIN_AGG_BIN)); if (i1 >= i0) { d[i0]++; d[i1 + 1]--; } }
    let c = 0, m = 0; for (let i = 0; i < n; i++) { c += d[i]; m = Math.max(m, c); }
    return m;
  })();
  const clinSvAll = clinSvForce === 'all' || clinSvStack <= CLIN_SV_STACK_MAX;
  const clinAgg = (() => {
    if (!clinSvShown || clinSvAll || !clinSvInView.length) return null;
    const n = Math.ceil(PW / CLIN_AGG_BIN), d = { loss: new Int32Array(n + 1), gain: new Int32Array(n + 1), other: new Int32Array(n + 1) };
    for (const v of clinSvInView) {
      const arr = v.kind === 'loss' ? d.loss : v.kind === 'gain' ? d.gain : d.other;
      const i0 = Math.max(0, Math.floor((refX(v.start) - X0) / CLIN_AGG_BIN)), i1 = Math.min(n - 1, Math.floor((refX(v.end) - X0) / CLIN_AGG_BIN));
      if (i1 < i0) continue;
      arr[i0]++; arr[i1 + 1]--;
    }
    const prof = (a: Int32Array) => { const o = new Int32Array(n); let c = 0; for (let i = 0; i < n; i++) { c += a[i]; o[i] = c; } return o; };
    const loss = prof(d.loss), gain = prof(d.gain), other = prof(d.other);
    let max = 1; for (let i = 0; i < n; i++) max = Math.max(max, loss[i], gain[i]);
    return { n, loss, gain, other, max, matched: clinSvInView.filter(clinMatch) };
  })();
  const clinSvItems: { v: ClinvarVariant; lane: number; a: number; b: number }[] = [];
  if (clinSvShown) {
    const ends: number[] = [];
    for (const v of [...(clinSvAll ? clinSvView.list : clinAgg?.matched ?? [])].sort((p, q) => p.start - q.start)) {
      const a = refX(v.start), b = Math.max(refX(v.end), a + 2);
      if (b < X0 || a > X1) continue;
      let lane = ends.findIndex(e => e + 2 < a);
      if (lane < 0) { if (ends.length < CLIN_SV_MAX_LANES) { lane = ends.length; ends.push(b); } else { lane = ends.indexOf(Math.min(...ends)); ends[lane] = Math.max(ends[lane], b); } }
      else ends[lane] = b;
      clinSvItems.push({ v, lane, a: cx(a), b: cx(b) });
    }
  }
  const yClinSv = clinSvShown ? y : -1;
  /** the bars' top: under the profile in the summary */
  const yClinBars = yClinSv + (clinAgg ? CLIN_AGG_H + 2 : 0);
  if (clinSvShown) y += keep(clinSvAll ? 'clinsv-all' : 'clinsv', Math.max(clinAgg ? 0 : 1, ...clinSvItems.map(d => d.lane + 1)) * CLIN_SV_LANE_H + (clinAgg ? CLIN_AGG_H + 2 : 0) + 12) + 2;
  // ClinVar SNVs and short indels, mirrored: loss of function rises above a line, missense hangs below it, in frame and
  // other sit on it; one mark per 3 px and class (its variants in its card), pathogenic filled over likely pathogenic hollow
  const clinSnvShown = showClinSnv && !!loaders.getClinvar;
  const snvH = CLIN_SNV_GENE_H;
  const isEsc = (v: ClinvarVariant) => clinSnvGene.nmd.get(v)?.call === 'escape';
  const clinSnvBins: { x: number; vs: ClinvarVariant[]; code: 'PG' | 'LP'; side: 'lof' | 'lofesc' | 'missense' | 'mid' }[] = [];
  // zoomed out: a mirrored histogram instead, per CLIN_SNV_BIN px: loss of function stacked above the line, missense below
  // (pathogenic at the base, likely pathogenic on it), in frame and other as a tick on the line
  const clinSnvSum = clinSnvShown && clinSnvMode === 'each' && span > CLIN_SNV_SUM_SPAN ? (() => {
    const bins = new Map<number, ClinvarVariant[]>();
    for (const v of clinSnvView.list) { const x = refX((v.start + v.end) / 2); if (x < X0 || x > X1) continue; const k = Math.floor((x - X0) / CLIN_SNV_BIN); bins.set(k, [...(bins.get(k) ?? []), v]); }
    const out = [...bins].map(([k, vs]) => {
      const n = (e: string, c?: string, esc?: boolean) => vs.filter(v => (e === 'mid' ? v.effect !== 'lof' && v.effect !== 'missense' : v.effect === e) && (!c || v.code === c) && (esc == null || isEsc(v) === esc)).length;
      return { k, x: X0 + k * CLIN_SNV_BIN, vs, lofP: n('lof', 'PG', false), lofL: n('lof', 'LP', false), escP: n('lof', 'PG', true), escL: n('lof', 'LP', true), misP: n('missense', 'PG'), misL: n('missense', 'LP'), mid: n('mid') };
    });
    let max = 1, midMax = 1; for (const b of out) { max = Math.max(max, b.lofP + b.lofL + b.escP + b.escL, b.misP + b.misL); midMax = Math.max(midMax, b.mid); }
    return { bins: out, max, midMax };
  })() : null;
  if (clinSnvShown && clinSnvMode === 'each' && !clinSnvSum) {
    const bins = new Map<string, ClinvarVariant[]>();
    for (const v of clinSnvView.list) {
      const x = refX((v.start + v.end) / 2); if (x < X0 || x > X1) continue;
      const side = v.effect === 'lof' ? (isEsc(v) ? 'lofesc' : 'lof') : v.effect === 'missense' ? v.effect : 'mid', k = `${side}:${Math.round((x - X0) / 3)}`;
      bins.set(k, [...(bins.get(k) ?? []), v]);
    }
    for (const [k, vs] of bins) { const [side, i] = k.split(':'); clinSnvBins.push({ x: X0 + Number(i) * 3, vs, code: vs.some(v => v.code === 'PG') ? 'PG' : 'LP', side: side as 'lof' | 'lofesc' | 'missense' | 'mid' }); }
    clinSnvBins.sort((p, q) => (p.code === q.code ? 0 : p.code === 'LP' ? -1 : 1));
  }
  // summed per gene: above the line a red circle (loss of function predicted to trigger NMD) and an orange one (predicted
  // to escape it), below a ▼ (missense), each with its count, centred on the gene's stretch in view; genes whose glyphs
  // would overlap merge into one (their counts summed)
  type SnvGroup = { genes: SvGene[]; vs: ClinvarVariant[]; nmd: number; esc: number; mis: number; rest: number; x: number; s: number; e: number; w: number };
  const clinSnvGroups: SnvGroup[] = [];
  if (clinSnvShown && clinSnvMode === 'gene') {
    const by = new Map<string, { g: SvGene | null; vs: ClinvarVariant[] }>();
    for (const v of clinSnvView.list) {
      const g = clinSnvGene.geneOf.get(v) ?? null, k = g ? `g:${g.symbol}` : `v:${v.id}:${v.start}`;
      if (!g && (v.end < vp.s || v.start > ve)) continue;
      if (g && (g.end < vp.s || g.start > ve)) continue;
      const cur = by.get(k); if (cur) cur.vs.push(v); else by.set(k, { g, vs: [v] });
    }
    const halfW = (q: Omit<SnvGroup, 'w'>) => {
      const lofW = (q.nmd ? 2 * snvR(q.nmd) : 0) + (q.esc ? 2 * snvR(q.esc) : 0) + (q.nmd && q.esc ? 2 : 0);
      return Math.max(lofW / 2, q.mis ? snvTri(q.mis) : 0, q.rest && !q.nmd && !q.esc && !q.mis ? 4 : 0) + (q.genes.length > 1 ? 9 : 0) + 1.5;
    };
    const make = (genes: SvGene[], vs: ClinvarVariant[], x: number, s0: number, e0: number): SnvGroup => {
      const q = { genes, vs, x, s: s0, e: e0, nmd: vs.filter(v => v.effect === 'lof' && !isEsc(v)).length, esc: vs.filter(v => v.effect === 'lof' && isEsc(v)).length, mis: vs.filter(v => v.effect === 'missense').length, rest: vs.filter(v => v.effect !== 'lof' && v.effect !== 'missense').length };
      return { ...q, w: halfW(q) };
    };
    const items = [...by.values()].map(({ g, vs }) => {
      const s0 = g ? Math.max(g.start, vp.s) : Math.min(...vs.map(v => v.start)), e0 = g ? Math.min(g.end, ve) : Math.max(...vs.map(v => v.end));
      return make(g ? [g] : [], vs, refX((s0 + e0) / 2), s0, e0);
    }).sort((a, b) => a.x - b.x);
    // a stack: merge the last two while their glyphs overlap (the merged glyph sits at their count-weighted centre)
    for (const it of items) {
      clinSnvGroups.push(it);
      while (clinSnvGroups.length > 1) {
        const b = clinSnvGroups[clinSnvGroups.length - 1], a = clinSnvGroups[clinSnvGroups.length - 2];
        if (b.x - a.x >= a.w + b.w + 2) break;
        clinSnvGroups.splice(-2, 2, make([...a.genes, ...b.genes], [...a.vs, ...b.vs], (a.x * a.vs.length + b.x * b.vs.length) / (a.vs.length + b.vs.length), Math.min(a.s, b.s), Math.max(a.e, b.e)));
      }
    }
  }
  const yClinSnv = clinSnvShown ? y : -1; if (clinSnvShown) y += keep('clinsnv', snvH) + 6;
  // AChro-Puce: one bar per region (upper half the deletion's class, lower half the duplication's), its name and classes after it
  const ACHRO_LANE_H = 30;
  const achroItems: { p: AchroPlaced; lane: number; a: number; b: number; text: string }[] = [];
  if (showAchro) {
    const ends: number[] = [];
    for (const p of [...achroPlaced].sort((u, v) => u.start - v.start)) {
      const a = cx(refX(p.start)), b = Math.max(cx(refX(p.end)), a + 3);
      if (b < X0 || a > X1) continue;
      const text = `${p.approx ? '≈ ' : ''}${p.r.name}`;
      // the labels follow each type's bar: the lane is free again after the farthest one
      const right = Math.max(...Object.values(p.types).map(sp => { const bb = Math.max(cx(refX(sp.end)), cx(refX(sp.start)) + 3); return Math.max(bb, X0) + 6; })) + textWidth(text) * 0.9 + 110;
      let lane = ends.findIndex(e => e + 6 < a);
      if (lane < 0) { lane = ends.length; ends.push(right); } else ends[lane] = right;
      achroItems.push({ p, lane, a, b, text });
    }
  }
  const yAchro = showAchro ? y : -1; if (showAchro) y += keep('achro', Math.max(1, ...achroItems.map(d => d.lane + 1)) * ACHRO_LANE_H) + 8;
  const yRefBand = rec ? y : -1; if (rec) y += BAND_H;
  const yRib = y; if (ders.length) y += RIB_H;
  const derY = ders.map((d, i) => {
    const band = y; y += BAND_H;
    const ruler = y; y += 24;
    const genesY = y; y += keep(`der${i}`, lod === 'density' ? derPacked[i].height + 40 : derPacked[i].height) + 4;
    const tad = domains.length ? y : -1; if (domains.length) y += TAD_H;
    y += 18;
    return { band, ruler, genes: genesY, tad };
  });
  const H = y + 6;

  // ---- cards: hover follows the pointer (moved directly, no redraw), a click keeps one open ----
  const boxPos = (e: React.MouseEvent) => { const r = boxRef.current!.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  /** below the pointer when it fits, else above; a card taller than either side is kept inside the figure */
  const fitY = (el: HTMLDivElement, yy: number, h: number) => {
    const ch = el.offsetHeight || 160;
    if (yy + 14 + ch <= h) { el.style.bottom = ''; el.style.top = `${yy + 14}px`; }
    else if (yy - 14 - ch >= 0) { el.style.top = ''; el.style.bottom = `${h - yy + 14}px`; }
    else { el.style.bottom = ''; el.style.top = `${Math.max(0, h - ch)}px`; }
  };
  const placeTip = (x: number, yy: number) => {
    const el = tipEl.current, box = boxRef.current; if (!el || !box) return;
    const w = box.clientWidth, h = box.clientHeight;
    // next to the pointer, kept inside the figure by the card's own width (left of the pointer near the right edge)
    const cw = el.offsetWidth || 400;
    el.style.left = `${Math.max(0, x + 14 + cw <= w ? x + 14 : Math.max(0, x - 14 - cw))}px`;
    fitY(el, yy, h);
  };
  const lastTipPos = useRef({ x: 0, y: 0 });
  const hoverOn = (key: string, dgv?: PopulationSv, clin?: ClinvarVariant[]) => (e: React.MouseEvent) => {
    if (dragRef.current?.moved) return;
    const p = boxPos(e); lastTipPos.current = p;
    const achro = key.startsWith('achro:') ? achroPlaced.find(x => x.r.id === key.split(':')[1]) : undefined, achroType = achro ? key.split(':')[2] as 'del' | 'dup' : undefined;
    if (tip?.gene !== key) setTip({ ...p, gene: key, dgv, achro, achroType, clin }); else placeTip(p.x, p.y);
  };
  const hoverOff = () => setTip(null);
  useLayoutEffect(() => { if (tip) placeTip(lastTipPos.current.x, lastTipPos.current.y); });
  const pinOn = (key: string, dgv?: PopulationSv, clin?: ClinvarVariant[]) => (e: React.MouseEvent) => {
    e.stopPropagation();
    if (suppressClick.current) return;
    const p = boxPos(e);
    const achro = key.startsWith('achro:') ? achroPlaced.find(x => x.r.id === key.split(':')[1]) : undefined, achroType = achro ? key.split(':')[2] as 'del' | 'dup' : undefined;
    if (!dgv && !achro && !clin) setLinked(key);
    setPinned(cur => (cur?.gene === key ? null : { ...p, gene: key, dgv, achro, achroType, clin }));
  };
  // the card points at its gene where it was: it closes when the view moves
  useEffect(() => { if (pinned && !anim.current.fly) setPinned(null); }, [Math.round(vp.s / Math.max(1, span) * 50), Math.round(Math.log(span) * 20)]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- gestures: drag to pan (a fling goes on), Ctrl / ⌘ + wheel or pinch to zoom, two-finger swipe to pan, double-click to zoom in,
  //      drag along the ruler to zoom on a stretch, click or drag on the ideogram to go there ----
  const dragRef = useRef<{ kind: 'pan' | 'band' | 'ideo' | 'select'; x0: number; vp0: Vp; moved: boolean; samples: { x: number; t: number }[]; id: number } | null>(null);
  /** the click that ends a drag opens no card */
  const suppressClick = useRef(false);
  const [band, setBand] = useState<{ a: number; b: number } | null>(null);
  /** an interval picked with Ctrl / ⌘ + drag: highlighted on the reference axis (it follows the view) until closed */
  const [sel, setSel] = useState<{ s: number; e: number } | null>(null);
  const [selCopied, setSelCopied] = useState(false);
  // an interval asked for (typed coordinates): highlighted as a Ctrl / ⌘ + drag would, with its bar
  useEffect(() => { if (highlight) { setSel({ s: highlight.start, e: highlight.end }); setSelCopied(false); } }, [highlight?.key]); // eslint-disable-line react-hooks/exhaustive-deps
  const toSvgX = (clientX: number) => { const r = svgRef.current!.getBoundingClientRect(); return ((clientX - r.left) / r.width) * VW; };
  const toSvgY = (clientY: number) => { const r = svgRef.current!.getBoundingClientRect(); return ((clientY - r.top) / r.height) * H; };
  const ideoPos = (x: number) => Math.max(0, Math.min(L, ((x - X0) / PW) * L));
  const onPointerDown = (e: React.PointerEvent) => {
    if (e.button !== 0) return;
    anim.current = {}; suppressClick.current = false;
    const x = toSvgX(e.clientX), yy = toSvgY(e.clientY);
    const kind = yy < IDEO_Y + IDEO_H + 4 ? 'ideo' : (e.ctrlKey || e.metaKey) && x >= X0 && x <= X1 ? 'select' : yy >= yRuler && yy < yRuler + RULER_H ? 'band' : 'pan';
    dragRef.current = { kind, x0: x, vp0: vpRef.current, moved: false, samples: [{ x: e.clientX, t: performance.now() }], id: e.pointerId };
    if (kind === 'ideo') { const sp = vpRef.current.e - vpRef.current.s, p = ideoPos(x); flyTo(p - sp / 2, p + sp / 2, 320); }
  };
  const onPointerMove = (e: React.PointerEvent) => {
    const x = toSvgX(e.clientX);
    const d = dragRef.current; if (!d || d.id !== e.pointerId) return;
    if (!d.moved && Math.abs(x - d.x0) < 4) return;
    if (!d.moved) { d.moved = true; svgRef.current?.setPointerCapture(e.pointerId); setTip(null); }
    d.samples.push({ x: e.clientX, t: performance.now() }); if (d.samples.length > 6) d.samples.shift();
    if (d.kind === 'pan') { const sh = -(x - d.x0) * bppRef.current; moveTo(d.vp0.s + sh, d.vp0.e + sh); }
    else if (d.kind === 'band' || d.kind === 'select') setBand({ a: Math.max(X0, Math.min(X1, d.x0)), b: Math.max(X0, Math.min(X1, x)) });
    else { const sp = vpRef.current.e - vpRef.current.s, p = ideoPos(x); moveTo(p - sp / 2, p + sp / 2, p); }
  };
  const onPointerUp = (e: React.PointerEvent) => {
    const d = dragRef.current; if (!d || d.id !== e.pointerId) return;
    dragRef.current = null; suppressClick.current = d.moved;
    if (d.kind === 'select') {
      const a = band ? Math.min(band.a, band.b) : d.x0, b = band ? Math.max(band.a, band.b) : d.x0; setBand(null);
      if (d.moved && b - a > 2) { const s0 = vpRef.current.s, k = bppRef.current; setSel({ s: Math.max(0, Math.round(s0 + (a - X0) * k)), e: Math.min(L, Math.round(s0 + (b - X0) * k)) }); setSelCopied(false); }
      return;
    }
    if (d.kind === 'band' && band) {
      const a = Math.min(band.a, band.b), b = Math.max(band.a, band.b); setBand(null);
      if (b - a > 4) { const s0 = vpRef.current.s, k = bppRef.current; flyTo(s0 + (a - X0) * k, s0 + (b - X0) * k); }
      return;
    }
    if (d.kind === 'pan' && d.moved && d.samples.length > 1) {
      const p = d.samples[0], q = d.samples[d.samples.length - 1], dt = q.t - p.t;
      if (dt > 0 && performance.now() - q.t < 60) {
        const r = svgRef.current!.getBoundingClientRect(), pxPerMs = (q.x - p.x) / dt * (VW / r.width);
        if (Math.abs(pxPerMs) > 0.25) { anim.current = { fling: { v: -pxPerMs * bppRef.current, t: performance.now() } }; kick(); }
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
        const f = Math.exp(e.deltaY * unit * 0.0022), at = v.s + (Math.max(X0, Math.min(X1, x)) - X0) * bppRef.current;
        moveTo(at - (at - v.s) * f, at + (v.e - at) * f, at);
      } else if (Math.abs(e.deltaX) > Math.abs(e.deltaY) || e.shiftKey) {
        e.preventDefault();
        const dx = (e.shiftKey && !e.deltaX ? e.deltaY : e.deltaX) * unit, sh = dx / r.width * VW * bppRef.current;
        moveTo(v.s + sh, v.e + sh);
      }
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, [txs.loading, L]); // eslint-disable-line react-hooks/exhaustive-deps
  const onDoubleClick = (e: React.MouseEvent) => {
    const x = toSvgX(e.clientX); if (x < X0 || x > X1) return;
    const v = vpRef.current, at = v.s + (x - X0) * bppRef.current, f = e.shiftKey ? 2 : 0.5;
    flyTo(at - (at - v.s) * f, at + (v.e - at) * f, 320, at);
  };
  const onKeyDown = (e: React.KeyboardEvent) => {
    const v = vpRef.current, sp = v.e - v.s;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') { e.preventDefault(); const sh = (e.key === 'ArrowLeft' ? -1 : 1) * sp * 0.3; flyTo(v.s + sh, v.e + sh, 260, 'keep'); }
    else if (e.key === '+' || e.key === '=') { e.preventDefault(); flyTo(v.s + sp * 0.25, v.e - sp * 0.25, 260); }
    else if (e.key === '-' || e.key === '_') { e.preventDefault(); flyTo(v.s - sp * 0.5, v.e + sp * 0.5, 260); }
    else if ((e.key === 'e' || e.key === 'E') && rec) { e.preventDefault(); eventView(); }
    else if (e.key === 'Escape' && sel) { e.preventDefault(); e.stopPropagation(); setSel(null); }
  };
  const zoomToSel = () => { if (sel) { const pad = Math.max(100, (sel.e - sel.s) * 0.04); flyTo(sel.s - pad, sel.e + pad); } };
  const selText = sel ? `${chrom}:${fmt(sel.s + 1)}-${fmt(sel.e)}` : '';

  // ---- ruler ----
  const step = niceStep(ve - vp.s, 9);
  const ticks: number[] = [];
  for (let t = Math.max(0, Math.ceil(vp.s / step) * step); t <= ve; t += step) ticks.push(t);
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

  /** the gene whose card is open: outlined, and joined to its copies on the rebuilt chromosome */
  const linkSym = linked;
  /** each copy of a gene's principal model on the rebuilt chromosome: whole or which part, inverted or not */
  const copiesOf = (gi: GeneImpact): string[] => story.chains.filter(ch => ch.parts.some(p => p.tx === gi.principal.tx)).map(ch => {
    const parts = ch.parts.filter(p => p.tx === gi.principal.tx), first = parts[0], last = parts[parts.length - 1];
    const five = first.tx.strand > 0 ? !first.cutLo || first.edgeLo : !first.cutHi || first.edgeHi;
    const three = last.tx.strand > 0 ? !last.cutHi || last.edgeHi : !last.cutLo || last.edgeLo;
    const what = five && three ? 'whole' : !five && !three ? 'inner part' : !five ? '3′ part' : '5′ part';
    const fused = ch.parts.some(p => p.tx.gene !== gi.gene.symbol) ? ` fused to ${[...new Set(ch.parts.filter(p => p.tx.gene !== gi.gene.symbol).map(p => p.tx.gene))].join(', ')}` : '';
    return `${what}${parts.some(p => p.dp.piece.reversed) ? ', inverted' : ''}${fused}`;
  });

  // ---- the SVG ----
  const vl = Math.max(2, ix(Math.min(L, ve)) - ix(Math.max(0, vp.s)));
  const svg = (
    <svg ref={svgRef} viewBox={`0 0 ${VW} ${H}`} width="100%" style={{ aspectRatio: `${VW} / ${H}`, cursor: dragRef.current?.moved && dragRef.current.kind === 'pan' ? 'grabbing' : 'grab', fontFamily: FONT, touchAction: 'pan-y' }} className="block select-none"
      onPointerDown={onPointerDown} onPointerMove={onPointerMove} onPointerUp={onPointerUp} onPointerCancel={onPointerUp}
      onContextMenu={e => { if (e.ctrlKey) e.preventDefault(); /* Ctrl + click on a Mac: the selection, not the menu */ }}
      onClick={() => { if (!suppressClick.current) { setPinned(null); setLinked(null); } }} onDoubleClick={onDoubleClick}>
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
      <path d={`M${ix(Math.max(0, vp.s))},${IDEO_Y + IDEO_H} L${X0},${IDEO_Y + IDEO_H + LENS_H} H${X1} L${ix(Math.max(0, vp.s)) + vl},${IDEO_Y + IDEO_H} Z`} fill="url(#svb-lens)" />
      <rect x={ix(Math.max(0, vp.s)) - 1} y={IDEO_Y - 2} width={vl + 2} height={IDEO_H + 4} rx={3} fill={C.accent} fillOpacity={0.12} stroke={C.accent} strokeWidth={1.4} />
      {/* track names */}
      <g fontSize={10.5} fill={C.muted}>
        {domains.length > 0 && <text x={14} y={yRefTad + TAD_H - 6}>TADs</text>}
        <text x={14} y={yRefGenes + 12} fontWeight={600} fill={C.text}>genes</text>
        <text x={14} y={yRefGenes + 25} fontSize={9.5}>{lod === 'density' ? 'density' : lod === 'bars' || mode === 'mane' ? 'MANE Select' : 'all NM_'}</text>
        {view.loading && <text x={14} y={yRefGenes + 38} fontSize={9} fill={C.faint}>loading…</text>}
        {yDgv >= 0 && <text x={14} y={yDgv + 8} fontWeight={600} fill={C.text}>DGV Gold<title>DGV Gold Standard: copy-number variants of the general population (curated, found in two studies and two samples at least)</title></text>}
        {yClinSv >= 0 && clinSvInView.length > 0 && clinSvStack > CLIN_SV_STACK_MAX && <text x={14} y={yClinSv + 20} fontSize={9.5} fill={C.accent} style={{ cursor: 'pointer' }}
          onClick={e => { e.stopPropagation(); const v = clinSvForce === 'all' ? 'auto' : 'all'; setClinSvForce(v); store.set('sashimi.svClinvarSvMode', v); }}>
          {clinSvAll ? `▾ summary` : `▸ all SVs (${clinSvInView.length})`}<title>{clinSvAll ? `Sum the SVs up (up to ${clinSvStack} lie over one point): how many cover each point (losses above, gains below), the ones matching the event as bars` : 'Every ClinVar SV in view, one bar each'}</title></text>}
        {yClinSv >= 0 && <text x={14} y={yClinSv + 8} fontWeight={600} fill={C.text}>ClinVar SV<title>ClinVar variants of 50 bp or more classified pathogenic or likely pathogenic (UCSC track clinvarCnv)</title></text>}
        {yClinSnv >= 0 && <g fontSize={9} pointerEvents="none">
          <text x={14} y={yClinSnv + 8}><tspan fill={CLIN_EFFECT_COLOR.lof}>● LoF</tspan><tspan fill={C.faint}> · </tspan><tspan fill={CLIN_LOF_ESC}>● NMD escape</tspan></text>
          <text x={14} y={yClinSnv + snvH - 2} fill={CLIN_EFFECT_COLOR.missense}>▼ missense</text></g>}
        {yClinSnv >= 0 && <text x={14} y={yClinSnv + snvH / 2 - 2} fontWeight={600} fill={C.text}>ClinVar SNV<title>ClinVar variants under 50 bp (SNVs, short indels) classified pathogenic or likely pathogenic (UCSC track clinvarMain)</title></text>}
        {yClinSnv >= 0 && clinSnvMode === 'gene' && <text x={14} y={yClinSnv + snvH / 2 + 11} fontSize={9.5} fill={C.accent} style={{ cursor: 'pointer' }}
          onClick={e => { e.stopPropagation(); setPinned(null); setTip(null); setClinSnvMode('each'); store.set('sashimi.svClinvarSnvMode', 'each'); }}>▸ each variant<title>Each variant at its place instead of the counts per gene</title></text>}
        {yClinSnv >= 0 && clinSnvMode === 'each' && <text x={14} y={yClinSnv + snvH / 2 + 11} fontSize={9.5} fill={C.accent} style={{ cursor: 'pointer' }}
          onClick={e => { e.stopPropagation(); setPinned(null); setTip(null); setClinSnvMode('gene'); store.set('sashimi.svClinvarSnvMode', 'gene'); }}>▾ per gene<title>Count the variants per gene: loss of function above (predicted NMD, predicted escape), missense below</title></text>}
        {yAchro >= 0 && <text x={14} y={yAchro + 10} fontWeight={600} fill={C.text}>AChro-Puce<title>{`${ACHRO_SOURCE.title}, ${ACHRO_SOURCE.version}`}</title></text>}
        {rec && <><text x={14} y={yRefBand + 9} fontSize={9.5}>reference</text><text x={14} y={yRefBand + 20} fontSize={9.5} fontFamily={MONO}>{kb(L)}</text></>}
        {ders.map((d, i) => { const len = derExtent(d).len, diff = len - L; return (
          <g key={d.name}><text x={14} y={derY[i].band + 9} fontWeight={700} fill={C.accent}>{d.name}</text>
            <text x={14} y={derY[i].band + 21} fontSize={9.5} fontFamily={MONO}>{kb(len)}</text>
            {Math.abs(diff) >= 1 && <text x={14} y={derY[i].band + 32} fontSize={9.5} fontWeight={700} fontFamily={MONO} fill={diff > 0 ? '#15803d' : C.lost}>{diff > 0 ? '+' : '−'}{kb(Math.abs(diff))}</text>}
          </g>); })}
      </g>
      {/* ruler: drag along it to zoom on a stretch */}
      <rect x={X0} y={yRuler} width={PW} height={RULER_H - 4} fill="#eef2f6" rx={4} style={{ cursor: 'col-resize' }}><title>Drag along the ruler to zoom on a stretch</title></rect>
      <g clipPath="url(#svb-clip)">
        {ticks.map(t => <g key={`tk${t}`}><line x1={refX(t)} x2={refX(t)} y1={yRuler + RULER_H - 9} y2={yRuler + RULER_H - 4} stroke={C.faint} />{refX(t) > X0 + 28 && refX(t) < X1 - 28 && <text x={refX(t)} y={yRuler + 11} fontSize={9} textAnchor="middle" fill={C.muted} fontFamily={MONO}>{tickText(t, step)}</text>}</g>)}
        {/* grid lines through the tracks */}
        {ticks.map(t => <line key={`gl${t}`} x1={refX(t)} x2={refX(t)} y1={yRuler + RULER_H} y2={rec ? yRefBand : H - 8} stroke={C.grid} strokeWidth={0.6} />)}
        {/* reference TADs */}
        {domains.filter(d => d.end > vp.s && d.start < ve).map((d, i) => {
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
              {labeled && badge(g.symbol, lx + (dot ? 10 : 0) + textWidth(g.symbol) + 4, yy + 9.5)}
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
        {/* ClinVar SVs: red losses, blue gains, purple inversions and complex; pathogenic solid, likely pathogenic paler; outlined when one matches a loss or gain of the event */}
        {yClinSv >= 0 && clinAgg && (() => {
          // the profile: losses above the line, gains below, inversions and complex on it; hover a point for the SVs over it
          const mid = yClinSv + CLIN_AGG_H / 2, hh = CLIN_AGG_H / 2 - 1.5, path = (a: Int32Array, dir: 1 | -1) => {
            let d = `M${X0},${mid}`;
            for (let i = 0; i < clinAgg.n; i++) { const yv = mid - dir * (a[i] / clinAgg.max) * hh; d += `L${X0 + i * CLIN_AGG_BIN},${yv.toFixed(1)}L${X0 + (i + 1) * CLIN_AGG_BIN},${yv.toFixed(1)}`; }
            return `${d}L${X1},${mid}Z`;
          };
          const others: JSX.Element[] = [];
          for (let i = 0; i < clinAgg.n; i++) if (clinAgg.other[i]) { let j = i; while (j + 1 < clinAgg.n && clinAgg.other[j + 1]) j++; others.push(<rect key={`o${i}`} x={X0 + i * CLIN_AGG_BIN} y={mid - 1} width={(j - i + 1) * CLIN_AGG_BIN} height={2} fill={CLIN_KIND_COLOR.inv} />); i = j; }
          const at = (e: React.MouseEvent) => {
            const x = toSvgX(e.clientX), pos = vp.s + (x - X0) * bppRef.current, tol = 2 * bppRef.current;
            const vs = clinSvInView.filter(v => v.start <= pos + tol && v.end >= pos - tol).sort((p, q) => (p.end - p.start) - (q.end - q.start));
            return { key: `cvagg:${Math.round((x - X0) / 4)}`, vs, pos: Math.round(pos) };
          };
          // ClinVar cards open on a click only (no card following the pointer over these tracks)
          const onClick = (e: React.MouseEvent) => { e.stopPropagation(); if (suppressClick.current) return; const r = at(e); if (!r.vs.length) return; setPinned(cur => (cur?.gene === r.key ? null : { ...boxPos(e), gene: r.key, clin: r.vs, clinAt: r.pos })); };
          return (
            <g style={{ cursor: 'pointer' }} onClick={onClick}>
              <rect x={X0} y={yClinSv} width={PW} height={CLIN_AGG_H} fill="transparent" />
              <path d={path(clinAgg.loss, 1)} fill={CLIN_KIND_COLOR.loss} fillOpacity={0.55} />
              <path d={path(clinAgg.gain, -1)} fill={CLIN_KIND_COLOR.gain} fillOpacity={0.55} />
              <line x1={X0} x2={X1} y1={mid} y2={mid} stroke={C.line} strokeWidth={0.8} />
              {others}
              <text x={X1 - 3} y={yClinSv + 8} fontSize={8.5} textAnchor="end" fill={C.faint} pointerEvents="none">{`up to ${clinSvStack} SVs over a point · axis ${clinAgg.max} per side`}</text>
            </g>
          );
        })()}
        {yClinSv >= 0 && clinSvItems.map(({ v, lane, a, b }) => {
          const yy = yClinBars + lane * CLIN_SV_LANE_H, col = CLIN_KIND_COLOR[v.kind] ?? CLIN_KIND_COLOR.other, key = `cv:${v.id}:${v.start}`;
          const same = clinMatch(v);
          return (
            <g key={key} style={{ cursor: 'pointer' }} onClick={pinOn(key, undefined, [v])}>
              <rect x={a} y={yy} width={b - a} height={CLIN_SV_LANE_H - 1} fill="transparent" />
              <rect x={a} y={yy + 1} width={Math.max(1.5, b - a)} height={6} rx={1} fill={col} fillOpacity={v.code === 'PG' ? 0.9 : 0.4} stroke={same ? '#111827' : 'none'} strokeWidth={0.8} />
            </g>
          );
        })}
        {yClinSv >= 0 && !clinSvItems.length && !clinAgg && <text x={X0 + 6} y={yClinSv + 8} fontSize={9.5} fill={C.faint}>{span > CLIN_SV_MAX_SPAN ? `zoom in below ${kb(CLIN_SV_MAX_SPAN)} to see the ClinVar SVs` : clinSvView.loading ? 'loading ClinVar…' : clinSvView.error ? `ClinVar unavailable: ${clinSvView.error}` : 'no pathogenic or likely pathogenic ClinVar SV here'}</text>}
        {/* ClinVar SNVs, mirrored: ▲ loss of function above the line, ● missense below, ◆ in frame and other on it */}
        {yClinSnv >= 0 && <line x1={X0} x2={X1} y1={yClinSnv + snvH / 2} y2={yClinSnv + snvH / 2} stroke={C.line} strokeWidth={1} />}
        {yClinSnv >= 0 && clinSnvSum && clinSnvSum.bins.map(b => {
          const mid = yClinSnv + snvH / 2, hh = snvH / 2 - 4, k = hh / clinSnvSum.max, w = CLIN_SNV_BIN - 0.8;
          const range: [number, number] = [Math.floor(vp.s + (b.x - X0) * bppRef.current), Math.ceil(vp.s + (b.x + CLIN_SNV_BIN - X0) * bppRef.current)];
          const vs = [...b.vs].sort((p, q) => ['lof', 'missense', 'inframe', 'other'].indexOf(p.effect) - ['lof', 'missense', 'inframe', 'other'].indexOf(q.effect) || (p.code === q.code ? 0 : p.code === 'PG' ? -1 : 1));
          const key = `cvsum:${b.k}`;
          return (
            <g key={key} style={{ cursor: 'pointer' }}
              onClick={e => { e.stopPropagation(); if (suppressClick.current) return; const pt = boxPos(e); setPinned(cur => (cur?.gene === key ? null : { ...pt, gene: key, clin: vs, clinRange: range })); }}>
              <rect x={b.x} y={yClinSnv} width={CLIN_SNV_BIN} height={snvH} fill="transparent" />
              {b.lofP > 0 && <rect x={b.x} y={mid - b.lofP * k} width={w} height={b.lofP * k} fill={CLIN_EFFECT_COLOR.lof} />}
              {b.lofL > 0 && <rect x={b.x} y={mid - (b.lofP + b.lofL) * k} width={w} height={b.lofL * k} fill={CLIN_EFFECT_COLOR.lof} fillOpacity={0.4} />}
              {b.escP > 0 && <rect x={b.x} y={mid - (b.lofP + b.lofL + b.escP) * k} width={w} height={b.escP * k} fill={CLIN_LOF_ESC} />}
              {b.escL > 0 && <rect x={b.x} y={mid - (b.lofP + b.lofL + b.escP + b.escL) * k} width={w} height={b.escL * k} fill={CLIN_LOF_ESC} fillOpacity={0.4} />}
              {b.misP > 0 && <rect x={b.x} y={mid} width={w} height={b.misP * k} fill={CLIN_EFFECT_COLOR.missense} />}
              {b.misL > 0 && <rect x={b.x} y={mid + b.misP * k} width={w} height={b.misL * k} fill={CLIN_EFFECT_COLOR.missense} fillOpacity={0.4} />}
              {b.mid > 0 && <rect x={b.x} y={mid - 1.5} width={w} height={3} fill={CLIN_EFFECT_COLOR.inframe} fillOpacity={0.25 + 0.6 * b.mid / clinSnvSum.midMax} />}
            </g>
          );
        })}
        {yClinSnv >= 0 && clinSnvSum && <text x={X1 - 3} y={yClinSnv + 8} fontSize={8.5} textAnchor="end" fill={C.faint} pointerEvents="none">{`up to ${clinSnvSum.max} per ${kb(CLIN_SNV_BIN * bppRef.current)} · zoom in below ${kb(CLIN_SNV_SUM_SPAN)} for each variant`}</text>}
        {yClinSnv >= 0 && clinSnvBins.map(({ x, vs, code, side }) => {
          const key = `cvs:${side}:${vs[0].id}:${vs[0].start}`, mid = yClinSnv + snvH / 2;
          const up = side === 'lof' || side === 'lofesc';
          const col = side === 'lofesc' ? CLIN_LOF_ESC : CLIN_EFFECT_COLOR[side === 'mid' ? (vs.some(v => v.effect === 'inframe') ? 'inframe' : 'other') : side];
          const filled = code === 'PG', fill = filled ? col : '#ffffff', sw = filled ? 0.9 : 1.4, stroke = filled ? '#ffffff' : col;
          // a stem a little longer for a mark gathering several variants (at most the lane's half)
          const stem = Math.min(snvH / 2 - 6, 13 + 2.4 * Math.log2(vs.length));
          const r = (vs.length > 1 ? 4.1 : 3.5) + (filled ? 0.3 : 0);
          const head = up
            ? <circle cx={x} cy={mid - stem} r={r} fill={fill} stroke={stroke} strokeWidth={sw} />
            : side === 'missense'
              ? <path d={`M${x},${mid + stem + r * 1.25}L${x + r * 1.15},${mid + stem - r * 0.75}L${x - r * 1.15},${mid + stem - r * 0.75}Z`} fill={fill} stroke={stroke} strokeWidth={sw} strokeLinejoin="round" />
              : <path d={`M${x},${mid - r}L${x + r},${mid}L${x},${mid + r}L${x - r},${mid}Z`} fill={fill} stroke={stroke} strokeWidth={sw} />;
          return (
            <g key={key} style={{ cursor: 'pointer' }} onClick={pinOn(key, undefined, vs)}>
              <rect x={x - 3} y={up ? yClinSnv : side === 'missense' ? mid : mid - 5} width={6} height={side === 'mid' ? 10 : snvH / 2} fill="transparent" />
              {side !== 'mid' && <line x1={x} x2={x} y1={mid} y2={up ? mid - stem : mid + stem} stroke={col} strokeOpacity={0.55} strokeWidth={1.1} />}
              {head}
            </g>
          );
        })}
        {yClinSnv >= 0 && clinSnvGroups.map(gp => {
          const midY = yClinSnv + snvH / 2, key = `cvg:${gp.genes.map(g => g.symbol).join(',') || `${gp.vs[0].id}:${gp.vs[0].start}`}`;
          const names = gp.genes.map(g => g.symbol), tx = gp.genes.map(g => `${g.principal.id}${g.principal.status === 'mane' ? ' (MANE Select)' : ''}`);
          const order = (v: ClinvarVariant) => (v.effect === 'lof' ? (isEsc(v) ? 1 : 0) : v.effect === 'missense' ? 2 : 3) * 2 + (v.code === 'PG' ? 0 : 1);
          const vs = [...gp.vs].sort((a, b) => order(a) - order(b));
          const info = { names, tx, s: gp.s, e: gp.e };
          const rN = snvR(gp.nmd), rE = snvR(gp.esc), both = gp.nmd > 0 && gp.esc > 0;
          const cxN = both ? gp.x - rN - 1 : gp.x, cxE = both ? gp.x + rE + 1 : gp.x;
          const t = snvTri(gp.mis), th = t * 1.3, triFs = gp.mis > 99 ? 7.5 : 9;
          const num = (n: number, x: number, yy: number, big: number) => <text x={x} y={yy} fontSize={n > 99 ? 7.5 : big > 9 ? 9.5 : 8.5} fontWeight={700} textAnchor="middle" fill="#fff" pointerEvents="none">{n}</text>;
          const gx0 = Math.max(X0, refX(gp.s)), gx1 = Math.min(X1, refX(gp.e));
          const lw = Math.max(gp.w * 2, 10);
          return (
            <g key={key} style={{ cursor: 'pointer' }}
              onClick={e => { e.stopPropagation(); if (suppressClick.current) return; const pt = boxPos(e); setPinned(cur => (cur?.gene === key ? null : { ...pt, gene: key, clin: vs, clinGenes: info })); }}>
              {/* the gene's stretch on the line */}
              {gx1 > gx0 && <rect x={gx0} y={midY - 1.5} width={Math.max(gx1 - gx0, 1)} height={3} rx={1.5} fill={C.faint} fillOpacity={0.35} />}
              <rect x={gp.x - lw / 2} y={yClinSnv} width={lw} height={snvH} fill="transparent" />
              {gp.nmd > 0 && <><circle cx={cxN} cy={midY - 4 - rN} r={rN} fill={CLIN_EFFECT_COLOR.lof} stroke="#fff" strokeWidth={1} />{num(gp.nmd, cxN, midY - 4 - rN + 3.3, rN)}</>}
              {gp.esc > 0 && <><circle cx={cxE} cy={midY - 4 - rE} r={rE} fill={CLIN_LOF_ESC} stroke="#fff" strokeWidth={1} />{num(gp.esc, cxE, midY - 4 - rE + 3.3, rE)}</>}
              {gp.mis > 0 && <><path d={`M${gp.x - t},${midY + 4}L${gp.x + t},${midY + 4}L${gp.x},${midY + 4 + th}Z`} fill={CLIN_EFFECT_COLOR.missense} stroke="#fff" strokeWidth={1} strokeLinejoin="round" /><text x={gp.x} y={midY + 4 + 1.5 + triFs * 0.74} fontSize={triFs} fontWeight={700} textAnchor="middle" fill="#fff" stroke={CLIN_EFFECT_COLOR.missense} strokeWidth={2} paintOrder="stroke" pointerEvents="none">{gp.mis}</text></>}
              {!gp.nmd && !gp.esc && !gp.mis && <path d={`M${gp.x},${midY - 3.5}L${gp.x + 3.5},${midY}L${gp.x},${midY + 3.5}L${gp.x - 3.5},${midY}Z`} fill={CLIN_EFFECT_COLOR.inframe} />}
              {gp.genes.length > 1 && <text x={gp.x + gp.w - 1} y={midY + 3} fontSize={8} textAnchor="end" fill={C.faint} pointerEvents="none">{`×${gp.genes.length}`}</text>}
            </g>
          );
        })}
        {yClinSnv >= 0 && !clinSnvBins.length && !clinSnvSum?.bins.length && !clinSnvGroups.length && <text x={X0 + 6} y={yClinSnv + snvH / 2 - 4} fontSize={9.5} fill={C.faint}>{span > CLIN_SNV_MAX_SPAN ? `zoom in below ${kb(CLIN_SNV_MAX_SPAN)} to see the ClinVar SNVs` : clinSnvView.loading ? 'loading ClinVar…' : clinSnvView.error ? `ClinVar unavailable: ${clinSnvView.error}` : 'no pathogenic or likely pathogenic ClinVar SNV here'}</text>}
        {/* AChro-Puce recurrent CNVs */}
        {yAchro >= 0 && achroItems.map(({ p, lane, text }) => {
          // two bars, each its own target: the deletion (red) above, the duplication (blue) below, its class after it
          const yy = yAchro + lane * ACHRO_LANE_H;
          const bar = (type: 'del' | 'dup', y0: number, first: boolean) => {
            const row = p.r.rows.find(r2 => r2.type === type && r2.sex !== 'female'), sp = p.types[type];
            if (!row || !sp) return null;
            // each type on its own stretch (a deletion's critical gene, a duplication's smallest overlap)
            const a = cx(refX(sp.start)), b = Math.max(cx(refX(sp.end)), a + 3);
            const lx = Math.max(a, X0) + Math.max(0, Math.min(b - Math.max(a, X0), PW)) + 6;
            const key = `achro:${p.r.id}:${type}`, col = ACHRO_TYPE_COLOR[type];
            const cls = `${row.cls === 'PIEV' ? 'PIEV' : `class ${row.cls}`}${row.also ? ` + ${row.also.cls === 'PIEV' ? 'PIEV' : `class ${row.also.cls}`} (${row.also.label})` : ''}${p.r.rows.some(r2 => r2.sex) ? ' (male)' : ''}`;
            const labelW = (first ? textWidth(text) * 0.88 + 8 : 0) + 30 + cls.length * 5.4;
            return (
              <g key={key} style={{ cursor: 'pointer' }} onMouseEnter={hoverOn(key)} onMouseMove={hoverOn(key)} onMouseLeave={hoverOff} onClick={pinOn(key)}>
                <rect x={Math.min(a, lx - 6)} y={y0 - 2} width={Math.max(b - a, 6) + 6 + labelW} height={13} fill="transparent" />
                <rect x={a} y={y0} width={Math.max(b - a, 3)} height={9} rx={2} fill={col} fillOpacity={0.85} stroke={col} strokeWidth={0.8} strokeDasharray={sp.approx ? '3 2' : undefined} />
                {lx + 30 < X1 && <text x={lx} y={y0 + 8} fontSize={9.5} fill={C.text}>
                  {first && <tspan fontWeight={700}>{text}  </tspan>}
                  <tspan fontWeight={700} fill={col}>{type === 'del' ? 'del' : 'dup'}</tspan>
                  <tspan dx={4} fontWeight={600}>{cls}</tspan></text>}
              </g>
            );
          };
          const hasDel = !!p.types.del && p.r.rows.some(r2 => r2.type === 'del');
          return <g key={`achro${p.r.id}`}>{bar('del', yy + 1, true)}{bar('dup', yy + 14, !hasDel)}</g>;
        })}
        {yAchro >= 0 && !achroItems.length && <text x={X0 + 6} y={yAchro + 10} fontSize={9.5} fill={C.faint}>{ACHRO_REGIONS.some(r => achroChrom(r) === chrom) ? 'no AChro-Puce recurrent CNV in view' : `no AChro-Puce sheet on ${chrom}`}</text>}
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
          // the ribbon joins each reference end to that end's place on the rebuilt chromosome: crossed when inverted
          const da = cx(derX(d, toDer(dp, lo))), db = cx(derX(d, toDer(dp, hi)));
          if (Math.max(ra, rb, da, db) < X0 - 2 || Math.min(ra, rb, da, db) > X1 + 2) return null;
          const y0 = yRefBand + BAND_H, y1 = derY[di].band, m = (y0 + y1) / 2;
          return <path key={`rib${di}-${dp.index}`} d={`M${ra},${y0} C${ra},${m} ${da},${m} ${da},${y1} L${db},${y1} C${db},${m} ${rb},${m} ${rb},${y0} Z`} fill={col} fillOpacity={0.16} stroke={col} strokeOpacity={0.45} strokeWidth={0.9}><title>{`${dp.seg.letter}${dp.piece.reversed ? '′ (inverted)' : ''}: ${chrom}:${fmt(dp.seg.start + 1)}-${fmt(dp.seg.end)} → ${d.name}${dp.piece.copy ? ` (copy ${dp.piece.copy + 1})` : ''}`}</title></path>;
        }))}
        {/* rebuilt chromosomes */}
        {ders.map((d, di) => {
          const Y = derY[di], packed = derPacked[di], dd = derDensity[di], yDd = Y.genes + packed.height + 2;
          return (
            <g key={`der${d.name}`}>
              {dd && <path d={dd.d} transform={`translate(0 ${yDd + 34}) scale(1 34)`} fill="url(#svb-dens)" stroke="none"><title>{`genes per ${kb(3 / ppb)} on ${d.name}: up to ${dd.max}`}</title></path>}
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
                    {b.labeled && !isFusion && badge(b.genes[0], lx + (dot ? 10 : 0) + textWidth(b.label) + 4, yy + 9.5)}
                    {b.labeled && isFusion && b.fusion && <g>
                      <rect x={lx + labelW(b.label) - 4} y={yy + 1} width={b.fusion.kind === 'in-frame' ? 52 : b.fusion.kind === 'promoter-swap' ? 74 : 66} height={12} rx={6} fill={b.fusion.kind === 'in-frame' ? '#7c3aed' : '#ede9fe'} />
                      <text x={lx + labelW(b.label) + 2} y={yy + 10} fontSize={9} fontWeight={700} fill={b.fusion.kind === 'in-frame' ? '#fff' : '#6d28d9'}>{b.fusion.kind === 'in-frame' ? 'in frame' : b.fusion.kind === 'promoter-swap' ? 'promoter swap' : b.fusion.kind === 'out-of-frame' ? 'out of frame' : b.fusion.kind === 'five-cds-whole' ? '3′ UTR join' : 'frame ?'}</text>
                    </g>}
                    {lod === 'density' && <line x1={lx} x2={lx} y1={yy + 12} y2={yDd + 6} stroke={head} strokeWidth={0.8} strokeOpacity={0.6} />}
                    {lod !== 'density' && b.flank && (() => {
                      const { g, dp } = b.flank, X = (pos: number) => derX(d, toDer(dp, pos));
                      return lod === 'full' ? glyph(g.principal, X, yl, head, `fl${g.symbol}`, g.principal.start, g.principal.end, (dp.piece.reversed ? -g.strand : g.strand) as 1 | -1)
                        : <rect x={cx(Math.min(X(g.start), X(g.end)))} y={yy + LABEL_H + 3} width={Math.max(1.5, Math.abs(cx(X(g.end)) - cx(X(g.start))))} height={6} rx={1.5} fill={head} />;
                    })()}
                    {lod !== 'density' && b.chains.map((ch, ci) => {
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
      {/* the chromosomes' ends (qter) in view: their lengths, and the difference between them */}
      {rec && ders.map((d, di) => {
        const ext = derExtent(d), xr = refX(L), xd = derX(d, ext.hi), diff = ext.len - L;
        const refEnd = xr <= X1 + 2, derEnd = xd <= X1 + 2, yR = yRefBand, yD = derY[di].band;
        return (
          <g key={`ends${di}`} pointerEvents="none">
            {refEnd && <text x={Math.min(xr, X1)} y={yR - 3} fontSize={9.5} fontWeight={700} textAnchor="end" fill={C.muted}>qter · {kb(L)}</text>}
            {derEnd && <text x={Math.min(xd, X1)} y={yD + BAND_H + 21} fontSize={9.5} fontWeight={700} textAnchor="end" fill={C.accent}>qter · {kb(ext.len)}</text>}
            {refEnd && derEnd && Math.abs(xd - xr) > 6 && (() => {
              const a = Math.min(xr, xd), b = Math.max(xr, xd), yb = diff > 0 ? yD - 9 : yR + BAND_H + 9, col = diff > 0 ? '#15803d' : C.lost;
              return <g>
                <line x1={xr} x2={xr} y1={yR + BAND_H} y2={yD} stroke={C.muted} strokeWidth={0.8} strokeDasharray="3 3" />
                <path d={`M${a},${yb - 4}v4H${b}v-4`} fill="none" stroke={col} strokeWidth={1.4} transform={diff > 0 ? undefined : `translate(0 0) scale(1 1)`} />
                <text x={b - a < 170 ? b : (a + b) / 2} y={yb - 6} fontSize={10} fontWeight={700} textAnchor={b - a < 170 ? 'end' : 'middle'} fill={col}>{diff > 0 ? `${d.name} ${kb(diff)} longer` : `${d.name} ${kb(-diff)} shorter`}</text>
              </g>;
            })()}
          </g>
        );
      })}
      {/* the gene whose card is open, joined to each of its copies on the rebuilt chromosome */}
      {linkSym && (() => {
        const ri = refItems.items.find(it => it.g.symbol === linkSym);
        const ref = ri ? { x0: Math.max(X0, ri.x0) - 4, x1: Math.min(X1, Math.max(ri.xg, ri.x0 + ri.lw)) + 4, y: yRefGenes + refItems.laneY[ri.lane] - 2, h: ri.h + 4 } : null;
        const copies = ders.flatMap((d, di) => derPacked[di].items.filter(b => b.genes.includes(linkSym)).map(b => ({ x0: Math.max(X0, b.x0) - 4, x1: Math.min(X1, Math.max(b.xg, b.x0 + b.lw)) + 4, y: derY[di].genes + derPacked[di].laneY[b.lane] - 2, h: b.h + 4, xg: (Math.max(X0, b.x0) + Math.min(X1, b.xg)) / 2 })));
        const rx = ri ? (Math.max(X0, ri.x0) + Math.min(X1, ri.xg)) / 2 : 0;
        return (
          <g pointerEvents="none" clipPath="url(#svb-clip)">
            {ref && <rect x={ref.x0} y={ref.y} width={ref.x1 - ref.x0} height={ref.h} rx={6} fill={C.accent} fillOpacity={0.08} stroke={C.accent} strokeWidth={1.6} />}
            {copies.map((k, i) => <rect key={`cp${i}`} x={k.x0} y={k.y} width={k.x1 - k.x0} height={k.h} rx={6} fill={C.accent} fillOpacity={0.08} stroke={C.accent} strokeWidth={1.6} />)}
            {ref && copies.map((k, i) => {
              const y0 = ref.y + ref.h, y1 = k.y, m = (y0 + y1) / 2;
              return <g key={`ln${i}`}>
                <path d={`M${rx},${y0} C${rx},${m} ${k.xg},${m} ${k.xg},${y1}`} fill="none" stroke="#fff" strokeWidth={4} strokeOpacity={0.85} />
                <path d={`M${rx},${y0} C${rx},${m} ${k.xg},${m} ${k.xg},${y1}`} fill="none" stroke={C.accent} strokeWidth={1.8} />
                <circle cx={k.xg} cy={y1} r={3} fill={C.accent} />
              </g>;
            })}
            {ref && <circle cx={rx} cy={ref.y + ref.h} r={3} fill={C.accent} />}
          </g>
        );
      })()}
      {domains.length > 0 && ders.map((d, i) => <text key={`tadl${i}`} x={14} y={derY[i].tad + 16} fontSize={10.5} fill={C.muted}>TADs</text>)}
      {/* the stretch being chosen on the ruler */}
      {sel && (() => {
        const a = Math.max(X0, refX(sel.s)), b = Math.min(X1, refX(sel.e));
        if (b <= a) return null;
        return <g pointerEvents="none" data-export="skip">
          <rect x={a} y={yRuler} width={Math.max(1, b - a)} height={H - yRuler - 6} fill="#f59e0b" fillOpacity={0.12} />
          {refX(sel.s) >= X0 && <line x1={a} x2={a} y1={yRuler} y2={H - 6} stroke="#d97706" strokeWidth={1.2} strokeDasharray="4 3" />}
          {refX(sel.e) <= X1 && <line x1={b} x2={b} y1={yRuler} y2={H - 6} stroke="#d97706" strokeWidth={1.2} strokeDasharray="4 3" />}
        </g>;
      })()}
      {/* the small variant searched: a line (a band when wider) through every track, its name on the ruler */}
      {variant && (() => {
        const a = refX(variant.start), b = refX(variant.end), mid = (a + b) / 2;
        if (b < X0 - 2 || a > X1 + 2) return null;
        const lw = Math.min(260, textWidth(variant.label) * 0.98 + 26), lx = Math.max(X0, Math.min(X1 - lw, mid - lw / 2)), py = yRuler + RULER_H / 2 - 8;
        const show = (e: React.MouseEvent) => { if (dragRef.current?.moved) return; const pt = boxPos(e); lastTipPos.current = pt; if (tip?.gene !== 'variant:') setTip({ ...pt, gene: 'variant:', variant: true }); else placeTip(pt.x, pt.y); };
        return <g>
          <g pointerEvents="none">
            {b - a > 2 && <rect x={Math.max(a, X0)} y={yRuler + RULER_H} width={Math.min(b, X1) - Math.max(a, X0)} height={H - yRuler - RULER_H - 6} fill={VARIANT_COLOR} fillOpacity={0.07} />}
            <line x1={mid} x2={mid} y1={py + 16} y2={H - 6} stroke={VARIANT_COLOR} strokeWidth={1.1} strokeDasharray="5 3" opacity={0.75} />
          </g>
          <g style={{ cursor: 'pointer' }} onMouseEnter={show} onMouseMove={show} onMouseLeave={hoverOff}
            onClick={e => { e.stopPropagation(); if (suppressClick.current) return; const pt = boxPos(e); setPinned(cur => (cur?.gene === 'variant:' ? null : { ...pt, gene: 'variant:', variant: true })); }}>
            <rect x={lx} y={py} width={lw} height={16} rx={8} fill={VARIANT_COLOR} />
            <path d={`M${Math.max(lx + 6, Math.min(lx + lw - 6, mid)) - 4},${py + 16}L${Math.max(lx + 6, Math.min(lx + lw - 6, mid)) + 4},${py + 16}L${Math.max(lx + 6, Math.min(lx + lw - 6, mid))},${py + 21}Z`} fill={VARIANT_COLOR} />
            <text x={lx + lw / 2} y={py + 11.5} fontSize={10} fontWeight={700} textAnchor="middle" fill="#fff" fontFamily={MONO}>{variant.label.length > 40 ? `${variant.label.slice(0, 39)}…` : variant.label}</text>
          </g>
        </g>;
      })()}
      {band && <g pointerEvents="none">
        <rect x={Math.min(band.a, band.b)} y={yRuler} width={Math.abs(band.b - band.a)} height={H - yRuler - 6} fill={C.accent} fillOpacity={0.08} stroke={C.accent} strokeDasharray="4 3" />
        <text x={(band.a + band.b) / 2} y={yRuler - 4} fontSize={10} fontWeight={700} textAnchor="middle" fill={C.accent}>{kb(Math.abs(band.b - band.a) / ppb)}{dragRef.current?.kind === 'select' ? ' · release to highlight' : ''}</text>
      </g>}
    </svg>
  );

  // ---- the card of the gene or population CNV under the pointer, or kept open ----
  useEffect(() => { const isGene = (t: TipAt | null) => !!t && !t.dgv && !t.achro && !t.clin && !t.variant; const sym = (isGene(pinned) ? pinned!.gene : null) ?? (isGene(tip) ? tip!.gene : null); if (sym) wantModes([sym]); }, [pinned?.gene, tip?.gene]); // eslint-disable-line react-hooks/exhaustive-deps
  const card = (at: TipAt, isPin: boolean) => {
    const place: React.CSSProperties = { left: Math.max(0, Math.min(at.x + 14, (boxRef.current?.clientWidth ?? 900) - 400)), ...(at.y > (boxRef.current?.clientHeight ?? 400) - 160 ? { bottom: (boxRef.current?.clientHeight ?? 400) - at.y + 14 } : { top: at.y + 14 }) };
    const ref = isPin ? (el: HTMLDivElement | null) => { if (el && boxRef.current) fitY(el, at.y, boxRef.current.clientHeight); } : (el: HTMLDivElement | null) => { tipEl.current = el; };
    if (at.achro) return <AchroCard key={isPin ? 'pin' : 'tip'} cardRef={ref} p={at.achro} type={at.achroType} chrom={chrom} cnv={cnv} pinned={isPin} place={isPin ? place : {}} onClose={() => setPinned(null)} onZoom={() => zoomTo(at.achro!.start, at.achro!.end)} />;
    if (at.clin) return <ClinvarCard key={isPin ? 'pin' : 'tip'} cardRef={ref} vs={at.clin} at={at.clinAt} range={at.clinRange} genes={at.clinGenes} nmdOf={nmdOf} chrom={chrom} cnv={cnv} pinned={isPin} place={isPin ? place : {}} onClose={() => setPinned(null)}
      onZoom={() => { if (at.clinGenes) { zoomTo(at.clinGenes.s, at.clinGenes.e); return; } const a = Math.min(...at.clin!.map(v => v.start)), b = Math.max(...at.clin!.map(v => v.end)); zoomTo(a, Math.max(b, a + 50)); }} />;
    if (at.variant && variant) {
      const pos = Math.floor((variant.start + variant.end) / 2);
      const inGenes = view.genes.filter(g => g.start <= pos && g.end > pos);
      const clin = clinSnvView.list.filter(v => v.start < variant.end && v.end > variant.start);
      return (
        <div key={isPin ? 'pin' : 'tip'} ref={ref} className={`absolute z-20 rounded-lg bg-slate-900 text-white text-[11px] px-3 py-2 shadow-xl w-max max-w-[460px] ${isPin ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={isPin ? place : {}}
          onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
          <div className="font-semibold flex items-start gap-2"><span className="flex-1 font-mono">{variant.label}</span>
            {isPin && <button onClick={() => setPinned(null)} className="-mr-1 -mt-0.5 px-1 text-slate-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
          {variant.lines.map((l, i) => <div key={i} className={i ? 'text-slate-400' : 'text-slate-200'}>{l}</div>)}
          {inGenes.length > 0
            ? inGenes.map(g => <div key={g.symbol} className="text-slate-200"><i>{g.symbol}</i>: {whereInTx(g.principal, pos)} <span className="text-slate-400">({g.principal.id}{g.principal.status === 'mane' ? ', MANE Select' : ''})</span></div>)
            : <div className="text-slate-400">{view.loading ? 'gene models loading…' : 'in no RefSeq gene'}</div>}
          <div className="mt-1 pt-1 border-t border-white/10">{clin.length
            ? <>ClinVar (pathogenic / likely pathogenic) here: {clin.slice(0, 4).map(v => <div key={`${v.id}:${v.start}`} className="text-slate-300">
                <span className="px-1 rounded text-[10px] font-bold text-white" style={{ background: v.code === 'PG' ? '#b91c1c' : '#f97316' }}>{v.code === 'PG' ? 'P' : 'LP'}</span> {v.start === variant.start && v.end === variant.end ? <b>same position</b> : 'overlapping'} · {v.conseq || v.type} <span className="font-mono text-[10px]">{v.protein ?? v.hgvs ?? ''}</span>
                {isPin && <> · <a href={clinvarUrl(v)} target="_blank" rel="noreferrer" className="text-accent-300 hover:underline">{v.id} ↗</a></>}</div>)}{clin.length > 4 ? <div className="text-slate-400">+{clin.length - 4} more</div> : null}</>
            : <span className="text-slate-400">{clinSnvView.loading ? 'ClinVar loading…' : 'no pathogenic or likely pathogenic ClinVar variant at this place'}</span>}</div>
          {isPin
            ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={() => zoomTo(variant.start, variant.end)} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on it</button></div>
            : <div className="mt-0.5 text-slate-400">click to keep this card open</div>}
        </div>
      );
    }
    if (at.dgv) return <DgvCard key={isPin ? 'pin' : 'tip'} cardRef={ref} v={at.dgv} chrom={chrom} cnv={cnv} pinned={isPin} place={isPin ? place : {}} onClose={() => setPinned(null)} onZoom={() => zoomTo(at.dgv!.outer[0], at.dgv!.outer[1])} />;
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
      <div key={isPin ? 'pin' : 'tip'} ref={ref} className={`absolute z-20 rounded-lg bg-slate-900 text-white text-[11px] px-3 py-2 shadow-xl w-max max-w-[390px] ${isPin ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={isPin ? place : undefined}
        onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
        <div className="font-semibold flex items-start gap-2"><span className="flex-1"><i>{g.symbol}</i>{gi ? ` · ${IMPACT[gi.kind].label}` : ''}</span>
          {isPin && <button onClick={() => setPinned(null)} className="-mr-1 -mt-0.5 px-1 text-slate-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
        <div className="text-slate-300">{tx.id}{tx.status === 'mane' ? ' (MANE Select)' : tx.status === 'select' ? ' (RefSeq Select)' : ''} · {chrom}:{fmt(g.start + 1)}-{fmt(g.end)} ({kb(g.end - g.start)}, {g.strand > 0 ? '+' : '−'} strand)</div>
        {gi && <div className="mt-0.5">{gi.principal.text}</div>}
        {(() => {
          const mo = modeOf(g.symbol);
          if (!mo.m && !mo.fromPanel && mo.cls === '?') return <div className="mt-0.5 text-slate-400">Inheritance (PanelApp): loading…</div>;
          const parts = mo.m?.parts.filter(p => p.cls !== '?' || mo.m!.parts.length === 1) ?? [];
          return <>
            <div className="mt-0.5"><span style={{ color: MOI[mo.cls].color === '#cbd5e1' ? '#94a3b8' : MOI[mo.cls].color }}>■</span> Inheritance: <b>{MOI[mo.cls].label}</b>
              <span className="text-slate-400">{mo.fromPanel ? ` (in ${mo.fromPanel})` : ''}{parts.length ? ` · ${mo.m!.source}: ${parts.map(p => `${MOI[p.cls].short} in ${p.panels} panel${p.panels > 1 ? 's' : ''}${p.green ? ` (${p.green} green)` : ''}`).join(', ')}` : mo.cls === 'none' ? ` · ${moiSrc.label}` : ''}</span></div>
            {(() => { const ds = diseasesOf(g.symbol); if (!ds.length) return null; const more = ds.length - 4;
              return <div className="mt-0.5" title={ds.map(d => `${d.name}${d.panels ? ` (${d.panels} panel${d.panels > 1 ? 's' : ''}${d.green ? `, ${d.green} green` : ''})` : ' (the loaded panel)'}`).join('\n')}>
                <span className="text-rose-300">✚</span> Disease (PanelApp): <b>{ds.slice(0, 4).map(d => d.name).join(' · ')}</b>{more > 0 ? <span className="text-slate-400"> +{more} more</span> : null}</div>; })()}
          </>;
        })()}
        {gi && ders.length > 0 && (() => { const cp = copiesOf(gi); return <div className="mt-0.5 text-teal-200">On {ders.map(d => d.name).join(', ')}: {cp.length ? `${cp.length} cop${cp.length > 1 ? 'ies' : 'y'} (${cp.join('; ')})` : 'none'}{isPin && cp.length ? ', joined in the figure (the links stay when this card is closed; click the background to clear them)' : ''}</div>; })()}
        {r?.text && <div className="mt-0.5 text-amber-200">{r.text}</div>}
        <div className="mt-0.5 text-slate-300"><span style={{ color: dot }}>●</span> {parts.join(' · ')}</div>
        {isPin
          ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={() => zoomTo(g.start, g.end)} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on <i>{g.symbol}</i></button></div>
          : <div className="mt-0.5 text-slate-400">click to keep this card open</div>}
      </div>
    );
  };

  // the Sashimi viewer's zoom group: one bordered strip of buttons
  const btn = 'min-w-[28px] h-[26px] grid place-items-center text-slate-600 hover:bg-slate-50';
  return (
    <div className="mt-2 border-t border-slate-200 pt-3 outline-none" style={{ fontFamily: FONT }} tabIndex={0} onKeyDown={onKeyDown}>
      {/* header */}
      <div className="px-5 flex flex-wrap items-center gap-2 text-xs">
        <div className="text-sm font-bold text-slate-900 mr-1">{rec ? 'Genes on the rebuilt chromosome' : 'Genes'}</div>
        <Segmented size="sm" label="Transcripts" title="The transcripts drawn" value={mode} onChange={setMode}
          options={(['mane', 'all'] as const).map(m => ({ value: m, label: m === 'mane' ? 'MANE Select' : 'all NM_ transcripts', hint: m === 'mane' ? 'One model per gene: MANE Select (else RefSeq Select, else the longest coding sequence)' : `Every RefSeq NM_ model of each gene (no XM_ prediction, no NR_ non-coding): the consequence can differ between isoforms. Drawn below ${kb(LOD_FULL)} in view.` }))} />
        <Segmented size="sm" label="Gene colour" title="What the genes are coloured by" value={colorBy} onChange={m => { setColorBy(m); store.set('sashimi.svColorBy', m); }}
          options={(['loeuf', 'impact', 'moi'] as const).map(m => ({ value: m, label: m === 'loeuf' ? 'colour: LOEUF' : m === 'moi' ? 'colour: inheritance' : 'colour: impact', hint: m === 'loeuf' ? `Genes coloured by their gnomAD LoF constraint (LOEUF): below ${LOEUF_INTOLERANT.v4} red (Mendelian threshold), ${LOEUF_INTOLERANT.v4}–${LOEUF_CONSTRAINED.v4} orange (gnomAD v4 general cut-off), above grey; v2.1.1 scores: below ${LOEUF_INTOLERANT['v2.1.1']} red` : m === 'moi' ? `Genes coloured by their mode of inheritance in PanelApp (${moiSrc.label}): the loaded panel's own value, else the commonest among the gene's panels, green ratings first` : 'Genes coloured by what the rearrangement does to them' }))} />
        <PanelPicker panel={panelData} onPanel={p => { setPanelData(p); store.set('sashimi.svPanel', p); }} />
        <GeneFilterControl filter={geneFilter} onChange={setGeneFilter} hasPanel={!!panel} panelName={panel?.name ?? null} hasEvent={!!rec}
          counts={(() => { const inV = visGenes.filter(g => g.end > vp.s && g.start < ve); const shownN = filtering ? inV.filter(g => passes(g.symbol)).length : inV.length; const waiting = geneFilter.moi.length ? inV.filter(g => { const m = modes.current.get(g.symbol); return m === undefined || m === 'loading'; }).length : 0; return { shown: shownN, total: inV.length, waiting }; })()} />
        {loaders.getPopulationSvs && (
          <Pill size="sm" icon="users" label="DGV Gold" on={showDgv} onChange={v => { setShowDgv(v); store.set('sashimi.svDgv', v); }}
            title="DGV Gold Standard (UCSC track dgvGold): curated copy-number variants of the general population, each found in two studies and two samples at least; drawn under the genes and compared with the event's losses and gains" />
        )}
        {loaders.getClinvar && <Pill size="sm" icon="star" label="ClinVar SV" on={showClinSv} onChange={v => { setShowClinSv(v); store.set('sashimi.svClinvarSv', v); }}
          title="ClinVar variants of 50 bp or more classified pathogenic or likely pathogenic (UCSC track clinvarCnv): losses red, gains blue, inversions and complex purple; pathogenic solid, likely pathogenic paler" />}
        {loaders.getClinvar && <Pill size="sm" icon="marker" label="ClinVar SNV" on={showClinSnv} onChange={v => { setShowClinSnv(v); store.set('sashimi.svClinvarSnv', v); }}
          title={`ClinVar SNVs and short indels (under 50 bp) classified pathogenic or likely pathogenic, UCSC track clinvarMain: counted per gene (loss of function by predicted NMD, missense) or each at its place; shown below ${kb(CLIN_SNV_MAX_SPAN)} in view`} />}
        <Pill size="sm" icon="grid" label="AChro-Puce" on={showAchro} onChange={v => { setShowAchro(v); store.set('sashimi.svAchro', v); }}
          title={`${ACHRO_SOURCE.title} (${ACHRO_SOURCE.version}): the recurrent CNVs the French network classifies, each with its class for the deletion and the duplication`} />
        {rec && <label className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium cursor-pointer whitespace-nowrap ${domainsFile ? 'bg-accent-50 border-accent-200 text-accent-600' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'}`} title="Chromatin domains (TADs) as a BED file (chrom, start, end, name; 0-based, the same genome build), e.g. from Hi-C of a relevant tissue. They show which genes change domain: an intact gene joined to another domain may be misregulated (position effect).">
          <Icon name="upload" size={13} />{domainsFile ? `TADs: ${domainsFile.name} (${domains.length})` : 'Add chromatin domains (BED)…'}
          <input type="file" accept=".bed,.txt,.tsv" className="hidden" onChange={async e => { const f = e.target.files?.[0]; if (f) { const v = { name: f.name, text: await f.text() }; setDomainsFile(v); store.set('sashimi.svDomains', v.text.length < 3e6 ? v : null); } e.target.value = ''; }} />
        </label>}
        {rec && domainsFile && <button onClick={() => { setDomainsFile(null); store.set('sashimi.svDomains', null); }} className="text-slate-400 hover:text-red-600" title="Remove the domains">×</button>}
        {rec && <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-xs text-slate-500" title={`Without chromatin domains, an intact gene within this distance of a breakpoint is reported as a possible position effect. TADs average about 1 Mb; position effects have been documented 1.3 Mb from a breakpoint (SOX9) and further.${domains.length ? ' With domains loaded, the domains are used instead.' : ''}`}>
          near a breakpoint ≤
          <select value={near} onChange={e => { const v = Number(e.target.value); setNear(v); store.set('sashimi.svNear', v); }} className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            {[0, 250e3, 500e3, 1e6, 2e6, 3e6].map(v => <option key={v} value={v}>{v === 0 ? 'off' : v >= 1e6 ? `${v / 1e6} Mb` : `${v / 1e3} kb`}</option>)}
          </select>
        </label>}
        <span className="ml-auto flex items-center gap-1">
          <span className="font-mono text-[11px] text-slate-600 mr-1 tabular-nums" title="The stretch in view">{chrom}:{fmt(vp.s + 1)}-{fmt(vp.e)} <span className="text-slate-400">· {kb(span)}</span></span>
          <span className="inline-flex border border-slate-200 rounded-[6px] overflow-hidden bg-white" role="group" aria-label="Zoom">
            <button className={btn} onClick={() => flyTo(vp.s - span * 0.5, vp.e + span * 0.5, 300)} title="Zoom out (−, Ctrl + wheel, Shift + double-click)" aria-label="Zoom out"><Icon name="minus" size={14} /></button>
            <button className={`${btn} border-l border-slate-200`} onClick={() => flyTo(vp.s + span * 0.25, vp.e - span * 0.25, 300)} title="Zoom in (+, Ctrl + wheel, double-click)" aria-label="Zoom in"><Icon name="plus" size={14} /></button>
            {rec && <button className={`${btn} border-l border-slate-200 px-2 text-xs font-medium`} onClick={eventView} title="Fit the event (E)">event</button>}
            <button className={`${btn} border-l border-slate-200 px-2 text-xs font-medium`} onClick={() => { const w = wholeView(); flyTo(w.s, w.e, 700); }} title="The whole chromosome (and the rebuilt one's ends)">{chrom.replace(/^chr/, 'chr ')}</button>
          </span>
        </span>
      </div>
      {rec && showDgv && loaders.getPopulationSvs && <DgvSummary cnv={cnv} loading={pop.loading} error={pop.error} />}
      {rec && showAchro && <AchroSummary cnv={cnv} placed={achroPlaced} />}
      <div className="px-3 mt-2"><div ref={widthBox} className="w-full" /></div>
      <div className="px-3 mt-0 relative" ref={boxRef}>
        {svg}
        {sel && (() => {
          // the bar over the highlighted interval, in the box's pixels (centred on its visible part)
          const r = svgRef.current?.getBoundingClientRect(), br = boxRef.current?.getBoundingClientRect();
          if (!r || !br) return null;
          const a = Math.max(X0, refX(sel.s)), b = Math.min(X1, refX(sel.e)), k = r.width / VW;
          const inView = b > a, mid = inView ? (a + b) / 2 : refX(sel.s) < X0 ? X0 : X1;
          const left = r.left - br.left + mid * k, top = r.top - br.top + (yRuler + RULER_H / 2) * k;
          return (
            <div className="absolute z-30 -translate-x-1/2 -translate-y-1/2 flex items-center gap-1 rounded-md border border-amber-300 bg-white/95 shadow-md px-1.5 py-0.5 text-[11px] whitespace-nowrap"
              style={{ left: Math.max(170, Math.min(br.width - 170, left)), top }} onMouseDown={e => e.stopPropagation()} onPointerDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
              {!inView && refX(sel.e) <= X0 ? <span className="text-slate-400">◂</span> : null}
              <span className="font-mono text-slate-700" title="the highlighted interval (1-based, inclusive)">{selText}</span>
              <span className="text-slate-400">· {kb(sel.e - sel.s)}</span>
              <button onClick={zoomToSel} className="ml-1 px-1.5 py-0.5 rounded bg-amber-500 text-white font-semibold hover:bg-amber-600" title="Zoom to the highlighted interval">🔍 zoom to region</button>
              <button onClick={() => { void navigator.clipboard?.writeText(selText).then(() => setSelCopied(true), () => {}); }} className="px-1.5 py-0.5 rounded border border-slate-300 text-slate-600 hover:bg-slate-50" title="Copy the coordinates">{selCopied ? '✓ copied' : 'copy'}</button>
              <button onClick={() => setSel(null)} className="px-1 text-slate-400 hover:text-slate-700 text-sm leading-none" title="Remove the highlight (Esc)">×</button>
              {!inView && refX(sel.s) >= X1 ? <span className="text-slate-400">▸</span> : null}
            </div>
          );
        })()}
        {pinned && card(pinned, true)}
        {tip && tip.gene !== pinned?.gene && card(tip, false)}
      </div>
      {/* legend */}
      <div className="px-5 mt-1 flex flex-wrap items-center gap-x-3 gap-y-1 text-[10.5px] text-slate-500">
        {colorBy === 'moi' ? (Object.keys(MOI) as (MoiClass | 'none')[]).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: MOI[k2].color }} />{MOI[k2].label}</span>)
          : colorBy === 'impact'
          ? (['fusion', 'disrupted', 'deleted', 'intragenic-del', 'intragenic-dup', 'partial-copy', 'duplicated', 'domain', 'moved', 'nearby', 'unaffected'] as ImpactKind[]).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: IMPACT[k2].color }} />{IMPACT[k2].label}</span>)
          : (['intolerant', 'constrained', 'tolerant', 'unknown'] as const).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: LOEUF_CLASS[k2].color }} />{LOEUF_CLASS[k2].label('v4')}</span>)}
        {colorBy === 'loeuf' && <span>(gnomAD v4; v2.1.1 scores: &lt; {LOEUF_INTOLERANT['v2.1.1']} intolerant)</span>}
        {rec && <span className="inline-flex items-center gap-1"><svg width="10" height="14"><path d="M6,1 l-3,3 l3,3 l-3,3 l3,3" fill="none" stroke={C.lost} strokeWidth="1.4" /></svg>cut by a breakpoint</span>}
        {showAchro && <span className="inline-flex items-center gap-1">AChro-Puce: <span className="inline-block w-4 h-2.5 rounded-sm" style={{ background: ACHRO_TYPE_COLOR.del }} /> deletion <span className="inline-block w-4 h-2.5 rounded-sm" style={{ background: ACHRO_TYPE_COLOR.dup }} /> duplication, the network's class after each bar (click a bar for its sheet); dashed: placed approximately</span>}
        {(showClinSv || showClinSnv) && loaders.getClinvar && <span className="inline-flex items-center gap-1">ClinVar (pathogenic / likely pathogenic):
          {showClinSv && !clinSvAll && <> SVs summed up (more than {CLIN_SV_STACK_MAX} overlap here): how many cover each point, <span style={{ color: CLIN_KIND_COLOR.loss }}>losses</span> above the line, <span style={{ color: CLIN_KIND_COLOR.gain }}>gains</span> below, <span style={{ color: CLIN_KIND_COLOR.inv }}>inversions and complex</span> on it{rec ? ', the ones matching a loss or gain of the event as outlined bars under it' : ''} (“▸ all SVs” lists them);</>}
          {showClinSv && clinSvAll && <><span className="inline-block w-4 h-2 rounded-sm" style={{ background: CLIN_KIND_COLOR.loss }} /> loss <span className="inline-block w-4 h-2 rounded-sm" style={{ background: CLIN_KIND_COLOR.gain }} /> gain <span className="inline-block w-4 h-2 rounded-sm" style={{ background: CLIN_KIND_COLOR.inv }} /> inversion, complex (pale: likely pathogenic){rec ? '; outlined when it matches a loss or gain of the event' : ''};</>}
          {showClinSnv && <span> SNV{clinSnvMode === 'gene' ? ' counted per gene (neighbours merged when crowded, ×n genes)' : clinSnvSum ? ` (summed per ${kb(CLIN_SNV_BIN * bppRef.current)}, zoom in below ${kb(CLIN_SNV_SUM_SPAN)} for each)` : ''}: <span style={{ color: CLIN_EFFECT_COLOR.lof }}>● loss of function predicted to trigger NMD</span> and <span style={{ color: CLIN_LOF_ESC }}>● predicted to escape it</span> above the line (nonsense, frameshift, canonical splice site, start lost; on the MANE Select transcript, a stop more than {NMD_BOUNDARY_NT} nt before the last exon junction triggers NMD, a splice site is assumed to skip its exon), <span style={{ color: CLIN_EFFECT_COLOR.missense }}>▼ missense</span> below{clinSnvMode === 'gene' ? '; pathogenic and likely pathogenic counted together' : <>, <span style={{ color: CLIN_EFFECT_COLOR.inframe }}>◆ in frame</span> and <span style={{ color: CLIN_EFFECT_COLOR.other }}>other</span> on it; filled pathogenic, hollow likely pathogenic</>}</span>}</span>}
        {showDgv && loaders.getPopulationSvs && <span className="inline-flex items-center gap-1"><svg width="26" height="10"><line x1="1" x2="25" y1="5" y2="5" stroke={DGV_COLOR.loss} /><line x1="1" x2="1" y1="2" y2="8" stroke={DGV_COLOR.loss} /><line x1="25" x2="25" y1="2" y2="8" stroke={DGV_COLOR.loss} /><rect x="8" y="2.5" width="10" height="5" fill={DGV_COLOR.loss} /></svg>DGV Gold CNV: outer extent (line), inner (bar); <span style={{ color: DGV_COLOR.loss }}>loss</span> <span style={{ color: DGV_COLOR.gain }}>gain</span> <span style={{ color: DGV_COLOR.both }}>both</span>{rec ? '; solid and outlined when it matches a loss or gain of the event' : ''}</span>}
        <span className="inline-flex items-center gap-1">badge after a name: mode of inheritance (PanelApp) — {(['AD', 'AR', 'AD+AR', 'XLR', 'XLD', 'MT'] as MoiClass[]).map(k2 => <span key={k2} className="inline-flex items-center rounded-full px-1 text-[9px] font-bold border" style={{ color: MOI[k2].color, borderColor: `${MOI[k2].color}88`, background: `${MOI[k2].color}21` }} title={MOI[k2].label}>{MOI[k2].short}</span>)}</span>
        {panel && <span className="inline-flex items-center gap-1">● panel rating: <span style={{ color: CONF[3].color }}>green</span> <span style={{ color: CONF[2].color }}>amber</span> <span style={{ color: CONF[1].color }}>red</span></span>}
        <span>· drag to pan (it glides on), Ctrl / ⌘ + scroll or pinch to zoom, two-finger swipe or Shift + scroll to pan, double-click to zoom in (Shift: out), Ctrl / ⌘ + drag to highlight an interval (then zoom to it), drag along the ruler to zoom on a stretch, click the ideogram to go there; ← → + − {rec ? 'E ' : ''}on the keyboard; hover a gene for its card, click to keep it open</span>
      </div>
      <div className="px-5 mt-2 pb-1">
        {rec && !story.genes.some(g => g.kind !== 'unaffected') && !txs.loading && <div className="text-slate-400 text-xs py-1">No gene of the window is changed by the selected junctions.</div>}
        <div className="text-[10.5px] text-slate-500 mt-2 leading-4">
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

// ======================== Gene filter ========================

/** The gene filter: a button with the number of criteria set, opening the criteria; every criterion set must hold. */
function GeneFilterControl({ filter, onChange, hasPanel, panelName, hasEvent, counts }: { filter: GeneFilter; onChange: (f: GeneFilter) => void; hasPanel: boolean; panelName: string | null; hasEvent: boolean; counts: { shown: number; total: number; waiting: number } }) {
  const [open, setOpen] = useState(false);
  const [custom, setCustom] = useState<string>(filter.loeuf != null && ![0.35, 0.45, 0.6].includes(filter.loeuf) ? String(filter.loeuf) : '');
  const n = filterCount(filter, hasPanel, hasEvent);
  const seg = <T extends string | number | null>(value: T, options: { v: T; label: string; title?: string }[], set: (v: T) => void, disabled = false) => (
    <span className={`inline-flex rounded-md border border-slate-300 overflow-hidden ${disabled ? 'opacity-40 pointer-events-none' : ''}`}>
      {options.map(o => <button key={String(o.v)} onClick={() => set(o.v)} title={o.title} className={`px-2 py-0.5 whitespace-nowrap ${value === o.v ? 'bg-accent-600 text-white' : 'bg-white text-slate-600 hover:bg-slate-50'}`}>{o.label}</button>)}
    </span>
  );
  const MOIS: (MoiClass | 'none')[] = ['AD', 'AR', 'AD+AR', 'XLR', 'XLD', 'MT', '?', 'none'];
  return (
    <span className="relative">
      <button onClick={() => setOpen(o => !o)} className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium whitespace-nowrap ${n ? 'bg-accent-50 border-accent-200 text-accent-600' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'}`}
        title="Show only some genes: in the panel, of some modes of inheritance, under a LOEUF, changed by the SV">
        <Icon name="filter" size={13} />Genes{n ? <span className="text-[10.5px] font-bold leading-4 px-1.5 rounded-full bg-accent-600 text-white">{n}</span> : null}{n ? <span className="text-slate-500">· {counts.shown}/{counts.total}</span> : null}
      </button>
      {open && (
        <div className="absolute z-40 mt-1 w-[470px] rounded-lg border border-slate-300 bg-white shadow-2xl p-3 text-xs space-y-2.5" onMouseDown={e => e.stopPropagation()}>
          <div className="flex items-center justify-between"><b>Show the genes</b><span className="flex items-center gap-2">
            <button onClick={() => { onChange(NO_FILTER); setCustom(''); }} className="text-accent-700 hover:underline">all genes</button>
            <button onClick={() => setOpen(false)} className="px-1.5 text-slate-400 hover:text-slate-700">×</button></span></div>
          <div className="flex items-center gap-2"><span className="w-36 text-slate-600">in the gene panel</span>
            {seg(filter.panel, [{ v: 'any', label: 'any' }, { v: 'in', label: 'in the panel' }, { v: 'green', label: 'green only' }], v => onChange({ ...filter, panel: v }), !hasPanel)}
            {!hasPanel && <span className="text-slate-400">load a panel first</span>}{hasPanel && <span className="text-slate-400 truncate max-w-[110px]" title={panelName ?? ''}>{panelName}</span>}</div>
          <div className="flex items-start gap-2"><span className="w-36 text-slate-600 pt-0.5">mode of inheritance</span>
            <span className="flex flex-wrap gap-1 flex-1">
              {MOIS.map(m => { const on = filter.moi.includes(m); return (
                <button key={m} onClick={() => onChange({ ...filter, moi: on ? filter.moi.filter(x => x !== m) : [...filter.moi, m] })} title={MOI[m].label}
                  className="px-1.5 py-0.5 rounded-full border text-[10.5px] font-bold" style={on ? { background: MOI[m].color, color: '#fff', borderColor: MOI[m].color } : { color: MOI[m].color === '#cbd5e1' ? '#64748b' : MOI[m].color, borderColor: `${MOI[m].color}88` }}>
                  {m === 'none' ? 'not in PanelApp' : m === '?' ? 'unknown' : MOI[m].short}</button>); })}
              <span className="text-slate-400 w-full">{filter.moi.length ? 'genes of any of these modes' : 'any mode'}</span>
            </span></div>
          <div className="flex items-center gap-2"><span className="w-36 text-slate-600">LOEUF below</span>
            {seg(filter.loeuf, [{ v: null, label: 'off' }, { v: 0.35, label: '0.35', title: 'gnomAD v2.1.1 LoF-intolerance cut-off' }, { v: 0.45, label: '0.45', title: 'gnomAD v4 browser threshold for Mendelian disorders' }, { v: 0.6, label: '0.6', title: 'gnomAD v4.0 general cut-off' }], v => { setCustom(''); onChange({ ...filter, loeuf: v }); })}
            <input value={custom} onChange={e => { setCustom(e.target.value); const v = Number(e.target.value.replace(',', '.')); if (e.target.value.trim() && Number.isFinite(v) && v > 0) onChange({ ...filter, loeuf: v }); }} placeholder="other" className="w-14 border rounded px-1 py-0.5" /></div>
          <div className="flex items-center gap-2"><span className="w-36 text-slate-600">changed by the SV</span>
            {seg(filter.changed, [{ v: 'any', label: 'any' }, { v: 'changed', label: 'changed' }, { v: 'near', label: 'changed or near a breakpoint' }], v => onChange({ ...filter, changed: v }), !hasEvent)}
            {!hasEvent && <span className="text-slate-400">no SV drawn</span>}</div>
          <div className="text-slate-500 border-t border-slate-100 pt-1.5">{counts.shown} of {counts.total} genes in view shown{counts.waiting ? ` · ${counts.waiting} waiting for their PanelApp modes` : ''}. Every criterion set must hold; a gene without the data a criterion needs (no LOEUF, not yet looked up) is hidden.</div>
        </div>
      )}
    </span>
  );
}

// ======================== AChro-Puce ========================

const typeWord = (r: AchroRow) => `${r.type === 'del' ? 'Deletion' : 'Duplication'}${r.sex ? (r.sex === 'male' ? ' (male)' : ' (female, carrier)') : ''}`;
/** reciprocal overlap of an event's loss or gain with a region */
/** reciprocal overlap of an event's loss or gain with a ClinVar variant */
const clinOverlap = (ev: CnvEvent, v: ClinvarVariant) => Math.max(0, Math.min(ev.end, v.end) - Math.max(ev.start, v.start)) / Math.max(ev.end - ev.start, v.end - v.start, 1);
const CLIN_KIND_WORD: Record<string, string> = { loss: 'deletion / loss', gain: 'duplication / gain', inv: 'inversion', complex: 'complex', snv: 'SNV', indel: 'indel', other: 'variant' };
/** The card of ClinVar variants (one SV, or the SNVs under one lollipop): significance, type, place, genes, conditions, review, the event's overlap. */
function ClinvarCard({ vs, at, range, genes, nmdOf, chrom, cnv, pinned, place, onClose, onZoom, cardRef }: { vs: ClinvarVariant[]; at?: number; range?: [number, number]; genes?: { names: string[]; tx: string[]; s: number; e: number }; nmdOf?: (v: ClinvarVariant) => NmdCall | undefined; chrom: string; cnv: { ev: CnvEvent }[]; pinned: boolean; place: React.CSSProperties; onClose: () => void; onZoom: () => void; cardRef?: (el: HTMLDivElement | null) => void }) {
  // over a point of the summary: the smallest first (the most specific); else pathogenic first, then the best reviewed
  const shown = (at != null || range || genes ? [...vs] : [...vs].sort((a, b) => (a.code === b.code ? b.stars - a.stars : a.code === 'PG' ? -1 : 1))).slice(0, genes ? 8 : 6);
  const nOf = (k: (v: ClinvarVariant) => boolean) => vs.filter(k).length;
  const lofN = (call: 'nmd' | 'escape') => nOf(v => v.effect === 'lof' && nmdOf?.(v)?.call === call), lofA = (call: 'nmd' | 'escape') => nOf(v => v.effect === 'lof' && nmdOf?.(v)?.call === call && !!nmdOf?.(v)?.assumed);
  const split = (k: (v: ClinvarVariant) => boolean) => { const pg = nOf(v => k(v) && v.code === 'PG'), lp = nOf(v => k(v) && v.code === 'LP'); return [pg && `${pg} P`, lp && `${lp} LP`].filter(Boolean).join(', '); };
  return (
    <div ref={cardRef} className={`absolute z-20 rounded-lg bg-slate-900 text-white text-[11px] px-3 py-2 shadow-xl w-max max-w-[460px] ${pinned ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={place}
      onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
      <div className="font-semibold flex items-start gap-2"><span className="flex-1">{genes ? <><i>{genes.names.length ? genes.names.slice(0, 6).join(', ') + (genes.names.length > 6 ? ` +${genes.names.length - 6}` : '') : 'outside genes'}</i> · ClinVar SNVs</> : 'ClinVar'}{genes ? '' : range ? ` · ${vs.length} variant${vs.length > 1 ? 's' : ''} in ${chrom}:${fmt(range[0] + 1)}-${fmt(range[1])}` : at != null ? ` · ${vs.length} SV${vs.length > 1 ? 's' : ''} over ${chrom}:${fmt(at + 1)}` : vs.length > 1 ? ` · ${vs.length} variants here` : ''} <span className="font-normal text-slate-400">· pathogenic / likely pathogenic</span></span>
        {pinned && <button onClick={onClose} className="-mr-1 -mt-0.5 px-1 text-slate-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
      {genes && <div className="text-slate-300 leading-snug">
        {lofN('nmd') > 0 && <div><span style={{ color: '#fb7185' }}>●</span> {lofN('nmd')} loss of function predicted to trigger NMD ({split(v => v.effect === 'lof' && nmdOf?.(v)?.call === 'nmd')}){lofA('nmd') ? `, ${lofA('nmd')} of them assumed` : ''}</div>}
        {lofN('escape') > 0 && <div><span style={{ color: '#fdba74' }}>●</span> {lofN('escape')} loss of function predicted to escape NMD ({split(v => v.effect === 'lof' && nmdOf?.(v)?.call === 'escape')}){lofA('escape') ? `, ${lofA('escape')} of them assumed` : ''}</div>}
        {nOf(v => v.effect === 'missense') > 0 && <div><span style={{ color: '#c4b5fd' }}>▼</span> {nOf(v => v.effect === 'missense')} missense ({split(v => v.effect === 'missense')})</div>}
        {nOf(v => v.effect !== 'lof' && v.effect !== 'missense') > 0 && <div><span className="text-slate-400">◆</span> {nOf(v => v.effect !== 'lof' && v.effect !== 'missense')} in frame or other ({split(v => v.effect !== 'lof' && v.effect !== 'missense')}), not drawn</div>}
        {genes.tx.length > 0 && <div className="text-slate-400">NMD predicted on {genes.tx.slice(0, 3).join(', ')}{genes.tx.length > 3 ? '…' : ''}: a stop more than {NMD_BOUNDARY_NT} nt before the last exon junction triggers it; splice sites assume the exon skipped</div>}
        {vs.length > shown.length && <div className="text-slate-400">the first {shown.length} below, loss of function first</div>}
      </div>}
      {range && <div className="text-slate-300">{(['lof', 'missense', 'inframe', 'other'] as const).map(e => [e, nOf(v => v.effect === e)] as const).filter(x => x[1]).map(([e, n]) => <span key={e} className="mr-2"><span style={{ color: CLIN_EFFECT_COLOR[e] === '#94a3b8' ? '#cbd5e1' : CLIN_EFFECT_COLOR[e] }}>■</span> {n} {CLIN_EFFECT_WORD[e]}</span>)}· {nOf(v => v.code === 'PG')} pathogenic, {nOf(v => v.code === 'LP')} likely pathogenic{vs.length > 6 ? ' · the first 6 (loss of function first)' : ''}</div>}
      {at != null && <div className="text-slate-300">{[['loss', nOf(v => v.kind === 'loss')], ['gain', nOf(v => v.kind === 'gain')], ['inversion or complex', nOf(v => v.kind !== 'loss' && v.kind !== 'gain')]].filter(x => x[1]).map(([w, n]) => `${n} ${w}${Number(n) > 1 ? (w === 'loss' ? 'es' : 's') : ''}`).join(' · ')} · {nOf(v => v.code === 'PG')} pathogenic, {nOf(v => v.code === 'LP')} likely pathogenic{vs.length > 6 ? ' · the 6 smallest below' : ''}</div>}
      {genes && <div className="mt-1 pt-1 border-t border-white/10">{shown.map(v => {
        const n = v.effect === 'lof' ? nmdOf?.(v) : undefined;
        return (
          <div key={`${v.id}:${v.start}`} className="truncate max-w-[440px]" title={[v.conditions, n?.why].filter(Boolean).join(' — ')}>
            <span className="px-1 rounded text-[10px] font-bold" style={{ background: v.code === 'PG' ? '#b91c1c' : '#f97316' }}>{v.code === 'PG' ? 'P' : 'LP'}</span>{' '}
            <span style={{ color: v.effect === 'lof' ? (n?.call === 'escape' ? '#fdba74' : '#fda4af') : v.effect === 'missense' ? '#c4b5fd' : '#cbd5e1' }}>{v.conseq || CLIN_EFFECT_WORD[v.effect]}</span>
            <span className="font-mono text-[10px] text-slate-300"> {v.protein ?? v.hgvs ?? ''}</span>
            {n && <span className="text-slate-400"> · {n.call === 'escape' ? 'escape' : 'NMD'}{n.assumed ? '?' : ''}: {n.why.replace(/: (escapes NMD|NMD assumed|NMD)$/, '')}</span>}
            {pinned && <> · <a href={clinvarUrl(v)} target="_blank" rel="noreferrer" className="text-accent-300 hover:underline">ClinVar ↗</a></>}
          </div>
        );
      })}</div>}
      {!genes && shown.map(v => {
        const big = v.end - v.start >= 50;
        const hits = big ? cnv.filter(c => (c.ev.kind === 'loss' ? v.kind === 'loss' : v.kind === 'gain') && clinOverlap(c.ev, v) > 0) : [];
        return (
          <div key={`${v.id}:${v.start}`} className="mt-1.5 pt-1 border-t border-white/10 first:border-0 first:pt-0">
            <div><span className="px-1 rounded text-[10px] font-bold" style={{ background: v.code === 'PG' ? '#b91c1c' : '#f97316' }}>{v.code === 'PG' ? 'P' : 'LP'}</span> <b>{v.sign}</b>
              {!big && <span className="ml-1 px-1 rounded text-[10px] font-semibold" style={{ background: CLIN_EFFECT_COLOR[v.effect], color: '#fff' }} title={v.conseq || 'no consequence given'}>{CLIN_EFFECT_WORD[v.effect]}{v.conseq && CLIN_EFFECT_WORD[v.effect] !== v.conseq ? ` · ${v.conseq}` : ''}</span>}
              <span className="text-slate-300"> · {v.type || CLIN_KIND_WORD[v.kind]}{big ? ` · ${kb(v.end - v.start)}` : ''}</span>
              <span className="text-amber-300" title={v.review}> {'★'.repeat(Math.min(4, v.stars))}{'☆'.repeat(Math.max(0, 4 - Math.min(4, v.stars)))}</span></div>
            {v.genes.length > 0 && <div className="text-slate-300">gene{v.genes.length > 1 ? 's' : ''}: <i>{v.genes.slice(0, 8).join(', ')}{v.genes.length > 8 ? ` +${v.genes.length - 8}` : ''}</i></div>}
            {v.conditions && <div className="text-slate-100">{v.conditions.length > 220 ? `${v.conditions.slice(0, 220)}…` : v.conditions}</div>}
            {v.effect === 'lof' && nmdOf?.(v) && (() => { const n = nmdOf(v)!; return <div style={{ color: n.call === 'escape' ? '#fdba74' : '#fda4af' }}>{n.call === 'escape' ? 'NMD escape' : 'NMD'}{n.assumed ? ' (assumed)' : ''}: {n.why}</div>; })()}
            {(v.hgvs || v.protein) && <div className="text-slate-400 font-mono text-[10px] break-all">{v.hgvs}{v.protein ? ` · ${v.protein}` : ''}</div>}
            <div className="text-slate-400">{chrom}:{fmt(v.start + 1)}{v.end - v.start > 1 ? `-${fmt(v.end)}` : ''} · {v.id}{v.origin ? ` · ${v.origin}` : ''}{v.review ? ` · ${v.review}` : ''}</div>
            {hits.map(h => <div key={h.ev.start} className="text-emerald-300">The event's {h.ev.kind} {h.ev.letters} overlaps it ({Math.round(clinOverlap(h.ev, v) * 100)} % reciprocal).</div>)}
            {pinned && <a href={clinvarUrl(v)} target="_blank" rel="noreferrer" className="text-accent-300 hover:underline">open in ClinVar ↗</a>}
          </div>
        );
      })}
      {vs.length > shown.length && <div className="mt-1 text-slate-400">+{vs.length - shown.length} more {genes ? 'in this gene: zoom on it and show each variant' : 'here'}{genes ? '' : at != null ? ' (larger ones): “▸ all SVs” lists them one by one' : ': zoom in to separate them'}</div>}
      {pinned
        ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={onZoom} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on it</button></div>
        : <div className="mt-0.5 text-slate-400">click to keep this card open (ClinVar links)</div>}
      <div className="mt-0.5 text-slate-500">ClinVar through the UCSC Genome Browser (tracks clinvarCnv, clinvarMain); conflicting interpretations left out.</div>
    </div>
  );
}
/**
 * How a loss (a deletion) or a gain (a duplication) of the event matches a region, in words, or null: by the rule the
 * sheet gives that type (touching its critical gene, spanning its smallest region of overlap), else a reciprocal overlap.
 */
const achroMatch = (ev: CnvEvent, p: AchroPlaced): string | null => {
  const t = ev.kind === 'loss' ? 'del' : 'dup', sp = p.types[t], rule = p.r.typed?.[t];
  if (!sp) return null;
  const ov = Math.max(0, Math.min(ev.end, sp.end) - Math.max(ev.start, sp.start)), len = Math.max(sp.end - sp.start, 1);
  if (rule?.match === 'overlaps') return ov > 0 ? `takes ${ov >= len ? 'all' : `${Math.round((ov / len) * 100)} %`} of ${rule.genes?.join(', ') ?? 'the critical region'}` : null;
  if (rule?.match === 'contains') return ov >= len ? 'spans the smallest region of overlap' : null;
  const ro = ov / Math.max(ev.end - ev.start, len);
  return ro >= RO_MATCH ? `overlaps it (${Math.round(ro * 100)} % reciprocal)` : null;
};

/** The card of an AChro-Puce region: its classes, a summary, how it was placed, the event's matching loss or gain, the sheet. */
function AchroCard({ p, type, chrom, cnv, pinned, place, onClose, onZoom, cardRef }: { p: AchroPlaced; type?: 'del' | 'dup'; chrom: string; cnv: { ev: CnvEvent }[]; pinned: boolean; place: React.CSSProperties; onClose: () => void; onZoom: () => void; cardRef?: (el: HTMLDivElement | null) => void }) {
  const hits = cnv.flatMap(c => { const m = achroMatch(c.ev, p); return m ? [{ ev: c.ev, m }] : []; });
  const spans = (['del', 'dup'] as const).filter(t => p.types[t]);
  const differ = spans.length === 2 && (p.types.del!.start !== p.types.dup!.start || p.types.del!.end !== p.types.dup!.end);
  const where = (sp: AchroSpan) => `${chrom}:${(sp.start + 1).toLocaleString('en-US')}-${sp.end.toLocaleString('en-US')} ${sp.approx ? '(approximate) ' : ''}· placed from the ${sp.how}.`;
  return (
    <div ref={cardRef} className={`absolute z-20 rounded-lg bg-slate-900 text-white text-[11px] px-3 py-2 shadow-xl w-max max-w-[440px] ${pinned ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={place}
      onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
      <div className="font-semibold flex items-start gap-2"><span className="flex-1">{p.r.name} <span className="font-normal text-slate-400">· AChro-Puce</span></span>
        {pinned && <button onClick={onClose} className="-mr-1 -mt-0.5 px-1 text-slate-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
      <table className="mt-1 w-full"><tbody>
        {p.r.rows.map((r, i) => (
          <tr key={i} className={`align-top ${type && r.type !== type ? 'opacity-45' : ''}`}><td className="pr-2 py-0.5 whitespace-nowrap font-bold" style={{ color: r.type === 'del' ? '#fca5a5' : '#93c5fd' }}>{typeWord(r)}</td>
            <td className="py-0.5"><span className="font-bold" style={{ color: ACHRO_COLOR[r.cls] === '#64748b' ? '#cbd5e1' : ACHRO_COLOR[r.cls] }}>{r.cls === 'PIEV' ? 'PIEV' : `Class ${r.cls}`}</span> <span className="text-slate-300">{r.label}</span>
              {r.also && <div><span className="font-bold" style={{ color: ACHRO_COLOR[r.also.cls] }}>{r.also.cls === 'PIEV' ? 'PIEV' : `Class ${r.also.cls}`}</span> <span className="text-slate-300">{r.also.label}</span></div>}</td></tr>
        ))}
      </tbody></table>
      <div className="mt-1 text-slate-200 leading-snug">{p.r.summary}</div>
      {hits.map(h => <div key={h.ev.start} className="mt-1 text-emerald-300">The event's {h.ev.kind} {h.ev.letters} {h.m}: {p.r.rows.filter(r => r.type === (h.ev.kind === 'loss' ? 'del' : 'dup')).map(r => `${r.cls === 'PIEV' ? 'PIEV' : `class ${r.cls}`}${r.sex ? ` (${r.sex})` : ''}`).join(', ') || 'no class given for this type'}.</div>)}
      {differ
        ? spans.map(t => <div key={t} className={`mt-1 text-slate-400 ${type && t !== type ? 'opacity-60' : ''}`}><b style={{ color: t === 'del' ? '#fca5a5' : '#93c5fd' }}>{t === 'del' ? 'Deletion' : 'Duplication'}</b> {where(p.types[t]!)}</div>)
        : <div className="mt-1 text-slate-400">{where(p.types[spans[0]] ?? p)}</div>}
      <div className="mt-0.5 text-slate-400">Sheet “{p.r.sheet}”, {p.r.authors} — {ACHRO_SOURCE.title}, {ACHRO_SOURCE.version}. The summary is written for this viewer: read the sheet.</div>
      {pinned
        ? <div className="mt-1.5 pt-1.5 border-t border-white/15 flex gap-2"><button onClick={onZoom} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on it</button>
            <a href={ACHRO_SOURCE.url} target="_blank" rel="noreferrer" className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">📄 the catalogue (PDF)</a>
            {p.r.coords?.some(c => c.from === 'table2') && <a href={ACHRO_GUIDELINES.url} target="_blank" rel="noreferrer" title={ACHRO_GUIDELINES.cite} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">📄 guidelines, Table 2</a>}</div>
        : <div className="mt-0.5 text-slate-400">click to keep this card open</div>}
    </div>
  );
}

/** One chip per loss or gain of the event that matches an AChro-Puce region (achroMatch), with the region's class for that type. */
function AchroSummary({ cnv, placed }: { cnv: { ev: CnvEvent }[]; placed: AchroPlaced[] }) {
  const hits = cnv.flatMap(c => placed.filter(p => achroMatch(c.ev, p)).map(p => ({ ev: c.ev, p, rows: p.r.rows.filter(r => r.type === (c.ev.kind === 'loss' ? 'del' : 'dup')) })));
  if (!hits.length) return null;
  return (
    <div className="px-5 mt-1.5 flex flex-wrap gap-1.5 text-[11.5px]" title={`${ACHRO_SOURCE.title}, ${ACHRO_SOURCE.version}: a loss or gain of the event overlapping a region reciprocally by ${Math.round(RO_MATCH * 100)} % or more; for a non-recurrent CNV, a loss touching the critical gene or a gain spanning the smallest region of overlap`}>
      <span className="text-slate-500 self-center">AChro-Puce:</span>
      {hits.map(({ ev, p, rows }) => (
        <span key={`${ev.start}${p.r.id}`} className="inline-flex items-center gap-1 rounded-full px-2.5 py-0.5 border bg-white" style={{ borderColor: `${ev.kind === 'loss' ? ACHRO_TYPE_COLOR.del : ACHRO_TYPE_COLOR.dup}66` }}>
          <b style={{ color: ev.kind === 'loss' ? ACHRO_TYPE_COLOR.del : ACHRO_TYPE_COLOR.dup }}>{ev.kind} {ev.letters}</b> ≈ {p.r.name} {ev.kind === 'loss' ? 'deletion' : 'duplication'}:
          {rows.length ? rows.map((r, i) => <b key={i} style={{ color: ACHRO_COLOR[r.cls] }}>{r.cls === 'PIEV' ? 'PIEV' : `class ${r.cls}`}{r.sex ? ` (${r.sex})` : ''}{r.also ? ` + ${r.also.cls === 'PIEV' ? 'PIEV' : `class ${r.also.cls}`} (${r.also.label})` : ''}</b>) : <span className="text-slate-500">not classified for this type</span>}
        </span>
      ))}
    </div>
  );
}

// ======================== PanelApp picker ========================

const TYPE_LABEL = { loss: 'loss', gain: 'gain', both: 'gain and loss', unknown: 'type not given' } as const;
const pct = (v: number) => `${Math.round(v * 100)} %`;
const evName = (ev: CnvEvent) => `${ev.kind} ${ev.letters} (${kb(ev.end - ev.start)})`;

/** The card of a population CNV: its extents, its frequency, its overlap with the event's losses and gains, and every field the source gives. */
function DgvCard({ v, chrom, cnv, pinned, place, onClose, onZoom, cardRef }: { v: PopulationSv; chrom: string; cnv: { ev: CnvEvent; matches: PopMatch[] }[]; pinned: boolean; place: React.CSSProperties; onClose: () => void; onZoom: () => void; cardRef?: (el: HTMLDivElement | null) => void }) {
  const span = (r: [number, number]) => `${chrom}:${fmt(r[0] + 1)}-${fmt(r[1])} (${kb(r[1] - r[0])})`;
  const sameInner = v.inner[0] === v.outer[0] && v.inner[1] === v.outer[1];
  const rows = cnv.map(c => ({ ev: c.ev, m: c.matches.find(m => m.sv.id === v.id && m.sv.outer[0] === v.outer[0]) })).filter(x => x.m);
  return (
    <div ref={cardRef} className={`absolute z-20 rounded-md bg-slate-900 text-white text-[11px] px-2.5 py-1.5 shadow-lg w-max max-w-[400px] ${pinned ? 'ring-1 ring-white/30' : 'pointer-events-none opacity-95'}`} style={place}
      onMouseDown={e => e.stopPropagation()} onClick={e => e.stopPropagation()}>
      <div className="font-semibold flex items-start gap-2"><span className="flex-1"><span style={{ color: { loss: '#fca5a5', gain: '#93c5fd', both: '#c4b5fd', unknown: '#cbd5e1' }[v.type] }}>●</span> {v.id} · population {TYPE_LABEL[v.type]} (DGV Gold Standard)</span>
        {pinned && <button onClick={onClose} className="-mr-1 -mt-0.5 px-1 text-slate-400 hover:text-white text-sm leading-none" title="Close">×</button>}</div>
      <div className="text-slate-300">{sameInner ? span(v.outer) : <>outer {span(v.outer)}<br />inner {span(v.inner)}</>}</div>
      <div className="mt-0.5">{[v.freqText != null ? `frequency ${v.freqText}` : 'frequency not given', v.samples != null ? `${fmt(v.samples)} samples` : '', v.studies != null ? `${fmt(v.studies)} studies` : ''].filter(Boolean).join(' · ')}</div>
      {rows.map(({ ev, m }) => (
        <div key={ev.start} className={`mt-0.5 ${m!.sameType && m!.ro >= RO_MATCH ? 'text-emerald-300' : 'text-slate-300'}`}>
          vs the {evName(ev)}: reciprocal overlap {pct(m!.ro)}{sameInner ? '' : ` (inner ${pct(m!.roInner)}, outer ${pct(m!.roOuter)})`}{m!.covers ? ', its inner extent holds the whole event' : ''}{m!.sameType ? '' : ', other type'}
        </div>
      ))}
      {pinned && v.fields.length > 0 && (
        <details className="mt-1 text-slate-400"><summary className="cursor-pointer">all fields ({v.fields.length})</summary>
          <div className="mt-0.5 max-h-40 overflow-auto font-mono text-[10px] leading-snug">{v.fields.map(([k2, val]) => <div key={k2}><span className="text-slate-500">{k2}</span> {val.length > 160 ? `${val.slice(0, 160)}…` : val}</div>)}</div>
        </details>
      )}
      {pinned
        ? <div className="mt-1.5 pt-1.5 border-t border-white/15"><button onClick={onZoom} className="px-2 py-0.5 rounded bg-white/10 hover:bg-white/20 text-white">🔍 zoom on it</button></div>
        : <div className="mt-0.5 text-slate-400">click to keep this card open</div>}
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
      <span className="text-slate-500 self-center">DGV Gold:</span>
      {loading ? <span className="text-slate-500">loading…</span> : error ? <span className="text-slate-500">unavailable ({error})</span> : cnv.map(({ ev, matches }) => {
        const same = matches.filter(m => m.sameType && m.ro >= RO_MATCH);
        const best = matches.find(m => m.sameType);
        const seen = same.length > 0, top = same.sort((a, b) => (b.sv.freq ?? 0) - (a.sv.freq ?? 0))[0];
        return (
          <span key={ev.start} className={`inline-flex items-center gap-1 rounded-full px-2.5 py-0.5 border ${seen ? 'bg-emerald-50 border-emerald-300 text-emerald-800' : 'bg-slate-50 border-slate-300 text-slate-700'}`}>
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
      <button onClick={() => setOpen(o => !o)} className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium whitespace-nowrap ${panel ? 'bg-accent-50 border-accent-200 text-accent-600' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'}`}
        title="Rate the genes against a PanelApp gene panel (green / amber / red)">
        <Icon name="layers" size={13} />{panel ? `${panel.name}${panel.version ? ` v${panel.version}` : ''} · ${panel.genes.length} genes` : 'Gene panel (PanelApp)…'}
      </button>
      {panel && <button onClick={() => onPanel(null)} className="ml-1 text-slate-400 hover:text-red-600" title="Remove the panel">×</button>}
      {open && (
        <div className="absolute z-40 mt-1 w-[460px] rounded-lg border border-slate-300 bg-white shadow-2xl p-2 text-xs" onMouseDown={e => e.stopPropagation()}>
          <div className="flex gap-1.5">
            <select value={src} onChange={e => setSrc(e.target.value)} className="border rounded px-1 py-1">{PANEL_SOURCES.map(s => <option key={s.id} value={s.id}>{s.label}</option>)}</select>
            <input autoFocus value={q} onChange={e => setQ(e.target.value)} placeholder="search a panel or a disorder" className="flex-1 border rounded px-2 py-1 outline-none focus:border-accent-300" />
            <button onClick={() => setOpen(false)} className="px-1.5 text-slate-400 hover:text-slate-700">×</button>
          </div>
          <div className="mt-1.5 max-h-[260px] overflow-y-auto">
            {list.loading && <div className="text-slate-500 p-2">loading the panels of {source.label}…</div>}
            {list.error && <div className="text-amber-700 p-2">{list.error}</div>}
            {hits.map(p => (
              <button key={p.id} disabled={!!busy} onClick={async () => { setBusy(p.id); try { const lp = await loadPanel(source, p.id); onPanel(lp); setOpen(false); } catch (e: any) { setList(l => ({ ...l, error: `${p.name}: ${e?.message || e}` })); } finally { setBusy(null); } }}
                className="w-full text-left px-2 py-1 rounded hover:bg-accent-50 flex justify-between gap-2">
                <span className="truncate">{p.name}</span><span className="text-slate-400 whitespace-nowrap">{busy === p.id ? 'loading…' : `v${p.version} · ${p.genes} genes`}</span>
              </button>
            ))}
            {!list.loading && !list.error && !hits.length && <div className="text-slate-400 p-2">no panel matches</div>}
          </div>
          <label className="mt-1.5 block border-t border-slate-100 pt-1.5 text-slate-600 cursor-pointer hover:text-accent-700">
            or load a file: a panel downloaded from PanelApp (TSV) or a list of gene symbols
            <input type="file" accept=".tsv,.txt,.csv" className="hidden" onChange={async e => { const f = e.target.files?.[0]; if (f) { onPanel(parsePanelFile(f.name, await f.text())); setOpen(false); } e.target.value = ''; }} />
          </label>
        </div>
      )}
    </span>
  );
}
