/**
 * The rearranged genome browser of the SV cartoon: the reference chromosome with its genes on top, the rebuilt
 * chromosome with the same genes in their new coordinates underneath, and a ribbon from every segment to wherever it
 * went (twisted when inverted, doubled when duplicated, none when lost). Genes are coloured by what happened to them
 * (svGenes.ts); a cut gene shows torn ends, a fusion its two partners joined under one label. Chromatin domains from a
 * BED the user loads are drawn as triangles over the reference and under the rebuilt chromosome, new domains outlined.
 * Underneath, every touched gene with its consequence, ClinGen dosage scores, gnomAD constraint and its rating in the
 * chosen PanelApp panel.
 */
import React, { useEffect, useMemo, useRef, useState } from 'react';
import { reconstruct, type SvReconstruction } from './svModel';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';
import { geneStory, dosageReading, dosageScoreText, ranksText, toDer, LOEUF_INTOLERANT, PLI_INTOLERANT, NEAR_BREAKPOINT_DEFAULT,
  type SvGene, type SvTranscript, type GeneAnnotations, type GenePanel, type Domain, type GeneImpact, type Chain, type Fragment, type DerLayout, type DosageScore, type Constraint, type ImpactKind } from './svGenes';
import { IMPACT, KIND_COLOR } from './svColors';
import { PANEL_SOURCES, listPanels, loadPanel, parsePanelFile, type PanelSummary, type LoadedPanel } from '../../standalone/panelapp';

const FONT = 'Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const MONO = 'ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
const VW = 1200, X0 = 112, X1 = 1186;
const TX_H = 13, LABEL_H = 13, LANE_GAP = 5, MAX_LANES = 8, TAD_H = 58, BAND_H = 20, RULER_H = 18, RIB_H = 118;
const C = { text: '#1f2937', muted: '#5b6b7e', faint: '#8796a8', grid: '#e3e8ee', line: '#d3dbe5', lost: '#dc2626', accent: '#0f766e', fusion5: '#7c3aed', fusion3: '#a855f7' };
const CONF = { 3: { label: 'Green', color: '#16a34a', bg: '#dcfce7' }, 2: { label: 'Amber', color: '#d97706', bg: '#fef3c7' }, 1: { label: 'Red', color: '#dc2626', bg: '#fee2e2' } } as const;
const fmt = (n: number) => n.toLocaleString('en-US');
const kb = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} Mb` : n >= 1000 ? `${(n / 1000).toFixed(n >= 1e4 ? 0 : 1)} kb` : `${n} bp`);

export interface SvGeneLoaders {
  getTranscripts?: (chrom: string, start: number, end: number) => Promise<SvTranscript[]>;
  getAnnotations?: (chrom: string, start: number, end: number, symbols: string[]) => Promise<{ dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[] }>;
  /** the models the viewer already has (the gene in view): used when the source cannot list a window's transcripts */
  fallback: SvTranscript[];
}

interface Props {
  rec: SvReconstruction; chrom: string; window: { start: number; end: number };
  segColors: Map<number, string>; flankCopies: number; zygosity: 'het' | 'hom'; build: string;
  loaders: SvGeneLoaders;
}

// ======================== Genes of the window ========================

/** One gene per symbol, its principal model the MANE Select (else RefSeq Select, else the longest coding sequence). */
function groupGenes(txs: SvTranscript[]): SvGene[] {
  const by = new Map<string, SvTranscript[]>();
  for (const t of txs) by.set(t.gene, [...(by.get(t.gene) ?? []), t]);
  const cds = (t: SvTranscript) => (t.cdsStart != null && t.cdsEnd != null ? t.exons.reduce((a, e) => a + Math.max(0, Math.min(e.end, t.cdsEnd!) - Math.max(e.start, t.cdsStart!)), 0) : 0);
  const rank = (t: SvTranscript) => (t.status === 'mane' ? 0 : t.status === 'select' ? 1 : t.status === 'mane-plus' ? 2 : 3);
  return [...by.entries()].map(([symbol, list]) => {
    const sorted = [...list].sort((a, b) => rank(a) - rank(b) || cds(b) - cds(a));
    return { symbol, strand: sorted[0].strand, start: Math.min(...list.map(t => t.start)), end: Math.max(...list.map(t => t.end)), transcripts: sorted, principal: sorted[0] };
  }).sort((a, b) => a.start - b.start);
}

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

// ======================== Component ========================

export default function SvGenomeBrowser({ rec: eventRec, chrom, window: eventWin, segColors, flankCopies, zygosity, build, loaders }: Props) {
  const [mode, setMode] = useState<'mane' | 'all'>('mane');
  const [txs, setTxs] = useState<{ list: SvTranscript[]; loading: boolean; error?: string; source: string }>({ list: [], loading: true, source: '' });
  const [ann, setAnn] = useState<{ dosage: Map<string, DosageScore>; constraint: Map<string, Constraint>; notes: string[]; loading: boolean }>({ dosage: new Map(), constraint: new Map(), notes: [], loading: false });
  const [panelData, setPanelData] = useState<LoadedPanel | null>(() => store.get<LoadedPanel>('sashimi.svPanel'));
  const [domainsFile, setDomainsFile] = useState<{ name: string; text: string } | null>(() => store.get('sashimi.svDomains'));
  const [hover, setHover] = useState<string | null>(null);
  const [tip, setTip] = useState<{ x: number; y: number; gene: string } | null>(null);
  const [showUnchanged, setShowUnchanged] = useState(false);
  const [showNear, setShowNear] = useState(false);
  const [view, setView] = useState({ z: 1, pan: 0 });
  const svgRef = useRef<SVGSVGElement>(null);
  const boxRef = useRef<HTMLDivElement>(null);
  const drag = useRef<{ x: number; pan: number } | null>(null);
  const [near, setNear] = useState<number>(() => store.get<number>('sashimi.svNear') ?? NEAR_BREAKPOINT_DEFAULT);

  // the region shown: the event's window, widened to `near` on each side of the outermost breakpoints (clamped to the
  // chromosome), so the genes a position effect could reach are there; the same junctions cut it (same segments, longer flanks)
  const bpKey = eventRec.junctions.flatMap(j => [j.a, j.b]).filter(e => e.chrom === chrom).map(e => e.pos).sort((a, b) => a - b).join(',');
  const win = useMemo(() => {
    const bps = bpKey ? bpKey.split(',').map(Number) : [];
    if (!bps.length || near <= 0) return eventWin;
    const len = chromLength(build as GenomeBuild, chrom) ?? Infinity;
    return { start: Math.max(0, Math.min(eventWin.start, bps[0] - near)), end: Math.min(Math.max(len, eventWin.end), Math.max(eventWin.end, bps[bps.length - 1] + near)) };
  }, [eventWin.start, eventWin.end, bpKey, near, build, chrom]); // eslint-disable-line react-hooks/exhaustive-deps
  const rec = useMemo(() => (win === eventWin ? eventRec : reconstruct([{ chrom, start: win.start, end: win.end }], eventRec.junctions, eventRec.choice)), [eventRec, win]); // eslint-disable-line react-hooks/exhaustive-deps

  // the NM_ models of the window
  useEffect(() => {
    let alive = true;
    if (!loaders.getTranscripts) { setTxs({ list: loaders.fallback, loading: false, source: 'the gene in view only (this source cannot list the window\'s genes)' }); return; }
    setTxs(t => ({ ...t, loading: true, error: undefined }));
    loaders.getTranscripts(chrom, win.start, win.end)
      .then(list => { if (alive) setTxs({ list, loading: false, source: 'RefSeq NM_ models (UCSC ncbiRefSeqCurated), MANE status from the UCSC mane track' }); })
      .catch(e => { if (alive) setTxs({ list: loaders.fallback, loading: false, error: `gene models of the window unavailable (${e?.message || e}): the gene in view only`, source: '' }); });
    return () => { alive = false; };
  }, [chrom, win.start, win.end]); // eslint-disable-line react-hooks/exhaustive-deps

  const genes = useMemo(() => groupGenes(txs.list), [txs.list]);
  const symbols = useMemo(() => genes.map(g => g.symbol), [genes]);
  useEffect(() => {
    if (!loaders.getAnnotations || !symbols.length) return;
    let alive = true;
    setAnn(a => ({ ...a, loading: true }));
    loaders.getAnnotations(chrom, win.start, win.end, symbols).then(r => { if (alive) setAnn({ ...r, loading: false }); }).catch(e => { if (alive) setAnn({ dosage: new Map(), constraint: new Map(), notes: [`annotations unavailable (${e?.message || e})`], loading: false }); });
    return () => { alive = false; };
  }, [symbols.join(','), chrom]); // eslint-disable-line react-hooks/exhaustive-deps

  const domains = useMemo(() => (domainsFile ? parseDomainsBed(domainsFile.text).get(chrom.startsWith('chr') ? chrom : `chr${chrom}`) ?? [] : []), [domainsFile, chrom]);
  const story = useMemo(() => geneStory(rec, genes, win, { flankCopies, zygosity, domains, near }), [rec, genes, win, flankCopies, zygosity, domains, near]);
  const panel = useMemo(() => (panelData ? toPanel(panelData) : null), [panelData]);
  const annotations: GeneAnnotations = useMemo(() => ({ dosage: ann.dosage, constraint: ann.constraint, panel }), [ann, panel]);
  const impactOf = useMemo(() => new Map(story.genes.map(g => [g.gene.symbol, g])), [story]);
  /** in a crowded view, the genes that keep their label: every changed gene, and the nearby ones worth a look */
  const worthLabel = (gi: GeneImpact | undefined) => !!gi && gi.kind !== 'unaffected' && (gi.kind !== 'nearby' || notable(gi, annotations));

  // ---- geometry ----
  const ders = story.ders.filter(d => d.rearranged);
  const winLen = win.end - win.start;
  const k = (X1 - X0) / Math.max(winLen, ...ders.map(d => d.length), 1);
  const z = view.z, pan = view.pan;
  const zx = (x: number) => (x - X0) * z + X0 + pan;
  const refX = (pos: number) => zx(X0 + ((X1 - X0) - winLen * k) / 2 + (pos - win.start) * k);
  const derX = (d: DerLayout, dpos: number) => zx(X0 + ((X1 - X0) - d.length * k) / 2 + dpos * k);

  // wheel zoom at the pointer, drag to pan, double-click to fit
  useEffect(() => {
    const el = svgRef.current; if (!el) return;
    const onWheel = (e: WheelEvent) => {
      if (!e.ctrlKey && !e.metaKey && Math.abs(e.deltaY) < Math.abs(e.deltaX)) return;
      e.preventDefault();
      const r = el.getBoundingClientRect(), xc = ((e.clientX - r.left) / r.width) * VW;
      setView(v => { const z2 = Math.max(1, Math.min(200, v.z * (e.deltaY < 0 ? 1.2 : 1 / 1.2))); return { z: z2, pan: z2 === 1 ? 0 : xc - X0 - (xc - X0 - v.pan) * (z2 / v.z) }; });
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, [txs.loading]);
  /** The view showing [lo, hi) of the reference over `frac` of the width. */
  const fitView = (lo: number, hi: number, frac = 1) => {
    const base = (pos: number) => X0 + ((X1 - X0) - winLen * k) / 2 + (pos - win.start) * k;
    const a0 = base(lo), b0 = base(hi), z2 = Math.max(1, Math.min(200, (X1 - X0) * frac / Math.max(1, b0 - a0)));
    return z2 <= 1 ? { z: 1, pan: 0 } : { z: z2, pan: (X0 + X1) / 2 - X0 - ((a0 + b0) / 2 - X0) * z2 };
  };
  const eventView = () => fitView(eventWin.start, eventWin.end);
  // open on the event; zooming out shows the genes around it
  useEffect(() => { setView(eventView()); }, [win.start, win.end, ders.map(d => d.length).join()]); // eslint-disable-line react-hooks/exhaustive-deps
  const zoomTo = (lo: number, hi: number) => {
    const a = refX(lo), b = refX(hi); const base = (x: number) => (x - X0 - pan) / z + X0;
    const a0 = base(a), b0 = base(b), z2 = Math.max(1, Math.min(200, (X1 - X0) * 0.7 / Math.max(1, b0 - a0)));
    setView({ z: z2, pan: (X0 + X1) / 2 - X0 - ((a0 + b0) / 2 - X0) * z2 });
  };

  // ---- layout of the reference genes ----
  const shown = (g: SvGene) => (mode === 'mane' ? [g.principal] : g.transcripts);
  const labelW = (s: string) => s.length * 6.4 + 16;
  /**
   * Genes in lanes, each with its label; when that would take more than MAX_LANES lanes (a gene-dense region seen whole),
   * glyphs alone are packed (sharing a lane past the cap) and only the genes worth a label keep it (`compact`).
   */
  function pack<T extends { x0: number; x1: number; xg: number; h: number }>(items: T[]): { items: (T & { lane: number; compact: boolean })[]; laneY: number[]; height: number } {
    const run = (ext: (it: T) => number, cap: number) => {
      const lanes: number[] = [], laneH: number[] = [];
      const placed = [...items].sort((a, b) => a.x0 - b.x0).map(it => {
        let lane = lanes.findIndex(end => end + 8 < it.x0);
        if (lane < 0 && lanes.length < cap) { lane = lanes.length; lanes.push(-Infinity); laneH.push(0); }
        if (lane < 0) lane = lanes.indexOf(Math.min(...lanes));
        lanes[lane] = Math.max(lanes[lane], ext(it)); laneH[lane] = Math.max(laneH[lane], it.h);
        return { ...it, lane };
      });
      return { placed, laneH };
    };
    let r = run(it => it.x1, Infinity), compact = false;
    if (r.laneH.length > MAX_LANES) { r = run(it => it.xg, MAX_LANES); compact = true; }
    const laneY: number[] = []; let y = 0;
    for (const h of r.laneH) { laneY.push(y); y += h + LANE_GAP; }
    return { items: r.placed.map(it => ({ ...it, compact })), laneY, height: Math.max(26, y) };
  }
  const inView = (a: number, b: number) => Math.max(a, b) >= X0 - 2 && Math.min(a, b) <= X1 + 2;
  const refItems = pack(genes.filter(g => inView(refX(g.start), refX(g.end))).map(g => { const a = Math.max(X0 + 2, refX(Math.max(g.start, win.start))), b = refX(Math.min(g.end, win.end)); return { g, x0: a, x1: Math.max(b, a + labelW(g.symbol)), xg: Math.max(b, a + 3), h: LABEL_H + shown(g).length * TX_H }; }));

  // ---- layout of the rebuilt chromosome's genes: one block per chain (a fusion's two partners in one block) ----
  interface DerBlock { der: DerLayout; chains: Chain[]; genes: string[]; label: string; x0: number; x1: number; xg: number; h: number; fusion: GeneImpact['principal']['fusions'][number] | null }
  const derBlocks = useMemo(() => {
    const keep = (c: Chain) => c.parts.every(p => (mode === 'mane' ? p.tx === impactOf.get(p.tx.gene)?.gene.principal : true));
    const out = new Map<DerLayout, DerBlock[]>();
    for (const d of ders) out.set(d, []);
    const copyNo = new Map<string, number>();
    for (const c of story.chains) {
      if (!keep(c)) continue;
      const d = story.ders[c.der]; if (!d.rearranged) continue;
      const geneList = [...new Set(c.parts.map(p => p.tx.gene))];
      const xs = c.parts.flatMap(p => [derX(d, toDer(p.dp, p.lo)), derX(d, toDer(p.dp, p.hi))]);
      if (!inView(Math.min(...xs), Math.max(...xs))) continue;
      const x0 = Math.max(X0 + 2, Math.min(...xs)), x1 = Math.max(...xs);
      const first = c.parts[0], last = c.parts[c.parts.length - 1];
      const gi = impactOf.get(geneList[0]);
      const fusion = geneList.length > 1 ? gi?.principal.fusions.find(f => f.five.gene === geneList[0] && f.three.gene === geneList[1]) ?? impactOf.get(geneList[1])?.principal.fusions.find(f => f.five.gene === geneList[0]) ?? null : null;
      const five = first.tx.strand > 0 ? !first.cutLo || first.edgeLo : !first.cutHi || first.edgeHi;
      const three = last.tx.strand > 0 ? !last.cutHi || last.edgeHi : !last.cutLo || last.edgeLo;
      let label = geneList.join('::');
      if (geneList.length === 1) {
        if (!five && !three) label += ' · inner part'; else if (!five) label += ' · 3′ part'; else if (!three) label += ' · 5′ part';
        else { const key = `${c.der}|${geneList[0]}|${first.tx.id}`; const n = (copyNo.get(key) ?? 0) + 1; copyNo.set(key, n); if ((gi?.principal.intact ?? 0) > 1 && n > 1) label += ` · copy ${n}`; }
      }
      const list = out.get(d)!;
      // in "all NM" mode the transcripts of one gene over the same span share a block
      const same = list.find(b => b.genes.join() === geneList.join() && Math.abs(b.x0 - x0) < 2 && Math.abs(b.x1 - x1) < 2);
      if (same) { same.chains.push(c); same.h += TX_H; continue; }
      list.push({ der: d, chains: [c], genes: geneList, label, x0, x1: Math.max(x1, x0 + labelW(label) + (fusion ? 70 : 0)), xg: Math.max(x1, x0 + 3), h: LABEL_H + TX_H, fusion });
    }
    return out;
  }, [story, mode, view, k]); // eslint-disable-line react-hooks/exhaustive-deps
  const derPacked = ders.map(d => pack(derBlocks.get(d) ?? []));

  // ---- vertical layout ----
  let y = 10;
  const yRefTad = domains.length ? y : -1; if (domains.length) y += TAD_H + 6;
  const yRefGenes = y; y += refItems.height + 4;
  const yRefRuler = y; y += RULER_H;
  const yRefBand = y; y += BAND_H;
  const yRib = y; y += RIB_H;
  const derY = ders.map((d, i) => {
    const band = y; y += BAND_H;
    const ruler = y; y += RULER_H + 6;
    const genesY = y; y += derPacked[i].height + 4;
    const tad = domains.length ? y : -1; if (domains.length) y += TAD_H;
    y += 18;
    return { band, ruler, genes: genesY, tad };
  });
  const H = y + 4;

  // ---- drawing helpers ----
  const dim = (sym: string | string[]) => (hover && !(Array.isArray(sym) ? sym.includes(hover) : sym === hover) ? 0.28 : 1);
  const enter = (sym: string) => (e: React.MouseEvent) => { setHover(sym); const r = boxRef.current!.getBoundingClientRect(); setTip({ x: e.clientX - r.left, y: e.clientY - r.top, gene: sym }); };
  const leave = () => { setHover(null); setTip(null); };
  const chevrons = (x0: number, x1: number, yy: number, dir: 1 | -1, color: string, key: string) => {
    const out: JSX.Element[] = [];
    for (let x = Math.min(x0, x1) + 10; x < Math.max(x0, x1) - 6; x += 16) out.push(<path key={`${key}c${x}`} d={`M${x - 2.5 * dir},${yy - 3} L${x + 1.5 * dir},${yy} L${x - 2.5 * dir},${yy + 3}`} fill="none" stroke={color} strokeWidth={1} opacity={0.7} />);
    return out;
  };
  /** A transcript (or its part [lo, hi]) drawn with `X` mapping reference positions; `dir` its direction on the drawing. */
  const txGlyph = (t: SvTranscript, X: (p: number) => number, yy: number, color: string, key: string, lo = t.start, hi = t.end, dir: 1 | -1 = t.strand) => {
    const a = X(lo), b = X(hi), els: JSX.Element[] = [];
    els.push(<line key={`${key}i`} x1={a} x2={b} y1={yy} y2={yy} stroke={color} strokeWidth={1.1} />);
    els.push(...chevrons(a, b, yy, dir, color, key));
    for (const e of t.exons) {
      const s = Math.max(e.start, lo), en = Math.min(e.end, hi);
      if (en <= s) continue;
      const xa = Math.min(X(s), X(en)), w = Math.max(1.2, Math.abs(X(en) - X(s)));
      els.push(<rect key={`${key}u${e.rank}`} x={xa} y={yy - 3.5} width={w} height={7} fill={color} opacity={0.55} />);
      if (t.cdsStart != null && t.cdsEnd != null) {
        const cs = Math.max(s, t.cdsStart), ce = Math.min(en, t.cdsEnd);
        if (ce > cs) els.push(<rect key={`${key}d${e.rank}`} x={Math.min(X(cs), X(ce))} y={yy - 5.5} width={Math.max(1.2, Math.abs(X(ce) - X(cs)))} height={11} fill={color} rx={0.8} />);
      }
    }
    return els;
  };
  const torn = (x: number, yy: number, color: string, key: string) => <path key={key} d={`M${x},${yy - 8} l-3,3 l3,3 l-3,3 l3,3 l-3,3`} fill="none" stroke={color} strokeWidth={1.4} />;
  const panelDot = (sym: string) => { const pg = panel?.genes.get(sym.toUpperCase()); return pg ? CONF[pg.confidence].color : null; };

  // ticks of the reference ruler
  const ticks = useMemo(() => {
    const span = winLen / z, raw = span / 8, p = Math.pow(10, Math.floor(Math.log10(raw))), step = [1, 2, 5, 10].map(m => m * p).find(s => s >= raw) ?? p * 10;
    const out: number[] = []; for (let t = Math.ceil(win.start / step) * step; t <= win.end; t += step) out.push(t);
    return out;
  }, [winLen, z, win.start, win.end]);

  // ---- the SVG ----
  const svg = (
    <svg ref={svgRef} viewBox={`0 0 ${VW} ${H}`} width="100%" style={{ aspectRatio: `${VW} / ${H}`, cursor: z > 1 ? (drag.current ? 'grabbing' : 'grab') : 'default', fontFamily: FONT }} className="block select-none"
      onMouseDown={e => { if (e.button === 0) drag.current = { x: e.clientX, pan: view.pan }; }}
      onMouseMove={e => { const d = drag.current; if (!d || z <= 1) return; const r = svgRef.current!.getBoundingClientRect(); setView(v => ({ ...v, pan: d.pan + ((e.clientX - d.x) / r.width) * VW })); }}
      onMouseUp={() => { drag.current = null; }} onMouseLeave={() => { drag.current = null; leave(); }} onDoubleClick={() => setView(eventView())}>
      <defs>
        <clipPath id="svb-clip"><rect x={X0 - 2} y={0} width={X1 - X0 + 4} height={H} /></clipPath>
        <pattern id="svb-lost" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#fef2f2" /><line x1="0" y1="0" x2="0" y2="6" stroke={C.lost} strokeWidth="1.6" /></pattern>
        <linearGradient id="svb-bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#fbfcfd" /><stop offset="1" stopColor="#f5f8fa" /></linearGradient>
      </defs>
      <rect width={VW} height={H} fill="url(#svb-bg)" rx={12} />
      {/* track names */}
      <g fontSize={10.5} fill={C.muted}>
        {domains.length > 0 && <text x={14} y={yRefTad + TAD_H - 6}>TADs</text>}
        <text x={14} y={yRefGenes + 12} fontWeight={600} fill={C.text}>genes</text>
        <text x={14} y={yRefGenes + 25} fontSize={9.5}>{mode === 'mane' ? 'MANE Select' : 'all NM_'}</text>
        <text x={14} y={yRefBand + 14} fontWeight={700} fill={C.text}>{chrom}</text>
        <text x={14} y={yRefBand + 26} fontSize={9.5}>reference</text>
        {ders.map((d, i) => <g key={d.name}><text x={14} y={derY[i].band + 14} fontWeight={700} fill={C.accent}>{d.name}</text><text x={14} y={derY[i].band + 26} fontSize={9.5}>rebuilt</text></g>)}
      </g>
      <g clipPath="url(#svb-clip)">
        {/* reference TADs */}
        {domains.map((d, i) => {
          const a = refX(d.start), b = refX(d.end), h = Math.min(TAD_H - 6, (b - a) / 2);
          return <path key={`rt${i}`} d={`M${a},${yRefTad + TAD_H} L${(a + b) / 2},${yRefTad + TAD_H - h} L${b},${yRefTad + TAD_H} Z`} fill={i % 2 ? '#dfe5ee' : '#e9eef4'} stroke="#a9b6c6" strokeWidth={0.8}><title>{`domain ${d.name ?? i + 1}: ${chrom}:${fmt(d.start + 1)}-${fmt(d.end)} (${kb(d.end - d.start)})`}</title></path>;
        })}
        {/* breakpoints through the reference */}
        {rec.segs.filter(s => !s.first).map(s => <line key={`bpl${s.id}`} x1={refX(s.start)} x2={refX(s.start)} y1={domains.length ? yRefTad : yRefGenes} y2={yRefBand + BAND_H} stroke={C.lost} strokeWidth={1} strokeDasharray="3 3" opacity={0.55} />)}
        {/* reference genes */}
        {refItems.items.map(({ g, lane, compact }) => {
          const gi = impactOf.get(g.symbol); const col = IMPACT[gi?.kind ?? 'unaffected'].color;
          const yy = yRefGenes + refItems.laneY[lane];
          const lostGene = gi?.kind === 'deleted';
          const dot = panelDot(g.symbol), lx = Math.max(X0 + 2, refX(Math.max(g.start, win.start)));
          return (
            <g key={`rg${g.symbol}`} opacity={dim(g.symbol)} onMouseEnter={enter(g.symbol)} onMouseMove={enter(g.symbol)} onMouseLeave={leave} style={{ cursor: 'pointer' }} onClick={() => zoomTo(g.start, g.end)}>
              {dot && <circle cx={lx + 3.5} cy={yy + 6.5} r={3.5} fill={dot} />}
              {(!compact || worthLabel(gi)) && <text x={lx + (dot ? 10 : 0)} y={yy + 10} fontSize={11} fontWeight={gi && gi.kind !== 'unaffected' ? 700 : 500} fontStyle="italic" fill={col}>{g.symbol}</text>}
              {shown(g).map((t, i) => <g key={t.id}>{txGlyph(t, refX, yy + LABEL_H + 6 + i * TX_H, col, `r${t.id}`)}</g>)}
              {lostGene && <line x1={refX(g.start)} x2={refX(g.end)} y1={yy + LABEL_H + 6} y2={yy + LABEL_H + 6} stroke={C.lost} strokeWidth={2.2} opacity={0.6} />}
            </g>
          );
        })}
        {/* reference ruler and band */}
        {ticks.map(t => <g key={`tk${t}`}><line x1={refX(t)} x2={refX(t)} y1={yRefRuler + 12} y2={yRefRuler + 16} stroke={C.faint} /><text x={refX(t)} y={yRefRuler + 9} fontSize={8.5} textAnchor="middle" fill={C.faint} fontFamily={MONO}>{fmt(t)}</text></g>)}
        {rec.segs.map(s => {
          const a = refX(s.start), b = refX(s.end), col = segColors.get(s.id)!, lost = rec.lost.includes(s);
          return (
            <g key={`rb${s.id}`}>
              <rect x={a} y={yRefBand} width={Math.max(1, b - a)} height={BAND_H} fill={lost ? 'url(#svb-lost)' : col} fillOpacity={lost ? 1 : 0.85} stroke="#fff" strokeWidth={1} />
              {b - a > 14 && <text x={(a + b) / 2} y={yRefBand + 14} fontSize={11.5} fontWeight={700} textAnchor="middle" fill={lost ? C.lost : '#fff'}>{s.letter}</text>}
              {lost && b - a > 30 && <text x={(a + b) / 2} y={yRefBand + BAND_H + 14} fontSize={10} fontWeight={700} textAnchor="middle" fill={C.lost}>lost · {kb(s.end - s.start)}</text>}
            </g>
          );
        })}
        {/* ribbons */}
        {ders.map((d, di) => d.pieces.map(dp => {
          const col = segColors.get(dp.seg.id)!;
          const ra = refX(dp.seg.start), rb = refX(dp.seg.end);
          let da = derX(d, dp.derStart), db = derX(d, dp.derEnd);
          if (dp.piece.reversed) [da, db] = [db, da];
          const y0 = yRefBand + BAND_H, y1 = derY[di].band, m = (y0 + y1) / 2;
          const path = `M${ra},${y0} C${ra},${m} ${da},${m} ${da},${y1} L${db},${y1} C${db},${m} ${rb},${m} ${rb},${y0} Z`;
          return <path key={`rib${di}-${dp.index}`} d={path} fill={col} fillOpacity={0.16} stroke={col} strokeOpacity={0.45} strokeWidth={0.9}><title>{`${dp.seg.letter}${dp.piece.reversed ? '′ (inverted)' : ''}: ${chrom}:${fmt(dp.seg.start + 1)}-${fmt(dp.seg.end)} → ${d.name} ${kb(dp.derStart)}–${kb(dp.derEnd)}${dp.piece.copy ? ` (copy ${dp.piece.copy + 1})` : ''}`}</title></path>;
        }))}
        {/* rebuilt chromosomes */}
        {ders.map((d, di) => {
          const Y = derY[di];
          return (
            <g key={`der${d.name}`}>
              {d.pieces.map(dp => {
                const a = derX(d, dp.derStart), b = derX(d, dp.derEnd), col = segColors.get(dp.seg.id)!;
                return (
                  <g key={`db${dp.index}`}>
                    <rect x={a} y={Y.band} width={Math.max(1, b - a)} height={BAND_H} fill={col} fillOpacity={0.85} stroke="#fff" strokeWidth={1} />
                    {b - a > 18 && <text x={(a + b) / 2} y={Y.band + 14} fontSize={11.5} fontWeight={700} textAnchor="middle" fill="#fff">{dp.piece.reversed ? `◂ ${dp.seg.letter}′` : dp.seg.letter}</text>}
                  </g>
                );
              })}
              {d.pieces.slice(1).map((dp, i) => {
                const x = derX(d, dp.derStart), via = dp.piece.via, alt = i % 2;
                if (!via || via === 'ref') return <line key={`js${dp.index}`} x1={x} x2={x} y1={Y.band} y2={Y.band + BAND_H} stroke="#fff" strokeWidth={1.5} />;
                const col = KIND_COLOR[via.kind];
                return (
                  <g key={`js${dp.index}`}>
                    <line x1={x} x2={x} y1={Y.band - 3} y2={Y.genes + derPacked[di].height} stroke={col} strokeWidth={1} strokeDasharray="2 3" opacity={0.6} />
                    <circle cx={x} cy={Y.band + BAND_H / 2} r={5} fill={col} stroke="#fff" strokeWidth={1.8} />
                    <text x={x} y={Y.ruler + (alt ? 16 : 8)} fontSize={9.5} fontWeight={700} textAnchor="middle" fill={col}>{via.id} · {fmt(via.reads)}</text>
                  </g>
                );
              })}
              {/* rebuilt genes */}
              {derPacked[di].items.map((b, bi) => {
                const yy = Y.genes + derPacked[di].laneY[b.lane];
                const gi = impactOf.get(b.genes[0]);
                const isFusion = b.genes.length > 1;
                const head = isFusion ? C.fusion5 : IMPACT[gi?.kind ?? 'unaffected'].color;
                const lx = Math.max(X0 + 2, b.x0);
                const dot = !isFusion ? panelDot(b.genes[0]) : null;
                return (
                  <g key={`dg${bi}`} opacity={dim(b.genes)} onMouseEnter={enter(b.genes[0])} onMouseMove={enter(b.genes[0])} onMouseLeave={leave} style={{ cursor: 'pointer' }}>
                    {dot && <circle cx={lx + 3.5} cy={yy + 6.5} r={3.5} fill={dot} />}
                    {(!b.compact || isFusion || worthLabel(gi)) && <text x={lx + (dot ? 10 : 0)} y={yy + 10} fontSize={11} fontWeight={700} fontStyle="italic" fill={head}>{b.label}</text>}
                    {isFusion && b.fusion && <g>
                      <rect x={lx + labelW(b.label) - 4} y={yy + 1} width={b.fusion.kind === 'in-frame' ? 52 : b.fusion.kind === 'promoter-swap' ? 74 : 66} height={12} rx={6} fill={b.fusion.kind === 'in-frame' ? '#7c3aed' : '#ede9fe'} />
                      <text x={lx + labelW(b.label) + 2} y={yy + 10} fontSize={9} fontWeight={700} fill={b.fusion.kind === 'in-frame' ? '#fff' : '#6d28d9'}>{b.fusion.kind === 'in-frame' ? 'in frame' : b.fusion.kind === 'promoter-swap' ? 'promoter swap' : b.fusion.kind === 'out-of-frame' ? 'out of frame' : b.fusion.kind === 'five-cds-whole' ? '3′ UTR join' : 'frame ?'}</text>
                    </g>}
                    {b.chains.map((c, ci) => {
                      const yl = yy + LABEL_H + 6 + ci * TX_H;
                      return (
                        <g key={`ch${ci}`}>
                          {c.parts.map((p: Fragment, pi) => {
                            const col = isFusion ? (p.tx.gene === b.genes[0] ? C.fusion5 : C.fusion3) : head;
                            const X = (pos: number) => derX(d, toDer(p.dp, pos));
                            const els = txGlyph(p.tx, X, yl, col, `d${bi}${ci}${pi}`, p.lo, p.hi, p.strand);
                            // torn ends: where the chain stops at a breakpoint
                            const leftX = Math.min(X(p.lo), X(p.hi)), rightX = Math.max(X(p.lo), X(p.hi));
                            const cutL = p.dp.piece.reversed ? p.cutHi && !p.edgeHi : p.cutLo && !p.edgeLo, cutR = p.dp.piece.reversed ? p.cutLo && !p.edgeLo : p.cutHi && !p.edgeHi;
                            const isFirst = (c.parts[0].strand > 0 ? pi === 0 : pi === c.parts.length - 1), isLast = (c.parts[0].strand > 0 ? pi === c.parts.length - 1 : pi === 0);
                            if (cutL && (c.parts[0].strand > 0 ? isFirst : isLast)) els.push(torn(leftX, yl, C.lost, `tl${bi}${ci}${pi}`));
                            if (cutR && (c.parts[0].strand > 0 ? isLast : isFirst)) els.push(torn(rightX, yl, C.lost, `tr${bi}${ci}${pi}`));
                            return <g key={`p${pi}`}>{els}</g>;
                          })}
                          {isFusion && (() => {
                            const xs = c.parts.flatMap(p => [derX(d, toDer(p.dp, p.lo)), derX(d, toDer(p.dp, p.hi))]);
                            return <rect x={Math.min(...xs) - 3} y={yl - 6.5} width={Math.max(...xs) - Math.min(...xs) + 6} height={13} rx={6.5} fill="none" stroke={C.fusion5} strokeWidth={1.1} strokeDasharray="4 2" />;
                          })()}
                        </g>
                      );
                    })}
                  </g>
                );
              })}
              {/* rebuilt TADs, hanging under the chromosome; new domains outlined */}
              {domains.length > 0 && story.derDomains.filter(x => story.ders[x.der] === d).map((x, i) => {
                const a = derX(d, x.derStart), b = derX(d, x.derEnd), h = Math.min(TAD_H - 8, (b - a) / 2);
                const top = Y.tad + 2;
                return (
                  <g key={`dt${i}`}>
                    <path d={`M${a},${top} L${(a + b) / 2},${top + h} L${b},${top} Z`} fill={x.neo ? '#ede9fe' : i % 2 ? '#dfe5ee' : '#e9eef4'} stroke={x.neo ? '#7c3aed' : '#a9b6c6'} strokeWidth={x.neo ? 1.4 : 0.8} strokeDasharray={x.neo ? '4 2' : undefined}>
                      <title>{x.neo ? `new domain: ${x.parts.map(p => `${kb(p.hi - p.lo)} of domain ${p.domain + 1}`).join(' + ')}` : `domain ${x.parts[0].domain + 1}${x.parts[0].hi - x.parts[0].lo < domains[x.parts[0].domain].end - domains[x.parts[0].domain].start ? ' (part)' : ''}`}</title>
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
    </svg>
  );

  // ---- gene table ----
  const nearby = story.genes.filter(g => g.kind === 'nearby');
  const touched = [...story.genes.filter(g => g.kind !== 'unaffected' && g.kind !== 'nearby'), ...nearby.filter(g => notable(g, annotations))];
  const nearOthers = nearby.filter(g => !notable(g, annotations));
  const unchanged = story.genes.filter(g => g.kind === 'unaffected');
  const rowsOf = (list: GeneImpact[]) => [...list].sort((a, b) => b.rank - a.rank || (panel?.genes.get(b.gene.symbol.toUpperCase())?.confidence ?? 0) - (panel?.genes.get(a.gene.symbol.toUpperCase())?.confidence ?? 0) || (ann.dosage.get(b.gene.symbol)?.hi ?? -1) - (ann.dosage.get(a.gene.symbol)?.hi ?? -1));
  const findings = useMemo(() => summary(story.genes, annotations, near), [story, annotations, near]);

  return (
    <div className="mt-2 border-t border-gray-200 pt-3" style={{ fontFamily: FONT }}>
      {/* header */}
      <div className="px-5 flex flex-wrap items-center gap-2 text-xs">
        <div className="text-sm font-bold mr-1">Genes on the rebuilt chromosome</div>
        <span className="inline-flex rounded-md border border-gray-300 overflow-hidden" role="group" aria-label="Transcripts">
          {(['mane', 'all'] as const).map(m => <button key={m} onClick={() => setMode(m)} className={`px-2 py-0.5 ${mode === m ? 'bg-accent-600 text-white' : 'bg-white text-gray-600 hover:bg-gray-50'}`} title={m === 'mane' ? 'One model per gene: MANE Select (else RefSeq Select, else the longest coding sequence)' : 'Every RefSeq NM_ model of each gene (no XM_ prediction, no NR_ non-coding): the consequence can differ between isoforms'}>{m === 'mane' ? 'MANE Select' : 'all NM_ transcripts'}</button>)}
        </span>
        <PanelPicker panel={panelData} onPanel={p => { setPanelData(p); store.set('sashimi.svPanel', p); }} />
        <label className="inline-flex items-center gap-1 px-2 py-0.5 rounded-md border border-gray-300 bg-white hover:bg-gray-50 cursor-pointer" title="Chromatin domains (TADs) as a BED file (chrom, start, end, name; 0-based, the same genome build as the alignments), e.g. from Hi-C of a relevant tissue. They show which genes change domain: an intact gene joined to another domain may be misregulated (position effect).">
          ◭ {domainsFile ? `TADs: ${domainsFile.name} (${domains.length} here)` : 'Add chromatin domains (BED)…'}
          <input type="file" accept=".bed,.txt,.tsv" className="hidden" onChange={async e => { const f = e.target.files?.[0]; if (f) { const v = { name: f.name, text: await f.text() }; setDomainsFile(v); store.set('sashimi.svDomains', v.text.length < 3e6 ? v : null); } e.target.value = ''; }} />
        </label>
        {domainsFile && <button onClick={() => { setDomainsFile(null); store.set('sashimi.svDomains', null); }} className="text-gray-400 hover:text-red-600" title="Remove the domains">×</button>}
        <label className="inline-flex items-center gap-1 text-gray-600" title={`Without chromatin domains, an intact gene within this distance of a breakpoint is reported as a possible position effect, and the genes that far around the event are fetched (zoom out to see them). TADs average about 1 Mb; position effects have been documented 1.3 Mb from a breakpoint (SOX9) and further.${domains.length ? ' With domains loaded, the domains are used instead.' : ''}`}>
          near a breakpoint ≤
          <select value={near} onChange={e => { const v = Number(e.target.value); setNear(v); store.set('sashimi.svNear', v); }} className="border rounded px-1 py-0.5">
            {[0, 250e3, 500e3, 1e6, 2e6, 3e6].map(v => <option key={v} value={v}>{v === 0 ? 'off' : v >= 1e6 ? `${v / 1e6} Mb` : `${v / 1e3} kb`}</option>)}
          </select>
          {domains.length > 0 && <span className="text-gray-400">(TADs used instead)</span>}
        </label>
        <span className="ml-auto flex items-center gap-1">
          {[['−', () => setView(v => ({ z: Math.max(1, v.z / 1.5), pan: v.z / 1.5 <= 1 ? 0 : v.pan / 1.5 }))], ['+', () => setView(v => ({ z: Math.min(200, v.z * 1.5), pan: v.pan * 1.5 - ((X0 + X1) / 2 - X0) * 0.5 }))], ['event', () => setView(eventView())], [near > 0 ? `±${near >= 1e6 ? `${near / 1e6} Mb` : `${near / 1e3} kb`}` : 'all', () => setView({ z: 1, pan: 0 })]].map(([l, f]) => (
            <button key={l as string} onClick={f as () => void} className="min-w-[28px] px-1.5 h-6 rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 font-semibold" title={l === 'event' ? 'Fit the event (double-click on the figure)' : l === '−' || l === '+' ? 'Zoom (or Ctrl + wheel)' : 'The whole region: the event and the genes around it'}>{l as string}</button>
          ))}
        </span>
      </div>
      {/* key findings */}
      {findings.length > 0 && (
        <div className="px-5 mt-2 flex flex-wrap gap-1.5">
          {findings.map((f, i) => (
            <span key={i} className="inline-flex items-center gap-1.5 rounded-full px-2.5 py-1 text-[11.5px] font-medium" style={{ background: `${f.color}14`, color: f.color, border: `1px solid ${f.color}40` }}>
              <b>{f.title}</b>{f.genes.map(g => <span key={g.symbol} className="italic font-semibold cursor-pointer hover:underline" onMouseEnter={() => setHover(g.symbol)} onMouseLeave={() => setHover(null)}>{g.symbol}{g.note ? <span className="not-italic font-normal opacity-80"> ({g.note})</span> : null}</span>).reduce((acc: React.ReactNode[], el, j) => (j ? [...acc, ', ', el] : [el]), [])}
            </span>
          ))}
        </div>
      )}
      <div className="px-3 mt-2 relative" ref={boxRef}>
        {txs.loading ? <div className="h-[220px] flex items-center justify-center text-sm text-gray-500">loading the gene models of the window…</div> : svg}
        {tip && (() => {
          const gi = impactOf.get(tip.gene); if (!gi) return null;
          const r = dosageReading(gi, annotations);
          return (
            <div className="pointer-events-none absolute z-10 rounded-md bg-gray-900/92 text-white text-[11px] px-2.5 py-1.5 shadow-lg max-w-[380px]" style={{ left: Math.min(tip.x + 14, (boxRef.current?.clientWidth ?? 900) - 390), top: tip.y + 14 }}>
              <div className="font-semibold"><i>{gi.gene.symbol}</i> · {IMPACT[gi.kind].label}</div>
              <div className="text-gray-300">{gi.principal.tx.id}{gi.principal.tx.status === 'mane' ? ' (MANE Select)' : gi.principal.tx.status === 'select' ? ' (RefSeq Select)' : ''}</div>
              <div className="mt-0.5">{gi.principal.text}</div>
              {r.text && <div className="mt-0.5 text-amber-200">{r.text}</div>}
            </div>
          );
        })()}
      </div>
      {/* legend */}
      <div className="px-5 mt-1 flex flex-wrap items-center gap-x-3 gap-y-1 text-[10.5px] text-gray-500">
        {(['fusion', 'disrupted', 'deleted', 'intragenic-del', 'intragenic-dup', 'partial-copy', 'duplicated', 'domain', 'moved', 'nearby', 'unaffected'] as ImpactKind[]).map(k2 => <span key={k2} className="inline-flex items-center gap-1"><span className="inline-block w-3 h-2.5 rounded-sm" style={{ background: IMPACT[k2].color }} />{IMPACT[k2].label}</span>)}
        <span className="inline-flex items-center gap-1"><svg width="10" height="14"><path d="M6,1 l-3,3 l3,3 l-3,3 l3,3" fill="none" stroke={C.lost} strokeWidth="1.4" /></svg>cut by a breakpoint</span>
        {panel && <span className="inline-flex items-center gap-1">● panel rating: <span style={{ color: CONF[3].color }}>green</span> <span style={{ color: CONF[2].color }}>amber</span> <span style={{ color: CONF[1].color }}>red</span></span>}
        <span>· scroll with Ctrl to zoom, drag to pan, click a gene to zoom on it, double-click to fit</span>
      </div>
      {/* table */}
      <div className="px-5 mt-3 pb-1">
        <table className="w-full border-collapse text-xs">
          <thead><tr className="text-[10px] uppercase tracking-wide text-gray-500 border-b border-gray-200">
            <th className="text-left font-semibold py-1 pr-2">gene</th><th className="text-left font-semibold pr-2">impact</th><th className="text-left font-semibold pr-2">consequence</th>
            <th className="text-left font-semibold pr-2">ClinGen dosage</th><th className="text-left font-semibold pr-2">gnomAD LoF constraint</th><th className="text-left font-semibold">{panel ? 'panel' : 'PanelApp'}</th>
          </tr></thead>
          <tbody>
            {rowsOf(touched).map(g => <GeneRow key={g.gene.symbol} gi={g} ann={annotations} mode={mode} hover={hover === g.gene.symbol} onHover={setHover} onZoom={() => zoomTo(g.gene.start, g.gene.end)} />)}
            {showNear && rowsOf(nearOthers).map(g => <GeneRow key={g.gene.symbol} gi={g} ann={annotations} mode={mode} hover={hover === g.gene.symbol} onHover={setHover} onZoom={() => zoomTo(g.gene.start, g.gene.end)} />)}
            {showUnchanged && rowsOf(unchanged).map(g => <GeneRow key={g.gene.symbol} gi={g} ann={annotations} mode={mode} hover={hover === g.gene.symbol} onHover={setHover} onZoom={() => zoomTo(g.gene.start, g.gene.end)} />)}
          </tbody>
        </table>
        <div className="flex flex-wrap gap-x-4">
          {nearOthers.length > 0 && <button onClick={() => setShowNear(v => !v)} className="mt-1 text-[11px] text-accent-700 hover:underline">{showNear ? 'hide' : 'show'} the {nearOthers.length} other gene{nearOthers.length > 1 ? 's' : ''} near a breakpoint (not in the panel, not dosage-sensitive, not LoF-intolerant)</button>}
          {unchanged.length > 0 && <button onClick={() => setShowUnchanged(v => !v)} className="mt-1 text-[11px] text-accent-700 hover:underline">{showUnchanged ? 'hide' : 'show'} the {unchanged.length} unchanged gene{unchanged.length > 1 ? 's' : ''}</button>}
        </div>
        {!touched.length && !txs.loading && <div className="text-gray-400 text-xs py-2">No gene of the window is changed by the selected junctions.</div>}
        <div className="text-[10.5px] text-gray-500 mt-2 leading-4">
          Genes: {txs.source || '—'}{txs.error ? <span className="text-amber-700"> · {txs.error}</span> : null}. {ann.loading ? 'Annotations loading… ' : ann.notes.length ? `${ann.notes.join(' · ')}. ` : ''}
          Copies count the {zygosity === 'het' ? `${flankCopies - 1} intact homologue${flankCopies - 1 === 1 ? '' : 's'} and the rearranged one` : 'rearranged homologues'}; a fusion needs both partners on the same strand of the rebuilt chromosome, its frame is read from the coding phases of the two exons joined (breakpoints in introns); a breakpoint inside an exon leaves the frame undetermined.
          LoF-intolerant: LOEUF &lt; {LOEUF_INTOLERANT.v4} (gnomAD v4; &lt; {LOEUF_INTOLERANT['v2.1.1']} for v2.1.1) or pLI ≥ {PLI_INTOLERANT}, as gnomAD suggests for Mendelian cases when a hard threshold is needed. ClinGen scores: 3 sufficient, 2 some, 1 little, 0 no evidence, 30 autosomal recessive, 40 dosage sensitivity unlikely.
          {domains.length ? ' Domains: a junction joining two parts of domains makes one new domain (boundaries between them lost); genes whose domain gains foreign sequence or loses part of its own are flagged, whatever the size.' : near > 0 ? ` No domains loaded: intact genes within ${near >= 1e6 ? `${near / 1e6} Mb` : `${near / 1e3} kb`} of a breakpoint are listed as a possible position effect; those in the panel, with a ClinGen score of 2 or 3, or LoF-intolerant are shown first, the others folded.` : ''}
        </div>
      </div>
    </div>
  );
}

// ======================== Genes worth a look ========================

/** In the chosen panel, a ClinGen haploinsufficiency or triplosensitivity score of 2 or 3, or LoF-intolerant in gnomAD. */
function notable(g: GeneImpact, a: GeneAnnotations): boolean {
  const d = a.dosage.get(g.gene.symbol), c = a.constraint.get(g.gene.symbol);
  return !!a.panel?.genes.has(g.gene.symbol.toUpperCase()) || [2, 3].includes(d?.hi ?? -1) || [2, 3].includes(d?.ts ?? -1)
    || (c?.loeuf != null && c.loeuf < (LOEUF_INTOLERANT[c.version] ?? LOEUF_INTOLERANT.v4)) || (c?.pli != null && c.pli >= PLI_INTOLERANT);
}

// ======================== Key findings ========================

function summary(genes: GeneImpact[], a: GeneAnnotations, near: number): { title: string; color: string; genes: { symbol: string; note?: string }[] }[] {
  const note = (g: GeneImpact) => {
    const d = a.dosage.get(g.gene.symbol), c = a.constraint.get(g.gene.symbol), p = a.panel?.genes.get(g.gene.symbol.toUpperCase());
    const positional = g.kind === 'nearby' || g.kind === 'domain' || g.kind === 'moved';
    return [p ? `${CONF[p.confidence].label} in panel` : '',
      d?.hi != null && [2, 3].includes(d.hi) && (g.principal.copies < g.principal.refCopies || positional) ? `HI ${d.hi}` : '',
      d?.ts != null && [2, 3].includes(d.ts) && (g.kind === 'duplicated' || positional) ? `TS ${d.ts}` : '',
      c?.loeuf != null && c.loeuf < (LOEUF_INTOLERANT[c.version] ?? 0.45) ? `LOEUF ${c.loeuf.toFixed(2)}` : ''].filter(Boolean).join(', ') || undefined;
  };
  const out: { title: string; color: string; genes: { symbol: string; note?: string }[] }[] = [];
  const add = (title: string, color: string, list: GeneImpact[]) => { if (list.length) out.push({ title, color, genes: list.map(g => ({ symbol: g.gene.symbol, note: note(g) })) }); };
  add('Loss of a functional copy:', '#b91c1c', genes.filter(g => g.kind === 'deleted' || g.kind === 'disrupted' || (g.kind === 'intragenic-del' && g.principal.frame !== 'kept') || (g.kind === 'intragenic-dup' && g.principal.frame === 'shifted')));
  const fz = new Map<string, string>();
  for (const g of genes) for (const f of g.principal.fusions) fz.set(`${f.five.gene}::${f.three.gene}`, f.kind === 'in-frame' ? 'in frame' : f.kind === 'promoter-swap' ? 'promoter swap' : f.kind === 'out-of-frame' ? 'out of frame' : f.kind === 'five-cds-whole' ? '3′ UTR' : 'frame undetermined');
  if (fz.size) out.push({ title: 'Fusion:', color: '#7c3aed', genes: [...fz].map(([symbol, n]) => ({ symbol, note: n })) });
  add('In-frame change:', '#ea580c', genes.filter(g => (g.kind === 'intragenic-del' || g.kind === 'intragenic-dup') && g.principal.frame === 'kept'));
  add('Gain of a whole copy:', '#15803d', genes.filter(g => g.kind === 'duplicated'));
  add('Chromatin domain changed:', '#b45309', genes.filter(g => g.kind === 'domain'));
  add('Intact, inverted or moved:', '#2563eb', genes.filter(g => g.kind === 'moved'));
  const nb = genes.filter(g => g.kind === 'nearby');
  if (nb.length) {
    const hot = nb.filter(g => notable(g, a));
    out.push({ title: `Near a breakpoint (≤ ${near >= 1e6 ? `${near / 1e6} Mb` : `${near / 1e3} kb`}), ${nb.length} gene${nb.length > 1 ? 's' : ''}${hot.length ? ', of which:' : ': none in the panel, dosage-sensitive or LoF-intolerant'}`, color: IMPACT.nearby.color, genes: hot.map(g => ({ symbol: g.gene.symbol, note: [note(g), g.nearest ? `${Math.round(g.nearest.distance / 1000)} kb` : ''].filter(Boolean).join(', ') })) });
  }
  return out;
}

// ======================== One table row ========================

function GeneRow({ gi, ann, mode, hover, onHover, onZoom }: { gi: GeneImpact; ann: GeneAnnotations; mode: 'mane' | 'all'; hover: boolean; onHover: (s: string | null) => void; onZoom: () => void }) {
  const sym = gi.gene.symbol, d = ann.dosage.get(sym), c = ann.constraint.get(sym), p = ann.panel?.genes.get(sym.toUpperCase());
  const r = dosageReading(gi, ann);
  const thr = c ? LOEUF_INTOLERANT[c.version] ?? 0.45 : 0.45;
  const scoreChip = (label: string, v: number | null | undefined, hot: string) => (
    <span className="inline-block px-1.5 rounded text-[10.5px] font-semibold mr-1" title={`${label}: ${dosageScoreText(v)}`}
      style={v === 3 ? { background: `${hot}1f`, color: hot } : v === 2 ? { background: '#fef3c7', color: '#b45309' } : { background: '#f1f4f8', color: '#5b6b7e' }}>{label} {v ?? '—'}</span>
  );
  return (
    <tr className={`border-b border-gray-100 align-top ${hover ? 'bg-accent-50' : ''}`} onMouseEnter={() => onHover(sym)} onMouseLeave={() => onHover(null)}>
      <td className="py-1.5 pr-2 whitespace-nowrap">
        <button onClick={onZoom} className="italic font-bold text-[12.5px] hover:underline" style={{ color: IMPACT[gi.kind].color }} title="Zoom on this gene">{sym}</button>
        <div className="text-[10px] text-gray-400 font-mono">{gi.principal.tx.id}{gi.principal.tx.status === 'mane' ? ' · MANE' : gi.principal.tx.status === 'select' ? ' · Select' : ''}</div>
      </td>
      <td className="py-1.5 pr-2"><span className="inline-block px-2 py-0.5 rounded-full text-[10.5px] font-semibold whitespace-nowrap" style={{ background: `${IMPACT[gi.kind].color}18`, color: IMPACT[gi.kind].color }}>{IMPACT[gi.kind].short}</span>
        <div className="text-[10px] text-gray-500 mt-0.5 font-mono">{gi.principal.copies} / {gi.principal.refCopies} copies</div></td>
      <td className="py-1.5 pr-2 text-gray-700 leading-4 max-w-[430px]">
        {gi.principal.text}
        {mode === 'all' && gi.others.length > 0 && <div className="mt-0.5 text-[10.5px] text-gray-500">other NM_: {gi.others.map(o => `${o.tx.id} ${IMPACT[o.kind].short}${o.frame === 'kept' || o.frame === 'shifted' ? `, frame ${o.frame}` : ''}${o.lost.length ? ` (${ranksText(o.lost)} lost)` : ''}`).join(' · ')}</div>}
        {mode === 'mane' && gi.others.length > 0 && <div className="mt-0.5 text-[10.5px] text-amber-700">{gi.others.length} other NM_ model{gi.others.length > 1 ? 's' : ''} read{gi.others.length > 1 ? '' : 's'} differently: switch to all NM_ transcripts</div>}
      </td>
      <td className="py-1.5 pr-2 whitespace-nowrap">
        {scoreChip('HI', d?.hi, '#b91c1c')}{scoreChip('TS', d?.ts, '#1d4ed8')}
        {r.text && <div className={`text-[10.5px] mt-0.5 whitespace-normal max-w-[220px] ${r.level === 'high' ? 'text-red-700 font-semibold' : r.level === 'moderate' ? 'text-amber-700' : 'text-gray-500'}`}>{r.text}</div>}
      </td>
      <td className="py-1.5 pr-2 whitespace-nowrap">
        {c ? (
          <div className="flex items-center gap-1.5" title={`gnomAD ${c.version}: LOEUF (upper bound of the observed / expected LoF ratio) ${c.loeuf ?? 'n/a'}, pLI ${c.pli ?? 'n/a'}; LoF-intolerant below LOEUF ${thr} (or pLI ≥ ${PLI_INTOLERANT})`}>
            <svg width="62" height="10"><rect x={0} y={3} width={62} height={4} rx={2} fill="#e3e8ee" /><line x1={(thr / 2) * 62} x2={(thr / 2) * 62} y1={0} y2={10} stroke="#94a3b8" strokeWidth={1} />{c.loeuf != null && <circle cx={Math.min(62, (c.loeuf / 2) * 62)} cy={5} r={3.5} fill={c.loeuf < thr ? '#b91c1c' : '#64748b'} />}</svg>
            <span className="font-mono text-[10.5px]">LOEUF {c.loeuf != null ? c.loeuf.toFixed(2) : '—'}</span>
            <span className="font-mono text-[10.5px] text-gray-500">pLI {c.pli != null ? c.pli.toFixed(2) : '—'}</span>
          </div>
        ) : <span className="text-gray-400">—</span>}
        {c && <div className="text-[10px] text-gray-400">gnomAD {c.version}</div>}
      </td>
      <td className="py-1.5 whitespace-nowrap">
        {p ? <span className="inline-block px-1.5 py-0.5 rounded text-[10.5px] font-semibold" style={{ background: CONF[p.confidence].bg, color: CONF[p.confidence].color }} title={p.phenotypes?.join('; ')}>{CONF[p.confidence].label}{p.moi ? ` · ${p.moi.replace(/^(\w+ ?\w*).*/, '$1').toLowerCase()}` : ''}</span> : <span className="text-gray-300">{ann.panel ? 'not in panel' : '—'}</span>}
      </td>
    </tr>
  );
}

// ======================== PanelApp picker ========================

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
