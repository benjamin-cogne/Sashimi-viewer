/**
 * Structural-variant cartoon: the junction arcs the user picked on a DNA track are turned into the rearranged
 * chromosome they imply (svModel.ts), drawn as an animated SVG scene: the reference cut at the breakpoints, every
 * segment flying to its place on the rebuilt chromosome (flipped when inverted, copied when duplicated), the segments
 * on no chromosome struck out as lost, then the copy number each segment should have against the read depth. Pure SVG
 * on requestAnimationFrame, as the splicing cartoon.
 */
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { CoverageRun } from './types';
import type { TxModel } from './geometry';
import { formatBp } from './geometry';
import { reconstruct, reconstructAll, rankPaths, SV_BREAKEND_CLUSTER_BP, pathText, type RankedPath, checkDepth, segmentContent, breakpointSite, exonRun, SV_DEPTH_MARGIN_BP, SV_DEPTH_MIN_BINS, SV_FIT_MIN_COPIES, SV_FIT_Z, SV_MOSAIC_MIN, SV_MOSAIC_MAX, SV_FLANK_MIN_BINS,
  type SvJunction, type SvSegment, type SvPiece, type SegmentDepth } from './svModel';
import { KIND_COLOR, segmentColors } from './svColors';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import SvNomenclatureBox, { useCytobands } from './SvNomenclatureBox';
import type { Cytoband } from './svNomenclature';

// ======================== Timeline ========================

const T_BREAK = 1.0, T_OPEN = 1.0, T_J0 = 2.5, STEP = 1.45, MOVE = 1.15;
const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
const easeInOut = (p: number) => (p < 0.5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2);
const phase = (t: number, a: number, b: number) => easeInOut(clamp01((t - a) / (b - a)));
const lerp = (a: number, b: number, p: number) => a + (b - a) * p;

// ======================== Palette and layout ========================

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const MONO = 'ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
const C = { text: '#1f2937', muted: '#5b6b7e', faint: '#8796a8', grid: '#e3e8ee', lost: '#dc2626', ok: '#16a34a', accent: '#0f766e', dup: '#15803d', inv: '#2563eb', bg: '#ffffff' };
const VW = 1200, LX = 120, RX = 1170;
const SEG_H = 30, GAP_OPEN = 26, GAP_JOIN = 8;
const fmt = (n: number) => n.toLocaleString('en-US');
const luminance = (hex: string) => { const n = parseInt(hex.slice(1), 16); return (0.299 * (n >> 16) + 0.587 * ((n >> 8) & 255) + 0.114 * (n & 255)) / 255; };
/** Widths grow with the log of the length so a 2 kb and a 2 Mb segment both stay readable. */
const widthOf = (s: SvSegment) => (s.first || s.last ? 92 : 46 + 44 * Math.log10(1 + (s.end - s.start) / 8000));

interface PlacedSeg extends SvSegment { color: string; w: number; xc: number; xo: number; y: number }
interface PlacedPiece extends SvPiece { x: number; y: number; t0: number }
interface Scene {
  segs: PlacedSeg[]; byChr: Record<string, PlacedSeg[]>; refRows: { chrom: string; y: number }[];
  derRows: { name: string; ends: string; rearranged: boolean; pieces: PlacedPiece[]; y: number }[];
  lost: PlacedSeg[]; H: number; tLost: number; tCn: number; tEnd: number;
}

function layout(rec: ReturnType<typeof reconstruct>): Scene {
  const segs = rec.segs as PlacedSeg[];
  const colors = segmentColors(rec);
  for (const s of segs) s.color = colors.get(s.id)!;
  const refRows = rec.chrs.map((chrom, k) => ({ chrom, y: 84 + k * 100 }));
  const refBottom = refRows[refRows.length - 1].y;
  const derRows = rec.ders.map((d, k) => ({ name: d.name, ends: d.ends, rearranged: d.rearranged, pieces: d.pieces as PlacedPiece[], y: refBottom + 134 + k * 96 }));
  const H = (derRows.length ? derRows[derRows.length - 1].y : refBottom + 40) + 70;
  const w0 = new Map(segs.map(s => [s.id, widthOf(s)]));
  const need = Math.max(
    ...refRows.map(rr => rec.byChr[rr.chrom].reduce((t, s) => t + w0.get(s.id)!, 0) + GAP_OPEN * (rec.byChr[rr.chrom].length - 1)),
    ...derRows.map(dr => dr.pieces.reduce((t, p) => t + w0.get(p.seg.id)!, 0) + GAP_JOIN * (dr.pieces.length - 1)));
  const k = Math.min(1, (RX - LX) / Math.max(1, need));
  for (const s of segs) s.w = w0.get(s.id)! * k;
  for (const rr of refRows) {
    const list = rec.byChr[rr.chrom] as PlacedSeg[];
    const closed = list.reduce((t, s) => t + s.w, 0), open = closed + GAP_OPEN * k * (list.length - 1);
    let xc = LX + (RX - LX - open) / 2 + (open - closed) / 2, xo = LX + (RX - LX - open) / 2;
    for (const s of list) { s.xc = xc; s.xo = xo; s.y = rr.y; xc += s.w; xo += s.w + GAP_OPEN * k; }
  }
  let order = 0;
  for (const dr of derRows) {
    const tot = dr.pieces.reduce((t, p) => t + (p.seg as PlacedSeg).w, 0) + GAP_JOIN * (dr.pieces.length - 1);
    let x = LX + (RX - LX - tot) / 2;
    for (const p of dr.pieces) { p.x = x; p.y = dr.y; x += (p.seg as PlacedSeg).w + GAP_JOIN; p.t0 = T_J0 + order * STEP; order++; }
  }
  const tLost = T_J0 + order * STEP + 0.3, tCn = tLost + 1.3, tEnd = tCn + 1.6;
  return { segs, byChr: rec.byChr as Record<string, PlacedSeg[]>, refRows, derRows, lost: rec.lost as PlacedSeg[], H, tLost, tCn, tEnd };
}

// ======================== Component ========================

export interface SvCartoonProps {
  chrom: string;
  /** the stretch the sample's coverage covers (0-based half-open): the chromosome is drawn from it */
  window: { start: number; end: number };
  /** every structural arc of the sample in the window, as junctions; the user ticks the ones of the event */
  candidates: SvJunction[];
  initialSelected: string[];
  sample: { name: string; color: string; coverage: CoverageRun[]; /** depths are scaled estimates (sampled window) */ sampled?: boolean };
  /** the other DNA samples shown, for the reference level */
  controls: { name: string; coverage: CoverageRun[] }[];
  models: TxModel[];
  build: string;
  /** the genes of the window and their annotations (the data source's), for the rearranged genome browser */
  geneLoaders: Omit<SvGeneLoaders, 'fallback'>;
  onClose: () => void;
}

export default function SvCartoon({ chrom, window: win, candidates, initialSelected, sample, controls, models, build, geneLoaders, onClose }: SvCartoonProps) {
  const [selected, setSelected] = useState<Set<string>>(() => new Set(initialSelected));
  const [zygosity, setZygosity] = useState<'het' | 'hom'>('het');
  const [copies, setCopies] = useState<number>(/^(chr)?[XY]$/i.test(chrom) ? 1 : 2);
  const [useControls, setUseControls] = useState(true);
  const [t, setT] = useState(0);
  const [playing, setPlaying] = useState(true);
  const [speed, setSpeed] = useState(1);
  const [explore, setExplore] = useState(false);
  const [view, setView] = useState({ k: 1, x: 0, y: 0 });
  const svgRef = useRef<SVGSVGElement>(null);
  const sceneRef = useRef<HTMLDivElement>(null);
  const panRef = useRef<{ x: number; y: number; vx: number; vy: number; moved: boolean } | null>(null);
  const raf = useRef<number | null>(null);
  const last = useRef<number | null>(null);
  const tRef = useRef(0);
  tRef.current = t;

  const chosen = useMemo(() => candidates.filter(j => selected.has(j.id)), [candidates, selected]);
  // several junctions leaving one breakend give several paths: ranked by the depth, the best shown unless one is picked
  const [pathKey, setPathKey] = useState<string | null>(null);
  useEffect(() => setPathKey(null), [chosen]);
  const ranked = useMemo(() => rankPaths(reconstructAll([{ chrom, start: win.start, end: win.end }], chosen), { name: sample.name, coverage: sample.coverage }, useControls ? controls : [], { zygosity, flankCopies: copies }),
    [chosen, chrom, win.start, win.end, sample, controls, useControls, zygosity, copies]);
  const story = useMemo(() => {
    const { rec, depth } = ranked.find(r => r.rec.choice.join(',') === pathKey) ?? ranked[0];
    const scene = layout(rec);
    const depthOf = new Map(depth.segments.map(d => [d.seg.id, d]));
    return { rec, depth, scene, depthOf };
  }, [ranked, pathKey, build]);
  const { scene, rec, depth } = story;
  const bands = useCytobands(chrom, geneLoaders.getCytobands);
  const segColorMap = useMemo(() => segmentColors(rec), [rec]);
  const loaders = useMemo<SvGeneLoaders>(() => ({ ...geneLoaders, fallback: models.map(m => ({
    id: m.transcriptId, gene: m.geneName, strand: m.strand, start: m.start, end: m.end, exons: m.exons.map(e => ({ start: e.start, end: e.end, rank: e.rank })),
    cdsStart: m.cdsStart, cdsEnd: m.cdsEnd, status: m.isMane ? 'mane' as const : null })) }), [geneLoaders, models]);
  const T_END = scene.tEnd;

  // ---- playback (one rAF, delta-time based) ----
  useEffect(() => {
    if (!playing) { last.current = null; return; }
    const step = (now: number) => {
      if (last.current != null) {
        const next = Math.min(T_END, tRef.current + ((now - last.current) / 1000) * speed);
        setT(next);
        if (next >= T_END) { setPlaying(false); setExplore(true); return; }
      }
      last.current = now;
      raf.current = requestAnimationFrame(step);
    };
    raf.current = requestAnimationFrame(step);
    return () => { if (raf.current != null) cancelAnimationFrame(raf.current); last.current = null; };
  }, [playing, speed, T_END]);
  const togglePlay = useCallback(() => {
    setExplore(false);
    setPlaying(p => { if (!p && tRef.current >= T_END) setT(0); return !p; });
  }, [T_END]);
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onClose();
      if (e.key === ' ' && !/INPUT|SELECT|TEXTAREA/.test((e.target as HTMLElement)?.tagName ?? '')) { e.preventDefault(); togglePlay(); }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [onClose, togglePlay]);
  const seek = (v: number) => { setT(v); last.current = null; setExplore(v >= T_END); };
  const replay = () => { setT(0); last.current = null; setExplore(false); setPlaying(true); };
  const toggleJunction = (id: string) => { setSelected(prev => { const n = new Set(prev); if (n.has(id)) n.delete(id); else n.add(id); return n; }); replay(); };

  // ---- explore: zoom and pan the still figure ----
  const resetView = () => setView({ k: 1, x: 0, y: 0 });
  const toScene = (e: { clientX: number; clientY: number }) => { const r = svgRef.current!.getBoundingClientRect(); return { x: ((e.clientX - r.left) / r.width) * VW, y: ((e.clientY - r.top) / r.height) * scene.H }; };
  const zoomAt = (factor: number, cx = VW / 2, cy = scene.H / 2) => setView(v => { const k = Math.max(0.5, Math.min(8, v.k * factor)); return { k, x: cx - ((cx - v.x) / v.k) * k, y: cy - ((cy - v.y) / v.k) * k }; });
  useEffect(() => {
    const el = svgRef.current;
    if (!el || !explore) return;
    const onWheel = (e: WheelEvent) => { if (!e.ctrlKey && !e.metaKey) return; e.preventDefault(); const p = toScene(e); zoomAt(e.deltaY < 0 ? 1.15 : 1 / 1.15, p.x, p.y); };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, [explore]); // eslint-disable-line react-hooks/exhaustive-deps
  const onPanStart = (e: React.MouseEvent) => { if (!explore || e.button !== 0) return; panRef.current = { x: e.clientX, y: e.clientY, vx: view.x, vy: view.y, moved: false }; };
  const onPanMove = (e: React.MouseEvent) => {
    const p = panRef.current; if (!p) return;
    const r = svgRef.current!.getBoundingClientRect();
    const dx = ((e.clientX - p.x) / r.width) * VW, dy = ((e.clientY - p.y) / r.height) * scene.H;
    if (Math.abs(dx) + Math.abs(dy) > 2) p.moved = true;
    setView(v => ({ ...v, x: p.vx + dx, y: p.vy + dy }));
  };
  const onPanEnd = () => { panRef.current = null; };
  const fileStem = `sv_${chrom}_${chosen.map(j => j.id).join('-') || 'none'}_${sample.name}`.replace(/[^\w.-]+/g, '_');
  const exportSvg = () => {
    const svg = svgRef.current; if (!svg) return;
    const clone = svg.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); clone.setAttribute('width', String(VW)); clone.setAttribute('height', String(scene.H));
    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = `${fileStem}.svg`; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  };
  const exportPng = async () => {
    const svg = svgRef.current; if (!svg) return;
    const clone = svg.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); clone.setAttribute('width', String(VW)); clone.setAttribute('height', String(scene.H));
    const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' }));
    const img = new Image();
    await new Promise<void>((res, rej) => { img.onload = () => res(); img.onerror = () => rej(new Error('render')); img.src = url; });
    const c = document.createElement('canvas'); c.width = VW * 2; c.height = scene.H * 2;
    const g = c.getContext('2d')!; g.scale(2, 2); g.fillStyle = '#fff'; g.fillRect(0, 0, VW, scene.H); g.drawImage(img, 0, 0, VW, scene.H);
    URL.revokeObjectURL(url);
    c.toBlob(b => { if (!b) return; const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = `${fileStem}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, 'image/png');
  };

  const tt = explore ? T_END : t;
  const chapters = [{ label: 'reference', at: 0 }, { label: 'breaks', at: T_BREAK }, { label: 'joins', at: T_J0 }, { label: 'lost', at: scene.tLost }, { label: 'copy number', at: scene.tCn }];
  const verdictCls = depth.verdict === 'explained' ? 'bg-emerald-600 text-white' : depth.verdict === 'balanced' ? 'bg-sky-600 text-white' : depth.verdict === 'unexplained' ? 'bg-red-600 text-white' : 'bg-gray-400 text-white';
  const verdictText = !chosen.length ? 'no junction selected' : depth.verdict === 'explained' ? '✓ depth fits the rebuilt chromosome' : depth.verdict === 'balanced' ? '≈ balanced: depth cannot confirm it' : depth.verdict === 'unexplained' ? `✗ depth does not fit ${depth.segments.filter(s => s.fits === false).map(s => s.seg.letter).join(', ')}` : 'depth not measurable';

  return (
    <div className="fixed inset-0 z-[70] bg-slate-900/80 backdrop-blur-sm flex items-start justify-center overflow-y-auto p-4" onMouseDown={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-[1240px] text-gray-900" onMouseDown={e => e.stopPropagation()} style={{ fontFamily: FONT }}>
        <div className="flex items-start justify-between px-5 pt-4 pb-2">
          <div>
            <div className="text-base font-bold">Structural variant cartoon <span className="font-normal text-gray-500">· {chrom}:{fmt(win.start + 1)}-{fmt(win.end)} · {chosen.length} junction{chosen.length === 1 ? '' : 's'}</span></div>
            <div className="text-xs text-gray-500">observed in {sample.name}{sample.sampled ? ' · sampled window: depths are estimates' : ''} · rebuilt from the junction reads, checked against the read depth · evidence, not a call</div>
          </div>
          <div className="flex items-center gap-2">
            <span className={`px-3 py-1 rounded-full text-xs font-semibold transition-opacity duration-700 ${verdictCls}`} style={{ opacity: tt >= scene.tCn + 0.8 ? 1 : 0.15 }}>{verdictText}</span>
            <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-2xl leading-none px-1" title="Close (Esc)">×</button>
          </div>
        </div>

        <div className="px-3 relative" ref={sceneRef}>
          <svg ref={svgRef} viewBox={`0 0 ${VW} ${scene.H}`} width="100%" className="block select-none" style={{ aspectRatio: `${VW} / ${scene.H}`, cursor: explore ? (panRef.current ? 'grabbing' : 'grab') : 'default' }} fontFamily={FONT}
            onMouseDown={onPanStart} onMouseMove={onPanMove} onMouseUp={onPanEnd} onMouseLeave={onPanEnd} onDoubleClick={resetView}>
            <defs>
              <linearGradient id="svc-bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#f8fafc" /><stop offset="1" stopColor="#ddf1ee" /></linearGradient>
            </defs>
            <rect x={0} y={0} width={VW} height={scene.H} fill="url(#svc-bg)" rx={14} />
            <g transform={explore ? `translate(${view.x}, ${view.y}) scale(${view.k})` : undefined}>
              <SceneDrawing scene={scene} t={tt} depthOf={story.depthOf} chosen={chosen} zygosity={zygosity} copies={copies} />
            </g>
            <Caption t={tt} scene={scene} rec={rec} depth={depth} chosen={chosen} zygosity={zygosity} />
          </svg>
          {explore && (
            <div className="absolute top-3 right-6 flex items-center gap-1 text-xs" onMouseDown={e => e.stopPropagation()}>
              <span className="mr-1 px-2 py-0.5 rounded-full bg-white/90 border border-gray-200 text-gray-600">final picture · Ctrl / ⌘ + scroll to zoom, drag to pan, double-click to reset</span>
              {([['−', () => zoomAt(1 / 1.3)], ['+', () => zoomAt(1.3)], ['fit', resetView]] as [string, () => void][]).map(([lab, fn]) => (
                <button key={lab} onClick={fn} className="w-7 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50 font-semibold">{lab}</button>
              ))}
              <button onClick={exportSvg} className="px-2 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50">SVG</button>
              <button onClick={exportPng} className="px-2 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50">PNG</button>
            </div>
          )}
        </div>

        {/* controls */}
        <div className="px-5 pt-3 pb-2 flex items-center gap-3 text-xs">
          <button onClick={togglePlay} disabled={!chosen.length} className="w-9 h-9 rounded-full bg-accent-600 text-white text-base flex items-center justify-center shadow hover:bg-accent-700 disabled:opacity-40" title="Play / pause (space)">{playing ? '❚❚' : '▶'}</button>
          <input type="range" min={0} max={T_END} step={0.01} value={t} onChange={e => seek(parseFloat(e.target.value))} className="flex-1 accent-accent-600" />
          <span className="font-mono w-16 text-right text-gray-500">{t.toFixed(1)} / {T_END.toFixed(1)} s</span>
          <select value={speed} onChange={e => setSpeed(parseFloat(e.target.value))} className="border rounded px-1 py-0.5 text-xs">
            {[0.5, 1, 1.5, 2].map(v => <option key={v} value={v}>{v}×</option>)}
          </select>
          <button onClick={() => { setPlaying(false); setT(T_END); setExplore(true); }} disabled={!chosen.length}
            className={`px-2.5 py-1 rounded-full border text-xs font-semibold ${explore ? 'bg-accent-600 border-accent-600 text-white' : 'border-accent-300 text-accent-700 hover:bg-accent-50'}`} title="Skip to the final picture and explore it">🔍 Final picture</button>
        </div>
        <div className="px-5 pb-3 flex flex-wrap items-center gap-1.5 text-[11px]">
          {chapters.map(c => (
            <button key={c.label} onClick={() => { seek(c.at + 0.001); setPlaying(true); }}
              className={`px-2 py-0.5 rounded-full border ${tt >= c.at && (chapters.find(o => o.at > c.at)?.at ?? Infinity) > tt ? 'bg-accent-600 border-accent-600 text-white' : 'border-gray-300 text-gray-600 hover:bg-gray-50'}`}>{c.label}</button>
          ))}
          <span className="ml-auto flex items-center gap-3 text-gray-600">
            <label className="flex items-center gap-1" title="One homologue rearranged and the others intact (a constitutional, heterozygous event), or every homologue (homozygous). On a chromosome with one copy the two coincide.">
              zygosity
              <select value={zygosity} onChange={e => { setZygosity(e.target.value as 'het' | 'hom'); }} className="border rounded px-1 py-0.5 text-xs" disabled={copies === 1}>
                <option value="het">heterozygous</option><option value="hom">homozygous</option>
              </select>
            </label>
            <label className="flex items-center gap-1" title="Copies of this chromosome in the sample outside the event: 2 for an autosome, 1 for X or Y in a male (the default on X and Y), 1 for a female's X after X-chromosome loss, etc. The depth of the flanks is taken as this many copies.">
              chromosome copies
              <select value={copies} onChange={e => setCopies(parseInt(e.target.value))} className="border rounded px-1 py-0.5 text-xs">
                {[1, 2, 3, 4].map(n => <option key={n} value={n}>{n}</option>)}
              </select>
            </label>
            {controls.length > 0 && (
              <label className="flex items-center gap-1" title={`Depth as a ratio to the median of the other DNA samples shown (${controls.map(c => c.name).join(', ')}) bin by bin, which cancels GC and mappability waves; off, the sample's own flanks are the reference`}>
                <input type="checkbox" checked={useControls} onChange={e => setUseControls(e.target.checked)} className="accent-accent-600" />normalise on the other samples
              </label>
            )}
          </span>
        </div>

        <Details rec={rec} depth={depth} ranked={ranked} onPath={setPathKey} candidates={candidates} selected={selected} onToggle={toggleJunction} models={models} build={build} bands={bands} chrom={chrom} />
        {chosen.length > 0 && <SvGenomeBrowser rec={rec} chrom={chrom} window={win} segColors={segColorMap} flankCopies={copies} zygosity={zygosity} build={build} loaders={loaders} />}
      </div>
    </div>
  );
}

// ======================== Scene ========================

/** A segment of chromosome: a bar with chevrons for its direction and its letter; `sx` < 0 reads it backwards. */
function SegShape({ x, y, w, sx, seg, opacity = 1, dashed, flankFade }: { x: number; y: number; w: number; sx: number; seg: PlacedSeg; opacity?: number; dashed?: boolean; flankFade?: boolean }) {
  const aw = Math.max(1, w * Math.abs(sx)), cx = x + w / 2, left = cx - aw / 2;
  const ink = luminance(seg.color) > 0.62 ? '#0f1b2b' : '#ffffff';
  const chevrons: JSX.Element[] = [];
  if (!dashed && aw > 18) {
    const dir = sx >= 0 ? 1 : -1, n = Math.max(1, Math.floor(aw / 26));
    for (let i = 0; i < n; i++) {
      const px = left + (i + 0.5) * aw / n;
      if (Math.abs(px - cx) < 9) continue;
      chevrons.push(<path key={i} d={`M${px - 3 * dir},${y - 5} L${px + 2 * dir},${y} L${px - 3 * dir},${y + 5}`} fill="none" stroke={ink} strokeOpacity={0.55} strokeWidth={1.6} strokeLinecap="round" strokeLinejoin="round" />);
    }
  }
  return (
    <g opacity={opacity}>
      <rect x={left} y={y - SEG_H / 2} width={aw} height={SEG_H} rx={3} fill={seg.color} fillOpacity={dashed ? 0.08 : 1} stroke={dashed ? seg.color : 'none'} strokeDasharray={dashed ? '4 3' : undefined} />
      {chevrons}
      {!dashed && aw > 18 && <text x={cx} y={y + 5} fontSize={14} fontWeight={700} textAnchor="middle" fill={ink}>{seg.letter}{sx < 0 ? '′' : ''}</text>}
      {dashed && <text x={cx} y={y + 5} fontSize={13} fontWeight={600} textAnchor="middle" fill={seg.color} fillOpacity={0.6}>{seg.letter}</text>}
      {flankFade && <rect x={seg.first ? left : left + aw - 26} y={y - SEG_H / 2} width={26} height={SEG_H} fill={`url(#svc-fade${seg.first ? 'L' : 'R'})`} />}
    </g>
  );
}

function SceneDrawing({ scene, t, depthOf, chosen, zygosity, copies }: { scene: Scene; t: number; depthOf: Map<number, SegmentDepth>; chosen: SvJunction[]; zygosity: 'het' | 'hom'; copies: number }) {
  const open = phase(t, T_BREAK, T_BREAK + T_OPEN);
  const els: JSX.Element[] = [];
  const allPieces = scene.derRows.flatMap(dr => dr.pieces);
  els.push(<defs key="defs">
    <linearGradient id="svc-fadeL" x1="0" x2="1"><stop offset="0" stopColor="#f8fafc" stopOpacity={1} /><stop offset="1" stopColor="#f8fafc" stopOpacity={0} /></linearGradient>
    <linearGradient id="svc-fadeR" x1="0" x2="1"><stop offset="0" stopColor="#f8fafc" stopOpacity={0} /><stop offset="1" stopColor="#f8fafc" stopOpacity={1} /></linearGradient>
  </defs>);
  for (const rr of scene.refRows) {
    const list = scene.byChr[rr.chrom];
    els.push(<g key={`rl${rr.chrom}`}><text x={16} y={rr.y - 4} fontSize={11} fontWeight={700} fill={C.text}>{rr.chrom}</text><text x={16} y={rr.y + 10} fontSize={10} fill={C.muted}>reference</text></g>);
    for (const seg of list) {
      const x = lerp(seg.xc, seg.xo, open);
      const firstUse = allPieces.find(p => p.seg === seg && p.copy === 0);
      const leaving = firstUse ? phase(t, firstUse.t0, firstUse.t0 + 0.15) : 0;
      const isLost = scene.lost.includes(seg);
      if (firstUse && t >= firstUse.t0) els.push(<SegShape key={`slot${seg.id}`} x={x} y={rr.y} w={seg.w} sx={1} seg={seg} dashed />);
      if (!firstUse || t < firstUse.t0) {
        const fall = isLost ? phase(t, scene.tLost, scene.tLost + 0.9) : 0;
        els.push(<SegShape key={`ref${seg.id}`} x={x} y={rr.y + fall * 6} w={seg.w} sx={1} seg={seg} opacity={1 - 0.6 * fall - leaving} flankFade={seg.first || seg.last} />);
        if (isLost && fall > 0) els.push(<g key={`lost${seg.id}`} opacity={fall}>
          <line x1={x + 4} x2={x + seg.w - 4} y1={rr.y + 6} y2={rr.y + 6} stroke={C.lost} strokeWidth={2} />
          <text x={x + seg.w / 2} y={rr.y + 34} fontSize={10} fontWeight={600} textAnchor="middle" fill={C.lost}>lost</text>
        </g>);
      }
      if (!seg.first && !seg.last) els.push(<text key={`len${seg.id}`} x={x + seg.w / 2} y={rr.y - 21} fontSize={9} textAnchor="middle" fill={C.faint}>{formatBp(seg.end - seg.start)}</text>);
    }
    for (let i = 1; i < list.length; i++) {
      const a = list[i - 1], xb = lerp(a.xc + a.w, (a.xo + a.w + list[i].xo) / 2, open), alpha = clamp01((t - 0.3) / 0.6);
      els.push(<g key={`bp${rr.chrom}${i}`} opacity={alpha}>
        <path d={`M${xb - 4},${rr.y - 24} l5,5 l-5,4 l5,5`} fill="none" stroke={C.lost} strokeWidth={1.6} strokeLinejoin="round" />
        <text x={xb} y={rr.y + (i % 2 ? 30 : 21)} fontSize={9} textAnchor="middle" fill={C.muted} fontFamily={MONO} opacity={open < 1 && i % 2 ? 0.5 : 1}>{fmt(list[i].start)}</text>
      </g>);
    }
    const cnA = phase(t, scene.tCn, scene.tCn + 0.8);
    if (cnA > 0) for (const seg of list) {
      const d = depthOf.get(seg.id);
      if (!d) continue;
      const col = d.fits == null ? C.faint : d.fits ? C.ok : C.lost;
      const l1 = `${d.fits == null ? '?' : d.fits ? '✓' : '✗'} ${d.expected} cop${d.expected === 1 ? 'y' : 'ies'}`;
      const l2 = d.observed == null ? 'too short' : `${d.observed.toFixed(2)} ± ${(2 * (d.se ?? 0)).toFixed(2)}`;
      els.push(<g key={`cn${seg.id}`} opacity={cnA}>
        <text x={seg.xo + seg.w / 2} y={rr.y + 46} fontSize={10} fontWeight={700} textAnchor="middle" fill={col}>{l1}</text>
        <text x={seg.xo + seg.w / 2} y={rr.y + 58} fontSize={9} textAnchor="middle" fill={col} fontFamily={MONO}>{l2}</text>
      </g>);
    }
  }
  els.push(<text key="homo" x={VW - 14} y={scene.refRows[0].y - 40} fontSize={10} textAnchor="end" fill={C.faint}>
    {copies === 1 ? 'one copy of this chromosome: no other homologue' : zygosity === 'het' ? `the other homologue${copies > 2 ? 's keep' : ' keeps'} the reference order` : 'every homologue carries the rearrangement'}
  </text>);
  for (const dr of scene.derRows) {
    const appear = clamp01((t - (dr.pieces[0]?.t0 ?? 0) + 0.3) / 0.4);
    els.push(<g key={`dl${dr.name}`} opacity={appear}>
      <text x={16} y={dr.y - 4} fontSize={11} fontWeight={700} fill={C.accent}>{dr.name}</text><text x={16} y={dr.y + 10} fontSize={10} fill={C.muted}>{dr.rearranged ? 'rebuilt' : 'unchanged'}</text>
      <line x1={LX} x2={RX} y1={dr.y + 22} y2={dr.y + 22} stroke={C.grid} strokeDasharray="2 4" />
    </g>);
    dr.pieces.forEach((p, i) => {
      if (i === 0 || t < p.t0 + MOVE * 0.85) return;
      const prev = dr.pieces[i - 1], xj = p.x - 4, a = clamp01((t - p.t0 - MOVE * 0.85) / 0.3);
      if (p.via === 'ref' || !p.via) { els.push(<line key={`rl${dr.name}${i}`} x1={prev.x + (prev.seg as PlacedSeg).w} x2={p.x} y1={dr.y} y2={dr.y} stroke={C.faint} strokeWidth={2} opacity={a} />); return; }
      const col = KIND_COLOR[p.via.kind], pulse = 1 + 0.6 * Math.max(0, 1 - (t - p.t0 - MOVE) / 0.6);
      els.push(<g key={`j${dr.name}${i}`} opacity={a}>
        <circle cx={xj} cy={dr.y} r={5 * pulse} fill={col} stroke="#fff" strokeWidth={2} />
        <line x1={xj} x2={xj} y1={dr.y - 18} y2={dr.y - 30} stroke={col} />
        <text x={xj} y={dr.y - 34} fontSize={10} fontWeight={700} textAnchor="middle" fill={col}>{p.via.id} · {fmt(p.via.reads)}</text>
      </g>);
    });
  }
  for (const dr of scene.derRows) for (const p of dr.pieces) {
    if (t < p.t0) continue;
    const seg = p.seg as PlacedSeg, sx0 = lerp(seg.xc, seg.xo, open), m = phase(t, p.t0, p.t0 + MOVE);
    const x = lerp(sx0, p.x, m), y = lerp(seg.y, p.y, m) - Math.sin(Math.PI * m) * 46;
    const sx = p.reversed ? Math.cos(Math.PI * m) : 1;
    const fadeIn = p.copy > 0 ? clamp01((t - p.t0) / 0.25) : 1;
    els.push(<SegShape key={`pc${dr.name}${seg.id}${p.copy}`} x={x} y={y} w={seg.w} sx={sx} seg={seg} opacity={fadeIn} />);
    if (p.copy > 0 && m < 1) els.push(<text key={`cp${dr.name}${seg.id}${p.copy}`} x={x + seg.w / 2} y={y - 22} fontSize={10} fontWeight={700} textAnchor="middle" fill={C.dup}>copy</text>);
    if (p.copy > 0 && m >= 1) els.push(<text key={`cp2${dr.name}${seg.id}${p.copy}`} x={p.x + seg.w / 2} y={p.y + 30} fontSize={9.5} fontWeight={700} textAnchor="middle" fill={C.dup}>copy {p.copy + 1}</text>);
    if (p.reversed && m >= 1) els.push(<text key={`inv${dr.name}${seg.id}${p.copy}`} x={p.x + seg.w / 2} y={p.y + 30} fontSize={9.5} fontWeight={700} textAnchor="middle" fill={C.inv}>inverted</text>);
  }
  if (!chosen.length) els.push(<text key="none" x={VW / 2} y={scene.H - 30} fontSize={12} textAnchor="middle" fill={C.faint}>No junction selected: tick the arcs of one event below, or shift-click them on the track.</text>);
  return <>{els}</>;
}

// ======================== Caption ========================

function Caption({ t, scene, rec, depth, chosen, zygosity }: { t: number; scene: Scene; rec: ReturnType<typeof reconstruct>; depth: ReturnType<typeof checkDepth>; chosen: SvJunction[]; zygosity: 'het' | 'hom' }) {
  const letters = (list: SvSegment[]) => list.map(x => x.letter).join(', ');
  let text: string;
  if (!chosen.length) text = '';
  else if (t < T_BREAK) text = `The reference ${rec.chrs.join(' and ')}, cut by the ${chosen.length} selected junction${chosen.length > 1 ? 's' : ''} into ${rec.segs.length} segments (${letters(rec.segs)}).`;
  else if (t < T_J0) text = 'Breaks: each junction joins two breakpoints; the chromosome opens at every one of them.';
  else if (t < scene.tLost) {
    const all = scene.derRows.flatMap(dr => dr.pieces.map(p => ({ p, dr })));
    const cur = [...all].reverse().find(o => t >= o.p.t0) ?? all[0];
    if (!cur) text = '';
    else {
      const { p, dr } = cur, i = dr.pieces.indexOf(p);
      if (i === 0) text = `${dr.name}, from its start: it begins with ${p.seg.letter}.`;
      else {
        const prev = dr.pieces[i - 1];
        const how = p.via === 'ref' || !p.via ? `the reference continues: ${p.seg.letter} follows ${prev.seg.letter}${prev.reversed ? '′' : ''} as usual`
          : `junction ${p.via.id} (${fmt(p.via.reads)} reads: ${p.via.evidence}) joins the end of ${prev.seg.letter}${prev.reversed ? '′' : ''} to ${p.seg.letter}`;
        const what = p.copy > 0 ? `: a second copy of ${p.seg.letter} (tandem duplication)` : p.reversed ? `, read backwards (${p.seg.letter}′, inverted)` : '';
        text = `${dr.name}, step ${i + 1} of ${dr.pieces.length}: ${how}${what}.`;
      }
    }
  } else if (t < scene.tCn) text = rec.lost.length ? `${letters(rec.lost)} ${rec.lost.length > 1 ? 'are' : 'is'} on no rebuilt chromosome: lost from this homologue${zygosity === 'hom' ? ' (and from the other)' : ', the other homologue keeps it'}.` : 'Every segment is kept: the event is balanced on this chromosome.';
  else {
    const bad = depth.segments.filter(s => s.fits === false);
    text = depth.verdict === 'none' ? 'The segments are too short for the depth to say anything (fewer than 3 bins each).'
      : bad.length ? `The read depth does not fit ${letters(bad.map(b => b.seg))}: a junction is missing, one does not belong to this event, or the zygosity or chromosome copies differ${bad.some(b => b.mosaic) ? ' (a mosaic fits better, see the table)' : ''}.`
      : depth.verdict === 'balanced' ? 'Balanced event: no segment changes copy number, so the depth can neither confirm nor refute it; the junction reads are the evidence.'
      : `The read depth fits every segment: the selected junctions explain the coverage (reference: ${depth.reference.how}).`;
  }
  return (
    <g>
      <rect x={16} y={scene.H - 48} width={VW - 32} height={36} rx={8} fill={C.bg} opacity={0.96} stroke={C.grid} />
      <text x={32} y={scene.H - 26} fontSize={13} fill={C.text}>{text}</text>
    </g>
  );
}

// ======================== Details ========================

function Details({ rec, depth, ranked, onPath, candidates, selected, onToggle, models, build, bands, chrom }: {
  rec: ReturnType<typeof reconstruct>; depth: ReturnType<typeof checkDepth>; ranked: RankedPath[]; onPath: (key: string) => void; candidates: SvJunction[]; selected: Set<string>; onToggle: (id: string) => void; models: TxModel[];
  build: string; bands?: Cytoband[]; chrom: string;
}) {
  const side = (e: SvJunction['a']) => (e.side === 'tail' ? '▸' : '◂');
  const breakpoints = [...new Set(rec.junctions.flatMap(j => [j.a.pos, j.b.pos]))].sort((a, b) => a - b);
  return (
    <div className="px-5 pb-4 grid grid-cols-[320px_1fr_1fr] gap-4 text-xs">
      <div>
        <div className="text-[10.5px] font-semibold text-gray-700 mb-1">Junctions in the window ({candidates.length}) · tick the ones of one event</div>
        <div className="space-y-1 max-h-[260px] overflow-y-auto pr-1">
          {candidates.map(j => (
            <label key={j.id} className={`grid grid-cols-[16px_26px_1fr] gap-1.5 items-start p-1.5 rounded-md cursor-pointer ${selected.has(j.id) ? 'bg-accent-50' : 'hover:bg-gray-50'}`}>
              <input type="checkbox" checked={selected.has(j.id)} onChange={() => onToggle(j.id)} className="mt-0.5 accent-accent-600" />
              <span className="font-bold">{j.id}</span>
              <span className="min-w-0">
                <span className="inline-block text-[9px] font-bold tracking-wide uppercase px-1 rounded text-white mr-1" style={{ background: KIND_COLOR[j.kind] }}>{j.kind}</span>
                <span className="font-mono text-[11px] text-gray-700">{fmt(j.a.pos)} {side(j.a)}—{side(j.b)} {fmt(j.b.pos)}</span>
                <div className="text-[11px] text-gray-500">{fmt(j.reads)} reads · {j.evidence} · ±{j.precision} bp{j.assumed ? ` · ${j.assumed}` : ''}</div>
              </span>
            </label>
          ))}
          {!candidates.length && <div className="text-gray-400">no structural arc in this window above the minimum reads</div>}
        </div>
        <div className="text-[10.5px] text-gray-500 mt-1">▸ keeps the sequence left of the breakpoint, ◂ the sequence right of it. An inversion arc gives its two junctions.</div>
      </div>
      <div>
        <div className="text-[10.5px] font-semibold text-gray-700 mb-1">Copy number against the read depth · reference: {depth.reference.note}</div>
        <table className="w-full border-collapse tabular-nums">
          <thead><tr className="text-[10px] text-gray-500 border-b border-gray-200"><th className="text-left font-semibold py-0.5">segment</th><th className="text-left font-semibold">content</th><th className="text-right font-semibold">expected</th><th className="text-right font-semibold">observed ± 2 SE</th><th className="text-right font-semibold">bins</th><th></th></tr></thead>
          <tbody>
            {depth.segments.map(d => {
              const seg = d.seg as PlacedSeg;
              const changes = d.expected !== depth.flankCopies;
              const content = segmentContent(seg, models).genes.map(gn => {
                const parts = [gn.whole.length ? `${exonRun(gn.whole)}${changes && gn.frame === 'shifted' ? ' (frame shifts)' : changes && gn.frame === 'kept' ? ' (frame kept)' : ''}` : '', gn.cut.length ? `cuts ${exonRun(gn.cut)}` : ''].filter(Boolean).join(', ');
                return `${gn.name}: ${parts || 'intronic'}`;
              }).join('; ') || 'no gene model';
              const cls = d.fits == null ? 'text-gray-400' : d.fits ? 'text-emerald-600' : 'text-red-600';
              return (
                <tr key={seg.id} className="border-b border-gray-100 align-top">
                  <td className="py-0.5 whitespace-nowrap"><span className="inline-block w-2.5 h-2.5 rounded-sm mr-1 align-[-1px]" style={{ background: seg.color }} /><b>{seg.letter}</b> <span className="text-gray-400 font-mono text-[10px]">{seg.first ? '…' : fmt(seg.start + 1)}–{seg.last ? '…' : fmt(seg.end)}</span></td>
                  <td className="py-0.5 text-[10.5px] text-gray-600">{content}</td>
                  <td className="py-0.5 text-right font-mono">{d.expected}</td>
                  <td className="py-0.5 text-right font-mono">{d.observed == null ? '—' : `${d.observed.toFixed(2)} ± ${(2 * (d.se ?? 0)).toFixed(2)}`}</td>
                  <td className="py-0.5 text-right font-mono text-gray-500">{d.bins}</td>
                  <td className={`py-0.5 font-bold ${cls}`} title={d.fits == null ? `fewer than ${SV_DEPTH_MIN_BINS} bins inside the margins: too short to assess` : d.mosaic ? `the depth lies between ${depth.flankCopies} and ${d.expected} copies: it fits about ${Math.round(d.mosaic.fraction * 100)} % of cells carrying the event (${Math.round(d.mosaic.lo * 100)}–${Math.round(d.mosaic.hi * 100)} %)` : d.fits ? 'within max(0.25 copies, 3 SE) of the expected copies' : `nearest integer: ${d.nearest} cop${d.nearest === 1 ? 'y' : 'ies'}`}>
                    {d.fits == null ? '?' : d.fits ? '✓' : d.mosaic ? `≈${Math.round(d.mosaic.fraction * 100)} % mosaic` : `✗ ~${d.nearest}`}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
        <div className="text-[10.5px] text-gray-500 mt-1">
          Bins of {formatBp(depth.binBp)} inside each segment, {SV_DEPTH_MARGIN_BP} bp kept away from each breakpoint; a segment fits within max({SV_FIT_MIN_COPIES} copies, {SV_FIT_Z} SE), the SE taking the correlation of neighbouring bins into account; a mosaic is reported between {Math.round(SV_MOSAIC_MIN * 100)} and {Math.round(SV_MOSAIC_MAX * 100)} % of cells; the flanks give the reference when they hold {SV_FLANK_MIN_BINS} bins or more.
        </div>
      </div>
      <div>
        {ranked.length > 1 && (
          <div className="mb-2">
            <div className="text-[10.5px] font-semibold text-gray-700 mb-0.5">{ranked.length} possible paths · the walk has a choice at some breakends</div>
            <div className="space-y-0.5">
              {ranked.slice(0, 5).map((r, i) => {
                const key = r.rec.choice.join(','), on = r.rec === rec;
                return (
                  <button key={key} onClick={() => onPath(key)} className={`w-full text-left rounded px-1.5 py-1 border ${on ? 'border-accent-400 bg-accent-50' : 'border-gray-200 hover:bg-gray-50'}`}>
                    <div className="font-mono text-[11px]">{pathText(r.rec) || 'no rearranged chromosome'}</div>
                    <div className="text-[10.5px] text-gray-600">
                      {i === 0 && <span className="text-emerald-700 font-semibold">best fit · </span>}
                      <span className={r.misfits ? 'text-red-600' : 'text-emerald-700'}>{r.assessed ? (r.misfits ? `depth misfits ${r.misfits} of ${r.assessed} segments (${r.depth.segments.filter(d => d.fits === false).map(d => d.seg.letter).join(', ')})` : `depth fits all ${r.assessed} segments`) : 'depth not assessable'}</span>
                      {r.rec.unplaced.length > 0 && <> · leaves out {r.rec.unplaced.map(j => `${j.id} (${j.reads})`).join(', ')}</>}
                    </div>
                  </button>
                );
              })}
            </div>
            {ranked.length > 5 && <div className="text-[10.5px] text-gray-500">{ranked.length - 5} less likely paths not listed</div>}
            <div className="text-[10.5px] text-gray-500 mt-0.5">Ranked by the segments whose depth does not fit, then the reads left out, then the summed deviation. The junctions left out by the best path may belong to a second version of the rearrangement (another derivative, or a subclone): their share of the reads at the shared breakend is a rough fraction.</div>
          </div>
        )}
        <div className="text-[10.5px] font-semibold text-gray-700 mb-1">Rebuilt chromosome{rec.ders.length > 1 ? 's' : ''}</div>
        {rec.ders.filter(d => d.rearranged || rec.ders.length === 1).map(d => (
          <div key={d.name} className="mb-1.5">
            <span className="font-semibold text-accent-700">{d.name}</span>{d.ends && <span className="text-gray-400"> ({d.ends})</span>}
            <div className="font-mono text-[11px] bg-gray-50 border border-gray-200 rounded px-2 py-1 mt-0.5">{d.pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ') || '—'}</div>
          </div>
        ))}
        {rec.lost.length > 0 && <div className="text-[11px] text-red-700 mb-1">Lost from the rebuilt homologue: {rec.lost.map(x => `${x.letter} (${formatBp(x.end - x.start)})`).join(', ')}.</div>}
        <SvNomenclatureBox rec={rec} build={build} bands={bands} />
        {breakpoints.length > 0 && (
          <div className="mt-1.5">
            <div className="text-[10.5px] font-semibold text-gray-700">Breakpoints ({chrom})</div>
            <ul className="text-[11px] text-gray-600">
              {breakpoints.map(p => <li key={p}><span className="font-mono">{fmt(p + 1)}</span> · {breakpointSite(p, models)}</li>)}
            </ul>
          </div>
        )}
        {rec.merged.length > 0 && <div className="text-[10.5px] text-gray-600 mt-1">Breakends within {SV_BREAKEND_CLUSTER_BP} bp taken as one breakpoint: {rec.merged.join('; ')}.</div>}
        {rec.warnings.map((w, i) => <div key={i} className="text-[11px] text-amber-800 mt-1">⚠ {w}</div>)}
        <div className="text-[10.5px] text-gray-500 mt-1.5">Breakpoints come from the reads, rounded to 5 bp (pairs alone: within an insert size); check them in the alignments before anything goes into a report. The scene draws segment widths on a log scale.</div>
      </div>
    </div>
  );
}
