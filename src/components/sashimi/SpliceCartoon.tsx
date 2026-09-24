/**
 * Splicing cartoon: an animated SVG scene showing how the canonical and the observed transcript
 * are spliced, translated, and either degraded by nonsense-mediated decay or turned into a
 * protein whose domains are kept, disrupted or lost. Pure SVG driven by requestAnimationFrame
 * (no animation library), so it runs at display refresh in the application and in the
 * single-file standalone viewer alike.
 */
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { cdnaPosition, type TxModel } from './geometry';
import type { Block, SpliceStory, Isoform, DomainStatus } from './spliceModel';

// ======================== Timeline ========================

const T_PRE = 1.0;      // pre-mRNA appears
const T_LOOP = 3.2;     // introns loop out
const T_CUT = 4.0;      // lariats released
const T_JOIN = 5.2;     // exons slide together
const T_MRNA = 6.2;     // cap, poly(A), EJCs
const T_TRANS = 11.2;   // translation
const T_END = 14.5;     // verdict
export const CARTOON_DURATION = T_END;

const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
const easeInOut = (p: number) => (p < 0.5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2);
const easeOut = (p: number) => 1 - Math.pow(1 - p, 3);
/** progress of [a, b] at time t, eased */
const phase = (t: number, a: number, b: number, ease = easeInOut) => ease(clamp01((t - a) / (b - a)));
const lerp = (a: number, b: number, p: number) => a + (b - a) * p;

// ======================== Palette ========================

const FONT = 'Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const C = {
  exon: '#334155', utr: '#94a3b8', intron: '#94a3b8', text: '#1f2937', muted: '#6b7280', faint: '#9ca3af',
  cryptic: '#b45309', lost: '#dc2626', ejc: '#7c3aed', ribosomeA: '#64748b', ribosomeB: '#94a3b8', upf1: '#dc2626',
  chain: '#cbd5e1', neo: '#f59e0b', ok: '#16a34a', bg: '#ffffff',
};
const DOMAIN_COLORS = ['#2a78d6', '#1baf7a', '#eda100', '#e87ba4', '#4a3aa7', '#008300', '#eb6834'];

// ======================== Layout ========================

const VW = 1200, VH = 640;
const X0 = 150, X1 = 1180;                 // drawable band
const ROW = [{ mrna: 128, protein: 218 }, { mrna: 378, protein: 468 }];   // canonical, observed
const EXON_H = 26, UTR_H = 14;
const GAP_PRE = 44;                        // intron width in the pre-mRNA
const CHAIN_Y_OFF = 0;

interface Placed extends Block {
  w: number;            // width (px)
  xPre: number;         // left x in the pre-mRNA
  xPost: number;        // left x once spliced
  len: number;
  /** genomic sub-ranges [start, end) that are coding */
  cds: [number, number] | null;
}

/** Widths grow with the square root of the length so short exons stay visible and long ones do not dominate. */
function widthOf(len: number, k: number) { return Math.max(14, Math.min(130, (10 + 3.4 * Math.sqrt(len)) * k)); }

function place(iso: Isoform, tx: TxModel, k: number, gapK: number): Placed[] {
  let xPre = X0, xPost = X0;
  return iso.blocks.map((b, i) => {
    const len = b.end - b.start, w = widthOf(len, k);
    if (i > 0 && !b.joinsPrev) xPre += GAP_PRE * gapK;
    const cds = tx.cdsStart != null && tx.cdsEnd != null ? [Math.max(b.start, tx.cdsStart), Math.min(b.end, tx.cdsEnd)] as [number, number] : null;
    const p: Placed = { ...b, w, xPre, xPost, len, cds: cds && cds[1] > cds[0] ? cds : null };
    xPre += w; xPost += w;
    return p;
  });
}

/** Scale factor so that the widest pre-mRNA row fits the band. */
function fitScale(isos: Isoform[]): { k: number; gapK: number } {
  const need = Math.max(...isos.map(iso => {
    let x = 0;
    iso.blocks.forEach((b, i) => { if (i > 0 && !b.joinsPrev) x += GAP_PRE; x += widthOf(b.end - b.start, 1); });
    return x;
  }));
  const avail = X1 - X0;
  if (need <= avail) return { k: 1, gapK: 1 };
  const gapK = Math.max(0.45, Math.min(1, avail / need));
  // shrink widths only as far as needed after the gaps were squeezed
  const needGapped = Math.max(...isos.map(iso => {
    let x = 0;
    iso.blocks.forEach((b, i) => { if (i > 0 && !b.joinsPrev) x += GAP_PRE * gapK; x += widthOf(b.end - b.start, 1); });
    return x;
  }));
  return { k: Math.min(1, avail / needGapped), gapK };
}

// ======================== Component ========================

export interface SpliceCartoonProps {
  story: SpliceStory | null;
  loading: string | null;
  error?: string;
  tx: TxModel;
  sampleName: string;
  sampleColor: string;
  junctionLabel: string;
  onClose: () => void;
}

export default function SpliceCartoon({ story, loading, error, tx, sampleName, sampleColor, junctionLabel, onClose }: SpliceCartoonProps) {
  const [t, setT] = useState(0);
  const [playing, setPlaying] = useState(true);
  const [speed, setSpeed] = useState(1);
  /** After the last frame the scene becomes a still figure to explore: zoom, pan, hover, select, export. */
  const [explore, setExplore] = useState(false);
  const [view, setView] = useState({ k: 1, x: 0, y: 0 });
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null);
  const [selected, setSelected] = useState<number | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const sceneRef = useRef<HTMLDivElement>(null);
  const panRef = useRef<{ x: number; y: number; vx: number; vy: number; moved: boolean } | null>(null);
  const raf = useRef<number | null>(null);
  const last = useRef<number | null>(null);
  const tRef = useRef(0);
  tRef.current = t;

  // ---- playback loop (one rAF, delta-time based, so it is fluid whatever the refresh rate) ----
  useEffect(() => {
    if (!playing || !story) { last.current = null; return; }
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
  }, [playing, speed, story]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onClose();
      if (e.key === ' ') { e.preventDefault(); togglePlay(); }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [onClose]); // eslint-disable-line react-hooks/exhaustive-deps

  const togglePlay = useCallback(() => {
    setExplore(false);
    setPlaying(p => { if (!p && tRef.current >= T_END) setT(0); return !p; });
  }, []);
  const seek = (v: number) => { setT(v); last.current = null; setExplore(v >= T_END); };
  const resetView = () => setView({ k: 1, x: 0, y: 0 });
  /** pointer position in viewBox units */
  const toScene = (e: { clientX: number; clientY: number }) => {
    const r = svgRef.current!.getBoundingClientRect();
    return { x: ((e.clientX - r.left) / r.width) * VW, y: ((e.clientY - r.top) / r.height) * VH };
  };
  const zoomAt = (factor: number, cx = VW / 2, cy = VH / 2) => setView(v => {
    const k = Math.max(0.5, Math.min(8, v.k * factor));
    // keep the scene point under (cx, cy) fixed
    return { k, x: cx - ((cx - v.x) / v.k) * k, y: cy - ((cy - v.y) / v.k) * k };
  });
  // non-passive wheel listener so the page does not scroll while zooming the figure
  useEffect(() => {
    const el = svgRef.current;
    if (!el || !explore) return;
    const onWheel = (e: WheelEvent) => {
      e.preventDefault();
      const p = toScene(e);
      zoomAt(e.deltaY < 0 ? 1.15 : 1 / 1.15, p.x, p.y);
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, [explore]); // eslint-disable-line react-hooks/exhaustive-deps
  const onPanStart = (e: React.MouseEvent) => {
    if (!explore || e.button !== 0) return;
    panRef.current = { x: e.clientX, y: e.clientY, vx: view.x, vy: view.y, moved: false };
  };
  const onPanMove = (e: React.MouseEvent) => {
    const p = panRef.current;
    if (!p) return;
    const r = svgRef.current!.getBoundingClientRect();
    const dx = ((e.clientX - p.x) / r.width) * VW, dy = ((e.clientY - p.y) / r.height) * VH;
    if (Math.abs(dx) + Math.abs(dy) > 2) p.moved = true;
    setView(v => ({ ...v, x: p.vx + dx, y: p.vy + dy }));
  };
  const onPanEnd = () => { const p = panRef.current; panRef.current = null; if (p && !p.moved) setSelected(null); };
  const hover = useCallback((lines: string[] | null, e?: React.MouseEvent) => {
    if (!lines || !e || !sceneRef.current) { setTip(null); return; }
    const r = sceneRef.current.getBoundingClientRect();
    setTip({ x: e.clientX - r.left, y: e.clientY - r.top, lines });
  }, []);
  const exportSvg = () => {
    const svg = svgRef.current; if (!svg) return;
    const clone = svg.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('width', String(VW)); clone.setAttribute('height', String(VH));
    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = `splicing_${tx.geneName}_${sampleName}.svg`; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  };
  const exportPng = async () => {
    const svg = svgRef.current; if (!svg) return;
    const clone = svg.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); clone.setAttribute('width', String(VW)); clone.setAttribute('height', String(VH));
    const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' }));
    const img = new Image();
    await new Promise<void>((res, rej) => { img.onload = () => res(); img.onerror = () => rej(new Error('render')); img.src = url; });
    const c = document.createElement('canvas'); c.width = VW * 2; c.height = VH * 2;
    const g = c.getContext('2d')!; g.scale(2, 2); g.fillStyle = '#fff'; g.fillRect(0, 0, VW, VH); g.drawImage(img, 0, 0, VW, VH);
    URL.revokeObjectURL(url);
    c.toBlob(b => { if (!b) return; const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = `splicing_${tx.geneName}_${sampleName}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, 'image/png');
  };
  const chapters = [
    { label: 'pre-mRNA', at: 0 }, { label: 'splicing', at: T_PRE }, { label: 'mRNA', at: T_JOIN }, { label: 'translation', at: T_MRNA }, { label: 'outcome', at: T_TRANS },
  ];

  const scene = useMemo(() => {
    if (!story) return null;
    const { k, gapK } = fitScale([story.canonical, story.aberrant]);
    const canon = place(story.canonical, tx, k, gapK), ab = place(story.aberrant, tx, k, gapK);
    // the spliced mRNA stretches to fill the band (same factor for both rows so exon sizes stay comparable)
    const widest = Math.max(canon.reduce((a, p) => a + p.w, 0), ab.reduce((a, p) => a + p.w, 0));
    const expand = Math.max(1, Math.min(1.9, (X1 - X0 - 90) / Math.max(1, widest)));
    return { canon, ab, expand };
  }, [story, tx]);

  return (
    <div className="fixed inset-0 z-[70] bg-slate-900/80 backdrop-blur-sm flex items-start justify-center overflow-y-auto p-4" onMouseDown={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-[1240px] text-gray-900" onMouseDown={e => e.stopPropagation()} style={{ fontFamily: FONT }}>
        {/* header */}
        <div className="flex items-start justify-between px-5 pt-4 pb-2">
          <div>
            <div className="text-base font-bold">Splicing cartoon <span className="font-normal text-gray-500">· {tx.geneName} · {junctionLabel}</span></div>
            <div className="text-xs text-gray-500">{story ? story.event.text : loading || error || ''} · observed in {sampleName}</div>
          </div>
          <div className="flex items-center gap-2">
            {story && <VerdictBadge story={story} reveal={t >= T_TRANS + 1.2} />}
            <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-2xl leading-none px-1" title="Close (Esc)">×</button>
          </div>
        </div>

        {/* scene */}
        <div className="px-3 relative" ref={sceneRef}>
          {story && scene ? (
            <svg ref={svgRef} viewBox={`0 0 ${VW} ${VH}`} width="100%" className="block select-none" style={{ aspectRatio: `${VW} / ${VH}`, cursor: explore ? (panRef.current ? 'grabbing' : 'grab') : 'default' }} fontFamily={FONT}
              onMouseDown={onPanStart} onMouseMove={onPanMove} onMouseUp={onPanEnd} onMouseLeave={() => { onPanEnd(); setTip(null); }} onDoubleClick={resetView}>
              <defs>
                <linearGradient id="sc-bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#f8fafc" /><stop offset="1" stopColor="#eef2ff" /></linearGradient>
                <radialGradient id="sc-ribo" cx="35%" cy="30%" r="80%"><stop offset="0" stopColor="#cbd5e1" /><stop offset="1" stopColor="#64748b" /></radialGradient>
                <radialGradient id="sc-upf1" cx="35%" cy="30%" r="80%"><stop offset="0" stopColor="#fca5a5" /><stop offset="1" stopColor="#b91c1c" /></radialGradient>
                <pattern id="sc-hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#fff7ed" /><line x1="0" y1="0" x2="0" y2="6" stroke={C.cryptic} strokeWidth="2" /></pattern>
                <pattern id="sc-hatch-lost" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#fef2f2" /><line x1="0" y1="0" x2="0" y2="6" stroke={C.lost} strokeWidth="1.5" /></pattern>
                <filter id="sc-shadow" x="-20%" y="-50%" width="140%" height="220%"><feDropShadow dx="0" dy="1.5" stdDeviation="1.5" floodColor="#0f172a" floodOpacity="0.18" /></filter>
              </defs>
              <rect x={0} y={0} width={VW} height={VH} fill="url(#sc-bg)" rx={14} />
              <g transform={explore ? `translate(${view.x}, ${view.y}) scale(${view.k})` : undefined}>
                <Row iso={story.canonical} placed={scene.canon} expand={scene.expand} row={0} t={explore ? T_END : t} tx={tx} story={story} kind="canonical" color={C.exon} label={['MANE Select', tx.transcriptId]}
                  explore={explore} onHover={hover} selected={selected} onSelect={setSelected} />
                <Row iso={story.aberrant} placed={scene.ab} expand={scene.expand} row={1} t={explore ? T_END : t} tx={tx} story={story} kind="observed" color={sampleColor} label={['Observed', sampleName]}
                  explore={explore} onHover={hover} selected={selected} onSelect={setSelected} />
              </g>
              {explore ? <Legend story={story} color={sampleColor} /> : <Caption t={t} story={story} />}
            </svg>
          ) : (
            <div className="h-[420px] flex items-center justify-center text-sm text-gray-500">{error ? <span className="text-red-600">{error}</span> : loading || 'preparing the scene…'}</div>
          )}
          {explore && story && (
            <div className="absolute top-3 right-6 flex items-center gap-1 text-xs" onMouseDown={e => e.stopPropagation()}>
              <span className="mr-1 px-2 py-0.5 rounded-full bg-white/90 border border-gray-200 text-gray-600">final picture · scroll to zoom, drag to pan, hover for details, click to highlight, double-click to reset</span>
              {[['−', () => zoomAt(1 / 1.3)], ['+', () => zoomAt(1.3)], ['fit', resetView]].map(([lab, fn]) => (
                <button key={lab as string} onClick={fn as () => void} className="w-7 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50 font-semibold">{lab as string}</button>
              ))}
              <button onClick={exportSvg} className="px-2 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50">SVG</button>
              <button onClick={exportPng} className="px-2 h-7 rounded bg-white border border-gray-300 text-gray-700 hover:bg-gray-50">PNG</button>
            </div>
          )}
          {tip && explore && (
            <div className="pointer-events-none absolute z-10 rounded-md bg-gray-900/90 text-white text-[11px] px-2.5 py-1.5 shadow-lg max-w-[320px]" style={{ left: Math.min(tip.x + 14, (sceneRef.current?.clientWidth ?? 800) - 330), top: tip.y + 12 }}>
              {tip.lines.map((l, i) => <div key={i} className={i === 0 ? 'font-semibold' : 'text-gray-200'}>{l}</div>)}
            </div>
          )}
        </div>

        {/* controls */}
        <div className="px-5 pt-3 pb-2 flex items-center gap-3 text-xs">
          <button onClick={togglePlay} className="w-9 h-9 rounded-full bg-indigo-600 text-white text-base flex items-center justify-center shadow hover:bg-indigo-700" title="Play / pause (space)">{playing ? '❚❚' : '▶'}</button>
          <input type="range" min={0} max={T_END} step={0.01} value={t} onChange={e => seek(parseFloat(e.target.value))} className="flex-1 accent-indigo-600" />
          <span className="font-mono w-16 text-right text-gray-500">{t.toFixed(1)} / {T_END.toFixed(1)} s</span>
          <select value={speed} onChange={e => setSpeed(parseFloat(e.target.value))} className="border rounded px-1 py-0.5 text-xs">
            {[0.5, 1, 1.5, 2].map(v => <option key={v} value={v}>{v}×</option>)}
          </select>
          <button onClick={() => { setPlaying(false); setT(T_END); setExplore(true); }} disabled={!story}
            className={`px-2.5 py-1 rounded-full border text-xs font-semibold ${explore ? 'bg-indigo-600 border-indigo-600 text-white' : 'border-indigo-300 text-indigo-700 hover:bg-indigo-50'}`} title="Skip to the final picture and explore it">🔍 Final picture</button>
        </div>
        <div className="px-5 pb-3 flex flex-wrap gap-1.5 text-[11px]">
          {chapters.map(c => (
            <button key={c.label} onClick={() => { seek(c.at); setPlaying(true); }}
              className={`px-2 py-0.5 rounded-full border ${t >= c.at && (chapters.find(o => o.at > c.at)?.at ?? Infinity) > t ? 'bg-indigo-600 border-indigo-600 text-white' : 'border-gray-300 text-gray-600 hover:bg-gray-50'}`}>{c.label}</button>
          ))}
        </div>

        {/* amino-acid strip + domains */}
        {story && <Details story={story} />}
      </div>
    </div>
  );
}

// ======================== Verdict badge ========================

function VerdictBadge({ story, reveal }: { story: SpliceStory; reveal: boolean }) {
  const v = story.verdict;
  const cls = v.degraded ? 'bg-red-600 text-white' : v.verdict === 'no_ptc' ? (story.diff.kind === 'identical' ? 'bg-emerald-600 text-white' : 'bg-amber-500 text-white') : v.verdict === 'escape_start_proximal' ? 'bg-amber-500 text-white' : 'bg-orange-600 text-white';
  return (
    <span className={`px-3 py-1 rounded-full text-xs font-semibold transition-opacity duration-700 ${cls}`} style={{ opacity: reveal ? 1 : 0.15 }} title={v.text}>
      {v.degraded ? '⛔ degraded by NMD' : v.verdict === 'no_ptc' ? (story.diff.kind === 'identical' ? '✓ normal protein' : `✎ ${story.diff.text}`) : `✂ ${v.headline}`}
    </span>
  );
}

// ======================== One transcript row ========================

interface RowProps {
  iso: Isoform; placed: Placed[]; expand: number; row: 0 | 1; t: number; tx: TxModel; story: SpliceStory; kind: 'canonical' | 'observed'; color: string; label: [string, string];
  /** still-picture mode: no dissolve, hover tooltips, click to highlight an exon in both rows */
  explore?: boolean;
  onHover?: (lines: string[] | null, e?: React.MouseEvent) => void;
  selected?: number | null;
  onSelect?: (rank: number | null) => void;
}

function Row({ iso, placed: placed0, expand, row, t, tx, story, kind, color, label, explore, onHover, selected, onSelect }: RowProps) {
  const y = ROW[row].mrna, yP = ROW[row].protein;
  const plus = tx.strand > 0;
  const appear = phase(t, 0.1 + row * 0.25, T_PRE, easeOut);
  const loop = phase(t, T_PRE, T_LOOP);
  const cut = phase(t, T_LOOP, T_CUT, easeOut);
  const join = phase(t, T_CUT, T_JOIN);
  const mrna = phase(t, T_JOIN, T_MRNA, easeOut);
  const trans = phase(t, T_MRNA, T_TRANS);
  const outcome = phase(t, T_TRANS, T_END, easeOut);
  const isObserved = kind === 'observed';
  const degraded = isObserved && story.verdict.degraded;
  const dissolve = degraded && !explore ? phase(t, T_TRANS + 1.3, T_TRANS + 2.6, easeOut) : 0;
  const hov = (lines: string[]) => (explore && onHover ? { onMouseEnter: (e: React.MouseEvent) => onHover(lines, e), onMouseMove: (e: React.MouseEvent) => onHover(lines, e), onMouseLeave: () => onHover(null) } : {});
  const cpos = (g: number) => cdnaPosition(g, tx).label;

  // blocks widen as they join (the mature mRNA fills the band); x of every block at time t (pre → post)
  const wk = lerp(1, expand, join);
  const placed: Placed[] = [];
  let xPostAcc = X0;
  for (const p of placed0) { placed.push({ ...p, w: p.w * wk, xPost: xPostAcc }); xPostAcc += p.w * wk; }
  const xs = placed.map((p, i) => lerp(placed0[i].xPre, p.xPost, join));
  const right = placed.length ? xs[placed.length - 1] + placed[placed.length - 1].w : X0;

  // ---- mRNA index → x (after joining) ----
  const xAt = (m: number): number => {
    for (let i = 0; i < placed.length; i++) {
      const s = iso.blockStart[i], len = placed[i].len;
      if (m >= s && m < s + len) return xs[i] + ((m - s) / len) * placed[i].w;
    }
    return m < 0 ? xs[0] : right;
  };

  // ---- translation state ----
  const cds = iso.cdsStart, stop = iso.stopIndex ?? iso.mrna.length;
  // both ribosomes advance at the same speed (nt per second); the observed one stops early at a PTC
  const spanOf = (x: Isoform) => (x.cdsStart != null ? (x.stopIndex ?? x.mrna.length) - x.cdsStart : 0);
  const span = Math.max(spanOf(story.canonical), spanOf(story.aberrant), 1);
  const readPos = cds != null ? Math.min(stop, cds + span * trans) : null;   // mRNA index under the ribosome
  const tArrive = T_MRNA + (T_TRANS - T_MRNA) * ((stop - (cds ?? 0)) / span);   // when this ribosome reaches its stop
  const aaMade = cds != null && readPos != null ? Math.max(0, Math.floor((readPos - cds) / 3)) : 0;
  const protLen = iso.protein.length;
  const aaPx = Math.min(2.4, (X1 - X0 - 80) / Math.max(1, story.canonical.protein.length));
  const chainX0 = X0;

  // ---- the skipped / affected exons in the pre-mRNA of the observed row (drawn as ghosts that leave with the lariat) ----
  const removed = isObserved ? story.event.canonical.filter(cb => !story.event.blocks.some(b => b.rank === cb.rank && b.kind === 'exon')) : [];
  const affectedRank = new Set(story.event.affected);

  // ---- introns between consecutive non-contiguous blocks: lariat loops ----
  const introns: JSX.Element[] = [];
  for (let i = 1; i < placed.length; i++) {
    if (placed[i].joinsPrev) continue;
    const xa = xs[i - 1] + placed[i - 1].w, xb = xs[i];
    const gapNow = xb - xa;
    const gapPre = placed[i].xPre - (placed[i - 1].xPre + placed[i - 1].w);
    const h = (34 + gapPre * 0.9) * loop;
    const mid = (xa + xb) / 2;
    // loop: teardrop from xa to xb rising by h; after the cut it detaches and floats away
    const d = `M${xa},${y} C${xa + gapNow * 0.05},${y - h} ${xb - gapNow * 0.05},${y - h} ${xb},${y}`;
    const removedHere = removed.filter(r => r.start >= placed[i - 1].end && r.end <= placed[i].start || (!plus && r.start >= placed[i].end && r.end <= placed[i - 1].start));
    introns.push(
      <g key={`in${i}`} opacity={(1 - cut) * appear} transform={`translate(0, ${-40 * cut})`}>
        <path d={d} fill="none" stroke={C.intron} strokeWidth={2} strokeLinecap="round" />
        {loop > 0.15 && <circle cx={mid} cy={y - h * 0.75} r={3 + 2 * loop} fill="#fff" stroke={C.intron} strokeWidth={1.5} opacity={loop} />}
        {removedHere.map(r => {
          const w = widthOf(r.end - r.start, 1) * 0.8;
          const yy = y - h * 0.62;
          return (
            <g key={`rm${r.rank}`} transform={`translate(${mid - w / 2}, ${yy - EXON_H / 2})`}>
              <rect width={w} height={EXON_H} rx={5} fill={C.lost} opacity={0.55 + 0.35 * (1 - loop)} />
              <text x={w / 2} y={EXON_H / 2 + 4} textAnchor="middle" fill="#fff" fontSize={11} fontWeight={700}>{r.rank}</text>
              {loop > 0.5 && <line x1={w * 0.15} y1={EXON_H * 0.15} x2={w * 0.85} y2={EXON_H * 0.85} stroke="#fff" strokeWidth={2} opacity={(loop - 0.5) * 2} />}
            </g>
          );
        })}
      </g>,
    );
  }

  // ---- blocks ----
  const blocks = placed.map((p, i) => {
    const x = xs[i];
    const isCryptic = p.kind === 'cryptic', isExt = p.kind === 'extension';
    const fill = isCryptic || isExt ? 'url(#sc-hatch)' : color;
    const affected = isObserved && p.rank != null && affectedRank.has(p.rank) && !isCryptic && !isExt;
    // portions: coding at full height, UTR at half height
    const parts: JSX.Element[] = [];
    if (p.cds && !isCryptic && !isExt) {
      const gs = plus ? p.start : p.end, ge = plus ? p.end : p.start;   // transcription order
      const f = (g: number) => (Math.abs(g - gs) / p.len) * p.w;
      const c0 = plus ? f(p.cds[0]) : f(p.cds[1]), c1 = plus ? f(p.cds[1]) : f(p.cds[0]);
      const lo = Math.min(c0, c1), hi = Math.max(c0, c1);
      if (lo > 0.5) parts.push(<rect key="u5" x={x} y={y - UTR_H / 2} width={lo} height={UTR_H} rx={3} fill={C.utr} />);
      parts.push(<rect key="c" x={x + lo} y={y - EXON_H / 2} width={Math.max(1, hi - lo)} height={EXON_H} rx={4} fill={fill} />);
      if (p.w - hi > 0.5) parts.push(<rect key="u3" x={x + hi} y={y - UTR_H / 2} width={p.w - hi} height={UTR_H} rx={3} fill={C.utr} />);
      void ge;
    } else if (!p.cds && !isCryptic && !isExt && tx.cdsStart != null) {
      parts.push(<rect key="u" x={x} y={y - UTR_H / 2} width={p.w} height={UTR_H} rx={3} fill={C.utr} />);
    } else {
      parts.push(<rect key="x" x={x} y={y - EXON_H / 2} width={p.w} height={EXON_H} rx={4} fill={fill} stroke={isCryptic || isExt ? C.cryptic : 'none'} strokeWidth={1.2} />);
    }
    // cryptic pieces are invisible in the pre-mRNA (they are intron) and materialise as the loop forms
    const vis = isCryptic || isExt ? Math.max(0.12, loop) : 1;
    const labelTxt = isCryptic ? 'cryptic' : isExt ? `+${p.len}` : p.rank != null ? String(p.rank) : '';
    const first = plus ? p.start : p.end - 1, lastB = plus ? p.end - 1 : p.start;
    const lines = isCryptic
      ? [`Cryptic exon · ${p.len} nt`, `${tx.chrom}:${(p.start + 1).toLocaleString('en-US')}-${p.end.toLocaleString('en-US')}`, `${cpos(first)} → ${cpos(lastB)} (intronic offsets)`]
      : isExt ? [`Intron retained into exon ${p.rank} · ${p.len} nt`, `${tx.chrom}:${(p.start + 1).toLocaleString('en-US')}-${p.end.toLocaleString('en-US')}`, `${cpos(first)} → ${cpos(lastB)}`]
      : [`Exon ${p.rank}${affected ? ' · affected' : ''} · ${p.len} nt`, `${tx.chrom}:${(p.start + 1).toLocaleString('en-US')}-${p.end.toLocaleString('en-US')}`, `${cpos(first)} → ${cpos(lastB)}`, p.cds ? `coding: ${p.cds[1] - p.cds[0]} nt` : 'non-coding (UTR)'];
    const isSel = explore && selected != null && p.rank === selected && !isCryptic && !isExt;
    return (
      <g key={`b${i}`} opacity={appear * vis * (1 - dissolve * 0.85)} filter="url(#sc-shadow)" {...hov(lines)}
        style={explore ? { cursor: 'pointer' } : undefined}
        onClick={explore && onSelect && p.rank != null && !isCryptic && !isExt ? (e => { e.stopPropagation(); onSelect(selected === p.rank ? null : p.rank!); }) : undefined}
        transform={dissolve ? `translate(${(i % 2 ? 1 : -1) * 18 * dissolve}, ${(i % 3 - 1) * 22 * dissolve}) rotate(${(i % 2 ? 8 : -8) * dissolve}, ${x + p.w / 2}, ${y})` : undefined}>
        {isSel && <rect x={x - 4} y={y - EXON_H / 2 - 4} width={p.w + 8} height={EXON_H + 8} rx={7} fill="none" stroke="#6366f1" strokeWidth={2.5} opacity={0.9} />}
        {parts}
        {affected && <rect x={x - 1.5} y={y - EXON_H / 2 - 1.5} width={p.w + 3} height={EXON_H + 3} rx={5} fill="none" stroke={C.lost} strokeWidth={1.5} strokeDasharray="3 2" />}
        {p.w >= 16 && <text x={x + p.w / 2} y={y + 4} textAnchor="middle" fill={isCryptic || isExt ? C.cryptic : '#fff'} fontSize={isCryptic || isExt ? 9 : 11} fontWeight={700}>{labelTxt}</text>}
      </g>
    );
  });

  // ---- mRNA decorations: cap, poly(A), EJCs ----
  const deco: JSX.Element[] = [];
  if (mrna > 0) {
    deco.push(<g key="cap" opacity={mrna}><circle cx={X0 - 12} cy={y} r={6} fill="#fbbf24" stroke="#b45309" strokeWidth={1} /><text x={X0 - 12} y={y - 11} textAnchor="middle" fontSize={8.5} fill={C.muted}>cap</text></g>);
    deco.push(<g key="pa" opacity={mrna * (1 - dissolve)}><text x={right + 6} y={y + 4} fontSize={10} fill={C.muted} fontFamily="ui-monospace, monospace">AAAAAAA</text></g>);
    iso.junctions.forEach((jm, i) => {
      const ex = xAt(jm) - 6;
      const passed = readPos != null && readPos > jm + 4;
      const downstream = iso.stopIndex != null && jm > iso.stopIndex;
      const pulse = degraded && downstream && readPos != null && readPos >= stop ? 1 + 0.25 * Math.sin(t * 9 + i) : 1;
      deco.push(
        <g key={`ejc${i}`} opacity={mrna * (passed ? 0 : 1) * (1 - dissolve)} transform={`translate(${ex}, ${y - EXON_H / 2 - 9}) scale(${pulse})`}
          {...hov([`Exon junction complex ${i + 1}`, `deposited ~20–24 nt upstream of the junction at mRNA position ${jm.toLocaleString('en-US')}`, passed ? 'displaced by the ribosome' : downstream ? 'still bound downstream of the stop: NMD signal' : 'not yet reached'])}>
          <circle r={5.5} fill={C.ejc} stroke="#fff" strokeWidth={1.2} />
          {i === 0 && <text y={-9} textAnchor="middle" fontSize={8.5} fill={C.ejc} fontWeight={600}>EJC</text>}
        </g>,
      );
    });
  }

  // ---- ribosome + nascent chain ----
  const ribo: JSX.Element[] = [];
  if (readPos != null && mrna >= 1 && trans > 0) {
    const rx = xAt(readPos);
    const stalled = degraded && readPos >= stop;
    const leave = !degraded ? outcome : 0;
    const chainTipX = chainX0 + aaMade * aaPx;
    ribo.push(
      <g key="ribo" opacity={1 - Math.max(leave, dissolve) } transform={`translate(${rx}, ${y})`}
        {...hov(stalled ? ['Ribosome stalled at the premature termination codon', `codon ${story.verdict.ptcCodon} · ${story.verdict.distanceFromStart} nt after the AUG`, story.verdict.headline] : ['Ribosome', `reading mRNA position ${Math.round(readPos).toLocaleString('en-US')} · ${aaMade} aa made`])}>
        <ellipse cx={0} cy={-15} rx={30} ry={17} fill="url(#sc-ribo)" opacity={0.85} />
        <ellipse cx={0} cy={13} rx={24} ry={10} fill="url(#sc-ribo)" opacity={0.7} />
        {stalled && <text y={-38} textAnchor="middle" fontSize={9.5} fill={C.lost} fontWeight={700}>PTC</text>}
        <path d={`M-6,-30 c -4,-10 4,-16 0,-26 c -3,-8 5,-12 2,-20`} fill="none" stroke="#94a3b8" strokeWidth={4} strokeLinecap="round" strokeDasharray="0.1 6" opacity={0.9} />
        <line x1={chainTipX - rx} y1={yP - y - 12} x2={chainTipX - rx} y2={yP - y - 4} stroke={C.neo} strokeWidth={2} opacity={trans < 1 ? 0.8 : 0} />
      </g>,
    );
    if (stalled) {
      const upf = phase(t, tArrive + 0.3, Math.max(tArrive + 1.3, T_TRANS + 0.4), easeOut);
      ribo.push(
        <g key="upf1" opacity={upf * (1 - dissolve)} transform={`translate(${lerp(rx + 160, rx + 44, upf)}, ${y - 34 + 16 * (1 - upf)})`}>
          <circle r={11} fill="url(#sc-upf1)" />
          <text y={-16} textAnchor="middle" fontSize={9.5} fill={C.upf1} fontWeight={700}>UPF1</text>
        </g>,
      );
    }
  }

  // ---- protein chain with domains ----
  const chain: JSX.Element[] = [];
  const shown = trans > 0 ? Math.min(protLen, aaMade) : 0;
  if (shown > 0 && iso.cdsStart != null) {
    const L = shown * aaPx;
    chain.push(<line key="bb" x1={chainX0} y1={yP} x2={chainX0 + L} y2={yP} stroke={C.chain} strokeWidth={9} strokeLinecap="round" opacity={1 - dissolve} />);
    chain.push(<line key="beads" x1={chainX0 + 2} y1={yP} x2={chainX0 + L} y2={yP} stroke="#94a3b8" strokeWidth={5} strokeLinecap="round" strokeDasharray="0.1 6" opacity={0.9 * (1 - dissolve)} />);
    // domains: canonical coordinates; on the observed row map through the diff
    const diff = story.diff;
    story.domains.forEach((ds, di) => {
      const colr = DOMAIN_COLORS[di % DOMAIN_COLORS.length];
      const segs = domainSegments(ds, isObserved, diff);
      segs.forEach((sg, k) => {
        const a = Math.min(sg[0], shown), b = Math.min(sg[1], shown);
        if (b <= a) return;
        const w = (b - a) * aaPx;
        const state = isObserved ? ds.state : 'intact';
        chain.push(
          <g key={`d${di}-${k}`} opacity={1 - dissolve} {...hov([`${ds.domain.description || ds.domain.id}`, `${ds.domain.type} ${ds.domain.id} · aa ${ds.domain.start}–${ds.domain.end}`, isObserved ? `${ds.state}: ${ds.note}` : 'MANE protein'])}>
            <title>{`${ds.domain.description || ds.domain.id} · aa ${ds.domain.start}–${ds.domain.end}${isObserved ? ` · ${ds.state}: ${ds.note}` : ''}`}</title>
            <rect x={chainX0 + a * aaPx} y={yP - 9} width={w} height={18} rx={6} fill={state === 'lost' ? 'url(#sc-hatch-lost)' : colr} opacity={state === 'disrupted' ? 0.55 : 0.95} stroke={state === 'intact' ? 'none' : C.lost} strokeWidth={1.2} strokeDasharray={state === 'disrupted' ? '3 2' : undefined} />
            {k === 0 && w > 28 && <text x={chainX0 + a * aaPx + 4} y={yP - 13} fontSize={9} fill={colr} fontWeight={600}>{(ds.domain.description || ds.domain.id).slice(0, 26)}</text>}
          </g>,
        );
      });
    });
    // frameshift neopeptide / inserted residues on the observed row
    if (isObserved && (diff.kind === 'frameshift' || diff.kind.startsWith('in_frame_ins') || diff.kind === 'in_frame_indel') && diff.inserted.length) {
      const a = Math.min(diff.prefix, shown), b = Math.min(diff.prefix + diff.inserted.length, shown);
      if (b > a) chain.push(<line key="neo" x1={chainX0 + a * aaPx} y1={yP} x2={chainX0 + b * aaPx} y2={yP} stroke={C.neo} strokeWidth={9} strokeLinecap="round" strokeDasharray="4 3" opacity={1 - dissolve}
        {...hov([diff.kind === 'frameshift' ? 'Frameshifted neopeptide' : 'Inserted residues', diff.inserted.length > 40 ? diff.inserted.slice(0, 40) + '…' : diff.inserted, diff.text])} />);
    }
    // in-frame deletion: a notch on the observed chain, a hatched ghost on the canonical one
    if (isObserved && (diff.kind === 'in_frame_deletion' || diff.kind === 'in_frame_indel') && shown > diff.prefix) {
      const nx = chainX0 + diff.prefix * aaPx + (diff.inserted.length * aaPx) / 2;
      chain.push(<path key="notch" d={`M${nx - 6},${yP - 12} L${nx},${yP - 3} L${nx + 6},${yP - 12}`} fill="none" stroke={C.lost} strokeWidth={2} strokeLinecap="round" opacity={1 - dissolve} />);
    }
    if (!isObserved && diff.lost[1] > diff.lost[0] && shown >= protLen) {
      const a = diff.lost[0], b = Math.min(diff.lost[1], protLen);
      chain.push(<rect key="ghost" x={chainX0 + a * aaPx} y={yP - 12} width={(b - a) * aaPx} height={24} rx={6} fill="url(#sc-hatch-lost)" opacity={0.6 * outcome} stroke={C.lost} strokeWidth={1}
        {...hov([`Residues ${a + 1}–${b} of the MANE protein are ${diff.kind === 'frameshift' || diff.kind === 'truncation' ? 'not made' : 'deleted'}`, diff.text])} />);
    }
    // stop marker
    if (shown >= protLen) {
      const endX = chainX0 + protLen * aaPx;
      chain.push(<g key="stop" opacity={1 - dissolve} {...hov(isObserved && story.verdict.ptc != null ? [`Premature termination codon`, `codon ${story.verdict.ptcCodon}, ${story.verdict.distanceFromStart} nt after the AUG`, story.verdict.headline] : ['Stop codon', `${protLen} aa protein`])}><circle cx={endX + 6} cy={yP} r={5} fill={isObserved && story.verdict.verdict !== 'no_ptc' ? C.lost : C.ok} /><text x={endX + 6} y={yP + 3.5} textAnchor="middle" fontSize={8} fill="#fff" fontWeight={700}>*</text></g>);
      chain.push(<text key="len" x={endX + 16} y={yP + 4} fontSize={10} fill={C.muted} opacity={1 - dissolve}>{protLen} aa{isObserved && degraded ? ' fragment' : ''}</text>);
    }
  }

  return (
    <g>
      <text x={18} y={y - 10} fontSize={12} fontWeight={700} fill={C.text}>{label[0]}</text>
      <text x={18} y={y + 6} fontSize={10} fill={C.muted}>{label[1]}</text>
      <text x={18} y={yP + 4} fontSize={10} fill={C.muted}>protein</text>
      <line x1={X0} y1={y} x2={right} y2={y} stroke={C.intron} strokeWidth={1.2} opacity={appear * (1 - dissolve) * (mrna > 0 ? 0.35 : 0)} />
      {introns}
      {blocks}
      {deco}
      {chain}
      {ribo}
      {degraded && (dissolve > 0.6 || explore) && <text x={X0} y={y - 32} fontSize={12} fontWeight={700} fill={C.lost} opacity={explore ? 1 : (dissolve - 0.6) / 0.4}>⛔ transcript degraded by NMD</text>}
    </g>
  );
}

/** Amino-acid ranges (0-based half-open, in the row's own protein coordinates) drawn for a domain. */
function domainSegments(ds: DomainStatus, observed: boolean, diff: SpliceStory['diff']): [number, number][] {
  const s = ds.domain.start - 1, e = ds.domain.end;
  if (!observed || diff.kind === 'identical') return [[s, e]];
  const [a, b] = diff.lost;
  const shift = diff.inserted.length - (b - a);
  if (diff.kind === 'frameshift' || diff.kind === 'truncation' || diff.kind === 'no_protein') {
    return e <= a ? [[s, e]] : s >= a ? [] : [[s, a]];
  }
  // in-frame: keep the parts outside [a, b), shifting what lies after b
  const out: [number, number][] = [];
  if (s < a) out.push([s, Math.min(e, a)]);
  if (e > b) out.push([Math.max(s, b) + shift, e + shift]);
  return out;
}

// ======================== Caption ========================

function Caption({ t, story }: { t: number; story: SpliceStory }) {
  const lines: [number, string][] = [
    [0, 'Pre-mRNA: exons (numbered) and introns of the MANE transcript; the observed transcript starts from the same gene.'],
    [T_PRE, 'Splicing: the spliceosome loops each intron into a lariat and joins the flanking exons.' + (story.event.kind === 'skip' ? ' The skipped exon leaves with the lariat.' : story.event.kind === 'cryptic_exon' ? ' A cryptic exon is recognised inside the intron.' : story.event.kind === 'exonic_site' ? ' A cryptic splice site inside the exon trims it.' : story.event.kind === 'intronic_site' ? ' A cryptic site in the intron retains part of it.' : '')],
    [T_JOIN, 'Mature mRNA: 5′ cap, poly(A) tail, and an exon junction complex (EJC) deposited ~20 nt upstream of every splice junction.'],
    [T_MRNA, 'Translation: the ribosome reads from the AUG and displaces each EJC it passes; the protein grows with its domains.'],
    [T_TRANS, story.verdict.text],
  ];
  const cur = [...lines].reverse().find(l => t >= l[0])!;
  const fade = phase(t, cur[0], cur[0] + 0.5, easeOut);
  return (
    <g opacity={fade}>
      <rect x={20} y={VH - 62} width={VW - 40} height={46} rx={10} fill="#ffffff" opacity={0.85} />
      <foreignObject x={30} y={VH - 58} width={VW - 60} height={40}>
        <div style={{ font: `12.5px ${FONT}`, color: '#1f2937', lineHeight: '19px' }}>{cur[1]}</div>
      </foreignObject>
    </g>
  );
}

// ======================== Legend of the still picture ========================

function Legend({ story, color }: { story: SpliceStory; color: string }) {
  const y = VH - 30;
  const items: { w: number; el: JSX.Element }[] = [];
  const add = (w: number, el: JSX.Element) => items.push({ w, el });
  add(158, <g><rect x={0} y={y - 8} width={16} height={16} rx={3} fill={C.exon} /><rect x={18} y={y - 8} width={16} height={16} rx={3} fill={color} /><text x={40} y={y + 4} fontSize={10} fill={C.muted}>exon (MANE / observed)</text></g>);
  add(60, <g><rect x={0} y={y - 4} width={16} height={8} rx={2} fill={C.utr} /><text x={22} y={y + 4} fontSize={10} fill={C.muted}>UTR</text></g>);
  if (story.event.kind === 'cryptic_exon' || story.event.kind === 'intronic_site') add(120, <g><rect x={0} y={y - 8} width={16} height={16} rx={3} fill="url(#sc-hatch)" stroke={C.cryptic} /><text x={22} y={y + 4} fontSize={10} fill={C.muted}>cryptic / retained intron</text></g>);
  if (story.event.affected.length) add(92, <g><rect x={0} y={y - 8} width={16} height={16} rx={3} fill="none" stroke={C.lost} strokeDasharray="3 2" /><text x={22} y={y + 4} fontSize={10} fill={C.muted}>affected exon</text></g>);
  add(58, <g><circle cx={6} cy={y} r={5.5} fill={C.ejc} /><text x={16} y={y + 4} fontSize={10} fill={C.muted}>EJC</text></g>);
  add(70, <g><rect x={0} y={y - 7} width={20} height={14} rx={5} fill={DOMAIN_COLORS[0]} /><text x={26} y={y + 4} fontSize={10} fill={C.muted}>domain</text></g>);
  add(78, <g><rect x={0} y={y - 7} width={20} height={14} rx={5} fill={DOMAIN_COLORS[0]} opacity={0.55} stroke={C.lost} strokeDasharray="3 2" /><text x={26} y={y + 4} fontSize={10} fill={C.muted}>disrupted</text></g>);
  add(60, <g><rect x={0} y={y - 7} width={20} height={14} rx={5} fill="url(#sc-hatch-lost)" stroke={C.lost} /><text x={26} y={y + 4} fontSize={10} fill={C.muted}>lost</text></g>);
  if (story.diff.kind === 'frameshift' || story.diff.inserted.length) add(88, <g><line x1={0} y1={y} x2={20} y2={y} stroke={C.neo} strokeWidth={8} strokeLinecap="round" strokeDasharray="4 3" /><text x={26} y={y + 4} fontSize={10} fill={C.muted}>{story.diff.kind === 'frameshift' ? 'neopeptide' : 'inserted aa'}</text></g>);
  add(90, <g><circle cx={6} cy={y} r={5} fill={C.lost} /><text x={6} y={y + 3.5} textAnchor="middle" fontSize={8} fill="#fff" fontWeight={700}>*</text><text x={16} y={y + 4} fontSize={10} fill={C.muted}>stop / PTC</text></g>);
  let x = 24;
  return (
    <g>
      <rect x={16} y={VH - 48} width={VW - 32} height={36} rx={8} fill="#ffffff" opacity={0.96} stroke="#e5e7eb" />
      {items.map((it, i) => { const el = <g key={i} transform={`translate(${x}, 0)`}>{it.el}</g>; x += it.w + 18; return el; })}
    </g>
  );
}

// ======================== Amino-acid strip and domain list ========================

function Details({ story }: { story: SpliceStory }) {
  const { diff, canonical, aberrant, verdict } = story;
  const from = Math.max(0, diff.prefix - 15), to = diff.prefix + 30;
  const canon = canonical.protein, ab = aberrant.protein;
  const cell = (ch: string, cls: string, key: number) => <span key={key} className={`inline-block w-[13px] text-center ${cls}`}>{ch}</span>;
  const rowC: JSX.Element[] = [], rowA: JSX.Element[] = [];
  for (let i = from; i < to; i++) {
    const inLost = i >= diff.lost[0] && i < diff.lost[1];
    if (i < canon.length) rowC.push(cell(canon[i], inLost ? 'text-red-600 line-through' : 'text-gray-700', i));
    else if (i === canon.length) rowC.push(cell('*', 'text-emerald-600 font-bold', i));
    else rowC.push(cell(' ', '', i));
    const inNew = i >= diff.prefix && i < diff.prefix + diff.inserted.length;
    if (i < ab.length) rowA.push(cell(ab[i], inNew ? 'text-amber-600 font-bold' : 'text-gray-700', i));
    else if (i === ab.length) rowA.push(cell(aberrant.stopIndex != null ? '*' : '…', verdict.verdict === 'no_ptc' ? 'text-emerald-600 font-bold' : 'text-red-600 font-bold', i));
    else rowA.push(cell(' ', '', i));
  }
  return (
    <div className="px-5 pb-4 grid grid-cols-[1fr_320px] gap-4 text-xs">
      <div>
        <div className="text-[10.5px] font-semibold text-gray-700 mb-1">Amino acids around the change (residues {from + 1}–{to}) · {diff.text}</div>
        <div className="font-mono text-[12px] leading-5 overflow-x-auto">
          <div><span className="inline-block w-20 text-gray-500 text-[10px]">MANE</span>{rowC}</div>
          <div><span className="inline-block w-20 text-gray-500 text-[10px]">observed</span>{rowA}</div>
        </div>
        <div className="text-[10.5px] text-gray-500 mt-1">
          {verdict.ptc != null ? `PTC at codon ${verdict.ptcCodon}, ${verdict.distanceFromStart} nt after the AUG` : 'no premature stop'}
          {verdict.distanceToLastJunction != null ? ` · ${verdict.distanceToLastJunction} nt upstream of the last exon–exon junction` : ''}
          {verdict.utr3Length != null ? ` · 3′ UTR ${verdict.utr3Length.toLocaleString('en-US')} nt${verdict.longUtr ? ' (long: EJC-independent NMD possible)' : ''}` : ''}
          {' · '}canonical protein {canonical.protein.length} aa, observed {aberrant.protein.length} aa
        </div>
      </div>
      <div>
        <div className="text-[10.5px] font-semibold text-gray-700 mb-1">Protein domains (UniProt / Pfam via UCSC)</div>
        {story.domains.length === 0 && <div className="text-gray-400">no domain annotation</div>}
        <ul className="space-y-0.5">
          {story.domains.map((d, i) => (
            <li key={i} className="flex items-center gap-2">
              <span className="inline-block w-3 h-3 rounded-sm" style={{ background: DOMAIN_COLORS[i % DOMAIN_COLORS.length] }} />
              <span className="flex-1 truncate" title={`${d.domain.id} · aa ${d.domain.start}–${d.domain.end}`}>{d.domain.description || d.domain.id} <span className="text-gray-400">{d.domain.start}–{d.domain.end}</span></span>
              <span className={`px-1.5 rounded-full text-[10px] font-semibold ${d.state === 'intact' ? 'bg-emerald-100 text-emerald-700' : d.state === 'lost' ? 'bg-red-100 text-red-700' : 'bg-amber-100 text-amber-700'}`} title={d.note}>{d.state}</span>
            </li>
          ))}
        </ul>
      </div>
    </div>
  );
}
