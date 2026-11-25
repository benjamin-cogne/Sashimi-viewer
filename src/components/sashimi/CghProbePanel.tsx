/**
 * The probes' log ratios around the call drawn (CGH_Log_Ratio.txt of its patient), above the browser: one dot per probe
 * of the visible window (black within ±0.2, red above, green below), a moving median, the call shaded with its mean log
 * ratio. Its own window (zoom / pan buttons, slider, drag, Ctrl / ⌘ + wheel): the SV browser's window is internal to it,
 * so the two are not synced. Without a call (`at`), the same drawing is the whole genome panel's chromosome view, with the
 * chromosome's call bars above (`calls`, a click opens one through `onCall`). Under the plot, the chromosome's ideogram
 * (`getBands`: the cytogenetic bands of the build, drawn with the SV browser's stains) with the window boxed; a click or a
 * drag on it moves the window.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { cgh, useCgh } from './cghStore';
import { STAIN } from './SvGenomeBrowser';
import type { Cytoband } from './svNomenclature';
import { callLines, CghTip, clipMark, lrAxis, lrLabel, LrScale, BTN } from './CghGenome';
import type { CghCnv } from './cghIntervals';
import { lowerBound, meanIn, type CghProbeChrom } from './cghProbes';
import type { ProbeSeg } from './cghSegment';

const H0 = 190, L = 34, R = 8, T0 = 8, B = 20, THR = 0.2, MIN_SPAN = 20_000;
const fmtPos = (bp: number, step: number) => (step >= 1e6 ? `${(bp / 1e6).toFixed(0)} Mb` : step >= 1e5 ? `${(bp / 1e6).toFixed(1)} Mb` : step >= 1e4 ? `${(bp / 1e6).toFixed(2)} Mb` : `${(bp / 1e3).toFixed(0)} kb`);

/** the margin each side of a call: half its size, at least 100 kb (this panel and the browser below open on it) */
export const callMargin = (c: CghCnv) => Math.max((c.end - c.start) / 2, 100_000);
const around = (c: CghCnv): [number, number] => { const pad = callMargin(c); return [Math.max(0, c.start - pad), c.end + pad]; };

/** the non-NaN probes of [from, to], with k/2 more on each side for the median's edges */
function visible(p: CghProbeChrom, from: number, to: number, k: number) {
  let i0 = lowerBound(p.start, from), extra = 0;
  while (i0 > 0 && extra < k >> 1) { i0--; if (!Number.isNaN(p.lr[i0])) extra++; }
  const idx: number[] = [];
  let after = 0;
  for (let i = i0; i < p.start.length; i++) {
    if (Number.isNaN(p.lr[i])) continue;
    idx.push(i);
    if (p.start[i] > to && ++after > k >> 1) break;
  }
  return idx;
}
/** the running median of k probes (k odd), centred */
function runningMedian(v: number[], k: number): number[] {
  const h = k >> 1, out = new Array<number>(v.length), w: number[] = [];
  for (let i = 0; i < v.length; i++) {
    w.length = 0;
    for (let j = Math.max(0, i - h); j <= Math.min(v.length - 1, i + h); j++) w.push(v[j]);
    w.sort((a, b) => a - b);
    out[i] = w.length % 2 ? w[w.length >> 1] : (w[w.length / 2 - 1] + w[w.length / 2]) / 2;
  }
  return out;
}

/** `segs`: the chromosome's segments from the probe calls (cghSegment.ts), drawn as a thin purple line */
type Bar = { x: number; y: number; w: number; h: number; c: CghCnv };
export function CghProbePanel({ c, at, chromLen, segs, calls, onCall, getBands }: { c?: CghCnv; at?: { sampleKey: string; chrom: string }; chromLen: number; segs?: ProbeSeg[]; calls?: CghCnv[]; onCall?: (c: CghCnv) => void; getBands?: (chrom: string) => Promise<Cytoband[]> }) {
  const sk = c?.sampleKey ?? at?.sampleKey ?? '', chrom = c?.chrom ?? at?.chrom ?? '';
  // call bars above the plot: file calls (loss red, gain blue), then probe calls (purple)
  const T = calls ? T0 + 22 : T0, H = H0 + T - T0;
  // read on first view (61k probes: some tens of ms), then kept
  const gen = cgh.probeGen(sk), p = useMemo(() => cgh.probes(sk)?.byChrom.get(chrom) ?? null, [sk, chrom, gen]); // eslint-disable-line react-hooks/exhaustive-deps
  const [win, setWin] = useState<[number, number]>(() => (c ? around(c) : [0, chromLen]));
  const bars = useRef<Bar[]>([]), drag = useRef<{ x: number; from: number; to: number; moved: boolean } | null>(null);
  const [k, setK] = useState(7);
  const [w, setW] = useState(0);
  const [hover, setHover] = useState<{ x: number; y: number; i: number } | null>(null);
  const [callTip, setCallTip] = useState<{ x: number; y: number } | null>(null);
  const [barTip, setBarTip] = useState<{ x: number; y: number; c: CghCnv } | null>(null);
  const s = useCgh(), CLIP = s.view.scale;
  const box = useRef<HTMLDivElement>(null), cv = useRef<HTMLCanvasElement>(null);
  const bands = useBands(getBands, chrom);
  useEffect(() => {
    const el = box.current; if (!el) return;
    const ro = new ResizeObserver(() => setW(el.clientWidth)); ro.observe(el); setW(el.clientWidth);
    return () => ro.disconnect();
  }, []);
  const len = Math.max(chromLen, p ? p.end[p.end.length - 1] ?? 0 : 0);
  const [from, to] = win, span = to - from;
  const xOf = (bp: number) => L + ((bp - from) / span) * (w - L - R);
  const yOf = (v: number) => T + ((CLIP - Math.max(-CLIP, Math.min(CLIP, v))) / (2 * CLIP)) * (H - T - B);
  const cs = c?.start ?? 0, ce = c?.end ?? 0;
  const fileMean = c?.logRatio ?? null, probeMean = useMemo(() => (p && c ? meanIn(p, cs, ce) : { mean: null, n: 0 }), [p, cs, ce]); // eslint-disable-line react-hooks/exhaustive-deps
  const idx = useMemo(() => (p ? visible(p, from, to, k) : []), [p, from, to, k]);
  const med = useMemo(() => (p ? runningMedian(idx.map(i => p.lr[i]), k) : []), [p, idx, k]);

  useEffect(() => {
    const el = cv.current; if (!el || !w || !p) return;
    const t0 = performance.now(), dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    g.font = '10px Inter, ui-sans-serif, system-ui, sans-serif'; g.textBaseline = 'middle';
    // the call: shaded, its mean log ratio a segment
    const a = c ? Math.max(L, xOf(cs)) : 0, b = c ? Math.min(w - R, xOf(ce)) : 0;
    if (b > a) { g.fillStyle = 'rgba(250, 204, 21, 0.16)'; g.fillRect(a, T, Math.max(1, b - a), H - T - B); }
    // grid: 0 and the axis' lines; the thresholds ±0.2 dashed
    const ax = lrAxis(CLIP);
    for (const v of ax.lines) {
      g.strokeStyle = v === 0 ? '#94a3b8' : '#e2e8f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(L, Math.round(yOf(v)) + 0.5); g.lineTo(w - R, Math.round(yOf(v)) + 0.5); g.stroke();
    }
    g.setLineDash([4, 3]); g.strokeStyle = '#64748b';
    for (const v of [THR, -THR]) { g.beginPath(); g.moveTo(L, Math.round(yOf(v)) + 0.5); g.lineTo(w - R, Math.round(yOf(v)) + 0.5); g.stroke(); }
    g.setLineDash([]);
    g.fillStyle = '#64748b'; g.textAlign = 'right';
    for (const v of ax.labels) g.fillText(lrLabel(v), L - 4, Math.min(H - B - 4, Math.max(T + 4, yOf(v))));
    // x ticks: a round step giving 4 to 10 ticks
    const raw = span / 6, mag = 10 ** Math.floor(Math.log10(raw)), step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= 8) ?? 10 * mag;
    g.textAlign = 'center'; g.fillStyle = '#64748b'; g.strokeStyle = '#cbd5e1';
    for (let t = Math.ceil(from / step) * step; t <= to; t += step) {
      const x = Math.round(xOf(t)) + 0.5; g.beginPath(); g.moveTo(x, H - B); g.lineTo(x, H - B + 3); g.stroke(); g.fillText(fmtPos(t, step), x, H - B + 10);
    }
    g.strokeStyle = '#cbd5e1'; g.strokeRect(L + 0.5, T + 0.5, w - L - R - 1, H - T - B - 1);
    g.save(); g.beginPath(); g.rect(L, T, w - L - R, H - T - B); g.clip();
    // probes: smaller dots when dense
    const inside = idx.filter(i => p.start[i] >= from && p.start[i] <= to).length, r = inside > (w - L - R) * 1.5 ? 1.2 : inside > (w - L - R) / 3 ? 1.8 : 2.4;
    for (const i of idx) {
      const v = p.lr[i], x = xOf((p.start[i] + p.end[i]) / 2);
      g.fillStyle = v > THR ? '#dc2626' : v < -THR ? '#16a34a' : '#111827';
      g.beginPath();
      // beyond the range: a triangle at the edge, pointing out
      if (v > CLIP) clipMark(g, x, T + 0.5, 'u'); else if (v < -CLIP) clipMark(g, x, H - B - 0.5, 'd'); else g.arc(x, yOf(v), r, 0, 2 * Math.PI);
      g.fill();
    }
    // the moving median (dark blue), then the call's mean on top
    g.strokeStyle = '#1e3a8a'; g.lineWidth = 1.6; g.beginPath();
    idx.forEach((i, j) => { const x = xOf((p.start[i] + p.end[i]) / 2), y = yOf(med[j]); if (j) g.lineTo(x, y); else g.moveTo(x, y); });
    g.stroke();
    const m = fileMean ?? probeMean.mean;
    if (m != null && b > a) { g.strokeStyle = '#ea580c'; g.lineWidth = 3; g.beginPath(); g.moveTo(a, yOf(m)); g.lineTo(b, yOf(m)); g.stroke(); }
    // the segments' means on top (over the median and the call's mean): a step where the probe calls cut
    if (segs) {
      g.strokeStyle = '#a21caf'; g.lineWidth = 1.5; g.beginPath();
      for (const sg of segs) { if (sg.end < from || sg.start > to) continue; const y = yOf(sg.mean); g.moveTo(Math.max(L, xOf(sg.start)), y); g.lineTo(Math.min(w - R, xOf(sg.end)), y); }
      g.stroke();
    }
    g.restore();
    // the call bars (chromosome view)
    const out: Bar[] = [];
    if (calls) for (const x of calls) {
      if (x.end < from || x.start > to) continue;
      const a0 = Math.max(L, xOf(x.start)), b0 = Math.min(w - R, Math.max(xOf(x.start) + 2, xOf(x.end)));
      const bar = x.fromProbes ? { x: a0, y: 16, w: b0 - a0, h: 7, c: x } : { x: a0, y: 4, w: b0 - a0, h: 9, c: x };
      g.fillStyle = x.fromProbes ? '#a21caf' : x.kind === 'loss' ? '#dc2626' : '#2563eb'; g.fillRect(bar.x, bar.y, bar.w, bar.h);
      if (x.id === s.shown) { g.strokeStyle = '#0f172a'; g.lineWidth = 1.5; g.strokeRect(bar.x - 1, bar.y - 1, bar.w + 2, bar.h + 2); }
      out.push(bar);
    }
    bars.current = out;
    el.setAttribute('data-drawn', String(inside)); el.setAttribute('data-ms', (performance.now() - t0).toFixed(0));
  }, [w, p, from, to, idx, med, cs, ce, fileMean, probeMean.mean, segs, CLIP, calls, s.shown]); // eslint-disable-line react-hooks/exhaustive-deps

  // Ctrl / ⌘ + wheel: zoom around the pointer (a plain wheel still scrolls the page)
  const live = useRef({ from, to, len, w }); live.current = { from, to, len, w };
  useEffect(() => {
    const el = box.current; if (!el) return;
    const f = (e: WheelEvent) => {
      if (!e.ctrlKey && !e.metaKey) return;
      e.preventDefault();
      const { from, to, len, w } = live.current, x = e.clientX - el.getBoundingClientRect().left;
      const bp = from + ((Math.max(L, Math.min(w - R, x)) - L) / (w - L - R)) * (to - from), k = e.deltaY > 0 ? 1.25 : 0.8;
      let a = bp - (bp - from) * k, b = bp + (to - bp) * k;
      if (b - a < MIN_SPAN) return; if (b - a > len) { a = 0; b = len; }
      const sh = a < 0 ? -a : b > len ? len - b : 0; setWin([Math.round(a + sh), Math.round(b + sh)]);
    };
    el.addEventListener('wheel', f, { passive: false });
    return () => el.removeEventListener('wheel', f);
  }, [!!p]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!p) return <div className="mx-5 mt-2 text-[11px] text-slate-500">Log ratio (probes): no probe on {chrom} in this patient's probe file.</div>;
  /** a window kept on the chromosome, at least MIN_SPAN wide */
  const go = (a: number, b: number) => { let s = Math.max(MIN_SPAN, b - a); s = Math.min(s, len); let x = Math.max(0, Math.min(len - s, (a + b) / 2 - s / 2)); x = Math.round(x); setWin([x, x + Math.round(s)]); };
  const mid = (from + to) / 2;
  const barAt = (x: number, y: number) => bars.current.find(b => x >= b.x - 1 && x <= b.x + b.w + 1 && y >= b.y - 1 && y <= b.y + b.h + 1);
  const onMove = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    setHover(null); setCallTip(null); setBarTip(null);
    // dragging: the window follows the pointer
    const d = drag.current;
    if (d && e.buttons & 1) { const dx = x - d.x; if (Math.abs(dx) > 3) d.moved = true; if (d.moved) { const sh = -dx * (d.to - d.from) / (w - L - R); go(d.from + sh, d.to + sh); } return; }
    const bb = barAt(x, y);
    if (bb) { setBarTip({ x, y, c: bb.c }); return; }
    if (x < L || x > w - R || y < T || y > H - B) return;
    // the nearest probe within 8 px (those of ± 8 px across); else, inside the call's shading, the call's tip
    const bpx = span / (w - L - R), bp = from + (x - L) * bpx;
    let lo = 0, hi = idx.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (p.start[idx[m]] < bp - 8 * bpx) lo = m + 1; else hi = m; }
    let best = -1, bd = 64;
    for (let j = lo; j < idx.length && p.start[idx[j]] <= bp + 8 * bpx; j++) {
      const i = idx[j], dx = xOf((p.start[i] + p.end[i]) / 2) - x, dy = yOf(p.lr[i]) - y;
      if (dx * dx + dy * dy < bd) { bd = dx * dx + dy * dy; best = i; }
    }
    if (best >= 0) setHover({ x: xOf((p.start[best] + p.end[best]) / 2), y: yOf(p.lr[best]), i: best });
    else if (c && x >= xOf(cs) && x <= Math.max(xOf(cs) + 2, xOf(ce))) setCallTip({ x, y });
  };
  const onDown = (e: React.MouseEvent) => { const rc = e.currentTarget.getBoundingClientRect(); drag.current = { x: e.clientX - rc.left, from, to, moved: false }; };
  const onClick = (e: React.MouseEvent) => {
    const moved = drag.current?.moved; drag.current = null; if (moved) return;
    const rc = e.currentTarget.getBoundingClientRect(), bb = barAt(e.clientX - rc.left, e.clientY - rc.top);
    if (bb && onCall) onCall(bb.c);
  };
  return (
    <div className="mx-5 mt-2" {...(c ? { 'data-cgh-probes': '' } : { 'data-cgh-chrom-probes': '' })}>
      <div className="flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] text-slate-600">
        {c && <span className="font-semibold text-slate-900 text-xs">Log ratio (probes)</span>}
        <span className="tabular-nums">{chrom}:{Math.round(from).toLocaleString('en-US')}-{Math.round(to).toLocaleString('en-US')}</span>
        <span className="inline-flex gap-1">
          <button className={BTN} title="Zoom out" aria-label="Zoom out" onClick={() => go(mid - span, mid + span)}>⊖</button>
          <button className={BTN} title="Zoom in" aria-label="Zoom in" onClick={() => go(mid - span / 4, mid + span / 4)}>⊕</button>
          <button className={BTN} title="Left by half the window" aria-label="Pan left" onClick={() => go(from - span / 2, to - span / 2)}>←</button>
          <button className={BTN} title="Right by half the window" aria-label="Pan right" onClick={() => go(from + span / 2, to + span / 2)}>→</button>
          <button className={BTN} onClick={() => setWin([0, len])}>whole chromosome</button>
          {c && <button className={BTN} title="Back to the call, centred, ± half its size (at least ± 100 kb)" onClick={() => setWin(around(c))}>call</button>}
        </span>
        <input type="range" min={0} max={Math.max(0, Math.round(len - span))} step={Math.max(1, Math.round(span / 50))} value={Math.round(from)} disabled={span >= len} onChange={e => { const a = Number(e.target.value); go(a, a + span); }}
          aria-label="Move along the chromosome" title="Move the window along the chromosome (also: drag the plot; Ctrl / ⌘ + wheel zooms)" className="w-32 accent-slate-500" data-cgh-pan />
        <label className="inline-flex items-center gap-1" title="Moving median over this many probes (dark blue line)">median of
          <select value={k} onChange={e => setK(Number(e.target.value))} aria-label="Probes in the moving median" className="h-[22px] rounded border border-slate-200 bg-white text-[11px] outline-none cursor-pointer">{[3, 7, 15].map(n => <option key={n}>{n}</option>)}</select> probes</label>
        {c && <span className="text-slate-500" title="The call's mean log ratio as the file gives it (orange segment), and the mean of the probes it spans">
          call: {fileMean != null ? <>{c.fromProbes ? 'segment' : 'file'} <b className="text-slate-800">{fileMean.toFixed(2)}</b></> : 'no log ratio in the file'}{probeMean.mean != null ? <> · probes' mean <b className="text-slate-800">{probeMean.mean.toFixed(2)}</b> ({probeMean.n})</> : ''}</span>}
        <LrScale />
        <span className="text-slate-400">· {c ? 'own window (not synced with the browser) · ' : ''}dashed: ±0.2 · clipped at ±{CLIP} (▲▼ beyond){segs ? <> · <span className="text-fuchsia-700">purple</span>: probe-call segments</> : ''}</span>
      </div>
      <div ref={box} className={`relative mt-1 select-none ${barTip ? 'cursor-pointer' : 'cursor-grab'}`} onMouseMove={onMove} onMouseDown={onDown} onClick={onClick} onMouseLeave={() => { setHover(null); setCallTip(null); setBarTip(null); drag.current = null; }} data-cgh-probe-plot>
        <canvas ref={cv} style={{ width: '100%', height: H, display: 'block' }} />
        {hover && <div className="absolute z-10 pointer-events-none rounded border border-slate-300 bg-white/95 shadow px-1.5 py-0.5 text-[10.5px] text-slate-800 whitespace-nowrap"
          style={{ left: Math.min(hover.x + 8, Math.max(0, w - 210)), top: Math.max(0, Math.min(hover.y - 30, H - 44)) }} data-cgh-probe-tip>
          <div className="font-semibold">{p.name[hover.i] || 'probe'}</div>
          <div className="tabular-nums">{chrom}:{p.start[hover.i].toLocaleString('en-US')}-{p.end[hover.i].toLocaleString('en-US')} · log ratio <b>{p.lr[hover.i].toFixed(3)}</b>{Math.abs(p.lr[hover.i]) > CLIP ? ` (beyond ±${CLIP})` : ''}</div>
        </div>}
        {barTip && <CghTip x={barTip.x} y={barTip.y} w={w} title={`${barTip.c.sample}${barTip.c.fromProbes ? ' · probe call' : ''}`} lines={callLines(s, barTip.c)} />}
        {callTip && c && <CghTip x={callTip.x} y={callTip.y} w={w} title={`${c.sample} · ${c.fromProbes ? 'probe call (from probes)' : 'call'}`} lines={callLines(s, c)} />}
        {hover && <div className="absolute pointer-events-none rounded-full border-2 border-slate-900" style={{ left: hover.x - 4, top: hover.y - 4, width: 8, height: 8 }} />}
      </div>
      {bands && w > 0 && <Ideogram bands={bands} chrom={chrom} len={len} w={w} from={from} to={to} onGo={(a, b) => go(a, b)} />}
    </div>
  );
}

/** the chromosome's cytogenetic bands, read once per chromosome (the loader caches them); undefined until read or when they cannot be */
function useBands(get: ((chrom: string) => Promise<Cytoband[]>) | undefined, chrom: string) {
  const [got, setGot] = useState<{ chrom: string; bands: Cytoband[] } | null>(null);
  const ref = useRef(get); ref.current = get;
  useEffect(() => {
    let alive = true;
    ref.current?.(chrom).then(b => { if (alive && b.length) setGot({ chrom, bands: b }); }).catch(() => { /* no ideogram */ });
    return () => { alive = false; };
  }, [chrom]);
  return got?.chrom === chrom ? got.bands : undefined;
}

const IY = 11, IH = 12;
/** the ideogram under the plot, drawn as the SV browser's (stains, centromere triangles, band names where they fit), the
 * window boxed; a click centres the window there, a drag moves it */
function Ideogram({ bands, chrom, len, w, from, to, onGo }: { bands: Cytoband[]; chrom: string; len: number; w: number; from: number; to: number; onGo: (a: number, b: number) => void }) {
  const PW = w - L - R, ix = (bp: number) => L + (bp / len) * PW, span = to - from;
  const drag = useRef<{ x: number; from: number } | null>(null);
  const bpAt = (e: React.PointerEvent) => ((e.clientX - e.currentTarget.getBoundingClientRect().left - L) / PW) * len;
  const a = ix(from), b = Math.max(a + 3, ix(to)), clip = `cgh-ideo-${chrom}`;
  const cen = bands.filter(x => x.stain === 'acen');
  return (
    <svg width={w} height={IY + IH + 3} className="block cursor-pointer select-none" data-cgh-ideogram aria-label={`${chrom} ideogram: click or drag to move the window`}
      onPointerDown={e => {
        const bp = bpAt(e), inside = bp >= from && bp <= to;
        e.currentTarget.setPointerCapture(e.pointerId);
        drag.current = { x: e.clientX, from: inside ? from : bp - span / 2 };
        if (!inside) onGo(bp - span / 2, bp + span / 2);
      }}
      onPointerMove={e => { const d = drag.current; if (!d || !(e.buttons & 1)) return; const sh = ((e.clientX - d.x) / PW) * len; onGo(d.from + sh, d.from + sh + span); }}
      onPointerUp={() => { drag.current = null; }}>
      <title>{`${chrom}: cytogenetic bands (UCSC), the window boxed. Click or drag to move it.`}</title>
      <defs><clipPath id={clip}><rect x={L} y={IY} width={PW} height={IH} rx={IH / 2} /></clipPath></defs>
      <text x={L - 4} y={IY + IH - 2} fontSize={9} textAnchor="end" fill="#64748b">{chrom.replace(/^chr/, '')}</text>
      <g clipPath={`url(#${clip})`}>
        <rect x={L} y={IY} width={PW} height={IH} fill={STAIN.gneg} />
        {bands.map(x => x.stain === 'acen'
          ? <path key={x.name + x.start} d={x.name.startsWith('p') ? `M${ix(x.start)},${IY} L${ix(x.end)},${IY + IH / 2} L${ix(x.start)},${IY + IH} Z` : `M${ix(x.end)},${IY} L${ix(x.start)},${IY + IH / 2} L${ix(x.end)},${IY + IH} Z`} fill={STAIN.acen} />
          : <rect key={x.name + x.start} x={ix(x.start)} y={IY} width={Math.max(0.5, ix(x.end) - ix(x.start))} height={IH} fill={STAIN[x.stain] ?? STAIN.gneg} />)}
      </g>
      <rect x={L} y={IY} width={PW} height={IH} rx={IH / 2} fill="none" stroke="#94a3b8" strokeWidth={0.8} />
      {cen.length > 0 && <text x={ix((cen[0].start + cen[cen.length - 1].end) / 2)} y={IY - 2} fontSize={8} textAnchor="middle" fill="#be123c">cen</text>}
      {bands.filter(x => x.stain !== 'acen' && ix(x.end) - ix(x.start) > x.name.length * 4.6 + 6).map(x => <text key={`n${x.name}${x.start}`} x={(ix(x.start) + ix(x.end)) / 2} y={IY - 2} fontSize={8} textAnchor="middle" fill="#94a3b8">{x.name}</text>)}
      <rect x={a} y={IY - 1.5} width={b - a} height={IH + 3} fill="rgba(37, 99, 235, 0.12)" stroke="#2563eb" strokeWidth={1.5} rx={2} data-cgh-ideo-window />
    </svg>
  );
}
