/**
 * The probes' log ratios around the call drawn (CGH_Log_Ratio.txt of its patient), above the browser: one dot per probe
 * of the visible window (black within ±0.2, red above, green below), a moving median, the call shaded with its mean log
 * ratio. Its own window (zoom / pan buttons): the SV browser's window is internal to it, so the two are not synced.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { cgh, useCgh } from './cghStore';
import { callLines, CghTip } from './CghGenome';
import type { CghCnv } from './cghIntervals';
import { lowerBound, meanIn, type CghProbeChrom } from './cghProbes';
import type { ProbeSeg } from './cghSegment';

const H = 190, L = 34, R = 8, T = 8, B = 20, CLIP = 2, THR = 0.2, MIN_SPAN = 20_000;
const BTN = 'h-[22px] min-w-[24px] px-1.5 rounded border border-slate-200 bg-white text-[11px] text-slate-700 hover:border-accent-300 hover:bg-accent-50';
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
export function CghProbePanel({ c, chromLen, segs }: { c: CghCnv; chromLen: number; segs?: ProbeSeg[] }) {
  // read on first view (61k probes: some tens of ms), then kept
  const p = useMemo(() => cgh.probes(c.sampleKey)?.byChrom.get(c.chrom) ?? null, [c.sampleKey, c.chrom]);
  const [win, setWin] = useState<[number, number]>(() => around(c));
  const [k, setK] = useState(7);
  const [w, setW] = useState(0);
  const [hover, setHover] = useState<{ x: number; y: number; i: number } | null>(null);
  const [callTip, setCallTip] = useState<{ x: number; y: number } | null>(null);
  const s = useCgh();
  const box = useRef<HTMLDivElement>(null), cv = useRef<HTMLCanvasElement>(null);
  useEffect(() => {
    const el = box.current; if (!el) return;
    const ro = new ResizeObserver(() => setW(el.clientWidth)); ro.observe(el); setW(el.clientWidth);
    return () => ro.disconnect();
  }, []);
  const len = Math.max(chromLen, p ? p.end[p.end.length - 1] ?? 0 : 0);
  const [from, to] = win, span = to - from;
  const xOf = (bp: number) => L + ((bp - from) / span) * (w - L - R);
  const yOf = (v: number) => T + ((CLIP - Math.max(-CLIP, Math.min(CLIP, v))) / (2 * CLIP)) * (H - T - B);
  const fileMean = c.logRatio, probeMean = useMemo(() => (p ? meanIn(p, c.start, c.end) : { mean: null, n: 0 }), [p, c.start, c.end]);
  const idx = useMemo(() => (p ? visible(p, from, to, k) : []), [p, from, to, k]);
  const med = useMemo(() => (p ? runningMedian(idx.map(i => p.lr[i]), k) : []), [p, idx, k]);

  useEffect(() => {
    const el = cv.current; if (!el || !w || !p) return;
    const dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    g.font = '10px Inter, ui-sans-serif, system-ui, sans-serif'; g.textBaseline = 'middle';
    // the call: shaded, its mean log ratio a segment
    const a = Math.max(L, xOf(c.start)), b = Math.min(w - R, xOf(c.end));
    if (b > a) { g.fillStyle = 'rgba(250, 204, 21, 0.16)'; g.fillRect(a, T, Math.max(1, b - a), H - T - B); }
    // grid: 0, ±0.25, ±1; the thresholds ±0.2 dashed
    for (const v of [1, 0.25, 0, -0.25, -1]) {
      g.strokeStyle = v === 0 ? '#94a3b8' : '#e2e8f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(L, Math.round(yOf(v)) + 0.5); g.lineTo(w - R, Math.round(yOf(v)) + 0.5); g.stroke();
    }
    g.setLineDash([4, 3]); g.strokeStyle = '#64748b';
    for (const v of [THR, -THR]) { g.beginPath(); g.moveTo(L, Math.round(yOf(v)) + 0.5); g.lineTo(w - R, Math.round(yOf(v)) + 0.5); g.stroke(); }
    g.setLineDash([]);
    g.fillStyle = '#64748b'; g.textAlign = 'right';
    for (const v of [2, 1, 0, -1, -2]) g.fillText(v > 0 ? `+${v}` : String(v), L - 4, Math.min(H - B - 4, Math.max(T + 4, yOf(v))));
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
      const v = p.lr[i];
      g.fillStyle = v > THR ? '#dc2626' : v < -THR ? '#16a34a' : '#111827';
      g.beginPath(); g.arc(xOf((p.start[i] + p.end[i]) / 2), yOf(v), r, 0, 2 * Math.PI); g.fill();
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
    cv.current?.setAttribute('data-drawn', String(inside));
  }, [w, p, from, to, idx, med, c.start, c.end, fileMean, probeMean.mean, segs]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!p) return <div className="mx-5 mt-2 text-[11px] text-slate-500">Log ratio (probes): no probe on {c.chrom} in this patient's probe file.</div>;
  /** a window kept on the chromosome, at least MIN_SPAN wide */
  const go = (a: number, b: number) => { let s = Math.max(MIN_SPAN, b - a); s = Math.min(s, len); let x = Math.max(0, Math.min(len - s, (a + b) / 2 - s / 2)); x = Math.round(x); setWin([x, x + Math.round(s)]); };
  const mid = (from + to) / 2;
  const onMove = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    setHover(null); setCallTip(null);
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
    else if (x >= xOf(c.start) && x <= Math.max(xOf(c.start) + 2, xOf(c.end))) setCallTip({ x, y });
  };
  return (
    <div className="mx-5 mt-2" data-cgh-probes>
      <div className="flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] text-slate-600">
        <span className="font-semibold text-slate-900 text-xs">Log ratio (probes)</span>
        <span className="tabular-nums">{c.chrom}:{Math.round(from).toLocaleString('en-US')}-{Math.round(to).toLocaleString('en-US')}</span>
        <span className="inline-flex gap-1">
          <button className={BTN} title="Zoom out" aria-label="Zoom out" onClick={() => go(mid - span, mid + span)}>⊖</button>
          <button className={BTN} title="Zoom in" aria-label="Zoom in" onClick={() => go(mid - span / 4, mid + span / 4)}>⊕</button>
          <button className={BTN} title="Left by half the window" aria-label="Pan left" onClick={() => go(from - span / 2, to - span / 2)}>←</button>
          <button className={BTN} title="Right by half the window" aria-label="Pan right" onClick={() => go(from + span / 2, to + span / 2)}>→</button>
          <button className={BTN} onClick={() => setWin([0, len])}>whole chromosome</button>
          <button className={BTN} title="Back to the call, centred, ± half its size (at least ± 100 kb)" onClick={() => setWin(around(c))}>call</button>
        </span>
        <label className="inline-flex items-center gap-1" title="Moving median over this many probes (dark blue line)">median of
          <select value={k} onChange={e => setK(Number(e.target.value))} aria-label="Probes in the moving median" className="h-[22px] rounded border border-slate-200 bg-white text-[11px] outline-none cursor-pointer">{[3, 7, 15].map(n => <option key={n}>{n}</option>)}</select> probes</label>
        <span className="text-slate-500" title="The call's mean log ratio as the file gives it (orange segment), and the mean of the probes it spans">
          call: {fileMean != null ? <>{c.fromProbes ? 'segment' : 'file'} <b className="text-slate-800">{fileMean.toFixed(2)}</b></> : 'no log ratio in the file'}{probeMean.mean != null ? <> · probes' mean <b className="text-slate-800">{probeMean.mean.toFixed(2)}</b> ({probeMean.n})</> : ''}</span>
        <span className="text-slate-400">· own window (not synced with the browser) · dashed: ±0.2 · clipped at ±2{segs ? <> · <span className="text-fuchsia-700">purple</span>: probe-call segments</> : ''}</span>
      </div>
      <div ref={box} className="relative mt-1 select-none" onMouseMove={onMove} onMouseLeave={() => { setHover(null); setCallTip(null); }} data-cgh-probe-plot>
        <canvas ref={cv} style={{ width: '100%', height: H, display: 'block' }} />
        {hover && <div className="absolute z-10 pointer-events-none rounded border border-slate-300 bg-white/95 shadow px-1.5 py-0.5 text-[10.5px] text-slate-800 whitespace-nowrap"
          style={{ left: Math.min(hover.x + 8, Math.max(0, w - 210)), top: Math.max(0, Math.min(hover.y - 30, H - 44)) }} data-cgh-probe-tip>
          <div className="font-semibold">{p.name[hover.i] || 'probe'}</div>
          <div className="tabular-nums">{c.chrom}:{p.start[hover.i].toLocaleString('en-US')}-{p.end[hover.i].toLocaleString('en-US')} · log ratio <b>{p.lr[hover.i].toFixed(3)}</b></div>
        </div>}
        {callTip && <CghTip x={callTip.x} y={callTip.y} w={w} title={`${c.sample} · ${c.fromProbes ? 'probe call (from probes)' : 'call'}`} lines={callLines(s, c)} />}
        {hover && <div className="absolute pointer-events-none rounded-full border-2 border-slate-900" style={{ left: hover.x - 4, top: hover.y - 4, width: 8, height: 8 }} />}
      </div>
    </div>
  );
}
