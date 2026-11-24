/**
 * The *Whole genome* panel of the array-CGH page: one patient's probes (CGH_Log_Ratio.txt) over chromosomes 1-22, X, Y
 * side by side, with its calls as bars above; a chromosome clicked opens it in the same panel, a call bar clicked there
 * opens the call. It stands where the call window does: expanded while no call is open, a strip once one is.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { cgh, useCgh, classOf, passes, probeCnv, LR_SCALES, LR_PRESETS, type CghState } from './cghStore';
import type { CghCnv, CghFile } from './cghIntervals';
import { lowerBound } from './cghProbes';
import { chromLength } from '../../standalone/assembly';
import { estCopies } from './cghMosaic';
import { viewOf } from './cghAnnotate';
import { CghProbePanel } from './CghProbePanel';
import { LEVEL_STYLE } from './cghScore';

const CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY'];
const L = 34, R = 8, PT = 24, PH = 170, PB = PT + PH, HH = PB + 18, THR = 0.2, GAP = 2;
/** vertical layout: a row's label band, its tallest column (chr1), the axis labels under it; the call bars beside each column */
const VT = 16, VB = 16, VBARS = 12;
export const BTN = 'h-[22px] min-w-[24px] px-1.5 rounded border border-slate-200 bg-white text-[11px] text-slate-700 hover:border-accent-300 hover:bg-accent-50 disabled:opacity-40 disabled:pointer-events-none';
const fmt = (n: number) => n.toLocaleString('en-US');
/** a chromosome's plot: its position axis p0 → p1 (x, or y going down in the vertical layout), its log-ratio axis q0 → q1
 * (y from +S down to −S; vertical: x from −S to +S) */
type Seg = { ch: string; to: number; vert: boolean; p0: number; p1: number; q0: number; q1: number };
type Bar = { x: number; y: number; w: number; h: number; c: CghCnv };
const posOf = (sg: Seg, bp: number) => sg.p0 + (bp / sg.to) * (sg.p1 - sg.p0);
const clampLr = (v: number, S: number) => Math.max(-S, Math.min(S, v));
const lrOf = (sg: Seg, v: number, S: number) => sg.q0 + ((sg.vert ? clampLr(v, S) + S : S - clampLr(v, S)) / (2 * S)) * (sg.q1 - sg.q0);
const xyOf = (sg: Seg, bp: number, v: number, S: number): [number, number] => (sg.vert ? [lrOf(sg, v, S), posOf(sg, bp)] : [posOf(sg, bp), lrOf(sg, v, S)]);

/** the log-ratio axis at ±S: the values labelled, the guide lines (0 included) */
export function lrAxis(S: number) {
  const step = S <= 0.5 ? 0.25 : S <= 1 ? 0.5 : S <= 3 ? 1 : S <= 4 ? 2 : 4, set = new Set([0, S, -S]);
  for (let v = step; v < S - 1e-9; v += step) { set.add(Math.round(v * 100) / 100); set.add(-Math.round(v * 100) / 100); }
  const labels = [...set].sort((a, b) => a - b);
  return { labels, lines: [...new Set([...labels, ...(S <= 2 ? [0.25, -0.25] : [])])].filter(v => Math.abs(v) < S) };
}
export const lrLabel = (v: number) => (v > 0 ? `+${v}` : String(v));
/** a small triangle at the edge for a point beyond the range, pointing out (u / d / l / r) */
export function clipMark(g: CanvasRenderingContext2D, x: number, y: number, dir: 'u' | 'd' | 'l' | 'r') {
  const [a, b] = dir === 'u' ? [0, 1] : dir === 'd' ? [0, -1] : dir === 'l' ? [1, 0] : [-1, 0];
  g.moveTo(x, y); g.lineTo(x + a * 4 - b * 2.5, y + b * 4 - a * 2.5); g.lineTo(x + a * 4 + b * 2.5, y + b * 4 + a * 2.5); g.closePath();
}
/** the log-ratio range (±), shared by the whole genome, chromosome and probe panels, kept as a setting */
export function LrScale() {
  const S = useCgh().view.scale, i = LR_SCALES.indexOf(S);
  return <span className="inline-flex items-center gap-1" data-cgh-scale title="Log-ratio axis range, shared by the whole genome, the chromosome and the call's probe panel (kept as a setting); points beyond it are drawn at the edge as small triangles; the colours stay at ±0.2">log ratio
    <span className="inline-flex rounded border border-slate-200 overflow-hidden" role="group" aria-label="Log-ratio range">
      {LR_PRESETS.map(v => <button key={v} onClick={() => cgh.setView({ scale: v })} aria-pressed={S === v} data-scale={v} className={`h-[22px] px-1.5 text-[11px] tabular-nums ${S === v ? 'bg-accent-50 text-accent-700 font-semibold' : 'bg-white text-slate-600 hover:bg-slate-50'}`}>±{v}</button>)}
      {!LR_PRESETS.includes(S) && <span className="h-[22px] px-1.5 text-[11px] leading-[22px] bg-accent-50 text-accent-700 font-semibold tabular-nums">±{S}</span>}</span>
    <button className={BTN} disabled={i <= 0} onClick={() => cgh.setView({ scale: LR_SCALES[i - 1] })} title="Narrower log-ratio range (zoom in)" aria-label="Zoom in on the log ratio">⊕</button>
    <button className={BTN} disabled={i >= LR_SCALES.length - 1} onClick={() => cgh.setView({ scale: LR_SCALES[i + 1] })} title="Wider log-ratio range (zoom out: ±3, ±4, ±6, ±8)" aria-label="Zoom out on the log ratio">⊖</button></span>;
}

/** a call's characteristics, for the hover tips (here and in the probe panel) */
export function callLines(s: CghState, c: CghCnv): string[] {
  const size = c.sizeKb >= 1000 ? `${(c.sizeKb / 1000).toFixed(2)} Mb` : `${c.sizeKb.toFixed(1)} kb`;
  const est = c.copies == null ? estCopies(c, s.sex[c.sampleKey]) : null, lv = c.fromProbes ? null : viewOf(s, c).score?.level;
  return [
    `${c.chrom}:${fmt(c.start)}-${fmt(c.end)} · ${size} · ${c.kind}`,
    `${c.probes ?? '?'} probes${c.logRatio != null ? ` · mean log ratio ${c.logRatio.toFixed(2)}` : ''}`,
    c.copies != null ? `copies x${c.copies}${c.copiesHi != null ? `~${c.copiesHi}` : ''} (file)` : est ? `copies ${est.text}` : 'copies not given',
    c.fromProbes ? 'from probes: not in the file, not exported' : `class: ${classOf(s, c).trim() || 'none'} · priority: ${lv ? LEVEL_STYLE[lv].label : 'not read yet'}`,
  ];
}
/** a small tip beside the pointer, kept inside the panel's width */
export function CghTip({ x, y, w, title, lines }: { x: number; y: number; w: number; title: string; lines: string[] }) {
  return <div className="absolute z-10 pointer-events-none rounded border border-slate-300 bg-white/95 shadow px-1.5 py-0.5 text-[10.5px] text-slate-800 whitespace-nowrap tabular-nums"
    style={{ left: Math.max(0, Math.min(x + 10, w - 280)), top: Math.max(0, y - 34) }} data-cgh-tip>
    <div className="font-semibold">{title}</div>{lines.map((l, i) => <div key={i}>{l}</div>)}</div>;
}

/** the median of a few values (sorted in place) */
const median = (v: number[]) => { v.sort((a, b) => a - b); const h = v.length >> 1; return v.length % 2 ? v[h] : (v[h - 1] + v[h]) / 2; };

export function CghGenomePanel({ open, onToggle, onOpen }: { open: boolean; onToggle: () => void; onOpen: (c: CghCnv, f: CghFile) => void }) {
  const s = useCgh();
  const all = useMemo(() => s.files.flatMap(f => f.samples.map(x => ({ x, f }))), [s.files]);
  const [pick, setPick] = useState<string | null>(null);
  // the active tab's patient; in All patients, the one picked here
  const cur = all.find(a => a.x.key === s.sample) ?? all.find(a => a.x.key === pick) ?? all[0];
  const key = cur?.x.key ?? '', build = cur?.f.build ?? 'GRCh37';
  const [chrom, setChrom] = useState<string | null>(null);
  const [w, setW] = useState(0);
  const [hover, setHover] = useState<{ x: number; y: number; title: string; lines: string[]; dot?: [number, number] } | null>(null);
  const box = useRef<HTMLDivElement>(null), base = useRef<HTMLCanvasElement>(null), top = useRef<HTMLCanvasElement>(null), bars = useRef<Bar[]>([]);
  const S = s.view.scale, vert = s.view.layout === 'vertical' && !chrom;
  const chromPanel = !!chrom && open && !!key && cgh.hasProbes(key);
  useEffect(() => {
    const el = box.current; if (!el) return;
    const ro = new ResizeObserver(() => setW(el.clientWidth)); ro.observe(el); setW(el.clientWidth);
    return () => ro.disconnect();
  }, [open, !!cur, chromPanel]); // eslint-disable-line react-hooks/exhaustive-deps
  // read on first view only (some tens of ms for 61k probes)
  const pr = useMemo(() => {
    if (!open || !key || !cgh.hasProbes(key)) return null;
    const t = performance.now(), p = cgh.probes(key);
    return p ? { p, ms: performance.now() - t } : null;
  }, [open, key]);
  /** the chromosomes' plots and the canvas height: side by side, one chromosome, or columns (one row, two when narrow) */
  const { segs, H, rows } = useMemo<{ segs: Seg[]; H: number; rows: number }>(() => {
    if (!w) return { segs: [], H: HH, rows: 1 };
    const lens = CHROMS.map(ch => chromLength(build, ch) ?? 1);
    if (chrom) return { segs: [{ ch: chrom, to: lens[CHROMS.indexOf(chrom)], vert: false, p0: L, p1: w - R, q0: PT, q1: PB }], H: HH, rows: 1 };
    if (vert) {
      const rows = w < 1000 ? 2 : 1, n = Math.ceil(CHROMS.length / rows), colW = (w - 8) / n, ph = rows === 1 ? 400 : 290, k = ph / Math.max(...lens);
      // each row as tall as its longest chromosome
      const tops = [0]; for (let r = 0; r < rows; r++) tops.push(tops[r] + VT + Math.max(...lens.slice(r * n, (r + 1) * n)) * k + VB);
      return {
        rows, H: Math.ceil(tops[rows]) + 2,
        segs: CHROMS.map((ch, i) => {
          const r = Math.floor(i / n), x = 4 + (i % n) * colW + 2, y = tops[r] + VT;
          return { ch, to: lens[i], vert: true, p0: y, p1: y + lens[i] * k, q0: x, q1: x + colW - 4 - VBARS };
        }),
      };
    }
    const tot = lens.reduce((a, b) => a + b, 0), k = (w - L - R - GAP * (CHROMS.length - 1)) / tot;
    let x = L;
    return { H: HH, rows: 1, segs: CHROMS.map((ch, i) => { const sg = { ch, to: lens[i], vert: false, p0: x, p1: x + lens[i] * k, q0: PT, q1: PB }; x = sg.p1 + GAP; return sg; }) };
  }, [w, chrom, build, vert]);
  const calls = cur ? cur.x.cnvs.filter(c => passes(s, c)) : [];
  const events = cur ? (cgh.probeCallsOf(key)?.events ?? []).map(e => probeCnv(e, key, cur.x.name)) : [];

  // the probes: one mark per pixel cell (the same picture, bounded by the pixels), a median line; redrawn only on resize / view / scale
  useEffect(() => {
    const el = base.current; if (!el || !w || !open) return;
    const t0 = performance.now(), dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    g.font = '10px Inter, ui-sans-serif, system-ui, sans-serif'; g.textBaseline = 'middle';
    const ax = lrAxis(S);
    if (vert) {
      g.textAlign = 'center';
      segs.forEach((sg, i) => {
        g.fillStyle = '#f8fafc'; g.fillRect(sg.q0, sg.p0, sg.q1 - sg.q0, sg.p1 - sg.p0);
        for (const v of ax.lines) { const x = Math.round(lrOf(sg, v, S)) + 0.5; g.strokeStyle = v === 0 ? '#94a3b8' : '#e2e8f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(x, sg.p0); g.lineTo(x, sg.p1); g.stroke(); }
        g.strokeStyle = '#cbd5e1'; g.strokeRect(Math.round(sg.q0) + 0.5, Math.round(sg.p0) + 0.5, Math.round(sg.q1 - sg.q0) - 1, Math.round(sg.p1 - sg.p0) - 1);
        g.fillStyle = '#334155'; g.fillText(sg.ch.slice(3), (sg.q0 + sg.q1) / 2, sg.p0 - 8);
        // the axis under the first column of each row
        if (i % Math.ceil(CHROMS.length / rows) === 0) {
          g.fillStyle = '#64748b'; g.font = '9px Inter, ui-sans-serif, system-ui, sans-serif';
          g.textAlign = 'left'; g.fillText(lrLabel(-S), sg.q0, sg.p1 + 8); g.textAlign = 'right'; g.fillText(lrLabel(S), sg.q1, sg.p1 + 8);
          g.textAlign = 'center'; g.font = '10px Inter, ui-sans-serif, system-ui, sans-serif';
        }
      });
    } else {
      segs.forEach((sg, i) => { if (i % 2) { g.fillStyle = '#f1f5f9'; g.fillRect(sg.p0, PT, sg.p1 - sg.p0, PH); } });
      for (const v of ax.lines) { const y = Math.round(lrOf(segs[0], v, S)) + 0.5; g.strokeStyle = v === 0 ? '#94a3b8' : '#e2e8f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(L, y); g.lineTo(w - R, y); g.stroke(); }
      g.fillStyle = '#64748b'; g.textAlign = 'right';
      for (const v of ax.labels) g.fillText(lrLabel(v), L - 4, Math.min(PB - 4, Math.max(PT + 4, lrOf(segs[0], v, S))));
      g.textAlign = 'center';
      if (chrom) {
        const sg = segs[0], step = [1e6, 2e6, 5e6, 1e7, 2e7, 5e7].find(st => sg.to / st <= 14) ?? 5e7;
        g.strokeStyle = '#cbd5e1';
        for (let t = 0; t <= sg.to; t += step) { const x = Math.round(posOf(sg, t)) + 0.5; g.beginPath(); g.moveTo(x, PB); g.lineTo(x, PB + 3); g.stroke(); g.fillText(`${t / 1e6} Mb`, x, PB + 10); }
      } else for (const sg of segs) g.fillText(sg.ch.slice(3), (sg.p0 + sg.p1) / 2, PB + 10);
      g.strokeStyle = '#cbd5e1'; g.strokeRect(L + 0.5, PT + 0.5, w - L - R - 1, PH - 1);
    }
    let n = 0, clipped = 0;
    if (pr) {
      const rws = H, cols = Math.ceil(w) + 1, grid = new Uint8Array(cols * rws), med: [number, number][][] = [];
      const clip: { sg: Seg; at: number; hi: boolean }[] = [], seen = new Set<string>();
      for (const sg of segs) {
        const q = pr.p.byChrom.get(sg.ch); if (!q) continue;
        const pts: [number, number][] = [], buf: number[] = [];
        let cell = -1;
        // the median line: per pixel along the position (whole genome), or over 15 probes (one chromosome)
        const flush = () => { if (buf.length) { const m = lrOf(sg, median(buf), S); pts.push(vert ? [m, cell + 0.5] : [cell + 0.5, m]); buf.length = 0; } };
        const lo = Math.min(sg.p0, sg.p1), hi = Math.max(sg.p0, sg.p1);
        for (let i = 0; i < q.lr.length; i++) {
          const v = q.lr[i]; if (Number.isNaN(v)) continue;
          const at = posOf(sg, (q.start[i] + q.end[i]) / 2); if (at < lo || at > hi) continue;
          const [x, y] = xyOf(sg, (q.start[i] + q.end[i]) / 2, v, S); n++;
          if (Math.abs(v) > S) {
            // beyond the range: one triangle per 3 px along the chromosome and side
            clipped++; const k = `${sg.ch}|${(at / 3) | 0}|${v > 0}`; if (!seen.has(k)) { seen.add(k); clip.push({ sg, at, hi: v > 0 }); }
          } else grid[(x | 0) * rws + (y | 0)] = v > THR ? 2 : v < -THR ? 3 : 1;
          const c = at | 0;
          if (chrom) { buf.push(v); if (buf.length === 15) { pts.push([x, lrOf(sg, median(buf.slice()), S)]); buf.shift(); } } else { if (c !== cell) { flush(); cell = c; } buf.push(v); }
        }
        if (!chrom) flush();
        med.push(pts);
      }
      const d = chrom ? 2 : 1.5;
      for (const [cls, col] of [[1, '#111827'], [2, '#dc2626'], [3, '#16a34a']] as const) {
        g.fillStyle = col; g.beginPath();
        for (let k = 0; k < grid.length; k++) if (grid[k] === cls) g.rect(Math.floor(k / rws) - d / 4, (k % rws) - d / 4, d, d);
        g.fill();
      }
      for (const up of [true, false]) {
        g.fillStyle = up ? '#dc2626' : '#16a34a'; g.beginPath();
        for (const c of clip) if (c.hi === up) {
          const e = lrOf(c.sg, up ? 1e9 : -1e9, S);
          if (c.sg.vert) clipMark(g, e, c.at, up ? 'r' : 'l'); else clipMark(g, c.at, e, up ? 'u' : 'd');
        }
        g.fill();
      }
      g.strokeStyle = '#1e3a8a'; g.lineWidth = vert ? 1 : 1.2;
      for (const pts of med) { g.beginPath(); pts.forEach(([x, y], j) => (j ? g.lineTo(x, y) : g.moveTo(x, y))); g.stroke(); }
    } else { g.fillStyle = '#64748b'; g.textAlign = 'center'; g.fillText('No probe file for this patient: its calls only', w / 2, vert ? H / 2 : PT + PH / 2 - 14); }
    el.setAttribute('data-ms', (performance.now() - t0).toFixed(0)); el.setAttribute('data-drawn', String(n)); el.setAttribute('data-clipped', String(clipped));
  }, [w, segs, H, pr, open, S]); // eslint-disable-line react-hooks/exhaustive-deps

  // the calls: file calls (loss red, gain blue) and probe calls (purple), cheap, redrawn with the state
  const barKey = `${calls.map(c => c.id).join()}|${events.length}|${s.shown}`;
  useEffect(() => {
    const el = top.current; if (!el || !w || !open) return;
    const dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    const out: Bar[] = [];
    // above the plot (horizontal), beside each column (vertical): file calls, then probe calls
    for (const [list, o, t] of [[calls, 3, 9], [events, 14, 7]] as const) for (const c of list) {
      const sg = segs.find(x => x.ch === c.chrom); if (!sg) continue;
      const a = posOf(sg, c.start), b = Math.max(a + 2, posOf(sg, c.end));
      const bar = sg.vert ? { x: sg.q1 + 2 + (o === 3 ? 0 : 6), y: a, w: o === 3 ? 5 : 4, h: b - a, c } : { x: a, y: o, w: b - a, h: t, c };
      g.fillStyle = c.fromProbes ? '#a21caf' : c.kind === 'loss' ? '#dc2626' : '#2563eb'; g.fillRect(bar.x, bar.y, bar.w, bar.h);
      if (c.id === s.shown) { g.strokeStyle = '#0f172a'; g.lineWidth = 1.5; g.strokeRect(bar.x - 1, bar.y - 1, bar.w + 2, bar.h + 2); }
      out.push(bar);
    }
    bars.current = out;
  }, [w, segs, H, open, barKey]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!cur) return null;
  const name = cur.x.name, male = s.sex[key] === 'male';
  if (!open) return (
    <button onClick={onToggle} data-cgh-genome-strip className="mx-4 mt-3 block w-[calc(100%-2rem)] text-left px-5 py-1.5 bg-white border border-slate-200 rounded-xl shadow-sm text-xs font-semibold text-slate-800 hover:border-accent-300"
      title="Show the whole genome of this patient (folds the call window)">▸ Whole genome — {name}</button>
  );
  /** the chromosome under the pointer: along its position axis (and, vertical, across its column and bars) */
  const segAt = (x: number, y: number) => segs.find(sg => sg.vert
    ? x >= sg.q0 - 3 && x <= sg.q1 + VBARS + 3 && y >= sg.p0 - 3 && y <= sg.p1 + 3
    : x >= sg.p0 - GAP / 2 && x <= sg.p1 + GAP / 2);
  const barAt = (x: number, y: number) => bars.current.find(b => x >= b.x - 1 && x <= b.x + b.w + 1 && y >= b.y - 1 && y <= b.y + b.h + 1);
  const onMove = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    const b = barAt(x, y);
    if (b) { setHover({ x, y, title: `${b.c.sample}${b.c.fromProbes ? ' · probe call' : ''}`, lines: callLines(s, b.c) }); return; }
    const sg = segAt(x, y), q = sg && pr?.p.byChrom.get(sg.ch);
    const across = sg ? (sg.vert ? x : y) : 0, along = sg ? (sg.vert ? y : x) : 0;
    if (!sg || !q || across < Math.min(sg.q0, sg.q1) - 2 || across > Math.max(sg.q0, sg.q1) + 2) { setHover(null); return; }
    // the nearest probe drawn within 6 px (the probes of ± 3 px of the pointer along the chromosome)
    const bpx = sg.to / (sg.p1 - sg.p0), bp = (along - sg.p0) * bpx;
    let best = -1, bd = 36;
    for (let i = lowerBound(q.start, bp - 3 * bpx); i < q.start.length && q.start[i] <= bp + 3 * bpx; i++) {
      const v = q.lr[i]; if (Number.isNaN(v)) continue;
      const [px, py] = xyOf(sg, (q.start[i] + q.end[i]) / 2, v, S), d2 = (px - x) ** 2 + (py - y) ** 2;
      if (d2 < bd) { bd = d2; best = i; }
    }
    if (best < 0) { setHover(null); return; }
    const [px, py] = xyOf(sg, (q.start[best] + q.end[best]) / 2, q.lr[best], S);
    setHover({ x: px, y: py, dot: [px, py], title: q.name[best] || 'probe', lines: [`${sg.ch}:${fmt(q.start[best])}-${fmt(q.end[best])} · log ratio ${q.lr[best].toFixed(3)}${Math.abs(q.lr[best]) > S ? ` (beyond ±${S})` : ''}`] });
  };
  const onClick = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    const b = barAt(x, y);
    if (chrom) { if (b) { cgh.setShown(b.c.id); onOpen(b.c, cur.f); } return; }
    const sg = segAt(x, y); if (sg) { setChrom(sg.ch); setHover(null); }
  };
  const step = (d: number) => setChrom(c => CHROMS[(CHROMS.indexOf(c!) + d + CHROMS.length) % CHROMS.length]);
  const lay = s.view.layout;
  return (
    <div className="mx-4 mt-3 bg-white border border-slate-200 rounded-xl shadow-sm pb-2" data-cgh-genome>
      <div className="px-5 pt-2 flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] text-slate-600">
        <button onClick={onToggle} className="text-xs font-semibold text-slate-900 hover:text-accent-700" title="Fold the whole genome">▾ Whole genome</button>
        {!s.sample && all.length > 1
          ? <select value={key} onChange={e => setPick(e.target.value)} aria-label="Patient" data-cgh-genome-pick className="h-[22px] rounded border border-slate-200 bg-white text-[12px] font-semibold outline-none cursor-pointer">{all.map(a => <option key={a.x.key} value={a.x.key}>{a.x.name}</option>)}</select>
          : <span className="text-sm font-bold text-slate-900">{name}</span>}
        {chrom && <span className="inline-flex items-center gap-1">
          <button className={BTN} onClick={() => setChrom(null)} title="Back to all chromosomes">whole genome</button>
          <button className={BTN} onClick={() => step(-1)} title="Previous chromosome" aria-label="Previous chromosome">←</button>
          <span className="font-semibold text-slate-900 text-xs">{chrom}</span>
          <button className={BTN} onClick={() => step(1)} title="Next chromosome" aria-label="Next chromosome">→</button></span>}
        {!chrom && <span className="inline-flex rounded border border-slate-200 overflow-hidden" role="group" aria-label="Layout" data-cgh-layout title="Chromosomes side by side (position across), or as columns (position down from pter, log ratio across); kept as a setting. A chromosome's own view stays horizontal.">
          {(['horizontal', 'vertical'] as const).map(v => <button key={v} onClick={() => cgh.setView({ layout: v })} aria-pressed={lay === v} className={`h-[22px] px-1.5 text-[11px] ${lay === v ? 'bg-accent-50 text-accent-700 font-semibold' : 'bg-white text-slate-600 hover:bg-slate-50'}`}>{v}</button>)}</span>}
        {!(chromPanel && pr) && <LrScale />}
        <span className="text-slate-400">· {chrom ? 'click a call bar: open it' : vert ? 'click a chromosome: its probes · position down from pter, log ratio across' : 'click a chromosome: its probes'} · dots black within ±0.2, red above, green below · bars: <span className="text-red-600">loss</span> / <span className="text-blue-600">gain</span> (calls listed), <span className="text-fuchsia-700">probe calls</span> · {chromPanel && pr ? 'every probe, as in the call\'s probe panel' : `clipped at ±${S} (▲▼ beyond)`}{male ? ' · male: X and Y near −1 (one copy, not normalised)' : ''}</span>
      </div>
      {chromPanel && pr ? <CghProbePanel key={`${key}|${chrom}`} at={{ sampleKey: key, chrom }} chromLen={chromLength(build, chrom) ?? 0} segs={cgh.probeCallsOf(key)?.segs.get(chrom)}
        calls={[...calls, ...events].filter(c => c.chrom === chrom)} onCall={c => { cgh.setShown(c.id); onOpen(c, cur.f); }} />
      : <div ref={box} className="relative mt-1 mx-1 select-none cursor-pointer" onMouseMove={onMove} onMouseLeave={() => setHover(null)} onClick={onClick}>
        <canvas ref={base} style={{ width: '100%', height: H, display: 'block' }} data-cgh-genome-probes data-layout={vert ? 'vertical' : 'horizontal'} />
        <canvas ref={top} style={{ width: '100%', height: H, display: 'block', position: 'absolute', inset: 0 }} />
        {hover?.dot && <div className="absolute pointer-events-none rounded-full border-2 border-slate-900" style={{ left: hover.dot[0] - 4, top: hover.dot[1] - 4, width: 8, height: 8 }} />}
        {hover && <CghTip x={hover.x} y={hover.y} w={w} title={hover.title} lines={hover.lines} />}
      </div>}
    </div>
  );
}
