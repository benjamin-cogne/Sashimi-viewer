/**
 * The *Whole genome* panel of the array-CGH page: one patient's probes (CGH_Log_Ratio.txt) over chromosomes 1-22, X, Y
 * side by side, with its calls as bars above; a chromosome clicked opens it in the same panel, a call bar clicked there
 * opens the call. It stands where the call window does: expanded while no call is open, a strip once one is.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { cgh, useCgh, classOf, passes, probeCnv, type CghState } from './cghStore';
import type { CghCnv, CghFile } from './cghIntervals';
import { lowerBound } from './cghProbes';
import { chromLength } from '../../standalone/assembly';
import { estCopies } from './cghMosaic';
import { viewOf } from './cghAnnotate';
import { LEVEL_STYLE } from './cghScore';

const CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY'];
const L = 34, R = 8, PT = 24, PH = 170, PB = PT + PH, H = PB + 18, CLIP = 2, THR = 0.2, GAP = 2;
const BTN = 'h-[22px] min-w-[24px] px-1.5 rounded border border-slate-200 bg-white text-[11px] text-slate-700 hover:border-accent-300 hover:bg-accent-50';
const yOf = (v: number) => PT + ((CLIP - Math.max(-CLIP, Math.min(CLIP, v))) / (2 * CLIP)) * PH;
const fmt = (n: number) => n.toLocaleString('en-US');
type Seg = { ch: string; x0: number; x1: number; to: number };
type Bar = { x: number; y: number; w: number; h: number; c: CghCnv };
const xOf = (sg: Seg, bp: number) => sg.x0 + (bp / sg.to) * (sg.x1 - sg.x0);

/** a call's characteristics, for the hover tips (here and in the probe panel) */
export function callLines(s: CghState, c: CghCnv): string[] {
  const size = c.sizeKb >= 1000 ? `${(c.sizeKb / 1000).toFixed(2)} Mb` : `${c.sizeKb.toFixed(1)} kb`;
  const est = c.copies == null ? estCopies(c, s.sex[c.sampleKey]) : null, lv = c.fromProbes ? null : viewOf(s, c).score?.level;
  return [
    `${c.chrom}:${fmt(c.start)}-${fmt(c.end)} · ${size} · ${c.kind}`,
    `${c.probes ?? '?'} probes${c.logRatio != null ? ` · mean log ratio ${c.logRatio.toFixed(2)}` : ''}`,
    c.copies != null ? `copies x${c.copies}${c.copiesHi != null ? `~${c.copiesHi}` : ''} (file)` : est ? `estimated ${est.text}` : 'copies not given',
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
  useEffect(() => {
    const el = box.current; if (!el) return;
    const ro = new ResizeObserver(() => setW(el.clientWidth)); ro.observe(el); setW(el.clientWidth);
    return () => ro.disconnect();
  }, [open, !!cur]); // eslint-disable-line react-hooks/exhaustive-deps
  // read on first view only (some tens of ms for 61k probes)
  const pr = useMemo(() => {
    if (!open || !key || !cgh.hasProbes(key)) return null;
    const t = performance.now(), p = cgh.probes(key);
    return p ? { p, ms: performance.now() - t } : null;
  }, [open, key]);
  const segs = useMemo<Seg[]>(() => {
    if (!w) return [];
    const lens = CHROMS.map(ch => chromLength(build, ch) ?? 1);
    if (chrom) return [{ ch: chrom, x0: L, x1: w - R, to: lens[CHROMS.indexOf(chrom)] }];
    const tot = lens.reduce((a, b) => a + b, 0), k = (w - L - R - GAP * (CHROMS.length - 1)) / tot;
    let x = L;
    return CHROMS.map((ch, i) => { const sg = { ch, x0: x, x1: x + lens[i] * k, to: lens[i] }; x = sg.x1 + GAP; return sg; });
  }, [w, chrom, build]);
  const calls = cur ? cur.x.cnvs.filter(c => passes(s, c)) : [];
  const events = cur ? (cgh.probeCallsOf(key)?.events ?? []).map(e => probeCnv(e, key, cur.x.name)) : [];

  // the probes: one mark per pixel cell (the same picture, bounded by the pixels), a median line; redrawn only on resize / view
  useEffect(() => {
    const el = base.current; if (!el || !w || !open) return;
    const t0 = performance.now(), dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    g.font = '10px Inter, ui-sans-serif, system-ui, sans-serif'; g.textBaseline = 'middle';
    segs.forEach((sg, i) => { if (i % 2) { g.fillStyle = '#f1f5f9'; g.fillRect(sg.x0, PT, sg.x1 - sg.x0, PH); } });
    for (const v of [1, 0.25, 0, -0.25, -1]) { g.strokeStyle = v === 0 ? '#94a3b8' : '#e2e8f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(L, Math.round(yOf(v)) + 0.5); g.lineTo(w - R, Math.round(yOf(v)) + 0.5); g.stroke(); }
    g.fillStyle = '#64748b'; g.textAlign = 'right';
    for (const v of [2, 1, 0, -1, -2]) g.fillText(v > 0 ? `+${v}` : String(v), L - 4, Math.min(PB - 4, Math.max(PT + 4, yOf(v))));
    g.textAlign = 'center';
    if (chrom) {
      const sg = segs[0], step = [1e6, 2e6, 5e6, 1e7, 2e7, 5e7].find(st => sg.to / st <= 14) ?? 5e7;
      g.strokeStyle = '#cbd5e1';
      for (let t = 0; t <= sg.to; t += step) { const x = Math.round(xOf(sg, t)) + 0.5; g.beginPath(); g.moveTo(x, PB); g.lineTo(x, PB + 3); g.stroke(); g.fillText(`${t / 1e6} Mb`, x, PB + 10); }
    } else for (const sg of segs) g.fillText(sg.ch.slice(3), (sg.x0 + sg.x1) / 2, PB + 10);
    g.strokeStyle = '#cbd5e1'; g.strokeRect(L + 0.5, PT + 0.5, w - L - R - 1, PH - 1);
    let n = 0;
    if (pr) {
      const rows = H, cols = Math.ceil(w) + 1, grid = new Uint8Array(cols * rows), med: { sg: Seg; pts: [number, number][] }[] = [];
      for (const sg of segs) {
        const q = pr.p.byChrom.get(sg.ch); if (!q) continue;
        const pts: [number, number][] = [], buf: number[] = [];
        let col = -1;
        const flush = () => { if (buf.length) { pts.push([col + 0.5, yOf(median(buf))]); buf.length = 0; } };
        for (let i = 0; i < q.lr.length; i++) {
          const v = q.lr[i]; if (Number.isNaN(v)) continue;
          const x = xOf(sg, (q.start[i] + q.end[i]) / 2); if (x < sg.x0 || x > sg.x1) continue;
          const cx = x | 0; n++;
          grid[cx * rows + (yOf(v) | 0)] = v > THR ? 2 : v < -THR ? 3 : 1;
          // the median line: per pixel column over the genome; over 15 probes for one chromosome (few probes a column)
          if (chrom) { buf.push(v); if (buf.length === 15) { pts.push([x, yOf(median(buf.slice()))]); buf.shift(); } } else { if (cx !== col) { flush(); col = cx; } buf.push(v); }
        }
        if (!chrom) flush();
        med.push({ sg, pts });
      }
      const d = chrom ? 2 : 1.5;
      for (const [cls, col] of [[1, '#111827'], [2, '#dc2626'], [3, '#16a34a']] as const) {
        g.fillStyle = col; g.beginPath();
        for (let k = 0; k < grid.length; k++) if (grid[k] === cls) g.rect(Math.floor(k / rows) - d / 4, (k % rows) - d / 4, d, d);
        g.fill();
      }
      g.strokeStyle = '#1e3a8a'; g.lineWidth = 1.2;
      for (const { pts } of med) { g.beginPath(); pts.forEach(([x, y], j) => (j ? g.lineTo(x, y) : g.moveTo(x, y))); g.stroke(); }
    } else { g.fillStyle = '#64748b'; g.fillText('No probe file for this patient: its calls only', (L + w - R) / 2, PT + PH / 2 - 14); }
    el.setAttribute('data-ms', (performance.now() - t0).toFixed(0)); el.setAttribute('data-drawn', String(n));
  }, [w, segs, pr, open]); // eslint-disable-line react-hooks/exhaustive-deps

  // the calls: file calls (loss red, gain blue) and probe calls (purple), cheap, redrawn with the state
  const barKey = `${calls.map(c => c.id).join()}|${events.length}|${s.shown}`;
  useEffect(() => {
    const el = top.current; if (!el || !w || !open) return;
    const dpr = window.devicePixelRatio || 1;
    el.width = Math.round(w * dpr); el.height = Math.round(H * dpr);
    const g = el.getContext('2d')!; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, H);
    const out: Bar[] = [];
    for (const [list, y, h] of [[calls, 3, 9], [events, 14, 7]] as const) for (const c of list) {
      const sg = segs.find(x => x.ch === c.chrom); if (!sg) continue;
      const a = xOf(sg, c.start), b = Math.max(a + 2, xOf(sg, c.end)), bar = { x: a, y, w: b - a, h, c };
      g.fillStyle = c.fromProbes ? '#a21caf' : c.kind === 'loss' ? '#dc2626' : '#2563eb'; g.fillRect(a, y, b - a, h);
      if (c.id === s.shown) { g.strokeStyle = '#0f172a'; g.lineWidth = 1.5; g.strokeRect(a - 1, y - 1, b - a + 2, h + 2); }
      out.push(bar);
    }
    bars.current = out;
  }, [w, segs, open, barKey]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!cur) return null;
  const name = cur.x.name, male = s.sex[key] === 'male';
  if (!open) return (
    <button onClick={onToggle} data-cgh-genome-strip className="mx-4 mt-3 block w-[calc(100%-2rem)] text-left px-5 py-1.5 bg-white border border-slate-200 rounded-xl shadow-sm text-xs font-semibold text-slate-800 hover:border-accent-300"
      title="Show the whole genome of this patient (folds the call window)">▸ Whole genome — {name}</button>
  );
  const segAt = (x: number) => segs.find(sg => x >= sg.x0 - GAP / 2 && x <= sg.x1 + GAP / 2);
  const barAt = (x: number, y: number) => bars.current.find(b => x >= b.x - 1 && x <= b.x + b.w + 1 && y >= b.y - 1 && y <= b.y + b.h + 1);
  const onMove = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    const b = barAt(x, y);
    if (b) { setHover({ x, y, title: `${b.c.sample}${b.c.fromProbes ? ' · probe call' : ''}`, lines: callLines(s, b.c) }); return; }
    const sg = segAt(x), q = sg && pr?.p.byChrom.get(sg.ch);
    if (!sg || !q || y < PT || y > PB) { setHover(null); return; }
    // the nearest probe drawn within 6 px (the probes of ± 3 px of the pointer)
    const bpx = sg.to / (sg.x1 - sg.x0), bp = (x - sg.x0) * bpx;
    let best = -1, bd = 36;
    for (let i = lowerBound(q.start, bp - 3 * bpx); i < q.start.length && q.start[i] <= bp + 3 * bpx; i++) {
      const v = q.lr[i]; if (Number.isNaN(v)) continue;
      const dx = xOf(sg, (q.start[i] + q.end[i]) / 2) - x, dy = yOf(v) - y, d2 = dx * dx + dy * dy;
      if (d2 < bd) { bd = d2; best = i; }
    }
    if (best < 0) { setHover(null); return; }
    const px = xOf(sg, (q.start[best] + q.end[best]) / 2), py = yOf(q.lr[best]);
    setHover({ x: px, y: py, dot: [px, py], title: q.name[best] || 'probe', lines: [`${sg.ch}:${fmt(q.start[best])}-${fmt(q.end[best])} · log ratio ${q.lr[best].toFixed(3)}`] });
  };
  const onClick = (e: React.MouseEvent) => {
    const rc = e.currentTarget.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top;
    const b = barAt(x, y);
    if (chrom) { if (b) { cgh.setShown(b.c.id); onOpen(b.c, cur.f); } return; }
    const sg = segAt(x); if (sg) { setChrom(sg.ch); setHover(null); }
  };
  const step = (d: number) => setChrom(c => CHROMS[(CHROMS.indexOf(c!) + d + CHROMS.length) % CHROMS.length]);
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
        <span className="text-slate-400">· {chrom ? 'click a call bar: open it' : 'click a chromosome: its probes'} · dots black within ±0.2, red above, green below · bars: <span className="text-red-600">loss</span> / <span className="text-blue-600">gain</span> (calls listed), <span className="text-fuchsia-700">probe calls</span> · clipped at ±2{male ? ' · male: X and Y near −1 (one copy, not normalised)' : ''}</span>
      </div>
      <div ref={box} className="relative mt-1 mx-1 select-none cursor-pointer" onMouseMove={onMove} onMouseLeave={() => setHover(null)} onClick={onClick}>
        <canvas ref={base} style={{ width: '100%', height: H, display: 'block' }} data-cgh-genome-probes />
        <canvas ref={top} style={{ width: '100%', height: H, display: 'block', position: 'absolute', inset: 0 }} />
        {hover?.dot && <div className="absolute pointer-events-none rounded-full border-2 border-slate-900" style={{ left: hover.dot[0] - 4, top: hover.dot[1] - 4, width: 8, height: 8 }} />}
        {hover && <CghTip x={hover.x} y={hover.y} w={w} title={hover.title} lines={hover.lines} />}
      </div>
    </div>
  );
}
