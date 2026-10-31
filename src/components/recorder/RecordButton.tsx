/**
 * The Record button: records what the viewer draws, then offers the recording as an animated SVG and / or an HTML
 * page with play / pause, light enough for a presentation.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { startRecording, type Recorder, type Recording } from './recorder';
import { recordingToSvg, svgToHtml, autoZooms, timeMap, invertTime, MAX_IDLE_DEFAULT_MS, ZOOM_LEVEL_DEFAULT, type Zoom } from './exportRecording';

const fmtTime = (ms: number) => { const s = Math.max(0, ms) / 1000; return `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`; };
const fmtSize = (b: number) => (b >= 1e6 ? `${(b / 1e6).toFixed(1)} MB` : `${Math.max(1, Math.round(b / 1000))} kB`);
const stamp = () => { const d = new Date(), p = (n: number) => String(n).padStart(2, '0'); return `${d.getFullYear()}-${p(d.getMonth() + 1)}-${p(d.getDate())}-${p(d.getHours())}${p(d.getMinutes())}`; };
function download(text: string, type: string, name: string) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = document.createElement('a');
  a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
/** The gzip size, as a web server would send it (when the browser can compress). */
async function gzipSize(text: string): Promise<number | null> {
  try {
    const CS = (globalThis as any).CompressionStream;
    if (!CS) return null;
    const stream = new Blob([text]).stream().pipeThrough(new CS('gzip'));
    return (await new Response(stream).blob()).size;
  } catch { return null; }
}

export default function RecordButton({ names, className }: { names: string[]; className: string }) {
  const [rec, setRec] = useState<Recorder | null>(null);
  const [elapsed, setElapsed] = useState(0);
  const [done, setDone] = useState<Recording | null>(null);
  useEffect(() => {
    if (!rec) return;
    const id = setInterval(() => setElapsed(rec.elapsed()), 250);
    return () => clearInterval(id);
  }, [rec]);
  const start = () => { setDone(null); setElapsed(0); setRec(startRecording()); };
  const stop = () => { if (!rec) return; setDone(rec.stop()); setRec(null); };
  return (
    <>
      <button onClick={rec ? stop : start} className={className} data-recorder="skip"
        title={rec ? 'Stop the recording' : 'Record what the viewer draws (tracks, cartoons, SV browser) with the pointer and clicks, as vectors: exported as an animated SVG or an HTML page for a presentation'}>
        {rec ? <><span className="inline-block w-2 h-2 rounded-full bg-red-600 animate-pulse" /> Stop · {fmtTime(elapsed)}</> : <><span className="inline-block w-2 h-2 rounded-full bg-red-600" /> Record</>}
      </button>
      {rec && createPortal(
        <div data-recorder="skip" className="fixed bottom-4 right-4 z-[90] flex items-center gap-2 h-9 pl-3 pr-1 rounded-full bg-gray-900/90 text-white text-xs shadow-lg">
          <span className="inline-block w-2.5 h-2.5 rounded-full bg-red-500 animate-pulse" /> REC {fmtTime(elapsed)}
          <button onClick={stop} className="h-7 px-3 rounded-full bg-white text-gray-900 font-semibold hover:bg-gray-100">Stop</button>
        </div>, document.body)}
      {done && createPortal(<ExportDialog rec={done} names={names} onClose={() => setDone(null)} onAgain={start} />, document.body)}
    </>
  );
}

type Rect = [number, number, number, number];

function ExportDialog({ rec, names, onClose, onAgain }: { rec: Recording; names: string[]; onClose: () => void; onAgain: () => void }) {
  const [format, setFormat] = useState<'svg' | 'html' | 'both'>('html');
  const [shorten, setShorten] = useState(true);
  const [loop, setLoop] = useState(true);
  const [anon, setAnon] = useState(false);
  const [title, setTitle] = useState('Sashimi viewer recording');
  const uniq = useMemo(() => [...new Set(names.filter(Boolean))], [names]);
  const [aliases, setAliases] = useState<string[]>(() => uniq.map((_, i) => `Sample ${i + 1}`));
  // zooms: automatic on the clicks (each can be turned off), and drawn on the preview
  const [autoOn, setAutoOn] = useState(true);
  const [level, setLevel] = useState(ZOOM_LEVEL_DEFAULT);
  const [autoOff, setAutoOff] = useState<Set<number>>(new Set());
  const [drawn, setDrawn] = useState<Zoom[]>([]);
  const [drawing, setDrawing] = useState(false);
  const [ptr, setPtr] = useState({ show: true, smooth: true, scale: 1.4, clicks: true });
  const autos = useMemo(() => autoZooms(rec, level), [rec, level]);
  const zooms = useMemo(() => [...(autoOn ? autos.filter((_, i) => !autoOff.has(i)) : []), ...drawn], [autos, autoOn, autoOff, drawn]);
  const map = useMemo(() => timeMap(rec, shorten ? MAX_IDLE_DEFAULT_MS : 0), [rec, shorten]);
  const out = useMemo(() => recordingToSvg(rec, {
    maxIdleMs: shorten ? MAX_IDLE_DEFAULT_MS : 0, loop, title, zooms, pointer: ptr,
    replace: anon ? uniq.map((n, i) => ({ from: n, to: aliases[i] ?? `Sample ${i + 1}` })) : [],
  }), [rec, shorten, loop, anon, uniq, aliases, title, zooms, ptr]);
  const html = useMemo(() => svgToHtml(out, title), [out, title]);
  const fmtOut = out.hasHtml ? 'html' : format;
  const [gz, setGz] = useState<number | null>(null);
  useEffect(() => { let alive = true; gzipSize(fmtOut === 'svg' ? out.svg : html).then(v => { if (alive) setGz(v); }); return () => { alive = false; }; }, [out, html, fmtOut]);
  const save = () => {
    const base = `sashimi-recording-${stamp()}`;
    if (fmtOut !== 'html') download(out.svg, 'image/svg+xml', `${base}.svg`);
    if (fmtOut !== 'svg') download(html, 'text/html', `${base}.html`);
  };
  const bytes = fmtOut === 'svg' ? out.bytes : fmtOut === 'html' ? new Blob([html]).size : out.bytes + new Blob([html]).size;

  // ---- preview: the HTML page itself, in a frame; the time is kept when it is rebuilt ----
  const frame = useRef<HTMLIFrameElement>(null);
  const keep = useRef<{ t: number; paused: boolean }>({ t: 0, paused: false });
  const player = () => (frame.current?.contentWindow as any)?.recPlayer as { pause(): void; play(): void; seek(ms: number): void; time(): number; paused(): boolean } | undefined;
  const ready = useRef(false);
  const onFrameLoad = () => { const p = player(); if (!p) return; if (keep.current.t > 0) p.seek(keep.current.t); if (keep.current.paused) p.pause(); ready.current = true; };
  useEffect(() => { ready.current = false; }, [html]);
  useEffect(() => {
    const id = setInterval(() => { const p = player(); if (p && ready.current) keep.current = { t: p.time(), paused: p.paused() }; }, 200);
    return () => clearInterval(id);
  }, []);
  // ---- drawing a zoom on the preview ----
  const [drag, setDrag] = useState<{ x0: number; y0: number; x1: number; y1: number } | null>(null);
  useEffect(() => {
    if (!drawing) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') { e.stopPropagation(); setDrawing(false); setDrag(null); } };
    addEventListener('keydown', onKey, true);
    return () => removeEventListener('keydown', onKey, true);
  }, [drawing]);
  const startDraw = () => { const p = player(); if (p) { p.pause(); keep.current = { t: p.time(), paused: true }; } setDrawing(true); };
  const finishDraw = () => {
    const d = drag; setDrag(null); setDrawing(false);
    const doc = frame.current?.contentDocument, p = player();
    if (!d || !doc || !p || Math.abs(d.x1 - d.x0) < 12 || Math.abs(d.y1 - d.y0) < 12) return;
    const g = doc.querySelector('#cam > g') as SVGGElement | null, svg = doc.querySelector('#stage > svg') as SVGSVGElement | null;
    const m = g?.getScreenCTM()?.inverse();
    if (!g || !svg || !m) return;
    const pt = (x: number, y: number) => { const q = svg.createSVGPoint(); q.x = x; q.y = y; return q.matrixTransform(m); };
    const a = pt(Math.min(d.x0, d.x1), Math.min(d.y0, d.y1)), b = pt(Math.max(d.x0, d.x1), Math.max(d.y0, d.y1));
    const played = p.time();
    const start = invertTime(map, played, rec.duration), end = invertTime(map, played + 2500, rec.duration);
    setDrawn(z => [...z, { start, end: Math.max(end, start + 500), rect: [a.x, a.y, b.x - a.x, b.y - a.y] as Rect, source: 'drawn' }]);
  };
  const pos = (e: React.MouseEvent) => { const r = (e.currentTarget as HTMLElement).getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const ts = (ms: number) => `${(map(ms) / 1000).toFixed(1)} s`;

  return (
    <div data-recorder="skip" className="fixed inset-0 z-[90] bg-slate-900/70 flex items-start justify-center overflow-y-auto p-4" onMouseDown={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-[1200px] text-gray-900 text-sm" onMouseDown={e => e.stopPropagation()}>
        <div className="flex items-start justify-between px-5 pt-4 pb-2">
          <div>
            <div className="text-base font-bold">Recording <span className="font-normal text-gray-500">· {fmtTime(out.duration)} · {out.frames} frames · {fmtSize(bytes)}{gz ? ` (${fmtSize(gz)} gzipped)` : ''}</span></div>
            <div className="text-xs text-gray-500">vectors only: sharp at any zoom, text kept as text · the drawings and the dialogs are recorded, not the toolbar</div>
          </div>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-2xl leading-none px-1" title="Close">×</button>
        </div>
        <div className="px-5">
          <div className="relative rounded-lg border border-gray-200 overflow-hidden" style={{ height: '58vh' }}>
            <iframe ref={frame} title="preview" srcDoc={html} onLoad={onFrameLoad} className="w-full h-full block" />
            {drawing && (
              <div className="absolute inset-0 cursor-crosshair bg-accent-500/5" onMouseDown={e => { const p = pos(e); setDrag({ x0: p.x, y0: p.y, x1: p.x, y1: p.y }); }}
                onMouseMove={e => { if (drag) { const p = pos(e); setDrag({ ...drag, x1: p.x, y1: p.y }); } }} onMouseUp={finishDraw}>
                {drag && <div className="absolute border-2 border-accent-600 bg-accent-500/10 rounded" style={{ left: Math.min(drag.x0, drag.x1), top: Math.min(drag.y0, drag.y1), width: Math.abs(drag.x1 - drag.x0), height: Math.abs(drag.y1 - drag.y0) }} />}
                {!drag && <div className="absolute top-2 left-1/2 -translate-x-1/2 px-3 py-1 rounded-full bg-gray-900/80 text-white text-xs">drag over the area to zoom on, from this moment (2.5 s) · Esc: cancel</div>}
              </div>
            )}
          </div>
        </div>
        <div className="px-5 py-3 grid grid-cols-[1.15fr_1fr_1fr] gap-5 text-xs">
          <div>
            <div className="font-semibold text-gray-700 mb-1">Zooms</div>
            <div className="flex items-center gap-3 mb-1">
              <label className="inline-flex items-center gap-1.5"><input type="checkbox" checked={autoOn} onChange={e => setAutoOn(e.target.checked)} className="accent-accent-600" />on the clicks</label>
              <select value={level} onChange={e => { setLevel(Number(e.target.value)); setAutoOff(new Set()); }} className="border border-gray-300 rounded px-1 py-0.5">
                {[1.5, 2, 3, 4].map(l => <option key={l} value={l}>{l}×</option>)}
              </select>
              <button onClick={startDraw} className="ml-auto h-7 px-2.5 rounded-md border border-accent-300 text-accent-700 hover:bg-accent-50">＋ draw a zoom</button>
            </div>
            <div className="space-y-0.5 max-h-28 overflow-y-auto">
              {autoOn && autos.map((z, i) => (
                <label key={`a${i}`} className="flex items-center gap-1.5 text-gray-700"><input type="checkbox" checked={!autoOff.has(i)} onChange={() => setAutoOff(s => { const n = new Set(s); if (n.has(i)) n.delete(i); else n.add(i); return n; })} className="accent-accent-600" />
                  click zoom · {ts(z.start)} → {ts(z.end)}</label>
              ))}
              {drawn.map((z, i) => (
                <div key={`d${i}`} className="flex items-center gap-1.5 text-gray-700"><span className="w-3.5 text-center text-accent-700">▣</span>drawn zoom · {ts(z.start)} →
                  <select value={Math.round((map(z.end) - map(z.start)) / 500) * 500} onChange={e => { const d = Number(e.target.value); setDrawn(list => list.map((x, k) => (k === i ? { ...x, end: invertTime(map, map(x.start) + d, rec.duration) } : x))); }} className="border border-gray-300 rounded px-1">
                    {[1500, 2500, 4000, 6000].map(d => <option key={d} value={d}>{d / 1000} s</option>)}
                  </select>
                  <button onClick={() => setDrawn(list => list.filter((_, k) => k !== i))} className="text-gray-400 hover:text-red-600" title="Remove">×</button></div>
              ))}
              {!zooms.length && <div className="text-gray-400">no zoom</div>}
            </div>
          </div>
          <div className="space-y-1.5">
            <div className="font-semibold text-gray-700">Pointer</div>
            <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
              <label className="inline-flex items-center gap-1.5"><input type="checkbox" checked={ptr.show} onChange={e => setPtr({ ...ptr, show: e.target.checked })} className="accent-accent-600" />shown</label>
              <label className="inline-flex items-center gap-1.5"><input type="checkbox" checked={ptr.smooth} onChange={e => setPtr({ ...ptr, smooth: e.target.checked })} className="accent-accent-600" />smoothed</label>
              <label className="inline-flex items-center gap-1.5"><input type="checkbox" checked={ptr.clicks} onChange={e => setPtr({ ...ptr, clicks: e.target.checked })} className="accent-accent-600" />click rings</label>
              <select value={ptr.scale} onChange={e => setPtr({ ...ptr, scale: Number(e.target.value) })} className="border border-gray-300 rounded px-1 py-0.5" title="Pointer size">
                {[1, 1.4, 2, 3].map(s => <option key={s} value={s}>size {s}×</option>)}
              </select>
            </div>
            <div className="font-semibold text-gray-700 pt-1">Playback</div>
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={shorten} onChange={e => setShorten(e.target.checked)} className="accent-accent-600" />shorten the pauses to {MAX_IDLE_DEFAULT_MS / 1000} s</label>
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={loop} onChange={e => setLoop(e.target.checked)} className="accent-accent-600" />loop</label>
          </div>
          <div className="space-y-1.5">
            <div className="font-semibold text-gray-700">File</div>
            <label className="block"><span className="text-gray-600 mr-2">Title</span><input value={title} onChange={e => setTitle(e.target.value)} className="border border-gray-300 rounded px-2 py-0.5 w-56" /></label>
            <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
              {(['html', 'svg', 'both'] as const).map(f => (
                <label key={f} className={`inline-flex items-center gap-1 ${out.hasHtml && f !== 'html' ? 'opacity-40' : ''}`}><input type="radio" disabled={out.hasHtml && f !== 'html'} checked={fmtOut === f} onChange={() => setFormat(f)} className="accent-accent-600" />
                  {f === 'html' ? 'HTML page' : f === 'svg' ? 'SVG' : 'both'}</label>
              ))}
            </div>
            {out.hasHtml && <div className="text-[11px] text-gray-500">a dialog is recorded with its text and tables: HTML page only</div>}
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={anon} onChange={e => setAnon(e.target.checked)} className="accent-accent-600" />anonymise the sample names</label>
            {anon && (
              <div className="space-y-1 max-h-24 overflow-y-auto">
                {uniq.map((n, i) => (
                  <div key={n} className="flex items-center gap-2"><span className="font-mono text-gray-500 w-32 truncate" title={n}>{n}</span>→
                    <input value={aliases[i] ?? ''} onChange={e => setAliases(a => { const b = [...a]; b[i] = e.target.value; return b; })} className="border border-gray-300 rounded px-2 py-0.5 w-28" /></div>
                ))}
                {!uniq.length && <div className="text-gray-400">no sample open</div>}
              </div>
            )}
          </div>
        </div>
        <div className="flex items-center justify-end gap-2 px-5 pb-4">
          <button onClick={() => { onClose(); onAgain(); }} className="h-8 px-3 rounded-md border border-gray-300 text-gray-700 hover:bg-gray-50">Record again</button>
          <button onClick={onClose} className="h-8 px-3 rounded-md border border-gray-300 text-gray-700 hover:bg-gray-50">Discard</button>
          <button onClick={save} className="h-8 px-4 rounded-md bg-accent-600 text-white font-semibold hover:bg-accent-700">Download</button>
        </div>
      </div>
    </div>
  );
}
