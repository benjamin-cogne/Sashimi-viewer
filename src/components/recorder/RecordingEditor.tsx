/**
 * The editor that opens when a recording stops: a live preview (the exported HTML page itself) and tabs to set the
 * zooms, the captions, the chapters, the pointer and the file, then download.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import type { Recording } from './recorder';
import { recordingToSvg, svgToHtml, autoZooms, timeMap, invertTime, MAX_IDLE_DEFAULT_MS, ZOOM_LEVEL_DEFAULT, type Zoom } from './exportRecording';
import { suggestCaptions, CAPTION_DEFAULT_MS, type Caption, type CaptionKind, type Chapter } from './captions';

type Rect = [number, number, number, number];
type Tab = 'zooms' | 'captions' | 'chapters' | 'pointer' | 'file';
/** What a drag (or a click) on the preview makes. */
type Pick = null | 'zoom' | 'highlight' | 'callout';

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
const KIND_LABEL: Record<CaptionKind, string> = { title: 'Title card', lower: 'Lower third', callout: 'Callout', highlight: 'Highlight' };
const KIND_HELP: Record<CaptionKind, string> = {
  title: 'a centred card over the whole view: a title and a subtitle',
  lower: 'a bar at the bottom left: one to three lines of narration',
  callout: 'a label pointing at a place of the drawing (click it on the preview); it follows the zooms',
  highlight: 'a ring around an area (drag it on the preview), optionally dimming the rest',
};

export default function RecordingEditor({ rec, names, onClose, onAgain }: { rec: Recording; names: string[]; onClose: () => void; onAgain: () => void }) {
  const [tab, setTab] = useState<Tab>('zooms');
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
  const [ptr, setPtr] = useState({ show: true, smooth: true, scale: 1.4, clicks: true });
  const [captions, setCaptions] = useState<Caption[]>([]);
  const [chapters, setChapters] = useState<Chapter[]>([]);
  const suggestions = useMemo(() => suggestCaptions(rec), [rec]);
  const autos = useMemo(() => autoZooms(rec, level), [rec, level]);
  const zooms = useMemo(() => [...(autoOn ? autos.filter((_, i) => !autoOff.has(i)) : []), ...drawn], [autos, autoOn, autoOff, drawn]);
  const map = useMemo(() => timeMap(rec, shorten ? MAX_IDLE_DEFAULT_MS : 0), [rec, shorten]);
  const out = useMemo(() => recordingToSvg(rec, {
    maxIdleMs: shorten ? MAX_IDLE_DEFAULT_MS : 0, loop, title, zooms, pointer: ptr, captions, chapters,
    replace: anon ? uniq.map((n, i) => ({ from: n, to: aliases[i] ?? `Sample ${i + 1}` })) : [],
  }), [rec, shorten, loop, anon, uniq, aliases, title, zooms, ptr, captions, chapters]);
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

  // ---- preview: the HTML page itself, in a frame; its time is kept when it is rebuilt ----
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
  /** The preview's moment, as recorded time (the preview pauses). */
  const nowRec = (): number => { const p = player(); if (p) { p.pause(); keep.current = { t: p.time(), paused: true }; } return invertTime(map, keep.current.t, rec.duration); };
  const span = (startRec: number, ms: number) => ({ start: startRec, end: invertTime(map, map(startRec) + ms, rec.duration) });
  const ts = (ms: number) => `${(map(ms) / 1000).toFixed(1)} s`;

  // ---- picking on the preview: a rectangle (zoom, highlight) or a point (callout) ----
  const [pick, setPick] = useState<{ what: Pick; onDone?: (r: Rect | [number, number]) => void }>({ what: null });
  const [drag, setDrag] = useState<{ x0: number; y0: number; x1: number; y1: number } | null>(null);
  useEffect(() => {
    if (!pick.what) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') { e.stopPropagation(); setPick({ what: null }); setDrag(null); } };
    addEventListener('keydown', onKey, true);
    return () => removeEventListener('keydown', onKey, true);
  }, [pick.what]);
  const startPick = (what: Pick, onDone: (r: Rect | [number, number]) => void) => { const p = player(); if (p) { p.pause(); keep.current = { t: p.time(), paused: true }; } setPick({ what, onDone }); };
  const toPage = (x: number, y: number): [number, number] | null => {
    const doc = frame.current?.contentDocument;
    const g = doc?.querySelector('#cam > g') as SVGGElement | null, svg = doc?.querySelector('#stage > svg') as SVGSVGElement | null;
    const m = g?.getScreenCTM()?.inverse();
    if (!svg || !m) return null;
    const q = svg.createSVGPoint(); q.x = x; q.y = y;
    const r = q.matrixTransform(m);
    return [Math.round(r.x), Math.round(r.y)];
  };
  const finishPick = () => {
    const d = drag, done = pick.onDone, what = pick.what;
    setDrag(null); setPick({ what: null });
    if (!d || !done) return;
    if (what === 'callout') { const p = toPage(d.x1, d.y1); if (p) done(p); return; }
    if (Math.abs(d.x1 - d.x0) < 12 || Math.abs(d.y1 - d.y0) < 12) return;
    const a = toPage(Math.min(d.x0, d.x1), Math.min(d.y0, d.y1)), b = toPage(Math.max(d.x0, d.x1), Math.max(d.y0, d.y1));
    if (a && b) done([a[0], a[1], b[0] - a[0], b[1] - a[1]]);
  };
  const pos = (e: React.MouseEvent) => { const r = (e.currentTarget as HTMLElement).getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };

  // ---- captions ----
  const addCaption = (kind: CaptionKind) => {
    const t = nowRec();
    const base: Caption = { kind, ...span(t, CAPTION_DEFAULT_MS[kind]), text: kind === 'title' ? title : kind === 'highlight' ? '' : 'Text' };
    if (kind === 'callout') startPick('callout', p => setCaptions(cs => [...cs, { ...base, at: p as [number, number] }]));
    else if (kind === 'highlight') startPick('highlight', r => setCaptions(cs => [...cs, { ...base, rect: r as Rect, dim: true }]));
    else setCaptions(cs => [...cs, base]);
  };
  const editCaption = (i: number, patch: Partial<Caption>) => setCaptions(cs => cs.map((c, k) => (k === i ? { ...c, ...patch } : c)));
  const durOf = (c: { start: number; end: number }) => Math.round((map(c.end) - map(c.start)) / 500) * 500;
  const DUR = [1500, 2500, 3000, 4000, 6000, 8000];
  const durSelect = (c: { start: number; end: number }, set: (end: number) => void) => (
    <select value={DUR.includes(durOf(c)) ? durOf(c) : 3000} onChange={e => set(invertTime(map, map(c.start) + Number(e.target.value), rec.duration))} className="border border-gray-300 rounded px-1 py-0.5">
      {DUR.map(d => <option key={d} value={d}>{d / 1000} s</option>)}
    </select>
  );
  const startInput = (c: { start: number; end: number }, set: (start: number, end: number) => void) => (
    <input type="number" step={0.1} min={0} value={(map(c.start) / 1000).toFixed(1)} className="border border-gray-300 rounded px-1 py-0.5 w-16" title="Start (s)"
      onChange={e => { const ps = Math.max(0, Number(e.target.value) * 1000), d = map(c.end) - map(c.start); set(invertTime(map, ps, rec.duration), invertTime(map, ps + d, rec.duration)); }} />
  );

  const tabBtn = (t: Tab, label: string, n?: number) => (
    <button onClick={() => setTab(t)} className={`h-8 px-3 rounded-md text-xs font-semibold ${tab === t ? 'bg-accent-600 text-white' : 'text-gray-600 hover:bg-gray-100'}`}>{label}{n ? ` · ${n}` : ''}</button>
  );
  const check = (on: boolean, set: (v: boolean) => void, label: string, title?: string) => (
    <label className="inline-flex items-center gap-1.5" title={title}><input type="checkbox" checked={on} onChange={e => set(e.target.checked)} className="accent-accent-600" />{label}</label>
  );

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
          <div className="relative rounded-lg border border-gray-200 overflow-hidden" style={{ height: '56vh' }}>
            <iframe ref={frame} title="preview" srcDoc={html} onLoad={onFrameLoad} className="w-full h-full block" />
            {pick.what && (
              <div className="absolute inset-0 cursor-crosshair bg-accent-500/5" onMouseDown={e => { const p = pos(e); setDrag({ x0: p.x, y0: p.y, x1: p.x, y1: p.y }); }}
                onMouseMove={e => { if (drag) { const p = pos(e); setDrag({ ...drag, x1: p.x, y1: p.y }); } }} onMouseUp={finishPick}>
                {drag && pick.what !== 'callout' && <div className="absolute border-2 border-accent-600 bg-accent-500/10 rounded" style={{ left: Math.min(drag.x0, drag.x1), top: Math.min(drag.y0, drag.y1), width: Math.abs(drag.x1 - drag.x0), height: Math.abs(drag.y1 - drag.y0) }} />}
                {!drag && <div className="absolute top-2 left-1/2 -translate-x-1/2 px-3 py-1 rounded-full bg-gray-900/80 text-white text-xs">
                  {pick.what === 'callout' ? 'click the place to point at' : pick.what === 'highlight' ? 'drag around the area to highlight' : 'drag over the area to zoom on'}, from this moment · Esc: cancel</div>}
              </div>
            )}
          </div>
        </div>
        <div className="px-5 pt-3 flex items-center gap-1 border-b border-gray-100 pb-2">
          {tabBtn('zooms', 'Zooms', zooms.length)}{tabBtn('captions', 'Captions', captions.length)}{tabBtn('chapters', 'Chapters', chapters.length)}{tabBtn('pointer', 'Pointer & playback')}{tabBtn('file', 'File')}
          <span className="ml-auto text-[11px] text-gray-500">new items start at the preview's moment</span>
        </div>
        <div className="px-5 py-3 text-xs min-h-[150px]">
          {tab === 'zooms' && (
            <div>
              <div className="flex items-center gap-3 mb-1.5">
                {check(autoOn, setAutoOn, 'zoom on the clicks')}
                <select value={level} onChange={e => { setLevel(Number(e.target.value)); setAutoOff(new Set()); }} className="border border-gray-300 rounded px-1 py-0.5">
                  {[1.5, 2, 3, 4].map(l => <option key={l} value={l}>{l}×</option>)}
                </select>
                <button onClick={() => startPick('zoom', r => { const t = nowRec(); setDrawn(z => [...z, { ...span(t, 2500), rect: r as Rect, source: 'drawn' }]); })} className="ml-auto h-7 px-2.5 rounded-md border border-accent-300 text-accent-700 hover:bg-accent-50">＋ draw a zoom</button>
              </div>
              <div className="space-y-1 max-h-32 overflow-y-auto">
                {autoOn && autos.map((z, i) => (
                  <label key={`a${i}`} className="flex items-center gap-1.5 text-gray-700"><input type="checkbox" checked={!autoOff.has(i)} onChange={() => setAutoOff(s => { const n = new Set(s); if (n.has(i)) n.delete(i); else n.add(i); return n; })} className="accent-accent-600" />
                    click zoom · {ts(z.start)} → {ts(z.end)}</label>
                ))}
                {drawn.map((z, i) => (
                  <div key={`d${i}`} className="flex items-center gap-1.5 text-gray-700"><span className="w-3.5 text-center text-accent-700">▣</span>drawn zoom · from
                    {startInput(z, (start, end) => setDrawn(l => l.map((x, k) => (k === i ? { ...x, start, end } : x))))} s for
                    {durSelect(z, end => setDrawn(l => l.map((x, k) => (k === i ? { ...x, end } : x))))}
                    <button onClick={() => setDrawn(list => list.filter((_, k) => k !== i))} className="text-gray-400 hover:text-red-600 text-base leading-none" title="Remove">×</button></div>
                ))}
                {!zooms.length && <div className="text-gray-400">no zoom</div>}
              </div>
            </div>
          )}
          {tab === 'captions' && (
            <div>
              <div className="flex flex-wrap items-center gap-1.5 mb-2">
                {(['title', 'lower', 'callout', 'highlight'] as CaptionKind[]).map(k => (
                  <button key={k} onClick={() => addCaption(k)} title={KIND_HELP[k]} className="h-7 px-2.5 rounded-md border border-accent-300 text-accent-700 hover:bg-accent-50">＋ {KIND_LABEL[k].toLowerCase()}</button>
                ))}
                {suggestions.filter(sg => !captions.some(c => c.text === sg.text)).slice(0, 3).map((sg, i) => (
                  <button key={i} onClick={() => setCaptions(cs => [...cs, sg])} className="h-7 px-2.5 rounded-md border border-dashed border-gray-300 text-gray-600 hover:bg-gray-50" title="Suggested from the recording: a dialog opened here">＋ “{sg.text}” at {ts(sg.start)}</button>
                ))}
              </div>
              <div className="space-y-1 max-h-36 overflow-y-auto pr-1">
                {captions.map((c, i) => (
                  <div key={i} className="flex items-center gap-1.5 text-gray-700">
                    <span className="w-20 font-semibold text-gray-600">{KIND_LABEL[c.kind]}</span>
                    {startInput(c, (start, end) => editCaption(i, { start, end }))} s
                    {durSelect(c, end => editCaption(i, { end }))}
                    <input value={c.text} onChange={e => editCaption(i, { text: e.target.value })} placeholder={c.kind === 'highlight' ? 'label (optional)' : 'text'} className="flex-1 border border-gray-300 rounded px-2 py-0.5" />
                    {c.kind === 'title' && <input value={c.sub ?? ''} onChange={e => editCaption(i, { sub: e.target.value })} placeholder="subtitle" className="w-56 border border-gray-300 rounded px-2 py-0.5" />}
                    {c.kind === 'highlight' && check(!!c.dim, v => editCaption(i, { dim: v }), 'dim the rest')}
                    <button onClick={() => setCaptions(cs => cs.filter((_, k) => k !== i))} className="text-gray-400 hover:text-red-600 text-base leading-none" title="Remove">×</button>
                  </div>
                ))}
                {!captions.length && <div className="text-gray-400">no caption · pause the preview where one should appear, then add it</div>}
              </div>
            </div>
          )}
          {tab === 'chapters' && (
            <div>
              <div className="flex items-center gap-2 mb-2">
                <button onClick={() => { const t = nowRec(); setChapters(cs => [...cs, { t, name: `Chapter ${cs.length + 1}`, card: true }].sort((a, b) => a.t - b.t)); }} className="h-7 px-2.5 rounded-md border border-accent-300 text-accent-700 hover:bg-accent-50">＋ chapter here</button>
                <span className="text-gray-500">the HTML page lists them, marks them on the time bar and stops at each one: → or PageDown (a presentation clicker) goes on, ← goes back</span>
              </div>
              <div className="space-y-1 max-h-36 overflow-y-auto pr-1">
                {chapters.map((c, i) => (
                  <div key={i} className="flex items-center gap-1.5 text-gray-700">
                    <span className="w-5 font-semibold text-gray-600">{i + 1}.</span>
                    <input type="number" step={0.1} min={0} value={(map(c.t) / 1000).toFixed(1)} className="border border-gray-300 rounded px-1 py-0.5 w-16" title="Time (s)"
                      onChange={e => { const t = invertTime(map, Math.max(0, Number(e.target.value) * 1000), rec.duration); setChapters(cs => cs.map((x, k) => (k === i ? { ...x, t } : x)).sort((a, b) => a.t - b.t)); }} /> s
                    <input value={c.name} onChange={e => setChapters(cs => cs.map((x, k) => (k === i ? { ...x, name: e.target.value } : x)))} className="flex-1 border border-gray-300 rounded px-2 py-0.5" />
                    {check(c.card, v => setChapters(cs => cs.map((x, k) => (k === i ? { ...x, card: v } : x))), 'title card', 'Show the chapter name as a title card when it starts')}
                    <button onClick={() => setChapters(cs => cs.filter((_, k) => k !== i))} className="text-gray-400 hover:text-red-600 text-base leading-none" title="Remove">×</button>
                  </div>
                ))}
                {!chapters.length && <div className="text-gray-400">no chapter · pause the preview where a part of the story starts, then add one</div>}
              </div>
            </div>
          )}
          {tab === 'pointer' && (
            <div className="grid grid-cols-2 gap-6">
              <div className="space-y-1.5">
                <div className="font-semibold text-gray-700">Pointer</div>
                <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
                  {check(ptr.show, v => setPtr({ ...ptr, show: v }), 'shown')}
                  {check(ptr.smooth, v => setPtr({ ...ptr, smooth: v }), 'smoothed')}
                  {check(ptr.clicks, v => setPtr({ ...ptr, clicks: v }), 'click rings')}
                  <select value={ptr.scale} onChange={e => setPtr({ ...ptr, scale: Number(e.target.value) })} className="border border-gray-300 rounded px-1 py-0.5" title="Pointer size">
                    {[1, 1.4, 2, 3].map(s => <option key={s} value={s}>size {s}×</option>)}
                  </select>
                </div>
              </div>
              <div className="space-y-1.5">
                <div className="font-semibold text-gray-700">Playback</div>
                <div>{check(shorten, setShorten, `shorten the pauses to ${MAX_IDLE_DEFAULT_MS / 1000} s`)}</div>
                <div>{check(loop, setLoop, 'loop')}</div>
              </div>
            </div>
          )}
          {tab === 'file' && (
            <div className="grid grid-cols-2 gap-6">
              <div className="space-y-1.5">
                <label className="block"><span className="text-gray-600 mr-2">Title</span><input value={title} onChange={e => setTitle(e.target.value)} className="border border-gray-300 rounded px-2 py-0.5 w-64" /></label>
                <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
                  {(['html', 'svg', 'both'] as const).map(f => (
                    <label key={f} className={`inline-flex items-center gap-1 ${out.hasHtml && f !== 'html' ? 'opacity-40' : ''}`}><input type="radio" disabled={out.hasHtml && f !== 'html'} checked={fmtOut === f} onChange={() => setFormat(f)} className="accent-accent-600" />
                      {f === 'html' ? 'HTML page' : f === 'svg' ? 'SVG' : 'both'}</label>
                  ))}
                </div>
                {out.hasHtml && <div className="text-[11px] text-gray-500">a dialog is recorded with its text and tables: HTML page only</div>}
                {!out.hasHtml && chapters.length > 0 && <div className="text-[11px] text-gray-500">chapters work in the HTML page; the SVG plays straight through</div>}
              </div>
              <div>
                {check(anon, setAnon, 'anonymise the sample names')}
                {anon && (
                  <div className="mt-1 space-y-1 max-h-24 overflow-y-auto">
                    {uniq.map((n, i) => (
                      <div key={n} className="flex items-center gap-2"><span className="font-mono text-gray-500 w-40 truncate" title={n}>{n}</span>→
                        <input value={aliases[i] ?? ''} onChange={e => setAliases(a => { const b = [...a]; b[i] = e.target.value; return b; })} className="border border-gray-300 rounded px-2 py-0.5 w-32" /></div>
                    ))}
                    {!uniq.length && <div className="text-gray-400">no sample open</div>}
                    <div className="text-[11px] text-gray-500">Replaced in every text drawn and in the captions; check the preview before sharing.</div>
                  </div>
                )}
              </div>
            </div>
          )}
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
