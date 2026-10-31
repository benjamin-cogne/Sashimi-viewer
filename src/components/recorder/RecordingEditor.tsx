/**
 * The editor that opens when a recording stops: a live preview (the exported HTML page itself) and tabs to set the
 * zooms, the captions, the chapters, the pointer and the file, then download.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import type { Recording } from './recorder';
import { recordingToSvg, svgToHtml, autoZooms, timeMap, invertTime, idleCuts, replaceText, MAX_IDLE_DEFAULT_MS, ZOOM_LEVEL_DEFAULT, type Zoom } from './exportRecording';
import Timeline, { type TimelineLane } from './Timeline';
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

/** Everything the editor sets: saved with the page when it is kept editable, and read back to edit it again. */
export interface Edits {
  title: string; format: 'svg' | 'html' | 'both'; shorten: boolean; loop: boolean; editable: boolean;
  anon: boolean; aliases: Record<string, string>;
  autoOn: boolean; level: number; autoOff: number[]; drawn: Zoom[];
  ptr: { show: boolean; smooth: boolean; scale: number; clicks: boolean };
  captions: Caption[]; chapters: Chapter[];
  cuts: { start: number; end: number }[]; speeds: { start: number; end: number; factor: number }[];
}
const DEFAULT_EDITS: Edits = {
  title: 'Sashimi viewer recording', format: 'html', shorten: true, loop: true, editable: false, anon: false, aliases: {},
  autoOn: true, level: ZOOM_LEVEL_DEFAULT, autoOff: [], drawn: [], ptr: { show: true, smooth: true, scale: 1.4, clicks: true },
  captions: [], chapters: [], cuts: [], speeds: [],
};
type Upd<T> = T | ((prev: T) => T);

/** The recording with the replacements made in every text it holds, for a page kept editable after anonymisation. */
function anonymised(rec: Recording, rep: { from: string; to: string }[]): Recording {
  if (!rep.length) return rec;
  return {
    ...rec,
    chunks: rec.chunks.map(c => ('leaf' in c ? { leaf: replaceText(c.leaf, rep) } : { ...c, open: replaceText(c.open, rep) })),
    html: rec.html.map(h => ({ ...h, html: replaceText(h.html, rep) })),
  };
}

export default function RecordingEditor({ rec, names, initial, onClose, onAgain }: { rec: Recording; names: string[]; initial?: Partial<Edits>; onClose: () => void; onAgain: () => void }) {
  const [tab, setTab] = useState<Tab>('zooms');
  const [ed, setEd] = useState<Edits>(() => ({ ...DEFAULT_EDITS, ...initial }));
  const field = <K extends keyof Edits>(k: K) => (v: Upd<Edits[K]>) => setEd(e => ({ ...e, [k]: typeof v === 'function' ? (v as (p: Edits[K]) => Edits[K])(e[k]) : v }));
  const { format, shorten, loop, anon, title, autoOn, level, drawn, ptr, captions, chapters, cuts, speeds, editable } = ed;
  const [setFormat, setShorten, setLoop, setAnon, setTitle, setAutoOn, setLevel, setDrawn, setPtr, setCaptions, setChapters, setCuts, setSpeeds, setEditable] =
    [field('format'), field('shorten'), field('loop'), field('anon'), field('title'), field('autoOn'), field('level'), field('drawn'), field('ptr'), field('captions'), field('chapters'), field('cuts'), field('speeds'), field('editable')];
  const autoOff = useMemo(() => new Set(ed.autoOff), [ed.autoOff]);
  const setAutoOff = (f: (s: Set<number>) => Set<number>) => setEd(e => ({ ...e, autoOff: [...f(new Set(e.autoOff))] }));
  const uniq = useMemo(() => [...new Set(names.filter(Boolean))], [names]);
  const aliases = uniq.map((n, i) => ed.aliases[n] ?? `Sample ${i + 1}`);
  const setAliases = (f: (a: string[]) => string[]) => { const next = f(aliases); setEd(e => ({ ...e, aliases: Object.fromEntries(uniq.map((n, i) => [n, next[i]])) })); };
  const replace = useMemo(() => (anon ? uniq.map((n, i) => ({ from: n, to: aliases[i] })) : []), [anon, uniq, ed.aliases]); // eslint-disable-line react-hooks/exhaustive-deps
  const suggestions = useMemo(() => suggestCaptions(rec), [rec]);
  const autos = useMemo(() => autoZooms(rec, level), [rec, level]);
  const zooms = useMemo(() => [...(autoOn ? autos.filter((_, i) => !autoOff.has(i)) : []), ...drawn], [autos, autoOn, autoOff, drawn]);
  const time = useMemo(() => ({ maxIdleMs: shorten ? MAX_IDLE_DEFAULT_MS : 0, cuts, speeds }), [shorten, cuts, speeds]);
  const map = useMemo(() => timeMap(rec, time), [rec, time]);
  const opts = useMemo(() => ({ time, loop, title, zooms, pointer: ptr, captions, chapters, replace }), [time, loop, title, zooms, ptr, captions, chapters, replace]);
  const out = useMemo(() => recordingToSvg(rec, opts), [rec, opts]);
  const project = useMemo(() => (editable ? JSON.stringify({ v: 1, rec: anonymised(rec, replace), edits: { ...ed, anon: false, aliases: {}, captions: ed.captions.map(c => ({ ...c, text: replaceText(`>${c.text}<`, replace).slice(1, -1) })) } }) : undefined), [editable, rec, replace, ed]);
  const html = useMemo(() => svgToHtml(out, title, { project }), [out, title, project]);
  const preview = useMemo(() => (project ? svgToHtml(out, title) : html), [out, title, project, html]);
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
  useEffect(() => { ready.current = false; }, [preview]);
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

  // ---- the timeline ----
  const [playhead, setPlayhead] = useState(0);
  useEffect(() => { const id = setInterval(() => setPlayhead(invertTime(map, keep.current.t, rec.duration)), 200); return () => clearInterval(id); }, [map, rec.duration]);
  const [sel, setSel] = useState<string | null>(null);
  const CAP_COLOR: Record<CaptionKind, string> = { title: '#0f1b2b', lower: '#475569', callout: '#14857a', highlight: '#7c3aed' };
  const lanes: TimelineLane[] = [
    { key: 'edit', label: 'Cut · speed', items: [
      ...cuts.map((c, i) => ({ id: `cut:${i}`, start: c.start, end: c.end, label: '✂ cut', color: '#dc2626' })),
      ...speeds.map((c, i) => ({ id: `sp:${i}`, start: c.start, end: c.end, label: `${c.factor < 1 ? '◀' : '⏩'} ×${c.factor}`, color: '#d97706' })),
    ] },
    { key: 'zoom', label: 'Zooms', items: [
      ...(autoOn ? autos.map((z, i) => ({ id: `az:${i}`, start: z.start, end: z.end, label: `click ${level}×`, color: '#4fb3a6', faded: autoOff.has(i) })) : []),
      ...drawn.map((z, i) => ({ id: `dz:${i}`, start: z.start, end: z.end, label: 'zoom', color: '#0f766e' })),
    ] },
    { key: 'cap', label: 'Captions', items: captions.map((c, i) => ({ id: `cp:${i}`, start: c.start, end: c.end, label: c.text || KIND_LABEL[c.kind], color: CAP_COLOR[c.kind] })) },
    { key: 'ch', label: 'Chapters', marks: true, items: chapters.map((c, i) => ({ id: `ch:${i}`, start: c.t, end: c.t + 1, label: `${i + 1}. ${c.name}`, color: '#0f766e' })) },
  ];
  const onTimelineChange = (_lane: string, id: string, start: number, end: number) => {
    const [kind, n] = id.split(':'), i = Number(n);
    if (kind === 'cut') setCuts(l => l.map((x, k) => (k === i ? { start, end } : x)));
    else if (kind === 'sp') setSpeeds(l => l.map((x, k) => (k === i ? { ...x, start, end } : x)));
    else if (kind === 'az') { const z = autos[i]; setAutoOff(s => new Set([...s, i])); setDrawn(l => [...l, { ...z, start, end, source: 'drawn' }]); setSel(`dz:${drawn.length}`); }
    else if (kind === 'dz') setDrawn(l => l.map((x, k) => (k === i ? { ...x, start, end } : x)));
    else if (kind === 'cp') editCaption(i, { start, end });
    else if (kind === 'ch') setChapters(l => l.map((x, k) => (k === i ? { ...x, t: start } : x)).sort((a, b) => a.t - b.t));
  };
  const removeSelected = () => {
    if (!sel) return;
    const [kind, n] = sel.split(':'), i = Number(n);
    const drop = <T,>(l: T[]) => l.filter((_, k) => k !== i);
    if (kind === 'cut') setCuts(drop); else if (kind === 'sp') setSpeeds(drop);
    else if (kind === 'az') setAutoOff(s => new Set([...s, i])); else if (kind === 'dz') setDrawn(drop);
    else if (kind === 'cp') setCaptions(drop); else if (kind === 'ch') setChapters(drop);
    setSel(null);
  };
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.key === 'Delete' || e.key === 'Backspace') && sel && !(e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement || e.target instanceof HTMLSelectElement)) { e.preventDefault(); removeSelected(); }
    };
    addEventListener('keydown', onKey);
    return () => removeEventListener('keydown', onKey);
  });
  const seekRec = (t: number) => { const p = player(); const pt = map(t); keep.current = { t: pt, paused: true }; if (p) { p.pause(); p.seek(pt); } setPlayhead(t); };
  const addAtPlayhead = (what: 'cut' | 'speed') => {
    const t = nowRec(), e = Math.min(rec.duration, t + 2000);
    if (e - t < 200) return;
    if (what === 'cut') { setCuts(l => [...l, { start: t, end: e }]); setSel(`cut:${cuts.length}`); }
    else { setSpeeds(l => [...l, { start: t, end: e, factor: 2 }]); setSel(`sp:${speeds.length}`); }
  };
  const selSpeed = sel?.startsWith('sp:') ? Number(sel.slice(3)) : -1;
  const idle = useMemo(() => idleCuts(rec, shorten ? MAX_IDLE_DEFAULT_MS : 0), [rec, shorten]);

  // ---- video: the drawing played full window, captured as a video of this tab ----
  const [video, setVideo] = useState<{ html: string; duration: number; stream: MediaStream } | null>(null);
  const [videoMsg, setVideoMsg] = useState<string | null>(null);
  const videoFrame = useRef<HTMLIFrameElement>(null);
  const startVideo = async () => {
    setVideoMsg(null);
    const md = navigator.mediaDevices as any;
    if (!md?.getDisplayMedia || typeof MediaRecorder === 'undefined') { setVideoMsg('This browser cannot record its tab: use Chrome or Edge on a computer.'); return; }
    let stream: MediaStream;
    try {
      // asked first, while the click still counts as the user's gesture
      stream = await md.getDisplayMedia({ video: { frameRate: 30 }, audio: false, preferCurrentTab: true, selfBrowserSurface: 'include', surfaceSwitching: 'exclude' });
    } catch (e) { setVideoMsg(`Video cancelled: the tab was not shared (${(e as Error).name}${(e as Error).message ? `: ${(e as Error).message}` : ''}).`); return; }
    const vOut = recordingToSvg(rec, { ...opts, loop: false });
    setVideo({ html: svgToHtml(vOut, title, { bare: true }), duration: vOut.duration, stream });
  };
  const runVideo = async () => {
    if (!video) return;
    const { stream, duration } = video;
    const track = stream.getVideoTracks()[0];
    const sleep = (ms: number) => new Promise(r => setTimeout(r, ms));
    try {
      const CT = (window as any).CropTarget;
      if (CT?.fromElement && videoFrame.current && (track as any).cropTo) await (track as any).cropTo(await CT.fromElement(videoFrame.current));
      await sleep(700);
      // H.264 in MP4 is what PowerPoint plays everywhere; else WebM (VP9); an MP4 of another codec last
      const mime = ['video/mp4;codecs=avc1.640028', 'video/mp4;codecs=avc1.42E01E', 'video/webm;codecs=vp9', 'video/webm', 'video/mp4'].find(m => MediaRecorder.isTypeSupported(m)) ?? '';
      const mr = new MediaRecorder(stream, { ...(mime ? { mimeType: mime } : {}), videoBitsPerSecond: 8_000_000 });
      const parts: Blob[] = [];
      mr.ondataavailable = e => { if (e.data.size) parts.push(e.data); };
      const stopped = new Promise(r => { mr.onstop = r; });
      mr.start(250);
      (videoFrame.current?.contentWindow as any)?.recPlayer?.play();
      await sleep(duration + 500);
      mr.stop(); await stopped;
      const type = mr.mimeType || mime || 'video/webm';
      const blob = new Blob(parts, { type });
      const url = URL.createObjectURL(blob), a = document.createElement('a');
      a.href = url; a.download = `sashimi-recording-${stamp()}.${type.includes('mp4') ? 'mp4' : 'webm'}`; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 4000);
      setVideoMsg(`Video saved: ${type.includes('mp4') ? (type.includes('avc1') ? 'MP4 (H.264)' : 'MP4') : 'WebM (VP9)'}, ${fmtSize(blob.size)}${CT?.fromElement ? '' : ' (this browser records the whole tab)'}.`);
    } catch (e) { setVideoMsg(`Video failed: ${(e as Error).message}`); }
    finally { stream.getTracks().forEach(t => t.stop()); setVideo(null); }
  };

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
            <iframe ref={frame} title="preview" srcDoc={preview} onLoad={onFrameLoad} className="w-full h-full block" />
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
        <div className="px-5 pt-2">
          <div className="flex items-center gap-2 text-xs mb-1">
            <span className="font-semibold text-gray-700">Timeline</span>
            <button onClick={() => addAtPlayhead('cut')} className="h-6 px-2 rounded border border-red-300 text-red-700 hover:bg-red-50" title="Remove 2 s from the playhead; drag its edges to fit">✂ cut</button>
            <button onClick={() => addAtPlayhead('speed')} className="h-6 px-2 rounded border border-amber-300 text-amber-700 hover:bg-amber-50" title="Play 2 s from the playhead faster or slower; drag its edges to fit">⏩ speed</button>
            {selSpeed >= 0 && speeds[selSpeed] && (
              <select value={speeds[selSpeed].factor} onChange={e => setSpeeds(l => l.map((x, k) => (k === selSpeed ? { ...x, factor: Number(e.target.value) } : x)))} className="border border-gray-300 rounded px-1 py-0.5">
                {[0.5, 1.5, 2, 4, 8].map(f => <option key={f} value={f}>×{f}{f < 1 ? ' (slower)' : ''}</option>)}
              </select>
            )}
            {sel && <button onClick={removeSelected} className="h-6 px-2 rounded border border-gray-300 text-gray-600 hover:bg-gray-50" title="Delete (Del)">delete</button>}
            <span className="ml-auto text-[11px] text-gray-500">recorded time · drag a block to move it, its edges to resize · click to place the playhead · Del removes the selection</span>
          </div>
          <Timeline duration={rec.duration} lanes={lanes} idle={idle} activity={rec.frames.map(f => f.t)} clicks={rec.clicks.map(c => c.t)}
            playhead={playhead} selected={sel} onSelect={setSel} onSeek={seekRec} onChange={onTimelineChange} />
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
                <select value={level} onChange={e => { setLevel(Number(e.target.value)); setAutoOff(() => new Set()); }} className="border border-gray-300 rounded px-1 py-0.5">
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
                <div>{check(editable, setEditable, 'keep the HTML page editable', 'The recording and these settings go into the page: open it again with ✎ next to Record to change it. About doubles the size.')}</div>
                <div className="flex items-center gap-2 pt-1">
                  <button onClick={startVideo} disabled={!!video} className="h-7 px-2.5 rounded-md border border-gray-300 text-gray-700 hover:bg-gray-50 disabled:opacity-50" title="For PowerPoint or Google Slides: the recording played full window and captured as a video of this tab (the browser asks to share it)">🎞 save as a video…</button>
                  <span className="text-[11px] text-gray-500">MP4 (H.264) where the browser can, else WebM · Chrome or Edge on a computer · plays in real time</span>
                </div>
                {videoMsg && <div className="text-[11px] text-gray-700">{videoMsg}</div>}
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
      {video && <iframe ref={videoFrame} title="video" srcDoc={video.html} onLoad={runVideo} className="fixed inset-0 w-screen h-screen z-[100] bg-white border-0" />}
    </div>
  );
}
