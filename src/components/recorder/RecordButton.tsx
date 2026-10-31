/**
 * The Record button: records what the viewer draws, then offers the recording as an animated SVG and / or an HTML
 * page with play / pause, light enough for a presentation.
 */
import { useEffect, useMemo, useState } from 'react';
import { createPortal } from 'react-dom';
import { startRecording, type Recorder, type Recording } from './recorder';
import { recordingToSvg, svgToHtml, MAX_IDLE_DEFAULT_MS } from './exportRecording';

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

function ExportDialog({ rec, names, onClose, onAgain }: { rec: Recording; names: string[]; onClose: () => void; onAgain: () => void }) {
  const [format, setFormat] = useState<'svg' | 'html' | 'both'>('html');
  const [shorten, setShorten] = useState(true);
  const [loop, setLoop] = useState(true);
  const [anon, setAnon] = useState(false);
  const [title, setTitle] = useState('Sashimi viewer recording');
  const uniq = useMemo(() => [...new Set(names.filter(Boolean))], [names]);
  const [aliases, setAliases] = useState<string[]>(() => uniq.map((_, i) => `Sample ${i + 1}`));
  const out = useMemo(() => recordingToSvg(rec, {
    maxIdleMs: shorten ? MAX_IDLE_DEFAULT_MS : 0, loop, title,
    replace: anon ? uniq.map((n, i) => ({ from: n, to: aliases[i] ?? `Sample ${i + 1}` })) : [],
  }), [rec, shorten, loop, anon, uniq, aliases, title]);
  const html = useMemo(() => (format === 'svg' ? '' : svgToHtml(out.svg, title, out.duration)), [out, format, title]);
  const [gz, setGz] = useState<number | null>(null);
  useEffect(() => { let alive = true; gzipSize(format === 'svg' ? out.svg : html).then(v => { if (alive) setGz(v); }); return () => { alive = false; }; }, [out, html, format]);
  const preview = useMemo(() => URL.createObjectURL(new Blob([out.svg], { type: 'image/svg+xml' })), [out]);
  useEffect(() => () => URL.revokeObjectURL(preview), [preview]);
  const save = () => {
    const base = `sashimi-recording-${stamp()}`;
    if (format !== 'html') download(out.svg, 'image/svg+xml', `${base}.svg`);
    if (format !== 'svg') download(html, 'text/html', `${base}.html`);
  };
  const bytes = format === 'svg' ? out.bytes : format === 'html' ? new Blob([html]).size : out.bytes + new Blob([html]).size;
  return (
    <div data-recorder="skip" className="fixed inset-0 z-[90] bg-slate-900/70 flex items-start justify-center overflow-y-auto p-4" onMouseDown={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-[1100px] text-gray-900 text-sm" onMouseDown={e => e.stopPropagation()}>
        <div className="flex items-start justify-between px-5 pt-4 pb-2">
          <div>
            <div className="text-base font-bold">Recording <span className="font-normal text-gray-500">· {fmtTime(out.duration)} · {out.frames} frames · {fmtSize(bytes)}{gz ? ` (${fmtSize(gz)} gzipped)` : ''}</span></div>
            <div className="text-xs text-gray-500">vectors only: sharp at any zoom, text kept as text · the toolbar and the tables are not recorded, dialogs show as their panel</div>
          </div>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-2xl leading-none px-1" title="Close">×</button>
        </div>
        <div className="px-5">
          <div className="rounded-lg border border-gray-200 bg-gray-50 overflow-hidden"><img src={preview} alt="preview of the recording" className="block w-full h-auto max-h-[55vh] object-contain" /></div>
        </div>
        <div className="px-5 py-3 grid grid-cols-[1fr_1fr] gap-6 text-xs">
          <div className="space-y-2">
            <label className="block"><span className="text-gray-600 mr-2">Title</span><input value={title} onChange={e => setTitle(e.target.value)} className="border border-gray-300 rounded px-2 py-1 w-72" /></label>
            <div className="flex items-center gap-3">
              <span className="text-gray-600">Format</span>
              {(['html', 'svg', 'both'] as const).map(f => (
                <label key={f} className="inline-flex items-center gap-1"><input type="radio" checked={format === f} onChange={() => setFormat(f)} className="accent-accent-600" />
                  {f === 'html' ? 'HTML page (play / pause, time bar)' : f === 'svg' ? 'SVG (plays on its own)' : 'both'}</label>
              ))}
            </div>
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={shorten} onChange={e => setShorten(e.target.checked)} className="accent-accent-600" />Shorten the pauses to {MAX_IDLE_DEFAULT_MS / 1000} s</label>
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={loop} onChange={e => setLoop(e.target.checked)} className="accent-accent-600" />Loop</label>
          </div>
          <div>
            <label className="flex items-center gap-1.5"><input type="checkbox" checked={anon} onChange={e => setAnon(e.target.checked)} className="accent-accent-600" />Anonymise the sample names</label>
            {anon && (
              <div className="mt-1 space-y-1 max-h-32 overflow-y-auto">
                {uniq.map((n, i) => (
                  <div key={n} className="flex items-center gap-2"><span className="font-mono text-gray-500 w-48 truncate" title={n}>{n}</span>→
                    <input value={aliases[i] ?? ''} onChange={e => setAliases(a => { const b = [...a]; b[i] = e.target.value; return b; })} className="border border-gray-300 rounded px-2 py-0.5 w-40" /></div>
                ))}
                {!uniq.length && <div className="text-gray-400">no sample open</div>}
                <div className="text-[11px] text-gray-500">Replaced in every text drawn; check the preview before sharing.</div>
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
