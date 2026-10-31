/**
 * The Record button: records what the viewer draws, then offers the recording as an animated SVG and / or an HTML
 * page with play / pause, light enough for a presentation.
 */
import { useEffect, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { startRecording, type Recorder, type Recording } from './recorder';
import RecordingEditor, { type Edits } from './RecordingEditor';
import { PROJECT_SCRIPT_ID } from './exportRecording';

const fmtTime = (ms: number) => { const s = Math.max(0, ms) / 1000; return `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`; };

export default function RecordButton({ names, className }: { names: string[]; className: string }) {
  const [rec, setRec] = useState<Recorder | null>(null);
  const [elapsed, setElapsed] = useState(0);
  const [done, setDone] = useState<{ rec: Recording; initial?: Partial<Edits> } | null>(null);
  const [openErr, setOpenErr] = useState<string | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);
  /** A page saved editable: its recording and settings, back into the editor. */
  const openFile = async (f: File) => {
    try {
      const doc = new DOMParser().parseFromString(await f.text(), 'text/html');
      const json = doc.getElementById(PROJECT_SCRIPT_ID)?.textContent;
      if (!json) throw new Error('this page was not saved editable (tick "keep the HTML page editable" before downloading)');
      const p = JSON.parse(json);
      if (!p?.rec?.frames) throw new Error('no recording inside');
      setOpenErr(null); setDone({ rec: p.rec as Recording, initial: p.edits as Partial<Edits> });
    } catch (e) { setOpenErr(`Cannot edit ${f.name}: ${(e as Error).message}`); setTimeout(() => setOpenErr(null), 6000); }
  };
  useEffect(() => {
    if (!rec) return;
    const id = setInterval(() => setElapsed(rec.elapsed()), 250);
    return () => clearInterval(id);
  }, [rec]);
  const start = () => { setDone(null); setElapsed(0); setRec(startRecording()); };
  const stop = () => { if (!rec) return; setDone({ rec: rec.stop() }); setRec(null); };
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
      {!rec && (
        <button onClick={() => fileRef.current?.click()} className={className} data-recorder="skip" title="Edit a recording again: open an HTML page saved editable">✎</button>
      )}
      <input ref={fileRef} type="file" accept=".html,text/html" className="hidden" onChange={e => { const f = e.target.files?.[0]; if (f) void openFile(f); e.target.value = ''; }} />
      {openErr && createPortal(<div data-recorder="skip" className="fixed bottom-4 right-4 z-[90] max-w-md px-3 py-2 rounded-lg bg-red-50 border border-red-200 text-red-800 text-xs shadow">{openErr}</div>, document.body)}
      {done && createPortal(<RecordingEditor rec={done.rec} initial={done.initial} names={done.initial ? [] : names} onClose={() => setDone(null)} onAgain={start} />, document.body)}
    </>
  );
}

