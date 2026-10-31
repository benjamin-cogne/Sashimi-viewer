/**
 * The recorder, driven from the Export menu: a 3-2-1 countdown, then the recording with a REC badge to stop it, then
 * the editor; it also opens a page saved editable back in the editor. Mounted once, so it outlives the menu.
 */
import { forwardRef, useEffect, useImperativeHandle, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { startRecording, type Recorder, type Recording } from './recorder';
import RecordingEditor, { type Edits } from './RecordingEditor';
import { PROJECT_SCRIPT_ID } from './exportRecording';

/** Seconds counted down before the recording starts. */
export const COUNTDOWN_S = 3;

export interface RecorderHandle {
  /** count down, then record */
  start(): void;
  stop(): void;
  /** pick an HTML page saved editable and open it in the editor */
  openFile(): void;
}

const fmtTime = (ms: number) => { const s = Math.max(0, ms) / 1000; return `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`; };

const RecorderHost = forwardRef<RecorderHandle, { names: string[]; onRecording?: (on: boolean) => void }>(function RecorderHost({ names, onRecording }, ref) {
  const [rec, setRec] = useState<Recorder | null>(null);
  const [count, setCount] = useState<number | null>(null);
  const [elapsed, setElapsed] = useState(0);
  const [done, setDone] = useState<{ rec: Recording; initial?: Partial<Edits> } | null>(null);
  const [openErr, setOpenErr] = useState<string | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);
  const recRef = useRef<Recorder | null>(null);
  recRef.current = rec;
  useEffect(() => { onRecording?.(!!rec || count != null); }, [rec, count]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    if (!rec) return;
    const id = setInterval(() => setElapsed(rec.elapsed()), 250);
    return () => clearInterval(id);
  }, [rec]);
  // the countdown: one step a second, then the recording starts (the countdown itself is not recorded)
  useEffect(() => {
    if (count == null) return;
    if (count === 0) { setCount(null); setElapsed(0); setRec(startRecording()); return; }
    const id = setTimeout(() => setCount(c => (c == null ? c : c - 1)), 1000);
    const esc = (e: KeyboardEvent) => { if (e.key === 'Escape') setCount(null); };
    addEventListener('keydown', esc);
    return () => { clearTimeout(id); removeEventListener('keydown', esc); };
  }, [count]);
  const start = () => { if (recRef.current || count != null) return; setDone(null); setCount(COUNTDOWN_S); };
  const stop = () => { const r = recRef.current; if (!r) return; setDone({ rec: r.stop() }); setRec(null); };
  /** A page saved editable: its recording and settings, back into the editor. */
  const openFile = async (f: File) => {
    try {
      const doc = new DOMParser().parseFromString(await f.text(), 'text/html');
      const json = doc.getElementById(PROJECT_SCRIPT_ID)?.textContent;
      if (!json) throw new Error('this page was not saved editable (tick "keep the HTML page editable", next to Download, before downloading)');
      const p = JSON.parse(json);
      if (!p?.rec?.frames) throw new Error('no recording inside');
      setOpenErr(null); setDone({ rec: p.rec as Recording, initial: p.edits as Partial<Edits> });
    } catch (e) { setOpenErr(`Cannot edit ${f.name}: ${(e as Error).message}`); setTimeout(() => setOpenErr(null), 6000); }
  };
  useImperativeHandle(ref, () => ({ start, stop, openFile: () => fileRef.current?.click() }));
  return (
    <>
      <input ref={fileRef} type="file" accept=".html,text/html" className="hidden" onChange={e => { const f = e.target.files?.[0]; if (f) void openFile(f); e.target.value = ''; }} />
      {count != null && count > 0 && createPortal(
        <div data-recorder="skip" className="fixed inset-0 z-[95] grid place-items-center bg-slate-900/25 pointer-events-none">
          <div key={count} className="w-40 h-40 rounded-full bg-gray-900/85 text-white grid place-items-center shadow-2xl animate-[pulse_1s_ease-out_1]">
            <span className="text-[88px] font-bold leading-none tabular-nums">{count}</span>
          </div>
          <div className="absolute bottom-10 left-1/2 -translate-x-1/2 text-white text-sm bg-gray-900/70 rounded-full px-3 py-1">recording starts… · Esc to cancel</div>
        </div>, document.body)}
      {rec && createPortal(
        <div data-recorder="skip" className="fixed bottom-4 right-4 z-[90] flex items-center gap-2 h-9 pl-3 pr-1 rounded-full bg-gray-900/90 text-white text-xs shadow-lg">
          <span className="inline-block w-2.5 h-2.5 rounded-full bg-red-500 animate-pulse" /> REC {fmtTime(elapsed)}
          <button onClick={stop} className="h-7 px-3 rounded-full bg-white text-gray-900 font-semibold hover:bg-gray-100">Stop</button>
        </div>, document.body)}
      {openErr && createPortal(<div data-recorder="skip" className="fixed bottom-4 right-4 z-[90] max-w-md px-3 py-2 rounded-lg bg-red-50 border border-red-200 text-red-800 text-xs shadow">{openErr}</div>, document.body)}
      {done && createPortal(<RecordingEditor rec={done.rec} initial={done.initial} names={done.initial ? [] : names} onClose={() => setDone(null)} onAgain={start} />, document.body)}
    </>
  );
});
export default RecorderHost;
