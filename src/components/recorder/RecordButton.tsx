/**
 * The Record button: records what the viewer draws, then offers the recording as an animated SVG and / or an HTML
 * page with play / pause, light enough for a presentation.
 */
import { useEffect, useState } from 'react';
import { createPortal } from 'react-dom';
import { startRecording, type Recorder, type Recording } from './recorder';
import RecordingEditor from './RecordingEditor';

const fmtTime = (ms: number) => { const s = Math.max(0, ms) / 1000; return `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`; };

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
      {done && createPortal(<RecordingEditor rec={done} names={names} onClose={() => setDone(null)} onAgain={start} />, document.body)}
    </>
  );
}

