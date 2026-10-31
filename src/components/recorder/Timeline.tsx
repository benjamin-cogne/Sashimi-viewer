/**
 * The editor's timeline, in recorded time: lanes of blocks (cuts and speed changes, zooms, captions) and of marks
 * (chapters), the stretches the pause shortening drops (hatched), the recording's activity (changes, clicks), and the
 * playhead. A block is moved by its body and resized by its edges, a mark moved; the change is applied when the mouse
 * is released. A click on an empty place moves the playhead there. For long recordings the timeline zooms (the
 * buttons, or Ctrl / ⌘ + wheel around the pointer) and pans (the wheel sideways, Shift + wheel, or the overview bar
 * under it); while zoomed in, it follows the playhead.
 */
import { useEffect, useRef, useState } from 'react';

export interface TimelineItem { id: string; start: number; end: number; label: string; color: string; faded?: boolean }
export interface TimelineLane { key: string; label: string; items: TimelineItem[]; marks?: boolean }

const fmt = (ms: number) => (ms >= 60000 ? `${Math.floor(ms / 60000)}:${String(Math.floor((ms % 60000) / 1000)).padStart(2, '0')}` : `${(ms / 1000).toFixed(ms < 10000 ? 1 : 0)} s`);

export default function Timeline({ duration, lanes, idle, activity, clicks, playhead, selected, onSelect, onSeek, onChange }: {
  duration: number; lanes: TimelineLane[];
  idle: { start: number; end: number }[]; activity: number[]; clicks: number[];
  playhead: number; selected: string | null;
  onSelect: (id: string | null) => void; onSeek: (t: number) => void;
  onChange: (lane: string, id: string, start: number, end: number) => void;
}) {
  const box = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(1000);
  useEffect(() => {
    const el = box.current;
    if (!el) return;
    const ro = new ResizeObserver(() => setWidth(el.clientWidth - LABEL_W));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);
  const LABEL_W = 92, LANE_H = 22;
  const D = Math.max(1, duration);
  // the stretch shown: [v0, v1] of recorded time
  const [view, setView] = useState<[number, number]>([0, D]);
  useEffect(() => setView([0, D]), [D]);
  const [v0, v1] = view, span = v1 - v0;
  const MIN_SPAN = Math.min(D, 1000);
  const clampView = (a: number, len: number): [number, number] => { len = Math.max(MIN_SPAN, Math.min(D, len)); a = Math.max(0, Math.min(D - len, a)); return [a, a + len]; };
  const zoomAt = (factor: number, at: number) => setView(([a, b]) => { const len = (b - a) / factor, f = (at - a) / (b - a || 1); return clampView(at - f * len, len); });
  const x = (t: number) => ((t - v0) / span) * width;
  const tOf = (px: number) => Math.max(0, Math.min(D, v0 + (px / width) * span));
  const lanesRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const el = lanesRef.current;
    if (!el) return;
    const wheel = (ev: WheelEvent) => {
      const r = el.getBoundingClientRect(), px = ev.clientX - r.left - LABEL_W;
      if (ev.ctrlKey || ev.metaKey) { ev.preventDefault(); zoomAt(Math.exp(-ev.deltaY * 0.002), tOf(px)); }
      else if (span < D && (Math.abs(ev.deltaX) > Math.abs(ev.deltaY) || ev.shiftKey)) { ev.preventDefault(); const d = (ev.shiftKey ? ev.deltaY : ev.deltaX) / width * span; setView(([a, b]) => clampView(a + d, b - a)); }
    };
    el.addEventListener('wheel', wheel, { passive: false });
    return () => el.removeEventListener('wheel', wheel);
  });
  // zoomed in: keep the playhead in view (a page further when it leaves on the right)
  useEffect(() => {
    if (span >= D) return;
    if (playhead > v1 || playhead < v0) setView(clampView(playhead - span * 0.1, span));
  }, [playhead]); // eslint-disable-line react-hooks/exhaustive-deps
  // the overview bar: drag the window to pan
  const [pan, setPan] = useState<{ x0: number; a0: number } | null>(null);
  useEffect(() => {
    if (!pan) return;
    const move = (ev: PointerEvent) => setView(([a, b]) => clampView(pan.a0 + ((ev.clientX - pan.x0) / width) * D, b - a));
    const up = () => setPan(null);
    addEventListener('pointermove', move); addEventListener('pointerup', up, { once: true });
    return () => { removeEventListener('pointermove', move); removeEventListener('pointerup', up); };
  }, [pan]); // eslint-disable-line react-hooks/exhaustive-deps
  // a drag in progress: the item and its provisional times
  const [drag, setDrag] = useState<{ lane: string; id: string; mode: 'move' | 'start' | 'end'; x0: number; s0: number; e0: number; s: number; e: number } | null>(null);
  useEffect(() => {
    if (!drag) return;
    const move = (ev: PointerEvent) => {
      const dt = ((ev.clientX - drag.x0) / width) * span;
      let s = drag.s0, e = drag.e0;
      if (drag.mode === 'move') { const len = e - s; s = Math.max(0, Math.min(D - len, s + dt)); e = s + len; }
      else if (drag.mode === 'start') s = Math.max(0, Math.min(e - 200, s + dt));
      else e = Math.min(D, Math.max(s + 200, e + dt));
      setDrag(d => (d ? { ...d, s, e } : d));
    };
    const up = () => { setDrag(d => { if (d && (d.s !== d.s0 || d.e !== d.e0)) onChange(d.lane, d.id, d.s, d.e); return null; }); };
    addEventListener('pointermove', move);
    addEventListener('pointerup', up, { once: true });
    return () => { removeEventListener('pointermove', move); removeEventListener('pointerup', up); };
  }, [drag?.id, drag?.mode]); // eslint-disable-line react-hooks/exhaustive-deps
  const grab = (lane: string, it: TimelineItem, mode: 'move' | 'start' | 'end') => (ev: React.PointerEvent) => {
    ev.stopPropagation(); ev.preventDefault();
    onSelect(it.id);
    setDrag({ lane, id: it.id, mode, x0: ev.clientX, s0: it.start, e0: it.end, s: it.start, e: it.end });
  };
  const seekAt = (ev: React.PointerEvent) => {
    const r = (ev.currentTarget as HTMLElement).getBoundingClientRect();
    onSelect(null); onSeek(tOf(ev.clientX - r.left));
  };
  // ruler ticks: about one every 80 px, on round seconds
  const step = [100, 200, 500, 1000, 2000, 5000, 10000, 15000, 30000, 60000, 120000, 300000].find(s => (s / span) * width >= 70) ?? 300000;
  const ticks: number[] = [];
  for (let t = Math.ceil(v0 / step) * step; t <= v1; t += step) ticks.push(t);
  const zoomed = span < D - 1;

  return (
    <div ref={box} className="select-none text-[10.5px]">
      <div className="flex">
        <div style={{ width: LABEL_W }} className="flex items-center gap-0.5">
          <button onClick={() => zoomAt(1 / 1.6, (v0 + v1) / 2)} disabled={!zoomed} className="w-5 h-5 rounded border border-gray-300 text-gray-600 hover:bg-gray-50 disabled:opacity-30" title="Zoom out">−</button>
          <button onClick={() => zoomAt(1.6, playhead >= v0 && playhead <= v1 ? playhead : (v0 + v1) / 2)} disabled={span <= MIN_SPAN} className="w-5 h-5 rounded border border-gray-300 text-gray-600 hover:bg-gray-50 disabled:opacity-30" title="Zoom in (Ctrl / ⌘ + wheel)">+</button>
          <button onClick={() => setView([0, D])} disabled={!zoomed} className="h-5 px-1 rounded border border-gray-300 text-gray-600 hover:bg-gray-50 disabled:opacity-30" title="The whole recording">fit</button>
        </div>
        <div className="relative h-5 cursor-pointer overflow-hidden" style={{ width }} onPointerDown={seekAt}>
          {ticks.map(t => <div key={t} className="absolute top-0 h-full border-l border-gray-300 text-gray-500 pl-0.5" style={{ left: x(t) }}>{fmt(t)}</div>)}
          {activity.map((t, i) => <div key={`a${i}`} className="absolute bottom-0 w-px h-1.5 bg-gray-400" style={{ left: x(t) }} />)}
          {clicks.map((t, i) => <div key={`c${i}`} className="absolute bottom-0 w-1.5 h-1.5 -ml-[3px] rounded-full bg-accent-600" style={{ left: x(t) }} title="click" />)}
        </div>
      </div>
      <div className="relative overflow-hidden" ref={lanesRef}>
        {lanes.map(lane => (
          <div key={lane.key} className="flex items-center border-t border-gray-100">
            <div className="text-gray-600 font-semibold truncate pr-1" style={{ width: LABEL_W }}>{lane.label}</div>
            <div className="relative cursor-pointer bg-gray-50 overflow-hidden" style={{ width, height: LANE_H }} onPointerDown={seekAt}>
              {lane.key === 'edit' && idle.map((c, i) => (
                <div key={`i${i}`} className="absolute top-0 h-full" title="dropped by the pause shortening"
                  style={{ left: x(c.start), width: Math.max(1, x(c.end) - x(c.start)), background: 'repeating-linear-gradient(45deg,#e2e8f0 0 4px,#f8fafc 4px 8px)' }} />
              ))}
              {lane.items.map(it => {
                const d = drag && drag.id === it.id ? drag : null;
                const s = d ? d.s : it.start, e = d ? d.e : it.end;
                const sel = selected === it.id;
                if (lane.marks) return (
                  <div key={it.id} className="absolute top-0 h-full cursor-ew-resize" style={{ left: x(s) - 6, width: 12 }} onPointerDown={grab(lane.key, { ...it, end: it.start + 1 }, 'move')} title={`${it.label} · ${fmt(s)}`}>
                    <div className="absolute left-[5px] top-0 h-full w-0.5" style={{ background: it.color }} />
                    <div className={`absolute left-[7px] top-0.5 px-1 rounded text-white whitespace-nowrap ${sel ? 'ring-2 ring-offset-1 ring-gray-800' : ''}`} style={{ background: it.color }}>{it.label}</div>
                  </div>
                );
                return (
                  <div key={it.id} className={`absolute top-0.5 rounded cursor-grab overflow-hidden text-white whitespace-nowrap ${sel ? 'ring-2 ring-offset-1 ring-gray-800 z-10' : ''}`}
                    style={{ left: x(s), width: Math.max(4, x(e) - x(s)), height: LANE_H - 4, background: it.color, opacity: it.faded ? 0.45 : 1 }}
                    onPointerDown={grab(lane.key, it, 'move')} title={`${it.label} · ${fmt(s)} → ${fmt(e)}`}>
                    <div className="absolute left-0 top-0 h-full w-1.5 cursor-ew-resize" onPointerDown={grab(lane.key, it, 'start')} />
                    <span className="px-2 leading-[18px]">{it.label}</span>
                    <div className="absolute right-0 top-0 h-full w-1.5 cursor-ew-resize" onPointerDown={grab(lane.key, it, 'end')} />
                  </div>
                );
              })}
            </div>
          </div>
        ))}
        {playhead >= v0 && playhead <= v1 && <div className="absolute top-0 bottom-0 w-0.5 bg-red-500 pointer-events-none" style={{ left: LABEL_W + x(playhead) }} />}
      </div>
      {zoomed && (
        <div className="flex items-center mt-1">
          <div className="text-gray-500" style={{ width: LABEL_W }}>{fmt(v0)} – {fmt(v1)}</div>
          <div className="relative h-2.5 rounded bg-gray-100" style={{ width }}>
            <div className="absolute top-0 h-full w-0.5 bg-red-400" style={{ left: (playhead / D) * width }} />
            <div className="absolute top-0 h-full rounded bg-accent-500/40 border border-accent-600 cursor-grab" style={{ left: (v0 / D) * width, width: Math.max(6, (span / D) * width) }}
              onPointerDown={ev => { ev.preventDefault(); setPan({ x0: ev.clientX, a0: v0 }); }} title="Drag to move along the recording" />
          </div>
        </div>
      )}
    </div>
  );
}
