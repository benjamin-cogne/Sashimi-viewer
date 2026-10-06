/**
 * The controls of the viewer's toolbars: one look for the whole page.
 *
 * - `Segmented`: one choice among a few (Reads | Usage, shared | own | % max), the active option a raised white chip.
 * - `Pill`: an on/off option, tinted indigo when on.
 * - `Stepper`: a number with − / + buttons, typed in place.
 * - `Select`: a native select dressed as a pill, with its key in grey ("window ≤ 100 kb").
 * - `Section`: an uppercase label and the controls it groups, on a pastel panel of its tone; `Sep` a thin rule.
 * - `Popover` / `MenuItem` / `SwitchRow`: a panel under a button, closed by a click outside or Escape.
 *
 * They live at module level, not inside the viewer's render: a component declared in a render is a new type on every
 * render, so React would remount it and a number being typed would lose its focus.
 */
import { useEffect, useRef, useState } from 'react';

const ICON_PATHS: Record<string, string> = {
  search: 'M11 4a7 7 0 1 0 0 14 7 7 0 0 0 0-14zM20 20l-3.5-3.5',
  plus: 'M12 5v14M5 12h14',
  minus: 'M5 12h14',
  reset: 'M3 12a9 9 0 1 0 3-6.7L3 8M3 3v5h5',
  chev: 'm6 9 6 6 6-6',
  layers: 'm12 2 9 5-9 5-9-5 9-5zM3 12l9 5 9-5M3 17l9 5 9-5',
  filter: 'M3 5h18l-7 8v6l-4 2v-8L3 5z',
  users: 'M9 4.5a3.5 3.5 0 1 0 0 7 3.5 3.5 0 0 0 0-7zM2.5 20a6.5 6.5 0 0 1 13 0M16 4.5a3.5 3.5 0 0 1 0 7M18 20a6 6 0 0 0-3-5.2',
  grid: 'M4 4h6v6H4zM14 4h6v6h-6zM4 14h6v6H4zM14 14h6v6h-6z',
  download: 'M12 3v12M7 10l5 5 5-5M4 21h16',
  image: 'M5 4h14a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V6a2 2 0 0 1 2-2zM9 8a2 2 0 1 0 0 4 2 2 0 0 0 0-4zM21 17l-5-5-9 9',
  file: 'M14 3H6a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V9zM14 3v6h6',
  folder: 'M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z',
  save: 'M5 3h11l5 5v11a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2zM7 3v5h8M7 21v-7h10v7',
  upload: 'M12 21V9M7 14l5-5 5 5M4 3h16',
  arc: 'M3 18c2-9 16-9 18 0',
  pct: 'M19 5 5 19M6.5 4a2.5 2.5 0 1 0 0 5 2.5 2.5 0 0 0 0-5zM17.5 15a2.5 2.5 0 1 0 0 5 2.5 2.5 0 0 0 0-5z',
  reads: 'M3 7h12M7 11h14M3 15h10M9 19h12',
  tx: 'M2 12h20M5 9h4v6H5zM15 9h4v6h-4z',
  snp: 'M12 11v10M12 5a3 3 0 1 0 0 6 3 3 0 0 0 0-6z',
  equal: 'M4 9h16M4 15h16',
  more: 'M5 11v2M12 11v2M19 11v2',
  info: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM12 11v6M12 7.5v.5',
  pair: 'M3 12h6M15 12h6M10 12h1M13 12h1',
  scissors: 'M6 3a3 3 0 1 0 0 6 3 3 0 0 0 0-6zM6 15a3 3 0 1 0 0 6 3 3 0 0 0 0-6zM8.6 7.6 20 18M8.6 16.4 20 6',
  insert: 'M12 4v16M8 4h8M8 20h8',
  consensus: 'M4 7h16M4 12h10M4 17h16',
  globe: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM3 12h18M12 3a14 14 0 0 1 0 18M12 3a14 14 0 0 0 0 18',
  camera: 'M4 8h3l2-3h6l2 3h3a1 1 0 0 1 1 1v10a1 1 0 0 1-1 1H4a1 1 0 0 1-1-1V9a1 1 0 0 1 1-1zM12 10a3.5 3.5 0 1 0 0 7 3.5 3.5 0 0 0 0-7z',
  x: 'M6 6l12 12M18 6 6 18',
  check: 'm5 12 5 5 9-10',
  star: 'm12 3 2.7 5.6 6.1.9-4.4 4.3 1 6.1L12 17l-5.4 2.9 1-6.1-4.4-4.3 6.1-.9z',
  marker: 'M12 21s-7-6.2-7-11a7 7 0 0 1 14 0c0 4.8-7 11-7 11zM12 7.5a2.5 2.5 0 1 0 0 5 2.5 2.5 0 0 0 0-5z',
  strands: 'M4 8h13l-3-3M20 16H7l3 3',
  panelUp: 'm6 15 6-6 6 6',
  panelDown: 'm6 9 6 6 6-6',
  eyeOff: 'M3 3l18 18M10.6 6.1A9.8 9.8 0 0 1 12 6c6.5 0 10 6 10 6a17 17 0 0 1-3.1 3.8M6.6 6.6C3.7 8.4 2 12 2 12s3.5 6 10 6c1.7 0 3.2-.4 4.5-1',
};

/** A 24-unit line icon (lucide style), drawn in the current text colour. */
export function Icon({ name, size = 16, className = '' }: { name: keyof typeof ICON_PATHS | string; size?: number; className?: string }) {
  return (
    <svg width={size} height={size} viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} strokeLinecap="round" strokeLinejoin="round"
      aria-hidden="true" className={`shrink-0 ${className}`}>
      <path d={ICON_PATHS[name] ?? ''} />
    </svg>
  );
}

/** One choice among a few. `size="sm"` for the secondary ones (depth axis, raw / collapsed). */
export function Segmented<T extends string>({ value, onChange, options, disabled, title, size = 'md', label, prefix }: {
  value: T; onChange: (v: T) => void; disabled?: boolean; title: string; size?: 'sm' | 'md'; label?: string;
  /** a key in grey inside the switch, before its choices ("AF ≥", "window ≤"): one control instead of a label and a menu */
  prefix?: string;
  /** `dot`: a colour dot (a sample's track colour) before the label, or overlapping dots for a choice covering several */
  options: { value: T; label: string; icon?: string; hint?: string; dot?: string | string[] }[];
}) {
  const sm = size === 'sm';
  return (
    <span role="radiogroup" aria-label={label ?? title} title={title}
      className={`inline-flex items-center gap-0.5 rounded-[10px] bg-slate-50 border border-slate-200 p-[3px] select-none ${disabled ? 'opacity-60' : ''}`}>
      {prefix && <span className={`pl-1.5 pr-1 whitespace-nowrap text-slate-500 ${sm ? 'text-[11px]' : 'text-xs'}`}>{prefix}</span>}
      {options.map(o => {
        const active = o.value === value;
        return (
          <button key={o.value} type="button" role="radio" aria-checked={active} disabled={disabled} onClick={() => onChange(o.value)} title={o.hint}
            className={`inline-flex items-center gap-1.5 whitespace-nowrap rounded-[7px] font-medium transition-all disabled:cursor-not-allowed ${sm ? 'h-[22px] px-2 text-[11.5px]' : 'h-[26px] px-2.5 text-[12.5px]'} ${active ? 'bg-white text-indigo-600 shadow-[0_1px_2px_rgba(15,23,42,.12),0_0_0_1px_rgba(15,23,42,.04)]' : 'text-slate-500 hover:text-slate-800'}`}>
            {o.icon && <Icon name={o.icon} size={sm ? 13 : 15} />}
            {o.dot && <Dots colors={Array.isArray(o.dot) ? o.dot : [o.dot]} faded={!active} />}
            <span className="truncate max-w-[150px]">{o.label}</span>
          </button>
        );
      })}
    </span>
  );
}

/** One colour dot, or several overlapping (a choice covering several samples); dimmed when not chosen. */
export function Dots({ colors, faded }: { colors: string[]; faded?: boolean }) {
  return (
    <span className={`inline-flex items-center shrink-0 ${faded ? 'opacity-60' : ''}`} aria-hidden="true">
      {colors.slice(0, 4).map((c, i) => (
        <span key={i} className="w-2 h-2 rounded-full ring-[1.5px] ring-white" style={{ background: c, marginLeft: i ? -3 : 0 }} />
      ))}
    </span>
  );
}

/** An on/off option: a pill, tinted indigo when on. */
export function Pill({ on, onChange, label, title, icon, disabled, badge, size = 'md' }: {
  on: boolean; onChange: (v: boolean) => void; label: string; title: string; icon?: string; disabled?: boolean; badge?: string | number; size?: 'sm' | 'md';
}) {
  return (
    <button type="button" aria-pressed={on} disabled={disabled} title={title} onClick={() => onChange(!on)}
      className={`inline-flex items-center gap-1.5 whitespace-nowrap border font-medium transition-colors disabled:opacity-40 disabled:cursor-not-allowed select-none ${size === 'sm' ? 'h-[26px] px-2 rounded-lg text-xs' : 'h-[30px] px-2.5 rounded-[9px] text-[12.5px]'} ${on ? 'bg-indigo-50 border-indigo-200 text-indigo-600' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'}`}>
      {icon && <Icon name={icon} size={size === 'sm' ? 13 : 15} />}
      {label}
      {badge != null && badge !== 0 && <span className="text-[10.5px] font-bold leading-4 px-1.5 rounded-full bg-indigo-600 text-white">{badge}</span>}
    </button>
  );
}

/** A number with − / + buttons; the value is typed in place (committed on each keystroke, clamped to [min, max]). */
export function Stepper({ value, onChange, min, max, step = 1, label, unit, title, width = 34, after, ariaLabel }: {
  value: number; onChange: (v: number) => void; min: number; max?: number; step?: number; label?: string; unit?: string; title: string;
  width?: number; after?: React.ReactNode; ariaLabel?: string;
}) {
  const clamp = (v: number) => Math.min(max ?? Infinity, Math.max(min, v));
  const [text, setText] = useState(String(value));
  useEffect(() => { if (parseFloat(text) !== value) setText(String(value)); }, [value]); // eslint-disable-line react-hooks/exhaustive-deps
  const bump = (d: number) => onChange(clamp(Math.round((value + d) * 1000) / 1000));
  return (
    <span title={title} className="inline-flex items-center h-[30px] border border-slate-200 rounded-[9px] bg-white overflow-hidden text-[12.5px]">
      {label && <span className="pl-2.5 pr-1.5 text-slate-500 text-xs whitespace-nowrap">{label}</span>}
      <input type="number" value={text} min={min} max={max} step={step} aria-label={ariaLabel ?? label ?? title}
        onChange={e => { setText(e.target.value); const v = parseFloat(e.target.value); if (Number.isFinite(v)) onChange(clamp(v)); }}
        onBlur={() => setText(String(value))}
        className="text-center font-semibold text-slate-900 bg-transparent outline-none focus:bg-indigo-50 [appearance:textfield] [&::-webkit-inner-spin-button]:appearance-none [&::-webkit-outer-spin-button]:appearance-none"
        style={{ width }} />
      {unit && <span className="pr-1.5 text-slate-500 text-xs whitespace-nowrap">{unit}</span>}
      {after}
      <button type="button" onClick={() => bump(-step)} disabled={value <= min} aria-label="decrease" className="w-6 h-full grid place-items-center text-slate-400 hover:text-slate-800 hover:bg-slate-50 border-l border-slate-100 disabled:opacity-40"><Icon name="minus" size={13} /></button>
      <button type="button" onClick={() => bump(step)} disabled={max != null && value >= max} aria-label="increase" className="w-6 h-full grid place-items-center text-slate-400 hover:text-slate-800 hover:bg-slate-50 border-l border-slate-100 disabled:opacity-40"><Icon name="plus" size={13} /></button>
    </span>
  );
}

/** A native select dressed as a pill: the key in grey, the value in ink, a chevron. */
export function Select({ label, value, onChange, children, title, ariaLabel, dot }: {
  label?: string; value: string | number; onChange: (v: string) => void; children: React.ReactNode; title: string; ariaLabel?: string; dot?: string;
}) {
  return (
    <label title={title} className="relative inline-flex items-center gap-1.5 h-[30px] pl-2.5 pr-7 border border-slate-200 rounded-[9px] bg-white text-[12.5px] text-slate-700 whitespace-nowrap hover:border-slate-300 cursor-pointer">
      {dot && <span className="w-2 h-2 rounded-full shrink-0" style={{ background: dot }} />}
      {label && <span className="text-slate-500 text-xs">{label}</span>}
      <select value={value} onChange={e => onChange(e.target.value)} aria-label={ariaLabel ?? label ?? title}
        className="appearance-none bg-transparent outline-none font-medium text-slate-800 cursor-pointer pr-0.5 [field-sizing:content]">
        {children}
      </select>
      <Icon name="chev" size={13} className="absolute right-2 text-slate-400 pointer-events-none" />
    </label>
  );
}

/**
 * The colour of a toolbar section: a pastel panel and a label in its hue, so the categories read apart at a glance
 * (literal class names: Tailwind keeps only those it sees).
 */
const TONES = {
  indigo: { panel: 'bg-indigo-50 border-indigo-200', label: 'text-indigo-600' },
  emerald: { panel: 'bg-emerald-50 border-emerald-200', label: 'text-emerald-700' },
  orange: { panel: 'bg-orange-50 border-orange-200', label: 'text-orange-700' },
  sky: { panel: 'bg-sky-50 border-sky-200', label: 'text-sky-700' },
} as const;
export type Tone = keyof typeof TONES;

/**
 * A group of the toolbar: its uppercase label and its controls on a pastel panel of its tone. A section without a
 * label (Filters) stands apart on the right of its line.
 */
export function Section({ label, title, tone, children }: { label?: string; title?: string; tone?: Tone; children: React.ReactNode }) {
  const t = tone ? TONES[tone] : null;
  return (
    <div className={`flex flex-wrap items-center gap-1.5 ${t ? `rounded-xl border px-2 py-1 ${t.panel}` : 'ml-auto'}`} role="group" aria-label={label}>
      {label && <span className={`text-[10px] font-bold tracking-[.08em] uppercase mr-0.5 pl-0.5 select-none ${t ? t.label : 'text-slate-400'}`} title={title}>{label}</span>}
      {children}
    </div>
  );
}

/** The toolbar: its sections side by side, wrapping on narrow windows. */
export function Bar({ children, className = '' }: { children: React.ReactNode; className?: string }) {
  return <div className={`flex flex-wrap items-center gap-x-2.5 gap-y-2 ${className}`}>{children}</div>;
}

/** The thin rule between two sections. */
export function Sep() {
  return <span className="w-px h-[26px] bg-slate-200 shrink-0" aria-hidden="true" />;
}

/**
 * A button and the panel it opens under itself. The panel closes on a click outside or Escape; it opens towards the
 * side that keeps it inside the window.
 */
export function Popover({ button, children, width = 320, title, open: openProp, onOpenChange, buttonClass, label, flush }: {
  button: React.ReactNode; children: React.ReactNode | ((close: () => void) => React.ReactNode); width?: number; title: string;
  open?: boolean; onOpenChange?: (v: boolean) => void; buttonClass?: string; label?: string;
  /** no padding: the panel holds a search box and a list edge to edge */
  flush?: boolean;
}) {
  const [own, setOwn] = useState(false);
  const open = openProp ?? own;
  const set = (v: boolean) => { setOwn(v); onOpenChange?.(v); };
  const ref = useRef<HTMLSpanElement>(null);
  const [right, setRight] = useState(false);
  useEffect(() => {
    if (!open) return;
    const down = (e: MouseEvent) => { if (ref.current && !ref.current.contains(e.target as Node)) set(false); };
    const key = (e: KeyboardEvent) => { if (e.key === 'Escape') set(false); };
    document.addEventListener('mousedown', down); document.addEventListener('keydown', key);
    return () => { document.removeEventListener('mousedown', down); document.removeEventListener('keydown', key); };
  }); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <span ref={ref} className="relative inline-flex">
      <button type="button" title={title} aria-haspopup="true" aria-expanded={open} aria-label={label}
        onClick={e => { const r = e.currentTarget.getBoundingClientRect(); setRight(r.left + width + 16 > window.innerWidth); set(!open); }}
        className={buttonClass ?? `inline-flex items-center gap-1.5 h-[30px] px-2.5 rounded-[9px] border text-[12.5px] font-medium whitespace-nowrap transition-colors ${open ? 'bg-indigo-50 border-indigo-200 text-indigo-600' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'}`}>
        {button}
      </button>
      {open && (
        <div role="dialog" className={`absolute top-full mt-2 ${right ? 'right-0' : 'left-0'} z-30 bg-white border border-slate-200 rounded-[14px] ${flush ? 'overflow-hidden' : 'p-2'} shadow-[0_18px_40px_-12px_rgba(15,23,42,.28),0_2px_6px_rgba(15,23,42,.08)] text-left`}
          style={{ width }}>
          {typeof children === 'function' ? children(() => set(false)) : children}
        </div>
      )}
    </span>
  );
}

/** A line of a menu: icon, label, optional hint under it. */
export function MenuItem({ icon, label, hint, onClick, disabled, title }: { icon?: string; label: string; hint?: string; onClick: () => void; disabled?: boolean; title?: string }) {
  return (
    <button type="button" onClick={onClick} disabled={disabled} title={title}
      className="w-full flex items-start gap-2.5 px-2.5 py-2 rounded-lg text-left hover:bg-slate-50 disabled:opacity-40 disabled:hover:bg-transparent">
      {icon && <span className="mt-px text-slate-500"><Icon name={icon} size={15} /></span>}
      <span className="flex-1 min-w-0">
        <span className="block text-[12.5px] font-medium text-slate-800">{label}</span>
        {hint && <span className="block text-[11.5px] text-slate-500 leading-snug mt-px">{hint}</span>}
      </span>
    </button>
  );
}

/** An option of a popover: its name and what it does, and a switch. */
export function SwitchRow({ on, onChange, label, hint, title, disabled }: { on: boolean; onChange: (v: boolean) => void; label: string; hint: string; title?: string; disabled?: boolean }) {
  return (
    <button type="button" role="switch" aria-checked={on} disabled={disabled} onClick={() => onChange(!on)} title={title}
      className="w-full flex items-center gap-3 px-2.5 py-2 rounded-lg text-left hover:bg-slate-50 disabled:opacity-45">
      <span className="flex-1 min-w-0">
        <span className="block text-[12.5px] font-semibold text-slate-900">{label}</span>
        <span className="block text-[11.5px] text-slate-500 leading-snug mt-px">{hint}</span>
      </span>
      <span className={`relative w-[30px] h-[18px] rounded-full shrink-0 transition-colors ${on ? 'bg-indigo-600' : 'bg-slate-300'}`}>
        <span className={`absolute top-[2px] w-[14px] h-[14px] rounded-full bg-white shadow transition-all ${on ? 'left-[14px]' : 'left-[2px]'}`} />
      </span>
    </button>
  );
}
