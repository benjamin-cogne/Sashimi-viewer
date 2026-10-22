/**
 * A search input that suggests gene names from the third letter typed (geneSuggest.ts): the genes the page knows first,
 * then HGNC's. Arrow keys move through the list, Enter or a click opens the gene (the input's form is submitted with the
 * symbol), Escape closes it. Text that does not read as a gene (coordinates, c. positions, exon numbers) gets no list.
 */
import { useEffect, useRef, useState } from 'react';
import { SUGGEST_MIN, looksLikeGene, suggestGenes, type GeneSuggestion } from '../../standalone/geneSuggest';

export function GeneSuggest({ value, onChange, local = [], className, placeholder, title, ariaLabel, width = 340 }: {
  value: string; onChange: (v: string) => void; local?: string[]; className?: string; placeholder?: string; title?: string; ariaLabel?: string;
  /** width of the list, px */
  width?: number;
}) {
  const [list, setList] = useState<GeneSuggestion[]>([]);
  const [open, setOpen] = useState(false);
  const [active, setActive] = useState(-1);
  const inputRef = useRef<HTMLInputElement>(null);
  /** the text was just set by a pick: no list for it */
  const picked = useRef<string | null>(null);
  const localKey = local.join('|');

  useEffect(() => {
    const q = value.trim();
    if (!looksLikeGene(q) || q.length < SUGGEST_MIN || picked.current === value) { setList([]); return; }
    const ctl = new AbortController();
    const timer = setTimeout(() => {
      suggestGenes(q, local, ctl.signal).then(l => { if (!ctl.signal.aborted) { setList(l); setActive(-1); } }).catch(() => {});
    }, 180);
    return () => { clearTimeout(timer); ctl.abort(); };
  }, [value, localKey]); // eslint-disable-line react-hooks/exhaustive-deps

  const pick = (s: GeneSuggestion) => {
    picked.current = s.symbol;
    onChange(s.symbol);
    setOpen(false); setList([]);
    // the form is submitted once the new text is in it (the handler reads it from the state)
    const form = inputRef.current?.form;
    setTimeout(() => form?.requestSubmit(), 0);
  };
  const shown = open && list.length > 0;

  return (
    <span className="relative flex-1 min-w-0 flex">
      <input ref={inputRef} value={value} placeholder={placeholder} title={title} aria-label={ariaLabel} className={className}
        role="combobox" aria-expanded={shown} aria-autocomplete="list" autoComplete="off" spellCheck={false}
        onChange={e => { picked.current = null; onChange(e.target.value); setOpen(true); }}
        onFocus={() => setOpen(true)}
        onBlur={() => setTimeout(() => setOpen(false), 120)}
        onKeyDown={e => {
          if (!shown) return;
          if (e.key === 'ArrowDown') { e.preventDefault(); setActive(a => (a + 1) % list.length); }
          else if (e.key === 'ArrowUp') { e.preventDefault(); setActive(a => (a <= 0 ? list.length - 1 : a - 1)); }
          else if (e.key === 'Escape') { setOpen(false); }
          else if (e.key === 'Enter' && active >= 0) { e.preventDefault(); pick(list[active]); }
        }} />
      {shown && (
        <div role="listbox" className="absolute left-0 top-full mt-2 z-40 bg-white border border-slate-200 rounded-xl py-1 shadow-[0_18px_40px_-12px_rgba(15,23,42,.28),0_2px_6px_rgba(15,23,42,.08)] text-left"
          style={{ width }}>
          {list.map((s, i) => (
            <button key={s.symbol} type="button" role="option" aria-selected={i === active}
              onMouseDown={e => e.preventDefault()} onClick={() => pick(s)} onMouseEnter={() => setActive(i)}
              className={`w-full flex items-baseline gap-2 px-3 py-1.5 text-left ${i === active ? 'bg-indigo-50' : 'hover:bg-slate-50'}`}>
              <span className="font-semibold text-[12.5px] text-slate-900 shrink-0">{s.symbol}</span>
              <span className="flex-1 min-w-0 truncate text-[11.5px] text-slate-500">
                {s.via ? <span className="text-amber-700">{s.via.kind === 'previous' ? 'previously' : 'alias'} {s.via.text} · </span> : null}
                {s.name ?? (s.local ? 'opened in this page before' : '')}
              </span>
              {s.location && <span className="shrink-0 text-[11px] font-mono text-slate-400">{s.location}</span>}
            </button>
          ))}
          <div className="px-3 pt-1 pb-0.5 text-[10.5px] text-slate-400 border-t border-slate-100 mt-1">↑ ↓ to choose · Enter to open · gene names from HGNC</div>
        </div>
      )}
    </span>
  );
}
