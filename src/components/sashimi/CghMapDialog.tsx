/**
 * The column mapping of an array-CGH table whose columns were not recognised (or are to be changed): delimiter and
 * header row, a preview, one select per field, what the Type values mean. The choice is remembered by the header.
 */
import { useMemo, useState } from 'react';
import { cgh, type CghState } from './cghStore';
import { parseCghFile, cghBuild, kindOfText, type CghFile } from './cghIntervals';
import { detectTable, guessCols, guessSizeUnit, rememberMap, splitDelimited, type ColMap, type Field, type MapOpts, type MapSpec, type TypeChoice } from './cghMapping';

const LABEL: Record<Field, string> = { chrom: 'Chromosome', start: 'Start', stop: 'End', type: 'Type (del / dup / loss / gain…)', copies: 'Copy number', iscn: 'ISCN', lr: 'Log ratio', probes: 'Probes', size: 'Size', band: 'Cytoband', genes: 'Genes', cls: 'Classification', note: 'Comment', sample: 'Sample name', sup: 'Suppress' };
const REQUIRED: Field[] = ['chrom', 'start', 'stop'], ONE_OF: Field[] = ['type', 'copies', 'iscn', 'lr'];
const OPTIONAL: Field[] = ['probes', 'size', 'band', 'genes', 'cls', 'note', 'sample', 'sup'];
const DELIMS: [string, string][] = [['\t', 'tab'], [';', 'semicolon ;'], [',', 'comma ,']];
const SEL = 'h-[26px] px-1.5 rounded border border-slate-200 bg-white text-[11.5px] outline-none cursor-pointer max-w-[210px]';

export function CghMapDialog({ item }: { item: CghState['toMap'][number] }) {
  const lines = useMemo(() => item.text.replace(/^﻿/, '').split(/\r?\n/), [item.text]);
  const first = useMemo(() => {
    const t = detectTable(lines);
    return item.spec?.opts ?? { delim: t?.delim ?? '\t', headerRow: t?.headerRow ?? Math.max(0, lines.findIndex(l => l.trim())), zeroBased: false, sizeUnit: 'kb', build: 'GRCh37', typeMap: {} } as MapOpts;
  }, [lines, item.spec]);
  const [delim, setDelim] = useState(first.delim), [headerRow, setHeaderRow] = useState(first.headerRow);
  const [cols, setCols] = useState<ColMap>(() => item.spec?.cols ?? guessCols(splitDelimited(lines[first.headerRow] ?? '', first.delim)));
  const [zeroBased, setZeroBased] = useState(first.zeroBased), [sizeUnit, setSizeUnit] = useState(first.sizeUnit);
  const [build, setBuild] = useState(first.build ?? 'GRCh37'), [typeMap, setTypeMap] = useState(first.typeMap);

  const header = splitDelimited(lines[headerRow] ?? '', delim);
  const data = useMemo(() => lines.slice(headerRow + 1).filter(l => l.trim()).map(l => splitDelimited(l, delim)), [lines, headerRow, delim]);
  const statesBuild = lines.slice(0, headerRow).some(l => cghBuild(l));
  // the Type values that say neither loss nor gain: each needs a choice (ignored until one is made)
  const odd = useMemo(() => {
    if (cols.type < 0) return [];
    const seen = new Map<string, string>();
    for (const r of data) { const v = (r[cols.type] ?? '').trim(); if (v && !kindOfText(v) && !seen.has(v.toLowerCase())) seen.set(v.toLowerCase(), v); }
    return [...seen.values()].slice(0, 30);
  }, [data, cols.type]);
  const spec: MapSpec = useMemo(() => ({ cols, opts: { delim, headerRow, zeroBased, sizeUnit, build: statesBuild ? null : build, typeMap: Object.fromEntries(odd.map(v => [v.toLowerCase(), typeMap[v.toLowerCase()] ?? 'ignore'])) } }), [cols, delim, headerRow, zeroBased, sizeUnit, build, statesBuild, odd, typeMap]);
  const ready = REQUIRED.every(f => cols[f] >= 0) && ONE_OF.some(f => cols[f] >= 0);
  const res = useMemo((): { file: CghFile | null; msg: string } => {
    if (!ready) return { file: null, msg: 'Choose Chromosome, Start, End and one of Type, Copy number, ISCN, Log ratio.' };
    try {
      const f = parseCghFile(item.text, item.name, item.fallback, spec), n = f.samples.reduce((k, x) => k + x.cnvs.length, 0);
      return { file: f, msg: `${n} row${n === 1 ? '' : 's'} read, ${f.skipped} skipped${f.firstSkip ? ` (first: ${f.firstSkip})` : ''}` };
    } catch (e) { return { file: null, msg: `0 rows read: ${(e as Error).message.replace(`${item.name}: `, '')}` }; }
  }, [ready, spec, item]);

  const guess = (d: string, h: number) => { setCols(guessCols(splitDelimited(lines[h] ?? '', d))); setTypeMap({}); };
  const setCol = (f: Field, v: number) => {
    const next = { ...cols, [f]: v };
    setCols(next);
    if ((f === 'size' || f === 'start' || f === 'stop') && next.size >= 0 && next.start >= 0 && next.stop >= 0) setSizeUnit(guessSizeUnit(data, next, delim !== ','));
  };
  const field = (f: Field, required: boolean) => <label key={f} className="flex items-center gap-2 py-0.5">
    <span className={`w-[170px] shrink-0 ${required ? 'font-semibold' : ''}`}>{LABEL[f]}{required ? ' *' : ''}</span>
    <select value={cols[f]} onChange={e => setCol(f, Number(e.target.value))} aria-label={LABEL[f]} className={SEL}>
      <option value={-1}>— none</option>
      {header.map((h, i) => <option key={i} value={i}>{h.trim() || `column ${i + 1}`}</option>)}
    </select>
    {f === 'size' && cols.size >= 0 && <select value={sizeUnit} onChange={e => setSizeUnit(e.target.value as MapOpts['sizeUnit'])} aria-label="Size unit" className={SEL}>{['bp', 'kb', 'Mb'].map(u => <option key={u}>{u}</option>)}</select>}
  </label>;
  const cancel = () => cgh.mapped(null);
  const use = () => { if (!res.file) return; rememberMap(header, spec); cgh.mapped(res.file); };
  return (
    <div className="fixed inset-0 z-[80] bg-black/40 flex items-start justify-center p-6 overflow-y-auto" role="dialog" aria-modal="true" aria-label="Column mapping"
      onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); cancel(); } }}>
      <div className="bg-white rounded-xl shadow-2xl border border-slate-200 w-full max-w-3xl text-slate-900 text-xs mt-6 outline-none" tabIndex={-1} ref={el => { if (el && !el.contains(document.activeElement)) el.focus(); }}>
        <div className="px-4 pt-3 pb-2 border-b border-slate-200">
          <div className="font-bold text-sm">Which column is which?</div>
          <div className="text-[11px] text-slate-500 mt-0.5 truncate" title={item.name}>{item.name}{item.spec ? '' : ': its columns were not recognised.'} The choice is remembered for tables with the same column names.</div>
        </div>
        <div className="px-4 py-3 space-y-3">
          <div className="flex flex-wrap items-center gap-x-4 gap-y-1">
            <label className="flex items-center gap-1.5"><span className="font-semibold">Delimiter</span>
              <select value={delim} onChange={e => { setDelim(e.target.value); guess(e.target.value, headerRow); }} aria-label="Delimiter" className={SEL}>{DELIMS.map(([v, t]) => <option key={v} value={v}>{t}</option>)}</select></label>
            <label className="flex items-center gap-1.5"><span className="font-semibold">Header row</span>
              <select value={headerRow} onChange={e => { const h = Number(e.target.value); setHeaderRow(h); guess(delim, h); }} aria-label="Header row" className={SEL}>
                {lines.slice(0, 30).map((l, i) => l.trim() ? <option key={i} value={i}>line {i + 1}: {l.slice(0, 28)}</option> : null)}</select></label>
          </div>
          <div className="overflow-x-auto border border-slate-200 rounded" data-map-preview>
            <table className="text-[11px] whitespace-nowrap">
              <thead><tr className="bg-slate-50 text-left">{header.map((h, i) => <th key={i} className="px-2 py-1 font-semibold border-b border-slate-200">{h.trim() || `column ${i + 1}`}</th>)}</tr></thead>
              <tbody>{data.slice(0, 8).map((r, j) => <tr key={j} className="border-b border-slate-100">{header.map((_, i) => <td key={i} className="px-2 py-0.5 text-slate-600 max-w-[180px] truncate">{r[i] ?? ''}</td>)}</tr>)}</tbody>
            </table>
          </div>
          <div className="grid sm:grid-cols-2 gap-x-6">
            <div>{REQUIRED.map(f => field(f, true))}
              <div className="text-[11px] text-slate-500 mt-1.5">At least one of these (in this order of preference):</div>{ONE_OF.map(f => field(f, false))}</div>
            <div><div className="text-[11px] text-slate-500">Optional</div>{OPTIONAL.map(f => field(f, false))}</div>
          </div>
          <div className="flex flex-wrap items-center gap-x-5 gap-y-1">
            <label className="flex items-center gap-1.5 cursor-pointer"><input type="checkbox" checked={zeroBased} onChange={e => setZeroBased(e.target.checked)} />Start is 0-based (BED-like)</label>
            {!statesBuild && <label className="flex items-center gap-1.5"><span className="font-semibold">Genome build</span>
              <select value={build} onChange={e => setBuild(e.target.value as 'GRCh37' | 'GRCh38')} aria-label="Genome build" className={SEL}><option>GRCh37</option><option>GRCh38</option></select><span className="text-slate-500">the file states none</span></label>}
          </div>
          {odd.length > 0 && <div className="rounded-lg border border-amber-300 bg-amber-50 px-3 py-2">
            <div className="font-semibold text-amber-900 mb-1">Type values that say neither loss nor gain</div>
            <div className="flex flex-wrap gap-x-4 gap-y-1">{odd.map(v => <label key={v} className="flex items-center gap-1.5"><span className="font-mono">{v}</span>
              <select value={typeMap[v.toLowerCase()] ?? 'ignore'} onChange={e => setTypeMap(m => ({ ...m, [v.toLowerCase()]: e.target.value as TypeChoice }))} aria-label={`Meaning of ${v}`} className={SEL}><option>ignore</option><option>loss</option><option>gain</option></select></label>)}</div>
          </div>}
          <div className={`text-[11.5px] ${res.file ? 'text-slate-700' : 'text-red-700'}`} data-map-status>{res.msg}</div>
        </div>
        <div className="px-4 py-2.5 border-t border-slate-200 flex justify-end gap-2">
          <button onClick={cancel} className="h-[30px] px-3 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50">Cancel</button>
          <button disabled={!res.file || res.file.samples.length === 0} onClick={use} className="h-[30px] px-4 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700 disabled:opacity-40 disabled:cursor-not-allowed">Use these columns</button>
        </div>
      </div>
    </div>
  );
}
