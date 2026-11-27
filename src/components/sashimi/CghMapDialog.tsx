/**
 * The column mapping of an array-CGH table whose columns were not recognised (or are to be changed): first what the table
 * is (calls, or probe log ratios: one row per probe), then delimiter and header row (or none: the columns are then named
 * column 1, 2…), a preview, one select per field, what the Type values mean. In probe mode: chromosome, position (or
 * start, + end), one log-ratio column per sample (named), probe name; the file is then read as a probe file (in the
 * worker), one patient per column. The choice is remembered by the header (a table without one is not remembered).
 */
import { useMemo, useState } from 'react';
import { PROBE_MIN_ROWS, guessProbeCols, rememberProbeMap, sampleOfColumn, type ProbeSpec } from './cghProbes';
import { cghChrom } from './cghIntervals';
import { cgh, type CghState } from './cghStore';
import { parseCghFile, cghBuild, kindOfText, type CghFile } from './cghIntervals';
import { detectTable, guessCols, guessLrMore, guessSizeUnit, looksLikeHeader, rememberMap, splitDelimited, splitRecords, type ColMap, type Field, type MapOpts, type MapSpec, type TypeChoice } from './cghMapping';

const LABEL: Record<Field, string> = { chrom: 'Chromosome', start: 'Start', stop: 'End', pos: 'Position (one column)', type: 'Type (del / dup / loss / gain…)', copies: 'Copy number', iscn: 'ISCN', lr: 'Log ratio (one or several)', probes: 'Probes', size: 'Size', band: 'Cytoband', genes: 'Genes', cls: 'Classification', note: 'Comment', sample: 'Sample name', sup: 'Suppress' };
const PLACE: Field[] = ['chrom', 'start', 'stop'], ONE_OF: Field[] = ['type', 'lr', 'copies', 'iscn'];
const OPTIONAL: Field[] = ['probes', 'size', 'band', 'genes', 'cls', 'note', 'sample', 'sup'];
const DELIMS: [string, string][] = [['\t', 'tab'], [';', 'semicolon ;'], [',', 'comma ,']];
const SEL = 'h-[26px] px-1.5 rounded border border-slate-200 bg-white text-[11.5px] outline-none cursor-pointer max-w-[210px]';

export function CghMapDialog({ item }: { item: CghState['toMap'][number] }) {
  const lines = useMemo(() => splitRecords(item.text), [item.text]);
  // a long table is likely a probe file: asked as one first (the dialog only reads its first rows then)
  const [mode, setMode] = useState<'calls' | 'probes'>(() => (!item.spec && lines.length > PROBE_MIN_ROWS ? 'probes' : 'calls'));
  /** a row's cells, and the header's (generated names when the table has no header row) */
  const rowsFrom = (d: string, h: number) => lines.slice(h + 1, mode === 'probes' ? h + 301 : undefined).filter(l => l.trim()).map(l => splitDelimited(l, d));
  const headerOf = (d: string, h: number) => (h >= 0 ? splitDelimited(lines[h] ?? '', d) : Array.from({ length: Math.max(0, ...rowsFrom(d, h).slice(0, 40).map(r => r.length)) }, () => ''));
  const first = useMemo(() => {
    const t = detectTable(lines);
    return item.spec?.opts ?? { delim: t?.delim ?? '\t', headerRow: t?.headerRow ?? Math.max(0, lines.findIndex(l => l.trim())), zeroBased: false, sizeUnit: 'kb', build: 'GRCh37', typeMap: {} } as MapOpts;
  }, [lines, item.spec]);
  const [delim, setDelim] = useState(first.delim), [headerRow, setHeaderRow] = useState(first.headerRow);
  /** the columns guessed by name and content, with every log-ratio column (gain / loss ones) */
  const guessAll = (d: string, h: number) => { const hd = headerOf(d, h), rs = rowsFrom(d, h); return guessLrMore(hd, rs, guessCols(hd, rs)); };
  const [init] = useState(() => (item.spec ? { cols: item.spec.cols, lrMore: item.spec.lrMore ?? [] } : guessAll(first.delim, first.headerRow)));
  const [cols, setCols] = useState<ColMap>(init.cols);
  // further log-ratio columns, merged with the first: the non-zero value of a row
  const [lrMore, setLrMore] = useState<number[]>(init.lrMore);
  const [zeroBased, setZeroBased] = useState(first.zeroBased), [sizeUnit, setSizeUnit] = useState(first.sizeUnit);
  const [build, setBuild] = useState(first.build ?? 'GRCh37'), [typeMap, setTypeMap] = useState(first.typeMap);

  const data = useMemo(() => rowsFrom(delim, headerRow), [lines, headerRow, delim, mode]); // eslint-disable-line react-hooks/exhaustive-deps
  const header = useMemo(() => headerOf(delim, headerRow), [lines, headerRow, delim]); // eslint-disable-line react-hooks/exhaustive-deps
  const statesBuild = lines.slice(0, Math.max(0, headerRow)).some(l => cghBuild(l));
  // the Type values that say neither loss nor gain: each needs a choice (ignored until one is made)
  const odd = useMemo(() => {
    if (cols.type < 0) return [];
    const seen = new Map<string, string>();
    for (const r of data) { const v = (r[cols.type] ?? '').trim(); if (v && !kindOfText(v) && !seen.has(v.toLowerCase())) seen.set(v.toLowerCase(), v); }
    return [...seen.values()].slice(0, 30);
  }, [data, cols.type]);
  const spec: MapSpec = useMemo(() => ({ cols, lrMore: lrMore.filter(j => j >= 0 && j !== cols.lr), opts: { delim, headerRow, zeroBased, sizeUnit, build: statesBuild ? null : build, ignoreOtherTypes: true,
    // only the values given a meaning: an ignored one is the default, and is not remembered
    typeMap: Object.fromEntries(odd.flatMap(v => { const t = typeMap[v.toLowerCase()]; return t && t !== 'ignore' ? [[v.toLowerCase(), t]] : []; })) } }), [cols, lrMore, delim, headerRow, zeroBased, sizeUnit, build, statesBuild, odd, typeMap]);
  // the place: Chromosome, Start and End, or one column holding it, or an ISCN (it holds the place and the copies)
  const placeOk = PLACE.every(f => cols[f] >= 0) || cols.pos >= 0 || cols.iscn >= 0;
  const typeOk = ONE_OF.some(f => cols[f] >= 0) || lrMore.some(j => j >= 0);
  // what is missing, said by name
  const missing = [!placeOk && `the place (${PLACE.filter(f => cols[f] < 0).map(f => LABEL[f]).join(', ')}, or a Position in one column)`,
    !typeOk && 'what the call is: a Type, a Log ratio or a Copy number column (or an ISCN)'].filter(Boolean);
  const ready = placeOk && typeOk;
  const res = useMemo((): { file: CghFile | null; msg: string } => {
    if (mode === 'probes') return { file: null, msg: '' };
    if (!ready) return { file: null, msg: `Missing: ${missing.join('; and ')}.` };
    try {
      const f = parseCghFile(item.text, item.name, item.fallback, spec), n = f.samples.reduce((k, x) => k + x.cnvs.length, 0);
      return { file: f, msg: `${n} row${n === 1 ? '' : 's'} read, ${f.skipped} skipped${f.firstSkip ? ` (first: ${f.firstSkip})` : ''}` };
    } catch (e) { return { file: null, msg: `0 rows read: ${(e as Error).message.replace(`${item.name}: `, '')}` }; }
  }, [ready, spec, item, mode]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- probe mode: chromosome, position (or start, + end), a log-ratio column per sample (named), probe name ----
  const [pc, setPc] = useState(() => guessProbeCols(header, data));
  const [names, setNames] = useState<string[]>([]);
  const fileBuild = cghBuild(lines.slice(0, Math.max(0, headerRow)).join(' '));
  const lrCols = pc.lr.filter(j => j >= 0);
  const probeSpecs: ProbeSpec[] = lrCols.map((j, k) => ({ delim, headerRow, chrom: pc.chrom, start: pc.start, end: pc.end, lr: j, name: pc.name, zeroBased, build: fileBuild ?? build, layout: 'columns chosen',
    sample: names[pc.lr.indexOf(j)]?.trim() || sampleOfColumn(header[j] ?? '', lrCols.length > 1, item.name) || `sample ${k + 1}` }));
  const probeMissing = [pc.chrom < 0 && 'Chromosome', pc.start < 0 && 'Position or start', !lrCols.length && 'a Log ratio column'].filter(Boolean);
  const dup = new Set(probeSpecs.map(x => x.sample)).size < probeSpecs.length;
  const probeRead = probeMissing.length ? 0 : data.filter(r => { const c = cghChrom((r[pc.chrom] ?? '').replace(/"/g, '')); return c && c !== 'chrM' && /^\d+$/.test((r[pc.start] ?? '').trim().replace(/"/g, '')); }).length;
  const probeMsg = probeMissing.length ? `Missing: ${probeMissing.join(', ')}.` : dup ? 'Two log-ratio columns have the same sample name: name them.'
    : `${probeRead} of the first ${data.length} rows read as probes (the whole file, ${Math.max(0, lines.length - headerRow - 1).toLocaleString('en')} lines, is read after); ${probeSpecs.length} patient${probeSpecs.length > 1 ? 's' : ''}: ${probeSpecs.map(x => x.sample).join(', ')}; positions read as ${zeroBased ? '0-based, moved to 1-based' : '1-based'}`;
  const probeOk = !probeMissing.length && !dup && probeRead > 0;
  const setP = (f: 'chrom' | 'start' | 'end' | 'name', v: number) => setPc(p => ({ ...p, [f]: v }));
  const pSel = (label: string, v: number, on: (v: number) => void, required: boolean) => <label key={label} className="flex items-center gap-2 py-0.5">
    <span className={`w-[170px] shrink-0 ${required ? 'font-semibold' : ''}`}>{label}{required ? ' *' : ''}</span>
    <select value={v} onChange={e => on(Number(e.target.value))} aria-label={label} className={SEL}>
      <option value={-1}>— none</option>{header.map((h, i) => <option key={i} value={i}>{h.trim() || `column ${i + 1}`}</option>)}</select></label>;

  const guess = (d: string, h: number) => { const g = guessAll(d, h); setCols(g.cols); setLrMore(g.lrMore); setTypeMap({}); setPc(guessProbeCols(headerOf(d, h), rowsFrom(d, h))); setNames([]); };
  const setCol = (f: Field, v: number) => {
    const next = { ...cols, [f]: v };
    setCols(next);
    if ((f === 'size' || f === 'start' || f === 'stop' || f === 'pos') && next.size >= 0 && ((next.start >= 0 && next.stop >= 0) || next.pos >= 0)) setSizeUnit(guessSizeUnit(data, next, delim !== ','));
  };
  const field = (f: Field, required: boolean) => <label key={f} className="flex items-center gap-2 py-0.5">
    <span className={`w-[170px] shrink-0 ${required ? 'font-semibold' : ''}`}>{LABEL[f]}{required ? ' *' : ''}</span>
    <select value={cols[f]} onChange={e => setCol(f, Number(e.target.value))} aria-label={LABEL[f]} className={SEL}>
      <option value={-1}>— none</option>
      {header.map((h, i) => <option key={i} value={i}>{h.trim() || `column ${i + 1}`}</option>)}
    </select>
    {f === 'lr' && cols.lr >= 0 && <button type="button" className="text-accent-700 hover:underline whitespace-nowrap" onClick={() => setLrMore(m => [...m, -1])} title="Gain and loss log ratios in separate columns (as Amplification / Gain / Loss / Deletion): a row takes the non-zero one">+ another</button>}
    {f === 'size' && cols.size >= 0 && <select value={sizeUnit} onChange={e => setSizeUnit(e.target.value as MapOpts['sizeUnit'])} aria-label="Size unit" className={SEL}>{['bp', 'kb', 'Mb'].map(u => <option key={u}>{u}</option>)}</select>}
  </label>;
  const cancel = () => cgh.mapped(null);
  const headerOk = headerRow >= 0 && looksLikeHeader(header);
  const use = () => {
    if (mode === 'probes') {
      if (!probeOk) return;
      // the columns only (never a row, nor a sample name typed)
      if (headerOk) rememberProbeMap(header, { delim, chrom: pc.chrom, start: pc.start, end: pc.end, lr: lrCols, name: pc.name, zeroBased, build: fileBuild ? null : build });
      cgh.mappedProbes(probeSpecs); return;
    }
    if (!res.file) return; if (headerOk) rememberMap(header, spec); cgh.mapped(res.file);
  };
  return (
    <div className="fixed inset-0 z-[80] bg-black/40 flex items-start justify-center p-6 overflow-y-auto" role="dialog" aria-modal="true" aria-label="Column mapping"
      onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); cancel(); } }}>
      <div className="bg-white rounded-xl shadow-2xl border border-slate-200 w-full max-w-3xl text-slate-900 text-xs mt-6 outline-none" tabIndex={-1} ref={el => { if (el && !el.contains(document.activeElement)) el.focus(); }}>
        <div className="px-4 pt-3 pb-2 border-b border-slate-200">
          <div className="font-bold text-sm">Which column is which?</div>
          <div className="text-[11px] text-slate-500 mt-0.5 truncate" title={item.name}>{item.name}{item.spec ? '' : ': its columns were not recognised.'} The choice is remembered for tables with the same column names.</div>
        </div>
        <div className="px-4 py-3 space-y-3">
          <div className="flex flex-wrap items-center gap-x-4 gap-y-1" data-map-mode={mode}><span className="font-semibold">This table holds</span>
            <label className="flex items-center gap-1.5 cursor-pointer"><input type="radio" name="cgh-map-mode" checked={mode === 'calls'} onChange={() => setMode('calls')} />Calls (one row per call)</label>
            <label className="flex items-center gap-1.5 cursor-pointer"><input type="radio" name="cgh-map-mode" checked={mode === 'probes'} onChange={() => setMode('probes')} />Probe log ratios (one row per probe)</label>
          </div>
          <div className="flex flex-wrap items-center gap-x-4 gap-y-1">
            <label className="flex items-center gap-1.5"><span className="font-semibold">Delimiter</span>
              <select value={delim} onChange={e => { setDelim(e.target.value); guess(e.target.value, headerRow); }} aria-label="Delimiter" className={SEL}>{DELIMS.map(([v, t]) => <option key={v} value={v}>{t}</option>)}</select></label>
            <label className="flex items-center gap-1.5"><span className="font-semibold">Header row</span>
              <select value={headerRow} onChange={e => { const h = Number(e.target.value); setHeaderRow(h); guess(delim, h); }} aria-label="Header row" className={SEL}>
                <option value={-1}>none (every line is a call)</option>
                {lines.slice(0, 30).map((l, i) => l.trim() ? <option key={i} value={i}>line {i + 1}: {l.slice(0, 28)}</option> : null)}</select></label>
          </div>
          <div className="overflow-x-auto border border-slate-200 rounded" data-map-preview>
            <table className="text-[11px] whitespace-nowrap">
              <thead><tr className="bg-slate-50 text-left">{header.map((h, i) => <th key={i} className="px-2 py-1 font-semibold border-b border-slate-200">{h.trim() || `column ${i + 1}`}</th>)}</tr></thead>
              <tbody>{data.slice(0, 8).map((r, j) => <tr key={j} className="border-b border-slate-100">{header.map((_, i) => <td key={i} className="px-2 py-0.5 text-slate-600 max-w-[180px] truncate">{r[i] ?? ''}</td>)}</tr>)}</tbody>
            </table>
          </div>
          {mode === 'probes' ? <div data-map-probes>
            {pSel('Chromosome', pc.chrom, v => setP('chrom', v), true)}
            {pSel('Position or start', pc.start, v => setP('start', v), true)}
            {pSel('End', pc.end, v => setP('end', v), false)}
            <div className="text-[11px] text-slate-500 mt-1.5">Log ratio: one column per sample, each a patient (attached to a loaded patient of that name, else a probes-only patient):</div>
            {(pc.lr.length ? pc.lr : [-1]).map((j, k) => <div key={k} className="flex items-center gap-2 py-0.5">
              {pSel(k ? `Log ratio ${k + 1}` : 'Log ratio', j, v => setPc(p => ({ ...p, lr: p.lr.length ? p.lr.map((y, x) => (x === k ? v : y)) : [v] })), k === 0)}
              {j >= 0 && <input value={names[k] ?? ''} placeholder={sampleOfColumn(header[j] ?? '', lrCols.length > 1, item.name)} onChange={e => setNames(n => { const c = [...n]; c[k] = e.target.value; return c; })} aria-label={`Sample name of log ratio ${k + 1}`} className="h-[26px] px-1.5 rounded border border-slate-200 text-[11.5px] w-[160px]" />}
              {k > 0 && <button type="button" className="text-slate-500 hover:text-red-700" onClick={() => { setPc(p => ({ ...p, lr: p.lr.filter((_, x) => x !== k) })); setNames(n => n.filter((_, x) => x !== k)); }} aria-label={`Remove log ratio ${k + 1}`}>×</button>}
            </div>)}
            {pc.lr[0] >= 0 && <button type="button" className="text-accent-700 hover:underline" onClick={() => setPc(p => ({ ...p, lr: [...p.lr, -1] }))}>+ another sample</button>}
            {pSel('Probe name', pc.name, v => setP('name', v), false)}
          </div> : <div className="grid sm:grid-cols-2 gap-x-6">
            <div>{PLACE.map(f => field(f, true))}
              <div className="text-[11px] text-slate-500 mt-1">or the place in one column (chr1:72523846-72667555, or an ISCN):</div>{field('pos', false)}
              <div className="text-[11px] text-slate-500 mt-1.5">At least one of these (in this order of preference). With no Type, the type comes from the copy number (below 2 loss, above gain; 1 is normal on a male's X and a Y), else from the log ratio's sign; a row with neither is skipped, counted:</div>
              {ONE_OF.map(f => <div key={f}>{field(f, false)}
                {f === 'lr' && lrMore.map((j, k) => <label key={k} className="flex items-center gap-2 py-0.5"><span className="w-[170px] shrink-0 pl-3 text-slate-500">and log ratio</span>
                  <select value={j} onChange={e => { const v = Number(e.target.value); setLrMore(m => (v < 0 ? m.filter((_, x) => x !== k) : m.map((y, x) => (x === k ? v : y)))); }} aria-label={`Log ratio column ${k + 2}`} className={SEL}>
                    <option value={-1}>— none</option>{header.map((h, i) => <option key={i} value={i}>{h.trim() || `column ${i + 1}`}</option>)}</select></label>)}</div>)}</div>
            <div><div className="text-[11px] text-slate-500">Optional</div>{OPTIONAL.map(f => field(f, false))}</div>
          </div>}
          <div className="flex flex-wrap items-center gap-x-5 gap-y-1">
            <label className="flex items-center gap-1.5 cursor-pointer"><input type="checkbox" checked={zeroBased} onChange={e => setZeroBased(e.target.checked)} />{mode === 'probes' ? 'Positions are 0-based (BED-like; else 1-based)' : 'Start is 0-based (BED-like)'}</label>
            {!statesBuild && <label className="flex items-center gap-1.5"><span className="font-semibold">Genome build</span>
              <select value={build} onChange={e => setBuild(e.target.value as 'GRCh37' | 'GRCh38')} aria-label="Genome build" className={SEL}><option>GRCh37</option><option>GRCh38</option></select><span className="text-slate-500">the file states none</span></label>}
          </div>
          {mode === 'calls' && odd.length > 0 && <div className="rounded-lg border border-amber-300 bg-amber-50 px-3 py-2">
            <div className="font-semibold text-amber-900 mb-1">Type values that say neither loss nor gain</div>
            <div className="flex flex-wrap gap-x-4 gap-y-1">{odd.map(v => <label key={v} className="flex items-center gap-1.5"><span className="font-mono">{v}</span>
              <select value={typeMap[v.toLowerCase()] ?? 'ignore'} onChange={e => setTypeMap(m => ({ ...m, [v.toLowerCase()]: e.target.value as TypeChoice }))} aria-label={`Meaning of ${v}`} className={SEL}><option>ignore</option><option>loss</option><option>gain</option></select></label>)}</div>
          </div>}
          {headerRow >= 0 && !headerOk && <div className="text-[11.5px] text-amber-800" data-map-noremember>This row does not look like a header: the choice is not remembered.</div>}
          <div className={`text-[11.5px] ${(mode === 'probes' ? probeOk : res.file) ? 'text-slate-700' : 'text-red-700'}`} data-map-status>{mode === 'probes' ? probeMsg : res.msg}</div>
        </div>
        <div className="px-4 py-2.5 border-t border-slate-200 flex justify-end gap-2">
          <button onClick={cancel} className="h-[30px] px-3 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50">Cancel</button>
          <button disabled={mode === 'probes' ? !probeOk : !res.file || res.file.samples.length === 0} onClick={use} className="h-[30px] px-4 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700 disabled:opacity-40 disabled:cursor-not-allowed">Use these columns</button>
        </div>
      </div>
    </div>
  );
}
