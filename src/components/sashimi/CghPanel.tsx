/**
 * The array-CGH page's frame: the button that adds samples, the questions asked once files are read (columns, minimum
 * number of probes), and the table of each patient's calls with their annotations. A call clicked is drawn by the page
 * in the SV browser's own browser, as a deletion or a duplication.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import { labHits, parseLabBenign, LAB_COVER, type LabHit } from './cghLabBenign';
import { topoLabel, topoWords, rearrangementHints } from './cghTopo';
import { seriesOf, duplicateNames, SERIES_OVERLAP, type SeriesInfo } from './cghSeries';
import { mosaicOf, MOSAIC_TAG_BELOW } from './cghMosaic';
import { CghMapDialog } from './CghMapDialog';
import { decodeText } from './cghFormats';
import { Icon, Stepper, Segmented, Pill, Popover } from './controls';
import { cgh, useCgh, passes, classOf, noteOf, isAnalysed, largestFirst, CGH_ACCEPT, CGH_FILES_TEXT, CGH_MIN_PROBES_DEFAULT, CGH_MAX_ANNOTATED, SERIES_OPS, type SeriesOp, type CghState } from './cghStore';
import { CGH_CLASSES, classColor, cnvHubUrl, copyText, exportCghFile, type CghCnv, type CghFile, type CghSex } from './cghIntervals';
import { annotate, retry, viewOf, CGH_SOURCES, type CghSource, bandOfCall, PANEL_SOURCE, PANELAPP_MAX_GENES, DGV_COVER, type DgvRelation } from './cghAnnotate';
import { LEVEL_RANK, LEVEL_STYLE, SCORE_RULES, type CnvGene, type GeneView, type Level, type RegionHit } from './cghScore';
import { MOI, type SvGeneLoaders } from './SvGenomeBrowser';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { MoiClass } from '../../standalone/panelapp';

const fmt = (n: number) => Math.round(n).toLocaleString('en-US');
const sizeText = (kb: number) => (kb >= 1000 ? `${(kb / 1000).toFixed(2)} Mb` : kb >= 1 ? `${kb.toFixed(kb >= 100 ? 0 : 1)} kb` : `${Math.round(kb * 1000)} bp`);
const KIND = { loss: { color: '#dc2626', bg: '#fee2e2', word: 'del' }, gain: { color: '#2563eb', bg: '#dbeafe', word: 'dup' } } as const;

/** "＋ Add samples": the files picked (several at once), read in the page. */
function CghAddButton({ compact }: { compact?: boolean }) {
  return (
    <label className={`inline-flex items-center gap-1.5 h-[30px] px-2.5 rounded-[6px] border text-[12.5px] font-medium cursor-pointer whitespace-nowrap bg-white border-slate-200 text-slate-700 hover:border-accent-300 hover:bg-accent-50 ${compact ? 'h-[26px] px-2 text-xs rounded-lg' : ''}`}
      title={`Array-CGH calls of one or more patients: ${CGH_FILES_TEXT}.\nRead in this page, never sent anywhere. Files can also be dropped on the page.`}>
      <Icon name="upload" size={compact ? 13 : 15} />{compact ? 'Add samples' : 'Add array-CGH samples'}
      <input data-cgh-input type="file" multiple accept={CGH_ACCEPT} className="hidden"
        onChange={e => { const fs = e.target.files; if (fs?.length) void cgh.add(fs); e.target.value = ''; }} />
    </label>
  );
}

/** The question asked once the files are read: the minimum number of probes, with what it keeps of each sample. */
export function CghProbeDialog() {
  const s = useCgh();
  if (s.toMap.length) return <CghMapDialog key={s.toMap[0].name} item={s.toMap[0]} />;
  return s.pending ? <ProbeQuestion p={s.pending} init={s.filter.minProbes} initKb={s.filter.minKb} loaded={s.files} /> : null;
}
/** the patient's sex: changes the priority of X-linked recessive genes and the copies drawn for chrX / chrY */
const SEX_MARK: Record<CghSex, string> = { male: '♂', female: '♀' };
function SexSelect({ value, onChange, title }: { value: CghSex | null; onChange: (v: CghSex | null) => void; title?: string }) {
  return <select value={value ?? ''} onChange={e => onChange((e.target.value || null) as CghSex | null)} aria-label="Sex of the patient"
    title={title ?? "The patient's sex: a male's chrX loss is hemizygous (a green X-linked recessive gene is then HIGH), and his chrX / chrY are drawn with one copy"}
    className={`h-[24px] px-1 rounded border border-slate-200 bg-white text-[11.5px] outline-none cursor-pointer ${value ? 'text-slate-800 font-medium' : 'text-slate-400'}`}>
    <option value="">sex not given</option><option value="male">♂ male</option><option value="female">♀ female</option>
  </select>;
}
/** a size in kb as the user reads it: 500 bp, 250 kb, 1.5 Mb (3 significant digits at most) */
const sig3 = (x: number) => Number(x.toPrecision(3));
const sizeWord = (v: number) => (v <= 0 ? 'any' : v >= 1000 ? `${sig3(v / 1000)} Mb` : v < 1 ? `${sig3(v * 1000)} bp` : `${sig3(v)} kb`);
/** the slider's span, kb, on a log scale; position 0 is "any" (no minimum) */
const SLIDE_LO = 10, SLIDE_HI = 5000, SLIDE_STEPS = 300;
const SIZE_UNITS = { bp: 0.001, kb: 1, Mb: 1000 } as const;
type SizeUnit = keyof typeof SIZE_UNITS;
const kbToPos = (kb: number) => (kb <= 0 ? 0 : 1 + Math.round(Math.min(1, Math.max(0, Math.log(kb / SLIDE_LO) / Math.log(SLIDE_HI / SLIDE_LO))) * (SLIDE_STEPS - 1)));
// rounded to 2 significant digits so that the slider lands on readable values (250 kb, 1.2 Mb)
const posToKb = (p: number) => (p <= 0 ? 0 : Number((SLIDE_LO * (SLIDE_HI / SLIDE_LO) ** ((p - 1) / (SLIDE_STEPS - 1))).toPrecision(2)));
const unitOf = (kb: number): SizeUnit => (kb >= 1000 ? 'Mb' : kb > 0 && kb < 1 ? 'bp' : 'kb');
/**
 * The minimum size: a log slider (any, 10 kb to 5 Mb) and the value typed with its unit, kept in sync; a typed value may
 * lie outside the slider's span (500 bp, 8 Mb). The value is committed on release, or `delay` ms after the last move,
 * so that the calls are not re-annotated at every pixel.
 */
function SizeSlider({ value, onChange, delay = 300, className = '' }: { value: number; onChange: (kb: number) => void; delay?: number; className?: string }) {
  const [kb, setKb] = useState(value);
  const [unit, setUnit] = useState<SizeUnit>(unitOf(value));
  const [text, setText] = useState(value > 0 ? String(sig3(value / SIZE_UNITS[unitOf(value)])) : '');
  const timer = useRef<number | undefined>(undefined);
  const sent = useRef(value);
  const show = (v: number) => { const u = unitOf(v); setUnit(u); setText(v > 0 ? String(sig3(v / SIZE_UNITS[u])) : ''); };
  // a value changed elsewhere (the other control, a stored filter)
  useEffect(() => { if (value !== sent.current) { sent.current = value; setKb(value); show(value); } }, [value]);
  useEffect(() => () => window.clearTimeout(timer.current), []);
  const commit = (v: number) => { window.clearTimeout(timer.current); if (v !== sent.current) { sent.current = v; onChange(v); } };
  const later = (v: number) => { setKb(v); window.clearTimeout(timer.current); if (delay <= 0) commit(v); else timer.current = window.setTimeout(() => commit(v), delay); };
  const typed = (t: string, u: SizeUnit) => { setText(t); setUnit(u); const n = Number(t.replace(',', '.')); if (t.trim() === '') later(0); else if (Number.isFinite(n) && n >= 0) later(n * SIZE_UNITS[u]); };
  return (
    <span className={`inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-xs text-slate-500 ${className}`}
      title={'Smaller calls are not analysed: neither listed nor annotated.\nSlider: any (left end), then 10 kb to 5 Mb on a log scale; a value outside can be typed (500 bp, 8 Mb).'}>
      size
      <input type="range" min={0} max={SLIDE_STEPS} step={1} value={kbToPos(kb)} aria-label="Minimum size (slider)" aria-valuetext={sizeWord(kb)}
        onChange={e => { const v = posToKb(Number(e.target.value)); show(v); later(v); }}
        onPointerUp={() => commit(kb)} onKeyUp={() => commit(kb)} onBlur={() => commit(kb)}
        className="w-[110px] h-[14px] accent-accent-600 cursor-pointer" />
      <span className="w-[62px] text-right font-medium text-slate-800 tabular-nums whitespace-nowrap">{kb > 0 ? `≥ ${sizeWord(kb)}` : 'any'}</span>
      <input type="text" inputMode="decimal" value={text} placeholder="–" aria-label="Minimum size" onKeyDown={e => { e.stopPropagation(); if (e.key === 'Enter') commit(kb); }}
        onChange={e => typed(e.target.value, unit)} onBlur={() => commit(kb)}
        className="w-[42px] h-[20px] px-1 rounded border border-slate-200 bg-white font-medium text-slate-800 outline-none focus:border-accent-300 tabular-nums" />
      <select value={unit} onChange={e => typed(text, e.target.value as SizeUnit)} aria-label="Unit of the minimum size" className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
        {(Object.keys(SIZE_UNITS) as SizeUnit[]).map(u => <option key={u} value={u}>{u}</option>)}
      </select>
    </span>
  );
}
const EXPORT_GAP_MS = 600;
const EXPORT_TITLE = 'Each file of the tab written back as it came, with the classes given here in its Classification column and the notes in a Comment column (added when the file has none; a BED gets both as two last columns, a VCF in INFO CGH_CLASS / CGH_NOTE; an .xlsx is written as tab-separated text, .classified.txt).\nClasses and notes are kept in this tab only: export them to keep them.\nSeveral files: one download each, and the browser may ask to allow several downloads.';
function ProbeQuestion({ p, init, initKb, loaded }: { p: NonNullable<CghState['pending']>; init: number; initKb: number; loaded: CghFile[] }) {
  const [n, setN] = useState<number>(init);
  const [kb, setKb] = useState<number>(initKb);
  // every call of the tab (loaded and added): past the cap, the largest are annotated now and the others on demand
  const every = [...loaded, ...p.files].flatMap(f => f.samples.flatMap(x => x.cnvs));
  const kept = every.filter(c => isAnalysed(c, n, kb)).length;
  const many = every.length > CGH_MAX_ANNOTATED, over = kept > CGH_MAX_ANNOTATED;
  const [sexes, setSexes] = useState<Record<string, CghSex | null>>(() => Object.fromEntries(p.files.flatMap(f => f.samples.map(x => [x.key, x.sex]))));
  const samples = p.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const noCount = p.files.some(f => f.format !== 'interval' || f.warnings.some(w => /no #Probes/.test(w)));
  return (
    <div className="fixed inset-0 z-[80] bg-black/40 flex items-start justify-center p-6" role="dialog" aria-modal="true" aria-label="Minimum number of probes"
      onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); cgh.cancel(); } if (e.key === 'Enter' && samples.length && !/^(SELECT|BUTTON|A)$/.test((e.target as HTMLElement).tagName)) cgh.accept(n, sexes, kb); }}>
      <div className="bg-white rounded-xl shadow-2xl border border-slate-200 w-full max-w-lg text-slate-900 text-xs mt-16 outline-none" tabIndex={-1} ref={el => { if (el && !el.contains(document.activeElement)) el.focus(); }}>
        <div className="px-4 pt-3 pb-2 border-b border-slate-200">
          <div className="font-bold text-sm">Array-CGH samples</div>
          <div className="text-[11px] text-slate-500 mt-0.5">Calls with fewer probes, or smaller, than this are not analysed: neither listed nor annotated (both can be changed in the list afterwards).</div>
        </div>
        <div className="px-4 py-3 space-y-3">
          {samples.length > 0 && <div className="flex items-center gap-2">
            <span className="font-semibold">Minimum number of probes</span>
            <Stepper value={n} onChange={setN} min={1} max={500} title="Calls with fewer probes are not analysed" ariaLabel="Minimum number of probes" />
            {n !== CGH_MIN_PROBES_DEFAULT && <button className="text-accent-700 hover:underline" onClick={() => setN(CGH_MIN_PROBES_DEFAULT)}>back to {CGH_MIN_PROBES_DEFAULT}</button>}
          </div>}
          {samples.length > 0 && <div className="flex items-center gap-2">
            <span className="font-semibold">Minimum size</span>
            <SizeSlider value={kb} onChange={setKb} delay={120} />
          </div>}
          {many && <div className="rounded-lg border px-3 py-2 text-[11.5px] border-amber-300 bg-amber-50 text-amber-900" role="status" data-cgh-count>
            <b>{every.length.toLocaleString('en-US')} calls in all</b> (this tab's patients and these). The calls with {n}+ probes{kb ? ` and of ${sizeWord(kb)} or more` : ''} are listed and annotated (genes, PanelApp, ClinGen, gnomAD, DGV Gold); the others are listed per patient on demand.
            {over ? <> <b>{kept.toLocaleString('en-US')} calls pass the minimums: the {CGH_MAX_ANNOTATED} largest are annotated now, the others on demand.</b></> : <> {kept.toLocaleString('en-US')} will be.</>}
          </div>}
          {samples.length > 0 && <table className="w-full text-[11.5px]">
            <thead><tr className="text-left text-slate-500 border-b border-slate-200"><th className="py-1 font-medium">Sample</th><th className="font-medium">File</th><th className="font-medium">Build</th><th className="font-medium">Sex</th><th className="font-medium text-right">Calls</th><th className="font-medium text-right">Analysed</th></tr></thead>
            <tbody>{samples.map(({ x, f }) => {
              const kept = x.cnvs.filter(c => isAnalysed(c, n, kb)).length;
              return <tr key={x.key} className="border-b border-slate-100">
                <td className="py-1 font-semibold">{x.name}</td><td className="text-slate-500 truncate max-w-[150px]" title={f.name}>{f.name}</td>
                <td className={f.build ? '' : 'text-amber-700'} title={f.buildText ? `the file says Genome ${f.buildText}` : 'the file states no build: read as GRCh37'}>{f.build ?? 'GRCh37 ?'}</td>
                <td><SexSelect value={sexes[x.key] ?? null} onChange={v => setSexes(m => ({ ...m, [x.key]: v }))} title={x.sex ? `the file says ${x.sex}` : undefined} /></td>
                <td className="text-right tabular-nums">{x.cnvs.length}</td><td className="text-right tabular-nums font-semibold">{kept}</td>
              </tr>;
            })}</tbody>
          </table>}
          {noCount && <div className="text-[11px] text-amber-700">A file gives no probe count (BED, VCF): its calls are all kept.</div>}
          {p.files.flatMap(f => f.warnings.filter(w => !/^(BED|VCF): no probe/.test(w)).map(w => <div key={f.id + w} className="text-[11px] text-amber-700">{f.name}: {w}</div>))}
          {p.files.filter(f => f.spec).map(f => <div key={f.id + 'cols'} className="text-[11px] text-slate-500">{f.name.split(/[›/\\]/).pop()}: {f.remembered ? 'columns as mapped before ' : ''}
            <button className="text-accent-700 hover:underline" onClick={() => cgh.remap(f)}>{f.remembered ? '(change…)' : 'columns…'}</button></div>)}
          {p.notes.map(e => <div key={e} className="text-[11px] text-slate-500">{e}</div>)}
          {duplicateNames([...loaded, ...p.files]).map(e => <div key={e} className="text-[11px] text-amber-700" data-cgh-dup>{e}</div>)}
          {p.errors.map(e => <div key={e} className="text-[11px] text-red-700">{e}</div>)}
        </div>
        <div className="px-4 py-2.5 border-t border-slate-200 flex justify-end gap-2">
          <button onClick={() => cgh.cancel()} className="h-[30px] px-3 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50">Cancel</button>
          {samples.length > 0 && <button autoFocus onClick={() => cgh.accept(n, sexes, kb)} className="h-[30px] px-4 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700 disabled:opacity-40 disabled:cursor-not-allowed">Analyse {samples.length} sample{samples.length > 1 ? 's' : ''}</button>}
        </div>
      </div>
    </div>
  );
}

/**
 * The frame of the array-CGH calls: one sample at a time (or all), filtered (probes, type, size, class, a gene or a
 * band), each call a row: clicked, it is drawn in the browser below; its class is chosen in the row; its CNV-Hub page
 * opens from it. ↑ ↓ go from one call to the next.
 */
export function CghFrame({ onShow, dropping, loadersFor }: { onShow: (c: CghCnv, f: CghFile) => void; dropping?: boolean; loadersFor?: (b: GenomeBuild) => Omit<SvGeneLoaders, 'fallback'> }) {
  const s = useCgh();
  /** a note opened from its row's ✎: its box takes the focus */
  const [noteFocus, setNoteFocus] = useState<string | null>(null);
  /** the detail under the call clicked folded (a second click on its row) */
  const [folded, setFolded] = useState(false);
  // the analysed calls of every patient, the CGH_MAX_ANNOTATED largest: their genes (RefSeq), dosage, constraint, PanelApp, AChro-Puce
  const buildOf = useMemo(() => new Map(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.map(c => [c.id, fl.build ?? 'GRCh37'] as const)))), [s.files]);
  const auto = useMemo(() => new Set(largestFirst(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.filter(c => isAnalysed(c, s.filter.minProbes, s.filter.minKb)))))
    .slice(0, CGH_MAX_ANNOTATED).map(c => c.id)), [s.files, s.filter.minProbes, s.filter.minKb]);
  useEffect(() => {
    if (!loadersFor) return;
    annotate(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.filter(c => auto.has(c.id)).map(c => ({ c, build: fl.build ?? 'GRCh37' })))), loadersFor);
  }, [auto, loadersFor]); // eslint-disable-line react-hooks/exhaustive-deps
  /** on demand: up to CGH_MAX_ANNOTATED of these, the largest first (the queue skips those annotated or queued) */
  const annotateNow = (cs: CghCnv[]) => { if (loadersFor) annotate(largestFirst(cs.filter(c => !s.ann[c.id])).slice(0, CGH_MAX_ANNOTATED).map(c => ({ c, build: buildOf.get(c.id) ?? 'GRCh37' })), loadersFor); };
  /** ask again, for these calls, only the sources given that failed (the queue skips the others) */
  const retryNow = (cs: CghCnv[], srcs: (CghSource | 'panelapp')[]) => { if (loadersFor) retry(cs.map(c => ({ c, build: buildOf.get(c.id) ?? 'GRCh37' })), srcs, loadersFor); };
  // the call clicked kept in the table's view (the table only, never the page), when it changes
  const tableBox = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const box = tableBox.current, row = box?.querySelector<HTMLElement>('tr[data-shown]');
    if (!box || !row) return;
    const r = row.getBoundingClientRect(), bb = box.getBoundingClientRect(), head = 26;
    if (r.top < bb.top + head) box.scrollTop -= bb.top + head - r.top;
    else if (r.bottom > bb.bottom) box.scrollTop += r.bottom - bb.bottom;
  }, [s.shown]);
  // each call's genes and priority, recomputed when an answer comes, not at each keystroke (viewOf reads ann, modes, sex)
  const views = useMemo(() => new Map(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs)).map(c => [c.id, viewOf(s, c)])), [s.files, s.ann, s.modes, s.sex]); // eslint-disable-line react-hooks/exhaustive-deps
  const hints = useMemo(() => rearrangementHints(s), [s.files, s.ann]); // eslint-disable-line react-hooks/exhaustive-deps
  if (!s.files.length) return null;
  const samples = s.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const listed = s.sample ? samples.filter(e => e.x.key === s.sample) : samples;
  const all = listed.flatMap(({ x, f }) => x.cnvs.map(c => ({ c, f })));
  const series = seriesOf(s.files);
  const lab = labHits(s.files, s.lab);
  const rank = (id: string) => { const v = views.get(id)?.score; return v ? LEVEL_RANK[v.level] : -1; };
  const byRank = (l: typeof all) => { if (s.filter.sort === 'priority') l.sort((a, b) => rank(b.c.id) - rank(a.c.id)); return l; };
  const isBelow = (c: CghCnv) => !isAnalysed(c, s.filter.minProbes, s.filter.minKb);
  // each patient's calls under the minimums (the other filters applied), listed after its own when revealed
  const below = new Map(listed.map(({ x, f }) => [x.key, byRank(x.cnvs.filter(c => isBelow(c) && passes(s, c, true)).map(c => ({ c, f })))]));
  /** a patient's line under its calls: N more under the minimums (to reveal), or hide them */
  type Item = { row: { c: CghCnv; f: CghFile } } | { more: { key: string; name: string; n: number; open: boolean } };
  const blockOf = (x: (typeof listed)[number]['x']): Item[] => {
    const b = below.get(x.key) ?? [], open = !!s.revealed[x.key];
    if (!b.length) return [];
    return open ? [{ more: { key: x.key, name: x.name, n: b.length, open } }, ...b.map(row => ({ row })), { more: { key: x.key, name: x.name, n: -1, open } }] : [{ more: { key: x.key, name: x.name, n: b.length, open } }];
  };
  // in genomic order the rows stay grouped by patient (each patient's line ends its block); by priority they mix, the lines at the bottom
  const grouped = s.filter.sort === 'position' || listed.length === 1;
  const items: Item[] = grouped
    ? listed.flatMap(({ x, f }) => [...byRank(x.cnvs.filter(c => passes(s, c)).map(c => ({ c, f }))).map(row => ({ row })), ...blockOf(x)])
    : [...byRank(all.filter(({ c }) => passes(s, c))).map(row => ({ row })), ...listed.flatMap(({ x }) => blockOf(x))];
  const rows = items.flatMap(it => ('row' in it ? [it.row] : []));
  const fewProbes = all.filter(({ c }) => isBelow(c)).length;
  // listed and not annotated: past the automatic set (or revealed under the minimums), never asked yet
  const notAnn = rows.filter(({ c }) => !s.ann[c.id] && !auto.has(c.id)).map(({ c }) => c);
  const minWords = `fewer than ${s.filter.minProbes} probe${s.filter.minProbes > 1 ? 's' : ''}${s.filter.minKb ? ` or smaller than ${sizeWord(s.filter.minKb)}` : ''}`;
  const anySuppressed = all.some(({ c }) => c.suppressed);
  // per column, each listed patient's calls (PanelApp: genes) to ask again; the control is hidden when none failed
  const retryOf = (col: RetryCol): RetryGroup[] => listed.map(({ x }) => {
    const cs = rows.filter(r => r.c.sampleKey === x.key).map(r => r.c), src = RETRY_SRC[col];
    const n = col === 'panelapp'
      ? new Set(cs.flatMap(c => (s.ann[c.id]?.genes ?? []).slice(0, PANELAPP_MAX_GENES).map(g => g.symbol).filter(g => s.modes[g] === 'error'))).size
      : cs.filter(c => src.some(k => s.ann[c.id]?.failed[k as CghSource] && !s.ann[c.id]!.asking.includes(k as CghSource))).length;
    const asking = col === 'panelapp' ? 0 : cs.filter(c => src.some(k => s.ann[c.id]?.asking.includes(k as CghSource))).length;
    return { key: x.key, name: x.name, n, asking, run: () => retryNow(cs, src) };
  });
  const rc = (col: RetryCol) => <RetryCtl col={col} groups={retryOf(col)} />;
  const f = s.filter;
  const moreN = (s.lab && f.lab !== 'any' ? 1 : 0) + (samples.length > 1 && f.series !== 'any' ? 1 : 0) + (anySuppressed && f.hideSuppressed ? 1 : 0);
  /** a call clicked: drawn below with its detail under its row; clicked again, its detail folds or unfolds */
  const show = (i: number, click = false) => {
    const r = rows[i]; if (!r) return;
    if (click && r.c.id === s.shown) { setFolded(x => !x); return; }
    cgh.setShown(r.c.id); onShow(r.c, r.f);
  };
  const at = rows.findIndex(r => r.c.id === s.shown);
  const exportFiles = [...new Set(listed.map(e => e.f))];
  /** one file written back; a file holds every sample it was read with: all of them are marked exported, only once its download is triggered */
  const exportOne = (fl: CghFile) => {
    const text = exportCghFile(fl, c => classOf(s, c), c => noteOf(s, c));
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain;charset=utf-8' }));
    a.download = fl.saveAs;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    cgh.exported(fl.samples.map(x => x.key));
  };
  // browsers block or ask before several downloads at once: one after another, each file marked as its own starts
  const exportAll = () => exportFiles.forEach((fl, i) => setTimeout(() => exportOne(fl), i * EXPORT_GAP_MS));
  /** classes or notes of the patients listed not exported yet */
  const listedDirty = listed.some(e => s.dirtyKeys[e.x.key]);
  const tab = (on: boolean) => `inline-flex items-center gap-1.5 h-8 pl-3 pr-2 -mb-px rounded-t-lg border text-xs whitespace-nowrap ${on ? 'bg-white border-slate-200 border-b-white text-slate-900 font-semibold shadow-[inset_0_2px_0_#0f766e]' : 'bg-slate-50 border-transparent text-slate-600 hover:text-slate-900 hover:bg-slate-100'}`;
  const filesBuilds = [...new Set(listed.map(e => e.f.build ?? 'GRCh37'))];
  return (
    <div className={`mx-4 mt-3 bg-white border rounded-xl shadow-sm ${dropping ? 'border-accent-400 ring-2 ring-accent-100' : 'border-slate-200'}`} data-cgh-frame>
      {/* one tab per patient, and all of them in one table */}
      <div className="px-5 pt-2 flex flex-wrap items-end gap-1 border-b border-slate-200">
        <div className="text-base font-bold mr-2 self-center">Array-CGH</div>
        <div className="flex flex-wrap items-end gap-1" role="tablist" aria-label="Patients">
        {samples.length > 1 && <button role="tab" aria-selected={!s.sample} className={tab(!s.sample)} onClick={() => cgh.setSample(null)}
          title="Every patient in one table, the sample as a column">All patients <span className="text-slate-400 font-normal tabular-nums">{samples.reduce((n, e) => n + e.x.cnvs.filter(c => passes(s, c)).length, 0)}</span></button>}
        {samples.map(({ x, f: fl }) => {
          const n = x.cnvs.filter(c => passes(s, c)).length, on = s.sample === x.key || samples.length === 1;
          return <span key={x.key} role="tab" aria-selected={on} className={tab(on)}>
            <button onClick={() => cgh.setSample(x.key)} title={`${fl.name}${fl.meta.length ? `\n${fl.meta.map(([k, v]) => `${k}: ${v}`).join('\n')}` : ''}`}>{x.name}{s.sex[x.key] ? <span className="ml-0.5 text-slate-500" title={s.sex[x.key]}>{SEX_MARK[s.sex[x.key]]}</span> : null} <span className="text-slate-400 font-normal tabular-nums">{n}/{x.cnvs.length}</span></button>
            <button className="text-slate-400 hover:text-red-600 leading-none text-sm px-0.5" title="Take this patient off the list"
              onClick={() => { if (!s.dirtyKeys[x.key] || confirm(`${x.name}: the classes and notes given here are not exported yet and will be lost.`)) cgh.remove(x.key); }}>×</button>
          </span>;
        })}
        </div>
        <span className="ml-auto self-center flex flex-wrap items-center justify-end gap-x-2 gap-y-1 py-1">
          <CghAddButton compact />
          <LabBenignButton />
          <span className="text-[11px] text-slate-500" title="The build the interval table states (Genome hg19); the calls are drawn on it">{filesBuilds.join(', ')}{listed.some(e => !e.f.build) ? ' (no build in a file: GRCh37 assumed)' : ''}</span>
{exportFiles.length > 1
            ? <Popover width={340} title={EXPORT_TITLE} label="Export classified tables"
                buttonClass={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium ${listedDirty ? 'bg-amber-50 border-amber-300 text-amber-800' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300'}`}
                button={<><Icon name="download" size={13} />Export classified tables{listedDirty ? ' •' : ''} ▾</>}>
                <div className="text-xs" data-cgh-export>
                  {exportFiles.map(fl => <div key={fl.id} className="flex items-center gap-2 py-1 border-b border-slate-100">
                    <span className="flex-1 min-w-0 truncate" title={fl.name}>{fl.name}{fl.samples.some(x => s.dirtyKeys[x.key]) ? <span className="text-amber-700"> •</span> : null}
                      <span className="block text-[10.5px] text-slate-500 truncate">{fl.samples.map(x => x.name).join(', ')} → {fl.saveAs}</span></span>
                    <button className="h-[24px] px-2 rounded border border-slate-200 hover:border-slate-300 font-medium" onClick={() => exportOne(fl)}>export</button>
                  </div>)}
                  <div className="flex items-center gap-2 pt-1.5">
                    <span className="flex-1 text-[10.5px] text-slate-500">The browser may ask to allow several downloads.</span>
                    <button className="h-[24px] px-2 rounded bg-accent-600 text-white font-semibold hover:bg-accent-700" title="Each file in turn, a moment apart; the browser may ask to allow several downloads" onClick={exportAll}>export all ({exportFiles.length})</button>
                  </div>
                </div>
              </Popover>
            : <button onClick={() => exportFiles[0] && exportOne(exportFiles[0])} className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium ${listedDirty ? 'bg-amber-50 border-amber-300 text-amber-800' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300'}`}
                title={EXPORT_TITLE}>
                <Icon name="download" size={13} />Export classified table{listedDirty ? ' •' : ''}</button>}
        </span>
      </div>
      <div className="pt-2" />
      <div className="px-5 pb-2 flex flex-wrap items-center gap-2 text-xs">
        {listed.length === 1 && <SexSelect value={s.sex[listed[0].x.key] ?? null} onChange={v => cgh.setSex(listed[0].x.key, v)} />}
        <span className="[&>span]:!h-[26px] [&>span]:!rounded-lg [&>span]:!text-xs"><Stepper value={f.minProbes} onChange={v => cgh.setFilter({ minProbes: v })} min={1} max={500} label="probes ≥" title="Calls with fewer probes are not analysed (3 by default)" ariaLabel="Minimum number of probes" /></span>
        <Segmented size="sm" label="Type" title="Losses, gains or both" value={f.kind} onChange={v => cgh.setFilter({ kind: v })}
          options={[{ value: 'all', label: 'all', hint: 'Losses and gains' }, { value: 'loss', label: 'losses', hint: 'Deletions only' }, { value: 'gain', label: 'gains', hint: 'Gains and amplifications only' }]} />
        <SizeSlider value={f.minKb} onChange={v => cgh.setFilter({ minKb: v })} />
        <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title="Classified or not (a class from the file counts)">class
          <select value={f.classed} onChange={e => cgh.setFilter({ classed: e.target.value as typeof f.classed })} className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            <option value="all">any</option><option value="unclassified">not classified</option><option value="classified">classified</option>
          </select></label>
        <Popover width={300} title="Less used filters: lab list, series, suppressed calls" label="More filters"
          buttonClass="inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border border-slate-200 bg-white text-xs font-medium text-slate-700 hover:border-slate-300 whitespace-nowrap"
          button={<>More filters ▾{moreN > 0 && <span className="text-[10.5px] font-bold leading-4 px-1.5 rounded-full bg-accent-600 text-white">{moreN}</span>}</>}>
          <div className="flex flex-col items-start gap-2 p-1 text-xs">
        {s.lab && <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title={`The laboratory's benign list (${s.lab.fileName}): calls sharing ${Math.round(LAB_COVER * 100)} % of the call or more with an entry of the same type`}>lab list
          <select value={f.lab} onChange={e => cgh.setFilter({ lab: e.target.value as typeof f.lab })} className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            <option value="any">any</option><option value="not">not in the list</option><option value="in">in the list</option>
          </select></label>}
        {samples.length > 1 && <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title={`Recurrence in the loaded series: the number of OTHER patients carrying the call, as the Series column counts it (same type, same chromosome, reciprocal overlap ≥ ${Math.round(SERIES_OVERLAP * 100)} %)`}>series
          <select value={f.series} onChange={e => cgh.setFilter({ series: e.target.value as SeriesOp })} aria-label="Series: comparison" className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            {SERIES_OPS.map(o => <option key={o} value={o}>{o}</option>)}
          </select>
          {f.series !== 'any' && <><input type="number" min={0} max={Math.max(0, samples.length - 1)} step={1} value={f.seriesN} aria-label="Series: number of other patients" onKeyDown={e => e.stopPropagation()}
            onChange={e => { const v = Math.max(0, Math.round(Number(e.target.value) || 0)); cgh.setFilter({ seriesN: v }); }} className="w-[40px] bg-transparent font-medium text-slate-800 outline-none tabular-nums" />
            <span>other{f.seriesN === 1 ? '' : 's'}</span></>}</label>}
        {anySuppressed && <Pill size="sm" label="hide suppressed" on={f.hideSuppressed} onChange={v => cgh.setFilter({ hideSuppressed: v })} title="Calls marked Suppress = true in CytoGenomics" />}
          </div>
        </Popover>
        <input value={f.query} onChange={e => cgh.setFilter({ query: e.target.value })} placeholder="gene, band (16p11), chromosome" aria-label="Filter the calls by gene, band or chromosome"
          className="h-[26px] w-[210px] px-2 rounded-lg border border-slate-200 text-xs placeholder:text-slate-400 focus:outline-none focus:border-accent-300" />
        <Segmented size="sm" label="Order" title="The calls in genomic order, or the highest priority first" value={f.sort} onChange={v => cgh.setFilter({ sort: v })}
          options={[{ value: 'position', label: 'position', hint: 'Genomic order, as in the file' }, { value: 'priority', label: 'priority', hint: 'HIGH first, then MODERATE, LOW, none' }]} />
        <span className="text-slate-500 ml-1">{rows.length} of {all.length} call{all.length > 1 ? 's' : ''} listed{fewProbes ? ` · ${fewProbes} under ${f.minProbes} probe${f.minProbes > 1 ? 's' : ''}${f.minKb ? ` or ${sizeWord(f.minKb)}` : ''} not analysed` : ''}
          {notAnn.length > 0 && <span data-cgh-notann> · {notAnn.length} not annotated · <button className="text-accent-700 hover:underline font-medium" onClick={() => annotateNow(notAnn)}
            title={`Annotate the listed calls not annotated yet (up to ${CGH_MAX_ANNOTATED} a click, the largest first)`}>annotate them</button></span>}</span>
      </div>
      <div ref={tableBox} style={{ maxHeight: 'max(250px, 40vh)' }} className="mx-5 mb-1.5 border border-slate-200 rounded-lg overflow-auto outline-none focus:ring-2 focus:ring-accent-100" tabIndex={0} aria-label="Array-CGH calls (↑ ↓ to go from one to the next)"
        onKeyDown={e => { if ((e.key === 'ArrowDown' || e.key === 'ArrowUp') && !/^(SELECT|INPUT|TEXTAREA)$/.test((e.target as HTMLElement).tagName)) { e.preventDefault(); e.stopPropagation(); show(Math.max(0, Math.min(rows.length - 1, at < 0 ? 0 : at + (e.key === 'ArrowDown' ? 1 : -1)))); } }}>
        <table className="w-full text-[11.5px] border-collapse">
          <thead className="sticky top-0 z-10 bg-slate-50 text-slate-500 text-left whitespace-nowrap">
            <tr className="border-b border-slate-200">
              {!s.sample && samples.length > 1 && <th className="px-1.5 py-0.5 font-medium">Sample</th>}
              <th className="px-1.5 py-0.5 font-medium">Band{rc('bands')}</th><th className="px-1.5 py-0.5 font-medium">Position (GRCh37)</th>
              <th className="px-1.5 py-0.5 font-medium text-right">Size</th><th className="px-1.5 py-0.5 font-medium">Type{rc('dgv')}</th>
              <th className="px-1.5 py-0.5 font-medium text-right" title="#Probes · mean log ratio">Probes · LogR</th>
              {samples.length > 1 && <th className="px-1.5 py-0.5 font-medium text-right" title={`Recurrence in the loaded series: other patients with a call of the same type on the same chromosome (reciprocal overlap ≥ ${Math.round(SERIES_OVERLAP * 100)} %) / other patients loaded. Every call counts, those under the minimums too.`}>Series</th>}
              <th className="px-1.5 py-0.5 font-medium" title="RefSeq NM_ genes of the interval (UCSC, the file's build); the file's Gene Name only when they cannot be read">Genes{rc('genes')}</th>
              <th className="px-1.5 py-0.5 font-medium" title={`${PANEL_SOURCE.label}, every panel: the genes by their best rating, green first (with their mode of inheritance), then amber, then red`}>PanelApp{rc('panelapp')}</th>
              <th className="px-1.5 py-0.5 font-medium" title="ClinGen dosage scores (haploinsufficiency for a loss, triplosensitivity for a gain; 3 sufficient, 2 some evidence) and the lowest gnomAD LOEUF">Dosage{rc('dosage')}</th>
              <th className="px-1.5 py-0.5 font-medium" title="A priority for the review from the genes' scores (the rules are listed under the table); not a classification">Priority</th>
              <th className="px-1.5 py-0.5 font-medium sticky right-[72px] z-[1] bg-slate-50 shadow-[-6px_0_6px_-6px_rgba(15,23,42,0.18)]">Classification</th><th className="px-1.5 py-0.5 font-medium sticky right-0 z-[1] bg-slate-50 w-[72px] min-w-[72px]">Links</th>
            </tr>
          </thead>
          <tbody>
            {items.map(it => {
              if ('more' in it) {
                const m = it.more;
                return <tr key={`more-${m.key}-${m.n}`} data-cgh-more className="border-b border-slate-100 bg-slate-50/60"><td colSpan={13} className="px-2 py-0.5 text-[11px]">
                  {m.n < 0 ? <button className="text-accent-700 hover:underline" onClick={() => cgh.setRevealed(m.key, false)}>− hide them</button>
                    : m.open ? <span className="text-slate-500">{m.n} call{m.n > 1 ? 's' : ''} of {m.name} under the minimums ({minWords}), annotated on demand</span>
                    : <button className="text-accent-700 hover:underline" title={`List them, and annotate up to ${CGH_MAX_ANNOTATED} (the largest first)`}
                        onClick={() => { cgh.setRevealed(m.key, true); annotateNow((below.get(m.key) ?? []).map(r => r.c)); }}>+ {m.n} more call{m.n > 1 ? 's' : ''} of {m.name} under the minimums ({minWords})</button>}
                </td></tr>;
              }
              const { c, f: fl } = it.row, i = rows.indexOf(it.row);
              const low = isBelow(c), na = !s.ann[c.id] && !auto.has(c.id);
              const cls = classOf(s, c), col = classColor(cls), on = s.shown === c.id, k = KIND[c.kind];
              const inList = CGH_CLASSES.some(o => o.value === cls), note = noteOf(s, c);
              const hub = cnvHubUrl(c, fl.build ?? 'GRCh37');
              const v = views.get(c.id)!;
              return [
                <tr key={c.id} data-cgh-row data-cgh-below={low || undefined} data-shown={on || undefined} onClick={() => show(i, true)} title={`${c.iscn ?? ''}${c.suppressed ? '\nSuppressed in CytoGenomics' : ''}\nClick: draw it in the browser below`}
                  className={`border-b border-slate-100 cursor-pointer ${on ? 'bg-accent-50 [&>td]:border-t-2 [&>td]:border-t-accent-600 [&>td:first-child]:border-l-2 [&>td:first-child]:border-l-accent-600 [&>td:last-child]:border-r-2 [&>td:last-child]:border-r-accent-600' : 'bg-white hover:bg-slate-50'} ${c.suppressed ? 'opacity-60' : ''} ${low ? 'text-slate-500' : ''}`}>
                  {!s.sample && samples.length > 1 && <td className="px-1.5 py-0.5 text-slate-600 whitespace-nowrap">{c.sample}{s.sex[c.sampleKey] ? ` ${SEX_MARK[s.sex[c.sampleKey]]}` : ''}</td>}
                  <td className="px-1.5 py-0.5 font-semibold whitespace-nowrap" title={!c.cytoband && v.ann?.band ? 'Band read from the cytogenetic bands (the file gives none)' : undefined}>{c.chrom.replace(/^chr/, '')}{bandOfCall(c, v.ann).replace(/\s+/g, '')}<TopoTag c={c} v={v} hints={hints.get(c.id)} /></td>
                  <td className="px-1.5 py-0.5 font-mono text-[10.5px] text-slate-700 whitespace-nowrap">{c.chrom}:{fmt(c.start)}-{fmt(c.end)}</td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap">{sizeText(c.sizeKb)}</td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap"><span className="inline-flex items-center rounded px-1.5 text-[10.5px] font-bold whitespace-nowrap" style={{ color: k.color, background: k.bg }} title={c.typeText}>{k.word} {copyText(c)}</span><DgvMedal v={v} /><LabTag hit={lab.get(c.id)} /></td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap"><span className={c.probes == null ? 'text-slate-400' : ''}>{c.probes ?? '?'}</span>{c.logRatio != null && <span className="text-slate-500"> · {c.logRatio.toFixed(2)}</span>}<MosaicTag c={c} sex={s.sex[c.sampleKey]} />{low && <BelowTag c={c} minProbes={f.minProbes} minKb={f.minKb} />}</td>
                  {samples.length > 1 && <SeriesCell info={series.get(c.id)} />}
                  {na ? <><NotAnnCell /><NotAnnCell /><NotAnnCell /><td className="px-1.5 py-0.5 whitespace-nowrap" onClick={e => e.stopPropagation()}>
                    <button className="px-1 rounded border border-slate-300 bg-white text-[10px] font-semibold text-slate-600 hover:bg-slate-50" onClick={() => annotateNow([c])} title="Annotate this call (genes, PanelApp, ClinGen, gnomAD, DGV Gold)">annotate</button></td></>
                  : <><GenesCell c={c} v={v} />
                  <PanelCell v={v} />
                  <DosageCell c={c} v={v} />
                  <td className="px-1.5 py-0.5 whitespace-nowrap"><PriorityBadge v={v} /></td></>}
                  <td className="px-1.5 py-0.5 whitespace-nowrap sticky right-[72px] z-[1] bg-inherit shadow-[-6px_0_6px_-6px_rgba(15,23,42,0.18)]" onClick={e => e.stopPropagation()}>
                    <select value={inList ? cls : ''} aria-label={`Classification of ${c.chrom}:${c.start}-${c.end}`} style={{ borderLeftColor: inList ? col ?? '#d3dbe5' : '#d3dbe5' }}
                      title={!inList && cls ? `The file says “${cls}”: kept in the export unless a class is chosen` : 'The class written in the export'}
                      onChange={e => cgh.setClass(c, e.target.value)}
                      className={`h-[22px] w-[128px] pl-1 pr-0.5 rounded border border-slate-200 border-l-[4px] bg-white text-[11px] outline-none cursor-pointer ${inList ? 'font-medium text-slate-800' : 'text-slate-400'}`}>
                      <option value="">— not classified</option>
                      {CGH_CLASSES.map(o => <option key={o.value} value={o.value}>{o.hint && o.value === 'NS' ? `NS (${o.hint})` : o.value}</option>)}
                    </select>
                    {!cls.trim() && lab.has(c.id) && <button onClick={() => cgh.setClass(c, 'Benign')} className="ml-1 px-1 rounded border border-slate-300 bg-white text-[10px] font-semibold text-slate-600 hover:bg-slate-50" title="Suggestion only, from the laboratory's benign list: sets the class Benign on this call. Nothing is classified by itself.">→ Benign?</button>}
                    <button onClick={() => { setNoteFocus(c.id); setFolded(false); if (!on) show(i); }} className={`ml-1 px-1 rounded ${note ? 'text-accent-700 bg-accent-50 font-bold' : 'text-slate-300 hover:text-slate-600'}`}
                      title={note ? `Note: ${note}` : 'Add a note on this call (exported in the Comment column)'} aria-label="Note">✎</button>
                  </td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap sticky right-0 z-[1] bg-inherit w-[72px] min-w-[72px]" onClick={e => e.stopPropagation()}>
                    {hub ? <a href={hub} target="_blank" rel="noopener noreferrer" className="text-accent-700 hover:underline font-medium" title={`CNV-Hub (cnvhub.net, Dijon University Hospital): this call's page, ${hub}. Opens in a new tab; the call's position and type are in the address.`}>CNV-Hub ↗</a>
                      : <span className="text-slate-300" title="CNV-Hub pages are known for hg19 (GRCh37) only">CNV-Hub</span>}
                  </td>
                </tr>,
                on && !folded && <tr key={`${c.id}-detail`}><td colSpan={13} className="p-0 border-x-2 border-b-2 border-accent-600"><CnvDetail c={c} v={v} onRetry={srcs => retryNow([c], srcs)} lab={lab.get(c.id)} hints={hints.get(c.id) ?? []} sex={s.sex[c.sampleKey]} series={series.get(c.id)} note={note} focusNote={noteFocus === c.id} onNoteFocused={() => setNoteFocus(null)} /></td></tr>,
              ];
            })}
            {!rows.length && <tr><td colSpan={13} className="px-3 py-3 text-slate-400">No call passes the filter.</td></tr>}
          </tbody>
        </table>
      </div>
      {s.files.some(fl => fl.warnings.some(w => !/^(BED|VCF): no probe/.test(w))) && <div className="px-5 text-[10.5px] text-amber-700">{s.files.flatMap(fl => fl.warnings.filter(w => !/^(BED|VCF): no probe/.test(w)).map(w => `${fl.name}: ${w}`)).join(' · ')}</div>}
      <details className="px-5 pb-2 text-[10.5px] text-slate-500 leading-4">
        <summary className="cursor-pointer text-accent-700 font-medium select-none">ⓘ How the calls are read, annotated, ranked and exported</summary>
        Click a call to draw it below as a deletion or a duplication (↑ ↓ for the next one), with its genes listed under its row (click it again to fold them): the genes, DGV Gold, ClinVar, AChro-Puce and PanelApp tracks of the browser are compared with it.
        Genes are RefSeq's NM_ models of the interval (UCSC ncbiRefSeqCurated, the file's build), the file's Gene Name only when they cannot be read (then marked <span className="text-amber-700">file</span>); ClinGen dosage through the UCSC Genome Browser, gnomAD v2.1.1 constraint on GRCh37, PanelApp from {PANEL_SOURCE.label} (every panel); DGV Gold marked when a record of the same type shares {Math.round(DGV_COVER * 100)} % of the call or more (the call inside it: 100 %), or lies inside the call whatever its share; one icon per way they lie, with the best share and how many (×n): <DgvIcon rel="within" /> the call inside a DGV CNV, <DgvIcon rel="contains" /> the call holding one, <DgvIcon rel="overlaps" /> overlapping on one side.
        A gain is drawn as a tandem duplication: array-CGH gives the copy number, not where the extra copy lies; a loss of both copies (x0) as homozygous.
        Positions and ISCN as in the file. Classes and notes are kept in this tab only, never in the browser's storage: <b>Export classified table</b> writes each file back with them (Classification column, and a Comment column; a VCF in INFO; an .xlsx as tab-separated text).
        <PriorityRules />
      </details>
    </div>
  );
}

// ---- the columns of the genes' evidence ----
const MoiBadge = ({ m }: { m: MoiClass }) => <span className="inline-flex items-center rounded-full px-1 text-[9px] font-bold border leading-[13px]" style={{ color: MOI[m].color, borderColor: `${MOI[m].color}88`, background: `${MOI[m].color}1a` }} title={MOI[m].label}>{MOI[m].short}</span>;
 /** a listed call never asked (past the automatic set, or under the minimums): grey, annotated on demand */
const NotAnnCell = () => <td className="px-1.5 py-0.5 text-slate-400 whitespace-nowrap" data-cgh-notann-cell>not annotated</td>;
/** the subtle mark of a call under the minimums, shown on demand: what it lacks */
function BelowTag({ c, minProbes, minKb }: { c: CghCnv; minProbes: number; minKb: number }) {
  const few = c.probes != null && c.probes < minProbes;
  return <span className="ml-1 inline-flex items-center rounded px-1 text-[9.5px] border border-slate-200 bg-slate-50 text-slate-500 leading-[14px] align-middle" title="Under the minimums: listed on demand">{few ? `${c.probes} probe${c.probes === 1 ? '' : 's'}` : `< ${sizeWord(minKb)}`}</span>;
}
const waitingDot = (n: number) => (n ? <span className="text-slate-400" title={`PanelApp: ${n} gene${n > 1 ? 's' : ''} still asked`}> …</span> : null);
type View = ReturnType<typeof viewOf>;

// ---- manual retry: per column, per patient, only what failed ----
type RetryCol = 'genes' | 'panelapp' | 'dosage' | 'dgv' | 'bands';
/** the sources each column asks again: Dosage holds ClinGen genes + gnomAD and the ClinGen regions; Type the DGV Gold medallion */
const RETRY_SRC: Record<RetryCol, (CghSource | 'panelapp')[]> = { genes: ['genes'], panelapp: ['panelapp'], dosage: ['dosage', 'regions'], dgv: ['dgv'], bands: ['bands'] };
const RETRY_WORD: Record<RetryCol, string> = { genes: 'RefSeq genes', panelapp: 'PanelApp', dosage: 'ClinGen dosage (genes and regions) and gnomAD constraint', dgv: 'DGV Gold', bands: 'cytogenetic bands' };
const SOURCE_WORD: Record<CghSource, string> = { genes: 'RefSeq genes', dosage: 'ClinGen dosage and gnomAD constraint', regions: 'ClinGen dosage regions', dgv: 'DGV Gold', bands: 'cytogenetic bands' };
interface RetryGroup { key: string; name: string; n: number; asking: number; run: () => void }
const RetryLink = ({ onClick }: { onClick: () => void }) => <button className="text-accent-700 hover:underline font-medium" onClick={onClick} title="Ask this source again for this call">retry</button>;
const RETRY_BTN = 'ml-1 px-1 rounded border border-amber-300 bg-amber-50 text-[10px] font-semibold text-amber-800 hover:bg-amber-100 align-middle';
/**
 * A column's retry: one patient listed, a button; several, a small menu with one button per patient (and all of them),
 * so that only that patient's failed calls (or genes) are asked again. Hidden when nothing failed.
 */
function RetryCtl({ col, groups }: { col: RetryCol; groups: RetryGroup[] }) {
  const [at, setAt] = useState<{ x: number; y: number } | null>(null);
  const failed = groups.filter(g => g.n > 0), n = failed.reduce((k, g) => k + g.n, 0), asking = groups.reduce((k, g) => k + g.asking, 0);
  useEffect(() => {
    if (!at) return;
    const off = (e: Event) => { if (!(e.target as HTMLElement).closest?.('[data-cgh-retry-menu]')) setAt(null); };
    document.addEventListener('mousedown', off); return () => document.removeEventListener('mousedown', off);
  }, [at]);
  if (!n) return asking ? <span className="ml-1 text-[10px] text-slate-400 font-normal" title={`${RETRY_WORD[col]}: asked again`}>↻ …</span> : null;
  const what = col === 'panelapp' ? `gene${n > 1 ? 's' : ''}` : `call${n > 1 ? 's' : ''}`;
  const title = `${RETRY_WORD[col]} did not answer for ${n} ${what}${groups.length > 1 ? ` (${failed.length} patient${failed.length > 1 ? 's' : ''})` : ''}: ask again only those`;
  if (groups.length === 1) return <button className={RETRY_BTN} data-cgh-retry={col} title={title} onClick={e => { e.stopPropagation(); failed[0].run(); }}>↻ retry ({n})</button>;
  return <>
    <button className={RETRY_BTN} data-cgh-retry={col} title={title} aria-haspopup="true" aria-expanded={!!at}
      onClick={e => { e.stopPropagation(); const r = e.currentTarget.getBoundingClientRect(); setAt(at ? null : { x: Math.min(r.left, window.innerWidth - 270), y: r.bottom + 4 }); }}>↻ retry ({n}) ▾</button>
    {/* fixed: the table's box scrolls and would clip it */}
    {at && <div data-cgh-retry-menu={col} role="dialog" className="fixed z-50 w-[260px] bg-white border border-slate-200 rounded-lg shadow-lg p-1.5 text-[11px] font-normal text-slate-700 whitespace-normal" style={{ left: at.x, top: at.y }}>
      <div className="px-1 pb-1 text-slate-500">{RETRY_WORD[col]}: ask again</div>
      {failed.map(g => <div key={g.key} className="flex items-center gap-2 px-1 py-0.5">
        <span className="flex-1 truncate font-semibold">{g.name}</span><span className="text-slate-500 tabular-nums">{g.n} {col === 'panelapp' ? `gene${g.n > 1 ? 's' : ''}` : `call${g.n > 1 ? 's' : ''}`}</span>
        <button className={RETRY_BTN} onClick={() => { g.run(); if (failed.length === 1) setAt(null); }}>↻ retry</button></div>)}
      {failed.length > 1 && <div className="flex justify-end px-1 pt-1 border-t border-slate-100 mt-1"><button className={RETRY_BTN} onClick={() => { failed.forEach(g => g.run()); setAt(null); }}>↻ all ({n})</button></div>}
    </div>}
  </>;
}

/** the small amber tag of a call whose log ratio says a share of the cells only carries it */
function MosaicTag({ c, sex }: { c: CghCnv; sex: CghSex | undefined }) {
  const m = mosaicOf(c, sex);
  if (!m?.ok || m.f >= MOSAIC_TAG_BELOW) return null;
  return <span className="ml-1 inline-flex items-center rounded-full px-1.5 text-[9.5px] font-bold border border-slate-300 bg-white text-slate-600 leading-[15px] align-middle" title="Fraction of cells estimated from the mean log ratio (AChro-Puce calculator, July 2025): indicative, see the call's detail">mosaic ~{Math.round(m.f * 100)}%</span>;
}
/** the call's detail line: expected and observed log ratio, and the fraction of cells it gives */
function MosaicLine({ c, sex }: { c: CghCnv; sex: CghSex | undefined }) {
  const m = mosaicOf(c, sex);
  if (!m) return null;
  const pc = (x: number) => `${Math.round(x * 100)} %`;
  return <div className="mb-1 text-[11px]"><span className="font-semibold text-slate-900">Log ratio</span>{' '}
    {m.ok ? <>observed <b>{c.logRatio!.toFixed(2)}</b>, expected for a one-copy {c.kind}: <b>{Number.isFinite(m.expected) ? `${m.expected > 0 ? '+' : ''}${m.expected.toFixed(2)}` : 'not finite'}</b> · cells carrying it: <b>{pc(m.f)}</b></>
      : <>observed <b>{c.logRatio!.toFixed(2)}</b> · fraction of cells: <span className="text-amber-700">sex needed (chrX / chrY outside the PARs)</span></>}
    <span className="text-slate-500"> · indicative (AChro-Puce calculator, July 2025): no biological variability</span></div>;
}
const seriesTip = (i: SeriesInfo) => i.hits.length ? i.hits.map(h => `${h.name}: ${Math.round(h.pct * 100)} % overlap`).join('\n') : 'No other patient has this call';
/** other patients carrying the call / other patients loaded; empty when only one patient is loaded */
function SeriesCell({ info }: { info: SeriesInfo | undefined }) {
  if (!info || !info.others) return <td />;
  return <td className={`px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap ${info.hits.length ? 'font-semibold text-slate-800' : 'text-slate-400'}`} title={seriesTip(info)}>{info.hits.length}/{info.others}</td>;
}
function GenesCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-1.5 py-0.5 text-slate-400">…</td>;
  const gs = v.ann.genes, file = v.ann.source === 'file';
  return (
    <td className="px-1.5 py-0.5 max-w-[120px] truncate" title={`${file ? `From the file (RefSeq not read): ${c.genes.join(', ')}` : `${gs.length} RefSeq gene${gs.length > 1 ? 's' : ''}: ${gs.map(g => `${g.symbol}${g.whole ? '' : ' (cut)'}`).join(', ')}`}`}>
      {file && <span className="text-amber-700 font-semibold mr-1">file</span>}
      {gs.length ? <><span className="tabular-nums text-slate-500 mr-1">{gs.length}</span><i>{gs.slice(0, 4).map(g => g.symbol).join(', ')}</i>{gs.length > 4 ? <span className="text-slate-400"> +{gs.length - 4}</span> : null}</>
        : <span className="text-slate-400">{file ? '—' : 'no gene'}</span>}
    </td>
  );
}
/** PanelApp's genes of the call, green first, then amber, then red (each with the modes of its best-rated entries) */
const RATING = { 3: { word: 'green', color: '#16a34a' }, 2: { word: 'amber', color: '#d97706' }, 1: { word: 'red', color: '#dc2626' } } as const;
const modesAt = (x: GeneView) => (x.pa!.best === 3 ? x.pa!.green : x.pa!.best === 2 ? x.pa!.amber : []).filter(m => m !== '?');
function PanelCell({ v }: { v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-1.5 py-0.5 text-slate-400">…</td>;
  const by = (r: 1 | 2 | 3) => v.genes.filter(x => x.pa?.best === r);
  const green = by(3), amber = by(2), red = by(1);
  const title = ([3, 2, 1] as const).flatMap(r => by(r).map(x => `${RATING[r].word}: ${x.g.symbol}${modesAt(x).length ? ` (${modesAt(x).map(m => MOI[m].short).join(', ')})` : ''}${x.pa!.diseases.length ? ` · ${x.pa!.diseases.slice(0, 2).join('; ')}` : ''}`)).join('\n');
  return (
    <td className="px-1.5 py-0.5 max-w-[150px] overflow-hidden" title={title || undefined}>
      <span className="flex items-center gap-x-1.5 whitespace-nowrap">
        {green.slice(0, 2).map(x => <span key={x.g.symbol} className="inline-flex items-center gap-0.5"><i>{x.g.symbol}</i>{modesAt(x).map(m => <MoiBadge key={m} m={m} />)}</span>)}
        {green.length > 2 && <span className="text-[10px] font-semibold" style={{ color: RATING[3].color }} title={`${green.length - 2} more green genes`}>+{green.length - 2}</span>}
        {amber.length > 0 && <span className="text-[10px] font-semibold" style={{ color: RATING[2].color }}>{green.length ? '' : <i className="font-normal">{amber[0].g.symbol} </i>}{amber.length} amber</span>}
        {red.length > 0 && <span className="text-[10px] font-semibold" style={{ color: RATING[1].color }}>{red.length} red</span>}
        {!green.length && !amber.length && !red.length && !v.waiting && <span className="text-slate-300">—</span>}
        {waitingDot(v.waiting)}
      </span>
    </td>
  );
}
/** the regions whose score for the call's type (HI for a loss, TS for a gain) is a grade, 1-3: not 0, 30 (AR) or 40 (unlikely), which is noise here */
const regionsScored = (c: CghCnv, regs: RegionHit[]) => regs.filter(r => { const x = c.kind === 'loss' ? r.hi : r.ts; return x === 1 || x === 2 || x === 3; });
const labText = (h: LabHit) => `${h.entry.name || 'unnamed entry'} · ${h.entry.chrom}:${h.entry.start.toLocaleString('en-US')}-${h.entry.end.toLocaleString('en-US')} (${h.entry.type === 'both' ? 'both types' : h.entry.type}) · ${Math.round(h.pct * 100)} % of the call`;
/** violet-grey tag: the call lies in the laboratory's benign list (distinct from DGV's violet pill) */
function LabTag({ hit }: { hit: LabHit | undefined }) {
  if (!hit) return null;
  return <span className="ml-1 inline-flex items-center rounded px-1 text-[9.5px] font-semibold align-middle bg-white text-slate-600 border border-slate-300" title={`Laboratory's benign list: ${labText(hit)}`}>lab benign</span>;
}
/** "Lab benign CNVs (BED)…": loads the laboratory's list; kept in this browser (localStorage), × removes it */
function LabBenignButton() {
  const s = useCgh(), lab = s.lab;
  const [msg, setMsg] = useState('');
  const load = async (f: File) => {
    try {
      const r = parseLabBenign(decodeText(new Uint8Array(await f.arrayBuffer())), f.name);
      if (!r.entries.length) { setMsg(`${f.name}: no interval read (${r.skipped} line${r.skipped === 1 ? '' : 's'} skipped)`); return; }
      setMsg(''); cgh.setLab(r);
    } catch (e) { setMsg(`${f.name}: ${(e as Error).message}`); }
  };
  return (
    <span className="inline-flex items-center gap-1">
      <label className="inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border border-slate-200 bg-white text-xs font-medium text-slate-700 cursor-pointer whitespace-nowrap hover:border-slate-300 hover:bg-slate-50"
        title={"The laboratory's list of benign CNVs: a BED (chrom, start 0-based, end, optional name that may say del / dup / loss / gain; without it the entry applies to both types). A call sharing " + Math.round(LAB_COVER * 100) + " % of itself with an entry of its type is marked and can be suggested as Benign.\nA laboratory reference, not patient data: kept in this browser for the next sessions."}>
        <Icon name="upload" size={13} />Lab benign CNVs (BED)…
        <input data-lab-input type="file" accept=".bed,.txt,.tsv,text/plain" className="hidden" onChange={e => { const f = e.target.files?.[0]; if (f) void load(f); e.target.value = ''; }} />
      </label>
      {lab && <span className="inline-flex items-center gap-1 h-[26px] px-2 rounded-lg bg-slate-100 border border-slate-300 text-[11px] text-slate-700" data-lab-chip
        title={`${lab.fileName}: ${lab.read} interval${lab.read === 1 ? '' : 's'} read, ${lab.skipped} line${lab.skipped === 1 ? '' : 's'} skipped${s.labStored ? '. Kept in this browser.' : '. Too large to keep: it is lost when the page closes.'}`}>
        {lab.fileName} · {lab.read} read{lab.skipped ? `, ${lab.skipped} skipped` : ''}{s.labStored ? '' : ' (not kept)'}
        <button onClick={() => cgh.setLab(null)} className="text-slate-500 hover:text-red-600" aria-label="Remove the lab benign list" title="Remove the list (also from this browser)">×</button></span>}
      {msg && <span className="text-[11px] text-amber-700">{msg}</span>}
    </span>
  );
}
/** a neutral tag in the Band cell: terminal, whole arm, whole chromosome; the rearrangement hints in its tooltip */
function TopoTag({ c, v, hints }: { c: CghCnv; v: View; hints: string[] | undefined }) {
  const t = v.ann?.topo, w = topoLabel(t);
  if (!w) return null;
  return <span className="ml-1 inline-flex items-center rounded border border-slate-300 px-1 text-[9.5px] font-medium text-slate-500" title={[...topoWords(t, c.chrom), ...(hints ?? [])].join('\n')}>{w}{hints?.length ? ' ⚠' : ''}</span>;
}
function DosageCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-1.5 py-0.5 text-slate-400">…</td>;
  const key = c.kind === 'loss' ? 'hi' : 'ts', word = c.kind === 'loss' ? 'HI' : 'TS';
  const scored = v.genes.filter(x => x[key] === 3 || x[key] === 2).sort((a, b) => (b[key] ?? 0) - (a[key] ?? 0));
  const regs = regionsScored(c, v.ann.regions);
  const low = v.genes.filter(x => x.loeuf != null).sort((a, b) => a.loeuf! - b.loeuf!)[0];
  return (
    <td className="px-1.5 py-0.5 whitespace-nowrap">
      {scored.slice(0, 1).map(x => <span key={x.g.symbol} className="mr-1.5 font-semibold" style={{ color: x[key] === 3 ? '#b91c1c' : '#c2410c' }} title={`ClinGen ${c.kind === 'loss' ? 'haploinsufficiency' : 'triplosensitivity'} score ${x[key]} (${x[key] === 3 ? 'sufficient' : 'some'} evidence)`}>{word}{x[key]} <i className="font-normal">{x.g.symbol}</i></span>)}
      {scored.length > 1 && <span className="text-slate-500 mr-1.5" title={scored.slice(1).map(x => `${word}${x[key]} ${x.g.symbol}`).join(', ')}>+{scored.length - 1}</span>}
      {regs.map(r => <span key={r.id} className="mr-1.5 font-semibold" style={{ color: r[key] === 3 ? '#b91c1c' : '#c2410c' }} title={`${r.label} (${r.id}): ClinGen region, ${word} ${r[key]}, ${Math.round(r.covered * 100)} % of it in the call`}>{word}{r[key]} region{r.covered < 0.999 ? ` ${Math.min(99, Math.round(r.covered * 100))} %` : ''}</span>)}
      {low && <span className="text-slate-600" title={`Lowest gnomAD LOEUF of the call (${low.version ?? ''}): ${low.g.symbol} ${low.loeuf!.toFixed(2)}${low.pli != null ? `, pLI ${low.pli.toFixed(2)}` : ''}`}>LOEUF {low.loeuf!.toFixed(2)}</span>}
      {!scored.length && !regs.length && !low && <span className="text-slate-300">—</span>}
    </td>
  );
}
function PriorityBadge({ v }: { v: View }) {
  if (!v.score) return <span className="text-slate-400">…</span>;
  const st = LEVEL_STYLE[v.score.level], incomplete = v.missing.length > 0 || v.waiting > 0 || v.asking > 0;
  const title = [...v.score.reasons.map(r => `${LEVEL_STYLE[r.level].label}: ${r.text}`), ...(v.waiting ? [`PanelApp still asked for ${v.waiting} gene${v.waiting > 1 ? 's' : ''}`] : []), ...v.missing.map(m => `Not available: ${m}`)].join('\n') || 'No protein-coding gene in the call';
  return <span className={`inline-flex items-center rounded px-1.5 text-[10.5px] font-bold ${incomplete && v.score.level !== 'high' ? 'border border-dashed' : ''}`} style={{ color: st.color, background: st.bg, borderColor: st.color }} title={title}>{st.label}{incomplete && v.score.level !== 'high' ? ' ?' : ''}</span>;
}

/** under the call's row: why its priority, the AChro-Puce regions, and each gene with its scores */
function CnvDetail({ c, v, onRetry, lab, hints, sex, series, note, focusNote, onNoteFocused }: { c: CghCnv; v: View; onRetry: (srcs: (CghSource | 'panelapp')[]) => void; lab: LabHit | undefined; hints: string[]; sex: CghSex | undefined; series: SeriesInfo | undefined; note: string; focusNote: boolean; onNoteFocused: () => void }) {
  const [allGenes, setAllGenes] = useState(false);
  const noteBox = <textarea value={note} onChange={e => cgh.setNote(c, e.target.value)} rows={1} placeholder="Note on this call (exported in the Comment column)"
    ref={el => { if (el && focusNote) { el.focus(); onNoteFocused(); } }}
    onKeyDown={e => e.stopPropagation()} aria-label="Note on this call"
    className="w-full mb-1 px-2 py-1 rounded border border-slate-300 bg-white text-[11.5px] text-slate-800 resize-y outline-none focus:border-accent-400 focus:ring-2 focus:ring-accent-100" />;
  const seen = series && series.others > 0 ? <div className="mb-1 text-[11px]"><span className="font-semibold text-slate-900">Series</span> {series.hits.length ? <>found in {series.hits.length} of {series.others} other patient{series.others > 1 ? 's' : ''}: {series.hits.map(h => `${h.name} (${Math.round(h.pct * 100)} %)`).join(', ')}</> : <>in none of the {series.others} other patient{series.others > 1 ? 's' : ''}</>}<span className="text-slate-500"> · same type, reciprocal overlap ≥ {Math.round(SERIES_OVERLAP * 100)} %</span></div> : null;
  const mosaic = <>{seen}<MosaicLine c={c} sex={sex} /></>;
  if (!v.ann || v.ann.state === 'loading' || !v.score) return <div className="px-4 py-2 bg-slate-50 text-slate-500 text-[11px]" onClick={e => e.stopPropagation()}>{noteBox}{mosaic}Reading the genes of the call…</div>;
  const ann = v.ann, key = c.kind === 'loss' ? 'hi' : 'ts';
  // PanelApp green first, then amber, then red, then the others; within each, the dosage score, then the position
  const genes = v.genes.slice().sort((a, b) => (b.pa?.best ?? 0) - (a.pa?.best ?? 0) || (b[key] ?? 0) % 10 - (a[key] ?? 0) % 10 || a.g.start - b.g.start);
  const shownGenes = allGenes ? genes : genes.slice(0, 8);
  const dot = (l: Level) => <span className="inline-block w-2 h-2 rounded-full mr-1.5 align-middle" style={{ background: LEVEL_STYLE[l].color }} />;
  return (
    <div className="sticky left-0 w-[calc(100vw-5.5rem)] max-w-full px-3 py-1.5 bg-slate-50 text-[11px] leading-4 text-slate-700 cursor-default" onClick={e => e.stopPropagation()}>
      {noteBox}
      {mosaic}
      {lab && <div className="mb-1 text-[11px]"><span className="font-semibold text-slate-900">Lab benign list</span> <span className="text-violet-800">{labText(lab)}</span> <span className="text-slate-500">· a suggestion: the class is yours to give</span></div>}
      {(topoWords(v.ann.topo, c.chrom).length > 0 || hints.length > 0) && <div className="mb-1 text-[11px]"><span className="font-semibold text-slate-900">Position</span> {topoWords(v.ann.topo, c.chrom).join(' · ')}{hints.map(h => <div key={h} className="text-violet-800">{h}</div>)}</div>}
      <div className="grid gap-x-6 gap-y-1 grid-cols-1 min-[1100px]:grid-cols-2">
        <div>
          <div className="font-semibold text-slate-900">Priority: <span style={{ color: LEVEL_STYLE[v.score.level].color }}>{LEVEL_STYLE[v.score.level].label}</span> <span className="font-normal text-slate-500">(a proposal for the review, not a class)</span></div>
          {v.score.reasons.length ? v.score.reasons.map(r => <div key={r.id}>{dot(r.level)}{r.text}</div>) : <div className="text-slate-500">No protein-coding gene in the call.</div>}
          {v.waiting > 0 && <div className="text-slate-500">PanelApp: {v.waiting} gene{v.waiting > 1 ? 's' : ''} still asked…</div>}
          {CGH_SOURCES.filter(x => x !== 'dgv' && (ann.failed[x] || ann.asking.includes(x))).map(x => <div key={x} className="text-amber-700" data-cgh-failed={x}>
            {ann.asking.includes(x) ? <span className="text-slate-500">{SOURCE_WORD[x]}: asked again…</span> : <>Not available: {ann.failed[x]} · <RetryLink onClick={() => onRetry([x])} /></>}</div>)}
          {v.paFailed > 0 && <div className="text-amber-700" data-cgh-failed="panelapp">Not available: PanelApp ({v.paFailed} gene{v.paFailed > 1 ? 's' : ''} not answered) · <RetryLink onClick={() => onRetry(['panelapp'])} /></div>}
          {v.missing.filter(m => /^PanelApp asked for/.test(m)).map(m => <div key={m} className="text-amber-700">Not available: {m}</div>)}
        </div>
        {(v.ann.dgv.length > 0 || v.ann.dgvError) && <div>
          <div className="font-semibold text-slate-900">DGV Gold (same type: ≥ {Math.round(DGV_COVER * 100)} % of the call shared, or inside the call whatever its share)</div>
          {dgvShown(v.ann.dgv).map(x => <div key={x.id} className="flex flex-wrap items-center gap-x-1"><DgvIcon rel={x.rel} />{x.id} <span className="text-slate-500">({x.type})</span>: <b>{pctText(x.cover)}</b> <span className="text-slate-500">· {DGV_REL[x.rel].word}</span>{x.freq ? ` · frequency ${x.freq}` : ''}{x.samples != null ? ` · ${x.samples} samples` : ''}</div>)}
          {v.ann.dgv.length > dgvShown(v.ann.dgv).length && <div className="text-slate-500">+{v.ann.dgv.length - dgvShown(v.ann.dgv).length} more DGV CNVs inside the call (in the medallion's tooltip)</div>}
          {v.ann.asking.includes('dgv') ? <div className="text-slate-500">asked again…</div>
            : v.ann.dgvError && <div className="text-amber-700" data-cgh-failed="dgv">not available ({v.ann.dgvError}) · <RetryLink onClick={() => onRetry(['dgv'])} /></div>}
        </div>}
        {regionsScored(c, v.ann.regions).length > 0 && <div>
          <div className="font-semibold text-slate-900">ClinGen dosage regions</div>
          {regionsScored(c, v.ann.regions).map(r => <div key={r.id}>{r.label} <span className="text-slate-500">({r.id})</span>: HI <b>{r.hi ?? '–'}</b>, TS <b>{r.ts ?? '–'}</b> <span className="text-slate-500">({Math.round(r.covered * 100)} % of the region in the call)</span></div>)}
        </div>}
        {v.ann.achro.length > 0 && <div>
          <div className="font-semibold text-slate-900">AChro-Puce</div>
          {v.ann.achro.map(a => <div key={a.id}>{a.name}: <b>{a.cls.join(' / ')}</b> <span className="text-slate-500">({a.how})</span></div>)}
        </div>}
      </div>
      {genes.length > 0 && <table className="mt-1.5 w-full border-collapse bg-white border border-slate-200 rounded">
        <thead className="text-slate-500 text-left"><tr className="border-b border-slate-200">
          <th className="px-1.5 py-0.5 font-medium">Gene</th><th className="px-1.5 py-0.5 font-medium">In the call</th>
          <th className="px-1.5 py-0.5 font-medium" title="ClinGen haploinsufficiency score">HI</th><th className="px-1.5 py-0.5 font-medium" title="ClinGen triplosensitivity score">TS</th>
          <th className="px-1.5 py-0.5 font-medium">LOEUF (pLI)</th><th className="px-1.5 py-0.5 font-medium">PanelApp</th><th className="px-1.5 py-0.5 font-medium">Diseases (PanelApp)</th>
        </tr></thead>
        <tbody>{shownGenes.map(x => (
          <tr key={x.g.symbol} className="border-b border-slate-100">
            <td className="px-1.5 py-0.5"><i className="font-semibold">{x.g.symbol}</i> <span className="text-slate-400">{x.g.tx}</span></td>
            <td className="px-1.5 py-0.5">{v.ann!.source === 'file' ? <span className="text-amber-700" title="RefSeq could not be read: the file names the gene, its extent is not known here">from the file</span> : <InCall g={x.g} />}</td>
            <td className="px-1.5 py-0.5 font-semibold" style={{ color: x.hi === 3 ? '#b91c1c' : x.hi === 2 ? '#c2410c' : undefined }}>{dosageText(x.hi)}</td>
            <td className="px-1.5 py-0.5 font-semibold" style={{ color: x.ts === 3 ? '#b91c1c' : x.ts === 2 ? '#c2410c' : undefined }}>{dosageText(x.ts)}</td>
            <td className="px-1.5 py-0.5 tabular-nums">{x.loeuf != null ? x.loeuf.toFixed(2) : ''}{x.pli != null ? <span className="text-slate-400"> ({x.pli.toFixed(2)})</span> : null}</td>
            <td className="px-1.5 py-0.5 whitespace-nowrap">{!x.pa ? <span className={x.paState === 'failed' ? 'text-amber-700' : 'text-slate-400'}>{x.paState === 'failed' ? 'not answered' : x.paState === 'not asked' ? `not asked (over ${PANELAPP_MAX_GENES} genes)` : '…'}</span>
              : x.pa.best === 0 ? <span className="text-slate-400">in no panel</span>
              : <span className="inline-flex items-center gap-1"><span className="font-semibold" style={{ color: RATING[x.pa.best].color }}>{RATING[x.pa.best].word}</span>
                {modesAt(x).map(m => <MoiBadge key={m} m={m} />)}</span>}</td>
            <td className="px-1.5 py-0.5 text-slate-600 max-w-[380px] truncate" title={x.pa?.diseases.join('\n')}>{x.pa?.diseases.slice(0, 2).join(' · ')}</td>
          </tr>))}</tbody>
      </table>}
      {genes.length > 8 && <button className="mt-1 text-accent-700 hover:underline" onClick={() => setAllGenes(!allGenes)}>{allGenes ? 'fewer genes' : `all ${genes.length} genes`}</button>}
    </div>
  );
}

/** the rules of the priority, as the score applies them (cghScore.ts) */
function PriorityRules() {
  return (
    <details className="mt-1">
      <summary className="cursor-pointer text-accent-700 font-medium">The priority's rules (HIGH, MODERATE, LOW)</summary>
      <div className="grid grid-cols-[1fr_1fr] gap-4 mt-1">
        {(['loss', 'gain'] as const).map(k => <div key={k}>
          <div className="font-semibold text-slate-700">{k === 'loss' ? 'Losses' : 'Gains'}</div>
          {SCORE_RULES.filter(r => r.kind === k).map(r => <div key={r.id}><span className="font-bold" style={{ color: LEVEL_STYLE[r.level].color }}>{LEVEL_STYLE[r.level].label}</span> <span className="text-slate-400">{r.id}</span> {r.label}</div>)}
          <div><span className="font-bold" style={{ color: LEVEL_STYLE.low.color }}>LOW</span> anything else, a call without a protein-coding gene included</div>
        </div>)}
      </div>
      <div className="mt-1">A call takes the highest level among the rules it meets. Evidence weighed as in the ACMG / ClinGen technical standards for CNVs (Riggs et al., Genet Med 2020;22:245-57: sections 2 and 3); a dashed badge with ? is incomplete (a source did not answer, PanelApp still asked).</div>
    </details>
  );
}

/** the DGV CNVs listed under a call: those holding or crossing it, then the 8 sharing most of those inside it */
const dgvShown = <T extends { rel: DgvRelation }>(d: T[]): T[] => [...d.filter(x => x.rel !== 'contains'), ...d.filter(x => x.rel === 'contains').slice(0, 4)];
/** a ClinGen dosage score: 30 and 40 are not grades (30: autosomal recessive gene, 40: dosage sensitivity unlikely) */
const dosageText = (v: number | null) => (v == null ? '' : v === 30 ? '30 (AR)' : v === 40 ? '40 (unlikely)' : String(v));
/** a share: under 1 % said so rather than 0 % */
const pctText = (f: number) => (f > 0 && f < 0.005 ? '<1%' : `${Math.round(f * 100)}%`);
/** how a DGV Gold CNV lies against the call, drawn: the call (dark) and the DGV CNV (violet) */
const DGV_REL: Record<DgvRelation, { word: string; call: [number, number]; dgv: [number, number] }> = {
  within: { word: 'the call lies inside a DGV Gold CNV', call: [7, 15], dgv: [1, 21] },
  contains: { word: 'the call holds a DGV Gold CNV', call: [1, 21], dgv: [6, 16] },
  overlaps: { word: 'the call and a DGV Gold CNV overlap on one side', call: [1, 14], dgv: [8, 21] },
};
const DgvIcon = ({ rel }: { rel: DgvRelation }) => (
  <svg width="22" height="11" viewBox="0 0 22 11" aria-hidden="true" className="shrink-0">
    <rect x={DGV_REL[rel].dgv[0]} y="1" width={DGV_REL[rel].dgv[1] - DGV_REL[rel].dgv[0]} height="3.5" rx="1" fill="#7c3aed" />
    <rect x={DGV_REL[rel].call[0]} y="6.5" width={DGV_REL[rel].call[1] - DGV_REL[rel].call[0]} height="3.5" rx="1" fill="#1f2937" />
  </svg>
);
/** the DGV Gold CNVs sharing 80 % of the call or more: one medallion per way they lie against it, with the best share */
function DgvMedal({ v }: { v: View }) {
  const d = v.ann?.dgv;
  if (!d?.length) return null;
  const rels = (['within', 'contains', 'overlaps'] as DgvRelation[]).map(rel => ({ rel, hs: d.filter(x => x.rel === rel) })).filter(r => r.hs.length);
  const title = rels.map(({ rel, hs }) => `${DGV_REL[rel].word} (${rel === 'contains' ? 'whatever its share of the call' : `≥ ${Math.round(DGV_COVER * 100)} % of the call shared`}): ${hs.length}\n${hs.slice(0, 12).map(x => `  ${x.id} (${x.type}): ${pctText(x.cover)} of the call${x.freq ? `, frequency ${x.freq}` : ''}${x.samples != null ? `, ${x.samples} samples` : ''}`).join('\n')}${hs.length > 12 ? `\n  +${hs.length - 12} more` : ''}`).join('\n');
  // one pill: each way they lie, with its best share
  return <span className="ml-1 whitespace-nowrap inline-flex items-center gap-1 rounded-full px-1.5 text-[9.5px] font-bold border border-slate-300 bg-slate-50 text-slate-600 leading-[15px] align-middle" title={`DGV Gold, same type:\n${title}`}>
    {rels.map(({ rel, hs }) => <span key={rel} className="inline-flex items-center gap-0.5"><DgvIcon rel={rel} />{pctText(Math.max(...hs.map(x => x.cover)))}</span>)}
  </span>;
}

const THEORETICAL = "theoretical: the call's first and last probes, the true breakpoints lie between probes";
/** what the call holds of a gene: whole, or the exons of its principal transcript it cuts and the frame */
function InCall({ g }: { g: CnvGene }) {
  if (g.whole) return <>whole</>;
  const c = g.cut;
  if (c?.effect === 'intronic' && c.intron) return <span className="text-slate-500" title={`${c.tx}\n${THEORETICAL}`}>intronic (between exons {c.intron[0]} and {c.intron[1]})</span>;
  if (!c || !c.text) return <>{g.exons ? 'cut (exons in)' : 'cut (no exon)'}</>;
  const tip = [`${c.tx} (${c.strand > 0 ? '+' : '−'} strand, ${c.n} exons)`,
    c.coding == null ? 'non-coding transcript' : `coding bases in the exons wholly inside: ${c.coding} (${c.coding % 3 === 0 ? 'a multiple of 3' : `${c.coding % 3} past a multiple of 3`})`,
    THEORETICAL].join('\n');
  return <span className={c.severe ? 'text-red-700' : 'text-slate-500'} title={tip}>{c.text}</span>;
}
