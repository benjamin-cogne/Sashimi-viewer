/**
 * Array-CGH in the SV browser: the button that adds samples (interval tables from CytoGenomics, or BED files), the
 * question of the minimum number of probes, and the frame listing each sample's calls between the search and the
 * browser. A call clicked is drawn in the browser as a deletion or a duplication, with the browser's own tracks.
 */
import { useEffect, useRef, useState } from 'react';
import { Icon, Stepper, Segmented, Pill } from './controls';
import { cgh, useCgh, passes, classOf, noteOf, isAnalysed, CGH_MIN_PROBES_DEFAULT, CGH_MAX_ANNOTATED, type CghState } from './cghStore';
import { CGH_CLASSES, classColor, cnvHubUrl, copyText, exportCghFile, type CghCnv, type CghFile, type CghSex } from './cghIntervals';
import { annotate, viewOf, bandOfCall, PANEL_SOURCE, PANELAPP_MAX_GENES, DGV_COVER, type DgvRelation } from './cghAnnotate';
import { LEVEL_RANK, LEVEL_STYLE, SCORE_RULES, type GeneView, type Level } from './cghScore';
import { MOI, type SvGeneLoaders } from './SvGenomeBrowser';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { MoiClass } from '../../standalone/panelapp';

const fmt = (n: number) => Math.round(n).toLocaleString('en-US');
const sizeText = (kb: number) => (kb >= 1000 ? `${(kb / 1000).toFixed(2)} Mb` : kb >= 1 ? `${kb.toFixed(kb >= 100 ? 0 : 1)} kb` : `${Math.round(kb * 1000)} bp`);
const KIND = { loss: { color: '#dc2626', bg: '#fee2e2', word: 'del' }, gain: { color: '#2563eb', bg: '#dbeafe', word: 'dup' } } as const;

/** "＋ Add samples": the files picked (several at once), read in the page. */
export function CghAddButton({ compact }: { compact?: boolean }) {
  const input = useRef<HTMLInputElement>(null);
  return (
    <label className={`inline-flex items-center gap-1.5 h-[30px] px-2.5 rounded-[6px] border text-[12.5px] font-medium cursor-pointer whitespace-nowrap bg-white border-slate-200 text-slate-700 hover:border-accent-300 hover:bg-accent-50 ${compact ? 'h-[26px] px-2 text-xs rounded-lg' : ''}`}
      title={'Array-CGH calls of one or more patients: the interval table of each sample exported by CytoGenomics (.txt: Chromosome, Start, Stop, Cytoband, Gene Name, Size(kb), Type, #Probes, Mean Log Ratio, ISCN, Suppress, Classification), a BED (chrom, start, end, a name saying del / dup or x1 / x3), or zip files holding them (one patient per zip).\nRead in this page, never sent anywhere. Files can also be dropped on the page.'}>
      <Icon name="upload" size={compact ? 13 : 15} />{compact ? 'Add samples' : 'Add array-CGH samples'}
      <input ref={input} data-cgh-input type="file" multiple accept=".txt,.tsv,.bed,.xls,.zip,text/plain,application/zip" className="hidden"
        onChange={e => { const fs = e.target.files; if (fs?.length) void cgh.add(fs); e.target.value = ''; }} />
    </label>
  );
}

/** The question asked once the files are read: the minimum number of probes, with what it keeps of each sample. */
export function CghProbeDialog() {
  const s = useCgh();
  return s.pending ? <ProbeQuestion p={s.pending} init={s.filter.minProbes} initKb={s.filter.minKb} loaded={s.files} /> : null;
}
/** the patient's sex: changes the priority of X-linked recessive genes and the copies drawn for chrX / chrY */
export const SEX_MARK: Record<CghSex, string> = { male: '♂', female: '♀' };
function SexSelect({ value, onChange, title }: { value: CghSex | null; onChange: (v: CghSex | null) => void; title?: string }) {
  return <select value={value ?? ''} onChange={e => onChange((e.target.value || null) as CghSex | null)} aria-label="Sex of the patient"
    title={title ?? "The patient's sex: a male's chrX loss is hemizygous (a green X-linked recessive gene is then HIGH), and his chrX / chrY are drawn with one copy"}
    className={`h-[24px] px-1 rounded border border-slate-200 bg-white text-[11.5px] outline-none cursor-pointer ${value ? 'text-slate-800 font-medium' : 'text-slate-400'}`}>
    <option value="">sex not given</option><option value="male">♂ male</option><option value="female">♀ female</option>
  </select>;
}
const SIZES = [0, 10, 50, 100, 200, 500, 1000];
const sizeWord = (v: number) => (v === 0 ? 'any' : v >= 1000 ? `${v / 1000} Mb` : `${v} kb`);
function ProbeQuestion({ p, init, initKb, loaded }: { p: NonNullable<CghState['pending']>; init: number; initKb: number; loaded: CghFile[] }) {
  const [n, setN] = useState<number>(init);
  const [kb, setKb] = useState<number>(initKb);
  // every call of the tab (loaded and added): too many, and a minimum must be raised before they are annotated
  const every = [...loaded, ...p.files].flatMap(f => f.samples.flatMap(x => x.cnvs));
  const kept = every.filter(c => isAnalysed(c, n, kb)).length;
  const many = every.length > CGH_MAX_ANNOTATED, over = kept > CGH_MAX_ANNOTATED;
  const [sexes, setSexes] = useState<Record<string, CghSex | null>>(() => Object.fromEntries(p.files.flatMap(f => f.samples.map(x => [x.key, x.sex]))));
  const samples = p.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const noCount = p.files.some(f => f.format === 'bed' || f.warnings.some(w => /no #Probes/.test(w)));
  return (
    <div className="fixed inset-0 z-[80] bg-black/40 flex items-start justify-center p-6" role="dialog" aria-modal="true" aria-label="Minimum number of probes"
      onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); cgh.cancel(); } if (e.key === 'Enter' && samples.length && !over && (e.target as HTMLElement).tagName !== 'SELECT') cgh.accept(n, sexes, kb); }}>
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
            <span className="font-semibold ml-3">size ≥</span>
            <select value={kb} onChange={e => setKb(Number(e.target.value))} aria-label="Minimum size" className="h-[30px] px-1.5 rounded-[6px] border border-slate-200 bg-white font-semibold outline-none cursor-pointer">
              {SIZES.map(v => <option key={v} value={v}>{sizeWord(v)}</option>)}
            </select>
          </div>}
          {many && <div className={`rounded-lg border px-3 py-2 text-[11.5px] ${over ? 'border-red-300 bg-red-50 text-red-800' : 'border-amber-300 bg-amber-50 text-amber-900'}`} role="alert">
            <b>{every.length.toLocaleString('en-US')} calls in all</b> (this tab's patients and these). Only the calls with {n}+ probes{kb ? ` and of ${sizeWord(kb)} or more` : ''} are annotated (genes, PanelApp, ClinGen, gnomAD, DGV Gold); the others are not.
            {over ? <> <b>{kept.toLocaleString('en-US')} would be: more than the {CGH_MAX_ANNOTATED} the page annotates.</b> Raise the minimum number of probes and / or the minimum size.</> : <> {kept.toLocaleString('en-US')} will be.</>}
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
          {noCount && <div className="text-[11px] text-amber-700">A file gives no probe count (BED): its calls are all kept.</div>}
          {p.files.flatMap(f => f.warnings.filter(w => !/^BED: no probe/.test(w)).map(w => <div key={f.id + w} className="text-[11px] text-amber-700">{f.name}: {w}</div>))}
          {p.notes.map(e => <div key={e} className="text-[11px] text-slate-500">{e}</div>)}
          {p.errors.map(e => <div key={e} className="text-[11px] text-red-700">{e}</div>)}
        </div>
        <div className="px-4 py-2.5 border-t border-slate-200 flex justify-end gap-2">
          <button onClick={() => cgh.cancel()} className="h-[30px] px-3 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50">Cancel</button>
          {samples.length > 0 && <button autoFocus disabled={over} onClick={() => cgh.accept(n, sexes, kb)} title={over ? `More than ${CGH_MAX_ANNOTATED} calls: raise a minimum first` : undefined} className="h-[30px] px-4 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700 disabled:opacity-40 disabled:cursor-not-allowed">Analyse {samples.length} sample{samples.length > 1 ? 's' : ''}</button>}
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
  // every analysed call of every patient: its genes (RefSeq, GRCh37), dosage, constraint, PanelApp, AChro-Puce
  useEffect(() => {
    if (!loadersFor) return;
    annotate(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.filter(c => isAnalysed(c, s.filter.minProbes, s.filter.minKb)).map(c => ({ c, build: fl.build ?? 'GRCh37' })))).slice(0, CGH_MAX_ANNOTATED), loadersFor);
  }, [s.files, s.filter.minProbes, s.filter.minKb, loadersFor]);
  // the call clicked kept in the table's view (the table only, never the page), when it changes
  const tableBox = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const box = tableBox.current, row = box?.querySelector<HTMLElement>('tr[data-shown]');
    if (!box || !row) return;
    const r = row.getBoundingClientRect(), bb = box.getBoundingClientRect(), head = 26;
    if (r.top < bb.top + head) box.scrollTop -= bb.top + head - r.top;
    else if (r.bottom > bb.bottom) box.scrollTop += r.bottom - bb.bottom;
  }, [s.shown]);
  if (!s.files.length) return null;
  const samples = s.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const listed = s.sample ? samples.filter(e => e.x.key === s.sample) : samples;
  const all = listed.flatMap(({ x, f }) => x.cnvs.map(c => ({ c, f })));
  const views = new Map(all.map(({ c }) => [c.id, viewOf(s, c)]));
  const rank = (id: string) => { const v = views.get(id)?.score; return v ? LEVEL_RANK[v.level] : -1; };
  const rows = all.filter(({ c }) => passes(s, c));
  if (s.filter.sort === 'priority') rows.sort((a, b) => rank(b.c.id) - rank(a.c.id));
  const fewProbes = all.filter(({ c }) => !isAnalysed(c, s.filter.minProbes, s.filter.minKb)).length;
  const analysedAll = s.files.reduce((k, fl) => k + fl.samples.reduce((m, x) => m + x.cnvs.filter(c => isAnalysed(c, s.filter.minProbes, s.filter.minKb)).length, 0), 0);
  const anySuppressed = all.some(({ c }) => c.suppressed);
  const f = s.filter;
  /** a call clicked: drawn below with its detail under its row; clicked again, its detail folds or unfolds */
  const show = (i: number, click = false) => {
    const r = rows[i]; if (!r) return;
    if (click && r.c.id === s.shown) { setFolded(f => !f); return; }
    cgh.setShown(r.c.id); onShow(r.c, r.f);
  };
  const at = rows.findIndex(r => r.c.id === s.shown);
  const exportAll = () => {
    for (const fl of new Set(listed.map(e => e.f))) {
      const text = exportCghFile(fl, c => classOf(s, c), c => noteOf(s, c));
      const a = document.createElement('a');
      a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain;charset=utf-8' }));
      a.download = fl.saveAs;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    }
    cgh.exported(listed.map(e => e.x.key));
  };
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
        <span className="ml-auto self-center flex items-center gap-2 py-1">
          <CghAddButton compact />
          <span className="text-[11px] text-slate-500" title="The build the interval table states (Genome hg19); the calls are drawn on it">{filesBuilds.join(', ')}{listed.some(e => !e.f.build) ? ' (no build in a file: GRCh37 assumed)' : ''}</span>
          <button onClick={exportAll} className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium ${listedDirty ? 'bg-amber-50 border-amber-300 text-amber-800' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300'}`}
            title={'Each file of the tab written back as it came, with the classes given here in its Classification column and the notes in a Comment column (added when the file has none; a BED gets both as two last columns).\nClasses and notes are kept in this tab only: export them to keep them.'}>
            <Icon name="download" size={13} />Export classified table{listedDirty ? ' •' : ''}</button>
        </span>
      </div>
      <div className="pt-2" />
      <div className="px-5 pb-2 flex flex-wrap items-center gap-2 text-xs">
        {listed.length === 1 && <SexSelect value={s.sex[listed[0].x.key] ?? null} onChange={v => cgh.setSex(listed[0].x.key, v)} />}
        <Stepper value={f.minProbes} onChange={v => cgh.setFilter({ minProbes: v })} min={1} max={500} label="probes ≥" title="Calls with fewer probes are not analysed (3 by default)" ariaLabel="Minimum number of probes" />
        <Segmented size="sm" label="Type" title="Losses, gains or both" value={f.kind} onChange={v => cgh.setFilter({ kind: v })}
          options={[{ value: 'all', label: 'all', hint: 'Losses and gains' }, { value: 'loss', label: 'losses', hint: 'Deletions only' }, { value: 'gain', label: 'gains', hint: 'Gains and amplifications only' }]} />
        <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title="Smaller calls are not analysed: neither listed nor annotated">size ≥
          <select value={f.minKb} onChange={e => cgh.setFilter({ minKb: Number(e.target.value) })} className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            {[0, 10, 50, 100, 200, 500, 1000].map(v => <option key={v} value={v}>{v === 0 ? 'any' : v >= 1000 ? `${v / 1000} Mb` : `${v} kb`}</option>)}
          </select></label>
        <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title="Classified or not (a class from the file counts)">class
          <select value={f.classed} onChange={e => cgh.setFilter({ classed: e.target.value as typeof f.classed })} className="bg-transparent font-medium text-slate-800 outline-none cursor-pointer">
            <option value="all">any</option><option value="unclassified">not classified</option><option value="classified">classified</option>
          </select></label>
        {anySuppressed && <Pill size="sm" label="hide suppressed" on={f.hideSuppressed} onChange={v => cgh.setFilter({ hideSuppressed: v })} title="Calls marked Suppress = true in CytoGenomics" />}
        <input value={f.query} onChange={e => cgh.setFilter({ query: e.target.value })} placeholder="gene, band (16p11), chromosome" aria-label="Filter the calls by gene, band or chromosome"
          className="h-[26px] w-[210px] px-2 rounded-lg border border-slate-200 text-xs placeholder:text-slate-400 focus:outline-none focus:border-accent-300" />
        <Segmented size="sm" label="Order" title="The calls in genomic order, or the highest priority first" value={f.sort} onChange={v => cgh.setFilter({ sort: v })}
          options={[{ value: 'position', label: 'position', hint: 'Genomic order, as in the file' }, { value: 'priority', label: 'priority', hint: 'HIGH first, then MODERATE, LOW, none' }]} />
        <span className="text-slate-500 ml-1">{rows.length} of {all.length} call{all.length > 1 ? 's' : ''} listed{fewProbes ? ` · ${fewProbes} under ${f.minProbes} probe${f.minProbes > 1 ? 's' : ''}${f.minKb ? ` or ${sizeWord(f.minKb)}` : ''} not analysed` : ''}</span>
      </div>
      {analysedAll > CGH_MAX_ANNOTATED && <div className="mx-5 mb-2 rounded-lg border border-red-300 bg-red-50 text-red-800 px-3 py-1.5 text-[11.5px]" role="alert">
        {analysedAll.toLocaleString('en-US')} calls analysed: only the first {CGH_MAX_ANNOTATED} are annotated. Raise the minimum number of probes and / or the minimum size.</div>}
      <div ref={tableBox} style={{ maxHeight: 'max(230px, 36vh)' }} className="mx-5 mb-1.5 border border-slate-200 rounded-lg overflow-auto outline-none focus:ring-2 focus:ring-accent-100" tabIndex={0} aria-label="Array-CGH calls (↑ ↓ to go from one to the next)"
        onKeyDown={e => { if ((e.key === 'ArrowDown' || e.key === 'ArrowUp') && !/^(SELECT|INPUT|TEXTAREA)$/.test((e.target as HTMLElement).tagName)) { e.preventDefault(); e.stopPropagation(); show(Math.max(0, Math.min(rows.length - 1, at < 0 ? 0 : at + (e.key === 'ArrowDown' ? 1 : -1)))); } }}>
        <table className="w-full text-[11.5px] border-collapse">
          <thead className="sticky top-0 z-10 bg-slate-50 text-slate-500 text-left whitespace-nowrap">
            <tr className="border-b border-slate-200">
              {!s.sample && samples.length > 1 && <th className="px-1.5 py-0.5 font-medium">Sample</th>}
              <th className="px-1.5 py-0.5 font-medium">Band</th><th className="px-1.5 py-0.5 font-medium">Position (GRCh37)</th>
              <th className="px-1.5 py-0.5 font-medium text-right">Size</th><th className="px-1.5 py-0.5 font-medium">Type</th>
              <th className="px-1.5 py-0.5 font-medium text-right" title="#Probes · mean log ratio">Probes · LogR</th>
              <th className="px-1.5 py-0.5 font-medium" title="RefSeq NM_ genes of the interval (UCSC, the file's build); the file's Gene Name only when they cannot be read">Genes</th>
              <th className="px-1.5 py-0.5 font-medium" title={`${PANEL_SOURCE.label}, every panel: the genes by their best rating, green first (with their mode of inheritance), then amber, then red`}>PanelApp</th>
              <th className="px-1.5 py-0.5 font-medium" title="ClinGen dosage scores (haploinsufficiency for a loss, triplosensitivity for a gain; 3 sufficient, 2 some evidence) and the lowest gnomAD LOEUF">Dosage</th>
              <th className="px-1.5 py-0.5 font-medium" title="A priority for the review from the genes' scores (the rules are listed under the table); not a classification">Priority</th>
              <th className="px-1.5 py-0.5 font-medium sticky right-[76px] z-[1] bg-slate-50 shadow-[-6px_0_6px_-6px_rgba(15,23,42,0.18)]">Classification</th><th className="px-1.5 py-0.5 font-medium sticky right-0 z-[1] bg-slate-50 w-[76px] min-w-[76px]">Links</th>
            </tr>
          </thead>
          <tbody>
            {rows.map(({ c, f: fl }, i) => {
              const cls = classOf(s, c), col = classColor(cls), on = s.shown === c.id, k = KIND[c.kind];
              const inList = CGH_CLASSES.some(o => o.value === cls), note = noteOf(s, c);
              const hub = cnvHubUrl(c, fl.build ?? 'GRCh37');
              const v = views.get(c.id)!;
              return [
                <tr key={c.id} data-cgh-row data-shown={on || undefined} onClick={() => show(i, true)} title={`${c.iscn ?? ''}${c.suppressed ? '\nSuppressed in CytoGenomics' : ''}\nClick: draw it in the browser below`}
                  className={`border-b border-slate-100 cursor-pointer ${on ? 'bg-accent-50 [&>td]:border-t-2 [&>td]:border-t-accent-600 [&>td:first-child]:border-l-2 [&>td:first-child]:border-l-accent-600 [&>td:last-child]:border-r-2 [&>td:last-child]:border-r-accent-600' : 'bg-white hover:bg-slate-50'} ${c.suppressed ? 'opacity-60' : ''}`}>
                  {!s.sample && samples.length > 1 && <td className="px-1.5 py-0.5 text-slate-600 whitespace-nowrap">{c.sample}{s.sex[c.sampleKey] ? ` ${SEX_MARK[s.sex[c.sampleKey]]}` : ''}</td>}
                  <td className="px-1.5 py-0.5 font-semibold whitespace-nowrap" title={!c.cytoband && v.ann?.band ? 'Band read from the cytogenetic bands (the file gives none)' : undefined}>{c.chrom.replace(/^chr/, '')}{bandOfCall(c, v.ann).replace(/\s+/g, '')}</td>
                  <td className="px-1.5 py-0.5 font-mono text-[10.5px] text-slate-700 whitespace-nowrap">{c.chrom}:{fmt(c.start)}-{fmt(c.end)}</td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap">{sizeText(c.sizeKb)}</td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap"><span className="inline-flex items-center rounded px-1.5 text-[10.5px] font-bold whitespace-nowrap" style={{ color: k.color, background: k.bg }} title={c.typeText}>{k.word} {copyText(c)}</span><DgvMedal v={v} /></td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap"><span className={c.probes == null ? 'text-slate-400' : ''}>{c.probes ?? '?'}</span>{c.logRatio != null && <span className="text-slate-500"> · {c.logRatio.toFixed(2)}</span>}</td>
                  <GenesCell c={c} v={v} />
                  <PanelCell v={v} />
                  <DosageCell c={c} v={v} />
                  <td className="px-1.5 py-0.5 whitespace-nowrap"><PriorityBadge v={v} /></td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap sticky right-[76px] z-[1] bg-inherit shadow-[-6px_0_6px_-6px_rgba(15,23,42,0.18)]" onClick={e => e.stopPropagation()}>
                    <select value={inList ? cls : ''} aria-label={`Classification of ${c.chrom}:${c.start}-${c.end}`} style={{ borderLeftColor: inList ? col ?? '#d3dbe5' : '#d3dbe5' }}
                      title={!inList && cls ? `The file says “${cls}”: kept in the export unless a class is chosen` : 'The class written in the export'}
                      onChange={e => cgh.setClass(c, e.target.value)}
                      className={`h-[20px] w-[128px] pl-1 pr-0.5 rounded border border-slate-200 border-l-[4px] bg-white text-[11px] outline-none cursor-pointer ${inList ? 'font-medium text-slate-800' : 'text-slate-400'}`}>
                      <option value="">— not classified</option>
                      {CGH_CLASSES.map(o => <option key={o.value} value={o.value}>{o.hint && o.value === 'NS' ? `NS (${o.hint})` : o.value}</option>)}
                    </select>
                    <button onClick={() => { setNoteFocus(c.id); setFolded(false); if (!on) show(i); }} className={`ml-1 px-1 rounded ${note ? 'text-accent-700 bg-accent-50 font-bold' : 'text-slate-300 hover:text-slate-600'}`}
                      title={note ? `Note: ${note}` : 'Add a note on this call (exported in the Comment column)'} aria-label="Note">✎</button>
                  </td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap sticky right-0 z-[1] bg-inherit w-[76px] min-w-[76px]" onClick={e => e.stopPropagation()}>
                    {hub ? <a href={hub} target="_blank" rel="noopener noreferrer" className="text-accent-700 hover:underline font-medium" title={`CNV-Hub (cnvhub.net, Dijon University Hospital): this call's page, ${hub}. Opens in a new tab; the call's position and type are in the address.`}>CNV-Hub ↗</a>
                      : <span className="text-slate-300" title="CNV-Hub pages are known for hg19 (GRCh37) only">CNV-Hub</span>}
                  </td>
                </tr>,
                on && !folded && <tr key={`${c.id}-detail`}><td colSpan={12} className="p-0 border-x-2 border-b-2 border-accent-600"><CnvDetail c={c} v={v} note={note} focusNote={noteFocus === c.id} onNoteFocused={() => setNoteFocus(null)} /></td></tr>,
              ];
            })}
            {!rows.length && <tr><td colSpan={12} className="px-3 py-3 text-slate-400">No call passes the filter.</td></tr>}
          </tbody>
        </table>
      </div>
      {s.files.some(fl => fl.warnings.some(w => !/^BED: no probe/.test(w))) && <div className="px-5 text-[10.5px] text-amber-700">{s.files.flatMap(fl => fl.warnings.filter(w => !/^BED: no probe/.test(w)).map(w => `${fl.name}: ${w}`)).join(' · ')}</div>}
      <details className="px-5 pb-2 text-[10.5px] text-slate-500 leading-4">
        <summary className="cursor-pointer text-accent-700 font-medium select-none">ⓘ How the calls are read, annotated, ranked and exported</summary>
        Click a call to draw it below as a deletion or a duplication (↑ ↓ for the next one), with its genes listed under its row (click it again to fold them): the genes, DGV Gold, ClinVar, AChro-Puce and PanelApp tracks of the browser are compared with it.
        Genes are RefSeq's NM_ models of the interval (UCSC ncbiRefSeqCurated, the file's build), the file's Gene Name only when they cannot be read (then marked <span className="text-amber-700">file</span>); ClinGen dosage through the UCSC Genome Browser, gnomAD v2.1.1 constraint on GRCh37, PanelApp from {PANEL_SOURCE.label} (every panel); DGV Gold marked when a record of the same type shares {Math.round(DGV_COVER * 100)} % of the call or more (the call inside it: 100 %), or lies inside the call whatever its share; one icon per way they lie, with the best share and how many (×n): <DgvIcon rel="within" /> the call inside a DGV CNV, <DgvIcon rel="contains" /> the call holding one, <DgvIcon rel="overlaps" /> overlapping on one side.
        A gain is drawn as a tandem duplication: array-CGH gives the copy number, not where the extra copy lies; a loss of both copies (x0) as homozygous.
        Positions and ISCN as in the file. Classes and notes are kept in this tab only, never in the browser's storage: <b>Export classified table</b> writes each file back with them (Classification column, and a Comment column).
        <PriorityRules />
      </details>
    </div>
  );
}

// ---- the columns of the genes' evidence ----
const MoiBadge = ({ m }: { m: MoiClass }) => <span className="inline-flex items-center rounded-full px-1 text-[9px] font-bold border leading-[13px]" style={{ color: MOI[m].color, borderColor: `${MOI[m].color}88`, background: `${MOI[m].color}1a` }} title={MOI[m].label}>{MOI[m].short}</span>;
const waitingDot = (n: number) => (n ? <span className="text-slate-400" title={`PanelApp: ${n} gene${n > 1 ? 's' : ''} still asked`}> …</span> : null);
type View = ReturnType<typeof viewOf>;

function GenesCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-1.5 py-0.5 text-slate-400">…</td>;
  const gs = v.ann.genes, file = v.ann.source === 'file';
  return (
    <td className="px-1.5 py-0.5 max-w-[140px] truncate" title={`${file ? `From the file (RefSeq not read): ${c.genes.join(', ')}` : `${gs.length} RefSeq gene${gs.length > 1 ? 's' : ''}: ${gs.map(g => `${g.symbol}${g.whole ? '' : ' (cut)'}`).join(', ')}`}`}>
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
    <td className="px-1.5 py-0.5 max-w-[180px] overflow-hidden" title={title || undefined}>
      <span className="flex items-center gap-x-1.5 whitespace-nowrap">
        {green.slice(0, 2).map(x => <span key={x.g.symbol} className="inline-flex items-center gap-0.5"><i>{x.g.symbol}</i>{modesAt(x).map(m => <MoiBadge key={m} m={m} />)}</span>)}
        {green.length > 2 && <span className="text-[10px] font-semibold" style={{ color: RATING[3].color }}>+{green.length - 2} green</span>}
        {amber.length > 0 && <span className="text-[10px] font-semibold" style={{ color: RATING[2].color }}>{green.length ? '' : <i className="font-normal">{amber[0].g.symbol} </i>}{amber.length} amber</span>}
        {red.length > 0 && <span className="text-[10px] font-semibold" style={{ color: RATING[1].color }}>{red.length} red</span>}
        {!green.length && !amber.length && !red.length && !v.waiting && <span className="text-slate-300">—</span>}
        {waitingDot(v.waiting)}
      </span>
    </td>
  );
}
function DosageCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-1.5 py-0.5 text-slate-400">…</td>;
  const key = c.kind === 'loss' ? 'hi' : 'ts', word = c.kind === 'loss' ? 'HI' : 'TS';
  const scored = v.genes.filter(x => x[key] === 3 || x[key] === 2).sort((a, b) => (b[key] ?? 0) - (a[key] ?? 0));
  const low = v.genes.filter(x => x.loeuf != null).sort((a, b) => a.loeuf! - b.loeuf!)[0];
  return (
    <td className="px-1.5 py-0.5 whitespace-nowrap">
      {scored.slice(0, 1).map(x => <span key={x.g.symbol} className="mr-1.5 font-semibold" style={{ color: x[key] === 3 ? '#b91c1c' : '#c2410c' }} title={`ClinGen ${c.kind === 'loss' ? 'haploinsufficiency' : 'triplosensitivity'} score ${x[key]} (${x[key] === 3 ? 'sufficient' : 'some'} evidence)`}>{word}{x[key]} <i className="font-normal">{x.g.symbol}</i></span>)}
      {scored.length > 1 && <span className="text-slate-500 mr-1.5" title={scored.slice(1).map(x => `${word}${x[key]} ${x.g.symbol}`).join(', ')}>+{scored.length - 1}</span>}
      {low && <span className="text-slate-600" title={`Lowest gnomAD LOEUF of the call (${low.version ?? ''}): ${low.g.symbol} ${low.loeuf!.toFixed(2)}${low.pli != null ? `, pLI ${low.pli.toFixed(2)}` : ''}`}>LOEUF {low.loeuf!.toFixed(2)}</span>}
      {!scored.length && !low && <span className="text-slate-300">—</span>}
    </td>
  );
}
function PriorityBadge({ v }: { v: View }) {
  if (!v.score) return <span className="text-slate-400">…</span>;
  const st = LEVEL_STYLE[v.score.level], incomplete = v.missing.length > 0 || v.waiting > 0;
  const title = [...v.score.reasons.map(r => `${LEVEL_STYLE[r.level].label}: ${r.text}`), ...(v.waiting ? [`PanelApp still asked for ${v.waiting} gene${v.waiting > 1 ? 's' : ''}`] : []), ...v.missing.map(m => `Not available: ${m}`)].join('\n') || 'No protein-coding gene in the call';
  return <span className={`inline-flex items-center rounded px-1.5 text-[10.5px] font-bold ${incomplete && v.score.level !== 'high' ? 'border border-dashed' : ''}`} style={{ color: st.color, background: st.bg, borderColor: st.color }} title={title}>{st.label}{incomplete && v.score.level !== 'high' ? ' ?' : ''}</span>;
}

/** under the call's row: why its priority, the AChro-Puce regions, and each gene with its scores */
function CnvDetail({ c, v, note, focusNote, onNoteFocused }: { c: CghCnv; v: View; note: string; focusNote: boolean; onNoteFocused: () => void }) {
  const [allGenes, setAllGenes] = useState(false);
  const noteBox = <textarea value={note} onChange={e => cgh.setNote(c, e.target.value)} rows={1} placeholder="Note on this call (exported in the Comment column)"
    ref={el => { if (el && focusNote) { el.focus(); onNoteFocused(); } }}
    onKeyDown={e => e.stopPropagation()} aria-label="Note on this call"
    className="w-full mb-1 px-2 py-1 rounded border border-slate-300 bg-white text-[11.5px] text-slate-800 resize-y outline-none focus:border-accent-400 focus:ring-2 focus:ring-accent-100" />;
  if (!v.ann || v.ann.state === 'loading' || !v.score) return <div className="px-4 py-2 bg-slate-50 text-slate-500 text-[11px]" onClick={e => e.stopPropagation()}>{noteBox}Reading the genes of the call…</div>;
  const key = c.kind === 'loss' ? 'hi' : 'ts';
  // PanelApp green first, then amber, then red, then the others; within each, the dosage score, then the position
  const genes = v.genes.slice().sort((a, b) => (b.pa?.best ?? 0) - (a.pa?.best ?? 0) || (b[key] ?? 0) % 10 - (a[key] ?? 0) % 10 || a.g.start - b.g.start);
  const shownGenes = allGenes ? genes : genes.slice(0, 8);
  const dot = (l: Level) => <span className="inline-block w-2 h-2 rounded-full mr-1.5 align-middle" style={{ background: LEVEL_STYLE[l].color }} />;
  return (
    <div className="px-4 py-2 bg-slate-50 text-[11px] text-slate-700 cursor-default" onClick={e => e.stopPropagation()}>
      {noteBox}
      <div className="flex flex-wrap gap-x-6 gap-y-1">
        <div className="min-w-[320px] flex-1">
          <div className="font-semibold text-slate-900">Priority: <span style={{ color: LEVEL_STYLE[v.score.level].color }}>{LEVEL_STYLE[v.score.level].label}</span> <span className="font-normal text-slate-500">(a proposal for the review, not a class)</span></div>
          {v.score.reasons.length ? v.score.reasons.map(r => <div key={r.id}>{dot(r.level)}{r.text}</div>) : <div className="text-slate-500">No protein-coding gene in the call.</div>}
          {v.waiting > 0 && <div className="text-slate-500">PanelApp: {v.waiting} gene{v.waiting > 1 ? 's' : ''} still asked…</div>}
          {v.missing.map(m => <div key={m} className="text-amber-700">Not available: {m}</div>)}
        </div>
        {(v.ann.dgv.length > 0 || v.ann.dgvError) && <div className="min-w-[240px] flex-1">
          <div className="font-semibold text-slate-900">DGV Gold (same type: ≥ {Math.round(DGV_COVER * 100)} % of the call shared, or inside the call whatever its share)</div>
          {dgvShown(v.ann.dgv).map(x => <div key={x.id} className="flex flex-wrap items-center gap-x-1"><DgvIcon rel={x.rel} />{x.id} <span className="text-slate-500">({x.type})</span>: <b>{pctText(x.cover)}</b> <span className="text-slate-500">· {DGV_REL[x.rel].word}</span>{x.freq ? ` · frequency ${x.freq}` : ''}{x.samples != null ? ` · ${x.samples} samples` : ''}</div>)}
          {v.ann.dgv.length > dgvShown(v.ann.dgv).length && <div className="text-slate-500">+{v.ann.dgv.length - dgvShown(v.ann.dgv).length} more DGV CNVs inside the call (in the medallion's tooltip)</div>}
          {v.ann.dgvError && <div className="text-amber-700">not available ({v.ann.dgvError})</div>}
        </div>}
        {v.ann.achro.length > 0 && <div className="min-w-[240px]">
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
            <td className="px-1.5 py-0.5">{v.ann!.source === 'file' ? <span className="text-amber-700" title="RefSeq could not be read: the file names the gene, its extent is not known here">from the file</span> : x.g.whole ? 'whole' : x.g.exons ? 'cut (exons in)' : 'cut (no exon)'}</td>
            <td className="px-1.5 py-0.5 font-semibold" style={{ color: x.hi === 3 ? '#b91c1c' : x.hi === 2 ? '#c2410c' : undefined }}>{x.hi ?? ''}</td>
            <td className="px-1.5 py-0.5 font-semibold" style={{ color: x.ts === 3 ? '#b91c1c' : x.ts === 2 ? '#c2410c' : undefined }}>{x.ts ?? ''}</td>
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
/** a share: under 1 % said so rather than 0 % */
const pctText = (f: number) => (f > 0 && f < 0.005 ? '<1%' : `${Math.round(f * 100)}%`);
/** how a DGV Gold CNV lies against the call, drawn: the call (dark) and the DGV CNV (violet) */
export const DGV_REL: Record<DgvRelation, { word: string; call: [number, number]; dgv: [number, number] }> = {
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
  return <span className="ml-1 whitespace-nowrap inline-flex items-center gap-1 rounded-full px-1.5 text-[9.5px] font-bold border border-violet-300 bg-violet-50 text-violet-700 leading-[15px] align-middle" title={`DGV Gold, same type:\n${title}`}>
    {rels.map(({ rel, hs }) => <span key={rel} className="inline-flex items-center gap-0.5"><DgvIcon rel={rel} />{pctText(Math.max(...hs.map(x => x.cover)))}</span>)}
  </span>;
}
