/**
 * Array-CGH in the SV browser: the button that adds samples (interval tables from CytoGenomics, or BED files), the
 * question of the minimum number of probes, and the frame listing each sample's calls between the search and the
 * browser. A call clicked is drawn in the browser as a deletion or a duplication, with the browser's own tracks.
 */
import { useEffect, useRef, useState } from 'react';
import { Icon, Stepper, Segmented, Pill } from './controls';
import { cgh, useCgh, passes, classOf, CGH_MIN_PROBES_DEFAULT, type CghState } from './cghStore';
import { CGH_CLASSES, classColor, cnvHubUrl, copyText, exportCghFile, type CghCnv, type CghFile } from './cghIntervals';
import { annotate, viewOf, PANEL_SOURCE, PANELAPP_MAX_GENES } from './cghAnnotate';
import { LEVEL_RANK, LEVEL_STYLE, SCORE_RULES, type GeneView, type Level } from './cghScore';
import { MOI, type SvGeneLoaders } from './SvGenomeBrowser';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { MoiClass } from '../../standalone/panelapp';

const fmt = (n: number) => Math.round(n).toLocaleString('en-US');
const sizeText = (kb: number) => (kb >= 1000 ? `${(kb / 1000).toFixed(2)} Mb` : kb >= 1 ? `${kb.toFixed(kb >= 100 ? 0 : 1)} kb` : `${Math.round(kb * 1000)} bp`);
const KIND = { loss: { color: '#dc2626', bg: '#fee2e2', word: 'del' }, gain: { color: '#2563eb', bg: '#dbeafe', word: 'dup' } } as const;
const OTHER = '\u0000other';

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
  return s.pending ? <ProbeQuestion p={s.pending} init={s.filter.minProbes} /> : null;
}
function ProbeQuestion({ p, init }: { p: NonNullable<CghState['pending']>; init: number }) {
  const [n, setN] = useState<number>(init);
  const samples = p.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const noCount = p.files.some(f => f.format === 'bed' || f.warnings.some(w => /no #Probes/.test(w)));
  return (
    <div className="fixed inset-0 z-[80] bg-black/40 flex items-start justify-center p-6" role="dialog" aria-modal="true" aria-label="Minimum number of probes"
      onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); cgh.cancel(); } if (e.key === 'Enter' && samples.length) cgh.accept(n); }}>
      <div className="bg-white rounded-xl shadow-2xl border border-slate-200 w-full max-w-lg text-slate-900 text-xs mt-16">
        <div className="px-4 pt-3 pb-2 border-b border-slate-200">
          <div className="font-bold text-sm">Array-CGH samples</div>
          <div className="text-[11px] text-slate-500 mt-0.5">Calls with fewer probes than this are not analysed (hidden from the list; the number can be changed in the list afterwards).</div>
        </div>
        <div className="px-4 py-3 space-y-3">
          {samples.length > 0 && <div className="flex items-center gap-2">
            <span className="font-semibold">Minimum number of probes</span>
            <Stepper value={n} onChange={setN} min={1} max={500} title="Calls with fewer probes are not analysed" ariaLabel="Minimum number of probes" />
            {n !== CGH_MIN_PROBES_DEFAULT && <button className="text-accent-700 hover:underline" onClick={() => setN(CGH_MIN_PROBES_DEFAULT)}>back to {CGH_MIN_PROBES_DEFAULT}</button>}
          </div>}
          {samples.length > 0 && <table className="w-full text-[11.5px]">
            <thead><tr className="text-left text-slate-500 border-b border-slate-200"><th className="py-1 font-medium">Sample</th><th className="font-medium">File</th><th className="font-medium">Build</th><th className="font-medium text-right">Calls</th><th className="font-medium text-right">Analysed</th></tr></thead>
            <tbody>{samples.map(({ x, f }) => {
              const kept = x.cnvs.filter(c => c.probes == null || c.probes >= n).length;
              return <tr key={x.key} className="border-b border-slate-100">
                <td className="py-1 font-semibold">{x.name}</td><td className="text-slate-500 truncate max-w-[150px]" title={f.name}>{f.name}</td>
                <td className={f.build ? '' : 'text-amber-700'} title={f.buildText ? `the file says Genome ${f.buildText}` : 'the file states no build: read as GRCh37'}>{f.build ?? 'GRCh37 ?'}</td>
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
          {samples.length > 0 && <button autoFocus onClick={() => cgh.accept(n)} className="h-[30px] px-4 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700">Analyse {samples.length} sample{samples.length > 1 ? 's' : ''}</button>}
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
  const [other, setOther] = useState<{ id: string; text: string } | null>(null);
  // every analysed call of every patient: its genes (RefSeq, GRCh37), dosage, constraint, PanelApp, AChro-Puce
  useEffect(() => {
    if (!loadersFor) return;
    annotate(s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.filter(c => c.probes == null || c.probes >= s.filter.minProbes).map(c => ({ c, build: fl.build ?? 'GRCh37' })))), loadersFor);
  }, [s.files, s.filter.minProbes, loadersFor]);
  if (!s.files.length) return null;
  const samples = s.files.flatMap(f => f.samples.map(x => ({ x, f })));
  const listed = s.sample ? samples.filter(e => e.x.key === s.sample) : samples;
  const all = listed.flatMap(({ x, f }) => x.cnvs.map(c => ({ c, f })));
  const views = new Map(all.map(({ c }) => [c.id, viewOf(s, c)]));
  const rank = (id: string) => { const v = views.get(id)?.score; return v ? LEVEL_RANK[v.level] : -1; };
  const rows = all.filter(({ c }) => passes(s, c));
  if (s.filter.sort === 'priority') rows.sort((a, b) => rank(b.c.id) - rank(a.c.id));
  const fewProbes = all.filter(({ c }) => c.probes != null && c.probes < s.filter.minProbes).length;
  const anySuppressed = all.some(({ c }) => c.suppressed);
  const f = s.filter;
  // the classes offered: the network's, then the words of the files and those typed here
  const own = [...new Set([...s.files.flatMap(fl => fl.samples.flatMap(x => x.cnvs.map(c => c.classification))), ...Object.values(s.classes)].map(v => v.trim()).filter(v => v && !CGH_CLASSES.some(c => c.value === v)))];
  const show = (i: number) => { const r = rows[i]; if (r) { cgh.setShown(r.c.id); onShow(r.c, r.f); } };
  const at = rows.findIndex(r => r.c.id === s.shown);
  const exportAll = () => {
    for (const fl of new Set(listed.map(e => e.f))) {
      const text = exportCghFile(fl, c => classOf(s, c));
      const a = document.createElement('a');
      a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain;charset=utf-8' }));
      a.download = fl.saveAs;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    }
    cgh.exported();
  };
  const tab = (on: boolean) => `inline-flex items-center gap-1.5 h-8 pl-3 pr-2 -mb-px rounded-t-lg border text-xs whitespace-nowrap ${on ? 'bg-white border-slate-200 border-b-white text-slate-900 font-semibold shadow-[inset_0_2px_0_#0f766e]' : 'bg-slate-50 border-transparent text-slate-600 hover:text-slate-900 hover:bg-slate-100'}`;
  const filesBuilds = [...new Set(listed.map(e => e.f.build ?? 'GRCh37'))];
  return (
    <div className={`mx-4 mt-3 bg-white border rounded-xl shadow-sm ${dropping ? 'border-accent-400 ring-2 ring-accent-100' : 'border-slate-200'}`} data-cgh-frame>
      <div className="px-5 pt-3 pb-2 flex flex-wrap items-center gap-2">
        <div className="text-base font-bold mr-1">Array-CGH</div>
        <CghAddButton compact />
        <span className="ml-auto flex items-center gap-2">
          <span className="text-[11px] text-slate-500" title="The build the interval table states (Genome hg19); the calls are drawn on it">{filesBuilds.join(', ')}{listed.some(e => !e.f.build) ? ' (no build in a file: GRCh37 assumed)' : ''}</span>
          <button onClick={exportAll} className={`inline-flex items-center gap-1.5 h-[26px] px-2 rounded-lg border text-xs font-medium ${s.dirty ? 'bg-amber-50 border-amber-300 text-amber-800' : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300'}`}
            title={'Each file of the tab written back as it came, its Classification column filled with the classes given here (a BED gets them as a last column).\nThe classes are kept in this tab only: export them to keep them.'}>
            <Icon name="download" size={13} />Export classified table{s.dirty ? ' •' : ''}</button>
        </span>
      </div>
      {/* one tab per patient, and all of them in one table */}
      <div className="px-5 flex flex-wrap items-end gap-1 border-b border-slate-200" role="tablist" aria-label="Patients">
        {samples.length > 1 && <button role="tab" aria-selected={!s.sample} className={tab(!s.sample)} onClick={() => cgh.setSample(null)}
          title="Every patient in one table, the sample as a column">All patients <span className="text-slate-400 font-normal tabular-nums">{samples.reduce((n, e) => n + e.x.cnvs.filter(c => passes(s, c)).length, 0)}</span></button>}
        {samples.map(({ x, f: fl }) => {
          const n = x.cnvs.filter(c => passes(s, c)).length, on = s.sample === x.key || samples.length === 1;
          return <span key={x.key} role="tab" aria-selected={on} className={tab(on)}>
            <button onClick={() => cgh.setSample(x.key)} title={`${fl.name}${fl.meta.length ? `\n${fl.meta.map(([k, v]) => `${k}: ${v}`).join('\n')}` : ''}`}>{x.name} <span className="text-slate-400 font-normal tabular-nums">{n}/{x.cnvs.length}</span></button>
            <button className="text-slate-400 hover:text-red-600 leading-none text-sm px-0.5" title="Take this patient off the list"
              onClick={() => { if (!s.dirty || !x.cnvs.some(c => s.classes[c.id] != null) || confirm(`${x.name}: the classes given here are not exported yet and will be lost.`)) cgh.remove(x.key); }}>×</button>
          </span>;
        })}
      </div>
      <div className="pt-2" />
      <div className="px-5 pb-2 flex flex-wrap items-center gap-2 text-xs">
        <Stepper value={f.minProbes} onChange={v => cgh.setFilter({ minProbes: v })} min={1} max={500} label="probes ≥" title="Calls with fewer probes are not analysed (3 by default)" ariaLabel="Minimum number of probes" />
        <Segmented size="sm" label="Type" title="Losses, gains or both" value={f.kind} onChange={v => cgh.setFilter({ kind: v })}
          options={[{ value: 'all', label: 'all', hint: 'Losses and gains' }, { value: 'loss', label: 'losses', hint: 'Deletions only' }, { value: 'gain', label: 'gains', hint: 'Gains and amplifications only' }]} />
        <label className="inline-flex items-center gap-1 h-[26px] pl-2 pr-1 rounded-lg border border-slate-200 bg-white text-slate-500" title="The smallest call listed">size ≥
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
        <span className="text-slate-500 ml-1">{rows.length} of {all.length} call{all.length > 1 ? 's' : ''} listed{fewProbes ? ` · ${fewProbes} under ${f.minProbes} probe${f.minProbes > 1 ? 's' : ''} not analysed` : ''}</span>
      </div>
      <div className="mx-5 mb-2 border border-slate-200 rounded-lg overflow-auto max-h-[460px] outline-none focus:ring-2 focus:ring-accent-100" tabIndex={0} aria-label="Array-CGH calls (↑ ↓ to go from one to the next)"
        onKeyDown={e => { if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); e.stopPropagation(); show(Math.max(0, Math.min(rows.length - 1, at < 0 ? 0 : at + (e.key === 'ArrowDown' ? 1 : -1)))); } }}>
        <table className="w-full text-[11.5px] border-collapse">
          <thead className="sticky top-0 z-10 bg-slate-50 text-slate-500 text-left">
            <tr className="border-b border-slate-200">
              {!s.sample && samples.length > 1 && <th className="px-2 py-1 font-medium">Sample</th>}
              <th className="px-2 py-1 font-medium">Band</th><th className="px-2 py-1 font-medium">Position (GRCh37)</th>
              <th className="px-2 py-1 font-medium text-right">Size</th><th className="px-2 py-1 font-medium">Type</th>
              <th className="px-2 py-1 font-medium text-right" title="#Probes">Probes</th><th className="px-2 py-1 font-medium text-right" title="Mean log ratio">Log ratio</th>
              <th className="px-2 py-1 font-medium" title="RefSeq NM_ genes of the interval (UCSC, the file's build); the file's Gene Name only when they cannot be read">Genes</th>
              <th className="px-2 py-1 font-medium" title={`${PANEL_SOURCE.label}, every panel: the genes rated green (diagnostic grade) with their mode of inheritance; amber genes counted`}>PanelApp green</th>
              <th className="px-2 py-1 font-medium" title="ClinGen dosage scores (haploinsufficiency for a loss, triplosensitivity for a gain; 3 sufficient, 2 some evidence) and the lowest gnomAD LOEUF">Dosage</th>
              <th className="px-2 py-1 font-medium" title="A priority for the review from the genes' scores (the rules are listed under the table); not a classification">Priority</th>
              <th className="px-2 py-1 font-medium">Classification</th><th className="px-2 py-1 font-medium">Links</th>
            </tr>
          </thead>
          <tbody>
            {rows.map(({ c, f: fl }, i) => {
              const cls = classOf(s, c), col = classColor(cls), on = s.shown === c.id, k = KIND[c.kind];
              const hub = cnvHubUrl(c, fl.build ?? 'GRCh37');
              const v = views.get(c.id)!;
              return [
                <tr key={c.id} data-cgh-row onClick={() => show(i)} title={`${c.iscn ?? ''}${c.suppressed ? '\nSuppressed in CytoGenomics' : ''}\nClick: draw it in the browser below`}
                  className={`border-b border-slate-100 cursor-pointer ${on ? 'bg-accent-50 shadow-[inset_3px_0_0_#0f766e]' : 'hover:bg-slate-50'} ${c.suppressed ? 'opacity-60' : ''}`}>
                  {!s.sample && samples.length > 1 && <td className="px-2 py-1 text-slate-600 whitespace-nowrap">{c.sample}</td>}
                  <td className="px-2 py-1 font-semibold whitespace-nowrap">{c.chrom.replace(/^chr/, '')}{c.cytoband.replace(/\s+/g, '')}</td>
                  <td className="px-2 py-1 font-mono text-[11px] text-slate-700 whitespace-nowrap">{c.chrom}:{fmt(c.start)}-{fmt(c.end)}</td>
                  <td className="px-2 py-1 text-right tabular-nums whitespace-nowrap">{sizeText(c.sizeKb)}</td>
                  <td className="px-2 py-1 whitespace-nowrap"><span className="inline-flex items-center rounded px-1.5 text-[10.5px] font-bold" style={{ color: k.color, background: k.bg }} title={c.typeText}>{k.word} {copyText(c)}</span></td>
                  <td className={`px-2 py-1 text-right tabular-nums ${c.probes == null ? 'text-slate-400' : ''}`}>{c.probes ?? '?'}</td>
                  <td className="px-2 py-1 text-right tabular-nums text-slate-600">{c.logRatio == null ? '' : c.logRatio.toFixed(2)}</td>
                  <GenesCell c={c} v={v} />
                  <PanelCell v={v} />
                  <DosageCell c={c} v={v} />
                  <td className="px-2 py-1 whitespace-nowrap"><PriorityBadge v={v} /></td>
                  <td className="px-2 py-1" onClick={e => e.stopPropagation()}>
                    {other?.id === c.id
                      ? <form className="flex gap-1" onSubmit={e => { e.preventDefault(); cgh.setClass(c, other.text); setOther(null); }}>
                          <input autoFocus value={other.text} onChange={e => setOther({ id: c.id, text: e.target.value })} onKeyDown={e => { if (e.key === 'Escape') { e.stopPropagation(); setOther(null); } }} placeholder="the class, as your lab writes it"
                            className="h-[22px] w-[170px] px-1.5 rounded border border-accent-300 text-[11px] outline-none" />
                          <button className="px-1.5 rounded bg-accent-600 text-white text-[11px]">OK</button></form>
                      : <select value={cls} aria-label={`Classification of ${c.chrom}:${c.start}-${c.end}`} style={{ borderLeftColor: col ?? '#d3dbe5' }}
                          onChange={e => { if (e.target.value === OTHER) setOther({ id: c.id, text: cls }); else cgh.setClass(c, e.target.value); }}
                          className={`h-[22px] max-w-[210px] pl-1.5 pr-1 rounded border border-slate-200 border-l-[4px] bg-white text-[11px] outline-none cursor-pointer ${cls ? 'font-medium text-slate-800' : 'text-slate-400'}`}>
                          <option value="">— not classified</option>
                          {CGH_CLASSES.map(o => <option key={o.value} value={o.value}>{o.value}</option>)}
                          {own.length > 0 && <option disabled>── from the files / typed ──</option>}
                          {own.map(v => <option key={v} value={v}>{v}</option>)}
                          <option value={OTHER}>Other (type it)…</option>
                        </select>}
                  </td>
                  <td className="px-2 py-1 whitespace-nowrap" onClick={e => e.stopPropagation()}>
                    {hub ? <a href={hub} target="_blank" rel="noopener noreferrer" className="text-accent-700 hover:underline font-medium" title={`CNV-Hub (cnvhub.net, Dijon University Hospital): this call's page, ${hub}. Opens in a new tab; the call's position and type are in the address.`}>CNV-Hub ↗</a>
                      : <span className="text-slate-300" title="CNV-Hub pages are known for hg19 (GRCh37) only">CNV-Hub</span>}
                  </td>
                </tr>,
                on && <tr key={`${c.id}-detail`}><td colSpan={13} className="p-0"><CnvDetail c={c} v={v} /></td></tr>,
              ];
            })}
            {!rows.length && <tr><td colSpan={13} className="px-3 py-3 text-slate-400">No call passes the filter.</td></tr>}
          </tbody>
        </table>
      </div>
      <div className="px-5 pb-3 text-[10.5px] text-slate-500 leading-4">
        Click a call to draw it below as a deletion or a duplication (↑ ↓ for the next one), with its genes listed under its row: the genes, DGV Gold, ClinVar, AChro-Puce and PanelApp tracks of the browser are compared with it.
        Genes are RefSeq's NM_ models of the interval (UCSC ncbiRefSeqCurated, the file's build), the file's Gene Name only when they cannot be read (then marked <span className="text-amber-700">file</span>); ClinGen dosage through the UCSC Genome Browser, gnomAD v2.1.1 constraint on GRCh37, PanelApp from {PANEL_SOURCE.label} (every panel).
        <PriorityRules />
        A gain is drawn as a tandem duplication: array-CGH gives the copy number, not where the extra copy lies; a loss of both copies (x0) as homozygous.
        Positions and ISCN as in the file (hover a row for its ISCN). The classes are kept in this tab only, never in the browser's storage: <b>Export classified table</b> writes each file back with them.
        {s.files.some(fl => fl.warnings.length) && <span className="text-amber-700"> {s.files.flatMap(fl => fl.warnings.filter(w => !/^BED: no probe/.test(w)).map(w => `${fl.name}: ${w}`)).join(' · ')}</span>}
      </div>
    </div>
  );
}

// ---- the columns of the genes' evidence ----
const MoiBadge = ({ m }: { m: MoiClass }) => <span className="inline-flex items-center rounded-full px-1 text-[9px] font-bold border leading-[13px]" style={{ color: MOI[m].color, borderColor: `${MOI[m].color}88`, background: `${MOI[m].color}1a` }} title={MOI[m].label}>{MOI[m].short}</span>;
const waitingDot = (n: number) => (n ? <span className="text-slate-400" title={`PanelApp: ${n} gene${n > 1 ? 's' : ''} still asked`}> …</span> : null);
type View = ReturnType<typeof viewOf>;

function GenesCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-2 py-1 text-slate-400">…</td>;
  const gs = v.ann.genes, file = v.ann.source === 'file';
  return (
    <td className="px-2 py-1 max-w-[230px] truncate" title={`${file ? `From the file (RefSeq not read): ${c.genes.join(', ')}` : `${gs.length} RefSeq gene${gs.length > 1 ? 's' : ''}: ${gs.map(g => `${g.symbol}${g.whole ? '' : ' (cut)'}`).join(', ')}`}`}>
      {file && <span className="text-amber-700 font-semibold mr-1">file</span>}
      {gs.length ? <><span className="tabular-nums text-slate-500 mr-1">{gs.length}</span><i>{gs.slice(0, 4).map(g => g.symbol).join(', ')}</i>{gs.length > 4 ? <span className="text-slate-400"> +{gs.length - 4}</span> : null}</>
        : <span className="text-slate-400">{file ? '—' : 'no gene'}</span>}
    </td>
  );
}
function PanelCell({ v }: { v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-2 py-1 text-slate-400">…</td>;
  const green = v.genes.filter(x => x.pa?.best === 3), amber = v.genes.filter(x => x.pa?.best === 2).length;
  return (
    <td className="px-2 py-1 max-w-[230px]" title={green.map(x => `${x.g.symbol}: ${x.pa!.green.map(m => MOI[m].label).join(', ')}${x.pa!.diseases.length ? ` · ${x.pa!.diseases.slice(0, 3).join('; ')}` : ''}`).join('\n') || undefined}>
      <span className="inline-flex flex-wrap items-center gap-x-1.5 gap-y-0.5">
        {green.slice(0, 4).map(x => <span key={x.g.symbol} className="inline-flex items-center gap-0.5 whitespace-nowrap"><i>{x.g.symbol}</i>{x.pa!.green.filter(m => m !== '?').map(m => <MoiBadge key={m} m={m} />)}</span>)}
        {green.length > 4 && <span className="text-slate-400">+{green.length - 4}</span>}
        {!green.length && !v.waiting && <span className="text-slate-300">—</span>}
        {amber > 0 && <span className="text-amber-700 text-[10.5px]" title="Genes rated amber in a panel at best">{amber} amber</span>}
        {waitingDot(v.waiting)}
      </span>
    </td>
  );
}
function DosageCell({ c, v }: { c: CghCnv; v: View }) {
  if (!v.ann || v.ann.state === 'loading') return <td className="px-2 py-1 text-slate-400">…</td>;
  const key = c.kind === 'loss' ? 'hi' : 'ts', word = c.kind === 'loss' ? 'HI' : 'TS';
  const scored = v.genes.filter(x => x[key] === 3 || x[key] === 2).sort((a, b) => (b[key] ?? 0) - (a[key] ?? 0));
  const low = v.genes.filter(x => x.loeuf != null).sort((a, b) => a.loeuf! - b.loeuf!)[0];
  return (
    <td className="px-2 py-1 whitespace-nowrap">
      {scored.slice(0, 2).map(x => <span key={x.g.symbol} className="mr-1.5 font-semibold" style={{ color: x[key] === 3 ? '#b91c1c' : '#c2410c' }} title={`ClinGen ${c.kind === 'loss' ? 'haploinsufficiency' : 'triplosensitivity'} score ${x[key]} (${x[key] === 3 ? 'sufficient' : 'some'} evidence)`}>{word}{x[key]} <i className="font-normal">{x.g.symbol}</i></span>)}
      {scored.length > 2 && <span className="text-slate-400 mr-1.5">+{scored.length - 2}</span>}
      {low && <span className="text-slate-600" title={`Lowest gnomAD LOEUF of the call (${low.version ?? ''}): ${low.g.symbol} ${low.loeuf!.toFixed(2)}${low.pli != null ? `, pLI ${low.pli.toFixed(2)}` : ''}`}>LOEUF {low.loeuf!.toFixed(2)} <i>{low.g.symbol}</i></span>}
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
function CnvDetail({ c, v }: { c: CghCnv; v: View }) {
  const [allGenes, setAllGenes] = useState(false);
  if (!v.ann || v.ann.state === 'loading' || !v.score) return <div className="px-4 py-2 bg-slate-50 text-slate-500 text-[11px]">Reading the genes of the call…</div>;
  const key = c.kind === 'loss' ? 'hi' : 'ts';
  const order = (x: GeneView) => (x[key] === 3 ? 0 : x.pa?.best === 3 ? 1 : x[key] === 2 ? 2 : (x.loeuf ?? 9) < 0.35 ? 3 : x.pa?.best === 2 ? 4 : 5);
  const genes = v.genes.slice().sort((a, b) => order(a) - order(b) || a.g.start - b.g.start);
  const shownGenes = allGenes ? genes : genes.slice(0, 30);
  const dot = (l: Level) => <span className="inline-block w-2 h-2 rounded-full mr-1.5 align-middle" style={{ background: LEVEL_STYLE[l].color }} />;
  return (
    <div className="px-4 py-2 bg-slate-50 border-t border-b border-accent-200 text-[11px] text-slate-700 cursor-default" onClick={e => e.stopPropagation()}>
      <div className="flex flex-wrap gap-x-6 gap-y-1">
        <div className="min-w-[320px] flex-1">
          <div className="font-semibold text-slate-900">Priority: <span style={{ color: LEVEL_STYLE[v.score.level].color }}>{LEVEL_STYLE[v.score.level].label}</span> <span className="font-normal text-slate-500">(a proposal for the review, not a class)</span></div>
          {v.score.reasons.length ? v.score.reasons.map(r => <div key={r.id}>{dot(r.level)}{r.text}</div>) : <div className="text-slate-500">No protein-coding gene in the call.</div>}
          {v.waiting > 0 && <div className="text-slate-500">PanelApp: {v.waiting} gene{v.waiting > 1 ? 's' : ''} still asked…</div>}
          {v.missing.map(m => <div key={m} className="text-amber-700">Not available: {m}</div>)}
        </div>
        {v.ann.achro.length > 0 && <div className="min-w-[240px]">
          <div className="font-semibold text-slate-900">AChro-Puce</div>
          {v.ann.achro.map(a => <div key={a.id}>{a.name}: <b>{a.cls.join(' / ')}</b> <span className="text-slate-500">({a.how})</span></div>)}
        </div>}
      </div>
      {genes.length > 0 && <table className="mt-1.5 w-full border-collapse bg-white border border-slate-200 rounded">
        <thead className="text-slate-500 text-left"><tr className="border-b border-slate-200">
          <th className="px-2 py-0.5 font-medium">Gene</th><th className="px-2 py-0.5 font-medium">In the call</th>
          <th className="px-2 py-0.5 font-medium" title="ClinGen haploinsufficiency score">HI</th><th className="px-2 py-0.5 font-medium" title="ClinGen triplosensitivity score">TS</th>
          <th className="px-2 py-0.5 font-medium">LOEUF (pLI)</th><th className="px-2 py-0.5 font-medium">PanelApp</th><th className="px-2 py-0.5 font-medium">Diseases (PanelApp)</th>
        </tr></thead>
        <tbody>{shownGenes.map(x => (
          <tr key={x.g.symbol} className="border-b border-slate-100">
            <td className="px-2 py-0.5"><i className="font-semibold">{x.g.symbol}</i> <span className="text-slate-400">{x.g.tx}</span></td>
            <td className="px-2 py-0.5">{x.g.whole ? 'whole' : x.g.exons ? 'cut (exons in)' : 'cut (no exon)'}</td>
            <td className="px-2 py-0.5 font-semibold" style={{ color: x.hi === 3 ? '#b91c1c' : x.hi === 2 ? '#c2410c' : undefined }}>{x.hi ?? ''}</td>
            <td className="px-2 py-0.5 font-semibold" style={{ color: x.ts === 3 ? '#b91c1c' : x.ts === 2 ? '#c2410c' : undefined }}>{x.ts ?? ''}</td>
            <td className="px-2 py-0.5 tabular-nums">{x.loeuf != null ? x.loeuf.toFixed(2) : ''}{x.pli != null ? <span className="text-slate-400"> ({x.pli.toFixed(2)})</span> : null}</td>
            <td className="px-2 py-0.5 whitespace-nowrap">{!x.pa ? <span className="text-slate-400">{v.ann!.genes.indexOf(x.g) < PANELAPP_MAX_GENES ? '…' : 'not asked'}</span>
              : x.pa.best === 0 ? <span className="text-slate-400">in no panel</span>
              : <span className="inline-flex items-center gap-1"><span className="font-semibold" style={{ color: x.pa.best === 3 ? '#16a34a' : x.pa.best === 2 ? '#d97706' : '#dc2626' }}>{x.pa.best === 3 ? 'green' : x.pa.best === 2 ? 'amber' : 'red'}</span>
                {(x.pa.best === 3 ? x.pa.green : x.pa.amber).filter(m => m !== '?').map(m => <MoiBadge key={m} m={m} />)}</span>}</td>
            <td className="px-2 py-0.5 text-slate-600 max-w-[380px] truncate" title={x.pa?.diseases.join('\n')}>{x.pa?.diseases.slice(0, 2).join(' · ')}</td>
          </tr>))}</tbody>
      </table>}
      {genes.length > 30 && <button className="mt-1 text-accent-700 hover:underline" onClick={() => setAllGenes(!allGenes)}>{allGenes ? 'fewer genes' : `all ${genes.length} genes`}</button>}
    </div>
  );
}

/** the rules of the priority, as the score applies them (cghScore.ts) */
function PriorityRules() {
  return (
    <details className="mt-1">
      <summary className="cursor-pointer text-accent-700 font-medium">The priority's rules (a proposal: HIGH, MODERATE, LOW, none)</summary>
      <div className="grid grid-cols-[1fr_1fr] gap-4 mt-1">
        {(['loss', 'gain'] as const).map(k => <div key={k}>
          <div className="font-semibold text-slate-700">{k === 'loss' ? 'Losses' : 'Gains'}</div>
          {SCORE_RULES.filter(r => r.kind === k).map(r => <div key={r.id}><span className="font-bold" style={{ color: LEVEL_STYLE[r.level].color }}>{LEVEL_STYLE[r.level].label}</span> <span className="text-slate-400">{r.id}</span> {r.label}</div>)}
          <div><span className="font-bold" style={{ color: LEVEL_STYLE.none.color }}>none</span>: no protein-coding gene in the call</div>
        </div>)}
      </div>
      <div className="mt-1">A call takes the highest level among the rules it meets. Evidence weighed as in the ACMG / ClinGen technical standards for CNVs (Riggs et al., Genet Med 2020;22:245-57: sections 2 and 3); a dashed badge with ? is incomplete (a source did not answer, PanelApp still asked).</div>
    </details>
  );
}
