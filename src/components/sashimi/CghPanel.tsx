/**
 * Array-CGH in the SV browser: the button that adds samples (interval tables from CytoGenomics, or BED files), the
 * question of the minimum number of probes, and the frame listing each sample's calls between the search and the
 * browser. A call clicked is drawn in the browser as a deletion or a duplication, with the browser's own tracks.
 */
import { useRef, useState } from 'react';
import { Icon, Stepper } from './controls';
import { cgh, useCgh, CGH_MIN_PROBES_DEFAULT, type CghState } from './cghStore';

/** "＋ Add samples": the files picked (several at once), read in the page. */
export function CghAddButton({ compact }: { compact?: boolean }) {
  const input = useRef<HTMLInputElement>(null);
  return (
    <label className={`inline-flex items-center gap-1.5 h-[30px] px-2.5 rounded-[6px] border text-[12.5px] font-medium cursor-pointer whitespace-nowrap bg-white border-slate-200 text-slate-700 hover:border-accent-300 hover:bg-accent-50 ${compact ? 'h-[26px] px-2 text-xs rounded-lg' : ''}`}
      title={'Array-CGH calls of one or more patients: the interval table of each sample exported by CytoGenomics (.txt: Chromosome, Start, Stop, Cytoband, Gene Name, Size(kb), Type, #Probes, Mean Log Ratio, ISCN, Suppress, Classification), or a BED (chrom, start, end, a name saying del / dup or x1 / x3).\nRead in this page, never sent anywhere.'}>
      <Icon name="upload" size={compact ? 13 : 15} />{compact ? 'Add samples' : 'Add array-CGH samples'}
      <input ref={input} data-cgh-input type="file" multiple accept=".txt,.tsv,.bed,.xls,text/plain" className="hidden"
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
