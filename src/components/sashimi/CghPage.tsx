/**
 * The array-CGH page, one level under the SV browser: the calls of the samples added (interval tables, BED, zips), one
 * tab per patient, and under them the SV browser's own browser drawing the call clicked (DEL / DUP on the file's build,
 * GRCh37), with its tracks. Esc folds what is open (back to the list, the whole genome folded); on that base view, Esc or
 * the bar's ← SV browser asks before going back to the SV browser.
 */
import { memo, useEffect, useMemo, useRef, useState } from 'react';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import { parseSv, type ParsedSv } from './svNomenclature';
import { segmentColors } from './svColors';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';
import { CghFrame, CghProbeDialog, CghExamplesMenu } from './CghPanel';
import { cgh, useCgh, CGH_ACCEPT, CGH_FILES_TEXT } from './cghStore';
import type { CghCnv, CghFile } from './cghIntervals';
import { bandOfCall } from './cghAnnotate';
import { clinGenTrack } from './CghClinGenTrack';
import { CghGuide, CghGuideButton } from './CghGuide';
import { CghProbePanel, callMargin } from './CghProbePanel';
import { DLRS_TIP, dlrsWord } from './cghProbes';
import { estCopies, estTip, parOf } from './cghMosaic';
import { CghGenomePanel } from './CghGenome';
import { CallIscnText, iscnOfCall } from './CghIscnView';


const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
type Loaders = Omit<SvGeneLoaders, 'fallback'>;
type Drawn = { cnv: CghCnv; build: GenomeBuild; parsed: ParsedSv };
/** a call as the SV browser reads it: DEL / DUP over its positions, on its file's build (GRCh37 when it states none) */
function drawOf(c: CghCnv, f: CghFile): Drawn {
  const b: GenomeBuild = f.build ?? 'GRCh37', len = chromLength(b, c.chrom) ?? c.end * 2;
  const parsed = parseSv(`${c.chrom}:g.${c.start}_${c.end}${c.kind === 'loss' ? 'del' : 'dup'}`, b, () => len);
  // centred on the call with M = max(50 % of its size, 100 kb) each side: the browser pads its window by 12 % more, so m + 0.12 (size + 2m) = M
  const size = c.end - c.start, m = Math.max(0, (callMargin(c) - 0.12 * size) / 1.24);
  const win = { start: Math.max(0, Math.round(c.start - m)), end: Math.min(len, Math.round(c.end + m)) };
  return { cnv: c, build: b, parsed: { ...parsed, window: win, rec: parsed.rebuild(win) } };
}
/** the browser redrawn only when its own props change, not at each answer of the annotations */
const Browser = memo(SvGenomeBrowser);
// one track object per build: the memoised browser redraws only when it changes

/** The SV browser's way in: the button (with the number of samples) and the page over it. Each listens to the array-CGH
 * state by itself, so that the SV browser is not redrawn when the calls are annotated. */
export function CghOpenButton() {
  const s = useCgh();
  const n = s.files.reduce((k, f) => k + f.samples.length, 0);
  return (
    <button onClick={() => cgh.setOpen(true)} className="ml-auto inline-flex items-center gap-1.5 h-[28px] px-2.5 rounded-[6px] border text-[12px] font-medium bg-white border-slate-200 text-slate-700 hover:border-accent-300 hover:bg-accent-50"
      title="Array-CGH: the calls of your patients (CytoGenomics interval tables, other tables, .xlsx, BED, VCF or zip), each call classified and drawn in this browser (GRCh37)">
      Array-CGH{n ? <span className="text-[10.5px] font-bold leading-4 px-1.5 rounded-full bg-accent-600 text-white">{n}</span> : null} ▸</button>
  );
}
export function CghPageSlot({ loadersFor }: { loadersFor: (b: GenomeBuild) => Loaders }) {
  const open = useCgh().open;
  return open ? <CghPage loadersFor={loadersFor} onBack={() => cgh.askLeave()} /> : null;
}

export default function CghPage({ loadersFor, onBack }: { loadersFor: (b: GenomeBuild) => Loaders; onBack: () => void }) {
  const s = useCgh();
  /** the call drawn: the one clicked, or, the page opened again, the one clicked before */
  const [drawn, setDrawn] = useState<Drawn | null>(() => {
    const id = cgh.get().shown, f = id ? s.files.find(fl => fl.samples.some(x => x.cnvs.some(c => c.id === id))) : undefined;
    const c = f?.samples.flatMap(x => x.cnvs).find(x => x.id === id);
    try { return c && f ? drawOf(c, f) : null; } catch { return null; }
  });
  // the whole genome stands where the call window does: open while no call is, folded to a strip once one is (and the other way round)
  const [genome, setGenome] = useState(() => !drawn);
  const [error, setError] = useState<string | null>(null);
  const [over, setOver] = useState(false);
  // the ClinGen track folded by default: HI 3 / TS 3 regions only, the recurrent list and its LCRs on demand
  const [recOpen, setRecOpen] = useState(false);
  const track = useMemo(() => clinGenTrack(drawn?.build ?? 'GRCh37', recOpen, () => setRecOpen(o => !o)), [drawn?.build, recOpen]);
  /** a call clicked: DEL / DUP over its positions, on the file's build (GRCh37 when it states none) */
  const show = (c: CghCnv, f: CghFile) => {
    setGenome(false);
    try { setDrawn(drawOf(c, f)); setError(null); } catch (e) { setError(`${c.chrom}:${c.start}-${c.end}: ${(e as Error).message}`); }
  };
  const loaders = useMemo(() => (drawn ? loadersFor(drawn.build) : null), [drawn?.build, loadersFor]); // eslint-disable-line react-hooks/exhaustive-deps
  const full = useMemo<SvGeneLoaders | null>(() => (loaders ? { ...loaders, fallback: [] } : null), [loaders]);
  const rec = drawn?.parsed.rec;
  const colors = useMemo(() => (rec ? segmentColors(rec) : new Map<number, string>()), [rec]);
  // the call drawn, while its patient is on the list
  const c = drawn && s.files.some(f => f.samples.some(x => x.key === drawn.cnv.sampleKey)) ? drawn.cnv : undefined;
  // a male's chrX and chrY outside the PARs (of the call's build): one copy, the call on it (no other homologue)
  const one = !!c && s.sex[c.sampleKey] === 'male' && (c.chrom === 'chrX' || c.chrom === 'chrY') && !(parOf(drawn?.build)[c.chrom] ?? []).some(([a, b]) => c.start >= a && c.end <= b);
  const band = c ? bandOfCall(c, s.ann[c.id]) : '', isc = c && drawn ? iscnOfCall(s, c, drawn.build) : null, iscn = isc?.text ?? isc?.file ?? null;
  // the probes of the call's patient (read on first view) and its copy number from the log ratio when the file gives none
  const probes = c && cgh.hasProbes(c.sampleKey) ? cgh.probes(c.sampleKey) : null, est = c ? estCopies(c, s.sex[c.sampleKey], drawn?.build) : null;
  // Esc: what is open beyond the base view (told to the store), folded when the store says so
  const busy = useRef(false);
  busy.current = (genome && s.files.length > 0) || !!c || !!error;
  useEffect(() => { cgh.setViewBusy(() => busy.current); return () => cgh.setViewBusy(null); }, []);
  const reset0 = useRef(s.viewReset);
  useEffect(() => {
    if (s.viewReset === reset0.current) return;
    reset0.current = s.viewReset; setGenome(false); setDrawn(null); setError(null);
  }, [s.viewReset]);
  return (
    <div className="fixed inset-0 z-[72] bg-page overflow-y-auto text-slate-900" style={{ fontFamily: FONT }} data-cgh-page
      onDragOver={e => { e.preventDefault(); e.stopPropagation(); if (!over) setOver(true); }}
      onDragLeave={e => { if (e.currentTarget === e.target) setOver(false); }}
      onDrop={e => { e.preventDefault(); e.stopPropagation(); setOver(false); if (e.dataTransfer.files.length) void cgh.add(e.dataTransfer.files); }}>
      <div className="w-full min-h-full pb-4">
        {/* the SV browser's bar, one level down */}
        <div className="sticky top-0 z-30 flex items-center gap-3 px-4 h-12 bg-bar border-b border-bar-line">
          <button onClick={onBack} className="inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="Back to the SV browser (asks first; the patients stay in this tab)">← SV browser</button>
          <span className="text-bar-muted">/</span>
          <span className="text-sm font-semibold text-bar-text">Array-CGH</span>
          <CghGuideButton />
          <span className="text-[11px] text-bar-muted" data-cgh-esc-hint title="Esc folds what is open (the call, its detail, the whole genome) back to the list; Esc again asks to go back to the SV browser">Esc: clear the view · twice: leave</span>
        </div>
        {!s.files.length && (
          <div className={`mx-4 mt-3 rounded-xl border-2 border-dashed ${over ? 'border-accent-400 bg-accent-50' : 'border-slate-300 bg-white'} px-6 py-10 text-center`}>
            <div className="text-base font-bold">Array-CGH samples</div>
            <div className="text-xs text-slate-500 mt-1 max-w-2xl mx-auto">Drop here, or add, {CGH_FILES_TEXT}. Read in this page, never sent anywhere.</div>
            <div className="mt-4 flex justify-center gap-2"><AddInline /><CghExamplesMenu /></div>
          </div>
        )}
        <CghFrame onShow={show} dropping={over} loadersFor={loadersFor} />
        {error && <div className="mx-4 mt-2 text-[12px] text-red-700">{error}</div>}
        {s.files.length > 0 && <CghGenomePanel open={genome} onToggle={() => setGenome(g => !g)} onOpen={show} loadersFor={loadersFor} />}
        {genome && drawn && rec && full && c && <button onClick={() => setGenome(false)} data-cgh-call-strip className="mx-4 mt-3 block w-[calc(100%-2rem)] text-left px-5 py-1.5 bg-white border border-slate-200 rounded-xl shadow-sm text-xs font-semibold text-slate-800 hover:border-accent-300"
          title="Show the call window again (folds the whole genome)">▸ Call — {c.sample} · {c.chrom}:{c.start.toLocaleString('en-US')}-{c.end.toLocaleString('en-US')} ({c.kind})</button>}
        {!genome && drawn && rec && full && c && (
          <div data-cgh-browser className="mx-4 mt-3 bg-white border border-slate-200 rounded-xl shadow-sm pb-2 overflow-hidden">
            <div className="px-5 pt-2 flex flex-wrap items-baseline gap-x-2 text-xs">
              <span className="text-sm font-bold">{c.sample}</span>
              {!iscn && band && <span className="font-semibold">{c.chrom.replace(/^chr/, '')}{band.replace(/\s+/g, '')}</span>}
              {iscn
                ? <span title="Array ISCN as in ISCN 2024, written from the band (the cytogenetic bands at sub-band precision, else the file's), the positions and the copy number (the file's, else from the log ratio), with the mosaic fraction when the log ratio says mosaic; the file's own ISCN beside when it differs"><CallIscnText i={isc!} /></span>
                : <span className="font-mono text-[12px] text-slate-800" title="The file gives no ISCN, and not both the band and the copy number to write one">{c.chrom}:{c.start.toLocaleString('en-US')}-{c.end.toLocaleString('en-US')} <span className="font-sans text-slate-500">{c.fromProbes ? <i>({c.kind}, from probes: not in the file, not exported)</i> : `(${c.kind}, no ISCN in the file)`}</span></span>}
              <span className="text-slate-500">· {c.sizeKb >= 1000 ? `${(c.sizeKb / 1000).toFixed(2)} Mb` : `${c.sizeKb.toFixed(1)} kb`} {c.probes != null ? ` · ${c.probes} probes` : ''}{c.logRatio != null ? ` · log ratio ${c.logRatio.toFixed(2)}` : ''}{c.pval != null ? ` · p ${c.pval.toExponential(1)}` : ''}{one ? ' · drawn with one copy (male)' : c.copies === 0 ? ' · drawn homozygous' : ''}{c.kind === 'gain' ? ' · extra copy drawn in tandem' : ''}</span>
              {est && <span className="font-semibold text-slate-700" title={estTip(est, c.logRatio!)} data-cgh-est>copies {est.text}</span>}
              {probes?.dlrs != null && <span className="text-slate-500" title={DLRS_TIP} data-cgh-dlrs>· DLRS <b className="text-slate-800">{probes.dlrs.toFixed(2)}</b> ({dlrsWord(probes.dlrs)})</span>}
            </div>
            {probes && <CghProbePanel key={c.id} c={c} chromLen={chromLength(drawn.build, c.chrom) ?? 0} segs={cgh.probeCallsOf(c.sampleKey)?.segs.get(c.chrom)} getBands={full.getCytobands} />}
            <div className="[&>div:first-child]:mt-2" data-cgh-svb>
              <Browser key={`${drawn.build}|${c.id}`} rec={rec} chrom={c.chrom} window={drawn.parsed.window} segColors={colors} flankCopies={one ? 1 : 2}
                zygosity={!one && c.copies === 0 ? 'hom' : 'het'} build={drawn.build} loaders={full} rebuild={drawn.parsed.rebuild} strongEvent extraTrack={track} />
            </div>
          </div>
        )}
      </div>
      <CghProbeDialog />
      <CghGuide />
      {s.leaving && <LeaveDialog />}
    </div>
  );
}

/** "Leave Array-CGH?": Cancel focused (Esc = Cancel, Enter = the focused button); the patients stay in the store */
function LeaveDialog() {
  const s = useCgh();
  const cancelBtn = useRef<HTMLButtonElement>(null);
  const n = Object.keys(s.dirtyKeys).length;
  const cancel = () => cgh.setLeaving(false);
  useEffect(() => {
    cancelBtn.current?.focus();
    // capture: before the page's own Esc, wherever the focus is
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') { e.stopPropagation(); e.preventDefault(); cancel(); } };
    window.addEventListener('keydown', onKey, true);
    return () => window.removeEventListener('keydown', onKey, true);
  }, []);
  const btn = 'h-[30px] px-4 rounded-[6px] text-[12.5px] font-semibold focus:outline-none focus-visible:ring-2 focus-visible:ring-accent-400 focus-visible:ring-offset-1';
  return (
    <div className="fixed inset-0 z-[90] bg-slate-900/30 flex items-center justify-center p-4" onClick={e => { if (e.target === e.currentTarget) cancel(); }} data-cgh-leave>
      <div role="alertdialog" aria-modal="true" aria-labelledby="cgh-leave-title" aria-describedby="cgh-leave-note" className="w-full max-w-[420px] bg-white rounded-xl shadow-xl border border-slate-200 px-5 py-4">
        <div id="cgh-leave-title" className="text-[14px] font-bold text-slate-900">Leave Array-CGH and go back to the SV browser?</div>
        <div id="cgh-leave-note" className="mt-1.5 text-[12px] text-slate-600 leading-5">
          The loaded patients, their classes and notes stay in this tab (Array-CGH ▸ brings them back) until the tab is closed or reloaded; export first to keep them in a file.
          {n > 0 && <div className="mt-1 font-semibold text-amber-700" data-cgh-leave-dirty>{n} patient{n > 1 ? 's' : ''} with unexported changes.</div>}
        </div>
        <div className="mt-4 flex justify-end gap-2">
          <button ref={cancelBtn} onClick={cancel} className={`${btn} border border-slate-300 bg-white text-slate-700 hover:bg-slate-50`} data-cgh-leave-cancel>Cancel</button>
          <button onClick={() => cgh.setOpen(false)} className={`${btn} bg-accent-600 text-white hover:bg-accent-700`} data-cgh-leave-ok>Leave</button>
        </div>
      </div>
    </div>
  );
}

/** the empty page's button */
function AddInline() {
  return (
    <label className="inline-flex items-center gap-1.5 h-[34px] px-4 rounded-[6px] bg-accent-600 text-white text-[12.5px] font-semibold hover:bg-accent-700 cursor-pointer">
      Add samples (.txt, .csv, .xlsx, .bed, .vcf, .zip)
      <input data-cgh-input type="file" multiple accept={CGH_ACCEPT} className="hidden"
        onChange={e => { const fs = e.target.files; if (fs?.length) void cgh.add(fs); e.target.value = ''; }} />
    </label>
  );
}
