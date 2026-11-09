/**
 * The array-CGH page, one level under the SV browser: the calls of the samples added (interval tables, BED, zips), one
 * tab per patient, and under them the SV browser's own browser drawing the call clicked (DEL / DUP on the file's build,
 * GRCh37), with its tracks. Back to the SV browser from the bar (or Esc).
 */
import { memo, useMemo, useState } from 'react';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import { parseSv, type ParsedSv } from './svNomenclature';
import { segmentColors } from './svColors';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';
import { CghFrame, CghProbeDialog } from './CghPanel';
import { cgh, useCgh } from './cghStore';
import type { CghCnv, CghFile } from './cghIntervals';
import { bandOfCall } from './cghAnnotate';

/** the call's array ISCN: the file's, else written from its band, positions and copies when the file gives both; else null */
const iscnOf = (c: CghCnv, build: GenomeBuild, band: string): string | null => c.iscn ?? (band && (c.copies != null || c.amp)
  ? `arr[${build}] ${c.chrom.replace(/^chr/, '')}${band.replace(/\s+/g, '').replace(/-/, '')}(${c.start}_${c.end})${c.copies != null ? `x${c.copies}${c.copiesHi != null ? `~${c.copiesHi}` : ''}` : 'amp'}`
  : null);

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
type Loaders = Omit<SvGeneLoaders, 'fallback'>;
type Drawn = { cnv: CghCnv; build: GenomeBuild; parsed: ParsedSv };
/** a call as the SV browser reads it: DEL / DUP over its positions, on its file's build (GRCh37 when it states none) */
function drawOf(c: CghCnv, f: CghFile): Drawn {
  const b: GenomeBuild = f.build ?? 'GRCh37';
  return { cnv: c, build: b, parsed: parseSv(`${c.chrom}:g.${c.start}_${c.end}${c.kind === 'loss' ? 'del' : 'dup'}`, b, ch => chromLength(b, ch)) };
}
/** the browser redrawn only when its own props change, not at each answer of the annotations */
const Browser = memo(SvGenomeBrowser);

/** The SV browser's way in: the button (with the number of samples) and the page over it. Each listens to the array-CGH
 * state by itself, so that the SV browser is not redrawn when the calls are annotated. */
export function CghOpenButton() {
  const s = useCgh();
  const n = s.files.reduce((k, f) => k + f.samples.length, 0);
  return (
    <button onClick={() => cgh.setOpen(true)} className="ml-auto inline-flex items-center gap-1.5 h-[28px] px-2.5 rounded-[6px] border text-[12px] font-medium bg-white border-slate-200 text-slate-700 hover:border-accent-300 hover:bg-accent-50"
      title="Array-CGH: the interval tables of your patients (CytoGenomics .txt, BED or zip), each call classified and drawn in this browser (GRCh37)">
      Array-CGH{n ? <span className="text-[10.5px] font-bold leading-4 px-1.5 rounded-full bg-accent-600 text-white">{n}</span> : null} ▸</button>
  );
}
export function CghPageSlot({ loadersFor }: { loadersFor: (b: GenomeBuild) => Loaders }) {
  const open = useCgh().open;
  return open ? <CghPage loadersFor={loadersFor} onBack={() => cgh.setOpen(false)} /> : null;
}

export default function CghPage({ loadersFor, onBack }: { loadersFor: (b: GenomeBuild) => Loaders; onBack: () => void }) {
  const s = useCgh();
  /** the call drawn: the one clicked, or, the page opened again, the one clicked before */
  const [drawn, setDrawn] = useState<Drawn | null>(() => {
    const id = cgh.get().shown, f = id ? s.files.find(fl => fl.samples.some(x => x.cnvs.some(c => c.id === id))) : undefined;
    const c = f?.samples.flatMap(x => x.cnvs).find(x => x.id === id);
    try { return c && f ? drawOf(c, f) : null; } catch { return null; }
  });
  const [error, setError] = useState<string | null>(null);
  const [over, setOver] = useState(false);
  /** a call clicked: DEL / DUP over its positions, on the file's build (GRCh37 when it states none) */
  const show = (c: CghCnv, f: CghFile) => {
    try { setDrawn(drawOf(c, f)); setError(null); } catch (e) { setError(`${c.chrom}:${c.start}-${c.end}: ${(e as Error).message}`); }
  };
  const loaders = useMemo(() => (drawn ? loadersFor(drawn.build) : null), [drawn?.build, loadersFor]); // eslint-disable-line react-hooks/exhaustive-deps
  const full = useMemo<SvGeneLoaders | null>(() => (loaders ? { ...loaders, fallback: [] } : null), [loaders]);
  const rec = drawn?.parsed.rec;
  const colors = useMemo(() => (rec ? segmentColors(rec) : new Map<number, string>()), [rec]);
  // the calls this page lists are in the store: a call drawn from a sample taken off stays until another is clicked
  // the call drawn, while its patient is on the list
  const c = drawn && s.files.some(f => f.samples.some(x => x.key === drawn.cnv.sampleKey)) ? drawn.cnv : undefined;
  // a male's chrX and chrY: one copy, the call on it (no other homologue)
  const one = !!c && s.sex[c.sampleKey] === 'male' && (c.chrom === 'chrX' || c.chrom === 'chrY');
  return (
    <div className="fixed inset-0 z-[72] bg-page overflow-y-auto text-slate-900" style={{ fontFamily: FONT }} data-cgh-page
      onDragOver={e => { e.preventDefault(); e.stopPropagation(); if (!over) setOver(true); }}
      onDragLeave={e => { if (e.currentTarget === e.target) setOver(false); }}
      onDrop={e => { e.preventDefault(); e.stopPropagation(); setOver(false); if (e.dataTransfer.files.length) void cgh.add(e.dataTransfer.files); }}>
      <div className="w-full min-h-full pb-4">
        {/* the SV browser's bar, one level down */}
        <div className="sticky top-0 z-30 flex items-center gap-3 px-4 h-12 bg-bar border-b border-bar-line">
          <button onClick={onBack} className="inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="Back to the SV browser (Esc)">← SV browser</button>
          <span className="text-bar-muted">/</span>
          <span className="text-sm font-semibold text-bar-text">Array-CGH</span>
          <span className="ml-auto text-[11px] text-bar-muted">Esc: back to the SV browser</span>
        </div>
        {!s.files.length && (
          <div className={`mx-4 mt-3 rounded-xl border-2 border-dashed ${over ? 'border-accent-400 bg-accent-50' : 'border-slate-300 bg-white'} px-6 py-10 text-center`}>
            <div className="text-base font-bold">Array-CGH samples</div>
            <div className="text-xs text-slate-500 mt-1 max-w-2xl mx-auto">Drop here, or add, the interval table of each patient exported by CytoGenomics (.txt), BED files, or zip files (one patient each: the interval tables inside are read, the other files skipped). Read in this page, never sent anywhere.</div>
            <div className="mt-4 flex justify-center"><AddInline /></div>
          </div>
        )}
        <CghFrame onShow={show} dropping={over} loadersFor={loadersFor} />
        {error && <div className="mx-4 mt-2 text-[12px] text-red-700">{error}</div>}
        {drawn && rec && full && c && (
          <div className="mx-4 mt-3 bg-white border border-slate-200 rounded-xl shadow-sm pb-2 overflow-hidden">
            <div className="px-5 pt-2 flex flex-wrap items-baseline gap-x-2 text-xs">
              <span className="text-sm font-bold">{c.sample}</span>
              {!iscnOf(c, drawn.build, bandOfCall(c, s.ann[c.id])) && bandOfCall(c, s.ann[c.id]) && <span className="font-semibold">{c.chrom.replace(/^chr/, '')}{bandOfCall(c, s.ann[c.id]).replace(/\s+/g, '')}</span>}
              {iscnOf(c, drawn.build, bandOfCall(c, s.ann[c.id]))
                ? <span className="font-mono text-[12px] text-slate-800" title={c.iscn ? 'ISCN as in the file' : 'ISCN written from the band (the file\'s, else the cytogenetic bands), the positions and the copies of the file (it gives no ISCN)'}>{iscnOf(c, drawn.build, bandOfCall(c, s.ann[c.id]))}</span>
                : <span className="font-mono text-[12px] text-slate-800" title="The file gives no ISCN, and not both the band and the copy number to write one">{c.chrom}:{c.start.toLocaleString('en-US')}-{c.end.toLocaleString('en-US')} <span className="font-sans text-slate-500">({c.kind}, no ISCN in the file)</span></span>}
              <span className="text-slate-500">· {c.sizeKb >= 1000 ? `${(c.sizeKb / 1000).toFixed(2)} Mb` : `${c.sizeKb.toFixed(1)} kb`} · {c.probes != null ? `${c.probes} probes` : 'probes not given'}{c.logRatio != null ? ` · log ratio ${c.logRatio.toFixed(2)}` : ''}{one ? ' · drawn with one copy (male)' : c.copies === 0 ? ' · drawn homozygous' : ''}{c.kind === 'gain' ? ' · extra copy drawn in tandem' : ''}</span>
            </div>
            <div className="[&>div:first-child]:mt-2">
              <Browser key={`${drawn.build}|${c.id}`} rec={rec} chrom={c.chrom} window={drawn.parsed.window} segColors={colors} flankCopies={one ? 1 : 2}
                zygosity={!one && c.copies === 0 ? 'hom' : 'het'} build={drawn.build} loaders={full} rebuild={drawn.parsed.rebuild} strongEvent />
            </div>
          </div>
        )}
      </div>
      <CghProbeDialog />
    </div>
  );
}

/** the empty page's button */
function AddInline() {
  return (
    <label className="inline-flex items-center gap-1.5 h-[34px] px-4 rounded-[6px] bg-accent-600 text-white text-[12.5px] font-semibold hover:bg-accent-700 cursor-pointer">
      Add samples (.txt, .bed, .zip)
      <input data-cgh-input type="file" multiple accept=".txt,.tsv,.bed,.xls,.zip,text/plain,application/zip" className="hidden"
        onChange={e => { const fs = e.target.files; if (fs?.length) void cgh.add(fs); e.target.value = ''; }} />
    </label>
  );
}
