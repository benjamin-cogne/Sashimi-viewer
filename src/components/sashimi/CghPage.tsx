/**
 * The array-CGH page, one level under the SV browser: the calls of the samples added (interval tables, BED, zips), one
 * tab per patient, and under them the SV browser's own browser drawing the call clicked (DEL / DUP on the file's build,
 * GRCh37), with its tracks. Back to the SV browser from the bar (or Esc).
 */
import { useMemo, useState } from 'react';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import SvNomenclatureBox, { useCytobands } from './SvNomenclatureBox';
import { parseSv, type ParsedSv } from './svNomenclature';
import { segmentColors } from './svColors';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';
import { CghFrame, CghProbeDialog } from './CghPanel';
import { cgh, useCgh } from './cghStore';
import { copyText, type CghCnv, type CghFile } from './cghIntervals';

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
type Loaders = Omit<SvGeneLoaders, 'fallback'>;

export default function CghPage({ loadersFor, onBack }: { loadersFor: (b: GenomeBuild) => Loaders; onBack: () => void }) {
  const s = useCgh();
  const [drawn, setDrawn] = useState<{ cnv: CghCnv; build: GenomeBuild; parsed: ParsedSv } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [over, setOver] = useState(false);
  /** a call clicked: DEL / DUP over its positions, on the file's build (GRCh37 when it states none) */
  const show = (c: CghCnv, f: CghFile) => {
    const b: GenomeBuild = f.build ?? 'GRCh37';
    try {
      const parsed = parseSv(`${c.chrom}:g.${c.start}_${c.end}${c.kind === 'loss' ? 'del' : 'dup'}`, b, ch => chromLength(b, ch));
      setDrawn({ cnv: c, build: b, parsed }); setError(null);
    } catch (e) { setError(`${c.chrom}:${c.start}-${c.end}: ${(e as Error).message}`); }
  };
  const loaders = useMemo(() => (drawn ? loadersFor(drawn.build) : null), [drawn?.build, loadersFor]); // eslint-disable-line react-hooks/exhaustive-deps
  const full = useMemo<SvGeneLoaders | null>(() => (loaders ? { ...loaders, fallback: [] } : null), [loaders]);
  const rec = drawn?.parsed.rec;
  const colors = useMemo(() => (rec ? segmentColors(rec) : new Map<number, string>()), [rec]);
  const bands = useCytobands(drawn?.cnv.chrom ?? 'chr1', loaders?.getCytobands);
  // the calls this page lists are in the store: a call drawn from a sample taken off stays until another is clicked
  const c = drawn?.cnv;
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
            <div className="px-5 pt-3 grid grid-cols-[1fr_1fr] gap-4 text-xs">
              <div>
                <div className="text-sm font-bold">{c.sample} · {c.chrom.replace(/^chr/, '')}{c.cytoband.replace(/\s+/g, '')} {c.kind === 'loss' ? 'loss' : 'gain'} {copyText(c)}</div>
                <div className="text-[11px] text-slate-600 mt-0.5">{c.chrom}:{c.start.toLocaleString('en-US')}-{c.end.toLocaleString('en-US')} ({drawn.build}) · {c.probes ?? '?'} probes{c.logRatio != null ? ` · log ratio ${c.logRatio.toFixed(2)}` : ''}{c.iscn ? <> · <span className="font-mono">{c.iscn}</span></> : null}</div>
                <div className="font-mono text-[11px] bg-slate-50 border border-slate-200 rounded px-2 py-1 mt-1">{rec.ders[0].name} = {rec.ders[0].pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ')}</div>
                <div className="text-[10.5px] text-slate-500 mt-0.5">{one ? `Drawn with one copy (${c.chrom} of a male)` : c.copies === 0 ? 'Drawn as homozygous (x0: both copies lost), 2 copies on the flanks' : 'Drawn as heterozygous, 2 copies on the flanks'}{c.kind === 'gain' ? '; the extra copy drawn in tandem (its place is not known from the array)' : ''}.</div>
              </div>
              <SvNomenclatureBox rec={rec} build={drawn.build} bands={bands} />
            </div>
            <div className="[&>div:first-child]:mt-2">
              <SvGenomeBrowser key={`${drawn.build}|${c.id}`} rec={rec} chrom={c.chrom} window={drawn.parsed.window} segColors={colors} flankCopies={one ? 1 : 2}
                zygosity={!one && c.copies === 0 ? 'hom' : 'het'} build={drawn.build} loaders={full} rebuild={drawn.parsed.rebuild} />
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
