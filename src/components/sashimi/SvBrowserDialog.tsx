/**
 * The rearranged genome browser opened from a typed description: an HGVS genomic description or the ISCN sequence
 * form is read into the rebuilt chromosome, which the browser draws with its genes and their annotations.
 */
import { useMemo, useState } from 'react';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import SvNomenclatureBox, { useCytobands } from './SvNomenclatureBox';
import { parseSv, type ParsedSv } from './svNomenclature';
import { segmentColors } from './svColors';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';

const FONT = 'Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const EXAMPLES = [
  { label: 'deletion', text: 'NC_000020.11:g.10489766_10783340del' },
  { label: 'inversion with deletions', text: 'NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]' },
  { label: 'inverted insertion', text: 'NC_000020.11:g.10807815_10807816ins10489766_10783340inv' },
  { label: 'ISCN sequence form', text: 'chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::chr20:g.10,808,806_cen_qter' },
];
const KEY = 'sashimi.svBrowserText';
const load = (): string => { try { return localStorage.getItem(KEY) ?? ''; } catch { return ''; } };
const save = (v: string) => { try { localStorage.setItem(KEY, v); } catch { /* storage blocked */ } };

export default function SvBrowserDialog({ build, loaders, onClose }: { build: string; loaders: Omit<SvGeneLoaders, 'fallback'>; onClose: () => void }) {
  const [text, setText] = useState(load);
  const [result, setResult] = useState<{ parsed: ParsedSv; text: string } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const draw = (t = text) => {
    try {
      const parsed = parseSv(t, build, c => chromLength(build as GenomeBuild, c));
      setResult({ parsed, text: t }); setError(null); save(t);
    } catch (e) { setResult(null); setError((e as Error).message); }
  };
  const chrom = result?.parsed.chrom ?? '';
  const bands = useCytobands(chrom || 'chr1', loaders.getCytobands);
  const colors = useMemo(() => (result ? segmentColors(result.parsed.rec) : new Map<number, string>()), [result]);
  const full = useMemo<SvGeneLoaders>(() => ({ ...loaders, fallback: [] }), [loaders]);
  const rec = result?.parsed.rec;
  return (
    <div className="fixed inset-0 z-[70] bg-slate-900/80 backdrop-blur-sm flex items-start justify-center overflow-y-auto p-4" onMouseDown={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-[1240px] text-gray-900" onMouseDown={e => e.stopPropagation()} style={{ fontFamily: FONT }}>
        <div className="flex items-start justify-between px-5 pt-4 pb-2">
          <div>
            <div className="text-base font-bold">SV browser <span className="font-normal text-gray-500">· {build}</span></div>
            <div className="text-xs text-gray-500">the genes of a rearrangement typed as HGVS (NC_… or chrN:g., del, dup, inv, ins, delins, allele form […;…]) or as the ISCN sequence form (chrN:g.…::chrN:g.…)</div>
          </div>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-2xl leading-none px-1" title="Close">×</button>
        </div>
        <div className="px-5 pb-3">
          <div className="flex gap-2">
            <input value={text} onChange={e => setText(e.target.value)} onKeyDown={e => { if (e.key === 'Enter') draw(); }} spellCheck={false} autoFocus
              placeholder="NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]"
              className="flex-1 font-mono text-[12.5px] h-9 px-3 rounded-md border border-gray-300 focus:outline-none focus:ring-2 focus:ring-accent-300" />
            <button onClick={() => draw()} className="h-9 px-4 rounded-md bg-accent-600 text-white text-sm font-semibold hover:bg-accent-700">Draw</button>
          </div>
          <div className="text-[11px] text-gray-500 mt-1">
            examples: {EXAMPLES.map((x, i) => <span key={x.label}>{i ? ' · ' : ''}<button className="text-accent-700 hover:underline" onClick={() => { setText(x.text); draw(x.text); }}>{x.label}</button></span>)}
            {' '}· band-level ISCN has no positions and is not drawn; one chromosome at a time (no translocation yet)
          </div>
          {error && <div className="text-[12px] text-red-700 mt-2">{error}</div>}
          {result && rec && (
            <div className="mt-2 grid grid-cols-[1fr_1fr] gap-4 text-xs">
              <div>
                <div className="text-[10.5px] font-semibold text-gray-700">Rebuilt chromosome</div>
                <div className="font-mono text-[11px] bg-gray-50 border border-gray-200 rounded px-2 py-1 mt-0.5">{rec.ders[0].name} = {rec.ders[0].pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ')}</div>
                {rec.lost.length > 0 && <div className="text-[11px] text-red-700 mt-0.5">Lost from the rearranged homologue: {rec.lost.map(s => s.letter).join(', ')}.</div>}
                {result.parsed.notes.map((n, i) => <div key={i} className="text-[10.5px] text-gray-500">{n}</div>)}
                <div className="text-[10.5px] text-gray-500 mt-0.5">Drawn as heterozygous, 2 copies on the flanks; no read depth (open the cartoon on a sample for that).</div>
              </div>
              <SvNomenclatureBox rec={rec} build={build} bands={bands} />
            </div>
          )}
        </div>
        {result && rec && (
          <SvGenomeBrowser key={result.text} rec={rec} chrom={chrom} window={result.parsed.window} segColors={colors} flankCopies={2} zygosity="het" build={build} loaders={full} rebuild={result.parsed.rebuild} />
        )}
      </div>
    </div>
  );
}
