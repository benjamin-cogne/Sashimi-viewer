/**
 * The SV browser opened from the toolbar: a rearrangement typed as HGVS or as the ISCN sequence form is rebuilt and
 * drawn with its genes; coordinates or a gene name go to that place of the genome (with the rearrangement when it is on
 * the same chromosome, else the reference alone). GRCh38 by default, GRCh37 on demand: an `NC_` accession's version
 * picks its build by itself.
 */
import { useEffect, useMemo, useState } from 'react';
import SvGenomeBrowser, { type SvGeneLoaders } from './SvGenomeBrowser';
import SvNomenclatureBox, { useCytobands } from './SvNomenclatureBox';
import { parseSv, type ParsedSv } from './svNomenclature';
import { segmentColors } from './svColors';
import { parseLocus } from './geometry';
import { buildOfNc } from './knownVariants';
import { chromLength } from '../../standalone/assembly';
import type { GenomeBuild } from '../../standalone/ensembl';

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const EXAMPLES = [
  { label: 'deletion', text: 'NC_000020.11:g.10489766_10783340del' },
  { label: 'inversion with deletions', text: 'NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]' },
  { label: 'inverted insertion', text: 'NC_000020.11:g.10807815_10807816ins10489766_10783340inv' },
  { label: 'ISCN sequence form', text: 'chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::chr20:g.10,808,806_cen_qter' },
  { label: 'coordinates', text: 'chr20:10,300,000-11,000,000' },
  { label: 'a gene', text: 'JAG1' },
];
const KEY = 'sashimi.svBrowserText', BUILD_KEY = 'sashimi.svBuild';
const load = (k: string): string => { try { return localStorage.getItem(k) ?? ''; } catch { return ''; } };
const save = (k: string, v: string) => { try { localStorage.setItem(k, v); } catch { /* storage blocked */ } };
/** A single position opens this much around it. */
const AROUND_POS = 50_000;

type Shown = { kind: 'sv'; parsed: ParsedSv; text: string } | { kind: 'region'; chrom: string; start: number; end: number };
type Loaders = Omit<SvGeneLoaders, 'fallback'>;

export default function SvBrowserDialog({ build: sessionBuild, loaders: sessionLoaders, loadersFor, onClose }: {
  build: string; loaders: Loaders;
  /** loaders for another build (the SV browser's own switch); without it, the session's build only */
  loadersFor?: (build: GenomeBuild) => Loaders;
  onClose: () => void;
}) {
  const [build, setBuild] = useState<GenomeBuild>(() => (loadersFor ? ((load(BUILD_KEY) as GenomeBuild) || 'GRCh38') : (sessionBuild as GenomeBuild)));
  const loaders = useMemo(() => (loadersFor ? loadersFor(build) : sessionLoaders), [build, loadersFor, sessionLoaders]);
  const [text, setText] = useState(() => load(KEY));
  const [shown, setShown] = useState<Shown | null>(null);
  const [focus, setFocus] = useState<{ start: number; end: number; key: number } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [note, setNote] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  /** Go to a place: with the rearrangement when it is on that chromosome (the view flies there), else the reference alone. */
  const goRegion = (chrom: string, start: number, end: number, b: GenomeBuild) => {
    const len = chromLength(b, chrom);
    if (!len) throw new Error(`Unknown chromosome ${chrom} on ${b}.`);
    if (start >= len) throw new Error(`${chrom} ends at ${len.toLocaleString('en-US')} on ${b}.`);
    const s = Math.max(0, start), e = Math.min(len, Math.max(end, s + 1));
    if (shown && (shown.kind === 'sv' ? shown.parsed.chrom : shown.chrom) === chrom && b === build) setFocus({ start: s, end: e, key: Date.now() });
    else { setShown({ kind: 'region', chrom, start: s, end: e }); setFocus(null); }
  };
  const draw = async (t = text, b: GenomeBuild = build) => {
    setError(null); setNote(null);
    const raw = t.trim();
    try {
      // coordinates: chr20:10,300,000-11,000,000, 20:10.3-11 Mb, chrX:31,100,000
      const loc = parseLocus(raw);
      if (loc) {
        const single = loc.start === loc.end;
        goRegion(loc.chrom, single ? loc.start - 1 - AROUND_POS : loc.start - 1, single ? loc.start + AROUND_POS : loc.end, b);
        save(KEY, t); return;
      }
      // a gene symbol
      if (/^[A-Za-z][A-Za-z0-9-]{0,19}(\.\d+)?$/.test(raw) && !/^NC_/i.test(raw) && !/:/.test(raw)) {
        if (!loaders.locateGene) throw new Error('Gene names need the UCSC API, not available here: type coordinates.');
        setBusy(true);
        const g = await loaders.locateGene(raw).finally(() => setBusy(false));
        const pad = Math.max(5000, (g.end - g.start) * 0.4);
        goRegion(g.chrom, g.start - pad, g.end + pad, b);
        setNote(`${g.name}: ${g.chrom}:${(g.start + 1).toLocaleString('en-US')}-${g.end.toLocaleString('en-US')} (${b})`);
        save(KEY, t); return;
      }
      // an HGVS or ISCN description: an NC_ accession's version decides the build
      const acc = raw.match(/NC_\d+\.\d+/i)?.[0];
      const accBuild = acc ? buildOfNc(acc) : null;
      if (accBuild && accBuild !== b) {
        if (!loadersFor) throw new Error(`${acc} is a ${accBuild} accession; this session reads ${b}.`);
        b = accBuild; setBuild(b); save(BUILD_KEY, b); setNote(`${acc} is a ${b} accession: switched to ${b}.`);
      }
      const parsed = parseSv(t, b, c => chromLength(b, c));
      setShown({ kind: 'sv', parsed, text: t }); setFocus(null); save(KEY, t);
    } catch (e) { setError((e as Error).message); }
  };
  const switchBuild = (b: GenomeBuild) => {
    if (b === build) return;
    setBuild(b); save(BUILD_KEY, b); setShown(null); setFocus(null);
    if (text.trim()) void draw(text, b);
  };

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onClose();
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [onClose]);
  const sv = shown?.kind === 'sv' ? shown.parsed : null;
  const chrom = sv ? sv.chrom : shown?.kind === 'region' ? shown.chrom : '';
  const bands = useCytobands(chrom || 'chr1', loaders.getCytobands);
  const colors = useMemo(() => (sv ? segmentColors(sv.rec) : new Map<number, string>()), [sv]);
  const full = useMemo<SvGeneLoaders>(() => ({ ...loaders, fallback: [] }), [loaders]);
  const rec = sv?.rec;
  return (
    <div className="fixed inset-0 z-[70] bg-white overflow-y-auto text-gray-900" style={{ fontFamily: FONT }}>
      <div className="w-full min-h-full">
        {/* a full-screen page: back to the Sashimi viewer from its bar (or Esc) */}
        <div className="sticky top-0 z-30 flex items-center gap-3 px-4 h-11 border-b border-gray-200 bg-white/95 shadow-sm">
          <button onClick={onClose} className="inline-flex items-center gap-1.5 h-8 px-3 rounded-md border border-gray-300 bg-white text-sm font-semibold text-gray-700 hover:bg-gray-50" title="Back to the Sashimi viewer (Esc)">← Sashimi viewer</button>
          <span className="text-gray-300">/</span>
          <span className="text-sm font-bold">SV browser</span>
          <span className="ml-auto text-[11px] text-gray-400">Esc: back to the viewer</span>
        </div>
        <div className="flex items-start justify-between px-5 pt-3 pb-2">
          <div>
            <div className="text-base font-bold flex items-center gap-2">SV browser
              {loadersFor
                ? <span className="inline-flex rounded-md border border-gray-300 overflow-hidden text-xs font-semibold" role="group" aria-label="Genome build">
                    {(['GRCh38', 'GRCh37'] as const).map(b => <button key={b} onClick={() => switchBuild(b)} className={`px-2 py-0.5 ${build === b ? 'bg-accent-600 text-white' : 'bg-white text-gray-600 hover:bg-gray-50'}`} title={b === 'GRCh38' ? 'GRCh38 / hg38 (default)' : 'GRCh37 / hg19: positions differ from GRCh38 (no liftover: type them in GRCh37)'}>{b}</button>)}
                  </span>
                : <span className="font-normal text-gray-500">· {build}</span>}
            </div>
            <div className="text-xs text-gray-500">a rearrangement typed as HGVS (NC_… or chrN:g., del, dup, inv, ins, delins, allele form […;…]) or as the ISCN sequence form (chrN:g.…::chrN:g.…); or a place to go to: coordinates (chr20:10.3-11 Mb) or a gene name</div>
          </div>
        </div>
        <div className="px-5 pb-3">
          <div className="flex gap-2">
            <input value={text} onChange={e => setText(e.target.value)} onKeyDown={e => { if (e.key === 'Enter') void draw(); }} spellCheck={false} autoFocus
              placeholder="NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]  ·  chr20:10,300,000-11,000,000  ·  JAG1"
              className="flex-1 font-mono text-[12.5px] h-9 px-3 rounded-md border border-gray-300 focus:outline-none focus:ring-2 focus:ring-accent-300" />
            <button onClick={() => void draw()} disabled={busy} className="h-9 px-4 rounded-md bg-accent-600 text-white text-sm font-semibold hover:bg-accent-700 disabled:opacity-60">{busy ? '…' : 'Go'}</button>
          </div>
          <div className="text-[11px] text-gray-500 mt-1">
            examples: {EXAMPLES.map((x, i) => <span key={x.label}>{i ? ' · ' : ''}<button className="text-accent-700 hover:underline" onClick={() => { setText(x.text); void draw(x.text); }}>{x.label}</button></span>)}
            {' '}· band-level ISCN has no positions and is not drawn; one chromosome at a time (no translocation yet)
          </div>
          {error && <div className="text-[12px] text-red-700 mt-2">{error}</div>}
          {note && <div className="text-[11.5px] text-accent-800 mt-1">{note}</div>}
          {sv && rec && (
            <div className="mt-2 grid grid-cols-[1fr_1fr] gap-4 text-xs">
              <div>
                <div className="text-[10.5px] font-semibold text-gray-700">Rebuilt chromosome</div>
                <div className="font-mono text-[11px] bg-gray-50 border border-gray-200 rounded px-2 py-1 mt-0.5">{rec.ders[0].name} = {rec.ders[0].pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ')}</div>
                {rec.lost.length > 0 && <div className="text-[11px] text-red-700 mt-0.5">Lost from the rearranged homologue: {rec.lost.map(s => s.letter).join(', ')}.</div>}
                {sv.notes.map((n, i) => <div key={i} className="text-[10.5px] text-gray-500">{n}</div>)}
                <div className="text-[10.5px] text-gray-500 mt-0.5">Drawn as heterozygous, 2 copies on the flanks; no read depth (open the cartoon on a sample for that).</div>
              </div>
              <SvNomenclatureBox rec={rec} build={build} bands={bands} />
            </div>
          )}
        </div>
        {shown && (
          <SvGenomeBrowser key={`${build}|${shown.kind === 'sv' ? shown.text : `region:${shown.chrom}`}`} rec={rec ?? null} chrom={chrom}
            window={sv ? sv.window : { start: (shown as { start: number }).start, end: (shown as { end: number }).end }}
            segColors={colors} flankCopies={2} zygosity="het" build={build} loaders={full} rebuild={sv?.rebuild} focus={focus} />
        )}
      </div>
    </div>
  );
}
