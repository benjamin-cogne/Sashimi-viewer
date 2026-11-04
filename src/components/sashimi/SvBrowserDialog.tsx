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
import { GeneSuggest } from './GeneSuggest';
import { Segmented } from './controls';
import { rememberGene } from '../../standalone/geneSuggest';
import { svLinkHref } from '../../standalone/link';
import { parseSmallGenomic, TRANSCRIPT_HGVS_RE, SMALL_KIND_WORD, SMALL_MAX, type SmallVariant } from './smallVariant';

const FONT = '"Inter Viewer", Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const EXAMPLES = [
  { label: 'deletion', text: 'NC_000020.11:g.10489766_10783340del' },
  { label: 'inversion with deletions', text: 'NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]' },
  { label: 'inverted insertion', text: 'NC_000020.11:g.10807815_10807816ins10489766_10783340inv' },
  { label: 'ISCN sequence form', text: 'chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::chr20:g.10,808,806_cen_qter' },
  { label: 'coordinates', text: 'chr20:10,300,000-11,000,000' },
  { label: 'DEL + coordinates', text: 'DEL 20:10489766-10783340' },
  { label: 'small variant', text: 'chr2:g.60461323del' },
  { label: 'transcript HGVS', text: 'NM_022893.4:c.1592del' },
  { label: 'a chromosome', text: 'chr20' },
  { label: 'a gene', text: 'JAG1' },
];
const KEY = 'sashimi.svBrowserText', BUILD_KEY = 'sashimi.svBuild';
const load = (k: string): string => { try { return localStorage.getItem(k) ?? ''; } catch { return ''; } };
const save = (k: string, v: string) => { try { localStorage.setItem(k, v); } catch { /* storage blocked */ } };
/** A single position opens this much around it. */
const AROUND_POS = 50_000;
/** A small variant opens with this much on each side (then on its gene, once the gene models are in). */
const AROUND_VAR = 10_000;
/** An interval opens with this share of its length on each side (and at least PAD_MIN), highlighted. */
const PAD_SHARE = 0.25, PAD_MIN = 2_000;
/** `DEL 20:399240-599326`, `dup chr20:399,240-599,326`, `INV: 20:…`: the rearrangement in the HGVS form the box reads. */
function typedSv(raw: string): string | null {
  const m = raw.match(/^(del(?:etion)?|dup(?:lication)?|inv(?:ersion)?)\s*[:\s]\s*(.+)$/i);
  if (!m) return null;
  const loc = parseLocus(m[2].trim());
  if (!loc || loc.end <= loc.start) return null;
  const kind = m[1].slice(0, 3).toLowerCase();
  return `${loc.chrom.startsWith('chr') ? loc.chrom : `chr${loc.chrom}`}:g.${loc.start}_${loc.end}${kind}`;
}

/** a small variant marked on the reference: as typed, and on the genome when typed on a transcript */
type ShownVariant = { kind: 'variant'; chrom: string; start: number; end: number; v: SmallVariant; typed: string; hgvsg: string | null; spdi: string | null; text: string };
type Shown = { kind: 'sv'; parsed: ParsedSv; text: string } | { kind: 'region'; chrom: string; start: number; end: number } | ShownVariant;
type Loaders = Omit<SvGeneLoaders, 'fallback'>;

export default function SvBrowserDialog({ build: sessionBuild, loaders: sessionLoaders, loadersFor, onClose, initial }: {
  build: string; loaders: Loaders;
  /** loaders for another build (the SV browser's own switch); without it, the session's build only */
  loadersFor?: (build: GenomeBuild) => Loaders;
  onClose: () => void;
  /** a search to run at once (a link: #sv=…), in its build when it gives one */
  initial?: { text: string; build?: GenomeBuild };
}) {
  const [build, setBuild] = useState<GenomeBuild>(() => (loadersFor ? (initial?.build ?? ((load(BUILD_KEY) as GenomeBuild) || 'GRCh38')) : (sessionBuild as GenomeBuild)));
  const loaders = useMemo(() => (loadersFor ? loadersFor(build) : sessionLoaders), [build, loadersFor, sessionLoaders]);
  const [text, setText] = useState(() => initial?.text ?? load(KEY));
  /** the interval asked for, highlighted in the view (coordinates typed) */
  const [mark, setMark] = useState<{ start: number; end: number; key: number } | null>(null);
  const [linkCopied, setLinkCopied] = useState(false);
  /** the address bar follows the search, so the page's address is a link to it */
  const setLink = (t: string | null, b: GenomeBuild) => { try { history.replaceState(history.state, '', svLinkHref(t, b)); } catch { /* a sandboxed frame */ } setLinkCopied(false); };
  const close = () => { setLink(null, build); onClose(); };
  const [shown, setShown] = useState<Shown | null>(null);
  const [focus, setFocus] = useState<{ start: number; end: number; key: number } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [note, setNote] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  /** Go to a place: with the rearrangement when it is on that chromosome (the view flies there), else the reference alone. */
  const goRegion = (chrom: string, start: number, end: number, b: GenomeBuild, highlight?: { start: number; end: number }) => {
    const len = chromLength(b, chrom);
    if (!len) throw new Error(`Unknown chromosome ${chrom} on ${b}.`);
    if (start >= len) throw new Error(`${chrom} ends at ${len.toLocaleString('en-US')} on ${b}.`);
    const s = Math.max(0, start), e = Math.min(len, Math.max(end, s + 1));
    setMark(highlight ? { start: Math.max(0, highlight.start), end: Math.min(len, highlight.end), key: Date.now() } : null);
    if (shown && shown.kind !== 'variant' && (shown.kind === 'sv' ? shown.parsed.chrom : shown.chrom) === chrom && b === build) setFocus({ start: s, end: e, key: Date.now() });
    else { setShown({ kind: 'region', chrom, start: s, end: e }); setFocus(null); }
  };
  const draw = async (t = text, b: GenomeBuild = build) => {
    setError(null); setNote(null);
    const raw = t.trim();
    try {
      // an NC_ accession's version decides the build (HGVS, ISCN, a small variant)
      const acc = raw.match(/NC_\d+\.\d+/i)?.[0];
      const accBuild = acc ? buildOfNc(acc) : null;
      if (accBuild && accBuild !== b) {
        if (!loadersFor) throw new Error(`${acc} is a ${accBuild} accession; this session reads ${b}.`);
        b = accBuild; setBuild(b); save(BUILD_KEY, b); setNote(`${acc} is a ${b} accession: switched to ${b}.`);
      }
      const showVariant = (v: SmallVariant, typed: string, hgvsg: string | null, spdi: string | null) => {
        const len = chromLength(b, v.chrom);
        if (!len) throw new Error(`Unknown chromosome ${v.chrom} on ${b}.`);
        if (v.end > len) throw new Error(`${v.chrom} ends at ${len.toLocaleString('en-US')} on ${b}.`);
        setShown({ kind: 'variant', chrom: v.chrom, start: Math.max(0, v.start - AROUND_VAR), end: Math.min(len, v.end + AROUND_VAR), v, typed, hgvsg, spdi, text: raw });
        setFocus(null); setMark(null); save(KEY, t); setLink(raw, b);
      };
      // a transcript description (NM_022893.4:c.1592del): moved to the genome by Ensembl, then marked (or rebuilt when large)
      if (TRANSCRIPT_HGVS_RE.test(raw.replace(/\s+/g, ''))) {
        const ld = loadersFor ? loadersFor(b) : loaders;
        if (!ld.recodeHgvs) throw new Error('Transcript descriptions need Ensembl, not available here: type the variant on the genome (chr2:g.…).');
        setBusy(true);
        const r = await ld.recodeHgvs(raw.replace(/\s+/g, '')).finally(() => setBusy(false));
        const small = parseSmallGenomic(r.hgvsg);
        const warn = r.warnings.length ? ` · Ensembl: ${r.warnings.join('; ')}` : '';
        if (small) { showVariant(small, raw, r.hgvsg, r.spdi); setNote(`${raw} = ${r.hgvsg} (${b}, Ensembl Variant Recoder)${warn}`); return; }
        const parsed = parseSv(r.hgvsg, b, c => chromLength(b, c));
        setShown({ kind: 'sv', parsed, text: r.hgvsg }); setFocus(null); setMark(null); save(KEY, t); setLink(raw, b);
        setNote(`${raw} = ${r.hgvsg} (${b}, Ensembl Variant Recoder): ${SMALL_MAX} bp or more, drawn as a rearrangement${warn}`);
        return;
      }
      // a small genomic variant (under SMALL_MAX bp): marked on the reference, no rebuilt chromosome
      const small = parseSmallGenomic(raw);
      if (small) { showVariant(small, raw, null, null); if (!accBuild) setNote(null); return; }
      // coordinates: chr20:10,300,000-11,000,000, 20:10.3-11 Mb, chrX:31,100,000
      // DEL / DUP / INV before coordinates: the rearrangement (read below as HGVS)
      const sv0 = typedSv(raw);
      if (sv0) {
        const parsed = parseSv(sv0, b, c => chromLength(b, c));
        setShown({ kind: 'sv', parsed, text: sv0 }); setFocus(null); setMark(null); save(KEY, t); setLink(raw, b);
        setNote(`${raw.split(/[\s:]/)[0].toUpperCase()} read as ${sv0}`);
        return;
      }
      // coordinates: the interval with some of the chromosome on each side, highlighted; a position, 50 kb around it
      const loc = parseLocus(raw);
      if (loc) {
        const single = loc.start === loc.end;
        const pad = single ? AROUND_POS : Math.max(PAD_MIN, Math.round((loc.end - loc.start + 1) * PAD_SHARE));
        goRegion(loc.chrom, loc.start - 1 - pad, loc.end + pad, b, { start: loc.start - 1, end: loc.end });
        save(KEY, t); setLink(raw, b); return;
      }
      // a whole chromosome: chr3, 3, X, chrM / MT
      const chr = raw.match(/^(?:chr)?(\d{1,2}|X|Y|M|MT)$/i);
      if (chr) {
        const name = `chr${chr[1].toUpperCase() === 'MT' ? 'M' : chr[1].toUpperCase()}`;
        const len = chromLength(b, name);
        if (!len) throw new Error(`Unknown chromosome ${name} on ${b}.`);
        goRegion(name, 0, len, b);
        save(KEY, t); setLink(raw, b); return;
      }
      // a gene symbol or an Ensembl gene id (looked up as in the Sashimi viewer's search)
      if (/^[A-Za-z][A-Za-z0-9-]{0,19}(\.\d+)?$/.test(raw) && !/^NC_/i.test(raw) && !/:/.test(raw)) {
        if (!loaders.locateGene) throw new Error('Gene names need the UCSC API, not available here: type coordinates.');
        setBusy(true);
        const g = await loaders.locateGene(raw).finally(() => setBusy(false));
        const pad = Math.max(5000, (g.end - g.start) * 0.4);
        goRegion(g.chrom, g.start - pad, g.end + pad, b);
        setNote(`${g.name}${g.name.toUpperCase() !== raw.toUpperCase() ? ` (${raw})` : ''}: ${g.chrom}:${(g.start + 1).toLocaleString('en-US')}-${g.end.toLocaleString('en-US')} (${b})`);
        rememberGene(g.name); setLink(raw, b);
        save(KEY, t); return;
      }
      // an HGVS or ISCN description (its NC_ accession's build already picked above)
      const parsed = parseSv(t, b, c => chromLength(b, c));
      setShown({ kind: 'sv', parsed, text: t }); setFocus(null); setMark(null); save(KEY, t); setLink(raw, b);
    } catch (e) { setError((e as Error).message); }
  };
  const switchBuild = (b: GenomeBuild) => {
    if (b === build) return;
    setBuild(b); save(BUILD_KEY, b); setShown(null); setFocus(null);
    if (text.trim()) void draw(text, b);
  };

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') close();
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [onClose, build]); // eslint-disable-line react-hooks/exhaustive-deps
  // a link's search runs once, at the start
  useEffect(() => { if (initial?.text.trim()) void draw(initial.text, build); }, []); // eslint-disable-line react-hooks/exhaustive-deps
  const sv = shown?.kind === 'sv' ? shown.parsed : null;
  const chrom = sv ? sv.chrom : shown?.kind === 'region' || shown?.kind === 'variant' ? shown.chrom : '';
  const variant = shown?.kind === 'variant' ? shown : null;
  const bands = useCytobands(chrom || 'chr1', loaders.getCytobands);
  const colors = useMemo(() => (sv ? segmentColors(sv.rec) : new Map<number, string>()), [sv]);
  const full = useMemo<SvGeneLoaders>(() => ({ ...loaders, fallback: [] }), [loaders]);
  const rec = sv?.rec;
  return (
    <div className="fixed inset-0 z-[70] bg-page overflow-y-auto text-slate-900" style={{ fontFamily: FONT }}>
      <div className="w-full min-h-full pb-4">
        {/* a full-screen page: back to the Sashimi viewer from its bar (or Esc) */}
        {/* the Sashimi viewer's app bar */}
        <div className="sticky top-0 z-30 flex items-center gap-3 px-4 h-12 bg-bar border-b border-bar-line">
          <button onClick={close} className="inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="Back to the Sashimi viewer (Esc)">← Sashimi viewer</button>
          <span className="text-bar-muted">/</span>
          <span className="text-sm font-semibold text-bar-text">SV browser</span>
          {shown && <button onClick={() => { void navigator.clipboard?.writeText(window.location.href).then(() => setLinkCopied(true), () => {}); }}
            className="ml-auto inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300"
            title={'A link to this search: the page opens on the SV browser with it (#sv=… after the page\'s address; add it to the viewer\'s address to share it)'}>{linkCopied ? '✓ link copied' : '🔗 copy link'}</button>}
          <span className={`${shown ? '' : 'ml-auto '}text-[11px] text-bar-muted`}>Esc: back to the viewer</span>
        </div>
        <div className="mx-4 mt-3 bg-white border border-slate-200 rounded-xl shadow-sm">
        <div className="flex items-start justify-between px-5 pt-3 pb-2">
          <div>
            <div className="text-base font-bold flex items-center gap-2">SV browser
              {loadersFor
                ? <Segmented size="sm" label="Genome build" title="The genome build of the positions" value={build} onChange={switchBuild}
                    options={[{ value: 'GRCh38', label: 'GRCh38', hint: 'GRCh38 / hg38 (default)' }, { value: 'GRCh37', label: 'GRCh37', hint: 'GRCh37 / hg19: positions differ from GRCh38 (no liftover: type them in GRCh37)' }]} />
                : <span className="font-normal text-slate-500">· {build}</span>}
            </div>
            <div className="text-xs text-slate-500">a rearrangement typed as HGVS (NC_… or chrN:g., del, dup, inv, ins, delins, allele form […;…]) or as the ISCN sequence form (chrN:g.…::chrN:g.…); or DEL / DUP / INV before coordinates (DEL 20:399240-599326); a small variant (under {SMALL_MAX} bp: chr2:g.60461323del, or on a transcript: NM_022893.4:c.1592del, placed by Ensembl), marked on the reference; or a place to go to: coordinates (chr20:10.3-11 Mb, opened with some flank and highlighted), a whole chromosome (chr3 or 3) or a gene name (suggested as you type, as in the Sashimi viewer)</div>
          </div>
        </div>
        <div className="px-5 pb-3">
          <form className="flex gap-2" onSubmit={e => { e.preventDefault(); void draw(); }}>
            {/* the Sashimi viewer's search box: gene names suggested from the third letter (HGNC, the recent genes first) */}
            <GeneSuggest value={text} onChange={setText} width={420} autoFocus ariaLabel="Rearrangement, place or gene"
              placeholder="NC_000020.11:g.[10489766_10783340del;10783341_10807815inv;10807816_10808805del]  ·  chr20:10,300,000-11,000,000  ·  chr20  ·  JAG1"
              title={'A rearrangement (HGVS or ISCN sequence form), coordinates (chr20:10.3-11 Mb), a chromosome (chr3 or 3), or a gene symbol or ENSG id.\nGene names are suggested from the third letter (HGNC), with their previous symbols and aliases.'}
              className="flex-1 min-w-0 font-mono text-[12.5px] h-[34px] px-3 rounded-[6px] border border-slate-200 bg-white text-slate-800 placeholder:text-slate-400 focus:outline-none focus:border-accent-300 focus:ring-2 focus:ring-accent-100" />
            <button type="submit" disabled={busy} className="h-[34px] px-4 rounded-[6px] bg-accent-600 text-white text-[12.5px] font-semibold hover:bg-accent-700 disabled:opacity-60">{busy ? '…' : 'Go'}</button>
          </form>
          <div className="text-[11px] text-slate-500 mt-1">
            examples: {EXAMPLES.map((x, i) => <span key={x.label}>{i ? ' · ' : ''}<button className="text-accent-700 hover:underline" onClick={() => { setText(x.text); void draw(x.text); }}>{x.label}</button></span>)}
            {' '}· band-level ISCN has no positions and is not drawn; one chromosome at a time (no translocation yet)
          </div>
          {error && <div className="text-[12px] text-red-700 mt-2">{error}</div>}
          {note && <div className="text-[11.5px] text-accent-800 mt-1">{note}</div>}
          {sv && rec && (
            <div className="mt-2 grid grid-cols-[1fr_1fr] gap-4 text-xs">
              <div>
                <div className="text-[10.5px] font-semibold text-slate-700">Rebuilt chromosome</div>
                <div className="font-mono text-[11px] bg-slate-50 border border-slate-200 rounded px-2 py-1 mt-0.5">{rec.ders[0].name} = {rec.ders[0].pieces.map(p => `${p.seg.letter}${p.reversed ? '′' : ''}`).join(' · ')}</div>
                {rec.lost.length > 0 && <div className="text-[11px] text-red-700 mt-0.5">Lost from the rearranged homologue: {rec.lost.map(s => s.letter).join(', ')}.</div>}
                {sv.notes.map((n, i) => <div key={i} className="text-[10.5px] text-slate-500">{n}</div>)}
                <div className="text-[10.5px] text-slate-500 mt-0.5">Drawn as heterozygous, 2 copies on the flanks; no read depth (open the cartoon on a sample for that).</div>
              </div>
              <SvNomenclatureBox rec={rec} build={build} bands={bands} />
            </div>
          )}
        </div>
        </div>
        {shown && (
          <div className="mx-4 mt-3 bg-white border border-slate-200 rounded-xl shadow-sm pb-2 overflow-hidden [&>div:first-child]:border-t-0 [&>div:first-child]:mt-0">
          <SvGenomeBrowser key={`${build}|${shown.kind === 'sv' ? shown.text : shown.kind === 'variant' ? `variant:${shown.text}` : `region:${shown.chrom}`}`} rec={rec ?? null} chrom={chrom}
            window={sv ? sv.window : { start: (shown as { start: number }).start, end: (shown as { end: number }).end }}
            segColors={colors} flankCopies={2} zygosity="het" build={build} loaders={full} rebuild={sv?.rebuild} focus={focus} highlight={mark}
            variant={variant ? {
              start: variant.v.start, end: variant.v.end, key: variant.text,
              label: variant.hgvsg ? variant.typed.replace(/\(([^)]*)\)/, '') : variant.v.hgvs,
              lines: [
                `${SMALL_KIND_WORD[variant.v.kind]}${variant.v.kind === 'snv' ? '' : ` of ${variant.v.size} bp`}`,
                ...(variant.hgvsg ? [`on the genome: ${variant.hgvsg} (${build}, Ensembl Variant Recoder)`, ...(variant.spdi ? [`SPDI: ${variant.spdi}`] : [])] : []),
                `${variant.v.chrom}:${(variant.v.start + 1).toLocaleString('en-US')}${variant.v.end - variant.v.start > 1 ? `-${variant.v.end.toLocaleString('en-US')}` : ''} (${build})`,
              ],
            } : undefined} />
          </div>
        )}
      </div>
    </div>
  );
}
