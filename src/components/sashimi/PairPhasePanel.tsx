/**
 * Phasing panel: a variant site against one or more positions typed (`chr:pos; chr:pos:A>G; …`), in cis or in trans
 * from the reads of one sample (pairPhase.ts). Each pair shows its two sites side by side, the verdict, the fragments
 * by allele and how they link the sites; several positions give a matrix of every pair. HTML, never exported.
 */
import React, { useCallback, useEffect, useRef, useState } from 'react';
import type { SashimiDataSource } from './datasource';
import type { AlignedRead, VariantSite } from './types';
import { READS_REGION_MAX_BP } from './datasource';
import { CHAIN_LONG_MIN_INDEL, CHAIN_READS_CAP, EPS, MAX_AGAINST, MIN_FRAGMENTS, PAIR_READS_CAP, chainPair, needsChain, parseTargets, phasePair, resolveTarget, targetWindow, type LinkKind, type PairCall, type PairFragment, type PairHaplotags, type PairSource, type SiteTagReads, type TagLane, type PairResult, type PairTarget } from '../../standalone/pairPhase';
import { HET_MAX, HET_MIN } from '../../standalone/phasing';
import { PHASE_HI, PHASE_LO } from './sitePhase';

export interface PairPhaseOrigin { sid: number; chrom: string; site: VariantSite }

interface Props {
  ds: SashimiDataSource;
  origin: PairPhaseOrigin;
  /** what was typed in the variant panel */
  text: string;
  samples: { id: number; name: string }[];
  /** the variant scan's sites of a sample on a chromosome, preferred to a call from the reads read here */
  knownSites: (sid: number, chrom: string) => VariantSite[];
  uniqueOnly: boolean;
  x: number; y: number;
  siteColor: (s: VariantSite) => string;
  siteShort: (s: VariantSite) => string;
  /** the hue of a haplotype of a phase set, as its reads and lanes are drawn */
  hpColor: (ps: number | null | undefined, hp: number) => string;
  /** the gene shown in the view, on the clicked site's chromosome: its exons (0-based half-open) placed on the map */
  gene?: { name: string; strand: number; exons: [number, number][] } | null;
  /** picking sites by clicking them in the variants track: on or off, and the last site picked there */
  picking: boolean;
  onPicking: (on: boolean) => void;
  picked: { n: number; chrom: string; site: VariantSite } | null;
  onClose: () => void;
}

interface SiteRow { label: string; target: PairTarget; site: VariantSite | null; note?: string; reads: AlignedRead[]; total: number; shown: number; error?: string }
type PairCell = { result: PairResult } | { reason: string };
interface Run { sid: number; sites: SiteRow[]; pairs: Map<string, PairCell>; errors: string[] }

export const PAIR_COLORS: Record<PairCall, string> = { cis: '#1d4ed8', trans: '#c2410c', conflict: '#b91c1c', few: '#9ca3af', hom: '#6b7280' };
const CALL_WORD: Record<PairCall, string> = { cis: 'in cis', trans: 'in trans', conflict: 'conflicting', few: 'undetermined', hom: 'homozygous' };
const SOURCE_WORD: Record<PairSource, string> = { reads: 'from the reads holding both sites', haplotags: "from the file's haplotags (HP / PS)", both: 'reads and haplotags agree', chain: 'from the reads phased through the sites between them' };
const REF_FILL = '#d1d5db';
const LINK_DASH: Record<LinkKind, string | undefined> = { read: undefined, spliced: '1.5 2', pair: '5 3', split: '1.5 2' };
const LINK_WORD: Record<LinkKind, string> = { read: 'one read', spliced: 'one spliced read', pair: 'two mates', split: 'a split read' };
const key = (i: number, j: number) => `${i}-${j}`;
const fmt = (n: number) => n.toLocaleString('en-US');

function originTarget(o: PairPhaseOrigin): PairTarget {
  const s = o.site;
  return { chrom: o.chrom, pos: s.pos, allele: s.kind === 'snv' ? { kind: 'snv', alt: s.alt } : { kind: s.kind, length: s.length }, text: `${o.chrom}:${fmt(s.pos + 1)}` };
}

export function PairPhasePanel({ ds, origin, text: text0, samples, knownSites, uniqueOnly, x, y, siteColor, siteShort, hpColor, gene, picking, onPicking, picked, onClose }: Props) {
  const [sid, setSid] = useState(origin.sid);
  const [text, setText] = useState(text0);
  const [run, setRun] = useState<Run | null>(null);
  const [busy, setBusy] = useState<string | false>(false);
  const [pick, setPick] = useState<[number, number]>([0, 1]);
  const abort = useRef<AbortController | null>(null);

  const go = useCallback(async (sampleId: number, typed: string) => {
    abort.current?.abort();
    const ac = new AbortController();
    abort.current = ac;
    setBusy('reading…');
    const { targets, errors } = parseTargets(typed, origin.chrom);
    const all = [originTarget(origin), ...targets];
    const sites: SiteRow[] = await Promise.all(all.map(async (t, i): Promise<SiteRow> => {
      const label = i ? (all.length > 2 ? `B${i}` : 'B') : 'A';
      try {
        const [s, e] = targetWindow(t);
        const r = await ds.getReads(sampleId, t.chrom, s, e, uniqueOnly, PAIR_READS_CAP, 'reads', 1, 0.01, { signal: ac.signal });
        const got = resolveTarget(t, r.reads, r.reference?.seq ?? null, r.reference?.start ?? 0, knownSites(sampleId, t.chrom));
        return { label, target: t, site: got.site, note: got.note, reads: r.reads, total: r.total, shown: r.shown };
      } catch (err: any) {
        return { label, target: t, site: null, reads: [], total: 0, shown: 0, error: err?.message ?? String(err) };
      }
    }));
    if (ac.signal.aborted) return;
    const pairs = new Map<string, PairCell>();
    for (let i = 0; i < sites.length; i++) for (let j = i + 1; j < sites.length; j++) {
      const a = sites[i], b = sites[j];
      if (!a.site || !b.site) pairs.set(key(i, j), { reason: `${!a.site ? a.label : b.label} not found in this sample` });
      else if (a.target.chrom !== b.target.chrom) pairs.set(key(i, j), { reason: 'on different chromosomes: no read joins them' });
      else if (a.site.pos === b.site.pos && a.site.kind === b.site.kind && a.site.alt === b.site.alt) pairs.set(key(i, j), { reason: 'the same variant' });
      else pairs.set(key(i, j), { result: phasePair(a.site, b.site, a.reads, b.reads, [a.label, b.label]) });
    }
    // pairs neither the reads holding both nor the haplotags decide: the reads of the whole stretch phased, once for
    // every such pair of the chromosome, linking them through the heterozygous sites between them
    const open = [...pairs].filter(([, c]) => 'result' in c && needsChain(c.result)).map(([k]) => k.split('-').map(Number) as [number, number]);
    if (open.length) {
      const chrom = sites[open[0][0]].target.chrom;
      const mine = open.filter(([i, j]) => sites[i].target.chrom === chrom && sites[j].target.chrom === chrom);
      const ps = mine.flatMap(([i, j]) => [sites[i].site!.pos, sites[j].site!.pos]);
      const lo = Math.max(0, Math.min(...ps) - 1000), hi = Math.max(...ps) + 1000;
      const note = (t: string) => { for (const [i, j] of mine) { const c = pairs.get(key(i, j)); if (c && 'result' in c) pairs.set(key(i, j), { result: { ...c.result, chainNote: t } }); } };
      if (hi - lo > READS_REGION_MAX_BP) note(`the stretch between them (${fmt(hi - lo)} bp) is too long to phase from the reads here (${fmt(READS_REGION_MAX_BP)} bp at most)`);
      else {
        setBusy('phasing the reads between the sites…');
        try {
          const r = await ds.getReads(sampleId, chrom, lo, hi, uniqueOnly, CHAIN_READS_CAP, 'collapsed', 3, 0.05, { haplotypes: 2, phaseSource: 'reads', longReadMinIndel: CHAIN_LONG_MIN_INDEL, signal: ac.signal });
          if (ac.signal.aborted) return;
          if (!r.phase) note('the reads between them could not be phased');
          else for (const [i, j] of mine) {
            const c = pairs.get(key(i, j));
            if (!c || !('result' in c)) continue;
            const ch = chainPair(r.phase, sites[i].site!, sites[j].site!, [sites[i].label, sites[j].label]);
            const sampledNote = r.total > r.shown ? ` (${fmt(r.shown)} of ${fmt(r.total)} reads, whole fragments)` : '';
            pairs.set(key(i, j), { result: ch.call
              ? { ...c.result, call: ch.call, source: 'chain', note: ch.note + sampledNote, chain: ch.chain }
              : { ...c.result, chainNote: ch.note + sampledNote } });
          }
        } catch (err: any) {
          if (ac.signal.aborted) return;
          note(`the reads between them could not be read: ${err?.message ?? err}`);
        }
      }
    }
    setRun({ sid: sampleId, sites, pairs, errors });
    setPick(p => (p[1] < sites.length ? p : [0, 1]));
    setBusy(false);
  }, [ds, origin, uniqueOnly, knownSites]);

  useEffect(() => { void go(sid, text); return () => abort.current?.abort(); }, [sid]); // eslint-disable-line react-hooks/exhaustive-deps

  // a site clicked in the variants track while picking: added to the list (once), the list phased again
  const seenPick = useRef(picked?.n ?? 0);
  useEffect(() => {
    if (!picked || picked.n === seenPick.current) return;
    seenPick.current = picked.n;
    const st = picked.site;
    const one = `${picked.chrom}:${st.pos + 1}:${st.kind === 'snv' ? st.alt : `${st.kind === 'ins' ? '+' : '-'}${st.length}`}`;
    const same = (t: { chrom: string; pos: number }) => t.chrom === picked.chrom && t.pos === st.pos;
    if (same({ chrom: origin.chrom, pos: origin.site.pos }) || parseTargets(text, origin.chrom).targets.some(same)) return;
    const next = text.trim() ? `${text.trim().replace(/;\s*$/, '')}; ${one}` : one;
    setText(next);
    void go(sid, next);
  }, [picked]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    if (!picking) return;
    const esc = (e: KeyboardEvent) => { if (e.key === 'Escape') onPicking(false); };
    window.addEventListener('keydown', esc);
    return () => window.removeEventListener('keydown', esc);
  }, [picking, onPicking]);

  const sites = run?.sites ?? [];
  const [pi, pj] = pick;
  const cell = run?.pairs.get(key(pi, pj));
  // where it opens, kept on screen; dragged by its title bar from there (it scrolls inside when taller than the room under its top)
  const [at, setAt] = useState(() => ({ left: Math.max(8, Math.min(x - 380, window.innerWidth - 772)), top: Math.max(8, Math.min(y + 12, window.innerHeight - 560)) }));
  const drag = useRef<{ dx: number; dy: number } | null>(null);
  const onDragStart = (e: React.PointerEvent<HTMLDivElement>) => {
    if ((e.target as HTMLElement).closest('select, button, input')) return;
    drag.current = { dx: e.clientX - at.left, dy: e.clientY - at.top };
    e.currentTarget.setPointerCapture(e.pointerId);
    e.preventDefault();
  };
  const onDragMove = (e: React.PointerEvent<HTMLDivElement>) => {
    const d = drag.current;
    if (!d) return;
    setAt({ left: Math.max(8 - 700, Math.min(e.clientX - d.dx, window.innerWidth - 60)), top: Math.max(0, Math.min(e.clientY - d.dy, window.innerHeight - 40)) });
  };
  const onDragEnd = () => { drag.current = null; };
  const top = at.top;
  const label = (r: SiteRow) => `${r.target.chrom}:${fmt((r.site?.pos ?? r.target.pos) + 1)}`;

  return (
    <div className="fixed z-50 w-[760px] max-w-[96vw] overflow-y-auto rounded-lg border border-gray-300 bg-white shadow-2xl text-xs" data-pair-panel
      style={{ left: at.left, top, maxHeight: window.innerHeight - top - 8 }}
      onMouseDown={e => e.stopPropagation()}>
      <div className="flex items-center justify-between gap-2 px-3 pt-2 pb-1 cursor-move select-none" title="Drag to move the panel"
        onPointerDown={onDragStart} onPointerMove={onDragMove} onPointerUp={onDragEnd} onPointerCancel={onDragEnd}>
        <div className="font-semibold text-gray-900">Phasing two variants</div>
        <div className="flex items-center gap-2">
          <select value={sid} onChange={e => setSid(Number(e.target.value))} className="border border-gray-300 rounded px-1 py-0.5 text-xs" title="The sample whose reads are read">
            {samples.map(s => <option key={s.id} value={s.id}>{s.name}</option>)}
          </select>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-700 text-base leading-none" title="Close">×</button>
        </div>
      </div>
      <form className="flex items-center gap-1.5 px-3 pt-1.5" onSubmit={e => { e.preventDefault(); void go(sid, text); }}>
        <span className="text-gray-500 shrink-0 font-mono">{origin.chrom}:{fmt(origin.site.pos + 1)} <span style={{ color: siteColor(origin.site) }}>{siteShort(origin.site)}</span> with</span>
        <input value={text} onChange={e => setText(e.target.value)} placeholder="chr:pos; chr:pos:A>G; …" spellCheck={false}
          className="min-w-0 flex-1 border border-gray-300 rounded px-1.5 py-0.5 font-mono text-xs outline-none focus:border-blue-400"
          title="Positions separated by ; — chr:pos (1-based), pos alone on this chromosome, an allele as chr:pos:A>G (VCF style, indels too), chr:pos:+2 / -3, or 1-55051215-G-GA" />
        <button type="button" onClick={() => onPicking(!picking)} title="Click sites in the variants track to add them to the list (click again or Esc to stop)"
          className={`rounded border px-2 py-0.5 shrink-0 ${picking ? 'border-blue-500 bg-blue-50 text-blue-700' : 'border-gray-300 hover:bg-gray-50'}`}>{picking ? 'picking…' : '⌖ Pick'}</button>
        <button type="submit" className="rounded border border-gray-300 px-2 py-0.5 hover:bg-gray-50 shrink-0" disabled={!!busy}>{busy || 'Phase'}</button>
      </form>
      {picking && <div className="px-3 pt-1 text-[10px] text-blue-700">Click a site in the variants track to add it (several in a row); ⌖ again or Esc to stop.</div>}
      {run && (run.errors.length > 0 || sites.some(s => s.note || s.error)) && (
        <div className="px-3 pt-1 text-[10px] leading-4">
          {run.errors.map((e, i) => <div key={`e${i}`} className="text-red-600">{e}</div>)}
          {sites.map((s, i) => (s.error || s.note) && <div key={`n${i}`} className={s.error || !s.site ? 'text-red-600' : 'text-gray-500'}>{s.label}: {s.error ?? s.note}</div>)}
        </div>
      )}
      {run && sites.length >= 2 && <PhaseMap sites={sites} pairs={run.pairs} siteColor={siteColor} siteShort={siteShort} hpColor={hpColor} gene={gene} />}
      {sites.length > 2 && (
        <Matrix sites={sites} pairs={run!.pairs} pick={pick} onPick={(i, j) => setPick([i, j])} siteColor={siteColor} siteShort={siteShort} />
      )}
      {run && sites.length >= 2 && cell && (
        <div className="px-3 pb-2 pt-2 mt-1.5 border-t border-gray-100">
          <div className="grid grid-cols-[1fr_auto_1fr] gap-3 items-start">
            <SiteCard row={sites[pi]} siteColor={siteColor} siteShort={siteShort} label={label(sites[pi])} />
            <Verdict cell={cell} a={sites[pi]} b={sites[pj]} />
            <SiteCard row={sites[pj]} siteColor={siteColor} siteShort={siteShort} label={label(sites[pj])} right />
          </div>
          {'result' in cell && <Evidence res={cell.result} a={sites[pi]} b={sites[pj]} siteColor={siteColor} hpColor={hpColor} />}
        </div>
      )}
      {!run && <div className="px-3 py-2 text-gray-400">reading the reads over each site…</div>}
      <div className="px-3 pb-2 text-[10px] text-gray-400 leading-4 border-t border-gray-100 pt-1">
        Fragments: the reads over each site joined by name (one read, a read spliced between them, two mates, a split read's parts). An allele under base quality 20, another allele, or two reads of a fragment disagreeing: unknown.
        A verdict needs {MIN_FRAGMENTS} informative fragments and at most {MAX_AGAINST * 100} % against; the LOD (log10 cis vs trans, {EPS * 100} % error per allele read) weighs them, it does not decide.
        A site under {HET_MIN * 100} % VAF counts only the fragments carrying it; one over {HET_MAX * 100} % is homozygous. Up to {fmt(PAIR_READS_CAP)} reads per site, whole fragments kept by name past it.
      </div>
    </div>
  );
}

/** Where each site sits: one of the two copies (0: A's, 1: the other), both (homozygous), or not placed. */
type Side = 0 | 1 | 'both' | null;

/**
 * The sites placed on the two copies of the chromosome: A on its haplotag lane when the tags place it (else on the first
 * copy), each other site from its verdict with a site already placed (cis: the same copy, trans: the other), walking
 * the decided pairs outwards from A; homozygous sites on both. A decided pair that contradicts the placement is kept.
 */
function placeSites(sites: SiteRow[], pairs: Map<string, PairCell>): { side: Side[]; via: (string | null)[]; ps: number | null; tagged: boolean; contradictions: [number, number][] } {
  const n = sites.length, side: Side[] = new Array(n).fill(null), via: (string | null)[] = new Array(n).fill(null);
  const res = (i: number, j: number) => { const c = pairs.get(key(Math.min(i, j), Math.max(i, j))); return c && 'result' in c ? c.result : null; };
  let ps: number | null = null, tagged = false;
  for (let i = 0; i < n; i++) if (sites[i].site && sites[i].site!.vaf > HET_MAX) { side[i] = 'both'; via[i] = 'homozygous'; }
  if (sites[0].site && side[0] == null) {
    side[0] = 0;
    // A's haplotag lane, from any of its pairs
    for (let j = 1; j < n && !tagged; j++) {
      const h = res(0, j)?.haplotags;
      if (h && (h.a.call === 'h1' || h.a.call === 'h2')) { side[0] = h.a.call === 'h1' ? 0 : 1; ps = h.a.ps; tagged = true; }
    }
    via[0] = tagged ? `on H${(side[0] as number) + 1} by its haplotags` : 'the reference copy of this drawing';
    const queue = [0];
    while (queue.length) {
      const i = queue.shift()!;
      for (let j = 0; j < n; j++) {
        if (side[j] != null || !sites[j].site || sites[j].target.chrom !== sites[0].target.chrom) continue;
        const r = res(i, j);
        if (!r || (r.call !== 'cis' && r.call !== 'trans')) continue;
        const si = side[i] as 0 | 1;
        side[j] = r.call === 'cis' ? si : ((1 - si) as 0 | 1);
        via[j] = `${CALL_WORD[r.call]} with ${sites[i].label} (${SOURCE_WORD[r.source]})`;
        queue.push(j);
      }
    }
  }
  const contradictions: [number, number][] = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const r = res(i, j), a = side[i], b = side[j];
    if (!r || (r.call !== 'cis' && r.call !== 'trans') || (a !== 0 && a !== 1) || (b !== 0 && b !== 1)) continue;
    if ((r.call === 'cis') !== (a === b)) contradictions.push([i, j]);
  }
  return { side, via, ps, tagged, contradictions };
}

/**
 * The two copies of the chromosome drawn as two lines, the sites at their positions (to scale) on the copy that carries
 * them, the view's gene and its exons above; a site not placed sits between the lines with a question mark.
 */
function PhaseMap({ sites, pairs, siteColor, siteShort, hpColor, gene }: { sites: SiteRow[]; pairs: Map<string, PairCell>; siteColor: (s: VariantSite) => string; siteShort: (s: VariantSite) => string; hpColor: (ps: number | null | undefined, hp: number) => string; gene?: Props['gene'] }) {
  const chrom = sites[0].target.chrom;
  const drawn = sites.map((r, i) => ({ r, i })).filter(({ r }) => r.site && r.target.chrom === chrom);
  if (drawn.length < 1 || !sites[0].site) return null;
  const { side, via, ps, tagged, contradictions } = placeSites(sites, pairs);
  const W = 712, X0 = 92, X1 = W - 14, GENE_Y = 16, L1 = 50, L2 = 92, AXIS = 122, H = 140;
  const pos = drawn.map(d => d.r.site!.pos);
  let lo = Math.min(...pos), hi = Math.max(...pos);
  const pad = Math.max(50, Math.round((hi - lo) * 0.06));
  lo = Math.max(0, lo - pad); hi += pad + 1;
  const xOf = (p: number) => X0 + ((p - lo) / (hi - lo)) * (X1 - X0);
  const lineColor = (h: 0 | 1) => (tagged ? hpColor(ps, h + 1) : '#9ca3af');
  const lineLabel = (h: 0 | 1) => (tagged ? `H${h + 1}` : h === 0 ? `copy with ${sites[0].label}` : 'other copy');
  // markers at least GAP apart: sites a few kb apart on a 100 kb map are spread, each joined to its true position on
  // the axis (and the gene) by a leader line
  const GAP = 26;
  const order = drawn.map(d => d.i).sort((p, q) => sites[p].site!.pos - sites[q].site!.pos);
  const xs = order.map(i => xOf(sites[i].site!.pos));
  for (let k = 1; k < xs.length; k++) xs[k] = Math.max(xs[k], xs[k - 1] + GAP);
  if (xs.length && xs[xs.length - 1] > X1) {
    xs[xs.length - 1] = X1;
    for (let k = xs.length - 2; k >= 0; k--) xs[k] = Math.min(xs[k], xs[k + 1] - GAP);
  }
  const xAt = new Map(order.map((i, k) => [i, xs[k]]));
  const ticks = [lo, Math.round((lo + hi) / 2), hi];
  const exons = gene ? gene.exons.filter(([s, e]) => e > lo && s < hi) : [];
  const geneSpan = gene && gene.exons.length ? [Math.min(...gene.exons.map(e => e[0])), Math.max(...gene.exons.map(e => e[1]))] : null;
  const showGene = !!geneSpan && geneSpan[1] > lo && geneSpan[0] < hi;
  const bad = new Set(contradictions.flat());
  // the heterozygous sites through which the reads of the stretch phased a pair with A, on their copy
  const relay = new Map<number, 0 | 1>();
  if (side[0] === 0 || side[0] === 1) for (let j = 1; j < sites.length; j++) {
    const c = pairs.get(key(0, j));
    if (c && 'result' in c && c.result.source === 'chain' && c.result.chain) for (const b of c.result.chain.between) relay.set(b.pos, (side[0] === 0 ? b.copy : 1 - b.copy) as 0 | 1);
  }
  const others = sites.filter(r => r.site && r.target.chrom !== chrom);
  return (
    <div className="px-3 pt-2">
      <svg width={W} height={H} className="block" style={{ maxWidth: '100%' }}>
        {showGene && (
          <g>
            <line x1={xOf(Math.max(lo, geneSpan![0]))} y1={GENE_Y} x2={xOf(Math.min(hi, geneSpan![1]))} y2={GENE_Y} stroke="#6b7280" strokeWidth={1} />
            {exons.map(([s, e], k) => <rect key={k} x={xOf(Math.max(lo, s))} y={GENE_Y - 5} width={Math.max(1.5, xOf(Math.min(hi, e)) - xOf(Math.max(lo, s)))} height={10} fill="#1f2937" />)}
            <text x={4} y={GENE_Y + 4} fontSize={10} fontWeight={700} fill="#374151">{gene!.name} <tspan fontWeight={400} fill="#9ca3af">{gene!.strand < 0 ? '←' : '→'}</tspan></text>
          </g>
        )}
        {([0, 1] as const).map(h => (
          <g key={h}>
            <line x1={X0} y1={h ? L2 : L1} x2={X1} y2={h ? L2 : L1} stroke={lineColor(h)} strokeWidth={4} strokeLinecap="round" />
            <text x={4} y={(h ? L2 : L1) + 4} fontSize={10} fontWeight={700} fill={tagged ? lineColor(h) : '#4b5563'}>{lineLabel(h)}</text>
            {tagged && h === 0 && <text x={4} y={L1 + 16} fontSize={8.5} fill="#9ca3af">PS {ps ?? 'none'}</text>}
          </g>
        ))}
        {[...relay].map(([p, c]) => <rect key={`r${p}`} x={xOf(p) - 0.75} y={(c ? L2 : L1) - 5} width={1.5} height={10} fill="#374151" opacity={0.6}><title>{`heterozygous site ${chrom}:${fmt(p + 1)}, a relay of the reads' phasing`}</title></rect>)}
        {drawn.map(({ r, i }) => {
          const st = r.site!, xt = xOf(st.pos), x = xAt.get(i)!, c = siteColor(st), sd = side[i];
          const tip = `${r.label} ${r.target.chrom}:${fmt(st.pos + 1)} ${siteShort(st)} · ${sd == null ? 'not placed: no decided pair links it to the others' : via[i]}${bad.has(i) ? ' · a decided pair contradicts this placement' : ''}`;
          const marks: React.ReactNode[] = [];
          const dot = (yy: number, k: string) => <circle key={k} cx={x} cy={yy} r={6} fill={c} stroke={bad.has(i) ? '#b91c1c' : '#fff'} strokeWidth={bad.has(i) ? 2 : 1.5} />;
          if (sd === 'both') { marks.push(dot(L1, 'a'), dot(L2, 'b')); }
          else if (sd === 0 || sd === 1) marks.push(dot(sd ? L2 : L1, 'a'));
          else marks.push(<g key="q"><circle cx={x} cy={(L1 + L2) / 2} r={7} fill="#fff" stroke={c} strokeWidth={1.5} strokeDasharray="2 2" /><text x={x} y={(L1 + L2) / 2 + 3.5} textAnchor="middle" fontSize={10} fontWeight={700} fill={c}>?</text></g>);
          // the label above the first copy, or below the second when the site sits there alone
          const below = sd === 1;
          const ly = below ? L2 + 17 : L1 - 10;
          return (
            <g key={i}>
              <title>{tip}</title>
              {showGene && <line x1={xt} y1={GENE_Y - 7} x2={xt} y2={GENE_Y + 7} stroke={c} strokeWidth={1.5} />}
              <line x1={x} y1={L1 - 6} x2={x} y2={L2 + 6} stroke={c} strokeWidth={0.8} strokeDasharray="2 3" opacity={0.7} />
              <line x1={x} y1={L2 + (below ? 22 : 6)} x2={xt} y2={AXIS} stroke={c} strokeWidth={0.8} opacity={0.7} />
              <circle cx={xt} cy={AXIS} r={2} fill={c} />
              {marks}
              <text x={x} y={ly} textAnchor="middle" fontSize={10} fontWeight={700} fill={c}>{r.label}</text>
            </g>
          );
        })}
        <line x1={X0} y1={AXIS} x2={X1} y2={AXIS} stroke="#d1d5db" />
        {ticks.map((t, k) => (
          <g key={k}>
            <line x1={xOf(t)} y1={AXIS} x2={xOf(t)} y2={AXIS + 4} stroke="#9ca3af" />
            <text x={xOf(t)} y={AXIS + 14} fontSize={9} fill="#6b7280" textAnchor={k === 0 ? 'start' : k === 2 ? 'end' : 'middle'}>{k === 0 ? `${chrom}:` : ''}{fmt(t + 1)}</text>
          </g>
        ))}
      </svg>
      <div className="text-[10px] text-gray-400 leading-4">
        The two copies of {chrom}{showGene ? `, ${gene!.name}'s exons above` : ''}; the axis and the gene to scale, sites close together spread apart, each joined to its position: each site on the copy carrying it ({tagged ? `H1 and H2 of phase set ${ps ?? '(no PS)'}, the file's haplotags` : `${sites[0].label}'s copy drawn first, the haplotypes being unnamed without haplotags`}), placed from the decided pairs; both copies when homozygous; ? when no decided pair places it.
        {contradictions.length > 0 && <span className="text-red-600"> Contradicting pairs: {contradictions.map(([i, j]) => `${sites[i].label}–${sites[j].label}`).join(', ')} (red outline).</span>}
        {relay.size > 0 && <span> Small dark ticks: the {relay.size} heterozygous sites the reads' phasing went through, on their copy.</span>}
        {others.length > 0 && <span> Not drawn (another chromosome): {others.map(r => r.label).join(', ')}.</span>}
      </div>
    </div>
  );
}

function Matrix({ sites, pairs, pick, onPick, siteColor, siteShort }: { sites: SiteRow[]; pairs: Map<string, PairCell>; pick: [number, number]; onPick: (i: number, j: number) => void; siteColor: (s: VariantSite) => string; siteShort: (s: VariantSite) => string }) {
  return (
    <div className="px-3 pt-2 overflow-x-auto">
      <table className="text-[10px] border-collapse">
        <tbody>
          <tr>
            <td />
            {sites.slice(1).map((s, j) => <td key={j} className="px-1 text-center font-semibold text-gray-700" title={s.target.text}>{s.label}</td>)}
          </tr>
          {sites.slice(0, -1).map((a, i) => (
            <tr key={i}>
              <td className="pr-2 font-semibold text-gray-700 whitespace-nowrap" title={a.target.text}>
                {a.label} <span className="font-normal text-gray-500 font-mono">{a.target.chrom}:{fmt((a.site?.pos ?? a.target.pos) + 1)}</span>
                {a.site && <span className="font-normal" style={{ color: siteColor(a.site) }}> {siteShort(a.site)}</span>}
              </td>
              {sites.slice(1).map((_, jj) => {
                const j = jj + 1;
                if (j <= i) return <td key={j} />;
                const c = pairs.get(key(i, j));
                const call = c && 'result' in c ? c.result.call : null;
                const on = pick[0] === i && pick[1] === j;
                return (
                  <td key={j} className="p-0.5">
                    <button onClick={() => onPick(i, j)} title={c ? ('result' in c ? `${a.label}–${sites[j].label}: ${CALL_WORD[c.result.call]}${c.result.call === 'cis' || c.result.call === 'trans' ? ` (${SOURCE_WORD[c.result.source]})` : ''} · ${c.result.note}` : c.reason) : ''}
                      className={`w-[84px] rounded px-1 py-0.5 text-white ${on ? 'ring-2 ring-offset-1 ring-gray-500' : ''}`} style={{ background: call ? PAIR_COLORS[call] : '#e5e7eb', color: call ? '#fff' : '#6b7280' }}>
                      {call ? CALL_WORD[call] : '—'}{c && 'result' in c && c.result.source === 'haplotags' ? ' · HP' : ''}
                    </button>
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function SiteCard({ row, label, right, siteColor, siteShort }: { row: SiteRow; label: string; right?: boolean; siteColor: (s: VariantSite) => string; siteShort: (s: VariantSite) => string }) {
  const s = row.site;
  return (
    <div className={`min-w-0 ${right ? 'text-right' : ''}`}>
      <div className="font-semibold text-gray-900"><span className="text-gray-400 mr-1">{row.label}</span><span className="font-mono">{label}</span></div>
      {s ? (
        <>
          <div className="font-semibold" style={{ color: siteColor(s) }}>{siteShort(s)}</div>
          <div className="text-gray-600">VAF {Math.round(s.vaf * 100)} % · {fmt(s.alt_count)} of {fmt(s.depth)} reads</div>
          <div className="text-[10px] text-gray-400">{s.vaf > HET_MAX ? 'homozygous' : s.vaf < HET_MIN ? `under ${HET_MIN * 100} %: mosaic or subclonal?` : 'heterozygous'}{row.total > row.shown ? ` · ${fmt(row.shown)} of ${fmt(row.total)} reads read (whole fragments)` : ''}</div>
        </>
      ) : <div className="text-red-600">{row.error ?? 'not found in this sample'}</div>}
    </div>
  );
}

function Verdict({ cell, a, b }: { cell: PairCell; a: SiteRow; b: SiteRow }) {
  const d = a.site && b.site && a.target.chrom === b.target.chrom ? Math.abs(b.site.pos - a.site.pos) : null;
  if (!('result' in cell)) return <div className="w-[230px] text-center text-gray-500 pt-1">{cell.reason}</div>;
  const r = cell.result;
  return (
    <div className="w-[230px] text-center">
      <div className="text-[15px] font-bold uppercase tracking-wide" style={{ color: PAIR_COLORS[r.call] }}>{CALL_WORD[r.call]}</div>
      {(r.call === 'cis' || r.call === 'trans') && <div className="text-[10px] font-semibold text-gray-500">{SOURCE_WORD[r.source]}</div>}
      <div className="text-gray-600 leading-4 mt-0.5">{r.note}</div>
      <div className="text-[10px] text-gray-400 mt-0.5">
        {r.cis + r.trans > 0 && <>LOD {Math.abs(r.lod) >= 99.95 ? '>99' : Math.abs(r.lod).toFixed(1)}{r.lod > 0 ? ' for cis' : r.lod < 0 ? ' for trans' : ''} · </>}
        {d != null && <>{fmt(d)} bp apart</>}
      </div>
    </div>
  );
}

function Evidence({ res, a, b, siteColor, hpColor }: { res: PairResult; a: SiteRow; b: SiteRow; siteColor: (s: VariantSite) => string; hpColor: (ps: number | null | undefined, hp: number) => string }) {
  const ca = siteColor(a.site!), cb = siteColor(b.site!);
  const counted = new Set(res.counted);
  const cellStyle = (i: number): React.CSSProperties => ({ color: counted.has(i) ? ([0, 3].includes(i) ? PAIR_COLORS.cis : PAIR_COLORS.trans) : '#9ca3af', fontWeight: counted.has(i) ? 700 : 400 });
  const links = (Object.keys(res.links) as LinkKind[]).filter(k => res.links[k]).map(k => `${fmt(res.links[k])} by ${LINK_WORD[k]}`);
  // the reads joined over both sites: their table and fragments; none (sites farther apart than the reads) leaves them out
  const joined = res.fragments.length + res.unknown > 0;
  return (
    <div className="mt-2">
      {joined && (
        <div className="flex flex-wrap items-start gap-x-6 gap-y-1">
          <table className="text-[11px] text-gray-700">
            <tbody>
              <tr><td /><td className="px-2 font-semibold" style={{ color: cb }}>{b.label} alt</td><td className="px-2 font-semibold text-gray-500">{b.label} ref</td></tr>
              <tr><td className="pr-2 font-semibold" style={{ color: ca }}>{a.label} alt</td><td className="px-2 text-center" style={cellStyle(0)}>{res.table[0]}</td><td className="px-2 text-center" style={cellStyle(1)}>{res.table[1]}</td></tr>
              <tr><td className="pr-2 font-semibold text-gray-500">{a.label} ref</td><td className="px-2 text-center" style={cellStyle(2)}>{res.table[2]}</td><td className="px-2 text-center" style={cellStyle(3)}>{res.table[3]}</td></tr>
            </tbody>
          </table>
          <div className="text-[10px] text-gray-500 leading-4 max-w-[430px]">
            <div><span style={{ color: PAIR_COLORS.cis }} className="font-semibold">cis</span> cells: alt–alt, ref–ref · <span style={{ color: PAIR_COLORS.trans }} className="font-semibold">trans</span> cells: alt–ref, ref–alt{res.counted.length < 4 ? ' · grey: not counted' : ''}</div>
            <div>Fragments over both sites: {links.length ? links.join(', ') : 'none'}{res.unknown ? `; ${res.unknown} more with an allele unknown` : ''}</div>
          </div>
        </div>
      )}
      {joined && <Fragments res={res} a={a} b={b} ca={ca} cb={cb} hpColor={hpColor} />}
      {(res.chain || res.chainNote) && (
        <div className="mt-2 text-[11px] leading-4">
          <span className="font-semibold text-gray-700">Phased through the reads between them</span>{' '}
          {res.chain
            ? <span className="text-gray-600">the block {fmt(res.chain.block.start + 1)}–{fmt(res.chain.block.end)}: {res.chain.block.sites.length} heterozygous sites ({res.chain.between.length} between the two), {fmt(res.chain.block.support[0])} fragments on its haplotype 1 and {fmt(res.chain.block.support[1])} on 2, {fmt(res.chain.block.ambiguous)} matching both, {fmt(res.chain.block.conflicting)} disagreeing at a site</span>
            : <span className="text-gray-500">{res.chainNote}</span>}
          {res.chain && <div className="text-[10px] text-gray-400">The page's read-based phasing of the stretch (the Collapsed reads' rule, WhatsHap's principle): heterozygous sites joined into a block when at least 2 fragments link them, at most 20 % disagreeing. The sites between the two are drawn on the map above.</div>}
        </div>
      )}
      {res.tagReads[0] && res.tagReads[1] && <TagReads res={res} a={a} b={b} ca={ca} cb={cb} hpColor={hpColor} />}
    </div>
  );
}

/**
 * The reads over each site by haplotag lane, one box per read in the allele's colour (grey: reference): H1 and H2 of the
 * site's phase set, then the reads of another set and the untagged ones. Lane by lane, A on the left and B on the right.
 */
function TagReads({ res, a, b, ca, cb, hpColor }: { res: PairResult; a: SiteRow; b: SiteRow; ca: string; cb: string; hpColor: (ps: number | null | undefined, hp: number) => string }) {
  const [ta, tb] = res.tagReads as [SiteTagReads, SiteTagReads];
  const h = res.haplotags;
  const W = 712, LABEL = 64, HALF = (W - LABEL) / 2, BOX = 4, STEP = 5, LINE = 9, AREA = HALF - 96, PER_LINE = Math.floor(AREA / STEP), MAX = PER_LINE * 6;
  const lanes: { key: TagLane; label: string; hp?: 1 | 2 }[] = [
    { key: 'h1', label: 'H1', hp: 1 }, { key: 'h2', label: 'H2', hp: 2 }, { key: 'other', label: 'other set' }, { key: 'none', label: 'untagged' },
  ];
  const shown = lanes.filter(l => ta.lanes[l.key].alt + ta.lanes[l.key].ref + tb.lanes[l.key].alt + tb.lanes[l.key].ref > 0);
  const lines = (t: SiteTagReads, k: TagLane) => Math.max(1, Math.ceil(Math.min(MAX, t.lanes[k].alt + t.lanes[k].ref) / PER_LINE));
  // where each site's allele sits: the lane the haplotags place it on, bold
  const on = (p: PairHaplotags['a'] | undefined): TagLane | null => (p?.call === 'h1' ? 'h1' : p?.call === 'h2' ? 'h2' : null);
  const onA = on(h?.a), onB = on(h?.b);
  let y = 24;
  const els: React.ReactNode[] = [];
  for (const l of shown) {
    const n = Math.max(lines(ta, l.key), lines(tb, l.key)), hgt = n * LINE;
    els.push(<line key={`g${l.key}`} x1={0} y1={y - 3} x2={W} y2={y - 3} stroke="#f3f4f6" />);
    els.push(<text key={`l${l.key}`} x={4} y={y + 8} fontSize={10} fontWeight={700} fill={l.hp ? '#374151' : '#9ca3af'}>{l.label}</text>);
    ([[ta, ca, LABEL, onA, a], [tb, cb, LABEL + HALF, onB, b]] as const).forEach(([t, color, x0, mine, row], side) => {
      const c = t.lanes[l.key], total = c.alt + c.ref;
      if (l.hp) els.push(<rect key={`h${l.key}${side}`} x={x0} y={y} width={6} height={hgt - 1} fill={hpColor(t.ps, l.hp)}><title>{`HP ${l.hp} of phase set ${t.ps ?? '(no PS)'}`}</title></rect>);
      const boxes: React.ReactNode[] = [];
      for (let i = 0; i < Math.min(total, MAX); i++) {
        const col = i % PER_LINE, ln = Math.floor(i / PER_LINE);
        boxes.push(<rect key={i} x={x0 + 10 + col * STEP} y={y + ln * LINE} width={BOX} height={LINE - 2} fill={i < c.alt ? color : REF_FILL} />);
      }
      els.push(<g key={`b${l.key}${side}`}><title>{`${row.label} · ${l.label}: ${c.alt} reads with the allele, ${c.ref} with the reference`}</title>{boxes}</g>);
      els.push(
        <text key={`n${l.key}${side}`} x={x0 + 10 + AREA + 6} y={y + 8} fontSize={10} fill={total ? '#374151' : '#9ca3af'} fontWeight={mine === l.key ? 700 : 400}>
          {total ? `${c.alt} / ${total} alt${total > MAX ? ` (+${total - MAX})` : ''}` : '–'}
        </text>,
      );
    });
    y += hgt + 6;
  }
  const head = (row: SiteRow, t: SiteTagReads, color: string, x: number) => (
    <text x={x + 10} y={12} fontSize={10} fontWeight={700} fill={color}>{row.label} <tspan fill="#6b7280" fontWeight={400}>{row.target.chrom}:{fmt(row.site!.pos + 1)} · PS {t.ps ?? 'none'}</tspan></text>
  );
  return (
    <div className="mt-2">
      <div className="text-[11px] text-gray-700">
        <span className="font-semibold">Reads by haplotag</span> <span className="text-gray-500">(the file's HP / PS){h ? ': ' : ''}</span>
        {h && <span style={{ color: h.call === 'cis' ? PAIR_COLORS.cis : h.call === 'trans' ? PAIR_COLORS.trans : undefined }} className={h.call === 'cis' || h.call === 'trans' ? 'font-semibold' : ''}>{h.call === 'cis' || h.call === 'trans' ? `${CALL_WORD[h.call]}, ` : ''}{h.note}</span>}
      </div>
      <div className="mt-1 max-h-[320px] overflow-y-auto border border-gray-100 rounded">
        <svg width={W} height={y} className="block" style={{ maxWidth: '100%' }}>
          {head(a, ta, ca, LABEL)}
          {head(b, tb, cb, LABEL + HALF)}
          <line x1={LABEL + HALF - 4} y1={2} x2={LABEL + HALF - 4} y2={y} stroke="#e5e7eb" />
          {els}
        </svg>
      </div>
      <div className="mt-0.5 text-[10px] text-gray-400">
        One box per read over the site, in the allele's colour (grey: the reference), every base quality as the variants track counts. A site is on a haplotype when at least {Math.round(PHASE_HI * 100)} % of its reads carry the allele and at most {Math.round(PHASE_LO * 100)} % of the other's; H1 and H2 are the phasing tool's, numbered afresh in each phase set.
      </div>
    </div>
  );
}

/** A band's fragments by haplotag, for its heading: " · H1 24, H2 1, untagged 3" */
function tagsOf(frags: PairFragment[]): string {
  const n = new Map<string, number>();
  for (const f of frags) { const k = f.hp != null ? `H${f.hp}` : 'untagged'; n.set(k, (n.get(k) ?? 0) + 1); }
  if (!n.size || (n.size === 1 && n.has('untagged'))) return '';
  return ' · ' + [...n].sort((p, q) => p[0].localeCompare(q[0])).map(([k, v]) => `${k} ${v}`).join(', ');
}

/** The fragments, one line each between a box at A (left) and at B (right), not to scale: cis cells on top, trans below. */
function Fragments({ res, a, b, ca, cb, hpColor }: { res: PairResult; a: SiteRow; b: SiteRow; ca: string; cb: string; hpColor: (ps: number | null | undefined, hp: number) => string }) {
  const W = 712, XA = 190, XB = W - 120, ROW = 5, PER_BAND = 80, HEAD = 13;
  const bands: { cell: number; title: string; frags: PairFragment[] }[] = [
    { cell: 0, title: `${a.label} alt · ${b.label} alt`, frags: [] }, { cell: 3, title: `${a.label} ref · ${b.label} ref`, frags: [] },
    { cell: 1, title: `${a.label} alt · ${b.label} ref`, frags: [] }, { cell: 2, title: `${a.label} ref · ${b.label} alt`, frags: [] },
  ];
  for (const f of res.fragments) bands.find(x => x.cell === (f.a ? 0 : 2) + (f.b ? 0 : 1))!.frags.push(f);
  const order: LinkKind[] = ['read', 'spliced', 'pair', 'split'];
  for (const band of bands) band.frags.sort((p, q) => order.indexOf(p.link) - order.indexOf(q.link) || (p.ps ?? 0) - (q.ps ?? 0) || (p.hp ?? 0) - (q.hp ?? 0));
  const shown = bands.filter(x => x.frags.length);
  if (!shown.length) return null;
  const counted = new Set(res.counted);
  let yy = 18;
  const els: React.ReactNode[] = [];
  for (const band of shown) {
    const isCis = band.cell === 0 || band.cell === 3, on = counted.has(band.cell);
    els.push(
      <text key={`h${band.cell}`} x={4} y={yy + 9} fontSize={10} fontWeight={700} fill={on ? (isCis ? PAIR_COLORS.cis : PAIR_COLORS.trans) : '#9ca3af'}>
        {band.title} — {band.frags.length} {on ? (isCis ? '(cis)' : '(trans)') : '(not counted)'}{tagsOf(band.frags)}
      </text>,
    );
    yy += HEAD;
    const rows = band.frags.slice(0, PER_BAND);
    rows.forEach((f, k) => {
      const y = yy + k * ROW + ROW / 2;
      const tip = `${f.name}\n${LINK_WORD[f.link]}: ${a.label} ${f.a ? 'alt' : 'ref'}, ${b.label} ${f.b ? 'alt' : 'ref'}${f.hp != null ? `\nHP ${f.hp}${f.ps != null ? `, PS ${f.ps}` : ''}` : ''}${f.conflict ? '\nits reads disagree at a site, the other site counted' : ''}`;
      els.push(
        <g key={`${band.cell}-${k}`} opacity={on ? 1 : 0.4}>
          <title>{tip}</title>
          <rect x={XA - 18} y={y - ROW / 2} width={XB - XA + 36} height={ROW} fill="transparent" />
          {f.hp != null && <rect x={XA - 24} y={y - 1.8} width={8} height={3.6} fill={hpColor(f.ps, f.hp)} />}
          <line x1={XA + 9} y1={y} x2={XB - 9} y2={y} stroke="#6b7280" strokeWidth={0.9} strokeDasharray={LINK_DASH[f.link]} />
          <rect x={XA - 9} y={y - 1.8} width={18} height={3.6} fill={f.a ? ca : REF_FILL} />
          <rect x={XB - 9} y={y - 1.8} width={18} height={3.6} fill={f.b ? cb : REF_FILL} />
        </g>,
      );
    });
    yy += rows.length * ROW;
    if (band.frags.length > rows.length) { els.push(<text key={`m${band.cell}`} x={XA} y={yy + 8} fontSize={9} fill="#6b7280">+ {band.frags.length - rows.length} more</text>); yy += 11; }
    yy += 6;
  }
  return (
    <div className="mt-2 max-h-[320px] overflow-y-auto border border-gray-100 rounded">
      <svg width={W} height={yy} className="block" style={{ maxWidth: '100%' }}>
        <text x={XA} y={11} textAnchor="middle" fontSize={10} fontWeight={700} fill={ca}>{a.label}</text>
        <text x={XB} y={11} textAnchor="middle" fontSize={10} fontWeight={700} fill={cb}>{b.label}</text>
        <g fontSize={8.5} fill="#6b7280">
          <line x1={XA + 40} y1={8} x2={XA + 62} y2={8} stroke="#6b7280" /><text x={XA + 66} y={11}>one read</text>
          <line x1={XA + 118} y1={8} x2={XA + 140} y2={8} stroke="#6b7280" strokeDasharray={LINK_DASH.pair} /><text x={XA + 144} y={11}>two mates</text>
          <line x1={XA + 200} y1={8} x2={XA + 222} y2={8} stroke="#6b7280" strokeDasharray={LINK_DASH.split} /><text x={XA + 226} y={11}>spliced or split read</text>
          {res.fragments.some(f => f.hp != null) && <text x={XA - 24} y={11} textAnchor="middle">HP</text>}
        </g>
        {els}
      </svg>
    </div>
  );
}
