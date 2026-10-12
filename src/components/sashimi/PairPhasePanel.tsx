/**
 * Phasing panel: a variant site against one or more positions typed (`chr:pos; chr:pos:A>G; …`), in cis or in trans
 * from the reads of one sample (pairPhase.ts). Each pair shows its two sites side by side, the verdict, the fragments
 * by allele and how they link the sites; several positions give a matrix of every pair. HTML, never exported.
 */
import React, { useCallback, useEffect, useRef, useState } from 'react';
import type { SashimiDataSource } from './datasource';
import type { AlignedRead, VariantSite } from './types';
import { EPS, MAX_AGAINST, MIN_FRAGMENTS, PAIR_READS_CAP, parseTargets, phasePair, resolveTarget, targetWindow, type LinkKind, type PairCall, type PairFragment, type PairResult, type PairTarget } from '../../standalone/pairPhase';
import { HET_MAX, HET_MIN } from '../../standalone/phasing';

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
  onClose: () => void;
}

interface SiteRow { label: string; target: PairTarget; site: VariantSite | null; note?: string; reads: AlignedRead[]; total: number; shown: number; error?: string }
type PairCell = { result: PairResult } | { reason: string };
interface Run { sid: number; sites: SiteRow[]; pairs: Map<string, PairCell>; errors: string[] }

export const PAIR_COLORS: Record<PairCall, string> = { cis: '#1d4ed8', trans: '#c2410c', conflict: '#b91c1c', few: '#9ca3af', hom: '#6b7280' };
const CALL_WORD: Record<PairCall, string> = { cis: 'in cis', trans: 'in trans', conflict: 'conflicting', few: 'undetermined', hom: 'homozygous' };
const REF_FILL = '#d1d5db';
const LINK_DASH: Record<LinkKind, string | undefined> = { read: undefined, spliced: '1.5 2', pair: '5 3', split: '1.5 2' };
const LINK_WORD: Record<LinkKind, string> = { read: 'one read', spliced: 'one spliced read', pair: 'two mates', split: 'a split read' };
const key = (i: number, j: number) => `${i}-${j}`;
const fmt = (n: number) => n.toLocaleString('en-US');

function originTarget(o: PairPhaseOrigin): PairTarget {
  const s = o.site;
  return { chrom: o.chrom, pos: s.pos, allele: s.kind === 'snv' ? { kind: 'snv', alt: s.alt } : { kind: s.kind, length: s.length }, text: `${o.chrom}:${fmt(s.pos + 1)}` };
}

export function PairPhasePanel({ ds, origin, text: text0, samples, knownSites, uniqueOnly, x, y, siteColor, siteShort, hpColor, onClose }: Props) {
  const [sid, setSid] = useState(origin.sid);
  const [text, setText] = useState(text0);
  const [run, setRun] = useState<Run | null>(null);
  const [busy, setBusy] = useState(false);
  const [pick, setPick] = useState<[number, number]>([0, 1]);
  const abort = useRef<AbortController | null>(null);

  const go = useCallback(async (sampleId: number, typed: string) => {
    abort.current?.abort();
    const ac = new AbortController();
    abort.current = ac;
    setBusy(true);
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
    setRun({ sid: sampleId, sites, pairs, errors });
    setPick(p => (p[1] < sites.length ? p : [0, 1]));
    setBusy(false);
  }, [ds, origin, uniqueOnly, knownSites]);

  useEffect(() => { void go(sid, text); return () => abort.current?.abort(); }, [sid]); // eslint-disable-line react-hooks/exhaustive-deps

  const sites = run?.sites ?? [];
  const [pi, pj] = pick;
  const cell = run?.pairs.get(key(pi, pj));
  // kept on screen: it scrolls inside when taller than the room under its top
  const top = Math.max(8, Math.min(y + 12, window.innerHeight - 560));
  const label = (r: SiteRow) => `${r.target.chrom}:${fmt((r.site?.pos ?? r.target.pos) + 1)}`;

  return (
    <div className="fixed z-50 w-[760px] max-w-[96vw] overflow-y-auto rounded-lg border border-gray-300 bg-white shadow-2xl text-xs" data-pair-panel
      style={{ left: Math.max(8, Math.min(x - 380, window.innerWidth - 772)), top, maxHeight: window.innerHeight - top - 8 }}
      onMouseDown={e => e.stopPropagation()}>
      <div className="flex items-center justify-between gap-2 px-3 pt-2">
        <div className="font-semibold text-gray-900">Phasing from the reads</div>
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
        <button type="submit" className="rounded border border-gray-300 px-2 py-0.5 hover:bg-gray-50" disabled={busy}>{busy ? 'reading…' : 'Phase'}</button>
      </form>
      {run && (run.errors.length > 0 || sites.some(s => s.note || s.error)) && (
        <div className="px-3 pt-1 text-[10px] leading-4">
          {run.errors.map((e, i) => <div key={`e${i}`} className="text-red-600">{e}</div>)}
          {sites.map((s, i) => (s.error || s.note) && <div key={`n${i}`} className={s.error || !s.site ? 'text-red-600' : 'text-gray-500'}>{s.label}: {s.error ?? s.note}</div>)}
        </div>
      )}
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
                    <button onClick={() => onPick(i, j)} title={c ? ('result' in c ? `${a.label}–${sites[j].label}: ${CALL_WORD[c.result.call]} · ${c.result.note}` : c.reason) : ''}
                      className={`w-[74px] rounded px-1 py-0.5 text-white ${on ? 'ring-2 ring-offset-1 ring-gray-500' : ''}`} style={{ background: call ? PAIR_COLORS[call] : '#e5e7eb', color: call ? '#fff' : '#6b7280' }}>
                      {call ? CALL_WORD[call] : '—'}
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
  const h = res.haplotags;
  return (
    <div className="mt-2">
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
          {h && (
            <div className="mt-0.5">
              <span className="font-semibold text-gray-700">Haplotags</span> (the file's HP / PS):{' '}
              <span style={{ color: h.call === 'cis' ? PAIR_COLORS.cis : h.call === 'trans' ? PAIR_COLORS.trans : undefined }} className={h.call === 'cis' || h.call === 'trans' ? 'font-semibold' : ''}>
                {h.call === 'cis' || h.call === 'trans' ? `${CALL_WORD[h.call]}, ` : ''}{h.note}
              </span>
              <div className="text-gray-400">{a.label}: {h.a.note.split(' · ')[0]} · {b.label}: {h.b.note.split(' · ')[0]}</div>
            </div>
          )}
        </div>
      </div>
      <Fragments res={res} a={a} b={b} ca={ca} cb={cb} hpColor={hpColor} />
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
