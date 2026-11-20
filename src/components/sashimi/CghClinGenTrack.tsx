/**
 * The "ClinGen CNV regions" track (SvGenomeBrowser's optional `extraTrack`), the same in every browser: the SV browser,
 * the SV cartoon's and the array-CGH page's (useClinGenTrack):
 * HI 3 / TS 3 dosage regions, the recurrent CNV list's regions and its breakpoint / LCR segments, from embedded data.
 */
import { useId, useMemo, useState } from 'react';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { SvExtraTrack } from './SvGenomeBrowser';
import { clinGenOn, recurrentOf, shortName, PIEV_RO, LCR_PAD, type ClinGenDosage, type ClinGenRecurrent } from './cghRecurrent';
import { CLINGEN_READ, CLINGEN_DOSAGE_EVALUATED } from './cghClinGenData';
import { ACHRO_COLOR } from './achropuce';

const COL = { hi: '#dc2626', ts: '#2563eb', rec: '#ea580c', recText: '#9a3412', lcr: '#94a3b8', text: '#1f2937', faint: '#8796a8', piev: ACHRO_COLOR.PIEV };
const ROW = 24, LCR_ROW = 19;
// the tag names the type when only one is PIEV (22q11.2 A-D: the duplication, not the DGS deletion)
const tagText = (t: ('del' | 'dup')[]) => (t.length === 2 ? 'PIEV' : `PIEV ${t[0]}`);
const tagW = (t: ('del' | 'dup')[]) => (t.length ? (t.length === 2 ? 27 : 44) : 0);
// rough widths (no canvas measure here): enough to keep labels of one lane apart
const tw = (s: string, px: number) => s.length * px * 0.56;
const fmt = (n: number) => n.toLocaleString('en-US');
const size = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} Mb` : `${Math.round(n / 1e3)} kb`);
const pievWord = (r: { piev: ('del' | 'dup')[]; pievName?: string }) => (r.piev.length ? `\nAChro-Puce: ${r.pievName}, PIEV (${r.piev.join(' and ')}): susceptibility, reduced penetrance` : '');

type Item = { r: ClinGenDosage | ClinGenRecurrent; a: number; b: number; lane: number; label: string };
/** greedy lanes: an item goes on the first lane whose last label ended left of it */
function lanes(rs: (ClinGenDosage | ClinGenRecurrent)[], x: (p: number) => number, x0: number, labelOf: (r: ClinGenDosage | ClinGenRecurrent) => string, px: number): Item[] {
  const ends: number[] = [];
  return rs.slice().sort((p, q) => p.start - q.start || q.end - p.end).map(r => {
    const a = x(r.start), b = Math.max(x(r.end), a + 2), left = Math.max(a, x0), label = labelOf(r);
    const right = Math.max(b, left + 2 + tagW(r.piev) + tw(label, px));
    let lane = ends.findIndex(e => e + 8 < left);
    if (lane < 0) { lane = ends.length; ends.push(right); } else ends[lane] = right;
    return { r, a, b, lane, label };
  });
}

function Tag({ x, y, t }: { x: number; y: number; t: ('del' | 'dup')[] }) {
  const w = tagW(t) - 3;
  return <g><rect x={x} y={y - 8} width={w} height={10} rx={2} fill="#fff" stroke={COL.piev} strokeWidth={0.8} /><text x={x + w / 2} y={y} fontSize={7.5} fontWeight={700} fill={COL.piev} textAnchor="middle">{tagText(t)}</text></g>;
}

const FOLD_ROW = 16;
/** the fold switch, on its own row under the dosage regions (inside the plot: the browser clips its tracks left of x0) */
function Fold({ open, n, onToggle, x, y }: { open: boolean; n: number; onToggle: () => void; x: number; y: number }) {
  // a box under the text: a click between its glyphs still hits
  return <g style={{ cursor: 'pointer' }} onClick={onToggle} data-cgh-fold>
    <title>{open ? 'Hide the recurrent CNV regions and breakpoint / LCR segments' : 'Show the recurrent CNV regions (list v2.0) and their breakpoint / LCR segments'}</title>
    <rect x={x} y={y} width={150} height={FOLD_ROW - 2} fill="transparent" />
    <text x={x + 2} y={y + 10} fontSize={9.5} fill="#0f766e">{open ? '▾ hide' : '▸ show'} recurrent CNVs / LCRs ({n})</text></g>;
}

function Track({ dos, rec, lcr, nd, nr, x0, fold }: { dos: Item[]; rec: Item[]; lcr: Item[]; nd: number; nr: number; x0: number; fold: { open: boolean; n: number; onToggle: () => void } | null }) {
  const hatch = `cgh-hatch${useId().replace(/:/g, '')}`;
  // no dosage region in view: a line saying so takes its place
  const yFold = nd ? nd * ROW : 16, yRec = yFold + (fold ? FOLD_ROW : 0), yLcr = yRec + nr * ROW;
  const sw = fold && <Fold {...fold} x={x0 + 2} y={yFold} />;
  const lx = (it: Item) => Math.max(it.a, x0) + 2;
  return (
    <g fontSize={9.5}>
      {!dos.length && <text x={x0 + 4} y={10} fontSize={9.5} fill={COL.faint}>no ClinGen HI 3 / TS 3 region in view</text>}
      {sw}
      <defs><pattern id={hatch} width={5} height={5} patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width={5} height={5} fill="#fff7ed" /><line x1={0} y1={0} x2={0} y2={5} stroke={COL.rec} strokeWidth={1.6} /></pattern></defs>
      {dos.map(it => { const r = it.r as ClinGenDosage, w = Math.max(it.b - it.a, 2), both = r.hi === 3 && r.ts === 3, t = lx(it) + tagW(r.piev);
        return <g key={r.id} transform={`translate(0 ${it.lane * ROW})`}><a href={r.url} target="_blank" rel="noreferrer">
          <title>{`${r.name}\n${r.id} · ${r.band} · ${size(r.end - r.start)}\nHI ${r.hi ?? '–'} · TS ${r.ts ?? '–'} (ClinGen dosage, last evaluated ${r.evaluated})\n${r.chrom}:${fmt(r.start + 1)}-${fmt(r.end)}${pievWord(r)}\nClick: the ClinGen page`}</title>
          <rect x={it.a} y={0} width={w} height={22} fill="transparent" />
          {both ? <><rect x={it.a} y={1} width={w} height={5} fill={COL.hi} opacity={0.85} /><rect x={it.a} y={6} width={w} height={5} fill={COL.ts} opacity={0.85} /></>
            : <rect x={it.a} y={1} width={w} height={10} rx={1.5} fill={r.hi === 3 ? COL.hi : COL.ts} opacity={0.85} />}
          {r.piev.length > 0 && <Tag x={lx(it)} y={21} t={r.piev} />}
          <text x={t} y={21} fill={COL.text}><tspan fontWeight={600}>{it.label}</tspan>
            {r.hi === 3 && <tspan dx={5} fontWeight={700} fill={COL.hi}>HI 3</tspan>}{r.ts === 3 && <tspan dx={5} fontWeight={700} fill={COL.ts}>TS 3</tspan>}</text>
        </a></g>; })}
      {rec.map(it => { const r = it.r as ClinGenRecurrent, y = yRec + it.lane * ROW;
        return <g key={`r${r.start}-${r.end}-${r.name}`}>
          <title>{`${r.name}\nClinGen recurrent CNV list v2.0 · ${size(r.end - r.start)}\n${r.chrom}:${fmt(r.start + 1)}-${fmt(r.end)}${pievWord(r)}`}</title>
          <rect x={it.a} y={y + 2} width={Math.max(it.b - it.a, 2)} height={8} rx={1} fill={`url(#${hatch})`} stroke={COL.rec} strokeWidth={1} />
          {r.piev.length > 0 && <Tag x={lx(it)} y={y + 21} t={r.piev} />}
          <text x={lx(it) + tagW(r.piev)} y={y + 21} fill={COL.recText}>{it.label}</text>
        </g>; })}
      {lcr.map(it => { const r = it.r, y = yLcr + it.lane * LCR_ROW;
        return <g key={`l${r.start}-${r.name}`}>
          <title>{`${r.name}: breakpoint / LCR segment (ClinGen recurrent CNV list v2.0) · ${size(r.end - r.start)}\n${r.chrom}:${fmt(r.start + 1)}-${fmt(r.end)}`}</title>
          <rect x={it.a} y={y + 2} width={Math.max(it.b - it.a, 2)} height={5} fill={COL.lcr} />
          <text x={lx(it)} y={y + 17} fontSize={9} fill={COL.faint}>{it.label}</text>
        </g>; })}
    </g>
  );
}

const sw = (style: React.CSSProperties) => <span className="inline-block w-4 h-2.5 rounded-sm" style={style} />;
/** the track for a build, folded (HI 3 / TS 3 regions only) or open (with the recurrent list and its LCRs): memoise
 * per build and fold (a new object redraws the browser) */
/** The track for a browser on `build`, folded at first: one per browser, so the SV browser, the SV cartoon and the array-CGH page draw the same. */
export function useClinGenTrack(build: GenomeBuild): SvExtraTrack {
  const [open, setOpen] = useState(false);
  return useMemo(() => clinGenTrack(build, open, () => setOpen(o => !o)), [build, open]);
}
export function clinGenTrack(build: GenomeBuild, open: boolean, onToggle: () => void): SvExtraTrack {
  return {
    id: 'cghClinGen', label: 'ClinGen CNV regions', side: ['ClinGen CNV', 'regions'],
    title: `ClinGen dosage regions scored HI 3 (red) or TS 3 (blue), and the ClinGen recurrent CNV list v2.0 (hatched orange) with its breakpoint / LCR segments (grey): embedded, read ${CLINGEN_READ}; PIEV: the AChro-Puce network classes the region a susceptibility CNV`,
    legend: <span className="inline-flex flex-wrap items-center gap-1">ClinGen (embedded, read {CLINGEN_READ}; dosage last evaluated {CLINGEN_DOSAGE_EVALUATED}): {sw({ background: COL.hi })} HI 3 {sw({ background: COL.ts })} TS 3 region (click: its page)
      {sw({ background: `repeating-linear-gradient(45deg, ${COL.rec} 0 1.5px, #fff7ed 1.5px 4px)`, border: `1px solid ${COL.rec}` })} recurrent CNV region (list v2.0)
      <span className="inline-block w-4 h-1 bg-slate-400" /> breakpoint / LCR segment
      <span className="rounded-sm border px-0.5 text-[8.5px] font-bold" style={{ color: COL.piev, borderColor: COL.piev }}>PIEV</span> AChro-Puce PIEV (for del / dup when one only) (≥ {Math.round(PIEV_RO * 100)} % reciprocal, same type): reduced penetrance</span>,
    layout: ({ chrom, s, e, x, x0 }) => {
      const all = clinGenOn(build, chrom).filter(r => r.end > s && r.start < e);
      const dos = lanes(all.filter(r => r.kind === 'dosage'), x, x0, r => shortName(r.name), 9.5);
      const extra = all.filter(r => r.kind !== 'dosage').length;
      const rec = open ? lanes(all.filter(r => r.kind === 'recurrent'), x, x0, r => shortName(r.name), 9.5) : [];
      const lcr = open ? lanes(all.filter(r => r.kind === 'lcr'), x, x0, r => r.name, 9) : [];
      const fold = extra ? { open, n: extra, onToggle } : null;
      const n = (its: Item[]) => (its.length ? Math.max(...its.map(i => i.lane)) + 1 : 0);
      const nd = n(dos), nr = n(rec), nl = n(lcr);
      return { height: (nd ? nd * ROW : 16) + (extra ? FOLD_ROW : 0) + nr * ROW + nl * LCR_ROW, node: <Track dos={dos} rec={rec} lcr={lcr} nd={nd} nr={nr} x0={x0} fold={fold} /> };
    },
  };
}

const pct = (f: number) => `${Math.round(f * 100)} %`;
/** a call's detail: the ClinGen recurrent regions it matches and whether its ends lie in the list's breakpoint / LCR segments */
export function CghRecurrentDetail({ chrom, start, end, kind, build }: { chrom: string; start: number; end: number; kind: 'loss' | 'gain'; build: GenomeBuild }) {
  const type = kind === 'loss' ? 'del' : 'dup';
  const { matches, ends, near } = recurrentOf(build, chrom, start, end);
  const s = ends.start.join(' / '), e = ends.end.join(' / ');
  if (!matches.length && !near && !s && !e) return null;
  const pad = `± ${LCR_PAD / 1000} kb`;
  const verdict = s && e && s !== e ? <>Ends in <b>{s}</b> and <b>{e}</b> ({pad}): consistent with an LCR-mediated recurrent CNV.</>
    : s && e ? <>Both ends in <b>{s}</b> ({pad}): within one breakpoint segment.</>
    : s || e ? <>{s ? <>Start in <b>{s}</b></> : 'Start in no breakpoint / LCR segment'}, {e ? <>end in <b>{e}</b></> : 'end in none'} ({pad}): one end only at a recurrent breakpoint.</>
    : <>Neither end in a breakpoint / LCR segment of the list ({pad}): not the recurrent breakpoints.</>;
  return <div data-cgh-recurrent>
    <div className="font-semibold text-slate-900">ClinGen recurrent CNVs (list v2.0)</div>
    {matches.map(m => <div key={m.name}>{shortName(m.name)}: <b>{pct(m.ofRegion)}</b> of the region in the call, <b>{pct(m.ofCall)}</b> of the call{m.piev.includes(type) && <span className="ml-1 rounded-sm border px-0.5 text-[9px] font-bold" style={{ color: COL.piev, borderColor: COL.piev }} title={`AChro-Puce classes this ${type === 'del' ? 'deletion' : 'duplication'} PIEV: susceptibility, reduced penetrance`}>PIEV</span>}</div>)}
    {!matches.length && <div className="text-slate-500">No recurrent region sharing half of the call and half of itself.</div>}
    <div>{verdict}</div>
  </div>;
}
