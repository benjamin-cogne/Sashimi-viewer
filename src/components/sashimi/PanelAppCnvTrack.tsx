/**
 * The "PanelApp CNV regions" track (SvGenomeBrowser's `extraTracks`): the known pathogenic CNVs of PanelApp Australia's
 * panel "Common deletion and duplication syndromes" (panel 3443), read live from its API (open to any origin). Its
 * regions are rated green, each a loss or a gain with its ClinGen ISCA id, HI / TS scores, mode of inheritance and
 * phenotypes; PanelApp gives their coordinates on GRCh38 only. On GRCh37 each region on the chromosome in view is
 * lifted over with Ensembl's assembly mapper (REST /map), from PanelApp's own GRCh38 coordinates (ClinGen's GRCh37
 * coordinates of the same ISCA id are not used: the panel's boundaries differ from ClinGen's for about a third of them).
 *
 * useCnvRegionTracks: this track (shown by default) and the ClinGen one (off, a toggle away), the same in the SV
 * browser, the SV cartoon and the array-CGH page.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import type { GenomeBuild } from '../../standalone/ensembl';
import type { SvExtraTrack } from './SvGenomeBrowser';
import { useClinGenTrack } from './CghClinGenTrack';

export const PANELAPP_CNV = { id: '3443', base: 'https://panelapp-aus.org', name: 'Common deletion and duplication syndromes' };
const API = `${PANELAPP_CNV.base}/api/v1/panels/${PANELAPP_CNV.id}/`;
const regionUrl = (entity: string) => `${PANELAPP_CNV.base}/panels/${PANELAPP_CNV.id}/region/${encodeURIComponent(entity)}/`;
const MAP_URL = 'https://rest.ensembl.org/map/human/GRCh38';

export interface PanelCnvRegion {
  entity: string; name: string; kind: 'loss' | 'gain';
  /** GRCh38, 0-based start, end exclusive (PanelApp gives 1-based inclusive) */
  chrom: string; start: number; end: number;
  hi: string | null; ts: string | null; moi: string; phenotypes: string[]; overlap: number | null; publications: string[];
}
export interface PanelCnvData { version: string; date: string; regions: PanelCnvRegion[] }

const score = (v: unknown) => (v == null || v === '' ? null : String(v));
/** The panel's green regions (confidence 3) with GRCh38 coordinates, as losses or gains (a "cnv_both" region is both). */
export function parsePanelCnv(p: any): PanelCnvData {
  const regions: PanelCnvRegion[] = [];
  for (const r of p?.regions ?? []) {
    const c = r?.grch38_coordinates;
    if (String(r?.confidence_level) !== '3' || !Array.isArray(c) || c.length !== 2 || !r?.chromosome) continue;
    const entity = String(r.entity_name ?? ''), t = String(r.type_of_variants ?? '');
    // the name says it (ISCA-37390-Loss is typed "small"); else the variant type
    const kinds: ('loss' | 'gain')[] = /-loss$/i.test(entity) ? ['loss'] : /-gain$/i.test(entity) ? ['gain']
      : t === 'cnv_loss' ? ['loss'] : t === 'cnv_gain' ? ['gain'] : t === 'cnv_both' ? ['loss', 'gain'] : [];
    for (const kind of kinds) regions.push({
      entity, name: String(r.verbose_name || entity).replace(/\s+/g, ' ').trim(), kind, chrom: String(r.chromosome).replace(/^chr/i, ''),
      start: Math.min(c[0], c[1]) - 1, end: Math.max(c[0], c[1]),
      hi: score(r.haploinsufficiency_score), ts: score(r.triplosensitivity_score), moi: String(r.mode_of_inheritance ?? ''),
      phenotypes: Array.isArray(r.phenotypes) ? r.phenotypes.map((x: unknown) => String(x).replace(/\s+/g, ' ').trim()) : [], overlap: r.required_overlap_percentage ?? null,
      publications: Array.isArray(r.publications) ? r.publications.map(String) : [],
    });
  }
  // the panel lists ISCA-37418-Loss twice (ISCA_37418-Loss too): one per place and type
  const seen = new Set<string>();
  return {
    version: String(p?.version ?? ''), date: String(p?.version_created ?? '').slice(0, 10),
    regions: regions.filter(r => { const k = `${r.chrom}:${r.start}-${r.end}:${r.kind}`; if (seen.has(k)) return false; seen.add(k); return true; }),
  };
}

let panelPromise: Promise<PanelCnvData> | null = null;
/** The panel, fetched once per page (again after a failure). */
export function loadPanelCnv(): Promise<PanelCnvData> {
  if (!panelPromise) panelPromise = (async () => {
    const r = await fetch(API, { headers: { Accept: 'application/json' } });
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    return parsePanelCnv(await r.json());
  })().catch(e => { panelPromise = null; throw e; });
  return panelPromise;
}

/** a GRCh38 stretch on GRCh37: the span of its mapped pieces on the same chromosome, null when none maps */
export interface Lifted { start: number; end: number; pieces: number; share: number }
const liftCache = new Map<string, Promise<Lifted | null>>();
export function liftToGrch37(chrom: string, start: number, end: number): Promise<Lifted | null> {
  const k = `${chrom}:${start}-${end}`;
  if (!liftCache.has(k)) liftCache.set(k, (async () => {
    const r = await fetch(`${MAP_URL}/${chrom}:${start + 1}..${end}:1/GRCh37?content-type=application/json`, { headers: { Accept: 'application/json' } });
    if (!r.ok) throw new Error(`Ensembl HTTP ${r.status}`);
    const ms = ((await r.json())?.mappings ?? []).filter((m: any) => String(m?.mapped?.seq_region_name) === chrom);
    if (!ms.length) return null;
    const covered = ms.reduce((a: number, m: any) => a + (m.original.end - m.original.start + 1), 0);
    return { start: Math.min(...ms.map((m: any) => m.mapped.start)) - 1, end: Math.max(...ms.map((m: any) => m.mapped.end)), pieces: ms.length, share: covered / (end - start) };
  })().catch(e => { liftCache.delete(k); throw e; }));
  return liftCache.get(k)!;
}

// ======================== Drawing ========================

const COL = { loss: '#dc2626', gain: '#2563eb', text: '#1f2937', faint: '#8796a8' };
const ROW = 24;
const fmt = (n: number) => n.toLocaleString('en-US');
const size = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} Mb` : `${Math.round(n / 1e3)} kb`);
const tw = (s: string, px: number) => s.length * px * 0.56;
/** "Chromosome 22q11.2 deletion syndrome, distal, MIM#611867" → "22q11.2 deletion syndrome" */
export const shortCnvName = (n: string) => n.split(',')[0].replace(/^Chromosome\s+/i, '').replace(/\s*MIM#?\d+/g, '').trim();

/** a place on the chromosome in view, with its loss and / or gain */
interface Place { chrom: string; start: number; end: number; loss?: PanelCnvRegion; gain?: PanelCnvRegion; lifted?: Lifted }
type Item = { p: Place; a: number; b: number; lane: number; label: string };

function lanes(ps: Place[], x: (pos: number) => number, x0: number): Item[] {
  const ends: number[] = [];
  return ps.slice().sort((p, q) => p.start - q.start || q.end - p.end).map(p => {
    const a = x(p.start), b = Math.max(x(p.end), a + 2), left = Math.max(a, x0);
    const names = [p.loss, p.gain].filter(Boolean).map(r => shortCnvName(r!.name));
    let label = [...new Set(names)].join(' / ');
    if (label.length > 64) label = `${label.slice(0, 62)}…`;
    const right = Math.max(b, left + 2 + tw(label, 9.5) + 34);
    let lane = ends.findIndex(e => e + 8 < left);
    if (lane < 0) { lane = ends.length; ends.push(right); } else ends[lane] = right;
    return { p, a, b, lane, label };
  });
}

function tip(r: PanelCnvRegion, p: Place, build: GenomeBuild): string {
  const ph = r.phenotypes.slice(0, 4).join('; ') + (r.phenotypes.length > 4 ? '; …' : '');
  const at = build === 'GRCh37' && p.lifted
    ? `GRCh37 ${p.chrom}:${fmt(p.lifted.start + 1)}-${fmt(p.lifted.end)} (lifted over by Ensembl${p.lifted.pieces > 1 ? `, ${p.lifted.pieces} pieces` : ''}${p.lifted.share < 0.98 ? `, ${Math.round(p.lifted.share * 100)} % mapped` : ''}; PanelApp gives GRCh38 ${p.chrom}:${fmt(r.start + 1)}-${fmt(r.end)})`
    : `GRCh38 ${r.chrom}:${fmt(r.start + 1)}-${fmt(r.end)}`;
  return `${r.name}\n${r.entity} · ${r.kind === 'loss' ? 'deletion' : 'duplication'} · ${size(r.end - r.start)} · green (PanelApp Australia)\n`
    + `HI ${r.hi ?? '–'} · TS ${r.ts ?? '–'}${r.overlap != null ? ` · required overlap ${r.overlap} %` : ''}\n${r.moi}${ph ? `\n${ph}` : ''}\n${at}\nClick: its PanelApp page`;
}

type State = { data: PanelCnvData | null; error: string | null };
function Track({ items, build, state, lifting, x0, none }: { items: Item[]; build: GenomeBuild; state: State; lifting: boolean; x0: number; none: string }) {
  const lx = (it: Item) => Math.max(it.a, x0) + 2;
  const note = state.error ? `PanelApp Australia unreachable (${state.error}): the ClinGen CNV regions track (embedded, no request) still shows the curated regions`
    : !state.data ? `Loading PanelApp Australia panel ${PANELAPP_CNV.id}…`
    : lifting && !items.length ? 'Lifting the regions over to GRCh37 (Ensembl)…'
    : !items.length ? none : '';
  return (
    <g fontSize={9.5}>
      {note && <text x={x0 + 4} y={10} fill={state.error ? '#b45309' : COL.faint}>{note}</text>}
      {items.map(it => { const { p } = it, w = Math.max(it.b - it.a, 2), both = !!(p.loss && p.gain), first = (p.loss ?? p.gain)!;
        const bar = (r: PanelCnvRegion, y: number, h: number) => <a key={r.kind} href={regionUrl(r.entity)} target="_blank" rel="noreferrer">
          <title>{tip(r, p, build)}</title><rect x={it.a} y={y} width={w} height={h} rx={both ? 0 : 1.5} fill={COL[r.kind]} opacity={0.85} /></a>;
        return <g key={`${p.start}-${p.end}`} transform={`translate(0 ${it.lane * ROW})`}>
          {both ? <>{bar(p.loss!, 1, 5)}{bar(p.gain!, 6, 5)}</> : bar(first, 1, 10)}
          <a href={regionUrl(first.entity)} target="_blank" rel="noreferrer"><title>{[p.loss, p.gain].filter(Boolean).map(r => tip(r!, p, build)).join('\n\n')}</title>
            <text x={lx(it)} y={21} fill={COL.text}><tspan fontWeight={600}>{it.label}</tspan>
              {p.loss && <tspan dx={5} fontWeight={700} fill={COL.loss}>del</tspan>}{p.gain && <tspan dx={5} fontWeight={700} fill={COL.gain}>dup</tspan>}</text></a>
        </g>; })}
    </g>
  );
}

const sw = (style: React.CSSProperties) => <span className="inline-block w-4 h-2.5 rounded-sm" style={style} />;

/**
 * The PanelApp Australia track for a browser on `build`: fetches the panel, and on GRCh37 lifts over the regions of the
 * chromosome drawn. A new object each time data arrives (the browser redraws).
 */
export function usePanelAppCnvTrack(build: GenomeBuild): SvExtraTrack {
  const [state, setState] = useState<State>({ data: null, error: null });
  const [lifts, setLifts] = useState<Map<string, Lifted | null>>(() => new Map());
  const pending = useRef(new Set<string>());
  const alive = useRef(true);
  useEffect(() => { alive.current = true; return () => { alive.current = false; }; }, []);
  useEffect(() => {
    loadPanelCnv().then(data => alive.current && setState({ data, error: null }), e => alive.current && setState({ data: null, error: (e as Error).message || 'network error' }));
  }, []);

  return useMemo(() => {
    const key = (r: PanelCnvRegion) => `${r.chrom}:${r.start}-${r.end}`;
    /** lift the chromosome's regions not yet asked for (outside the render: the browser lays tracks out while drawing) */
    const askLifts = (rs: PanelCnvRegion[]) => {
      const todo = [...new Set(rs.map(key))].filter(k => !lifts.has(k) && !pending.current.has(k));
      if (!todo.length) return;
      todo.forEach(k => pending.current.add(k));
      setTimeout(() => {
        for (const k of todo) {
          const r = rs.find(x => key(x) === k)!;
          liftToGrch37(r.chrom, r.start, r.end).then(v => v, () => null).then(v => {
            pending.current.delete(k);
            if (alive.current) setLifts(m => new Map(m).set(k, v));
          });
        }
      }, 0);
    };
    const v = state.data;
    return {
      id: 'panelAppCnv', label: 'PanelApp CNV regions', side: ['PanelApp', 'CNV regions'], defaultOn: true,
      title: `Known pathogenic CNVs: PanelApp Australia panel ${PANELAPP_CNV.id} "${PANELAPP_CNV.name}"${v ? ` (v${v.version}, ${v.date})` : ''}, green regions, deletions red and duplications blue; read live from PanelApp${build === 'GRCh37' ? '; GRCh37: lifted over from the panel\'s GRCh38 coordinates by Ensembl' : ''}`,
      legend: <span className="inline-flex flex-wrap items-center gap-1">PanelApp Australia, <a className="underline" href={`${PANELAPP_CNV.base}/panels/${PANELAPP_CNV.id}/`} target="_blank" rel="noreferrer">{PANELAPP_CNV.name}</a>{v ? ` (v${v.version}, ${v.date}, ${v.regions.length} green regions)` : ''}:
        {sw({ background: COL.loss })} deletion {sw({ background: COL.gain })} duplication (click: its PanelApp page){build === 'GRCh37' ? '; lifted over to GRCh37 (Ensembl)' : ''}</span>,
      layout: ({ chrom, s, e, x, x0 }) => {
        const c = chrom.replace(/^chr/i, '');
        const on = v ? v.regions.filter(r => r.chrom === c) : [];
        if (build === 'GRCh37' && on.length) askLifts(on);
        const places = new Map<string, Place>();
        let lifting = false;
        for (const r of on) {
          let start = r.start, end = r.end, lifted: Lifted | undefined;
          if (build === 'GRCh37') {
            const l = lifts.get(key(r));
            if (l === undefined) { lifting = true; continue; }
            if (l === null) continue;
            lifted = l; start = l.start; end = l.end;
          }
          if (end <= s || start >= e) continue;
          const k = `${start}-${end}`, p = places.get(k) ?? { chrom: c, start, end, lifted };
          p[r.kind] = r; places.set(k, p);
        }
        const items = lanes([...places.values()], x, x0);
        const n = items.length ? Math.max(...items.map(i => i.lane)) + 1 : 0;
        const none = on.length ? 'no PanelApp CNV region in view' : `no PanelApp CNV region on ${chrom}`;
        return { height: n ? n * ROW : 16, node: <Track items={items} build={build} state={state} lifting={lifting} x0={x0} none={none} /> };
      },
    };
  }, [build, state, lifts]);
}

/** The known-CNV tracks of every browser: PanelApp Australia and ClinGen, both shown by default (each a toggle away from
 * hidden). `clinGenOn`: the array-CGH page's browser, whose ClinGen on / off choice is kept under its own id, apart from
 * the other browsers' */
export function useCnvRegionTracks(build: GenomeBuild, { clinGenOn = false }: { clinGenOn?: boolean } = {}): SvExtraTrack[] {
  const panelApp = usePanelAppCnvTrack(build), clinGen = useClinGenTrack(build);
  return useMemo(() => [panelApp, clinGenOn ? { ...clinGen, id: `${clinGen.id}.cghPage`, defaultOn: true } : { ...clinGen, defaultOn: true }], [panelApp, clinGen, clinGenOn]);
}
