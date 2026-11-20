/**
 * The ClinGen CNV regions of the array-CGH browser's track and of a call's detail, from embedded data
 * (cghClinGenData.ts: HI 3 / TS 3 dosage regions and the recurrent CNV list v2.0 with its breakpoint / LCR segments).
 * Embedded so the track needs no request; the priority rules keep reading dosage regions live (cghAnnotate.ts).
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { CLINGEN_DOSAGE_GRCh37, CLINGEN_DOSAGE_GRCh38, CLINGEN_RECURRENT_GRCh37, CLINGEN_RECURRENT_GRCh38 } from './cghClinGenData';
import { ACHRO_REGIONS } from './achropuce';

/** 0-based half-open */
export interface ClinGenDosage { kind: 'dosage'; chrom: string; start: number; end: number; id: string; name: string; hi: number | null; ts: number | null; band: string; evaluated: string; url: string; piev: ('del' | 'dup')[]; pievName?: string }
export interface ClinGenRecurrent { kind: 'recurrent' | 'lcr'; chrom: string; start: number; end: number; name: string; piev: ('del' | 'dup')[]; pievName?: string }

/** AChro-Puce PIEV match: reciprocal overlap at least this share */
export const PIEV_RO = 0.8;
/** a call's end "in" a breakpoint segment when within this distance: arrays have few probes inside LCRs, so calls stop at their edge */
export const LCR_PAD = 50_000;

const ro = (a0: number, a1: number, b0: number, b1: number) => { const o = Math.min(a1, b1) - Math.max(a0, b0); return o <= 0 ? 0 : Math.min(o / (a1 - a0), o / (b1 - b0)); };
/** the AChro-Puce PIEV rows (of the given types) of a sheet lying on the region (its coordinates in this build only: no lift-over here) */
function pievOf(build: GenomeBuild, chrom: string, s: number, e: number, types: ('del' | 'dup')[]) {
  for (const r of ACHRO_REGIONS) {
    const co = r.coords?.find(x => x.build === build && x.chrom === chrom);
    if (!co || ro(s, e, co.start - 1, co.end) < PIEV_RO) continue;
    const t = types.filter(ty => r.rows.some(w => w.type === ty && w.cls === 'PIEV'));
    if (t.length) return { piev: t, pievName: r.name };
  }
  return { piev: [] as ('del' | 'dup')[] };
}

/** the 22q11.2 breakpoints as most reports name them (the list says BP-A … BP-H) */
export const lcrName = (chrom: string, n: string) => { const m = chrom === 'chr22' && n.match(/^BP-([A-H]) \(/); return m ? `LCR22${m[1]}` : n; };
/** a long ClinGen name shortened for a label: no "(includes …)", no "recurrent" / "region" */
export const shortName = (n: string) => { const t = n.replace(/\s*\(includes[^)]*\)/g, '').replace(/\s+(recurrent|region)\b/g, '').replace(/\s+/g, ' ').trim(); return t.length > 52 ? `${t.slice(0, 50)}…` : t; };

const cache = new Map<string, (ClinGenDosage | ClinGenRecurrent)[]>();
export function clinGenOn(build: GenomeBuild, chrom: string): (ClinGenDosage | ClinGenRecurrent)[] {
  const k = `${build}|${chrom}`, hit = cache.get(k); if (hit) return hit;
  const c = chrom.replace(/^chr/, ''), b38 = build === 'GRCh38';
  const dos: ClinGenDosage[] = (b38 ? CLINGEN_DOSAGE_GRCh38 : CLINGEN_DOSAGE_GRCh37).filter(r => r[0] === c).map(([, s, e, n, name, hi, ts, band, evaluated]) => ({
    kind: 'dosage', chrom, start: s, end: e, id: `ISCA-${n}`, name, hi, ts, band, evaluated, url: `https://search.clinicalgenome.org/kb/gene-dosage/region/ISCA-${n}`,
    ...pievOf(build, chrom, s, e, [...(hi === 3 ? ['del' as const] : []), ...(ts === 3 ? ['dup' as const] : [])]),
  }));
  const rec: ClinGenRecurrent[] = (b38 ? CLINGEN_RECURRENT_GRCh38 : CLINGEN_RECURRENT_GRCh37).filter(r => r[0] === c).map(([, s, e, name, t]) => ({
    kind: t === 'b' ? 'lcr' : 'recurrent', chrom, start: s, end: e, name: t === 'b' ? lcrName(chrom, name) : name,
    ...(t === 'b' ? { piev: [] } : pievOf(build, chrom, s, e, ['del', 'dup'])),
  }));
  const out = [...dos, ...rec];
  cache.set(k, out);
  return out;
}

export interface RecurrentMatch { name: string; ofRegion: number; ofCall: number; piev: ('del' | 'dup')[] }
export interface RecurrentEnds { start: string[]; end: string[] }
/**
 * A call (1-based inclusive) against the recurrent list: the regions sharing at least half of it and of themselves,
 * and the breakpoint / LCR segments (within LCR_PAD) each end falls in.
 */
export function recurrentOf(build: GenomeBuild, chrom: string, start1: number, end1: number): { matches: RecurrentMatch[]; ends: RecurrentEnds; near: boolean } {
  const s = start1 - 1, e = end1, all = clinGenOn(build, chrom);
  const matches = all.filter((r): r is ClinGenRecurrent => r.kind === 'recurrent').flatMap(r => {
    const o = Math.min(e, r.end) - Math.max(s, r.start);
    return o > 0 && o / (r.end - r.start) >= 0.5 && o / (e - s) >= 0.5 ? [{ name: r.name, ofRegion: o / (r.end - r.start), ofCall: o / (e - s), piev: r.piev }] : [];
  }).sort((a, b) => Math.min(b.ofRegion, b.ofCall) - Math.min(a.ofRegion, a.ofCall)).slice(0, 4);
  const lcr = all.filter(r => r.kind === 'lcr');
  const at = (p: number) => [...new Set(lcr.filter(r => p >= r.start - LCR_PAD && p < r.end + LCR_PAD).map(r => r.name))];
  // `near`: a recurrent region overlaps the call, so the ends are worth stating even when they miss every segment
  return { matches, ends: { start: at(s), end: at(e - 1) }, near: all.some(r => r.kind === 'recurrent' && r.start < e && r.end > s) };
}
