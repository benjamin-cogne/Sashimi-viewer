/**
 * Array ISCN (ISCN 2024, McGowan-Jordan, Hastings, Moore eds., Karger, doi:10.1159/000538512, chapter 8 "Microarray"),
 * as in its examples `arr[GRCh38] 1p36.11p35.3(25,889,422_28,493,687)x3mat,11p12p11.12(41,444,865_49,146,935)x1dn[0.5]`
 * and `arr[GRCh38] Xp22.33q28(168,546_155,233,730)x1[0.6]`: the build, the chromosome, the bands from pter to qter (the
 * band at the start first; one band written once), the positions with thousands commas, the copy number, a mosaic
 * fraction in brackets; several abnormalities of a patient separated by commas, by chromosome (1-22, X, Y) then position.
 * The copy number: the file's, else the one inferred from the log ratio (cghMosaic.ts `estCopies`); the mosaic fraction
 * when the AChro-Puce estimate classes the call as mosaic (one-copy change in < 80 % of cells).
 */
import type { CghCnv, CghSex } from './cghIntervals';
import { estCopies, mosaicOf, MOSAIC_TAG_BELOW } from './cghMosaic';

const ORDER = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
const num = (n: number) => n.toLocaleString('en-US');
const chr = (c: CghCnv) => c.chrom.replace(/^chr/i, '');

/** the band of a call from the cytogenetic bands (UCSC cytoBand, 0-based half-open), at their finest resolution:
 * `q11.23`, or `q11.21 - q11.23` (the band at the start first, so pter to qter); undefined when a position has none */
export function bandSpan(bands: { start: number; end: number; name: string }[], s0: number, e0: number): string | undefined {
  const at = (pos: number) => bands.find(b => b.start <= pos && pos < b.end)?.name;
  const a = at(s0), z = at(Math.max(s0, e0 - 1));
  return !a || !z ? undefined : a === z ? a : `${a} - ${z}`;
}
/** the band written in a call's ISCN: the one read from the cytogenetic bands (sub-band precision), else the file's */
export const iscnBandOf = (c: CghCnv, ucscBand: string | undefined) => ucscBand || c.cytoband || '';
/** the bands of a call as ISCN writes them: `q11.21q11.23` (first at the start), a single band once; '' when none */
export function iscnBands(band: string): string {
  const b = band.replace(/\s+/g, '').match(/[pq]\d+(?:\.\d+)?|cen/g);
  if (!b?.length) return '';
  return b[0] === b[b.length - 1] ? b[0] : `${b[0]}${b[b.length - 1]}`;
}
/** a fraction as ISCN brackets it: 0.5, 0.38 */
const frac = (f: number) => String(Math.round(f * 100) / 100);

export interface CallIscn {
  /** this call written alone (`arr[GRCh37] 22q11.21(18,894,835_21,505,417)x1`), null without a band nor a copy number */
  text: string | null;
  /** the part after `arr[build] `, to join a patient's calls */
  part: string | null;
  /** the file's own ISCN, and whether ours differs from it (build names, spaces and commas aside) */
  file: string | null;
  differs: boolean;
}
/** the copy number and the mosaic fraction written for a call: the file's copies, else inferred from the log ratio */
export function iscnCopies(c: CghCnv, sex: CghSex | undefined): { copies: string | null; mosaic: number | null } {
  const m = mosaicOf(c, sex), mos = m?.ok && m.f < MOSAIC_TAG_BELOW ? m.f : null;
  if (c.copies != null) return { copies: `x${c.copies}${c.copiesHi != null ? `~${c.copiesHi}` : ''}`, mosaic: c.copiesHi == null ? mos : null };
  if (c.amp) return { copies: 'amp', mosaic: null };
  const e = estCopies(c, sex);
  return e ? { copies: `x${e.n}`, mosaic: e.mosaic } : { copies: null, mosaic: null };
}
const norm = (t: string) => t.replace(/\b(hg19|GRCh37)\b/gi, 'GRCh37').replace(/\b(hg38|GRCh38)\b/gi, 'GRCh38').replace(/[\s,]/g, '').toLowerCase();

/** a call's ISCN; `whole`: the call spans the whole chromosome (written `(21)x3`) */
export function callIscn(c: CghCnv, build: string, band: string, sex: CghSex | undefined, whole = false): CallIscn {
  const { copies, mosaic } = iscnCopies(c, sex), b = iscnBands(band);
  const where = whole ? `(${chr(c)})` : b ? `${chr(c)}${b}(${num(c.start)}_${num(c.end)})` : null;
  const part = where && copies ? `${where}${copies}${mosaic != null ? `[${frac(mosaic)}]` : ''}` : null;
  const text = part ? `arr[${build}] ${part}` : null;
  return { text, part, file: c.iscn, differs: !!(c.iscn && text && norm(c.iscn) !== norm(text)) };
}

/** a patient's ISCN from its calls' parts, by chromosome then position; null when none */
export function patientIscn(build: string, calls: { c: CghCnv; part: string | null }[]): string | null {
  const ok = calls.filter(x => x.part).sort((a, b) => ORDER.indexOf(chr(a.c)) - ORDER.indexOf(chr(b.c)) || a.c.start - b.c.start || a.c.end - b.c.end);
  return ok.length ? `arr[${build}] ${ok.map(x => x.part).join(',')}` : null;
}
