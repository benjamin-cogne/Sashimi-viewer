/**
 * The reads that support an arc, for "Show supporting reads" (the arc's panel): the reads track then holds only them
 * (and their mates), downsampled when there are many.
 *
 * - **junction** (RNA): an intron of the read (CIGAR N) at the junction, to the base, or within `tol` at both ends
 *   (long reads, whose junctions a few bases off were counted with it: junctionSnap.ts).
 * - **deletion**: a deletion in the CIGAR (≥ SV_MIN_DELETION) with both ends within `tol`, or a read broken at either
 *   end of it (below).
 * - **split, duplication, inversion**: a read broken at either end of the arc: an alignment end followed by a clip of
 *   RESCUE_MIN_CLIP bases or more, or an end of one of its parts (SA tag), within `tol`. That is the split and clipped
 *   reads the arc counts, placed or rescued, and a few that only share a breakpoint with it.
 * - **discordant** pairs: mates on the chromosome of the arc, of the arc's class (→ ← apart, ← → facing away, one
 *   strand), the left one in the arc's first 500 bp bins and the right one in its last.
 *
 * Coordinates 0-based half-open, as the arcs'.
 */
import type { AlignedRead } from '../components/sashimi/types';

export type ArcSupportKind = 'junction' | 'deletion' | 'split' | 'duplication' | 'inversion' | 'discordant';
export interface ArcSupport {
  kind: ArcSupportKind; start: number; end: number;
  /** discordant pairs: their class */
  pairKind?: 'deletion' | 'duplication' | 'inversion';
  /** bases either end may be off (0: exact) */
  tol: number;
}
/** Soft or hard clip after which an alignment end counts as a break (as the rescue of clipped reads: alignments.ts). */
const MIN_CLIP = 8;
const MIN_DELETION = 50;
/**
 * A discordant arc's ends are where its pairs' reads point (or, in files exported before, the edges of 500 bp bins
 * joined within one bin): a supporting pair has a read starting within this of each end.
 */
const PAIR_SLACK = 1500;

/** What the tests need of a read, from a record (CIGAR string) or from a decoded read. */
interface ArcRead {
  start: number; end: number;
  /** introns (N) and deletions (D) as [start, end, 'N' | 'D'] */
  gaps: [number, number, 'N' | 'D'][];
  clipL: number; clipR: number;
  sa: string | null;
  flags: number;
  rev: boolean; mateRev: boolean;
  matePos?: number; mateChrom?: string;
}

export function arcReadFromCigar(r: { start: number; cigar: string; flags: number; sa?: string | null; matePos?: number; mateChrom?: string }): ArcRead {
  const gaps: [number, number, 'N' | 'D'][] = [];
  let pos = r.start, clipL = 0, clipR = 0, seen = false;
  const re = /(\d+)([MIDNSHP=X])/g;
  let m: RegExpExecArray | null;
  while ((m = re.exec(r.cigar))) {
    const len = +m[1], op = m[2];
    if (op === 'S' || op === 'H') { if (seen) clipR += len; else clipL += len; continue; }
    if (op === 'M' || op === '=' || op === 'X') { pos += len; seen = true; }
    else if (op === 'N' || op === 'D') { gaps.push([pos, pos + len, op]); pos += len; seen = true; }
    else if (op === 'I') seen = true;
  }
  return { start: r.start, end: pos, gaps, clipL, clipR, sa: r.sa ?? null, flags: r.flags, rev: (r.flags & 16) !== 0, mateRev: (r.flags & 32) !== 0, matePos: r.matePos, mateChrom: r.mateChrom };
}

export function arcReadFromAligned(r: AlignedRead, chrom: string): ArcRead {
  const gaps: [number, number, 'N' | 'D'][] = r.d.map(([a, b]) => [a, b, 'D']);
  for (let k = 0; k + 1 < r.b.length; k++) {
    const be = r.b[k][1], ns = r.b[k + 1][0];
    if (ns > be && !r.d.some(([a, b]) => a === be && b === ns)) gaps.push([be, ns, 'N']);
  }
  return { start: r.s, end: r.e, gaps, clipL: r.c[0] + (r.h?.[0] ?? 0), clipR: r.c[1] + (r.h?.[1] ?? 0), sa: r.sa ?? null, flags: r.f, rev: r.r === 1, mateRev: (r.f & 32) !== 0,
    matePos: r.mp, mateChrom: r.mp == null ? undefined : (r.mc ?? chrom) };
}

const near = (x: number, y: number, tol: number) => Math.abs(x - y) <= tol;
const sameChrom = (a: string, b: string) => a === b || a.replace(/^chr/i, '') === b.replace(/^chr/i, '');

/** Where a read breaks: its alignment ends followed by a clip, and the ends of its other parts (SA tag) on `chrom`. */
function breakpoints(r: ArcRead, chrom: string): number[] {
  const out: number[] = [];
  if (r.clipL >= MIN_CLIP) out.push(r.start);
  if (r.clipR >= MIN_CLIP) out.push(r.end);
  if (r.sa) for (const part of r.sa.split(';')) {
    const f = part.split(',');
    if (f.length < 4 || !sameChrom(f[0], chrom)) continue;
    const s = parseInt(f[1]) - 1;
    if (!Number.isFinite(s)) continue;
    let len = 0;
    for (const [, n, op] of f[3].matchAll(/(\d+)([MIDNSHP=X])/g)) if ('MDN=X'.includes(op)) len += +n;
    out.push(s, s + len);
  }
  return out;
}

export function supportsArc(r: ArcRead, chrom: string, a: ArcSupport): boolean {
  switch (a.kind) {
    case 'junction':
      return r.gaps.some(([s, e, op]) => op === 'N' && near(s, a.start, a.tol) && near(e, a.end, a.tol));
    case 'deletion':
      if (r.gaps.some(([s, e, op]) => op === 'D' && e - s >= MIN_DELETION && near(s, a.start, a.tol) && near(e, a.end, a.tol))) return true;
      return breakpoints(r, chrom).some(p => near(p, a.start, a.tol) || near(p, a.end, a.tol));
    case 'split': case 'duplication': case 'inversion':
      return breakpoints(r, chrom).some(p => near(p, a.start, a.tol) || near(p, a.end, a.tol));
    case 'discordant': {
      if (r.matePos == null || !r.mateChrom || !sameChrom(r.mateChrom, chrom) || !(r.flags & 1)) return false;
      const lo = Math.min(r.start, r.matePos), hi = Math.max(r.start, r.matePos);
      if (Math.abs(lo - a.start) > PAIR_SLACK || Math.abs(hi - a.end) > PAIR_SLACK) return false;
      // as far apart as the evidence asks (alignments.ts: 1 kb for mates facing each other, 300 bp facing away) and
      // spanning at least half the arc: an ordinary pair between the two ends is none of its
      const cls0 = r.rev === r.mateRev ? 'inversion' : (r.start <= r.matePos ? r.rev : r.mateRev) ? 'duplication' : 'deletion';
      if (hi - lo < Math.max(cls0 === 'duplication' ? 300 : 1000, 0.5 * (a.end - a.start))) return false;
      const leftRev = r.start <= r.matePos ? r.rev : r.mateRev, rightRev = r.start <= r.matePos ? r.mateRev : r.rev;
      const cls = r.rev === r.mateRev ? 'inversion' : leftRev && !rightRev ? 'duplication' : 'deletion';
      return !a.pairKind || cls === a.pairKind;
    }
  }
}
