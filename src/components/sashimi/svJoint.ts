/**
 * The base-level joint of a junction, from its split reads: the bases inserted between the two sides, or the bases
 * they share (microhomology), and whether an insertion is a copy of nearby sequence (templated).
 *
 * For each split read whose primary alignment and one part of its SA tag meet at the junction's two breakends:
 * - the read is turned so that it runs along the derivative, the part at breakend A first;
 * - bases of the read between the two parts are inserted at the joint; bases both parts aligned (they overlap on the
 *   read) are shared by the two sides;
 * the reads vote, and the majority joint is checked on the reference: the homology at the joint is how far the
 * breakpoint can slide without changing the derivative's sequence (HOMLEN / HOMSEQ in Manta and DELLY), whatever the
 * aligner did with the shared bases. An insertion of TEMPLATE_MIN bases or more is looked up within TEMPLATE_SPAN of
 * either breakend, both strands, at most one mismatch per 10 bases.
 *
 * The mechanism hint follows the usual reading of breakpoint junctions (Hastings et al., Nat Rev Genet 2009;
 * Ottaviani et al., EMBO J 2014; Lee, Carvalho & Lupski, Cell 2007): blunt or 1–4 shared bases, end joining (NHEJ or
 * alt-EJ); 5 shared bases or more, microhomology-mediated end joining; a few non-templated inserted bases, NHEJ; a
 * templated insertion, replication-based repair (FoSTeS / MMBIR). A hint, not a call.
 */
import type { AlignedRead } from './types';
import type { Breakend, SvJunction } from './svModel';

export type Strand = '+' | '-';
export interface SvJoint {
  /** split reads read at this joint, and those agreeing with the call */
  reads: number; agree: number;
  /** bases inserted at the joint, along the derivative from A to B */
  inserted: string;
  /** bases the two sides share at the joint (the breakpoint can sit anywhere along them) */
  homology: string;
  /** the joint without the shared bases counted twice: A keeps up to posA, B from posB (0-based, as the breakends) */
  posA: number; posB: number;
  /** how far the joint slides back and on along the shared bases */
  slideLeft: number; slideRight: number;
  /** the derivative around the joint, from A to B: A's flank, the inserted or shared bases, B's flank */
  left: string; mid: string; right: string;
  /** the stretch an insertion copies (0-based half-open) */
  templated: { chrom: string; start: number; end: number; strand: Strand; mismatches: number } | null;
  /** the joints the reads gave, most voted first */
  votes: { key: string; n: number }[];
}

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
export const revComp = (s: string) => { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] ?? 'N'; return o; };

/** Inserted bases looked up as a copy of nearby sequence from this length, within this distance of a breakend. */
export const TEMPLATE_MIN = 10, TEMPLATE_SPAN = 5000;
/** A part of a split read is taken at a breakend when it ends within this of it. */
const AT_BREAKEND_BP = 20;
const FLANK = 20;

/** What the joint needs of a read: its primary alignment's geometry, its whole sequence as stored, its SA tag. */
export interface JointRead { name: string; chrom: string; start: number; end: number; rev: boolean; lead: number; trail: number; seq: string; flags: number; sa: string }

/**
 * The split reads of a reads window, with their whole sequence rebuilt from the reference, their mismatches,
 * insertions and soft-clipped bases: primary records only (a supplementary one is hard-clipped), with an SA tag and
 * their clipped bases at hand.
 */
export function jointReads(reads: AlignedRead[], chrom: string, ref: { start: number; seq: string } | null): JointRead[] {
  if (!ref) return [];
  const out: JointRead[] = [];
  for (const r of reads) {
    if (!r.sa || r.f & 0x900 || (r.c[0] + r.c[1] > 0 && !r.cs) || r.h?.[0] || r.h?.[1]) continue;
    // clipped bases cut short (a container or page keeps at most so many per side): the read's coordinates would slip
    if (r.cs && (r.cs[0].length !== r.c[0] || r.cs[1].length !== r.c[1])) continue;
    const mism = new Map(r.m.map(([p, b]) => [p, b]));
    const ins = new Map(r.i.map(([p], k) => [p, r.is?.[k] ?? null]));
    if ([...ins.values()].some(v => v == null)) continue;
    // the aligned bases: the reference with the read's mismatches, its inserted bases before the base they precede
    let aligned = '', ok = true;
    const used = new Set<number>();
    const insertAt = (p: number) => { if (ins.has(p) && !used.has(p)) { used.add(p); aligned += ins.get(p)!; } };
    for (const [bs, be] of r.b) {
      for (let p = bs; p < be; p++) {
        const k = p - ref.start;
        if (k < 0 || k >= ref.seq.length) { ok = false; break; }
        insertAt(p);
        aligned += mism.get(p) ?? ref.seq[k];
      }
      if (!ok) break;
      insertAt(be);
    }
    if (!ok || used.size !== ins.size) continue;
    const seq = (r.cs?.[0] ?? '') + aligned + (r.cs?.[1] ?? '');
    out.push({ name: r.n, chrom, start: r.s, end: r.e, rev: r.r === 1, lead: r.c[0], trail: r.c[1], seq: seq.toUpperCase(), flags: r.f, sa: r.sa });
  }
  return out;
}

interface Part { chrom: string; start: number; end: number; rev: boolean; qs: number; qe: number }
/** A part of the read from its CIGAR (an SA entry), query coordinates in the read's original orientation. */
function partOf(chrom: string, start: number, rev: boolean, cigar: string): Part {
  let refLen = 0, qLen = 0, lead = 0, trail = 0, seen = false;
  for (const m of cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) {
    const n = +m[1], op = m[2];
    if (op === 'S' || op === 'H') { if (!seen) lead += n; else trail += n; qLen += n; continue; }
    if ('MI=X'.includes(op)) { qLen += n; seen = true; }
    if ('MDN=X'.includes(op)) { refLen += n; seen = true; }
  }
  return { chrom, start, end: start + refLen, rev, qs: rev ? trail : lead, qe: qLen - (rev ? lead : trail) };
}
const sameChrom = (a: string, b: string) => a === b || a.replace(/^chr/i, '') === b.replace(/^chr/i, '');

/**
 * The joint of the junction a → b from the split reads around its two breakends. `refAt` gives the reference (+
 * strand, upper case) of a stretch, null outside what was read. Null when no split read meets both breakends.
 */
export function jointOf(a: Breakend, b: Breakend, reads: JointRead[], refAt: (chrom: string, start: number, end: number) => string | null): SvJoint | null {
  const touches = (p: Part, e: Breakend) => sameChrom(p.chrom, e.chrom) && (e.side === 'tail' ? Math.abs(p.end - e.pos) <= AT_BREAKEND_BP : Math.abs(p.start - e.pos) <= AT_BREAKEND_BP);
  // the derivative reads A on its strand (a tail forward up to pos, a head backward down to pos), B the other way round
  const strandA: Strand = a.side === 'tail' ? '+' : '-', strandB: Strand = b.side === 'head' ? '+' : '-';
  const votes = new Map<string, { n: number; ins: string; ovl: number; posA: number; posB: number }>();
  const seen = new Set<string>();
  let reads_ = 0;
  for (const r of reads) {
    if (seen.has(r.name)) continue;
    const L = r.seq.length;
    const parts: Part[] = [{ chrom: r.chrom, start: r.start, end: r.end, rev: r.rev, qs: r.rev ? r.trail : r.lead, qe: L - (r.rev ? r.lead : r.trail) }];
    for (const x of r.sa.split(';').filter(Boolean)) {
      const f = x.split(',');
      if (f.length >= 4) parts.push(partOf(f[0], parseInt(f[1], 10) - 1, f[2] === '-', f[3]));
    }
    const pa = parts.find(p => touches(p, a)), pb = parts.find(p => p !== pa && touches(p, b));
    if (!pa || !pb) continue;
    // the read turned so that it runs along the derivative (A's part on A's strand), A's part first
    const same = (pa.rev ? '-' : '+') === strandA;
    const orig = r.rev ? revComp(r.seq) : r.seq, dread = same ? orig : revComp(orig);
    const flip = (p: Part) => (same ? p : { ...p, qs: L - p.qe, qe: L - p.qs });
    const A = flip(pa), B = flip(pb);
    if (B.qs < A.qs) continue;
    seen.add(r.name); reads_++;
    const gap = B.qs - A.qe;
    const ins = gap > 0 ? dread.slice(A.qe, B.qs) : '';
    const ovl = gap < 0 ? -gap : 0;
    const posA = strandA === '+' ? A.end : A.start, posB = strandB === '+' ? B.start : B.end;
    const key = `${ins}|${ovl}|${posA}|${posB}`;
    const v = votes.get(key) ?? { n: 0, ins, ovl, posA, posB };
    v.n++; votes.set(key, v);
  }
  if (!reads_) return null;
  const ranked = [...votes.values()].sort((x, y) => y.n - x.n);
  const best = ranked[0];
  const cutA = strandA === '+' ? best.posA - best.ovl : best.posA + best.ovl;
  const base = (chrom: string, p: number, strand: Strand) => { const s = refAt(chrom, p, p + 1); return !s ? '?' : strand === '+' ? s : (COMP[s] ?? 'N'); };
  // A's k-th base before the joint (k ≥ 1; k ≤ 0 runs past it), B's k-th base after it (k ≥ 0; k < 0 before it)
  const aBefore = (k: number) => (strandA === '+' ? base(a.chrom, cutA - k, '+') : base(a.chrom, cutA + k - 1, '-'));
  const bAfter = (k: number) => (strandB === '+' ? base(b.chrom, best.posB + k, '+') : base(b.chrom, best.posB - 1 - k, '-'));
  let slideLeft = 0, slideRight = 0;
  if (!best.ins) {
    while (slideLeft < 200 && aBefore(1 + slideLeft) !== '?' && aBefore(1 + slideLeft) === bAfter(-1 - slideLeft)) slideLeft++;
    while (slideRight < 200 && aBefore(-slideRight) !== '?' && aBefore(-slideRight) === bAfter(slideRight)) slideRight++;
  }
  const homology = !best.ins ? Array.from({ length: slideLeft + slideRight }, (_, k) => aBefore(slideLeft - k)).join('') : '';
  const left = Array.from({ length: FLANK }, (_, k) => aBefore(slideLeft + FLANK - k)).join('');
  const right = Array.from({ length: FLANK }, (_, k) => bAfter(slideRight + k)).join('');
  let templated: SvJoint['templated'] = null;
  if (best.ins.length >= TEMPLATE_MIN && !/N/.test(best.ins)) {
    const maxMis = Math.floor(best.ins.length / 10);
    for (const e of [a, b]) {
      const lo = Math.max(0, e.pos - TEMPLATE_SPAN), text = refAt(e.chrom, lo, e.pos + TEMPLATE_SPAN);
      if (!text) continue;
      for (const [strand, q] of [['+', best.ins], ['-', revComp(best.ins)]] as const) {
        for (let i = 0; i + q.length <= text.length; i++) {
          let mis = 0;
          for (let k = 0; k < q.length && mis <= maxMis; k++) if (text[i + k] !== q[k]) mis++;
          if (mis <= maxMis && (!templated || mis < templated.mismatches)) templated = { chrom: e.chrom, start: lo + i, end: lo + i + q.length, strand, mismatches: mis };
        }
      }
    }
  }
  return {
    reads: reads_, agree: best.n, inserted: best.ins, homology, posA: cutA, posB: best.posB, slideLeft, slideRight,
    left, mid: best.ins || homology, right, templated,
    votes: ranked.map(v => ({ key: v.ins ? `+${v.ins.length > 12 ? `${v.ins.length} bp` : v.ins}` : v.ovl ? `${v.ovl} shared` : 'blunt', n: v.n })),
  };
}

/** What the joint says of the repair that made it: one sentence, a hint (see the module note). */
export function jointHint(j: SvJoint): string {
  if (j.templated) return 'a templated insertion, as replication-based repair (FoSTeS / MMBIR) leaves';
  if (j.inserted) return j.inserted.length <= 20 ? 'a few non-templated bases, as end joining (NHEJ) leaves' : 'a longer non-templated insertion: its origin is not near the breakpoints';
  const h = j.homology.length;
  if (h >= 5) return `${h} bp of microhomology: microhomology-mediated end joining (MMEJ) is likely`;
  if (h > 0) return `${h} bp of microhomology: end joining (NHEJ or alt-EJ)`;
  return 'a blunt joint: end joining (NHEJ)';
}

/** The joint in a few words, for the junction list. */
export function jointSummary(j: SvJoint): string {
  if (j.templated) return `+${j.inserted.length} bp templated`;
  if (j.inserted) return `+${j.inserted.length} bp inserted`;
  if (j.homology) return `${j.homology.length} bp shared`;
  return 'blunt';
}

/**
 * The junction with its joint: the breakends where the reads put them, both sides kept forward (a deletion, a tandem
 * duplication) moved as far 3′ as the shared bases allow (the HGVS 3′ rule), the joint carried for the nomenclature.
 */
export function refineJunction(j: SvJunction, joint: SvJoint | null | undefined): SvJunction {
  if (!joint) return j;
  const forward = j.a.side === 'tail' && j.b.side === 'head' && !joint.inserted;
  const shift = forward ? joint.slideRight : 0;
  return { ...j, a: { ...j.a, pos: joint.posA + shift }, b: { ...j.b, pos: joint.posB + shift }, precision: 0, joint };
}
