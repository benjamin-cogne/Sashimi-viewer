/**
 * Browser-side port of backend/bam_reader.py: turn BAM/CRAM records into the compact
 * read encoding the viewer consumes, plus coverage runs and junction counts.
 * All coordinates are 0-based half-open.
 */
import type { AlignedRead, DiscordantArc, JunctionArc, StructuralEvidence, RealignedClip, Breakpoint, RescuedClips, ClipCluster, ElsewhereLink, SvArc } from '../components/sashimi/types';
import { clusterSv, type SvMember } from './svmerge';

/** Aligner-agnostic view of one record (BAM or CRAM). */
export interface RawRead {
  name: string;
  start: number;         // 0-based
  cigar: string;
  seq: string;           // read bases ('' when unavailable)
  qual: ArrayLike<number> | null;
  flags: number;
  mapq: number;
  nh: number | null;
  /** Mismatches already known from the record (CRAM substitution features); used when `seq` is empty. */
  mismatches?: [number, string, number][];
  /** pair and split-read fields, filled when structural evidence is wanted */
  tlen?: number;
  /** mate chromosome name ('' when unmapped or unknown) and 0-based mate start */
  mateChrom?: string; matePos?: number;
  /** SA tag: supplementary alignments "rname,pos,strand,CIGAR,mapQ,NM;" */
  sa?: string | null;
  /** haplotag of a phased file (WhatsHap / LongPhase haplotag, HiPhase, DRAGEN): HP haplotype, PS phase set, PC confidence */
  hp?: number | null; ps?: number | null; pc?: number | null;
  /** base-modification tags (MM / ML, MN), kept for the methylation colours of the reads; absent in light scans */
  mods?: { mm: string; ml: ArrayLike<number> | null; mn: number | null } | null;
}

const FLAG_PAIRED = 1, FLAG_UNMAPPED = 4, FLAG_REVERSE = 16, FLAG_READ2 = 128, FLAG_SECONDARY = 256, FLAG_QCFAIL = 512, FLAG_DUP = 1024;

/** Mismatches carried by CRAM read features: substitutions (X), base+quality (B) and base runs (b). */
export function cramMismatches(features: CramFeatureLike[] | undefined, qual: ArrayLike<number> | null): [number, string, number][] {
  const out: [number, string, number][] = [];
  for (const f of features || []) {
    const q = qual && qual.length > f.pos ? qual[f.pos] : 30;
    if (f.code === 'X' && f.sub) out.push([f.refPos, f.sub, q]);
    else if (f.code === 'B' && Array.isArray(f.data)) out.push([f.refPos, f.data[0], f.data[1] ?? q]);
    else if (f.code === 'b' && typeof f.data === 'string') for (let k = 0; k < f.data.length; k++) out.push([f.refPos + k, f.data[k], qual && qual.length > f.pos + k ? qual[f.pos + k] : 30]);
  }
  return out;
}

/** Same filter as pysam count_coverage(read_callback='all'): flags of the records left out. */
export const DROP_FLAGS = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_QCFAIL | FLAG_DUP;
export const keepFlags = (flags: number): boolean => (flags & DROP_FLAGS) === 0;
/** Uniquely mapped: NH:1 when the tag is there, else MAPQ ≥ 30. */
export const uniqueFrom = (nh: number | null, mapq: number): boolean => (nh != null ? nh === 1 : mapq >= 30);
export function keepRead(r: RawRead): boolean { return keepFlags(r.flags); }
export function isUnique(r: RawRead): boolean { return uniqueFrom(r.nh, r.mapq); }

const CIGAR_RE = /(\d+)([MIDNSHP=X])/g;

export function parseCigar(cigar: string): [number, string][] {
  const out: [number, string][] = [];
  for (const m of cigar.matchAll(CIGAR_RE)) out.push([parseInt(m[1]), m[2]]);
  return out;
}

/** Port of bam_reader._encode_read. `ref` covers [refStart, refStart+ref.length). */
export function encodeRead(r: RawRead, ref: string | null, refStart: number): AlignedRead {
  const blocks: [number, number][] = [], dels: [number, number][] = [], ins: [number, number][] = [];
  const mism: [number, string, number][] = [];
  const clips: [number, number] = [0, 0], hard: [number, number] = [0, 0], clipSeq: [string, string] = ['', ''], insSeq: string[] = [];
  let rpos = r.start, qpos = 0, first = true;   // `first`: no aligned base seen yet (leading clips, whatever their order, are the left ones)
  const refEnd = ref ? refStart + ref.length : 0;
  for (const [len, op] of parseCigar(r.cigar)) {
    if (op === 'M' || op === '=' || op === 'X') {
      blocks.push([rpos, rpos + len]);
      if (ref && r.seq) {
        const lo = Math.max(rpos, refStart), hi = Math.min(rpos + len, refEnd);
        if (hi > lo) {
          const qs = r.seq.substring(qpos + (lo - rpos), qpos + (hi - rpos));
          const rs = ref.substring(lo - refStart, hi - refStart);
          if (qs !== rs) {
            for (let k = 0; k < qs.length; k++) {
              const a = qs[k], b = rs[k];
              if (a !== b && b !== 'N') {
                const qi = qpos + (lo - rpos) + k;
                mism.push([lo + k, a, r.qual && r.qual.length > qi ? r.qual[qi] : 30]);
              }
            }
          }
        }
      }
      rpos += len; qpos += len; first = false;
    } else if (op === 'I') { ins.push([rpos, len]); insSeq.push(r.seq ? r.seq.substring(qpos, qpos + len) : ''); qpos += len; first = false; }
    else if (op === 'D') { dels.push([rpos, rpos + len]); rpos += len; first = false; }
    else if (op === 'N') { rpos += len; first = false; }
    else if (op === 'S') { const side = first ? 0 : 1; clips[side] += len; if (r.seq) clipSeq[side] += r.seq.substring(qpos, qpos + len); qpos += len; }
    else if (op === 'H') { hard[first ? 0 : 1] += len; }
  }
  if (!r.seq && r.mismatches && ref) {
    for (const m of r.mismatches) if (m[0] >= refStart && m[0] < refEnd && ref[m[0] - refStart] !== m[1] && ref[m[0] - refStart] !== 'N') mism.push(m);
  }
  const out: AlignedRead = {
    n: r.name, s: r.start, e: rpos, r: (r.flags & FLAG_REVERSE) ? 1 : 0, q: r.mapq, f: r.flags, nh: r.nh,
    b: blocks, d: dels, i: ins,
    m: mism,
    c: clips,
  };
  if (r.seq && (clips[0] || clips[1])) out.cs = clipSeq;
  if (hard[0] || hard[1]) out.h = hard;
  if (r.seq && ins.length) out.is = insSeq;
  if (r.sa) out.sa = r.sa;
  if (r.hp != null && r.hp > 0) { out.hp = r.hp; if (r.ps != null) out.ps = r.ps; if (r.pc != null) out.pc = r.pc; }
  return out;
}

/** The SA tag parsed: the other parts of a split read (0-based starts). */
export function parseSa(sa: string | undefined): { chrom: string; start: number; strand: '+' | '-'; cigar: string; mapq: number }[] {
  if (!sa) return [];
  const out: { chrom: string; start: number; strand: '+' | '-'; cigar: string; mapq: number }[] = [];
  for (const part of sa.split(';')) {
    const f = part.split(',');
    if (f.length < 4) continue;
    const start = parseInt(f[1]) - 1;
    if (!(start >= 0)) continue;
    out.push({ chrom: f[0], start, strand: f[2] === '-' ? '-' : '+', cigar: f[3], mapq: parseInt(f[4]) || 0 });
  }
  return out;
}

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', a: 't', c: 'g', g: 'c', t: 'a', n: 'n' };
export function reverseComplement(seq: string): string { let out = ''; for (let i = seq.length - 1; i >= 0; i--) out += COMP[seq[i]] ?? 'N'; return out; }

/**
 * The bases a supplementary record hard-clipped, taken from the read's primary record (which carries the whole
 * sequence, soft-clipped): the part's read-coordinate interval is cut out of the primary sequence, reverse-complemented
 * first when the two records are on opposite strands.
 */
export function hardClippedBases(part: { h?: [number, number]; c: [number, number]; cigar?: string; r: 0 | 1; queryLen: number }, primary: { seq: string; flags: number }): [string, string] | null {
  if (!primary.seq || !part.h) return null;
  const seq = ((primary.flags & FLAG_REVERSE) ? 1 : 0) === part.r ? primary.seq : reverseComplement(primary.seq);
  const total = part.h[0] + part.queryLen + part.h[1];
  if (seq.length !== total) return null;
  return [seq.substring(0, part.h[0]), seq.substring(total - part.h[1])];
}


/** Anchors of a boundary-spanning read: aligned bases required on the exon side and on the intron side of the boundary. */
export const SPAN_EXON_ANCHOR = 6, SPAN_INTRON_ANCHOR = 10;



/** Depth at a position from encoded reads (for site calling). */
export function depthArray(reads: AlignedRead[], start: number, end: number): Int32Array {
  const n = Math.max(0, end - start);
  const diff = new Int32Array(n + 1);
  for (const r of reads) for (const [bs, be] of r.b) {
    const lo = Math.max(bs, start), hi = Math.min(be, end);
    if (hi > lo) { diff[lo - start] += 1; diff[hi - start] -= 1; }
  }
  const out = new Int32Array(n);
  let d = 0;
  for (let i = 0; i < n; i++) { d += diff[i]; out[i] = d; }
  return out;
}

// ======================== CRAM read features → CIGAR + sequence ========================

export interface CramFeatureLike { code: string; pos: number; refPos: number; data?: any; sub?: string; ref?: string }

/**
 * Rebuild a CIGAR string from CRAM read features. `pos` is the read position of the
 * feature as exposed by @gmod/cram, which is 0-based in practice (a substitution at
 * pos 48 of a read starting at 115101 sits on reference base 115149, 0-based); the
 * walk below therefore treats `pos + 1` as the 1-based position of the specification.
 */
export function cramCigar(features: CramFeatureLike[] | undefined, readLength: number, lengthOnRef: number): string {
  if (!features || features.length === 0) return `${readLength}M`;
  features = features.map(f => ({ ...f, pos: f.pos + 1 }));
  const ops: [number, string][] = [];
  const push = (len: number, op: string) => {
    if (len <= 0) return;
    const last = ops[ops.length - 1];
    if (last && last[1] === op) last[0] += len; else ops.push([len, op]);
  };
  let readPos = 1; // next unaccounted 1-based read position
  for (const f of [...features].sort((a, b) => a.pos - b.pos)) {
    if (f.pos > readPos) push(f.pos - readPos, 'M');
    readPos = Math.max(readPos, f.pos);
    switch (f.code) {
      case 'S': push((f.data as string).length, 'S'); readPos += (f.data as string).length; break;
      case 'I': push((f.data as string).length, 'I'); readPos += (f.data as string).length; break;
      case 'i': push(1, 'I'); readPos += 1; break;
      case 'D': push(f.data as number, 'D'); break;
      case 'N': push(f.data as number, 'N'); break;
      case 'H': push(f.data as number, 'H'); break;
      case 'P': push(f.data as number, 'P'); break;
      case 'X': case 'B': case 'Q': push(1, 'M'); readPos += 1; break;   // substitution / base+qual / quality: consume one read base
      case 'b': push((f.data as string).length, 'M'); readPos += (f.data as string).length; break;
      default: break;
    }
  }
  if (readPos <= readLength) push(readLength - readPos + 1, 'M');
  // Sanity: reference consumption must match lengthOnRef; otherwise fall back to a plain match
  const refLen = ops.reduce((n, [len, op]) => n + ('MDN=X'.includes(op) ? len : 0), 0);
  if (lengthOnRef > 0 && refLen !== lengthOnRef) return `${readLength}M`;
  return ops.map(([len, op]) => `${len}${op}`).join('');
}

// ======================== Exon usage (port of backend/exon_usage.py) ========================

export type StrandnessCall = 'firststrand' | 'secondstrand' | 'unstranded' | 'unknown';

/** Transcript strand implied by a read under fr-firststrand (dUTP): read 2 on the transcript strand, read 1 / single-end opposite. */
export function txStrandFirststrand(r: RawRead): number {
  const mapped = (r.flags & FLAG_REVERSE) ? -1 : 1;
  return (r.flags & FLAG_PAIRED) && (r.flags & FLAG_READ2) ? mapped : -mapped;
}

export function detectStrandness(reads: RawRead[], geneStrand: number): { strandness: StrandnessCall; fraction: number | null } {
  const sample = reads.slice(0, 3000);
  if (sample.length < 100) return { strandness: 'unknown', fraction: null };
  const agree = sample.filter(r => txStrandFirststrand(r) === geneStrand).length;
  const f = agree / sample.length;
  return { strandness: f >= 0.9 ? 'firststrand' : f <= 0.1 ? 'secondstrand' : 'unstranded', fraction: f };
}

export function strandKeeper(strandness: StrandnessCall, geneStrand: number): (r: RawRead) => boolean {
  if (strandness === 'firststrand') return r => txStrandFirststrand(r) === geneStrand;
  if (strandness === 'secondstrand') return r => txStrandFirststrand(r) !== geneStrand;
  return () => true;
}

/** Median / mean depth of aligned bases over [start, end) and the number of reads with a block on it. */
export function exonDepth(reads: AlignedRead[], start: number, end: number): { median: number; mean: number; reads: number } {
  const depth = depthArray(reads, start, end);
  let n = 0;
  for (const r of reads) if (r.b.some(([bs, be]) => Math.min(be, end) > Math.max(bs, start))) n++;
  const sorted = Array.from(depth).sort((a, b) => a - b);
  const m = sorted.length >> 1;
  const med = sorted.length ? (sorted.length % 2 ? sorted[m] : (sorted[m - 1] + sorted[m]) / 2) : 0;
  const mean = sorted.length ? sorted.reduce((a, b) => a + b, 0) / sorted.length : 0;
  return { median: med, mean, reads: n };
}

// ======================== Structural evidence (genomic libraries) ========================

export const SV_MIN_DELETION = 50, SV_MIN_CLIP = 20, SV_MIN_SUPPORT = 3;
/** Mates facing away from each other (← →) at least this far apart (and twice the median insert) are duplication-type. */
export const OUTWARD_MIN_BP = 300;
const FLAG_PROPER = 2, FLAG_MATE_UNMAPPED = 8, FLAG_MATE_REVERSE = 32, FLAG_SUPPLEMENTARY = 2048;

const cigarRefLen = (cigar: string) => parseCigar(cigar).reduce((n, [len, op]) => n + ('MDN=X'.includes(op) ? len : 0), 0);
const sameChrom = (a: string, b: string) => a.replace(/^chr/i, '').toUpperCase() === b.replace(/^chr/i, '').toUpperCase();

/** One aligned part of a read: reference interval (0-based half-open), strand, and the read (query) interval it covers. */
interface Segment { chrom: string; start: number; end: number; rev: boolean; qs: number; qe: number }
/** Segment geometry from a CIGAR: leading and trailing clips (S or H) give the query interval, reversed on the minus strand. */
function segment(chrom: string, start: number, rev: boolean, cigar: string): Segment {
  const ops = parseCigar(cigar);
  let refLen = 0, qLen = 0, lead = 0, trail = 0, seenAligned = false;
  for (const [len, op] of ops) {
    if (op === 'S' || op === 'H') { if (!seenAligned) lead += len; else trail += len; qLen += len; continue; }
    if ('MI=X'.includes(op)) { qLen += len; seenAligned = true; }
    if ('MDN=X'.includes(op)) { refLen += len; seenAligned = true; }
  }
  const qs = rev ? trail : lead, qe = qLen - (rev ? lead : trail);
  return { chrom, start, end: start + refLen, rev, qs, qe };
}

/** Majority consensus of clipped sequences anchored at the breakpoint: `right` clips start there, `left` clips end there. */
export function clipConsensus(seqs: string[], side: 'left' | 'right'): { seq: string; depth: number[] } {
  const rows = side === 'left' ? seqs.map(x => x.split('').reverse().join('')) : seqs;
  const need = Math.min(2, rows.length);
  const out: string[] = [], depth: number[] = [];
  for (let k = 0; ; k++) {
    const counts: Record<string, number> = {}; let covering = 0;
    for (const r of rows) if (r.length > k) { covering++; counts[r[k]] = (counts[r[k]] ?? 0) + 1; }
    if (covering < need || covering === 0) break;
    const best = Object.entries(counts).sort((a, b) => b[1] - a[1])[0];
    out.push(best[1] / covering >= 0.6 ? best[0] : 'N'); depth.push(covering);
  }
  return side === 'left' ? { seq: out.reverse().join(''), depth: depth.reverse() } : { seq: out.join(''), depth };
}

/** Longest window the clipped consensus is checked over once its seed is found; mismatches allowed: one per 20 bases. */
const REALIGN_CHECK_BP = 60;
/**
 * Places a clipped consensus on the reference of the window: its 20 bases next to the breakpoint (the first ones of
 * a right clip, the last ones of a left clip) must occur exactly once over both strands, and the consensus must
 * then agree with the reference over up to 60 bases with at most one mismatch per 20. Returns the 0-based interval
 * covered on the forward strand, the strand, and the length matched; null when not placed.
 */
export function placeClip(consensus: string, side: 'left' | 'right', ref: { start: number; seq: string; rc?: string }): { start: number; end: number; strand: '+' | '-'; matched: number } | null {
  const seq = consensus.replace(/N+$/, '').replace(/^N+/, '');
  if (seq.length < SV_MIN_CLIP || seq.includes('N')) return null;
  const check = Math.min(seq.length, REALIGN_CHECK_BP);
  // the part checked is the one next to the breakpoint
  const probe = side === 'right' ? seq.slice(0, check) : seq.slice(seq.length - check);
  const seed = side === 'right' ? probe.slice(0, SV_MIN_CLIP) : probe.slice(probe.length - SV_MIN_CLIP);
  const rc = ref.rc ?? (ref.rc = reverseComplement(ref.seq));
  const hits: { strand: '+' | '-'; at: number }[] = [];
  for (const [strand, text] of [['+', ref.seq], ['-', rc]] as const) {
    let from = 0;
    while (hits.length < 2) { const at = text.indexOf(seed, from); if (at < 0) break; hits.push({ strand, at }); from = at + 1; }
    if (hits.length >= 2) return null;
  }
  if (hits.length !== 1) return null;
  const h = hits[0];
  const text = h.strand === '+' ? ref.seq : rc;
  // the probe starts `probeOff` bases before the seed (0 for a right clip, check - 20 for a left one)
  const probeOff = side === 'right' ? 0 : check - SV_MIN_CLIP;
  const at = h.at - probeOff;
  if (at < 0 || at + check > text.length) return null;
  let mism = 0;
  for (let k = 0; k < check; k++) if (text[at + k] !== probe[k]) mism++;
  if (mism > Math.floor(check / 20)) return null;
  // forward-strand interval of the matched probe
  const fwdStart = h.strand === '+' ? at : text.length - (at + check);
  return { start: ref.start + fwdStart, end: ref.start + fwdStart + check, strand: h.strand, matched: check };
}

/** True when some read of the scan is soft-clipped by 20 bases or more without an SA tag and carries its sequence: the window is worth a reference for realignment. */
export function hasRealignableClips(reads: RawRead[]): boolean {
  return reads.some(r => !r.sa && r.seq && parseCigar(r.cigar).some(([len, op]) => op === 'S' && len >= SV_MIN_CLIP));
}
/** True when some read is soft-clipped by RESCUE_MIN_CLIP bases or more without an SA tag and carries its sequence: a known breakpoint could rescue it. */
export function hasRescuableClips(reads: RawRead[]): boolean {
  return reads.some(r => !r.sa && r.seq && parseCigar(r.cigar).some(([len, op]) => op === 'S' && len >= RESCUE_MIN_CLIP));
}

/** The breakpoints of a sample's evidence, as rescue targets: every arc of the four kinds. */
export function breakpointsOf(sv: StructuralEvidence): Breakpoint[] {
  return [
    ...sv.deletions.map(j => ({ start: j.start, end: j.end, kind: 'deletion' as const })),
    ...sv.splits.map(j => ({ start: j.start, end: j.end, kind: 'split' as const })),
    ...(sv.duplications ?? []).map(j => ({ start: j.start, end: j.end, kind: 'duplication' as const })),
    ...(sv.inversions ?? []).map(j => ({ start: j.start, end: j.end, kind: 'inversion' as const })),
  ];
}
/** Shortest soft clip a known breakpoint can rescue: 8 bases with one candidate at the position, 12 with several. */
export const RESCUE_MIN_CLIP = 8, RESCUE_MIN_CLIP_MULTI = 12;
/** Arc ends are rounded to 5 bp: a clip within this distance of an end belongs to it. */
export const RESCUE_TOLERANCE_BP = 2;
const RESCUE_TOLERANCE = RESCUE_TOLERANCE_BP;

/** One soft- or hard-clipped read end without SA tag: where it is clipped, on which side, and the clipped bases (empty for a hard clip). */
interface ClipEnd { pos: number; side: 'left' | 'right'; seq: string; hard: boolean }
/** The clipped ends of the reads (both sides when both are clipped), soft clips of RESCUE_MIN_CLIP bases or more with their sequence, hard clips of SV_MIN_CLIP or more without. */
export function clipEnds(reads: RawRead[]): ClipEnd[] {
  const out: ClipEnd[] = [];
  for (const r of reads) {
    if (r.sa) continue;
    const ops = parseCigar(r.cigar);
    let pos = r.start, leftClip = 0, rightClip = 0, leftHard = 0, rightHard = 0, seenAligned = false, qLen = 0;
    for (const [len, op] of ops) {
      if (op === 'S') { if (!seenAligned) leftClip += len; else rightClip += len; qLen += len; }
      else if (op === 'H') { if (!seenAligned) leftHard += len; else rightHard += len; }
      else if ('MI=X'.includes(op)) { qLen += len; seenAligned = true; }
      if ('MDN=X'.includes(op)) { pos += len; seenAligned = true; }
    }
    const seqOk = !!r.seq && r.seq.length === qLen;
    if (leftClip >= RESCUE_MIN_CLIP && seqOk) out.push({ pos: r.start, side: 'left', seq: r.seq.substring(0, leftClip), hard: false });
    else if (!leftClip && leftHard >= SV_MIN_CLIP) out.push({ pos: r.start, side: 'left', seq: '', hard: true });
    if (rightClip >= RESCUE_MIN_CLIP && seqOk) out.push({ pos, side: 'right', seq: r.seq.substring(qLen - rightClip), hard: false });
    else if (!rightClip && rightHard >= SV_MIN_CLIP) out.push({ pos, side: 'right', seq: '', hard: true });
  }
  return out;
}

/** Bases of the forward reference from `from` (0-based) over `len`, or null when outside the reference window. */
const refSlice = (ref: { start: number; seq: string }, from: number, len: number): string | null => {
  if (from < ref.start || from + len > ref.start + ref.seq.length || len <= 0) return null;
  return ref.seq.substring(from - ref.start, from - ref.start + len);
};
/** A clipped sequence against the reference bases it should equal: no mismatch under 20 bases, one per 20 above. */
const clipMatches = (clip: string, target: string | null): boolean => {
  if (!target || target.length !== clip.length) return false;
  let mism = 0;
  for (let k = 0; k < clip.length; k++) if (clip[k] !== target[k]) { mism++; if (mism > Math.floor(clip.length / 20)) return false; }
  return true;
};
/**
 * What the clipped bases of a read clipped at `pos` on `side` should read if the read crossed breakpoint `b`,
 * for the end of `b` the clip sits at (`atStart`: the clip is at b.start, else at b.end). The offset of the clip from
 * the rounded end is carried to the other end (breakpoints shift together along a microhomology). Null when this side
 * of the read cannot cross `b` that way.
 *  deletion-type / CIGAR deletion (lo, hi): right clip at lo reads ref[hi…]; left clip at hi reads ref[…lo]
 *  duplication (lo, hi): right clip at hi reads ref[lo…]; left clip at lo reads ref[…hi]
 *  inversion (lo, hi): right clip at lo reads rc(ref[…hi]); right clip at hi reads rc(ref[…lo]);
 *                      left clip at hi reads rc(ref[lo…]); left clip at lo reads rc(ref[hi…])
 */
function expectedClip(b: Breakpoint, atStart: boolean, pos: number, side: 'left' | 'right', len: number, ref: { start: number; seq: string }): string | null {
  const d = pos - (atStart ? b.start : b.end);
  const lo = b.start + d, hi = b.end + d;
  const fwd = (from: number) => refSlice(ref, from, len);
  const rc = (from: number) => { const x = refSlice(ref, from, len); return x == null ? null : reverseComplement(x); };
  if (b.kind === 'deletion' || b.kind === 'split') {
    if (atStart && side === 'right') return fwd(hi);
    if (!atStart && side === 'left') return fwd(lo - len);
    return null;
  }
  if (b.kind === 'duplication') {
    if (!atStart && side === 'right') return fwd(lo);
    if (atStart && side === 'left') return fwd(hi - len);
    return null;
  }
  // inversion
  if (atStart && side === 'right') return rc(hi - len);
  if (!atStart && side === 'right') return rc(lo - len);
  if (!atStart && side === 'left') return rc(lo);
  if (atStart && side === 'left') return rc(hi);
  return null;
}
/**
 * Rescues clipped read ends at known breakpoints: a soft clip whose bases match the reference at the other end of a
 * breakpoint whose end it sits at (within the 5 bp rounding), a hard clip by position alone. A clip that fits several
 * breakpoints is dropped. Returns the counts per breakpoint, keyed `${kind}:${start}-${end}`.
 */
export function rescueClipEnds(ends: ClipEnd[], breakpoints: Breakpoint[], ref: { start: number; seq: string } | null): Map<string, { count: number; hard: number; /** rescued ends that were members of a clip cluster (hard, or soft of SV_MIN_CLIP bases or more) */ long: number }> {
  const out = new Map<string, { count: number; hard: number; long: number }>();
  if (!breakpoints.length) return out;
  const bump = (b: Breakpoint, e: ClipEnd) => { const k = `${b.kind}:${b.start}-${b.end}`; const long = e.hard || e.seq.length >= SV_MIN_CLIP ? 1 : 0; const x = out.get(k); if (x) { x.count++; x.long += long; if (e.hard) x.hard++; } else out.set(k, { count: 1, hard: e.hard ? 1 : 0, long }); };
  for (const e of ends) {
    // the breakpoints an end of which the clip sits at
    const near = breakpoints.flatMap(b => {
      const tol = b.kind === 'deletion' ? 0 : RESCUE_TOLERANCE;
      const at: { b: Breakpoint; atStart: boolean }[] = [];
      if (Math.abs(e.pos - b.start) <= tol) at.push({ b, atStart: true });
      if (Math.abs(e.pos - b.end) <= tol) at.push({ b, atStart: false });
      return at;
    }).filter(x => e.hard || expectedClip(x.b, x.atStart, e.pos, e.side, 1, ref ?? { start: 0, seq: 'N' }) !== undefined);
    if (!near.length) continue;
    if (e.hard) {
      // no sequence: attached by position when exactly one breakpoint has an end there, on a side that can hold a clip
      const fits = near.filter(x => (x.b.kind === 'deletion' || x.b.kind === 'split') ? (x.atStart ? e.side === 'right' : e.side === 'left')
        : x.b.kind === 'duplication' ? (x.atStart ? e.side === 'left' : e.side === 'right') : true);
      if (fits.length === 1) bump(fits[0].b, e);
      continue;
    }
    if (!ref) continue;
    if (e.seq.length < (near.length > 1 ? RESCUE_MIN_CLIP_MULTI : RESCUE_MIN_CLIP)) continue;
    const hits = near.filter(x => clipMatches(e.seq, expectedClip(x.b, x.atStart, e.pos, e.side, e.seq.length, ref)));
    if (hits.length === 1) bump(hits[0].b, e);
  }
  return out;
}

/**
 * Deletions inside reads, split reads, soft-clip clusters and discordant pairs of a window (0-based half-open),
 * from the light records of the scan. `rate` scales every count back when the window was sampled; support
 * thresholds apply to the scaled counts. `chrom` is the window's chromosome as named in the file.
 *
 * Split reads are read as chains: every part of a read (the record itself plus the parts its SA tag lists,
 * primary and supplementary alike) is placed on the read by its clips and the parts are ordered along the
 * read; each pair of adjacent parts is one breakpoint, classified by where the read continues: further on
 * the same strand (deletion-type), backwards (duplication-type), on the other strand (inversion), on another
 * chromosome (translocation), or after an unaligned stretch of the read (insertion). Each read counts once,
 * whichever of its parts fall in the window.
 *
 * Clipped reads without an SA tag form clip clusters. With the reference of the window (`ref`), the clipped
 * consensus of each cluster is placed by realignment (`placeClip`) and the cluster joins the arc a split read
 * would give between its clip position and the placed sequence; hard-clipped records without SA tag, which carry
 * no sequence, count in the cluster at their clip position and follow it into the arc.
 */
/**
 * `insertMedian`, when given, is the median |template length| of the proper pairs of the window, for callers that
 * pass only the records carrying structural evidence (the others would count for that median alone).
 */
export function structuralEvidence(reads: RawRead[], chrom: string, start: number, end: number, rate: number, ref?: { start: number; seq: string } | null, insertMedian?: number | null): StructuralEvidence {
  // every arc before merging, with its evidence units by source and the names of the reads behind them
  const dels = new Map<string, SvMember>(), splits = new Map<string, SvMember>(), dups = new Map<string, SvMember>(), invs = new Map<string, SvMember>();
  const discKind = { deletion: new Map<string, SvMember>(), duplication: new Map<string, SvMember>(), inversion: new Map<string, SvMember>() };
  /**
   * The breakpoints each discordant pair points to, from its reads (the bins only group the pairs): mates facing away
   * (duplication) lie inside the duplicated stretch, so its ends are at or before the left mates' starts and at or after
   * the right mates' ends; mates facing each other (deletion) lie outside, the ends at or after the left mates' ends and
   * at or before the right mates' starts; same-strand mates (inversion) end before (→ →) or start after (← ←) the
   * breakpoints.
   */
  const pairEnds = new Map<SvMember, [number, number][]>();
  const clips = new Map<string, ClipCluster>(), elsewhere = new Map<string, ElsewhereLink>();
  const clipSeqs = new Map<string, string[]>();
  const insertions = new Map<number, { pos: number; len: number; count: number }>();
  const add = (m: Map<string, SvMember>, s: number, e: number, n = 1, src = 'other', name?: string) => {
    const k = `${s}-${e}`;
    let j = m.get(k);
    if (!j) m.set(k, j = { start: s, end: e, count: 0, src: {}, names: new Set(), anon: 0 });
    j.count += n; j.src[src] = (j.src[src] ?? 0) + n;
    if (name) j.names.add(name); else j.anon += n;
  };
  const far = (m: Map<string, ElsewhereLink>, kind: 'split' | 'pair', pos: number, target: string) => { const k = `${kind}${target}@${pos}`; const x = m.get(k); if (x) x.count++; else m.set(k, { kind, pos, chrom: target, count: 1 }); };
  const r5 = (x: number) => Math.round(x / 5) * 5;
  const inWindow = (a: number, b: number) => b > start && a < end;
  let median = insertMedian;
  if (median === undefined) {
    const inserts: number[] = [];
    for (const r of reads) if (r.flags & FLAG_PAIRED && r.flags & FLAG_PROPER && r.tlen) inserts.push(Math.abs(r.tlen));
    median = inserts.length ? inserts.sort((a, b) => a - b)[inserts.length >> 1] : null;
  }
  const farInsert = median ? Math.max(5 * median, 1000) : Infinity;
  // mates facing away from each other (← →) further apart than this: the junction of a tandem duplication. Closer, they
  // are short fragments whose reads ran past each other (adapter read-through), not an event
  const outwardMin = Math.max(OUTWARD_MIN_BP, 2 * (median ?? 0));
  // primary records of each pair among `reads`: a pair is counted from its leftmost mate, or from the other one when the
  // leftmost is not among them (a mate beyond the window: the reads of a duplication's far side, 4 kb away, were all lost)
  const primaries = new Map<string, number>();
  // where each primary record's alignment ends, by name and start: a mate's real end, when its record is here (a mate
  // clipped at the junction ends there, not a read length after its start)
  const alignedEnd = new Map<string, number>();
  for (const r of reads) if (r.name && r.flags & FLAG_PAIRED && !(r.flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY))) {
    primaries.set(r.name, (primaries.get(r.name) ?? 0) + 1);
    let e = r.start;
    for (const [len, op] of parseCigar(r.cigar)) if ('MDN=X'.includes(op)) e += len;
    alignedEnd.set(`${r.name}:${r.start}`, e);
  }
  /** one record per split read, the primary when it is in the window */
  const chains = new Map<string, RawRead>();
  // Deleted and skipped bases (CIGAR D, N) inside the primary records of each pair, by read name. The template length
  // runs from one mate's outer end to the other's, across them: a pair whose gap is a deletion one of its reads already
  // carries has a long TLEN but is not a discordant pair (it is the same event as that D, counted with the deletions).
  const innerGap = (cigar: string) => { let g = 0; for (const [len, op] of parseCigar(cigar)) if (op === 'D' || op === 'N') g += len; return g; };
  const pairGap = new Map<string, number>();
  for (const r of reads) {
    if (!r.name || r.flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) || !(r.flags & FLAG_PAIRED)) continue;
    const g = innerGap(r.cigar);
    if (g) pairGap.set(r.name, (pairGap.get(r.name) ?? 0) + g);
  }
  const cluster = (key: string, pos: number, side: 'left' | 'right', hard: boolean, seq?: string) => {
    const c = clips.get(key);
    if (c) { c.count++; if (hard) c.hard = (c.hard ?? 0) + 1; } else clips.set(key, { pos, side, count: 1, hard: hard ? 1 : 0 });
    if (seq) { const l = clipSeqs.get(key); if (l) l.push(seq); else clipSeqs.set(key, [seq]); }
  };
  for (const r of reads) {
    const ops = parseCigar(r.cigar);
    let pos = r.start, leftClip = 0, rightClip = 0, leftHard = 0, rightHard = 0, seenAligned = false, qLen = 0;
    ops.forEach(([len, op]) => {
      if (op === 'S') { if (!seenAligned) leftClip += len; else rightClip += len; qLen += len; }
      else if (op === 'H') { if (!seenAligned) leftHard += len; else rightHard += len; }
      else if ('MI=X'.includes(op)) { qLen += len; seenAligned = true; }
      if (op === 'D' && len >= SV_MIN_DELETION && pos + len > start && pos < end) add(dels, pos, pos + len, 1, 'cigar', r.name);
      if ('MDN=X'.includes(op)) { pos += len; seenAligned = true; }
    });
    const alnEnd = pos;
    // clip clusters: reads whose clipped part is placed elsewhere (SA tag) are split reads, drawn from their chain instead;
    // soft clips bring their sequence when the scan decoded it, hard clips without SA tag only their position
    if (!r.sa) {
      const seqOk = r.seq && r.seq.length === qLen;
      if (leftClip >= SV_MIN_CLIP && r.start >= start && r.start < end) cluster(`L${r.start}`, r.start, 'left', false, seqOk ? r.seq.substring(0, leftClip) : undefined);
      else if (leftHard >= SV_MIN_CLIP && r.start >= start && r.start < end) cluster(`L${r.start}`, r.start, 'left', true);
      if (rightClip >= SV_MIN_CLIP && alnEnd > start && alnEnd <= end) cluster(`R${alnEnd}`, alnEnd, 'right', false, seqOk ? r.seq.substring(qLen - rightClip) : undefined);
      else if (rightHard >= SV_MIN_CLIP && alnEnd > start && alnEnd <= end) cluster(`R${alnEnd}`, alnEnd, 'right', true);
    }
    if (r.sa && r.name) { const prev = chains.get(r.name); if (!prev || ((prev.flags & FLAG_SUPPLEMENTARY) && !(r.flags & FLAG_SUPPLEMENTARY))) chains.set(r.name, r); }
    // discordant pairs, counted once per pair: from the leftmost mate (the first of the pair when both start at the same
    // base), and from primary records only (a supplementary record repeats its primary's mate fields)
    if (r.flags & FLAG_PAIRED && !(r.flags & (FLAG_MATE_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) && r.mateChrom != null && r.matePos != null) {
      if (!sameChrom(r.mateChrom, chrom)) far(elsewhere, 'pair', Math.floor(r.start / 500) * 500, r.mateChrom);   // mates elsewhere never share a start: binned like the discordant pairs
      else if (r.start < r.matePos || (r.start === r.matePos && !(r.flags & FLAG_READ2)) || (r.name && primaries.get(r.name) === 1)) {
        const rev = (r.flags & FLAG_REVERSE) !== 0, mateRev = (r.flags & FLAG_MATE_REVERSE) !== 0;
        const leftmost = r.start <= r.matePos;
        // the pair's extent: from the leftmost start to the rightmost end. From the positions, the mate's end taken as
        // its start plus this read's aligned length (its own deletions and introns left out): TLEN is no guide across an
        // event (BWA gave pairs of a 4 kb duplication, mates 3.9 kb apart, a TLEN of 110 and the proper-pair flag)
        const readLen = alnEnd - r.start - innerGap(r.cigar);
        const lo = Math.min(r.start, r.matePos);
        const hi = r.matePos >= r.start ? Math.max(alnEnd, r.matePos + readLen) : alnEnd;
        const span = hi - lo - (r.name ? pairGap.get(r.name) ?? 0 : innerGap(r.cigar));
        // orientation of the leftmost mate and of the other: → ← normal, ← → outward (duplication), same strand (inversion)
        const leftRev = leftmost ? rev : mateRev, rightRev = leftmost ? mateRev : rev;
        const kind: 'deletion' | 'duplication' | 'inversion' | null = rev === mateRev ? 'inversion'
          : leftRev && !rightRev ? (Math.abs(r.matePos - r.start) > outwardMin ? 'duplication' : null)
          : span > farInsert ? 'deletion' : null;
        if (kind) {
          const a = Math.floor(lo / 500) * 500, b = Math.ceil(hi / 500) * 500;
          if (b > a) {
            add(discKind[kind], a, b, 1, 'pair', r.name);
            // where the breakpoints lie from this pair's reads (see pairEnds): the arc is drawn there, not at the bins
            const mateEnd = (r.name ? alignedEnd.get(`${r.name}:${r.matePos}`) : undefined) ?? r.matePos + readLen;
            const leftEnd = leftmost ? alnEnd : mateEnd, rightStart = Math.max(r.start, r.matePos), rightEnd = leftmost ? mateEnd : alnEnd;
            const est: [number, number] = kind === 'duplication' ? [lo, rightEnd] : kind === 'deletion' ? [leftEnd, rightStart]
              : !rev ? [leftEnd, rightEnd] : [lo, rightStart];
            const m = discKind[kind].get(`${a}-${b}`)!;
            const e = pairEnds.get(m);
            if (e) e.push(est); else pairEnds.set(m, [est]);
          }
        }
      }
    }
  }
  /** One breakpoint between two adjacent parts of a read: classified and counted `n` times; returns the arc it joined, if any. */
  const breakpoint = (a: Segment, b: Segment, n = 1, src = 'split', name?: string): RealignedClip['arc'] | null => {
    const aOut = a.rev ? a.start : a.end;          // where the read leaves part a on the reference
    const bIn = b.rev ? b.end : b.start;           // where it enters part b
    const aHere = sameChrom(a.chrom, chrom), bHere = sameChrom(b.chrom, chrom);
    if (!aHere && !bHere) return null;
    if (!sameChrom(a.chrom, b.chrom)) { if (aHere && aOut >= start && aOut < end) far(elsewhere, 'split', aOut, b.chrom); else if (bHere && bIn >= start && bIn < end) far(elsewhere, 'split', bIn, a.chrom); return null; }
    const lo = r5(Math.min(aOut, bIn)), hi = r5(Math.max(aOut, bIn));
    if (a.rev !== b.rev) { if (inWindow(lo, hi) && hi > lo) { add(invs, lo, hi, n, src === 'split' ? (a.rev ? '-/+' : '+/-') : src, name); return { start: lo, end: hi, kind: 'inversion' }; } return null; }
    const refGap = a.rev ? a.start - b.end : b.start - a.end;
    const qGap = b.qs - a.qe;
    if (refGap >= SV_MIN_DELETION) { if (inWindow(lo, hi)) { add(splits, lo, hi, n, src, name); return { start: lo, end: hi, kind: 'split' }; } }
    else if (refGap <= -SV_MIN_DELETION) { if (inWindow(lo, hi) && hi > lo) { add(dups, lo, hi, n, src, name); return { start: lo, end: hi, kind: 'duplication' }; } }
    else if (qGap >= SV_MIN_DELETION && aOut >= start && aOut < end) { const k = r5(aOut); const x = insertions.get(k); if (x) { x.count += n; x.len = Math.round((x.len * (x.count - n) + qGap * n) / x.count); } else insertions.set(k, { pos: k, len: qGap, count: n }); }
    return null;
  };
  // split reads: the chain of every part of the read, ordered along the read
  for (const r of chains.values()) {
    const segs: Segment[] = [segment(chrom, r.start, (r.flags & FLAG_REVERSE) !== 0, r.cigar)];
    for (const part of (r.sa ?? '').split(';')) {
      const f = part.split(',');
      if (f.length < 4) continue;
      const saStart = parseInt(f[1]) - 1;
      if (!Number.isFinite(saStart)) continue;
      segs.push(segment(f[0], saStart, f[2] === '-', f[3]));
    }
    segs.sort((a, b) => a.qs - b.qs);
    for (let i = 0; i + 1 < segs.length; i++) breakpoint(segs[i], segs[i + 1], 1, 'split', r.name);
  }
  // clip clusters placed by realignment of their clipped consensus: the cluster becomes a breakpoint like a split read's
  const realigned: RealignedClip[] = [];
  const placedKeys = new Set<string>();
  if (ref && ref.seq) {
    const refWin = { start: ref.start, seq: ref.seq, rc: undefined as string | undefined };
    for (const [key, c] of [...clips]) {
      const seqs = clipSeqs.get(key);
      if (!seqs || c.count * rate < SV_MIN_SUPPORT) continue;
      const cons = clipConsensus(seqs, c.side).seq;
      const hit = placeClip(cons, c.side, refWin);
      if (!hit) continue;
      const L = hit.matched;
      // the aligned side as a one-base anchor at the clip position, the placed clip as the other part, in read order
      const anchor: Segment = c.side === 'right' ? { chrom, start: c.pos - 1, end: c.pos, rev: false, qs: 0, qe: 1 } : { chrom, start: c.pos, end: c.pos + 1, rev: false, qs: L, qe: L + 1 };
      const placed: Segment = { chrom, start: hit.start, end: hit.end, rev: hit.strand === '-', qs: c.side === 'right' ? 1 : 0, qe: c.side === 'right' ? 1 + L : L };
      const arc = c.side === 'right' ? breakpoint(anchor, placed, c.count, 'clip') : breakpoint(placed, anchor, c.count, 'clip');
      if (!arc) continue;
      realigned.push({ pos: c.pos, side: c.side, count: c.count * rate, hard: (c.hard ?? 0) * rate, target: hit.start, strand: hit.strand, matched: L, arc });
      placedKeys.add(key);
      clips.delete(key);
    }
  }
  // clipped reads rescued at the sample's own breakpoints (arcs from chains, deletions and placed clusters): their clipped
  // bases must match the reference at the other end; reads of a placed cluster are already counted and stay out
  const rescued: RescuedClips[] = [];
  const own: Breakpoint[] = [
    ...[...dels.values()].map(j => ({ start: j.start, end: j.end, kind: 'deletion' as const })),
    ...[...splits.values()].map(j => ({ start: j.start, end: j.end, kind: 'split' as const })),
    ...[...dups.values()].map(j => ({ start: j.start, end: j.end, kind: 'duplication' as const })),
    ...[...invs.values()].map(j => ({ start: j.start, end: j.end, kind: 'inversion' as const })),
  ];
  if (own.length) {
    // members of a placed cluster (hard, or soft clips of 20 bases or more at its position) are already counted with it
    const ends = clipEnds(reads).filter(e => e.pos >= start && e.pos <= end && !(placedKeys.has(`${e.side === 'left' ? 'L' : 'R'}${e.pos}`) && (e.hard || e.seq.length >= SV_MIN_CLIP)));
    for (const [k, n] of rescueClipEnds(ends, own, ref ?? null)) {
      const [kind, span] = k.split(':'); const [a, b] = span.split('-').map(Number);
      const m = kind === 'deletion' ? dels : kind === 'split' ? splits : kind === 'duplication' ? dups : invs;
      add(m, a, b, n.count, 'rescued');
      rescued.push({ start: a, end: b, kind: kind as Breakpoint['kind'], count: n.count * rate, hard: n.hard * rate, own: true });
      // rescued members of an unplaced cluster leave it (short clips were never in one)
      let left = n.long;
      for (const side of ['L', 'R'] as const) for (const pos of [a, b]) for (let d = -RESCUE_TOLERANCE; d <= RESCUE_TOLERANCE && left > 0; d++) {
        const c = clips.get(`${side}${pos + d}`);
        if (c) { const take = Math.min(c.count, left); c.count -= take; left -= take; if (!c.count) clips.delete(`${side}${pos + d}`); }
      }
    }
  }
  // one arc per event: arcs of a family whose breakpoints lie within the merge tolerance (svmerge.ts) are merged, the
  // deletions inside CIGARs with the deletion-type split reads and placed clips; the placed and rescued clips then
  // name the event they joined
  const eventOf = new Map<string, SvArc & { kind: Breakpoint['kind'] }>();
  const family = (kind: Breakpoint['kind'], maps: [string, Map<string, SvMember>][]) => {
    const members = maps.flatMap(([k, m]) => [...m.values()].map(x => ({ x, k })));
    const byMember = new Map(members.map(({ x, k }) => [x, k]));
    return clusterSv(members.map(({ x }) => x)).map(({ event, members: ms }) => {
      const ev = { ...event, count: event.count * rate };
      for (const m of ms) eventOf.set(`${byMember.get(m)}:${m.start}-${m.end}`, { ...ev, kind });
      return ev;
    }).sort((a, b) => a.start - b.start || a.end - b.end);
  };
  const deletionEvents = family('deletion', [['deletion', dels], ['split', splits]]);
  const duplicationEvents = family('duplication', [['duplication', dups]]);
  const inversionEvents = family('inversion', [['inversion', invs]]);
  for (const x of realigned) { const ev = eventOf.get(`${x.arc.kind}:${x.arc.start}-${x.arc.end}`); if (ev) x.arc = { start: ev.start, end: ev.end, kind: ev.kind }; }
  for (const x of rescued) { const ev = eventOf.get(`${x.kind}:${x.start}-${x.end}`); if (ev) { x.start = ev.start; x.end = ev.end; x.kind = ev.kind; } }
  // the pairs of one event fall in neighbouring bins (their insert sizes vary): bins of one class within one bin of the
  // strongest at both ends join it, one arc per event
  const discordant: DiscordantArc[] = (Object.keys(discKind) as (keyof typeof discKind)[]).flatMap(kind => {
    const bins = [...discKind[kind].values()].sort((a, b) => b.count - a.count || a.start - b.start);
    const out: DiscordantArc[] = [], used = new Set<SvMember>();
    for (const x of bins) {
      if (used.has(x)) continue;
      const arc: DiscordantArc = { start: x.start, end: x.end, count: 0, kind };
      const ests: [number, number][] = [];
      for (const y of bins) if (!used.has(y) && Math.abs(y.start - x.start) <= 500 && Math.abs(y.end - x.end) <= 500) {
        used.add(y); arc.count += y.count * rate; arc.start = Math.min(arc.start, y.start); arc.end = Math.max(arc.end, y.end);
        ests.push(...(pairEnds.get(y) ?? []));
      }
      // the arc at the breakpoints the pairs point to: the outermost reads of a duplication's pairs, the innermost of a
      // deletion's (the breakpoint lies beyond them, within an insert size), the median of an inversion's
      if (ests.length) {
        const as = ests.map(e => e[0]).sort((p, q) => p - q), bs = ests.map(e => e[1]).sort((p, q) => p - q);
        const [s0, e0] = kind === 'duplication' ? [as[0], bs[bs.length - 1]] : kind === 'deletion' ? [as[as.length - 1], bs[0]] : [as[as.length >> 1], bs[bs.length >> 1]];
        if (e0 > s0) { arc.start = s0; arc.end = e0; }
      }
      out.push(arc);
    }
    return out;
  }).sort((a, b) => a.start - b.start || a.end - b.end);
  return {
    deletions: deletionEvents, splits: [], duplications: duplicationEvents, inversions: inversionEvents, discordant,
    insertions: [...insertions.values()].map(x => ({ ...x, count: x.count * rate })).filter(x => x.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    elsewhere: [...elsewhere.values()].map(x => ({ ...x, count: x.count * rate })).filter(x => x.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    clips: [...clips.values()].map(c => ({ ...c, count: c.count * rate, hard: (c.hard ?? 0) * rate })).filter(c => c.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    realigned: realigned.length ? realigned.sort((a, b) => a.pos - b.pos) : undefined,
    rescued: rescued.length ? rescued.sort((a, b) => a.start - b.start || a.end - b.end) : undefined,
    insertMedian: median, reads: reads.length,
  };
}
