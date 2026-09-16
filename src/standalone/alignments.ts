/**
 * Browser-side port of backend/bam_reader.py: turn BAM/CRAM records into the compact
 * read encoding the viewer consumes, plus coverage runs and junction counts.
 * All coordinates are 0-based half-open.
 */
import type { BoundarySpanning, AlignedRead, CoverageRun, JunctionArc, StructuralEvidence, ClipCluster, ElsewhereLink } from '../components/sashimi/types';

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
  const clips: [number, number] = [0, 0];
  let rpos = r.start, qpos = 0, first = true;
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
      rpos += len; qpos += len;
    } else if (op === 'I') { ins.push([rpos, len]); qpos += len; }
    else if (op === 'D') { dels.push([rpos, rpos + len]); rpos += len; }
    else if (op === 'N') { rpos += len; }
    else if (op === 'S') { clips[first ? 0 : 1] = len; qpos += len; }
    first = false;
  }
  if (!r.seq && r.mismatches && ref) {
    for (const m of r.mismatches) if (m[0] >= refStart && m[0] < refEnd && ref[m[0] - refStart] !== m[1] && ref[m[0] - refStart] !== 'N') mism.push(m);
  }
  return {
    n: r.name, s: r.start, e: rpos, r: (r.flags & FLAG_REVERSE) ? 1 : 0, q: r.mapq, f: r.flags, nh: r.nh,
    b: blocks, d: dels, i: ins,
    m: mism,
    c: clips,
  };
}

/** Run-length coverage of [start, end) from the aligned blocks of encoded reads. */
export function coverageRuns(reads: AlignedRead[], start: number, end: number): CoverageRun[] {
  const n = end - start;
  if (n <= 0) return [];
  const diff = new Int32Array(n + 1);
  for (const r of reads) for (const [bs, be] of r.b) {
    const lo = Math.max(bs, start), hi = Math.min(be, end);
    if (hi > lo) { diff[lo - start] += 1; diff[hi - start] -= 1; }
  }
  const runs: CoverageRun[] = [];
  let depth = 0, runStart = 0;
  for (let i = 0; i < n; i++) {
    const d = depth + diff[i];
    if (i > 0 && d !== depth) { runs.push({ start: start + runStart, end: start + i, depth }); runStart = i; }
    depth = d;
  }
  runs.push({ start: start + runStart, end, depth });
  return runs;
}

/** Anchors of a boundary-spanning read: aligned bases required on the exon side and on the intron side of the boundary. */
export const SPAN_EXON_ANCHOR = 6, SPAN_INTRON_ANCHOR = 10;

/**
 * Unspliced reads through exon–intron boundaries: one aligned block covering SPAN_EXON_ANCHOR bases
 * on the exon side and SPAN_INTRON_ANCHOR on the intron side. `intronStarts` are boundaries with the
 * intron to the right (exon ends), `intronEnds` with the intron to the left (exon starts).
 */
export function boundarySpanning(reads: AlignedRead[], intronStarts: Iterable<number>, intronEnds: Iterable<number>): BoundarySpanning {
  const starts = [...new Set(intronStarts)].sort((a, b) => a - b), ends = [...new Set(intronEnds)].sort((a, b) => a - b);
  const startCount = new Map<number, number>(starts.map(p => [p, 0])), endCount = new Map<number, number>(ends.map(p => [p, 0]));
  const lowerBound = (arr: number[], v: number) => { let lo = 0, hi = arr.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; } return lo; };
  for (const r of reads) for (const [bs, be] of r.b) {
    // intron start p: block must cover [p - exonAnchor, p + intronAnchor)
    for (let i = lowerBound(starts, bs + SPAN_EXON_ANCHOR); i < starts.length && starts[i] + SPAN_INTRON_ANCHOR <= be; i++) startCount.set(starts[i], (startCount.get(starts[i]) || 0) + 1);
    // intron end q: block must cover [q - intronAnchor, q + exonAnchor)
    for (let i = lowerBound(ends, bs + SPAN_INTRON_ANCHOR); i < ends.length && ends[i] + SPAN_EXON_ANCHOR <= be; i++) endCount.set(ends[i], (endCount.get(ends[i]) || 0) + 1);
  }
  return { intronStart: Object.fromEntries(startCount), intronEnd: Object.fromEntries(endCount) };
}

/** Junctions (CIGAR N gaps) overlapping [start, end). */
export function junctionCounts(reads: AlignedRead[], start: number, end: number): JunctionArc[] {
  const counts = new Map<string, JunctionArc>();
  for (const r of reads) {
    const dels = new Set(r.d.map(([s, e]) => `${s}-${e}`));
    for (let k = 0; k + 1 < r.b.length; k++) {
      const s = r.b[k][1], e = r.b[k + 1][0];
      if (e <= s || dels.has(`${s}-${e}`)) continue;
      if (e > start && s < end) {
        const key = `${s}-${e}`;
        const j = counts.get(key);
        if (j) j.count++; else counts.set(key, { start: s, end: e, count: 1 });
      }
    }
  }
  return [...counts.values()].sort((a, b) => a.start - b.start || a.end - b.end);
}

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
 */
export function structuralEvidence(reads: RawRead[], chrom: string, start: number, end: number, rate: number): StructuralEvidence {
  const dels = new Map<string, JunctionArc>(), splits = new Map<string, JunctionArc>(), dups = new Map<string, JunctionArc>(), invs = new Map<string, JunctionArc>(), disc = new Map<string, JunctionArc>();
  const clips = new Map<string, ClipCluster>(), elsewhere = new Map<string, ElsewhereLink>();
  const insertions = new Map<number, { pos: number; len: number; count: number }>();
  const add = (m: Map<string, JunctionArc>, s: number, e: number) => { const k = `${s}-${e}`; const j = m.get(k); if (j) j.count++; else m.set(k, { start: s, end: e, count: 1 }); };
  const far = (m: Map<string, ElsewhereLink>, kind: 'split' | 'pair', pos: number, target: string) => { const k = `${kind}${target}@${pos}`; const x = m.get(k); if (x) x.count++; else m.set(k, { kind, pos, chrom: target, count: 1 }); };
  const r5 = (x: number) => Math.round(x / 5) * 5;
  const inWindow = (a: number, b: number) => b > start && a < end;
  const inserts: number[] = [];
  for (const r of reads) if (r.flags & FLAG_PAIRED && r.flags & FLAG_PROPER && r.tlen) inserts.push(Math.abs(r.tlen));
  const median = inserts.length ? inserts.sort((a, b) => a - b)[inserts.length >> 1] : null;
  const farInsert = median ? Math.max(5 * median, 1000) : Infinity;
  /** one record per split read, the primary when it is in the window */
  const chains = new Map<string, RawRead>();
  for (const r of reads) {
    const ops = parseCigar(r.cigar);
    let pos = r.start, leftClip = 0, rightClip = 0;
    ops.forEach(([len, op], i) => {
      if (op === 'S') { if (i === 0 || (i === 1 && ops[0][1] === 'H')) leftClip = len; else rightClip = len; }
      if (op === 'D' && len >= SV_MIN_DELETION && pos + len > start && pos < end) add(dels, pos, pos + len);
      if ('MDN=X'.includes(op)) pos += len;
    });
    const alnEnd = pos;
    // clip clusters: reads whose clipped part is placed elsewhere (SA tag) are split reads, drawn from their chain instead
    if (!r.sa && leftClip >= SV_MIN_CLIP && r.start >= start && r.start < end) { const k = `L${r.start}`; const c = clips.get(k); if (c) c.count++; else clips.set(k, { pos: r.start, side: 'left', count: 1 }); }
    if (!r.sa && rightClip >= SV_MIN_CLIP && alnEnd > start && alnEnd <= end) { const k = `R${alnEnd}`; const c = clips.get(k); if (c) c.count++; else clips.set(k, { pos: alnEnd, side: 'right', count: 1 }); }
    if (r.sa && r.name) { const prev = chains.get(r.name); if (!prev || ((prev.flags & FLAG_SUPPLEMENTARY) && !(r.flags & FLAG_SUPPLEMENTARY))) chains.set(r.name, r); }
    // discordant pairs, counted once from the leftmost mate
    if (r.flags & FLAG_PAIRED && !(r.flags & FLAG_MATE_UNMAPPED) && r.mateChrom != null && r.matePos != null) {
      if (!sameChrom(r.mateChrom, chrom)) far(elsewhere, 'pair', r.start, r.mateChrom);
      else if (r.start <= r.matePos) {
        const sameStrand = ((r.flags & FLAG_REVERSE) !== 0) === ((r.flags & FLAG_MATE_REVERSE) !== 0);
        const span = Math.abs(r.tlen ?? (r.matePos + (alnEnd - r.start) - r.start));
        if (sameStrand || span > farInsert) {
          const a = Math.floor(r.start / 500) * 500, b = Math.ceil((r.matePos + (alnEnd - r.start)) / 500) * 500;
          if (b > a) add(disc, a, b);
        }
      }
    }
  }
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
    for (let i = 0; i + 1 < segs.length; i++) {
      const a = segs[i], b = segs[i + 1];
      const aOut = a.rev ? a.start : a.end;          // where the read leaves part a on the reference
      const bIn = b.rev ? b.end : b.start;           // where it enters part b
      const aHere = sameChrom(a.chrom, chrom), bHere = sameChrom(b.chrom, chrom);
      if (!aHere && !bHere) continue;
      if (!sameChrom(a.chrom, b.chrom)) { if (aHere && aOut >= start && aOut < end) far(elsewhere, 'split', aOut, b.chrom); else if (bHere && bIn >= start && bIn < end) far(elsewhere, 'split', bIn, a.chrom); continue; }
      const lo = r5(Math.min(aOut, bIn)), hi = r5(Math.max(aOut, bIn));
      if (a.rev !== b.rev) { if (inWindow(lo, hi) && hi > lo) add(invs, lo, hi); continue; }
      const refGap = a.rev ? a.start - b.end : b.start - a.end;
      const qGap = b.qs - a.qe;
      if (refGap >= SV_MIN_DELETION) { if (inWindow(lo, hi)) add(splits, lo, hi); }
      else if (refGap <= -SV_MIN_DELETION) { if (inWindow(lo, hi) && hi > lo) add(dups, lo, hi); }
      else if (qGap >= SV_MIN_DELETION && aOut >= start && aOut < end) { const k = r5(aOut); const x = insertions.get(k); if (x) { x.count++; x.len = Math.round((x.len * (x.count - 1) + qGap) / x.count); } else insertions.set(k, { pos: k, len: qGap, count: 1 }); }
    }
  }
  const scaled = (m: Map<string, JunctionArc>) => [...m.values()].map(j => ({ ...j, count: j.count * rate })).sort((a, b) => a.start - b.start || a.end - b.end);
  return {
    deletions: scaled(dels), splits: scaled(splits), duplications: scaled(dups), inversions: scaled(invs), discordant: scaled(disc),
    insertions: [...insertions.values()].map(x => ({ ...x, count: x.count * rate })).filter(x => x.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    elsewhere: [...elsewhere.values()].map(x => ({ ...x, count: x.count * rate })).filter(x => x.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    clips: [...clips.values()].map(c => ({ ...c, count: c.count * rate })).filter(c => c.count >= SV_MIN_SUPPORT).sort((a, b) => a.pos - b.pos),
    insertMedian: median, reads: reads.length,
  };
}
