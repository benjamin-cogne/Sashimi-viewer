/**
 * Browser-side port of backend/bam_reader.py: turn BAM/CRAM records into the compact
 * read encoding the viewer consumes, plus coverage runs and junction counts.
 * All coordinates are 0-based half-open.
 */
import type { AlignedRead, CoverageRun, JunctionArc } from '../components/sashimi/types';

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

/** Same filter as pysam count_coverage(read_callback='all'). */
export function keepRead(r: RawRead): boolean {
  return (r.flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_QCFAIL | FLAG_DUP)) === 0;
}
export function isUnique(r: RawRead): boolean {
  return r.nh != null ? r.nh === 1 : r.mapq >= 30;
}

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
