/**
 * Compact storage of reads and coverage for exported pages.
 *
 * Three layers, each simple on its own:
 *  1. columns of small integers: every field of the reads is one array across all reads (starts as deltas from the
 *     previous read, everything inside a read relative to its start), written as variable-length integers;
 *  2. sections and a directory: the columns are grouped into named sections ("core", "reference", "pairs"), each
 *     compressed on its own, listed at the front with their sizes, so a reader decodes what it needs and skips
 *     what it does not know;
 *  3. deflate (the browser's own, with a JavaScript fallback for older browsers) and base64, so the bytes can live
 *     inside an HTML page.
 *
 * Reads are written in blocks of READS_PER_BLOCK, each block its own deflated column set with its genomic span in
 * the directory: a reader can decode only the blocks overlapping a window. The format is versioned; every field a
 * section holds is listed in docs/embedded-format.md.
 */
import { inflateSync, deflateSync } from 'fflate';
import type { AlignedRead, CoverageRun, JunctionArc } from '../components/sashimi/types';

export const COLUMNAR_VERSION = 5;
export const READS_PER_BLOCK = 1000;

// ======================== variable-length integers ========================

/** Signed → unsigned fold ("zigzag"): 0, -1, 1, -2, 2 … become 0, 1, 2, 3, 4 … so small negatives stay small. */
const zz = (v: number) => (v < 0 ? -2 * v - 1 : 2 * v);
const unzz = (u: number) => (u & 1 ? -(u + 1) / 2 : u / 2);

class Writer {
  private buf = new Uint8Array(1 << 16);
  private n = 0;
  private grow(extra: number) {
    if (this.n + extra <= this.buf.length) return;
    let size = this.buf.length * 2;
    while (size < this.n + extra) size *= 2;
    const b = new Uint8Array(size); b.set(this.buf); this.buf = b;
  }
  /** unsigned variable-length integer: 7 bits per byte, high bit = another byte follows */
  u(v: number) {
    this.grow(10);
    if (v < 0 || !Number.isFinite(v)) throw new Error(`unsigned varint out of range: ${v}`);
    while (v >= 128) { this.buf[this.n++] = (v % 128) | 128; v = Math.floor(v / 128); }
    this.buf[this.n++] = v;
  }
  /** signed variable-length integer */
  s(v: number) { this.u(zz(v)); }
  byte(b: number) { this.grow(1); this.buf[this.n++] = b & 255; }
  bytes(b: Uint8Array) { this.grow(b.length); this.buf.set(b, this.n); this.n += b.length; }
  /** length-prefixed UTF-8 string */
  str(s: string) { const b = new TextEncoder().encode(s); this.u(b.length); this.bytes(b); }
  done(): Uint8Array { return this.buf.slice(0, this.n); }
}

class Reader {
  n = 0;
  constructor(private buf: Uint8Array) {}
  u(): number {
    let v = 0, mul = 1, b: number;
    do {
      if (this.n >= this.buf.length) throw new Error('truncated column data');
      b = this.buf[this.n++]; v += (b & 127) * mul; mul *= 128;
    } while (b & 128);
    return v;
  }
  s(): number { return unzz(this.u()); }
  byte(): number { return this.buf[this.n++]; }
  bytes(len: number): Uint8Array { const b = this.buf.subarray(this.n, this.n + len); this.n += len; return b; }
  str(): string { const len = this.u(); return new TextDecoder().decode(this.bytes(len)); }
  get done() { return this.n >= this.buf.length; }
}

// ======================== deflate / inflate ========================

const hasStreams = typeof CompressionStream === 'function' && typeof DecompressionStream === 'function';

async function throughStream(data: Uint8Array, stream: { readable: ReadableStream<Uint8Array>; writable: WritableStream<BufferSource> }): Promise<Uint8Array> {
  const writer = stream.writable.getWriter();
  void writer.write(data as unknown as BufferSource).then(() => writer.close()).catch(() => { /* surfaced by the reader below */ });
  const reader = stream.readable.getReader();
  const chunks: Uint8Array[] = []; let total = 0;
  for (;;) { const { value, done } = await reader.read(); if (done) break; chunks.push(value); total += value.length; }
  const out = new Uint8Array(total); let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

/** Raw deflate: the browser's compressor when it has one (asynchronous, off the main thread), else fflate. */
export async function deflate(data: Uint8Array): Promise<Uint8Array> {
  if (hasStreams) {
    try { return await throughStream(data, new CompressionStream('deflate-raw')); } catch { /* fall back */ }
  }
  return deflateSync(data, { level: 6 });
}
/** Raw inflate, same choice of engine. */
export async function inflate(data: Uint8Array): Promise<Uint8Array> {
  if (hasStreams) {
    try { return await throughStream(data, new DecompressionStream('deflate-raw')); } catch { /* fall back */ }
  }
  return inflateSync(data);
}
/** True when the page can decode without the JavaScript fallback (which is slower but always available). */
export const nativeCodec = hasStreams;

// ======================== base64 ========================

const NodeBuffer: any = (globalThis as any).Buffer;
export function toBase64(bytes: Uint8Array): string {
  if (NodeBuffer) return NodeBuffer.from(bytes).toString('base64');
  let s = '';
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000) as unknown as number[]);
  return btoa(s);
}
export function fromBase64(text: string): Uint8Array {
  if (NodeBuffer) return new Uint8Array(NodeBuffer.from(text, 'base64'));
  const s = atob(text);
  const out = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}

// ======================== sections ========================

/** One compressed section: its name, what it holds and where it sits in the byte stream. */
export interface SectionEntry {
  name: string;
  /** deflated bytes */
  bytes: number;
  /** for read blocks: index of the block and the genomic span [start, end) of the reads it holds */
  block?: number; start?: number; end?: number;
  /** number of reads in a block */
  n?: number;
}
export interface Directory {
  v: number;
  /** what the stream holds */
  kind: 'reads' | 'coverage';
  sections: SectionEntry[];
  /** free-form facts of the whole stream (window, totals…) */
  meta: Record<string, unknown>;
}

/**
 * Section names of the binary directory, by code. Append only: a code never changes meaning. A name that is not here
 * is written as NAME_ESCAPE followed by the name itself, so a section added later still travels; a reader that meets
 * a code it does not know names the section `#code` and skips it (its length is in the directory).
 *
 * `mods` (5), `methylation` (9) and `seq` (12, whole read sequences) are companion sections written by other
 * producers of these streams; this viewer does not write them and skips them unless the caller of decodeReads decodes
 * them itself. Code 10 is the escape,
 * fixed for good, so growing the table never moves it; the empty name at 10 only keeps the codes in place.
 */
const SECTION_NAMES = ['core', 'pairs', 'clips', 'inserts', 'sa', 'mods', 'reference', 'runs', 'junctions', 'methylation', '', 'hap', 'seq'];
const NAME_ESCAPE = 10;
const F_BLOCK = 1, F_START = 2, F_END = 4, F_N = 8;
const M_WINDOW = 1, M_TOTAL = 2, M_READS = 4, M_SOURCE = 8;

/**
 * The directory as varints (version 4, docs/embedded-format.md 2). Version 1-3 wrote it as JSON, about 226 bytes a
 * stream, which is most of a small window's stream (a few dozen reads) and adds up over the many windows of a
 * page. This is about a tenth of it.
 */
function encodeDirectory(dir: Directory): Uint8Array {
  const w = new Writer();
  w.u(dir.v); w.u(dir.kind === 'reads' ? 0 : 1); w.u(dir.sections.length);
  let prev = 0;
  for (const s of dir.sections) {
    const code = s.name ? SECTION_NAMES.indexOf(s.name) : -1;
    w.u(code >= 0 ? code : NAME_ESCAPE);
    if (code < 0) w.str(s.name);
    const flags = (s.block != null ? F_BLOCK : 0) | (s.start != null ? F_START : 0) | (s.end != null ? F_END : 0) | (s.n != null ? F_N : 0);
    w.u(flags); w.u(s.bytes);
    if (s.block != null) w.u(s.block);
    if (s.start != null) { w.s(s.start - prev); prev = s.start; }
    if (s.end != null) w.s(s.end - (s.start ?? prev));
    if (s.n != null) w.u(s.n);
  }
  // meta: the fields every reads stream has, then whatever else as JSON (empty when nothing)
  const { window, total, reads, reference_source, ...rest } = dir.meta as { window?: { start: number; end: number }; total?: number; reads?: number; reference_source?: string | null };
  const mask = (window ? M_WINDOW : 0) | (total != null ? M_TOTAL : 0) | (reads != null ? M_READS : 0) | (reference_source !== undefined ? M_SOURCE : 0);
  w.u(mask);
  if (window) { w.u(window.start); w.s(window.end - window.start); }
  if (total != null) w.u(total);
  if (reads != null) w.u(reads);
  if (reference_source !== undefined) w.str(reference_source ?? '');
  w.str(Object.keys(rest).length ? JSON.stringify(rest) : '');
  return w.done();
}
function decodeDirectory(bytes: Uint8Array): Directory {
  const rd = new Reader(bytes);
  const v = rd.u(), kind = rd.u() === 0 ? 'reads' : 'coverage', n = rd.u();
  const sections: SectionEntry[] = [];
  let prev = 0;
  for (let i = 0; i < n; i++) {
    const code = rd.u();
    const name = code === NAME_ESCAPE ? rd.str() : SECTION_NAMES[code] ?? `#${code}`;
    const flags = rd.u();
    const s: SectionEntry = { name, bytes: rd.u() };
    if (flags & F_BLOCK) s.block = rd.u();
    if (flags & F_START) { s.start = prev + rd.s(); prev = s.start; }
    if (flags & F_END) s.end = (s.start ?? prev) + rd.s();
    if (flags & F_N) s.n = rd.u();
    sections.push(s);
  }
  const mask = rd.u();
  const meta: Record<string, unknown> = {};
  if (mask & M_WINDOW) { const start = rd.u(); meta.window = { start, end: start + rd.s() }; }
  if (mask & M_TOTAL) meta.total = rd.u();
  if (mask & M_READS) meta.reads = rd.u();
  if (mask & M_SOURCE) meta.reference_source = rd.str() || null;
  const rest = rd.str();
  if (rest) Object.assign(meta, JSON.parse(rest));
  return { v, kind, sections, meta };
}

/** Stream = [u32 little-endian directory length][directory][section bytes in directory order]. */
function packSections(dir: Directory, parts: Uint8Array[]): Uint8Array {
  const head = encodeDirectory(dir);
  const total = 4 + head.length + parts.reduce((n, p) => n + p.length, 0);
  const out = new Uint8Array(total);
  new DataView(out.buffer).setUint32(0, head.length, true);
  out.set(head, 4);
  let o = 4 + head.length;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}
/** The directory of a stream: JSON up to version 3 (it starts with "{"), varints from version 4. */
export function readDirectory(stream: Uint8Array): { dir: Directory; offsets: number[] } {
  const headLen = new DataView(stream.buffer, stream.byteOffset, stream.byteLength).getUint32(0, true);
  const head = stream.subarray(4, 4 + headLen);
  const dir = head[0] === 0x7b ? JSON.parse(new TextDecoder().decode(head)) as Directory : decodeDirectory(head);
  if (typeof dir.v !== 'number' || dir.v > COLUMNAR_VERSION) throw new Error(`unknown columnar version ${dir.v}`);
  const offsets: number[] = []; let o = 4 + headLen;
  for (const s of dir.sections) { offsets.push(o); o += s.bytes; }
  return { dir, offsets };
}

// ======================== reads ========================

/** The reads of one window, as the export stores them and the page reads them back. */
export interface ReadsPayload {
  window: { start: number; end: number };
  /** reads passing the filters in the window, before the cap */
  total: number;
  reads: AlignedRead[];
  reference: { start: number; seq: string } | null;
  reference_source: 'fasta' | 'ensembl' | 'browser' | null;
}

const BASE_CODE: Record<string, number> = { A: 0, C: 1, G: 2, T: 3, N: 4 };
const CODE_BASE = 'ACGTN';

/** Op codes of the alignment shape (version 5): the CIGAR operations the reads keep, =/X counted as M. */
const OP_M = 0, OP_I = 1, OP_D = 2, OP_N = 3;

/**
 * The alignment shape of a read as CIGAR-like ops, from its blocks, deletions and insertions: one M per block, an
 * insertion before the reference base at its position (inside a block, or at a block's edge), and the gap between
 * two blocks as its deletions and, around them, skipped bases (N). A trailing gap (a CIGAR ending in D or N) is kept
 * from the read's end. `encodeRead` makes one block per M op, so this gives back the CIGAR the read came from, clips
 * aside.
 */
function opsOf(r: AlignedRead): number[] {
  const out: number[] = [];                          // op, length, op, length…
  const ins = r.i.length > 1 ? [...r.i].sort((a, b) => a[0] - b[0]) : r.i;
  let ii = 0, di = 0;
  const dels = r.d;
  const gap = (from: number, to: number) => {        // the bases between two blocks: deletions, skips, and an insertion between them
    let pos = from;
    while (di < dels.length && dels[di][0] < from) di++;
    for (;;) {
      const a = di < dels.length && dels[di][0] < to ? dels[di][0] : Infinity;
      const ip = ii < ins.length && ins[ii][0] < to ? ins[ii][0] : Infinity;
      if (a === Infinity && ip === Infinity) break;
      if (ip <= a) { if (ip > pos) { out.push(OP_N, ip - pos); pos = ip; } out.push(OP_I, ins[ii++][1]); }
      else { if (a > pos) out.push(OP_N, a - pos); out.push(OP_D, dels[di][1] - a); pos = dels[di++][1]; }
    }
    if (to > pos) out.push(OP_N, to - pos);
  };
  for (let k = 0; k < r.b.length; k++) {
    const [bs, be] = r.b[k];
    let p = bs;
    for (; ii < ins.length && ins[ii][0] <= be; ii++) {
      const [ip, il] = ins[ii];
      if (ip > p) { out.push(OP_M, ip - p); p = ip; }
      out.push(OP_I, il);
    }
    if (be > p) out.push(OP_M, be - p);
    const next = k + 1 < r.b.length ? r.b[k + 1][0] : r.e;
    if (next > be) gap(be, next);
  }
  for (; ii < ins.length; ii++) out.push(OP_I, ins[ii][1]);   // a read with no block (unmapped) has none; kept anyway
  return out;
}

/**
 * One block of reads as columns (version 5). Reads are sorted by start; each column lists one field for every read
 * of the block. The alignment shape is a stream of ops, split into four columns so that each holds one kind of number:
 * how many ops a read has, their codes (one byte each), the lengths of the M ops, and those of the others (I, D, N).
 * Blocks, deletions, insertions and the read's end are all rebuilt from it; up to version 4 they were stored apart
 * (the end, the first block, each later block as gap and length, each deletion and insertion again with its offset),
 * which for a long read with hundreds of indels cost three times the ops (1,726 against 562 bytes per ONT read).
 * Mismatches are chained positions from the read start, then base and quality; soft-clip lengths close the block.
 */
function encodeReadBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  w.u(reads.length);
  let prev = 0;
  for (const r of reads) { w.u(r.s - prev); prev = r.s; }                 // starts (deltas)
  for (const r of reads) w.u(r.f);                                        // flags
  for (const r of reads) w.u(r.q);                                        // MAPQ
  for (const r of reads) w.s(r.nh == null ? -1 : r.nh);                   // NH (-1 = absent)
  const ops = reads.map(opsOf);
  for (const o of ops) w.u(o.length >> 1);                                // op counts
  for (const o of ops) for (let k = 0; k < o.length; k += 2) w.byte(o[k]);                     // op codes
  for (const o of ops) for (let k = 0; k < o.length; k += 2) if (o[k] === OP_M) w.u(o[k + 1]); // M lengths
  for (const o of ops) for (let k = 0; k < o.length; k += 2) if (o[k] !== OP_M) w.u(o[k + 1]); // I, D, N lengths
  for (const r of reads) {
    w.u(r.m.length);
    let mp = r.s;
    for (const [pos, base, qual] of r.m) { w.u(pos - mp); mp = pos; w.byte(BASE_CODE[base] ?? 4); w.byte(Math.min(255, Math.max(0, qual))); }
  }
  for (const r of reads) { w.u(r.c[0]); w.u(r.c[1]); }                    // soft clips
  return w.done();
}
function decodeReadBlock(bytes: Uint8Array, names: (i: number) => string, version: number): AlignedRead[] {
  if (version < 5) return decodeReadBlockV4(bytes, names);
  const rd = new Reader(bytes);
  const n = rd.u();
  const s = new Array<number>(n), f = new Array<number>(n), q = new Array<number>(n), nh = new Array<number | null>(n), cnt = new Array<number>(n);
  let prev = 0;
  for (let i = 0; i < n; i++) { prev += rd.u(); s[i] = prev; }
  for (let i = 0; i < n; i++) f[i] = rd.u();
  for (let i = 0; i < n; i++) q[i] = rd.u();
  for (let i = 0; i < n; i++) { const v = rd.s(); nh[i] = v < 0 ? null : v; }
  let total = 0;
  for (let i = 0; i < n; i++) { cnt[i] = rd.u(); total += cnt[i]; }
  const codes = new Uint8Array(total);
  for (let k = 0; k < total; k++) codes[k] = rd.byte();
  const lens = new Array<number>(total);
  for (let k = 0; k < total; k++) if (codes[k] === OP_M) lens[k] = rd.u();
  for (let k = 0; k < total; k++) if (codes[k] !== OP_M) lens[k] = rd.u();
  const out: AlignedRead[] = new Array(n);
  let k = 0;
  for (let i = 0; i < n; i++) {
    const b: [number, number][] = [], d: [number, number][] = [], ins: [number, number][] = [];
    let pos = s[i];
    for (const end = k + cnt[i]; k < end; k++) {
      const len = lens[k];
      switch (codes[k]) {
        case OP_M: b.push([pos, pos + len]); pos += len; break;
        case OP_I: ins.push([pos, len]); break;
        case OP_D: d.push([pos, pos + len]); pos += len; break;
        default: pos += len;
      }
    }
    out[i] = { n: names(i), s: s[i], e: pos, r: (f[i] & 16) ? 1 : 0, q: q[i], f: f[i], nh: nh[i], b, d, i: ins, m: [], c: [0, 0] };
  }
  for (let i = 0; i < n; i++) {
    const c = rd.u(); let mp = s[i];
    for (let j = 0; j < c; j++) { mp += rd.u(); const base = CODE_BASE[rd.byte()] ?? 'N'; out[i].m.push([mp, base, rd.byte()]); }
  }
  for (let i = 0; i < n; i++) out[i].c = [rd.u(), rd.u()];
  return out;
}
/** The core block of versions 1–4: the end, the blocks, then deletions and insertions with their offsets. */
function decodeReadBlockV4(bytes: Uint8Array, names: (i: number) => string): AlignedRead[] {
  const rd = new Reader(bytes);
  const n = rd.u();
  const s = new Array<number>(n), len = new Array<number>(n), f = new Array<number>(n), q = new Array<number>(n), nh = new Array<number | null>(n);
  let prev = 0;
  for (let i = 0; i < n; i++) { prev += rd.u(); s[i] = prev; }
  for (let i = 0; i < n; i++) len[i] = rd.u();
  for (let i = 0; i < n; i++) f[i] = rd.u();
  for (let i = 0; i < n; i++) q[i] = rd.u();
  for (let i = 0; i < n; i++) { const v = rd.s(); nh[i] = v < 0 ? null : v; }
  const extra: [number, number][][] = [];
  for (let i = 0; i < n; i++) { const k = rd.u(); const list: [number, number][] = []; for (let j = 0; j < k; j++) list.push([rd.u(), rd.u()]); extra.push(list); }
  const first = new Array<number>(n);
  for (let i = 0; i < n; i++) first[i] = rd.u();
  const out: AlignedRead[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const b: [number, number][] = [[s[i], s[i] + first[i]]];
    for (const [gap, l] of extra[i]) { const st = b[b.length - 1][1] + gap; b.push([st, st + l]); }
    out[i] = { n: names(i), s: s[i], e: s[i] + len[i], r: (f[i] & 16) ? 1 : 0, q: q[i], f: f[i], nh: nh[i], b, d: [], i: [], m: [], c: [0, 0] };
  }
  for (let i = 0; i < n; i++) { const k = rd.u(); for (let j = 0; j < k; j++) { const a = s[i] + rd.u(); out[i].d.push([a, a + rd.u()]); } }
  for (let i = 0; i < n; i++) { const k = rd.u(); for (let j = 0; j < k; j++) out[i].i.push([s[i] + rd.u(), rd.u()]); }
  for (let i = 0; i < n; i++) {
    const k = rd.u(); let mp = s[i];
    for (let j = 0; j < k; j++) { mp += rd.u(); const base = CODE_BASE[rd.byte()] ?? 'N'; out[i].m.push([mp, base, rd.byte()]); }
  }
  for (let i = 0; i < n; i++) out[i].c = [rd.u(), rd.u()];
  return out;
}

/** Span of a read pair as the aligner writes TLEN: leftmost to rightmost mapped base, negative for the right read. */
const derivedTlen = (r: { s: number; e: number }, m: { s: number; e: number }) =>
  (r.s <= m.s ? 1 : -1) * (Math.max(r.e, m.e) - Math.min(r.s, m.s));

/**
 * For each read of a stream, the signed distance to a read of the SAME stream that starts where its mate starts,
 * or 0 when there is none.
 *
 * Both halves of a pair are nearly always in the same window (98.2% measured on 100 kb windows of a paired RNA-seq
 * file), and the mate
 * then carries its own start and end: storing where it is costs a small signed distance instead of a position, and
 * the template length becomes a residual that is zero for all but a handful of reads. That is what makes `pairs`,
 * the bulkiest section of paired data, about a third of its former size.
 *
 * For the position and the template length, the link only has to land on a read whose start equals `mp`: that
 * reproduces the mate position exactly and the residual absorbs the rest. But a link that runs both ways is also
 * taken as *the* pairing of the two reads (`mk`, see decodeReads), which is what joins the right two when several
 * reads share a start — so the true mate is preferred: the read carrying the same mate key, or the same name, the
 * other pair bit and a primary record. Then the read whose residual vanishes, then the nearest, the shortest varint.
 */
function pairLinks(reads: AlignedRead[]): Int32Array {
  const byStart = new Map<number, number[]>();
  for (let i = 0; i < reads.length; i++) {
    const list = byStart.get(reads[i].s);
    if (list) list.push(i); else byStart.set(reads[i].s, [i]);
  }
  const links = new Int32Array(reads.length);
  for (let i = 0; i < reads.length; i++) {
    const r = reads[i];
    if (r.mp == null || r.mc) continue;             // unpaired, or a mate on another chromosome: not in this stream
    const cands = byStart.get(r.mp);
    if (!cands) continue;
    let best = 0, bestKey = Infinity;
    for (const j of cands) {
      if (j === i) continue;
      const m = reads[j];
      const mate = r.mk != null ? m.mk === r.mk
        : m.n === r.n && (m.f & 192) !== (r.f & 192) && !(m.f & 0x900);
      const key = (mate ? 0 : 2e9) + ((r.tl ?? 0) === derivedTlen(r, m) ? 0 : 1e9) + Math.abs(j - i);
      if (key < bestKey) { bestKey = key; best = j - i; }
    }
    links[i] = best;
  }
  return links;
}

/**
 * Pair fields of a block (docs/embedded-format.md 3.2, version 5): the links first, each with a bit saying whether the
 * read it points to, later in the same block, links back; such a read then has no entry of its own (its link is the
 * way back). Then, for the reads with no link, the mate chromosome (dictionary) and the mate start; and a template
 * length, residual or verbatim, for every paired read. A link to another block is written by both reads, so that a
 * block decodes without the ones before it.
 *
 * `all` is the whole stream and `base` the index of `block[0]` in it, because a link may cross a block boundary.
 */
function encodePairBlock(block: AlignedRead[], all: AlignedRead[], base: number, links: Int32Array): Uint8Array {
  const w = new Writer();
  const chroms: string[] = []; const idx = new Map<string, number>();
  const chromIdx = (c: string) => { let i = idx.get(c); if (i == null) { i = chroms.length; chroms.push(c); idx.set(c, i); } return i; };
  const n = block.length;
  const implied = new Uint8Array(n);                 // the second read of a mutual link inside the block: nothing written
  for (let i = 0; i < n; i++) {
    const d = links[base + i];
    if (d > 0 && i + d < n && links[base + i + d] === -d) implied[i + d] = 1;
  }
  const codes = block.map((r, i) => (links[base + i] ? 0 : r.mp == null ? -1 : r.mc ? chromIdx(r.mc) + 1 : 0));   // -1 no mate, 0 same chromosome, k+1 = chroms[k]
  w.u(chroms.length); for (const c of chroms) w.str(c);
  for (let i = 0; i < n; i++) {
    if (implied[i]) continue;
    const d = links[base + i];
    // the bit: the read linked to links back to THIS read. Another read may link to it too (look-alike reads at one
    // start); up to now that one got the bit as well, and the target decoded with its mate instead of its own
    w.u(d ? (zz(d) << 1) | (d > 0 && i + d < n && links[base + i + d] === -d ? 1 : 0) : 0);
  }
  for (let i = 0; i < n; i++) if (!links[base + i]) w.s(codes[i]);
  for (let i = 0; i < n; i++) if (!links[base + i] && codes[i] >= 0) w.s(block[i].mp! - block[i].s);
  for (let i = 0; i < n; i++) {
    if (codes[i] < 0) continue;
    const d = links[base + i];
    w.s(d ? (block[i].tl ?? 0) - derivedTlen(block[i], all[base + i + d]) : (block[i].tl ?? 0));
  }
  return w.done();
}
/** Where a linked mate is read from: the read at that index in the stream, decoding its block if it is not out yet. */
type MateLookup = (index: number) => Promise<{ s: number; e: number } | null>;

/** Decodes a block's pair fields; from version 3 also returns, per read, the stream index its link points to (-1 none). */
async function decodePairBlock(bytes: Uint8Array, reads: AlignedRead[], version: number, base: number, mateAt: MateLookup): Promise<Int32Array | null> {
  const rd = new Reader(bytes);
  const k = rd.u(); const chroms: string[] = []; for (let i = 0; i < k; i++) chroms.push(rd.str());
  const codes = version < 5 ? reads.map(() => rd.s()) : new Array<number>(reads.length);
  if (version < 3) {
    const mp = reads.map(() => rd.s());
    const tl = reads.map(() => rd.s());
    reads.forEach((r, i) => {
      if (codes[i] < 0) return;
      r.mp = r.s + mp[i]; r.tl = tl[i];
      if (codes[i] > 0) r.mc = chroms[codes[i] - 1];
    });
    return null;
  }
  let links: number[];
  if (version < 5) {
    links = reads.map(() => rd.s());
  } else {
    // the links first; a mutual one also gives the read it points to its link back, which that read then omits
    links = new Array(reads.length).fill(0);
    const implied = new Uint8Array(reads.length);
    for (let i = 0; i < reads.length; i++) {
      if (implied[i]) continue;
      const v = rd.u();
      if (!v) continue;
      const d = unzz(v >> 1);
      links[i] = d;
      if (v & 1) { links[i + d] = -d; implied[i + d] = 1; }
    }
    for (let i = 0; i < reads.length; i++) codes[i] = links[i] ? 0 : rd.s();
  }
  const mp = reads.map((_, i) => (codes[i] >= 0 && !links[i] ? rd.s() : 0));
  const tl = reads.map((_, i) => (codes[i] >= 0 ? rd.s() : 0));
  const targets = new Int32Array(reads.length).fill(-1);
  for (let i = 0; i < reads.length; i++) {
    const r = reads[i];
    if (codes[i] < 0) continue;
    if (codes[i] > 0) r.mc = chroms[codes[i] - 1];   // a linked mate is in this stream, so on this chromosome
    if (links[i]) {
      const m = await mateAt(base + i + links[i]);
      if (!m) continue;                              // the mate's block is gone: better no mate than a wrong one
      r.mp = m.s; r.tl = tl[i] + derivedTlen(r, m);
      targets[i] = base + i + links[i];
    } else {
      r.mp = r.s + mp[i]; r.tl = tl[i];
    }
  }
  return targets;
}

/** Clip fields of a block: soft-clipped bases at each end and hard-clipped lengths. */
function encodeClipBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  for (const r of reads) { w.str(r.cs?.[0] ?? ''); w.str(r.cs?.[1] ?? ''); w.u(r.h?.[0] ?? 0); w.u(r.h?.[1] ?? 0); }
  return w.done();
}
function decodeClipBlock(bytes: Uint8Array, reads: AlignedRead[]) {
  const rd = new Reader(bytes);
  for (const r of reads) {
    const l = rd.str(), rt = rd.str(), hl = rd.u(), hr = rd.u();
    if (l || rt) r.cs = [l, rt];
    if (hl || hr) r.h = [hl, hr];
  }
}
/** Inserted bases of a block: one string per insertion of each read, in the order of the core section's insertions. */
function encodeInsertBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  for (const r of reads) for (let k = 0; k < r.i.length; k++) w.str(r.is?.[k] ?? '');
  return w.done();
}
function decodeInsertBlock(bytes: Uint8Array, reads: AlignedRead[]) {
  const rd = new Reader(bytes);
  for (const r of reads) { if (!r.i.length) continue; const list = r.i.map(() => rd.str()); if (list.some(x => x)) r.is = list; }
}
/** SA tags of a block: one string per read (empty when the read is not split). */
function encodeSaBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  for (const r of reads) w.str(r.sa ?? '');
  return w.done();
}
function decodeSaBlock(bytes: Uint8Array, reads: AlignedRead[]) {
  const rd = new Reader(bytes);
  for (const r of reads) { const v = rd.str(); if (v) r.sa = v; }
}
/** Haplotags of a block: per read HP (0 when untagged), then for a tagged read PS + 1 and PC + 1 (0 when absent). */
function encodeHapBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  for (const r of reads) { w.u(r.hp ?? 0); if (r.hp) { w.u(r.ps != null ? r.ps + 1 : 0); w.u(r.pc != null ? r.pc + 1 : 0); } }
  return w.done();
}
function decodeHapBlock(bytes: Uint8Array, reads: AlignedRead[]) {
  const rd = new Reader(bytes);
  for (const r of reads) {
    const hp = rd.u();
    if (!hp) continue;
    r.hp = hp;
    const ps = rd.u(), pc = rd.u();
    if (ps) r.ps = ps - 1;
    if (pc) r.pc = pc - 1;
  }
}

/** Encodes the reads of a window: blocks of READS_PER_BLOCK, a pairs section per block when any read has a mate, the reference bases. */
export async function encodeReads(p: ReadsPayload): Promise<Uint8Array> {
  const sorted = [...p.reads].sort((a, b) => a.s - b.s || a.e - b.e);
  const sections: SectionEntry[] = []; const parts: Uint8Array[] = [];
  const anyPair = sorted.some(r => r.mp != null);
  const links = anyPair ? pairLinks(sorted) : new Int32Array(0);
  for (let b = 0; b * READS_PER_BLOCK < sorted.length; b++) {
    const base = b * READS_PER_BLOCK;
    const block = sorted.slice(base, base + READS_PER_BLOCK);
    const core = await deflate(encodeReadBlock(block));
    sections.push({ name: 'core', block: b, start: block[0].s, end: Math.max(...block.map(r => r.e)), n: block.length, bytes: core.length });
    parts.push(core);
    if (anyPair) {
      const pairs = await deflate(encodePairBlock(block, sorted, base, links));
      sections.push({ name: 'pairs', block: b, bytes: pairs.length }); parts.push(pairs);
    }
    if (block.some(r => r.cs || r.h)) {
      const clips = await deflate(encodeClipBlock(block));
      sections.push({ name: 'clips', block: b, bytes: clips.length }); parts.push(clips);
    }
    if (block.some(r => r.is)) {
      const ins = await deflate(encodeInsertBlock(block));
      sections.push({ name: 'inserts', block: b, bytes: ins.length }); parts.push(ins);
    }
    if (block.some(r => r.sa)) {
      const sa = await deflate(encodeSaBlock(block));
      sections.push({ name: 'sa', block: b, bytes: sa.length }); parts.push(sa);
    }
    if (block.some(r => r.hp)) {
      const hap = await deflate(encodeHapBlock(block));
      sections.push({ name: 'hap', block: b, bytes: hap.length }); parts.push(hap);
    }
  }
  if (p.reference) {
    const w = new Writer(); w.u(p.reference.start); w.str(p.reference.seq);
    const ref = await deflate(w.done());
    sections.push({ name: 'reference', bytes: ref.length }); parts.push(ref);
  }
  const dir: Directory = { v: COLUMNAR_VERSION, kind: 'reads', sections, meta: { window: p.window, total: p.total, reads: sorted.length, reference_source: p.reference_source } };
  return packSections(dir, parts);
}

/**
 * Walks the alignment shapes of a reads stream without building a read: for each read of the blocks overlapping
 * `range` (every block without one), `visit(start, flags, mapq, nh, ops)` with `ops` the read's CIGAR ops, clips aside,
 * packed as a BAM record packs them (length × 16 + op: M 0, I 1, D 2, N 3). `ops` is valid during the call only.
 *
 * This is what a coverage, a junction count or the reads through an exon boundary need of a read, and in version 5 it
 * is the first columns of the core section: the mismatches, the clips and every other section are not read. A reader
 * counting depth this way decodes about as many bytes as the reads' starts and shapes, a fraction of a full decode.
 * Returns false, having visited nothing, for a stream before version 5, whose shape a caller takes from decodeReads.
 */
export async function scanShapes(stream: Uint8Array, range: { start: number; end: number } | undefined,
  visit: (start: number, flags: number, mapq: number, nh: number | null, ops: Uint32Array) => void): Promise<boolean> {
  const { dir, offsets } = readDirectory(stream);
  if (dir.v < 5) return false;
  let ops = new Uint32Array(256);
  for (let i = 0; i < dir.sections.length; i++) {
    const s = dir.sections[i];
    if (s.name !== 'core' || (range && (s.end! <= range.start || s.start! >= range.end))) continue;
    const rd = new Reader(await inflate(stream.subarray(offsets[i], offsets[i] + s.bytes)));
    const n = rd.u();
    const st = new Float64Array(n), f = new Int32Array(n), q = new Int32Array(n), nh = new Int32Array(n), cnt = new Int32Array(n);
    let prev = 0, total = 0;
    for (let k = 0; k < n; k++) { prev += rd.u(); st[k] = prev; }
    for (let k = 0; k < n; k++) f[k] = rd.u();
    for (let k = 0; k < n; k++) q[k] = rd.u();
    for (let k = 0; k < n; k++) nh[k] = rd.s();
    for (let k = 0; k < n; k++) { cnt[k] = rd.u(); total += cnt[k]; }
    const codes = rd.bytes(total);
    const lens = new Float64Array(total);
    for (let k = 0; k < total; k++) if (codes[k] === OP_M) lens[k] = rd.u();
    for (let k = 0; k < total; k++) if (codes[k] !== OP_M) lens[k] = rd.u();
    for (let r = 0, k = 0; r < n; r++) {
      const c = cnt[r];
      if (c > ops.length) ops = new Uint32Array(Math.max(c, ops.length * 2));
      let end = st[r];
      for (let j = 0; j < c; j++, k++) {
        const len = lens[k], op = codes[k];
        ops[j] = len * 16 + op;
        if (op !== OP_I) end += len;
      }
      if (range && (end <= range.start || st[r] >= range.end) && !(end === st[r] && st[r] >= range.start && st[r] < range.end)) continue;
      visit(st[r], f[r], q[r], nh[r] < 0 ? null : nh[r], ops.subarray(0, c));
    }
  }
  return true;
}

/** Facts of an encoded reads stream, without decoding a read. */
export function readsInfo(stream: Uint8Array): { window: { start: number; end: number }; total: number; reads: number; blocks: { start: number; end: number; n: number }[] } {
  const { dir } = readDirectory(stream);
  const meta = dir.meta as { window: { start: number; end: number }; total: number; reads: number };
  return { window: meta.window, total: meta.total, reads: meta.reads, blocks: dir.sections.filter(s => s.name === 'core').map(s => ({ start: s.start!, end: s.end!, n: s.n! })) };
}

/** What a caller of decodeReads may ask for beyond the reads themselves. */
export interface DecodeReadsOptions {
  /** decode this block alone (its number in the directory), whatever the range */
  onlyBlock?: number;
  /**
   * Decoders of companion sections this module does not know, by name: each gets the inflated bytes, the reads of the
   * block the section belongs to (in stream order) and the stream's version, and adds what it reads to those reads.
   */
  sections?: Record<string, (bytes: Uint8Array, block: AlignedRead[], version: number) => void>;
}

/**
 * Decodes the reads of a stream, only the blocks overlapping [start, end) when a range is given. Read names are
 * synthetic ("read N", N counted over the whole stream) since the export drops them.
 *
 * From version 3 a read's mate is stored as a link to another read of the stream, which may sit in a block this
 * call is filtering out; such a block is then inflated on the side, for its starts and ends alone, and cached in
 * case the reader asks for it next. Decoding a whole stream, the usual case, never pays for it: every block it
 * needs is one it was going to decode anyway.
 *
 * Two returned reads whose links point at each other, one of each pair bit, get the same mate key `mk`, and the
 * viewer joins them on it rather than by position (`pairMates`). The key is the stream's window start and the lower
 * of the two indices, so decoding the same stream twice gives the same keys and reads of different streams never share
 * one. A one-way link is not a pairing — the writer may have linked a read to another that merely starts where its
 * mate does — and is left to the positional rule, as is every read of a stream before version 3.
 */
export async function decodeReads(stream: Uint8Array, range?: { start: number; end: number }, opts?: DecodeReadsOptions): Promise<ReadsPayload> {
  const onlyBlock = opts?.onlyBlock;
  const { dir, offsets } = readDirectory(stream);
  const meta = dir.meta as { window: { start: number; end: number }; total: number; reference_source: ReadsPayload['reference_source'] };
  const reads: AlignedRead[] = [];
  let reference: ReadsPayload['reference'] = null;
  let counted = 0;

  const coreAt = new Map<number, number>();          // block number → its section index
  for (let i = 0; i < dir.sections.length; i++) if (dir.sections[i].name === 'core') coreAt.set(dir.sections[i].block!, i);
  const decoded = new Map<number, AlignedRead[]>();  // block number → its reads, whether or not they are being returned
  const blockOf = async (b: number): Promise<AlignedRead[] | null> => {
    const have = decoded.get(b);
    if (have) return have;
    const si = coreAt.get(b);
    if (si == null) return null;
    const s = dir.sections[si];
    const out = decodeReadBlock(await inflate(stream.subarray(offsets[si], offsets[si] + s.bytes)), k => `read ${k + 1}`, dir.v);
    decoded.set(b, out);
    return out;
  };
  const mateAt: MateLookup = async index => {
    const block = await blockOf(Math.floor(index / READS_PER_BLOCK));
    return block?.[index % READS_PER_BLOCK] ?? null;
  };
  const returned = new Map<number, AlignedRead>();   // stream index → read, for the reads this call returns
  const linkOf = new Map<number, number>();          // stream index → the stream index its link points to

  for (let i = 0; i < dir.sections.length; i++) {
    const s = dir.sections[i];
    const bytes = () => stream.subarray(offsets[i], offsets[i] + s.bytes);
    if (s.name === 'core') {
      const base = counted; counted += s.n ?? 0;
      if (onlyBlock != null && s.block !== onlyBlock) continue;
      if (range && (s.end! <= range.start || s.start! >= range.end)) continue;
      const block = await blockOf(s.block!) ?? decodeReadBlock(await inflate(bytes()), k => `read ${k + 1}`, dir.v);
      block.forEach((r, k) => { r.n = `read ${base + k + 1}`; });   // named now: a mate lookup may have decoded it first
      // the block's companion sections follow it: pairs, clips, inserts, sa, hap, then any the caller decodes itself
      // (`opts.sections`); the rest are skipped
      for (let j = i + 1; j < dir.sections.length && dir.sections[j].name !== 'core' && dir.sections[j].block === s.block; j++) {
        const c = dir.sections[j];
        const data = () => inflate(stream.subarray(offsets[j], offsets[j] + c.bytes));
        if (c.name === 'pairs') {
          const targets = await decodePairBlock(await data(), block, dir.v, s.block! * READS_PER_BLOCK, mateAt);
          targets?.forEach((t, k) => { if (t >= 0) linkOf.set(s.block! * READS_PER_BLOCK + k, t); });
        }
        else if (c.name === 'clips') decodeClipBlock(await data(), block);
        else if (c.name === 'inserts') decodeInsertBlock(await data(), block);
        else if (c.name === 'sa') decodeSaBlock(await data(), block);
        else if (c.name === 'hap') decodeHapBlock(await data(), block);
        else if (opts?.sections?.[c.name]) opts.sections[c.name](await data(), block, dir.v);
      }
      block.forEach((r, k) => returned.set(s.block! * READS_PER_BLOCK + k, r));
      reads.push(...block);
    } else if (s.name === 'reference') {
      const rd = new Reader(await inflate(bytes()));
      reference = { start: rd.u(), seq: rd.str() };
    }
    // unknown sections (written by a newer version) are skipped
  }
  for (const [g, t] of linkOf) {
    if (t <= g || linkOf.get(t) !== g) continue;
    const a = returned.get(g), b = returned.get(t);
    if (!a || !b || (a.f & 192) === (b.f & 192)) continue;
    a.mk = b.mk = `${meta.window.start}:${g}`;
  }
  return { window: meta.window, total: meta.total, reads, reference, reference_source: meta.reference_source ?? null };
}

// ======================== coverage ========================

/** Coverage and junctions of one window, the part of SampleCoverage that is bulky. */
export interface CoveragePayload {
  runs: CoverageRun[];
  junctions: JunctionArc[];
}

/** Runs as (length, depth delta) columns, junctions as (start delta, length, count). Everything else of a coverage answer stays JSON. */
export async function encodeCoverage(p: CoveragePayload): Promise<Uint8Array> {
  const w = new Writer();
  w.u(p.runs.length);
  w.u(p.runs.length ? p.runs[0].start : 0);
  let prevDepth = 0;
  for (const r of p.runs) w.u(r.end - r.start);
  for (const r of p.runs) { w.s(r.depth - prevDepth); prevDepth = r.depth; }
  const runs = await deflate(w.done());
  const wj = new Writer();
  const js = [...p.junctions].sort((a, b) => a.start - b.start || a.end - b.end);
  wj.u(js.length);
  let prev = 0;
  for (const j of js) { wj.u(j.start - prev); prev = j.start; }
  for (const j of js) wj.u(j.end - j.start);
  for (const j of js) wj.u(Math.round(j.count));
  const junctions = await deflate(wj.done());
  const sections: SectionEntry[] = [{ name: 'runs', bytes: runs.length }, { name: 'junctions', bytes: junctions.length }];
  const parts = [runs, junctions];
  const dir: Directory = { v: COLUMNAR_VERSION, kind: 'coverage', sections, meta: {} };
  return packSections(dir, parts);
}
export async function decodeCoverage(stream: Uint8Array): Promise<CoveragePayload> {
  const { dir, offsets } = readDirectory(stream);
  const out: CoveragePayload = { runs: [], junctions: [] };
  for (let i = 0; i < dir.sections.length; i++) {
    const s = dir.sections[i];
    const bytes = await inflate(stream.subarray(offsets[i], offsets[i] + s.bytes));
    if (s.name === 'runs') {
      const rd = new Reader(bytes); const n = rd.u(); let pos = rd.u(); let depth = 0;
      const lens = new Array<number>(n); for (let k = 0; k < n; k++) lens[k] = rd.u();
      for (let k = 0; k < n; k++) { depth += rd.s(); out.runs.push({ start: pos, end: pos + lens[k], depth }); pos += lens[k]; }
    } else if (s.name === 'junctions') {
      const rd = new Reader(bytes); const n = rd.u();
      const starts = new Array<number>(n); let prev = 0; for (let k = 0; k < n; k++) { prev += rd.u(); starts[k] = prev; }
      const lens = new Array<number>(n); for (let k = 0; k < n; k++) lens[k] = rd.u();
      const counts = new Array<number>(n); for (let k = 0; k < n; k++) counts[k] = rd.u();
      for (let k = 0; k < n; k++) out.junctions.push({ start: starts[k], end: starts[k] + lens[k], count: counts[k] });
    }
    // other sections (a newer version's, or another producer's) are skipped
  }
  return out;
}
