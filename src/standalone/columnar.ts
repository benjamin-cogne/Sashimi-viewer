/**
 * Compact storage of reads and coverage for exported pages (and, later, converted files).
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

export const COLUMNAR_VERSION = 1;
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

/** Stream = [u32 little-endian directory length][directory JSON][section bytes in directory order]. */
function packSections(dir: Directory, parts: Uint8Array[]): Uint8Array {
  const head = new TextEncoder().encode(JSON.stringify(dir));
  const total = 4 + head.length + parts.reduce((n, p) => n + p.length, 0);
  const out = new Uint8Array(total);
  new DataView(out.buffer).setUint32(0, head.length, true);
  out.set(head, 4);
  let o = 4 + head.length;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}
export function readDirectory(stream: Uint8Array): { dir: Directory; offsets: number[] } {
  const headLen = new DataView(stream.buffer, stream.byteOffset, stream.byteLength).getUint32(0, true);
  const dir = JSON.parse(new TextDecoder().decode(stream.subarray(4, 4 + headLen))) as Directory;
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

/**
 * One block of reads as columns. Reads are sorted by start; each column lists one field for every read of the
 * block, the variable-length lists (blocks, deletions, insertions, mismatches) as a count per read followed by
 * their values, every position relative to the read start (blocks and mismatches chained: each relative to the
 * previous one).
 */
function encodeReadBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  w.u(reads.length);
  let prev = 0;
  for (const r of reads) { w.u(r.s - prev); prev = r.s; }                 // starts (deltas)
  for (const r of reads) w.u(r.e - r.s);                                  // aligned length
  for (const r of reads) w.u(r.f);                                        // flags
  for (const r of reads) w.u(r.q);                                        // MAPQ
  for (const r of reads) w.s(r.nh == null ? -1 : r.nh);                   // NH (-1 = absent)
  for (const r of reads) {                                                // aligned blocks after the first: gap then length
    w.u(r.b.length - 1);
    for (let k = 1; k < r.b.length; k++) { w.u(r.b[k][0] - r.b[k - 1][1]); w.u(r.b[k][1] - r.b[k][0]); }
  }
  for (const r of reads) { w.u(r.b.length ? r.b[0][1] - r.b[0][0] : 0); }  // first block length
  for (const r of reads) { w.u(r.d.length); for (const [a, b] of r.d) { w.u(a - r.s); w.u(b - a); } }
  for (const r of reads) { w.u(r.i.length); for (const [p, len] of r.i) { w.u(p - r.s); w.u(len); } }
  for (const r of reads) {
    w.u(r.m.length);
    let mp = r.s;
    for (const [pos, base, qual] of r.m) { w.u(pos - mp); mp = pos; w.byte(BASE_CODE[base] ?? 4); w.byte(Math.min(255, Math.max(0, qual))); }
  }
  for (const r of reads) { w.u(r.c[0]); w.u(r.c[1]); }                    // soft clips
  return w.done();
}
function decodeReadBlock(bytes: Uint8Array, names: (i: number) => string): AlignedRead[] {
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

/** Pair fields of a block: mate start as a delta from the read start, template length, mate chromosome (dictionary). */
function encodePairBlock(reads: AlignedRead[]): Uint8Array {
  const w = new Writer();
  const chroms: string[] = []; const idx = new Map<string, number>();
  const chromIdx = (c: string) => { let i = idx.get(c); if (i == null) { i = chroms.length; chroms.push(c); idx.set(c, i); } return i; };
  const codes = reads.map(r => (r.mp == null ? -1 : r.mc ? chromIdx(r.mc) + 1 : 0));   // -1 no mate, 0 same chromosome, k+1 = chroms[k]
  w.u(chroms.length); for (const c of chroms) w.str(c);
  for (let i = 0; i < reads.length; i++) w.s(codes[i]);
  for (const r of reads) w.s(r.mp == null ? 0 : r.mp - r.s);
  for (const r of reads) w.s(r.tl ?? 0);
  return w.done();
}
function decodePairBlock(bytes: Uint8Array, reads: AlignedRead[]) {
  const rd = new Reader(bytes);
  const k = rd.u(); const chroms: string[] = []; for (let i = 0; i < k; i++) chroms.push(rd.str());
  const codes = reads.map(() => rd.s());
  const mp = reads.map(() => rd.s());
  const tl = reads.map(() => rd.s());
  reads.forEach((r, i) => {
    if (codes[i] < 0) return;
    r.mp = r.s + mp[i]; r.tl = tl[i];
    if (codes[i] > 0) r.mc = chroms[codes[i] - 1];
  });
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

/** Encodes the reads of a window: blocks of READS_PER_BLOCK, a pairs section per block when any read has a mate, the reference bases. */
export async function encodeReads(p: ReadsPayload): Promise<Uint8Array> {
  const sorted = [...p.reads].sort((a, b) => a.s - b.s || a.e - b.e);
  const sections: SectionEntry[] = []; const parts: Uint8Array[] = [];
  const anyPair = sorted.some(r => r.mp != null);
  for (let b = 0; b * READS_PER_BLOCK < sorted.length; b++) {
    const block = sorted.slice(b * READS_PER_BLOCK, (b + 1) * READS_PER_BLOCK);
    const core = await deflate(encodeReadBlock(block));
    sections.push({ name: 'core', block: b, start: block[0].s, end: Math.max(...block.map(r => r.e)), n: block.length, bytes: core.length });
    parts.push(core);
    if (anyPair) {
      const pairs = await deflate(encodePairBlock(block));
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
  }
  if (p.reference) {
    const w = new Writer(); w.u(p.reference.start); w.str(p.reference.seq);
    const ref = await deflate(w.done());
    sections.push({ name: 'reference', bytes: ref.length }); parts.push(ref);
  }
  const dir: Directory = { v: COLUMNAR_VERSION, kind: 'reads', sections, meta: { window: p.window, total: p.total, reads: sorted.length, reference_source: p.reference_source } };
  return packSections(dir, parts);
}

/** Facts of an encoded reads stream, without decoding a read. */
export function readsInfo(stream: Uint8Array): { window: { start: number; end: number }; total: number; reads: number; blocks: { start: number; end: number; n: number }[] } {
  const { dir } = readDirectory(stream);
  const meta = dir.meta as { window: { start: number; end: number }; total: number; reads: number };
  return { window: meta.window, total: meta.total, reads: meta.reads, blocks: dir.sections.filter(s => s.name === 'core').map(s => ({ start: s.start!, end: s.end!, n: s.n! })) };
}

/**
 * Decodes the reads of a stream, only the blocks overlapping [start, end) when a range is given. Read names are
 * synthetic ("read N", N counted over the whole stream) since the export drops them.
 */
export async function decodeReads(stream: Uint8Array, range?: { start: number; end: number }): Promise<ReadsPayload> {
  const { dir, offsets } = readDirectory(stream);
  const meta = dir.meta as { window: { start: number; end: number }; total: number; reference_source: ReadsPayload['reference_source'] };
  const reads: AlignedRead[] = [];
  let reference: ReadsPayload['reference'] = null;
  let counted = 0;
  for (let i = 0; i < dir.sections.length; i++) {
    const s = dir.sections[i];
    const bytes = () => stream.subarray(offsets[i], offsets[i] + s.bytes);
    if (s.name === 'core') {
      const base = counted; counted += s.n ?? 0;
      if (range && (s.end! <= range.start || s.start! >= range.end)) continue;
      const block = decodeReadBlock(await inflate(bytes()), k => `read ${base + k + 1}`);
      // the block's companion sections follow it: pairs, clips, inserts, sa (any order; unknown names skipped)
      for (let j = i + 1; j < dir.sections.length && dir.sections[j].name !== 'core' && dir.sections[j].block === s.block; j++) {
        const c = dir.sections[j];
        const data = () => inflate(stream.subarray(offsets[j], offsets[j] + c.bytes));
        if (c.name === 'pairs') decodePairBlock(await data(), block);
        else if (c.name === 'clips') decodeClipBlock(await data(), block);
        else if (c.name === 'inserts') decodeInsertBlock(await data(), block);
        else if (c.name === 'sa') decodeSaBlock(await data(), block);
      }
      reads.push(...block);
    } else if (s.name === 'reference') {
      const rd = new Reader(await inflate(bytes()));
      reference = { start: rd.u(), seq: rd.str() };
    }
    // unknown sections (written by a newer version) are skipped
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
  const dir: Directory = { v: COLUMNAR_VERSION, kind: 'coverage', sections: [{ name: 'runs', bytes: runs.length }, { name: 'junctions', bytes: junctions.length }], meta: {} };
  return packSections(dir, [runs, junctions]);
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
  }
  return out;
}
