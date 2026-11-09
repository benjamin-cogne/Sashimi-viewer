/**
 * Alignment files read from a URL, next to the files read from this computer.
 *
 * A sample's file is either a File (picked, dropped, or found in a run folder: read on demand with File.slice, as
 * always) or a RemoteSource: plain data naming a URL, read with HTTP Range requests, so only the header, the index
 * and the blocks of the windows looked at travel, never the whole file. RemoteSource is plain data on purpose: it
 * crosses to the variant worker by postMessage like a File does, and the worker opens its own handle on it.
 *
 * What the server has to do (as for igv.js): answer Range requests with 206 Partial Content and, when the page is not
 * served from the same origin as the files, send CORS headers allowing the page's origin and the Range header and
 * exposing Content-Range. A page served from the same server as the BAMs needs no CORS header at all.
 */
import { BlobFile, RemoteFile, type FilehandleOptions, type GenericFilehandle } from 'generic-filehandle2';

/** A file at a URL: what a sample keeps of it (plain data, cloneable). */
export interface RemoteSource {
  readonly kind: 'url';
  /** the URL as given (resolved against the page), query string kept: pre-signed URLs carry their signature there */
  readonly url: string;
  /** the file name, from the URL's path */
  readonly name: string;
  /**
   * bytes, from the server's Content-Range, else the Content-Length of a HEAD; Infinity when the server tells neither
   * to this page (reads then stop at the end of the file by themselves: a short or empty answer)
   */
  readonly size: number;
}

/** What a sample's data or index file is: a file of this computer, or a file at a URL. */
export type SampleFile = File | RemoteSource;

export const isRemote = (f: SampleFile | null | undefined): f is RemoteSource =>
  !!f && !(typeof Blob !== 'undefined' && f instanceof Blob) && (f as RemoteSource).kind === 'url';

/** The same file: the same File object, or the same URL. */
export const sameSource = (a: SampleFile, b: SampleFile): boolean => a === b || (isRemote(a) && isRemote(b) && a.url === b.url);

/** File name of a URL: the last segment of its path, decoded. */
export function urlFileName(url: string): string {
  try {
    const seg = new URL(url, base()).pathname.split('/').filter(Boolean).pop() ?? '';
    return decodeURIComponent(seg) || url;
  } catch { return url; }
}

/** The URL with a suffix added to its path, before any query or fragment (x.bam?sig=… → x.bam.bai?sig=…). */
export function withPathSuffix(url: string, suffix: string, strip?: RegExp): string {
  const u = new URL(url, base());
  u.pathname = (strip ? u.pathname.replace(strip, '') : u.pathname) + suffix;
  return u.toString();
}

const base = () => (typeof location !== 'undefined' && /^https?:/.test(location.href) ? location.href : undefined);

// ---------------- Reading ----------------

/** Bytes per cached block: reads are rounded out to whole blocks, so neighbouring reads cost one request. */
const BLOCK = 128 * 1024;
/** Blocks kept per file (least recently used dropped first): 32 MB. */
const MAX_BLOCKS = 256;
/** A read this long or longer goes straight to the server, uncached (an index read whole, a long tile). */
const DIRECT = 8 * 1024 * 1024;

/**
 * RemoteFile with a block cache. The BAM and CRAM readers ask for many small, overlapping ranges (the header grown by
 * steps, a BGZF block per tile, tiles that start inside the previous one); over a network each would be a round trip.
 * Here reads are rounded out to 128 kB blocks, the missing blocks of a read are fetched in one request per run of
 * consecutive blocks, the blocks are kept (32 MB per file), and concurrent reads of the same block share one request.
 */
class CachedRemoteFile extends RemoteFile {
  private blocks = new Map<number, Uint8Array>();
  private inflight = new Map<number, Promise<Uint8Array>>();

  protected async fetchBytes(length: number, position: number, opts: FilehandleOptions): Promise<Uint8Array<ArrayBuffer>> {
    if (length >= DIRECT) return super.fetchBytes(length, position, opts);
    const first = Math.floor(position / BLOCK), last = Math.floor((position + length - 1) / BLOCK);
    const want: Promise<Uint8Array>[] = [];
    let run: number[] = [];
    const flush = () => {
      if (!run.length) return;
      const ids = run; run = [];
      // fetched without the caller's signal: other reads may be waiting on the same blocks
      const p = super.fetchBytes((ids.length) * BLOCK, ids[0] * BLOCK, { ...opts, signal: undefined });
      ids.forEach((id, k) => {
        const b = p.then(bytes => {
          const blk = bytes.subarray(k * BLOCK, Math.min(bytes.length, (k + 1) * BLOCK)).slice();
          this.keep(id, blk);
          return blk;
        }).finally(() => this.inflight.delete(id));
        this.inflight.set(id, b);
      });
    };
    for (let id = first; id <= last; id++) {
      const have = this.blocks.get(id);
      if (have) { this.blocks.delete(id); this.blocks.set(id, have); }   // most recently used last
      const pending = have ? null : this.inflight.get(id);
      if (have || pending) flush(); else run.push(id);
    }
    flush();
    for (let id = first; id <= last; id++) want.push(this.blocks.get(id) ? Promise.resolve(this.blocks.get(id)!) : this.inflight.get(id)!);
    const parts = await abortable(Promise.all(want), opts.signal);
    const out = new Uint8Array(length);
    let filled = 0;
    for (let k = 0; k < parts.length; k++) {
      const blkStart = (first + k) * BLOCK;
      const from = Math.max(0, position - blkStart), to = Math.min(parts[k].length, position + length - blkStart);
      if (to <= from) break;   // past the end of the file
      out.set(parts[k].subarray(from, to), filled);
      filled += to - from;
      if (parts[k].length < BLOCK) break;   // the last block of the file
    }
    return filled === length ? out : out.slice(0, filled);
  }

  private keep(id: number, blk: Uint8Array) {
    this.blocks.set(id, blk);
    while (this.blocks.size > MAX_BLOCKS) this.blocks.delete(this.blocks.keys().next().value!);
  }

  clear() { this.blocks.clear(); }
}

function abortable<T>(p: Promise<T>, signal?: AbortSignal): Promise<T> {
  if (!signal) return p;
  if (signal.aborted) return Promise.reject(signal.reason ?? new DOMException('Aborted', 'AbortError'));
  return new Promise<T>((resolve, reject) => {
    const onAbort = () => reject(signal.reason ?? new DOMException('Aborted', 'AbortError'));
    signal.addEventListener('abort', onAbort, { once: true });
    p.then(v => { signal.removeEventListener('abort', onAbort); resolve(v); }, e => { signal.removeEventListener('abort', onAbort); reject(e); });
  });
}

/** One handle per URL in this thread (page or worker), so the header read and the readers share its blocks. */
const remoteHandles = new Map<string, CachedRemoteFile>();

/** A file handle for the GMOD readers (BAM, CRAM, FASTA): a File read on demand, or a URL read by ranges. */
export function handleOf(f: SampleFile): GenericFilehandle {
  if (!isRemote(f)) return new BlobFile(f);
  let h = remoteHandles.get(f.url);
  if (!h) { h = new CachedRemoteFile(f.url); remoteHandles.set(f.url, h); }
  return h;
}

/** Bytes [offset, offset + length) of a file, fewer at its end. */
export async function readBytes(f: SampleFile, offset: number, length: number, signal?: AbortSignal): Promise<Uint8Array> {
  if (!isRemote(f)) return new Uint8Array(await f.slice(offset, Math.min(offset + length, f.size)).arrayBuffer());
  const len = Math.max(0, Math.min(length, f.size - offset));
  return len ? handleOf(f).read(len, offset, { signal }) : new Uint8Array(0);
}

/** Forget the cached blocks of a URL (the sample was removed). */
export function forgetRemote(f: SampleFile | null | undefined) {
  if (isRemote(f)) { remoteHandles.get(f.url)?.clear(); remoteHandles.delete(f.url); }
}

// ---------------- Opening ----------------

/**
 * Checks that a URL can be read by ranges and returns it as a RemoteSource. One request for its first byte must come
 * back as 206 Partial Content; everything else is turned into the reason a person can act on. A server that ignores
 * Range answers 200 with the whole file: the request is cancelled as soon as the headers arrive, so a 30 GB BAM is not
 * downloaded by a probe.
 *
 * The size comes from the 206's Content-Range. Across origins the page may read that header only when the server lists
 * it in Access-Control-Expose-Headers, which many servers configured for CORS do not; the size is then taken from the
 * Content-Length of a HEAD request (a header every page may read), and when even that is missing it is left unknown:
 * BAM and CRAM reads go by the index and do not need it.
 */
export async function probeUrl(url: string, signal?: AbortSignal): Promise<RemoteSource> {
  let abs: string;
  try { abs = new URL(url.trim(), base()).toString(); } catch { throw new Error(`${url} is not a URL`); }
  if (!/^https?:/i.test(abs)) throw new Error(`${url}: only http:// and https:// URLs can be read`);
  const name = urlFileName(abs);
  const ctl = new AbortController();
  const onAbort = () => ctl.abort();
  signal?.addEventListener('abort', onAbort, { once: true });
  let res: Response;
  try {
    res = await fetch(abs, { headers: { range: 'bytes=0-0' }, signal: ctl.signal, cache: 'no-store' });
  } catch (e: any) {
    signal?.removeEventListener('abort', onAbort);
    if (e?.name === 'AbortError') throw e;
    throw new Error(`${name}: the server could not be reached, or it refused this page (no CORS header allowing it: Access-Control-Allow-Origin, with Range in Access-Control-Allow-Headers). Serving this page from the same server as the files avoids CORS altogether.`);
  }
  signal?.removeEventListener('abort', onAbort);
  if (res.status === 200) { ctl.abort(); throw new Error(`${name}: the server ignores Range requests (it answered with the whole file), so the file cannot be read in parts`); }
  if (res.status === 401 || res.status === 403) { ctl.abort(); throw new Error(`${name}: access denied (HTTP ${res.status})`); }
  if (res.status === 404) { ctl.abort(); throw new Error(`${name}: not found (HTTP 404)`); }
  if (res.status !== 206) { ctl.abort(); throw new Error(`${name}: HTTP ${res.status}`); }
  const enc = res.headers.get('content-encoding');
  if (enc && enc !== 'identity') { ctl.abort(); throw new Error(`${name}: the server compresses the file on the fly (Content-Encoding: ${enc}), which breaks reading it in parts: serve .bam/.bai/.cram/.crai as application/octet-stream without gzip`); }
  const cr = res.headers.get('content-range');
  await res.arrayBuffer().catch(() => {});
  let size = cr ? Number(/\/(\d+)\s*$/.exec(cr)?.[1]) : NaN;
  if (!(size > 0)) size = await headLength(abs, signal);
  return { kind: 'url', url: abs, name, size: size > 0 ? size : Infinity };
}

/** Content-Length of a HEAD request (readable across origins without being exposed), or NaN. */
async function headLength(url: string, signal?: AbortSignal): Promise<number> {
  try {
    const res = await fetch(url, { method: 'HEAD', signal, cache: 'no-store' });
    if (!res.ok) return NaN;
    const enc = res.headers.get('content-encoding');
    if (enc && enc !== 'identity') return NaN;   // the length of a compressed body, not of the file
    return Number(res.headers.get('content-length'));
  } catch (e: any) {
    if (e?.name === 'AbortError') throw e;
    return NaN;
  }
}

/** A file size for a person: "1.18 GB", or "size not told by the server" for a URL whose size is unknown. */
export const sizeLabel = (bytes: number): string => (Number.isFinite(bytes) ? `${(bytes / 1e9).toFixed(2)} GB` : 'size not told by the server');

/** A probe that answers null instead of failing (index candidates). */
const tryProbe = (url: string, signal?: AbortSignal) => probeUrl(url, signal).catch(e => { if (e?.name === 'AbortError') throw e; return null; });

/**
 * The index of an alignment URL, as for files: `x.bam.bai` then `x.bai`, `x.cram.crai` then `x.crai` (the suffixes
 * given, for a registered kind), each beside the file. The first error of the data file itself is what is reported.
 */
export async function findIndexUrl(data: RemoteSource, suffixes: string[], stripExt: RegExp, signal?: AbortSignal): Promise<RemoteSource | null> {
  for (const s of suffixes) {
    for (const u of [withPathSuffix(data.url, s), withPathSuffix(data.url, s, stripExt)]) {
      if (u === data.url) continue;
      const got = await tryProbe(u, signal);
      if (got) return got;
    }
  }
  return null;
}
