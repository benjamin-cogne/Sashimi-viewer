/**
 * Reading annotation files as a stream of lines, and the name index of a tabix-indexed one (annotation.ts).
 *
 * Lines: plain text, or gzipped, inflated with fflate as the bytes come: a bgzipped file (tabix's BGZF, many gzip
 * members one after the other) is inflated member after member, which the browsers' DecompressionStream does not do.
 *
 * Name index: an indexed file is read by region (its genes, when a window is drawn), so finding a gene by its name
 * needs, once, a pass over the whole file. It keeps only what a lookup needs (each gene's name, contig and span, and
 * the names it is known by: gene ID, transcript IDs and names), a few MB for a whole human annotation, not the models.
 * It runs in a worker (annotationWorker.ts) or, without one, here; no import of the page, so that a worker can hold it.
 */
import { Gunzip } from 'fflate';

export type AnnotationFormat = 'gff3' | 'gtf' | 'bed';

/** Batches of the file's lines (no line ends), inflated as a stream when gzipped, decoded as UTF-8. */
export async function* readLineBatches(file: Blob, gz: boolean, onBytes?: (read: number) => void): AsyncGenerator<string[]> {
  const reader = file.stream().getReader();
  const dec = new TextDecoder();
  const text: string[] = [];
  const gunzip = gz ? new Gunzip(d => { text.push(dec.decode(d, { stream: true })); }) : null;
  let rest = '', read = 0;
  try {
    for (;;) {
      const { value, done } = await reader.read();
      if (value) {
        read += value.length; onBytes?.(read);
        if (gunzip) gunzip.push(value, false); else text.push(dec.decode(value, { stream: true }));
      }
      if (done) { if (gunzip) gunzip.push(new Uint8Array(0), true); text.push(dec.decode()); }
      if (text.length) {
        const lines = (rest + text.join('')).split('\n');
        text.length = 0;
        rest = lines.pop() ?? '';
        for (let i = 0; i < lines.length; i++) if (lines[i].endsWith('\r')) lines[i] = lines[i].slice(0, -1);
        if (lines.length) yield lines;
      }
      if (done) break;
    }
  } finally { reader.releaseLock(); }
  if (rest) yield [rest.endsWith('\r') ? rest.slice(0, -1) : rest];
}

/** The file's lines one by one (readLineBatches). */
export async function* readLines(file: Blob, gz: boolean): AsyncGenerator<string> {
  for await (const batch of readLineBatches(file, gz)) yield* batch;
}

export async function isGzipped(file: Blob, name: string): Promise<boolean> {
  if (/\.gz$/i.test(name)) return true;
  const b = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return b[0] === 0x1f && b[1] === 0x8b;
}

/** The format from the name; a `.gff` is sniffed (key=value attributes: GFF3; key "value"; : GTF). */
export async function formatOf(file: Blob, name: string, gz: boolean): Promise<AnnotationFormat> {
  const n = name.replace(/\.gz$/i, '').toLowerCase();
  if (n.endsWith('.gtf')) return 'gtf';
  if (n.endsWith('.bed')) return 'bed';
  if (n.endsWith('.gff3')) return 'gff3';
  let k = 0;
  for await (const line of readLines(file, gz)) {
    if (!line || line.startsWith('#')) continue;
    const col = line.split('\t')[8] ?? '';
    if (/^\s*\w+\s+"/.test(col)) return 'gtf';
    if (/=/.test(col)) return 'gff3';
    if (++k > 50) break;
  }
  return 'gff3';
}

// ======================== the name index ========================

const enc = new TextEncoder(), dec8 = new TextDecoder();
/**
 * A copy of a string that holds nothing else. A JavaScript engine may keep a piece cut from a longer string as a view
 * of it (V8 does from 13 characters): a name cut from a line, itself cut from the text inflated from a block of the
 * file, would keep that whole block alive, and a name index the whole file.
 */
export const own = (s: string): string => (s.length < 13 ? s : dec8.decode(enc.encode(s)));

/** Where a gene is: its contig (the file's spelling) and span, 0-based half-open. */
export interface GeneLoc { name: string; chrom: string; start: number; end: number }
export interface NameIndex {
  locs: GeneLoc[];
  /** every name a gene is known by, upper case, with the gene (index in `locs`) */
  keys: [string, number][];
  lines: number; genes: number; transcripts: number;
}

/** feature types a gene is not found by (GFF3): leaves of a transcript, or spans of whole sequences */
export const LEAF_TYPE = /^(exon|pseudogenic_exon|noncoding_exon|CDS|five_prime_UTR|three_prime_UTR|UTR|5'UTR|3'UTR|five_prime_utr|three_prime_utr|start_codon|stop_codon|intron|polyA_site|polyA_signal_sequence)$/;
export const SKIP_LONE = /^(region|chromosome|contig|supercontig|scaffold|source|databank_entry|sequence_assembly|match|match_part|cDNA_match|biological_region)$/i;

/** A GTF attribute's value (key "value";), or undefined. */
function gtfAttr(attrs: string, key: string): string | undefined {
  let from = 0;
  for (;;) {
    const i = attrs.indexOf(key, from);
    if (i < 0) return undefined;
    from = i + key.length;
    if (i > 0 && attrs[i - 1] !== ' ' && attrs[i - 1] !== ';' && attrs[i - 1] !== '\t') continue;
    let j = from;
    while (attrs[j] === ' ') j++;
    if (attrs[j] === '"') { const k = attrs.indexOf('"', j + 1); return k < 0 ? undefined : attrs.slice(j + 1, k); }
    if (j === from) continue;   // a longer key starting the same way
    const k = attrs.indexOf(';', j);
    return attrs.slice(j, k < 0 ? undefined : k).trim();
  }
}
const decode = (s: string) => { try { return decodeURIComponent(s); } catch { return s; } };
function gffAttrMap(col: string): Map<string, string> {
  const m = new Map<string, string>();
  for (const kv of col.split(';')) { const i = kv.indexOf('='); if (i > 0) m.set(kv.slice(0, i).trim(), decode(kv.slice(i + 1).trim())); }
  return m;
}
/** the 9 columns' starts: positions of the first 8 tabs, or null for a line with fewer */
function tabs(line: string): number[] | null {
  const t: number[] = [];
  let i = -1;
  for (let k = 0; k < 8; k++) { i = line.indexOf('\t', i + 1); if (i < 0) return null; t.push(i); }
  return t;
}

class Builder {
  locs: GeneLoc[] = [];
  keys: [string, number][] = [];
  private last = new Map<number, Set<string>>();
  private chroms = new Map<string, string>();
  add(name: string, chrom: string, start: number, end: number): number {
    let c = this.chroms.get(chrom);
    if (c == null) this.chroms.set(chrom, c = own(chrom));
    this.locs.push({ name: own(name), chrom: c, start, end }); this.alias(this.locs.length - 1, name); return this.locs.length - 1;
  }
  alias(i: number, k: string | undefined): void {
    if (!k) return;
    const u = k.toUpperCase();
    let s = this.last.get(i);
    if (!s) this.last.set(i, s = new Set());
    if (!s.has(u)) { s.add(u); this.keys.push([own(u), i]); }
  }
}

/**
 * The names of a file's genes and where each is, from one pass over it. `onProgress` gets the fraction of the file
 * read (its compressed bytes).
 */
export async function scanNames(file: Blob, gz: boolean, format: AnnotationFormat, onProgress?: (fraction: number) => void): Promise<NameIndex> {
  const b = new Builder();
  let lines = 0, transcripts = 0, lastTick = 0;
  const tick = (read: number) => { const f = read / Math.max(1, file.size); if (f - lastTick >= 0.01) { lastTick = f; onProgress?.(f); } };
  if (format === 'gtf') {
    const genes = new Map<string, number>(), txSeen = new Set<string>();   // keys only while reading
    let sawTx = false;
    for await (const batch of readLineBatches(file, gz, tick)) for (const line of batch) {
      lines++;
      if (!line || line.charCodeAt(0) === 35) continue;
      const t = tabs(line);
      if (!t) continue;
      const type = line.slice(t[1] + 1, t[2]);
      // files with transcript rows give each transcript's span there: their exon, CDS… rows are passed over
      if (type === 'transcript') sawTx = true;
      else if (type !== 'gene' && sawTx) continue;
      const chrom = line.slice(0, t[0]), s = +line.slice(t[2] + 1, t[3]) - 1, e = +line.slice(t[3] + 1, t[4]);
      if (!(s >= 0) || !(e > s)) continue;
      const attrs = line.slice(t[7] + 1);
      const tid = type === 'gene' ? undefined : gtfAttr(attrs, 'transcript_id');
      const gid = gtfAttr(attrs, 'gene_id') ?? gtfAttr(attrs, 'gene_name') ?? tid;
      if (!gid) continue;
      const key = `${gid}\t${chrom}`;
      let i = genes.get(key);
      const gname = gtfAttr(attrs, 'gene_name');
      if (i == null) { i = b.add(gname ?? gid, chrom, s, e); genes.set(key, i); b.alias(i, gid); }
      else { const l = b.locs[i]; if (s < l.start) l.start = s; if (e > l.end) l.end = e; if (gname && l.name === gid) { l.name = own(gname); b.alias(i, gname); } }
      if (tid) {
        b.alias(i, tid); b.alias(i, gtfAttr(attrs, 'transcript_name'));
        const tk = `${tid}\t${chrom}`;
        if (!txSeen.has(tk)) { txSeen.add(tk); transcripts++; }
      }
    }
  } else if (format === 'gff3') {
    interface F { id: string | null; parent: string | null; name: string; chrom: string; start: number; end: number; names: string[] }
    const byId = new Map<string, F>(), all: F[] = [];
    for await (const batch of readLineBatches(file, gz, tick)) {
      let stop = false;
      for (const line of batch) {
        lines++;
        if (!line) continue;
        if (line.charCodeAt(0) === 35) { if (line.startsWith('##FASTA')) { stop = true; break; } continue; }
        const t = tabs(line);
        if (!t) continue;
        const type = line.slice(t[1] + 1, t[2]);
        if (LEAF_TYPE.test(type) || SKIP_LONE.test(type)) continue;
        const s = +line.slice(t[2] + 1, t[3]) - 1, e = +line.slice(t[3] + 1, t[4]);
        if (!(s >= 0) || !(e > s)) continue;
        const a = gffAttrMap(line.slice(t[7] + 1));
        const id = a.get('ID') ?? null, chrom = decode(line.slice(0, t[0]));
        const name = a.get('Name') ?? a.get('gene_name') ?? a.get('gene') ?? a.get('locus_tag') ?? id ?? `${chrom}:${s + 1}-${e}`;
        const f: F = { id, parent: a.get('Parent')?.split(',')[0] ?? null, name, chrom, start: s, end: e, names: [name, id ?? '', ...['transcript_id', 'gene_id', 'locus_tag', 'gene', 'Alias', 'gene_synonym', 'old_locus_tag'].flatMap(k => (a.get(k) ?? '').split(','))].map(x => x.trim()).filter(Boolean) };
        all.push(f);
        if (id && !byId.has(id)) byId.set(id, f);
      }
      if (stop) break;
    }
    const rootOf = (f: F): F => { let r = f; for (let d = 0; d < 8 && r.parent; d++) { const p = byId.get(r.parent); if (!p) break; r = p; } return r; };
    const at = new Map<F, number>();
    for (const f of all) {
      const r = rootOf(f);
      let i = at.get(r);
      if (i == null) { i = b.add(r.name, r.chrom, r.start, r.end); at.set(r, i); for (const n of r.names) b.alias(i, n); }
      if (r !== f) { transcripts++; for (const n of f.names) b.alias(i, n); const l = b.locs[i]; if (f.start < l.start) l.start = f.start; if (f.end > l.end) l.end = f.end; }
    }
  } else {
    for await (const batch of readLineBatches(file, gz, tick)) for (const line of batch) {
      lines++;
      if (!line.trim() || /^(track|browser)\b/.test(line) || line.charCodeAt(0) === 35) continue;
      const c = line.split('\t');
      const s = +c[1], e = +c[2];
      if (!(s >= 0) || !(e > s)) continue;
      b.add(c[3] && c[3] !== '.' ? c[3] : `${c[0]}:${s + 1}-${e}`, c[0], s, e);
      transcripts++;
    }
  }
  onProgress?.(1);
  return { locs: b.locs, keys: b.keys, lines, genes: b.locs.length, transcripts };
}
