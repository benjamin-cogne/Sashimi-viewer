/**
 * File intake for the standalone page: which local files form a sample, and how to collect them
 * from a picker, a drop or a folder.
 *
 * A web page cannot open a file by path: the browser only hands over what the user selected or
 * dropped. So "add sample.bam and find sample.bam.bai next to it" works in two ways. (1) Everything
 * the user ever added goes into one pool, and alignments are paired with their index by name
 * whatever the order of arrival. (2) The user gives the folder (picker or drop) and the page reads
 * the directory listing itself, which is what a desktop tool would do with the path.
 */

export const ALIGN_EXT = /\.(bam|cram)$/i;
export const INDEX_EXT = /\.(bai|csi|crai)$/i;
export const FASTA_EXT = /\.(fa|fasta|fna)(\.gz)?$/i;
export const FASTA_INDEX_EXT = /\.(fai|gzi)$/i;
/** Every extension the page has a use for; other files in a folder are ignored without being read. */
export const WANTED_EXT = /\.(bam|cram|bai|csi|crai|fa|fasta|fna|fa\.gz|fasta\.gz|fna\.gz|fai|gzi)$/i;

export interface FileLike { name: string }
export type AlignKind = 'bam' | 'cram';

export interface PairedSample<F extends FileLike = File> { name: string; kind: AlignKind; file: F; index: F }
export interface PendingFile<F extends FileLike = File> {
  file: F;
  /** Index file names that would complete it, most conventional first (samtools default first). */
  wanted: string[];
}
export interface PairedFasta<F extends FileLike = File> { fa: F; fai: F; gzi?: F }

export function alignKind(name: string): AlignKind | null {
  const m = ALIGN_EXT.exec(name);
  return m ? (m[1].toLowerCase() as AlignKind) : null;
}

/** Candidate index names for an alignment: `x.bam.bai`, `x.bai`, `x.bam.csi`, `x.csi`; `x.cram.crai`, `x.crai`. */
export function indexCandidates(name: string): string[] {
  const kind = alignKind(name);
  if (!kind) return [];
  const stem = name.replace(ALIGN_EXT, '');
  const exts = kind === 'bam' ? ['bai', 'csi'] : ['crai'];
  return exts.flatMap(ext => [`${name}.${ext}`, `${stem}.${ext}`]);
}

/** Index names for a FASTA: `x.fa.fai` always, plus `x.fa.gz.gzi` when bgzipped. */
export function fastaIndexCandidates(name: string): { fai: string; gzi?: string } {
  return { fai: `${name}.fai`, gzi: /\.gz$/i.test(name) ? `${name}.gzi` : undefined };
}

/**
 * Resolve a pool of files (keyed by lower-cased name) into samples, alignments still waiting for an
 * index, and the reference FASTA. Pure: the caller decides which samples are new.
 */
export function pairPool<F extends FileLike>(pool: Map<string, F>): { samples: PairedSample<F>[]; pending: PendingFile<F>[]; fasta?: PairedFasta<F>; fastaPending?: PendingFile<F> } {
  const lookup = (n: string) => pool.get(n.toLowerCase());
  const samples: PairedSample<F>[] = [];
  const pending: PendingFile<F>[] = [];
  for (const f of pool.values()) {
    const kind = alignKind(f.name);
    if (!kind) continue;
    const wanted = indexCandidates(f.name);
    const idx = wanted.map(lookup).find(Boolean);
    if (idx) samples.push({ name: f.name.replace(ALIGN_EXT, ''), kind, file: f, index: idx });
    else pending.push({ file: f, wanted });
  }
  samples.sort((a, b) => a.name.localeCompare(b.name));
  pending.sort((a, b) => a.file.name.localeCompare(b.file.name));
  let fasta: PairedFasta<F> | undefined, fastaPending: PendingFile<F> | undefined;
  const fa = [...pool.values()].find(f => FASTA_EXT.test(f.name));
  if (fa) {
    const want = fastaIndexCandidates(fa.name);
    const fai = lookup(want.fai), gzi = want.gzi ? lookup(want.gzi) : undefined;
    if (fai && (!want.gzi || gzi)) fasta = { fa, fai, gzi };
    else fastaPending = { file: fa, wanted: [want.fai, want.gzi].filter((n): n is string => !!n && !lookup(n)) };
  }
  return { samples, pending, fasta, fastaPending };
}

// ---- folders ----

export const MAX_FOLDER_DEPTH = 3;
export const MAX_FOLDER_FILES = 5000;

/** Files of a folder chosen with the File System Access API (Chromium: Edge, Chrome). Only wanted extensions are materialised. */
export async function filesFromDirectoryHandle(dir: any, depth = 0, out: File[] = []): Promise<File[]> {
  for await (const [name, handle] of dir.entries()) {
    if (out.length >= MAX_FOLDER_FILES) break;
    if (handle.kind === 'file') { if (WANTED_EXT.test(name)) out.push(await handle.getFile()); }
    else if (handle.kind === 'directory' && depth < MAX_FOLDER_DEPTH && !name.startsWith('.')) await filesFromDirectoryHandle(handle, depth + 1, out);
  }
  return out;
}

/** Files behind the items of a drop: plain files, plus the contents of dropped folders (webkitGetAsEntry, all current browsers). */
export async function filesFromDrop(dt: DataTransfer): Promise<File[]> {
  const items = dt.items ? Array.from(dt.items) : [];
  const entries = items.map(it => (typeof it.webkitGetAsEntry === 'function' ? it.webkitGetAsEntry() : null));
  if (!entries.some(Boolean)) return Array.from(dt.files);
  const out: File[] = [];
  const readAll = (reader: any): Promise<any[]> => new Promise((resolve, reject) => {
    const acc: any[] = [];
    const step = () => reader.readEntries((batch: any[]) => { if (!batch.length) resolve(acc); else { acc.push(...batch); step(); } }, reject);
    step();
  });
  const walk = async (entry: any, depth: number) => {
    if (out.length >= MAX_FOLDER_FILES) return;
    if (entry.isFile) {
      if (depth === 0 || WANTED_EXT.test(entry.name)) out.push(await new Promise<File>((res, rej) => entry.file(res, rej)));
    } else if (entry.isDirectory && depth < MAX_FOLDER_DEPTH && !entry.name.startsWith('.')) {
      for (const child of await readAll(entry.createReader())) await walk(child, depth + 1);
    }
  };
  for (let i = 0; i < entries.length; i++) {
    const e = entries[i];
    if (e) await walk(e, 0);
    else { const f = items[i].getAsFile(); if (f) out.push(f); }
  }
  return out;
}

export function supportsDirectoryPicker(): boolean {
  return typeof (globalThis as any).showDirectoryPicker === 'function';
}
