/**
 * Getting the files of a session back without asking for each of them again.
 *
 * A page never sees a file's path, so a session names its files relative to a *run folder*. In
 * every browser the user can drop or pick that folder, and the page finds each file by its
 * relative path (only metadata is enumerated; a BAM is read on demand as always). Chromium
 * browsers (Chrome, Edge, Opera) also hand over FileSystemHandle objects for picked or dropped
 * files and folders: a handle is a small bookmark, not a copy, and kept in IndexedDB it reopens
 * the folder in a later page load after one permission click. Firefox and Safari expose no
 * handles, so there the folder is dropped or picked again.
 */
import { kindExtensions } from './fileKinds';

/** Minimal typing of the WICG permission methods, which lib.dom does not declare. */
type Permissioned = { queryPermission?(o: { mode: 'read' | 'readwrite' }): Promise<PermissionState>; requestPermission?(o: { mode: 'read' | 'readwrite' }): Promise<PermissionState> };
export type FSHandle = FileSystemFileHandle & Permissioned;
export type FSDirHandle = FileSystemDirectoryHandle & Permissioned;

export const hasFileSystemAccess = (): boolean => typeof (window as any).showDirectoryPicker === 'function' && typeof indexedDB !== 'undefined';

/** A file with its path relative to the run folder ("run42/patient.bam"), when it came from a folder. */
export interface PathedFile { file: File; path?: string }

// ---------------- IndexedDB bookmarks (Chromium) ----------------

const DB = 'sashimi-viewer', STORE = 'handles';

function openDb(): Promise<IDBDatabase> {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(DB, 1);
    req.onupgradeneeded = () => { if (!req.result.objectStoreNames.contains(STORE)) req.result.createObjectStore(STORE, { keyPath: 'key' }); };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}
async function put(rows: { key: string; handle: FileSystemHandle; size?: number }[]): Promise<void> {
  if (!rows.length) return;
  const db = await openDb();
  await new Promise<void>((resolve, reject) => {
    const tx = db.transaction(STORE, 'readwrite');
    for (const r of rows) tx.objectStore(STORE).put({ ...r, saved: Date.now() });
    tx.oncomplete = () => resolve(); tx.onerror = () => reject(tx.error);
  });
  db.close();
}
async function get(key: string): Promise<any | null> {
  try {
    const db = await openDb();
    const row = await new Promise<any>((resolve, reject) => {
      const req = db.transaction(STORE, 'readonly').objectStore(STORE).get(key);
      req.onsuccess = () => resolve(req.result ?? null); req.onerror = () => reject(req.error);
    });
    db.close();
    return row;
  } catch { return null; }
}

export const rememberFolder = (name: string, handle: FSDirHandle) => put([{ key: `folder:${name}`, handle }]).catch(() => {});
export const recallFolder = async (name: string): Promise<FSDirHandle | null> => (await get(`folder:${name}`))?.handle ?? null;
export const rememberFiles = (entries: { name: string; size: number; handle: FSHandle }[]) => put(entries.map(e => ({ key: `file:${e.name}`, handle: e.handle, size: e.size }))).catch(() => {});
export const recallFile = async (name: string): Promise<FSHandle | null> => (await get(`file:${name}`))?.handle ?? null;

/** Read permission on a handle; `ask` lets the browser prompt (needs a user gesture). */
export async function permitted(h: Permissioned, ask: boolean): Promise<boolean> {
  try {
    let state: PermissionState = 'granted';
    if (h.queryPermission) state = await h.queryPermission({ mode: 'read' });
    if (state !== 'granted' && ask && h.requestPermission) state = await h.requestPermission({ mode: 'read' });
    return state === 'granted';
  } catch { return false; }
}

// ---------------- Pickers (Chromium) ----------------

export async function pickFolder(): Promise<FSDirHandle> {
  return (window as any).showDirectoryPicker({ mode: 'read' });
}
export async function pickFiles(): Promise<{ file: File; handle: FSHandle }[]> {
  const handles: FSHandle[] = await (window as any).showOpenFilePicker({
    multiple: true,
    types: [{ description: 'Alignments, indexes, reference', accept: { 'application/octet-stream': ['.bam', '.bai', '.cram', '.crai', ...kindExtensions(), '.fa', '.fasta', '.fna', '.gz', '.fai', '.gzi'] } }],
  });
  return Promise.all(handles.map(async h => ({ file: await h.getFile(), handle: h })));
}

const KEEP = /\.(bam|bai|cram|crai|fa|fasta|fna|fai|gzi|gz)$/i;
/** a file the viewer reads: alignments, indexes, reference, and the files of registered kinds (asked each time: plugins register at load) */
const keep = (name: string) => KEEP.test(name) || kindExtensions().some(x => name.toLowerCase().endsWith(x));
const MAX_DEPTH = 4, MAX_FILES = 5000;

/** Alignment-related files of a folder handle with their relative paths (metadata only; bounded depth and count). */
export async function filesInFolder(dir: FSDirHandle): Promise<PathedFile[]> {
  const out: PathedFile[] = [];
  const walk = async (d: FileSystemDirectoryHandle, prefix: string, depth: number) => {
    for await (const [name, h] of (d as any).entries() as AsyncIterable<[string, FileSystemHandle]>) {
      if (out.length >= MAX_FILES) return;
      if (h.kind === 'directory') { if (depth < MAX_DEPTH) await walk(h as FileSystemDirectoryHandle, `${prefix}${name}/`, depth + 1); }
      else if (keep(name)) out.push({ file: await (h as FileSystemFileHandle).getFile(), path: `${prefix}${name}` });
    }
  };
  await walk(dir, '', 0);
  return out;
}

/** One file of a folder handle by relative path, or null. */
export async function fileInFolder(dir: FSDirHandle, relPath: string): Promise<File | null> {
  try {
    const parts = relPath.split('/').filter(Boolean);
    let d: FileSystemDirectoryHandle = dir;
    for (const p of parts.slice(0, -1)) d = await d.getDirectoryHandle(p);
    return await (await d.getFileHandle(parts[parts.length - 1])).getFile();
  } catch { return null; }
}

// ---------------- Folder input and drops (every browser) ----------------

/** Files of a `webkitdirectory` input: the folder name and paths inside it. */
export function filesFromFolderInput(list: FileList): { folder: string | null; files: PathedFile[] } {
  const files: PathedFile[] = [];
  let folder: string | null = null;
  for (const f of Array.from(list)) {
    const rel = (f as any).webkitRelativePath as string | undefined;
    if (rel && rel.includes('/')) {
      const [top, ...rest] = rel.split('/');
      folder = folder ?? top;
      if (keep(f.name)) files.push({ file: f, path: rest.join('/') });
    } else if (keep(f.name)) files.push({ file: f });
  }
  return { folder, files };
}

/**
 * Files of a drop: plain files as they are; a dropped folder is walked (every browser, through
 * webkitGetAsEntry) and its files carry their relative paths. Chromium also gives the handles.
 */
export async function filesFromDrop(dt: DataTransfer): Promise<{ folder: string | null; folderHandle: FSDirHandle | null; files: PathedFile[]; fileHandles: Map<string, FSHandle>; /** dropped files that are not alignments: session files (.json) */ others: File[] }> {
  const files: PathedFile[] = [];
  const fileHandles = new Map<string, FSHandle>();
  const others: File[] = [];
  let folder: string | null = null, folderHandle: FSDirHandle | null = null;
  // Everything the DataTransfer can give is taken now, synchronously: browsers (Chrome first) empty the item list
  // and `dt.files` as soon as the drop handler yields, so a getAsFile() after an await returns null.
  const items = Array.from(dt.items || []);
  const entries = items.map(it => (typeof it.webkitGetAsEntry === 'function' ? it.webkitGetAsEntry() : null));
  const plain = items.map(it => (it.kind === 'file' ? it.getAsFile() : null));
  const fallback = Array.from(dt.files || []);
  const handlePromises = items.map(it => (typeof (it as any).getAsFileSystemHandle === 'function' ? (it as any).getAsFileSystemHandle().catch(() => null) : Promise.resolve(null)));
  const readDir = (d: FileSystemDirectoryEntry): Promise<FileSystemEntry[]> => new Promise((resolve, reject) => {
    const reader = d.createReader(); const all: FileSystemEntry[] = [];
    const step = () => reader.readEntries(batch => { if (!batch.length) resolve(all); else { all.push(...batch); step(); } }, reject);
    step();
  });
  const toFile = (e: FileSystemFileEntry): Promise<File> => new Promise((resolve, reject) => e.file(resolve, reject));
  const walk = async (e: FileSystemEntry, prefix: string, depth: number) => {
    if (files.length >= MAX_FILES) return;
    if (e.isDirectory) { if (depth < MAX_DEPTH) for (const c of await readDir(e as FileSystemDirectoryEntry)) await walk(c, `${prefix}${e.name}/`, depth + 1); }
    else if (keep(e.name)) files.push({ file: await toFile(e as FileSystemFileEntry), path: `${prefix}${e.name}` });
  };
  const handles = await Promise.all(handlePromises);
  for (let i = 0; i < items.length; i++) {
    const e = entries[i], h = handles[i];
    if (e?.isDirectory) {
      folder = folder ?? e.name;
      if (h && h.kind === 'directory' && !folderHandle) folderHandle = h as FSDirHandle;
      // paths are relative to the dropped folder itself
      for (const c of await readDir(e as FileSystemDirectoryEntry)) await walk(c, '', 1);
    } else {
      const f = plain[i];
      if (!f) continue;
      if (keep(f.name)) { files.push({ file: f }); if (h && h.kind === 'file') fileHandles.set(f.name, h as FSHandle); }
      else others.push(f);
    }
  }
  if (!items.length) for (const f of fallback) { if (keep(f.name)) files.push({ file: f }); else others.push(f); }
  return { folder, folderHandle, files, fileHandles, others };
}
