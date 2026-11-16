/**
 * The array-CGH samples of the SV browser, kept for the tab's life (the SV browser can be closed and opened again):
 * the files read, the class given to each call, the filter, the call shown. Nothing of a patient is written to the
 * browser's storage: only the filter here (the lab's benign list and the column mappings in their own modules).
 */
import { useSyncExternalStore } from 'react';
import { unzipSync } from 'fflate';
import { parseCghFile, leafOf, toClass, type CghFile, type CghCnv, type CghSex } from './cghIntervals';
import { NeedsMapping, findMap, splitRecords, type MapSpec } from './cghMapping';
import { decodeText, gunzip, isGzip, isOle, isZip, xlsxText, OLE_MESSAGE } from './cghFormats';
import { isVcf, parseVcf } from './cghVcf';
import { seriesOf } from './cghSeries';
import { labHits, loadStoredLab, storeLab, type LabBenign } from './cghLabBenign';
import type { CnvAnn } from './cghAnnotate';
import type { GeneModes } from '../../standalone/panelapp';

export interface CghFilter {
  /** calls with fewer probes are not analysed (a BED's calls, without a count, are kept) */
  minProbes: number;
  kind: 'all' | 'loss' | 'gain';
  /** smallest size shown, kb */
  minKb: number;
  classed: 'all' | 'unclassified' | 'classified';
  hideSuppressed: boolean;
  /** a gene, a band or a chromosome typed */
  query: string;
  /** the calls in genomic order, or the highest priority first */
  sort: 'position' | 'priority';
  /** recurrence in the series: any call, in no other patient, or in at least SERIES_RECURRENT others */
  series: 'any' | 'unique' | 'recurrent';
  /** the laboratory's benign list: any call, those not in it, those in it */
  lab: 'any' | 'not' | 'in';
}
/** "recurrent": the call is found in at least this many other patients */
export const SERIES_RECURRENT = 2;
export interface CghState {
  files: CghFile[];
  /** the class given in the browser, by call id (else the file's own) */
  classes: Record<string, string>;
  /** the sample listed (its key), or null for every sample */
  sample: string | null;
  filter: CghFilter;
  /** the call drawn in the browser */
  shown: string | null;
  /** files read, waiting for the minimum number of probes; the files that could not be read */
  pending: { files: CghFile[]; errors: string[]; notes: string[] } | null;
  /** tables whose columns were not recognised (or whose columns are changed), waiting for the column mapping; `orig` is kept when it is cancelled */
  toMap: { name: string; text: string; fallback?: string; spec?: MapSpec; orig?: CghFile }[];
  /** classes changed since the last export */
  dirty: boolean;
  /** the patients (sample keys) with classes or notes not exported yet */
  dirtyKeys: Record<string, true>;
  /** the array-CGH page is open (over the SV browser) */
  open: boolean;
  /** what each call holds (genes, dosage, constraint, AChro-Puce), by call id; PanelApp's entries by gene */
  ann: Record<string, CnvAnn>;
  modes: Record<string, GeneModes | 'loading' | 'error'>;
  /** each patient's sex, by sample key (absent: not given) */
  sex: Record<string, CghSex>;
  /** each patient's DLRS (derivative log ratio spread), by sample key, for the mosaic interval (absent: not given; not persisted) */
  dlrs: Record<string, number>;
  /** a note on each call, by call id (else the file's) */
  notes: Record<string, string>;
  /** the laboratory's list of benign CNVs (a reference, kept in localStorage), and whether it could be stored */
  lab: LabBenign | null; labStored: boolean;
}

const FILTER_KEY = 'sashimi.cghFilter';
export const CGH_MIN_PROBES_DEFAULT = 3;
/** at most this many calls are annotated (genes, PanelApp, ClinGen, gnomAD, DGV): more, and a minimum must be raised */
export const CGH_MAX_ANNOTATED = 200;
/** a call is analysed (listed and annotated) when it has enough probes (a BED's calls have no count) and is large enough */
export const isAnalysed = (c: CghCnv, minProbes: number, minKb: number) => (c.probes == null || c.probes >= minProbes) && c.sizeKb >= minKb;
const loadFilter = (): CghFilter => {
  const d: CghFilter = { minProbes: CGH_MIN_PROBES_DEFAULT, kind: 'all', minKb: 0, classed: 'all', hideSuppressed: true, query: '', sort: 'position', series: 'any', lab: 'any' };
  try { const v = JSON.parse(localStorage.getItem(FILTER_KEY) ?? 'null') as Partial<CghFilter> | null; return v ? { ...d, ...v, query: '', series: 'any', lab: 'any' } : d; } catch { return d; }
};

let state: CghState = { files: [], classes: {}, sample: null, filter: loadFilter(), shown: null, pending: null, toMap: [], dirty: false, dirtyKeys: {}, open: false, ann: {}, modes: {}, sex: {}, dlrs: {}, notes: {}, lab: loadStoredLab(), labStored: true };
const subs = new Set<() => void>();
const set = (patch: Partial<CghState>) => { state = { ...state, ...patch }; subs.forEach(f => f()); };
const subscribe = (f: () => void) => { subs.add(f); return () => { subs.delete(f); }; };
export const useCgh = () => useSyncExternalStore(subscribe, () => state);

// a tab closed with classes typed and not exported asks first
const onUnload = (e: BeforeUnloadEvent) => { if (state.dirty) { e.preventDefault(); e.returnValue = ''; } };
if (typeof window !== 'undefined') window.addEventListener('beforeunload', onUnload);

/** the files a zip's entries are read from (a .xls may be a text export, or an old binary refused with a message) */
const TABLE = /\.(txt|tsv|csv|xlsx?|bed|vcf)(\.gz)?$/i;
/** the file pickers' list (kept in step with what `add` reads) */
export const CGH_ACCEPT = '.txt,.tsv,.csv,.bed,.xls,.xlsx,.vcf,.gz,.zip,text/plain,application/zip';
/** what the pickers take, in words (the page's help and the buttons' tooltips) */
export const CGH_FILES_TEXT = 'the interval table of each patient exported by CytoGenomics (.txt), any other table of calls (.txt, .tsv, .csv, .xlsx: its columns are asked when not recognised; a position in one cell such as chr1:72523846-72667555 or an ISCN is read), BED files, SV / CNV VCF files (.vcf, .vcf.gz: DEL, DUP, CNV), or zip files of them (one patient per zip). Text in UTF-8, UTF-16 or Windows-1252; an old binary .xls must be saved as .xlsx, .csv or tab-separated text first';
/** a file's text read with the mapping remembered for its header when there is one, else by its header's names */
function readTable(text: string, name: string, fallback?: string): CghFile {
  const hit = findMap(splitRecords(text));
  if (hit) { try { return { ...parseCghFile(text, name, fallback, { ...hit.spec, opts: { ...hit.spec.opts, headerRow: hit.row } }), remembered: true }; } catch { /* the remembered columns do not fit: by names */ } }
  return parseCghFile(text, name, fallback);
}
type ReadOne = { file: CghFile; note?: string } | { toMap: CghState['toMap'][number] };
/** one file's bytes: gzip unpacked, an .xlsx's first worksheet, a VCF, else a table or a BED (its columns asked when not recognised) */
function readBytes(bytes: Uint8Array, name: string, fallback?: string): ReadOne {
  const b = isGzip(bytes) ? gunzip(bytes) : bytes;
  if (isOle(b)) throw new Error(OLE_MESSAGE);
  // a zip here is a workbook (zips themselves are opened by `add`)
  const xlsx = isZip(b);
  const text = xlsx ? xlsxText(b) : decodeText(b);
  if (!xlsx && isVcf(text)) return { file: parseVcf(text, name, fallback) };
  try {
    const file = readTable(text, name, fallback);
    return { file, note: xlsx ? `${leafOf(name)}: first worksheet read; exported as tab-separated text (${file.saveAs})` : undefined };
  } catch (x) { if (x instanceof NeedsMapping) return { toMap: { name, text, fallback } }; throw x; }
}
/** the tables of a zip, and of the zips inside it (one level): macOS's resource files left out, the other files counted */
function zipEntries(data: Uint8Array, depth = 0): { entries: { path: string; bytes: Uint8Array }[]; skipped: number } {
  let skipped = 0;
  const files = unzipSync(data, { filter: f => {
    if (f.name.endsWith('/') || /(^|\/)__MACOSX\/|(^|\/)\._/.test(f.name)) return false;
    const keep = (TABLE.test(f.name) || (depth === 0 && /\.zip$/i.test(f.name))) && f.originalSize < 200e6;
    if (!keep) skipped++;
    return keep;
  } });
  const entries: { path: string; bytes: Uint8Array }[] = [];
  for (const [path, bytes] of Object.entries(files)) {
    if (/\.zip$/i.test(path)) { const inner = zipEntries(bytes, depth + 1); skipped += inner.skipped; entries.push(...inner.entries.map(e => ({ path: `${path} › ${e.path}`, bytes: e.bytes }))); }
    else entries.push({ path, bytes });
  }
  return { entries, skipped };
}

/** a patient's classes or notes changed / exported: the patients still to export, and whether any is */
const touched = (key: string) => { const dirtyKeys = { ...state.dirtyKeys, [key]: true as const }; return { dirtyKeys, dirty: true }; };
const exportedOf = (keys: string[]) => { const dirtyKeys = { ...state.dirtyKeys }; for (const k of keys) delete dirtyKeys[k]; return { dirtyKeys, dirty: Object.keys(dirtyKeys).length > 0 }; };

export const cgh = {
  get: () => state,
  setOpen(open: boolean) { set({ open }); },
  setAnn(id: string, a: CnvAnn) { set({ ann: { ...state.ann, [id]: a } }); },
  setModes(m: Record<string, GeneModes | 'loading' | 'error'>) { set({ modes: { ...state.modes, ...m } }); },
  /** read the files picked; they wait for the minimum number of probes */
  async add(list: FileList | File[]) {
    const files: CghFile[] = [], errors: string[] = [], notes: string[] = [], toMap: CghState['toMap'] = [];
    const seen = new Set([...state.files, ...(state.pending?.files ?? [])].map(f => f.lines.join('\n')));
    const take = (f: CghFile) => { const k = f.lines.join('\n'); if (seen.has(k)) { notes.push(`${f.name}: already added`); return; } seen.add(k); files.push(f); };
    for (const f of Array.from(list)) {
      try {
        const buf = new Uint8Array(await f.arrayBuffer());
        if (/\.zip$/i.test(f.name)) {
          // one patient per zip: its interval tables (a .txt must be one: other text files are skipped), its BED files
          const { entries, skipped } = zipEntries(buf);
          const patient = f.name.replace(/\.zip$/i, '');
          let n = 0, other = skipped, asked = 0;
          for (const e of entries) {
            try {
              const r = readBytes(e.bytes, `${f.name} › ${e.path}`, patient);
              if ('toMap' in r) { toMap.push(r.toMap); asked++; continue; }
              if (r.file.format === 'bed' && !/\.bed$/i.test(e.path)) { other++; continue; }
              take(r.file); n++;
              if (r.note) notes.push(r.note);
            } catch (x) { if ((x as Error).message === OLE_MESSAGE) errors.push(`${f.name} › ${e.path}: ${OLE_MESSAGE}`); else other++; }
          }
          if (!n && !asked) errors.push(`${f.name}: no interval table or BED inside (${entries.length + skipped} file${entries.length + skipped > 1 ? 's' : ''})`);
          else if (other) notes.push(`${f.name}: ${n} table${n > 1 ? 's' : ''} read, ${other} other file${other > 1 ? 's' : ''} skipped`);
        } else {
          const r = readBytes(buf, f.name);
          if ('toMap' in r) toMap.push(r.toMap); else { take(r.file); if (r.note) notes.push(r.note); }
        }
      } catch (e) { errors.push(`${f.name}: ${(e as Error).message.replace(`${f.name}: `, '')}`); }
    }
    // files dropped while the question is open join it
    const p = state.pending, next = { files: [...(p?.files ?? []), ...files], errors: [...(p?.errors ?? []), ...errors], notes: [...(p?.notes ?? []), ...notes] };
    set({ pending: next.files.length || next.errors.length || next.notes.length ? next : null, toMap: [...state.toMap, ...toMap] });
  },
  /** the table at the head of the queue read with the columns chosen (null: skipped, or kept as it was) */
  mapped(file: CghFile | null) {
    const [head, ...toMap] = state.toMap;
    const add = file ?? head?.orig;
    const pending = state.pending ?? { files: [], errors: [], notes: [] };
    const next = add ? { ...pending, files: [...pending.files, add] } : pending;
    set({ toMap, pending: toMap.length || next.files.length || next.errors.length || next.notes.length ? next : null });
  },
  /** a file read, to be read again with other columns */
  remap(f: CghFile) {
    const p = state.pending; if (!p || !f.spec) return;
    set({ pending: { ...p, files: p.files.filter(x => x !== f) }, toMap: [...state.toMap, { name: f.name, text: f.lines.join(f.eol), fallback: f.fallback, spec: f.spec, orig: f }] });
  },
  /** the files accepted with this minimum and the sexes chosen (by sample key; the file's own when not chosen) */
  accept(minProbes: number, sexes: Record<string, CghSex | null> = {}, minKb = state.filter.minKb) {
    const p = state.pending; if (!p) return;
    const sex = { ...state.sex };
    for (const x of p.files.flatMap(f => f.samples)) { const v = x.key in sexes ? sexes[x.key] : x.sex; if (v) sex[x.key] = v; else delete sex[x.key]; }
    const filter = { ...state.filter, minProbes, minKb };
    try { localStorage.setItem(FILTER_KEY, JSON.stringify({ ...filter, query: '' })); } catch { /* storage blocked */ }
    const files = [...state.files, ...p.files];
    set({ files, pending: null, filter, sex, sample: p.files[0]?.samples[0]?.key ?? state.sample });
  },
  cancel() { set({ pending: null }); },
  setLab(lab: LabBenign | null) { set({ lab, labStored: storeLab(lab), filter: lab ? state.filter : { ...state.filter, lab: 'any' } }); },
  setFilter(patch: Partial<CghFilter>) {
    const filter = { ...state.filter, ...patch };
    try { localStorage.setItem(FILTER_KEY, JSON.stringify({ ...filter, query: '' })); } catch { /* storage blocked */ }
    set({ filter });
  },
  setSample(key: string | null) { set({ sample: key }); },
  setClass(cnv: CghCnv, v: string) { set({ classes: { ...state.classes, [cnv.id]: v }, ...touched(cnv.sampleKey) }); },
  setNote(cnv: CghCnv, v: string) { set({ notes: { ...state.notes, [cnv.id]: v }, ...touched(cnv.sampleKey) }); },
  setShown(id: string | null) { set({ shown: id }); },
  setSex(key: string, v: CghSex | null) { const sex = { ...state.sex }; if (v) sex[key] = v; else delete sex[key]; set({ sex }); },
  setDlrs(key: string, v: number | null) { const dlrs = { ...state.dlrs }; if (v != null && v > 0) dlrs[key] = v; else delete dlrs[key]; set({ dlrs }); },
  /** a sample taken off the list (its file too when it was the file's only sample) */
  remove(sampleKey: string) {
    const files = state.files.map(f => ({ ...f, samples: f.samples.filter(s => s.key !== sampleKey) })).filter(f => f.samples.length);
    const left = files.flatMap(f => f.samples);
    set({ files, sample: state.sample === sampleKey ? left[0]?.key ?? null : state.sample, ...exportedOf([sampleKey]) });
  },
  /** these patients' classes and notes written to a file */
  exported(keys: string[]) { set(exportedOf(keys)); },
};

/** a call's class: the one given in the browser, else the file's (in the laboratory's words when recognised) */
export const classOf = (s: CghState, c: CghCnv) => s.classes[c.id] ?? toClass(c.classification) ?? c.classification;
/** a call's note: the one written in the browser, else the file's */
export const noteOf = (s: CghState, c: CghCnv) => s.notes[c.id] ?? c.note;

/** whether a call is listed */
export function passes(s: CghState, c: CghCnv): boolean {
  const f = s.filter;
  if (!isAnalysed(c, f.minProbes, f.minKb)) return false;
  if (f.kind !== 'all' && c.kind !== f.kind) return false;
  if (f.hideSuppressed && c.suppressed) return false;
  const cls = classOf(s, c).trim();
  if (f.classed === 'unclassified' && cls) return false;
  if (f.classed === 'classified' && !cls) return false;
  if (f.series !== 'any') { const n = seriesOf(s.files).get(c.id)?.hits.length ?? 0; if (f.series === 'unique' ? n > 0 : n < SERIES_RECURRENT) return false; }
  if (f.lab !== 'any' && (labHits(s.files, s.lab).has(c.id) !== (f.lab === 'in'))) return false;
  const q = f.query.trim().toLowerCase();
  if (q) {
    // a chromosome (3, chrX) exactly; a band from its start (1p3, 16p11); a gene by part of its name
    const qc = q.replace(/^chr/, '');
    const hit = /^(\d{1,2}|x|y)$/.test(qc) ? c.chrom.toLowerCase() === `chr${qc}`
      : `${c.chrom.replace(/^chr/, '')}${c.cytoband || s.ann[c.id]?.band || ''}`.toLowerCase().startsWith(qc) || c.genes.some(g => g.toLowerCase().includes(q));
    if (!hit) return false;
  }
  return true;
}
