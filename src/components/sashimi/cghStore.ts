/**
 * The array-CGH samples of the SV browser, kept for the tab's life (the SV browser can be closed and opened again):
 * the files read, the class given to each call, the filter, the call shown. Nothing of a patient is written to the
 * browser's storage: only the filter here (the lab's benign list and the column mappings in their own modules).
 */
import { useSyncExternalStore } from 'react';
import { unzipSync } from 'fflate';
import { parseCghFile, aberrationHeader, leafOf, toClass, type CghFile, type CghCnv, type CghSex } from './cghIntervals';
import { NeedsMapping, findMap, splitRecords, type MapSpec } from './cghMapping';
import { decodeText, gunzip, isGzip, isOle, isZip, xlsxText, OLE_MESSAGE } from './cghFormats';
import { isVcf, parseVcf } from './cghVcf';
import { parseProbes, probeSampleOf, type CghProbes, type CghProbeChrom } from './cghProbes';
import { seriesOf } from './cghSeries';
import { detect, SEG_DEFAULTS, type ProbeCalls, type ProbeEvent, type SegOpts } from './cghSegment';
import ProbeWorker from './cghProbeWorker?worker&inline';
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
  /** recurrence in the series: the number of OTHER patients carrying the call (the Series column) compared to seriesN */
  series: SeriesOp;
  seriesN: number;
  /** the laboratory's benign list: any call, those not in it, those in it */
  lab: 'any' | 'not' | 'in';
}
export type SeriesOp = 'any' | '=' | '≥' | '≤' | '>' | '<';
export const SERIES_OPS: SeriesOp[] = ['any', '=', '≥', '≤', '>', '<'];
const seriesOk = (op: SeriesOp, n: number, k: number) =>
  op === '=' ? n === k : op === '≥' ? n >= k : op === '≤' ? n <= k : op === '>' ? n > k : op === '<' ? n < k : true;
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
  /** a note on each call, by call id (else the file's) */
  notes: Record<string, string>;
  /** the laboratory's list of benign CNVs (a reference, kept in localStorage), and whether it could be stored */
  lab: LabBenign | null; labStored: boolean;
  /** the patients (sample keys) whose calls under the minimums are listed too (this tab only) */
  revealed: Record<string, true>;
  /** the quick guide is open (Esc closes it, not the page) */
  help: boolean;
  /** the settings of the probe calls (cghSegment.ts), kept in localStorage */
  probeOpts: SegOpts;
  /** events detected from each patient's probes, by sample key (computed in the background, for these settings and sex) */
  probeCalls: Record<string, ProbeCalls & { key: string }>;
  /** the patients whose probe calls not in the file are listed */
  probeOpen: Record<string, true>;
  /** how the probes are drawn (kept in localStorage): the whole genome's layout, the log-ratio range ±scale (shared by the genome, chromosome and probe panels) */
  view: CghView;
}
/** `iscnClasses`: the classes whose calls make the patient's ISCN */
export interface CghView { layout: 'horizontal' | 'vertical'; scale: number; iscnClasses: string[] }
/** the log-ratio ranges (±): presets ±0.5, ±1 (default), ±2; ⊖ goes further out */
export const LR_SCALES = [0.5, 1, 2, 3, 4, 6, 8];
export const LR_PRESETS = [0.5, 1, 2];
const VIEW_KEY = 'sashimi.cghView';
const loadView = (): CghView => {
  const d: CghView = { layout: 'horizontal', scale: 1, iscnClasses: ['Pathogenic', 'Likely pathogenic', 'Uncertain significance'] };
  try { const v = { ...d, ...(JSON.parse(localStorage.getItem(VIEW_KEY) ?? 'null') as Partial<CghView> | null) }; return { layout: v.layout === 'vertical' ? 'vertical' : 'horizontal', scale: LR_SCALES.includes(v.scale) ? v.scale : 1, iscnClasses: Array.isArray(v.iscnClasses) ? v.iscnClasses.filter(x => typeof x === 'string') : d.iscnClasses }; } catch { return d; }
};

const FILTER_KEY = 'sashimi.cghFilter';
export const CGH_MIN_PROBES_DEFAULT = 3;
/** at most this many calls are annotated by themselves (genes, PanelApp, ClinGen, gnomAD, DGV), the largest first; the others on demand, as many a click */
export const CGH_MAX_ANNOTATED = 200;
/** the largest first, ties in file order (Array.sort is stable) */
export const largestFirst = <T extends { sizeKb: number }>(cs: T[]) => cs.slice().sort((a, b) => b.sizeKb - a.sizeKb);
/** a call is analysed (listed and annotated) when it has enough probes (a BED's calls have no count) and is large enough */
export const isAnalysed = (c: CghCnv, minProbes: number, minKb: number) => (c.probes == null || c.probes >= minProbes) && c.sizeKb >= minKb;
const PROBE_OPTS_KEY = 'sashimi.cghProbeCalls';
const loadProbeOpts = (): SegOpts => { try { return { ...SEG_DEFAULTS, ...(JSON.parse(localStorage.getItem(PROBE_OPTS_KEY) ?? 'null') as Partial<SegOpts> | null) }; } catch { return SEG_DEFAULTS; } };
const loadFilter = (): CghFilter => {
  const d: CghFilter = { minProbes: CGH_MIN_PROBES_DEFAULT, kind: 'all', minKb: 0, classed: 'all', hideSuppressed: true, query: '', sort: 'position', series: 'any', seriesN: 1, lab: 'any' };
  try { const v = JSON.parse(localStorage.getItem(FILTER_KEY) ?? 'null') as Partial<CghFilter> | null; return v ? { ...d, ...v, query: '', series: 'any', lab: 'any' } : d; } catch { return d; }
};

let state: CghState = { files: [], classes: {}, sample: null, filter: loadFilter(), shown: null, pending: null, toMap: [], dirty: false, dirtyKeys: {}, open: false, ann: {}, modes: {}, sex: {}, notes: {}, lab: loadStoredLab(), labStored: true, revealed: {}, help: false, probeOpts: loadProbeOpts(), probeCalls: {}, probeOpen: {}, view: loadView() };
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
export const CGH_FILES_TEXT = 'the interval table of each patient exported by CytoGenomics (.txt, or its automated IntervalBasedReport .xls), any other table of calls (.txt, .tsv, .csv, .xlsx: its columns are asked when not recognised; a position in one cell such as chr1:72523846-72667555 or an ISCN is read), BED files, SV / CNV VCF files (.vcf, .vcf.gz: DEL, DUP, CNV), or zip files of them (one patient per zip). Text in UTF-8, UTF-16 or Windows-1252; an old binary .xls must be saved as .xlsx, .csv or tab-separated text first';
/** a file's text read with the mapping remembered for its header when there is one, else by its header's names */
function readTable(text: string, name: string, fallback?: string): CghFile {
  const recs = splitRecords(text);
  // CytoGenomics' interval-based report is read by its own layout, never by a remembered mapping
  if (aberrationHeader(recs) >= 0) return parseCghFile(text, name, fallback);
  const hit = findMap(recs);
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
function zipEntries(data: Uint8Array, depth = 0): { entries: { path: string; bytes: Uint8Array }[]; skipped: string[] } {
  const skipped: string[] = [];
  const files = unzipSync(data, { filter: f => {
    if (f.name.endsWith('/') || /(^|\/)__MACOSX\/|(^|\/)\._/.test(f.name)) return false;
    const keep = (TABLE.test(f.name) || (depth === 0 && /\.zip$/i.test(f.name))) && f.originalSize < 200e6;
    if (!keep) skipped.push(f.name);
    return keep;
  } });
  const entries: { path: string; bytes: Uint8Array }[] = [];
  for (const [path, bytes] of Object.entries(files)) {
    if (/\.zip$/i.test(path)) { const inner = zipEntries(bytes, depth + 1); skipped.push(...inner.skipped); entries.push(...inner.entries.map(e => ({ path: `${path} › ${e.path}`, bytes: e.bytes }))); }
    else entries.push({ path, bytes });
  }
  return { entries, skipped };
}

/** a patient's classes or notes changed / exported: the patients still to export, and whether any is */
const touched = (key: string) => { const dirtyKeys = { ...state.dirtyKeys, [key]: true as const }; return { dirtyKeys, dirty: true }; };
/**
 * The probes' log ratios of each patient (CGH_Log_Ratio.txt), by sample key: the file's bytes until first asked, then
 * read once into typed arrays. Kept out of the state (no redraw), in this tab only: drawn by the probe panel, their DLRS in the import question.
 */
const probeBytes = new Map<string, Uint8Array>(), probeRead = new Map<string, CghProbes | null>();
/** each patient's probe data version: bumped when its probes are replaced or dropped (memos over the probes, results computed on older ones) */
const probeGen = new Map<string, number>();
const bump = (k: string) => probeGen.set(k, (probeGen.get(k) ?? 0) + 1);
const dropProbes = (keys: string[]) => { for (const k of keys) { probeBytes.delete(k); probeRead.delete(k); bump(k); } callsQueue = callsQueue.filter(k => !keys.includes(k)); };
/** what a patient's probe calls depend on: the settings and the sex (chrY read for a male only) */
const callsKey = (k: string) => `${state.sex[k] ?? ''}|${state.probeOpts.minProbes}|${state.probeOpts.minLr}|${state.probeOpts.t}`;
let callsQueue: string[] = [], callsBusy = false;
/** the worker computing the probe calls (cghSegWorker.ts); null: none could be made, they are computed on the page */
let segWorker: Worker | null | undefined, segReq = 0;
const segWaiting = new Map<number, (r: any) => void>();
function worker(): Worker | null {
  if (segWorker !== undefined) return segWorker;
  try {
    segWorker = new ProbeWorker();
    segWorker.onmessage = (e: MessageEvent<{ id: number; result?: unknown; error?: string }>) => { const f = segWaiting.get(e.data.id); segWaiting.delete(e.data.id); f?.(e.data.result ?? null); };
    segWorker.onerror = () => { segWorker?.terminate(); segWorker = null; for (const f of segWaiting.values()) f(null); segWaiting.clear(); };
  } catch { segWorker = null; }
  return segWorker;
}
/** a patient's probe calls: in the worker (a 1M-probe array takes seconds), else on the page */
function detectOff(pr: CghProbes, o: SegOpts, male: boolean): Promise<ProbeCalls> {
  const w = worker();
  if (!w) return Promise.resolve(detect(pr, o, male));
  const id = ++segReq, byChrom = new Map([...pr.byChrom].map(([c, p]) => [c, { start: p.start, end: p.end, lr: p.lr, name: [] as string[] }]));
  return new Promise(res => { segWaiting.set(id, r => res(r ?? detect(pr, o, male))); w.postMessage({ id, byChrom, build: pr.build, opts: o, male }); });
}
/** probe files being read in the worker, by sample key */
const probeParsing = new Map<string, Promise<CghProbes | null>>();
const parseSync = (b: Uint8Array) => { try { return parseProbes(b); } catch { return null; } };
/** a patient's probes read (in the worker when there is one: 1M probes take over a second); the page redrawn when they are */
function probesAsync(k: string): Promise<CghProbes | null> {
  if (probeRead.has(k)) return Promise.resolve(probeRead.get(k) ?? null);
  const b = probeBytes.get(k); if (!b) return Promise.resolve(null);
  const had = probeParsing.get(k); if (had) return had;
  const w = worker(), gen = probeGen.get(k);
  const pr: Promise<CghProbes | null> = !w ? Promise.resolve(parseSync(b)) : new Promise(res => {
    const id = ++segReq;
    segWaiting.set(id, (r: (Omit<CghProbes, 'byChrom'> & { byChrom: Map<string, Omit<CghProbeChrom, 'name'> & { name: string }> }) | null) =>
      res(r ? { ...r, byChrom: new Map([...r.byChrom].map(([c, x]) => [c, { ...x, name: x.name.split('\n') }])) } : parseSync(b)));
    w.postMessage({ id, parse: b });
  });
  const done = pr.then(p => {
    if (probeParsing.get(k) === done) probeParsing.delete(k);
    // the probes replaced or dropped meanwhile: this reading is dropped
    if (probeGen.get(k) !== gen || probeBytes.get(k) !== b) return probeRead.get(k) ?? null;
    probeRead.set(k, p); probeBytes.delete(k); bump(k); set({});
    return p;
  });
  probeParsing.set(k, done);
  return done;
}
/** one patient at a time, in the worker, so the page stays usable while they are computed */
async function nextCalls() {
  if (callsBusy) return;
  const k = callsQueue.shift(); if (k == null) return;
  callsBusy = true;
  try {
    const pr = await probesAsync(k), key = callsKey(k), gen = probeGen.get(k);
    if (pr && state.probeCalls[k]?.key !== key) {
      const r = await detectOff(pr, state.probeOpts, state.sex[k] === 'male');
      // kept only when the settings, sex and probes are still those it was computed for (else asked again)
      if (callsKey(k) === key && probeGen.get(k) === gen && cgh.hasProbes(k)) set({ probeCalls: { ...state.probeCalls, [k]: { ...r, key } } });
      else if (cgh.hasProbes(k) && !callsQueue.includes(k)) callsQueue.push(k);
    }
  } finally { callsBusy = false; }
  if (callsQueue.length) setTimeout(nextCalls, 0);
}
/** a probe-only event as a call, to be drawn like one (never in a file, never exported) */
export const probeCnv = (e: ProbeEvent, sampleKey: string, sample: string): CghCnv => ({
  id: `probe|${sampleKey}|${e.chrom}:${e.start}-${e.end}`, sample, sampleKey, chrom: e.chrom, start: e.start, end: e.end, cytoband: '', genes: [],
  sizeKb: (e.end - e.start + 1) / 1000, kind: e.kind, typeText: 'from probes', probes: e.n, logRatio: e.mean, iscn: null, copies: null, copiesHi: null,
  amp: false, suppressed: false, classification: '', classRest: '', note: '', line: -1, fromProbes: true });
/** a few names of the files left aside in a zip */
const someNames = (names: string[]) => names.map(n => leafOf(n)).slice(0, 4).join(', ') + (names.length > 4 ? ', …' : '');
const exportedOf = (keys: string[]) => { const dirtyKeys = { ...state.dirtyKeys }; for (const k of keys) delete dirtyKeys[k]; return { dirtyKeys, dirty: Object.keys(dirtyKeys).length > 0 }; };

export const cgh = {
  get: () => state,
  setOpen(open: boolean) { set({ open }); },
  setHelp(help: boolean) { set({ help }); },
  setAnn(id: string, a: CnvAnn) { set({ ann: { ...state.ann, [id]: a } }); },
  /** many calls' annotations at once (the queue's placeholders: one redraw) */
  setAnns(m: Record<string, CnvAnn>) { set({ ann: { ...state.ann, ...m } }); },
  /** a patient's calls under the minimums listed, or hidden again */
  setRevealed(key: string, on: boolean) { const revealed = { ...state.revealed }; if (on) revealed[key] = true; else delete revealed[key]; set({ revealed }); },
  setModes(m: Record<string, GeneModes | 'loading' | 'error'>) { set({ modes: { ...state.modes, ...m } }); },
  /** whether a patient has its probes' log ratios (CGH_Log_Ratio.txt) */
  hasProbes: (sampleKey: string) => probeBytes.has(sampleKey) || !!probeRead.get(sampleKey),
  /** a patient's probes are being read (in the worker) */
  probesReading: (sampleKey: string) => probeParsing.has(sampleKey),
  /** a patient's probe data version (a memo over its probes depends on it) */
  probeGen: (sampleKey: string) => probeGen.get(sampleKey) ?? 0,
  /** a patient's probes' log ratios, read on first use in the worker (null: none, not readable, or being read: the page is redrawn when they are) */
  probes(sampleKey: string): CghProbes | null {
    if (!probeRead.has(sampleKey)) { void probesAsync(sampleKey); return probeRead.get(sampleKey) ?? null; }
    return probeRead.get(sampleKey) ?? null;
  },
  /** the probe calls of these patients, computed in the background when missing or out of date */
  runProbeCalls(keys: string[]) {
    const todo = keys.filter(k => cgh.hasProbes(k) && state.probeCalls[k]?.key !== callsKey(k) && !callsQueue.includes(k));
    if (!todo.length) return;
    callsQueue.push(...todo); if (!callsBusy) setTimeout(nextCalls, 0);
  },
  /** the current probe calls of a patient (null: none yet, or out of date) */
  probeCallsOf: (k: string): ProbeCalls | null => (state.probeCalls[k]?.key === callsKey(k) ? state.probeCalls[k] : null),
  setProbeOpts(patch: Partial<SegOpts>) {
    const probeOpts = { ...state.probeOpts, ...patch };
    try { localStorage.setItem(PROBE_OPTS_KEY, JSON.stringify(probeOpts)); } catch { /* storage blocked */ }
    set({ probeOpts });
  },
  setView(patch: Partial<CghView>) {
    const view = { ...state.view, ...patch };
    try { localStorage.setItem(VIEW_KEY, JSON.stringify(view)); } catch { /* storage blocked */ }
    set({ view });
  },
  setProbeOpen(key: string, on: boolean) { const probeOpen = { ...state.probeOpen }; if (on) probeOpen[key] = true; else delete probeOpen[key]; set({ probeOpen }); },
  /** read the files picked; they wait for the minimum number of probes */
  async add(list: FileList | File[]) {
    const files: CghFile[] = [], errors: string[] = [], notes: string[] = [], toMap: CghState['toMap'] = [];
    const seen = new Map([...state.files, ...(state.pending?.files ?? [])].map(f => [f.lines.join('\n'), f] as const));
    /** the file accepted, or the one already loaded with the same content (its samples get the zip's probes) */
    const take = (f: CghFile): CghFile => { const k = f.lines.join('\n'), had = seen.get(k); if (had) { notes.push(`${f.name}: already added`); return had; } seen.set(k, f); files.push(f); return f; };
    // probe log-ratio files given alone (or in a zip without a table): matched to a sample by name after the loop
    const loose: { name: string; sample: string; bytes: Uint8Array }[] = [];
    // new probes for a patient: its probe calls of the former ones dropped, out of the queue
    const stale: string[] = [];
    const keep = (keys: string[], bytes: Uint8Array) => { dropProbes(keys); for (const k of keys) { probeBytes.set(k, bytes); stale.push(k); } };
    for (const f of Array.from(list)) {
      try {
        const buf = new Uint8Array(await f.arrayBuffer());
        if (/\.zip$/i.test(f.name)) {
          // one patient per zip: its interval tables (a .txt must be one: other text files are skipped), its BED files
          const { entries, skipped } = zipEntries(buf);
          const patient = f.name.replace(/\.zip$/i, '');
          // the tables are found by their content; the other files (images, XML, QC tables…) are left aside, counted
          const other = [...skipped], asks: CghState['toMap'] = [], read: CghFile[] = [], probesIn: { name: string; sample: string; bytes: Uint8Array }[] = [];
          for (const e of entries) {
            const path = `${f.name} › ${e.path}`;
            try {
              const raw = isGzip(e.bytes) ? gunzip(e.bytes) : e.bytes;
              const ps = probeSampleOf(raw);
              if (ps != null) { probesIn.push({ name: path, sample: ps, bytes: raw }); continue; }
              const r = readBytes(raw, path, patient);
              if ('toMap' in r) { asks.push(r.toMap); continue; }
              if (r.file.format === 'bed' && !/\.bed$/i.test(e.path)) { other.push(e.path); continue; }
              read.push(take(r.file));
              if (r.note) notes.push(r.note);
            } catch (x) { if ((x as Error).message === OLE_MESSAGE) errors.push(`${path}: ${OLE_MESSAGE}`); else other.push(e.path); }
          }
          const n = read.length;
          // a table whose columns are not recognised is asked only when the zip holds no table read by itself
          if (n) other.push(...asks.map(a => a.name)); else toMap.push(...asks);
          // one patient per zip: its probe log ratios go with its sample(s) (by name when several)
          const ss = read.flatMap(x => x.samples);
          for (const p of probesIn) {
            const to = ss.length === 1 ? ss : ss.filter(x => x.name === p.sample);
            if (to.length) { keep(to.map(x => x.key), p.bytes); notes.push(`${f.name}: probe log ratios (${leafOf(p.name)}) kept for ${to[0].name}`); }
            else loose.push(p);
          }
          const total = entries.length + skipped.length;
          if (!n && !asks.length && !probesIn.length) errors.push(`${f.name}: no interval table or BED inside (${total} file${total > 1 ? 's' : ''})`);
          else if (other.length) notes.push(`${f.name}: ${n} table${n !== 1 ? 's' : ''} read; ignored: ${other.length} other file${other.length > 1 ? 's' : ''} (${someNames(other)})`);
        } else {
          const raw = isGzip(buf) ? gunzip(buf) : buf, ps = probeSampleOf(raw);
          if (ps != null) { loose.push({ name: f.name, sample: ps, bytes: raw }); continue; }
          const r = readBytes(buf, f.name);
          if ('toMap' in r) toMap.push(r.toMap); else { take(r.file); if (r.note) notes.push(r.note); }
        }
      } catch (e) { errors.push(`${f.name}: ${(e as Error).message.replace(`${f.name}: `, '')}`); }
    }
    // a probe log-ratio file alone: the sample(s) of that name, loaded, waiting or read now
    const all = [...state.files, ...(state.pending?.files ?? []), ...files].flatMap(x => x.samples);
    for (const p of loose) {
      const to = all.filter(x => x.name === p.sample);
      if (to.length) { keep(to.map(x => x.key), p.bytes); notes.push(`${leafOf(p.name)}: probe log ratios kept for ${p.sample}`); }
      else notes.push(`${leafOf(p.name)}: probe log ratios of "${p.sample}" not matched to a loaded sample (not kept)`);
    }
    // files dropped while the question is open join it
    const p = state.pending, next = { files: [...(p?.files ?? []), ...files], errors: [...(p?.errors ?? []), ...errors], notes: [...(p?.notes ?? []), ...notes] };
    const probeCalls = { ...state.probeCalls }; for (const k of stale) delete probeCalls[k];
    set({ pending: next.files.length || next.errors.length || next.notes.length ? next : null, toMap: [...state.toMap, ...toMap], probeCalls });
    // probes given to patients already listed: their probe calls now (the list's effect runs only when the files change)
    cgh.runProbeCalls(stale.filter(k => state.files.some(f => f.samples.some(x => x.key === k))));
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
  cancel() { dropProbes(state.pending?.files.flatMap(f => f.samples.map(x => x.key)) ?? []); set({ pending: null }); },
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
  /** a sample taken off the list (its file too when it was the file's only sample) */
  remove(sampleKey: string) {
    const files = state.files.map(f => ({ ...f, samples: f.samples.filter(s => s.key !== sampleKey) })).filter(f => f.samples.length);
    const left = files.flatMap(f => f.samples);
    dropProbes([sampleKey]);
    const probeCalls = { ...state.probeCalls }; delete probeCalls[sampleKey];
    set({ files, probeCalls, sample: state.sample === sampleKey ? left[0]?.key ?? null : state.sample, ...exportedOf([sampleKey]) });
  },
  /** these patients' classes and notes written to a file */
  exported(keys: string[]) { set(exportedOf(keys)); },
};

/** a call's class: the one given in the browser, else the file's (in the laboratory's words when recognised) */
export const classOf = (s: CghState, c: CghCnv) => s.classes[c.id] ?? toClass(c.classification) ?? c.classification;
/** the build of a patient's file (GRCh37 when it states none) */
export const buildOfKey = (s: CghState, key: string) => s.files.find(f => f.samples.some(x => x.key === key))?.build ?? 'GRCh37';
/** a call's note: the one written in the browser, else the file's */
export const noteOf = (s: CghState, c: CghCnv) => s.notes[c.id] ?? c.note;

/** whether a call is listed; `belowToo`: the minimums left aside (a patient's calls under them, revealed on demand) */
export function passes(s: CghState, c: CghCnv, belowToo = false): boolean {
  const f = s.filter;
  if (!belowToo && !isAnalysed(c, f.minProbes, f.minKb)) return false;
  if (f.kind !== 'all' && c.kind !== f.kind) return false;
  // a chrX row typed by copy number in a file without the sex: hidden when it is the normal one for the sex set (x1 male, x2 female)
  if (c.cnSex) { const sx = s.sex[c.sampleKey]; if (sx && c.copies === (sx === 'male' ? 1 : 2)) return false; }
  if (f.hideSuppressed && c.suppressed) return false;
  const cls = classOf(s, c).trim();
  if (f.classed === 'unclassified' && cls) return false;
  if (f.classed === 'classified' && !cls) return false;
  if (f.series !== 'any') { const n = seriesOf(s.files).get(c.id)?.hits.length ?? 0; if (!seriesOk(f.series, n, f.seriesN)) return false; }
  if (f.lab !== 'any' && (labHits(s.files, s.lab).has(c.id) !== (f.lab === 'in'))) return false;
  const q = f.query.trim().toLowerCase();
  if (q) {
    // a chromosome (3, chrX) exactly; a band from its start (1p3, 16p11); a gene by part of its name
    const qc = q.replace(/^chr/, '');
    const hit = /^(\d{1,2}|x|y)$/.test(qc) ? c.chrom.toLowerCase() === `chr${qc}`
      : `${c.chrom.replace(/^chr/, '')}${c.cytoband || s.ann[c.id]?.band || ''}`.toLowerCase().startsWith(qc)
        // the file's Gene Name and the RefSeq genes read for the call
        || c.genes.some(g => g.toLowerCase().includes(q)) || !!s.ann[c.id]?.genes.some(g => g.symbol.toLowerCase().includes(q));
    if (!hit) return false;
  }
  return true;
}
