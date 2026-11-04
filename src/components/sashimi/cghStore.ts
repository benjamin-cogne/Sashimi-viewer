/**
 * The array-CGH samples of the SV browser, kept for the tab's life (the SV browser can be closed and opened again):
 * the files read, the class given to each call, the filter, the call shown. Nothing of a patient is written to the
 * browser's storage: only the settings (the minimum number of probes, the filter) are remembered.
 */
import { useSyncExternalStore } from 'react';
import { parseCghFile, type CghFile, type CghCnv } from './cghIntervals';

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
}
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
  pending: { files: CghFile[]; errors: string[] } | null;
  /** classes changed since the last export */
  dirty: boolean;
  /** the array-CGH page is open (over the SV browser) */
  open: boolean;
}

const FILTER_KEY = 'sashimi.cghFilter';
export const CGH_MIN_PROBES_DEFAULT = 3;
const loadFilter = (): CghFilter => {
  const d: CghFilter = { minProbes: CGH_MIN_PROBES_DEFAULT, kind: 'all', minKb: 0, classed: 'all', hideSuppressed: true, query: '' };
  try { const v = JSON.parse(localStorage.getItem(FILTER_KEY) ?? 'null') as Partial<CghFilter> | null; return v ? { ...d, ...v, query: '' } : d; } catch { return d; }
};

let state: CghState = { files: [], classes: {}, sample: null, filter: loadFilter(), shown: null, pending: null, dirty: false, open: false };
const subs = new Set<() => void>();
const set = (patch: Partial<CghState>) => { state = { ...state, ...patch }; subs.forEach(f => f()); };
const subscribe = (f: () => void) => { subs.add(f); return () => { subs.delete(f); }; };
export const useCgh = () => useSyncExternalStore(subscribe, () => state);

// a tab closed with classes typed and not exported asks first
const onUnload = (e: BeforeUnloadEvent) => { if (state.dirty) { e.preventDefault(); e.returnValue = ''; } };
if (typeof window !== 'undefined') window.addEventListener('beforeunload', onUnload);

/** a file's text: UTF-8, else Windows-1252 (an export from a Windows PC) */
async function readText(f: File): Promise<string> {
  const buf = await f.arrayBuffer();
  try { return new TextDecoder('utf-8', { fatal: true }).decode(buf); } catch { return new TextDecoder('windows-1252').decode(buf); }
}

export const cgh = {
  get: () => state,
  setOpen(open: boolean) { set({ open }); },
  /** read the files picked; they wait for the minimum number of probes */
  async add(list: FileList | File[]) {
    const files: CghFile[] = [], errors: string[] = [];
    for (const f of Array.from(list)) {
      try { files.push(parseCghFile(await readText(f), f.name)); } catch (e) { errors.push((e as Error).message); }
    }
    set({ pending: { files, errors } });
  },
  /** the files accepted with this minimum: listed, the first new sample shown */
  accept(minProbes: number) {
    const p = state.pending; if (!p) return;
    const filter = { ...state.filter, minProbes };
    try { localStorage.setItem(FILTER_KEY, JSON.stringify({ ...filter, query: '' })); } catch { /* storage blocked */ }
    const files = [...state.files, ...p.files];
    set({ files, pending: null, filter, sample: p.files[0]?.samples[0]?.key ?? state.sample });
  },
  cancel() { set({ pending: null }); },
  setFilter(patch: Partial<CghFilter>) {
    const filter = { ...state.filter, ...patch };
    try { localStorage.setItem(FILTER_KEY, JSON.stringify({ ...filter, query: '' })); } catch { /* storage blocked */ }
    set({ filter });
  },
  setSample(key: string | null) { set({ sample: key }); },
  setClass(cnv: CghCnv, v: string) { set({ classes: { ...state.classes, [cnv.id]: v }, dirty: true }); },
  setShown(id: string | null) { set({ shown: id }); },
  /** a sample taken off the list (its file too when it was the file's only sample) */
  remove(sampleKey: string) {
    const files = state.files.map(f => ({ ...f, samples: f.samples.filter(s => s.key !== sampleKey) })).filter(f => f.samples.length);
    const left = files.flatMap(f => f.samples);
    set({ files, sample: state.sample === sampleKey ? left[0]?.key ?? null : state.sample, dirty: left.length ? state.dirty : false });
  },
  exported() { set({ dirty: false }); },
};

/** a call's class: the one given in the browser, else the file's */
export const classOf = (s: CghState, c: CghCnv) => s.classes[c.id] ?? c.classification;

/** whether a call is listed */
export function passes(s: CghState, c: CghCnv): boolean {
  const f = s.filter;
  if (c.probes != null && c.probes < f.minProbes) return false;
  if (f.kind !== 'all' && c.kind !== f.kind) return false;
  if (c.sizeKb < f.minKb) return false;
  if (f.hideSuppressed && c.suppressed) return false;
  const cls = classOf(s, c).trim();
  if (f.classed === 'unclassified' && cls) return false;
  if (f.classed === 'classified' && !cls) return false;
  const q = f.query.trim().toLowerCase();
  if (q) {
    // a chromosome (3, chrX) exactly; a band from its start (1p3, 16p11); a gene by part of its name
    const qc = q.replace(/^chr/, '');
    const hit = /^(\d{1,2}|x|y)$/.test(qc) ? c.chrom.toLowerCase() === `chr${qc}`
      : `${c.chrom.replace(/^chr/, '')}${c.cytoband}`.toLowerCase().startsWith(qc) || c.genes.some(g => g.toLowerCase().includes(q));
    if (!hit) return false;
  }
  return true;
}
