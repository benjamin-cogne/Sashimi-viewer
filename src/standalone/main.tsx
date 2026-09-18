/**
 * Standalone Sashimi viewer: a single HTML page. The user adds BAM/CRAM files (with
 * their .bai/.crai), optionally a reference FASTA, and types a gene. Everything is
 * decoded in the browser; only gene lookups (RefSeq models from the UCSC API, Ensembl REST as
 * fallback) and reference sequence (when no FASTA is given) are fetched from the network.
 */
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import ReactDOM from 'react-dom/client';
import SashimiViewer, { DEFAULT_VIEWER_SETTINGS, type ViewerSettings, type ViewerState } from '../components/SashimiViewer';
import { buildSession, defaultSessionName, matchSession, parseSession, viewerSettingsOf, type SessionFile } from './session';
import { fileInFolder, filesFromDrop, filesFromFolderInput, filesInFolder, hasFileSystemAccess, permitted, pickFiles, pickFolder, recallFile, recallFolder, rememberFiles, rememberFolder, type FSDirHandle, type FSHandle, type PathedFile } from './handles';
import { LocalDataSource, type LocalSample } from './localSource';
import { EMBEDDED_APP, EMBEDDED_VERSION, EmbeddedDataSource, buildExportHtml, embeddedSamples, encodeCoverageV2, encodeReadsV2, pageIsUnbuilt, readEmbedded, type EmbeddedExport, type EmbeddedView, type EncodedCoverage, type EncodedCoverageV2, type EncodedReadsV2 } from './embedded';
import type { GenomeBuild } from './ensembl';
import { parseLocus, toTxModel } from '../components/sashimi/geometry';
import type { KnownVariant, LibraryEvidence, LibraryType } from '../components/sashimi/types';
import { safeFileName, serializePlotSvg, stackSvgs } from '../components/sashimi/svgExport';
import { describeLink, parseLink, variantOfInterest } from './link';
import '../index.css';

/** Unreleased build (branch dev published under /dev/): banner, tab title and a red favicon, so it is never mistaken for the stable page. */
const DEV = import.meta.env.VITE_DEV_MODE === '1'
  ? { branch: import.meta.env.VITE_DEV_BRANCH || 'dev', sha: (import.meta.env.VITE_DEV_SHA || '').slice(0, 7), date: import.meta.env.VITE_DEV_DATE || '', stable: import.meta.env.VITE_STABLE_URL || '../' }
  : null;
if (DEV) {
  document.title = `[DEV] ${document.title}`;
  const icon = document.querySelector<HTMLLinkElement>('link[rel="icon"]');
  if (icon) icon.href = 'data:image/svg+xml,' + encodeURIComponent('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64"><rect width="64" height="64" rx="14" fill="#dc2626"/><path d="M11 44 Q32 -22 53 44" fill="none" stroke="#fff" stroke-width="6" stroke-linecap="round"/><rect x="2" y="44" width="16" height="14" rx="2" fill="#fff"/><rect x="24" y="44" width="16" height="14" rx="2" fill="#fff"/><rect x="46" y="44" width="16" height="14" rx="2" fill="#fff"/></svg>');
}

/** Request carried by the URL (deep link from another tool), read once at start-up. */
const LINK = parseLink(window.location.hash, window.location.search);
const LINK_TEXT = LINK ? `${LINK.mark.chrom}:${LINK.mark.start.toLocaleString('en-US')}${LINK.mark.end > LINK.mark.start ? `-${LINK.mark.end.toLocaleString('en-US')}` : ''}` : '';

/** The data an exported page carries (null for the normal viewer). */
const EMBEDDED = readEmbedded();

/** A region open in the viewer: the gene (1-based inclusive bounds) and, optionally, the window, the pinned locus and the reads track. */
type Opened = { geneName: string; geneId?: string; chrom: string; start: number; end: number; view?: { start: number; end: number }; /** undefined: the window itself is pinned when there is one; null: nothing pinned */ mark?: { start: number; end: number } | null; reads?: boolean };
/** A registered view (tab): its region, the last state the viewer reported for it, and the options it must start with when reopened. */
interface ViewTab { id: number; label: string; opened: Opened; state: ViewerState | null; settings?: Partial<ViewerSettings> }
const fmtLocus = (chrom: string, start: number, end: number) => `${chrom}:${start.toLocaleString('en-US')}${end > start ? `-${end.toLocaleString('en-US')}` : ''}`;
/** Tab label: the gene, with the searched locus when there is one (the window itself is in the tooltip). */
const labelOf = (o: Opened) => `${o.geneName}${o.mark ? ` · ${fmtLocus(o.chrom, o.mark.start, o.mark.end)}` : ''}`;
/** A tab's state: what its viewer last reported, else (a tab loaded from a session or an export and never opened yet) its region with the options it will start with. */
const stateOfTab = (t: ViewTab): ViewerState => t.state ?? {
  ...DEFAULT_VIEWER_SETTINGS, ...(t.settings ?? {}),
  gene: { name: t.opened.geneName, id: t.opened.geneId, chrom: t.opened.chrom, start: t.opened.start, end: t.opened.end },
  view: { chrom: t.opened.chrom, start: t.opened.view?.start ?? t.opened.start, end: t.opened.view?.end ?? t.opened.end },
  mark: t.opened.mark ?? null,
};
const openedOfState = (st: ViewerState, prev?: Opened): Opened => ({ geneName: st.gene.name, geneId: st.gene.id ?? prev?.geneId, chrom: st.gene.chrom, start: st.gene.start, end: st.gene.end, view: { ...st.view }, mark: st.mark ?? null });
/** Largest window fetched around a view (the viewer's own rule), for the export. */
const MAX_FETCH_BP = 2_000_000;
/** The reads track exists below this window size (the viewer's rule); the export follows it. */
const READS_MAX_VIEW_BP = 100_000;
/** Choices of the export dialog: the window exported around each view (coverage, junctions, retention and reads alike) and the reads per sample. */
interface ExportOptions { window: 'view' | 'margin' | 'max'; readsCap: 'shown' | 'dense' | 'all' }
const READS_CAPS: Record<ExportOptions['readsCap'], number> = { shown: 20000, dense: 100000, all: Number.MAX_SAFE_INTEGER };

const ALIGN_EXT = /\.(bam|cram)$/i;
const INDEX_EXT = /\.(bai|crai)$/i;
const FASTA_EXT = /\.(fa|fasta|fna)(\.gz)?$/i;

/** Pair alignment files with their indexes by name (case-insensitive, Windows tools often shout): x.bam + x.bam.bai or x.bai. */
type PairedSample = { name: string; kind: 'bam' | 'cram'; file: File; index: File; path?: string; indexPath?: string };
/** Pairs alignments with their index (same directory when the files carry relative paths) and finds the reference FASTA. */
function pairFiles(input: (File | PathedFile)[]): { samples: PairedSample[]; unmatched: string[]; fasta?: { fa: File; fai: File; gzi?: File }; fastaMissing?: string } {
  const items: PathedFile[] = input.map(x => (x instanceof File ? { file: x } : x));
  const files = items.map(x => x.file);
  const dirOf = (x: PathedFile) => (x.path ? x.path.slice(0, x.path.lastIndexOf('/') + 1) : '');
  const byKey = new Map(items.map(x => [`${dirOf(x)}|${x.file.name.toLowerCase()}`, x]));
  const lookupIn = (dir: string, n: string) => byKey.get(`${dir}|${n.toLowerCase()}`);
  const lookup = (n: string) => lookupIn('', n)?.file;
  const samples: PairedSample[] = [];
  const unmatched: string[] = [];
  for (const x of items) {
    const f = x.file;
    if (!ALIGN_EXT.test(f.name)) continue;
    const kind = f.name.toLowerCase().endsWith('.cram') ? 'cram' : 'bam';
    const stem = f.name.replace(ALIGN_EXT, '');
    const dir = dirOf(x);
    const idx = lookupIn(dir, `${f.name}.${kind === 'bam' ? 'bai' : 'crai'}`) || lookupIn(dir, `${stem}.${kind === 'bam' ? 'bai' : 'crai'}`);
    if (idx) samples.push({ name: stem, kind, file: f, index: idx.file, path: x.path, indexPath: idx.path });
    else unmatched.push(f.name);
  }
  let fasta: { fa: File; fai: File; gzi?: File } | undefined, fastaMissing: string | undefined;
  const fa = files.find(f => FASTA_EXT.test(f.name));
  if (fa) {
    const fai = lookup(`${fa.name}.fai`);
    const gzi = lookup(`${fa.name}.gzi`);
    if (fai && (!fa.name.toLowerCase().endsWith('.gz') || gzi)) fasta = { fa, fai, gzi };
    else fastaMissing = fa.name;
  }
  return { samples, unmatched, fasta, fastaMissing };
}

/** The project mark (docs/logo): three exons, canonical arcs, and a dashed red arc skipping the middle exon. */
function Logo({ size = 28 }: { size?: number }) {
  return (
    <svg width={size} height={size} viewBox="0 0 64 64" aria-hidden="true" className="shrink-0">
      <path d="M11 46 Q21.5 14 32 46 M32 46 Q42.5 14 53 46" fill="none" stroke="#3730a3" strokeWidth="3.5" strokeLinecap="round" />
      <path d="M11 46 Q32 -26 53 46" fill="none" stroke="#dc2626" strokeWidth="3.5" strokeLinecap="round" strokeDasharray="5 4" />
      <path d="M18 51 H25 M39 51 H46" stroke="#3730a3" strokeWidth="2.5" strokeLinecap="round" />
      <rect x="4" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
      <rect x="25" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
      <rect x="46" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
    </svg>
  );
}

function App() {
  const [build, setBuild] = useState<GenomeBuild>(EMBEDDED?.build ?? LINK?.build ?? 'GRCh38');
  const [samples, setSamples] = useState<LocalSample[]>(() => (EMBEDDED ? embeddedSamples(EMBEDDED) : []));
  const [renaming, setRenaming] = useState<{ id: number; value: string } | null>(null);
  const [fasta, setFasta] = useState<{ fa: File; fai: File; gzi?: File } | undefined>();
  const [notes, setNotes] = useState<string[]>([]);
  /** What the page is doing with a folder or files that were just given, before their chips exist (listing a folder, reading a drop, reopening a session's folder). */
  const [intake, setIntake] = useState<string | null>(null);
  const [gene, setGene] = useState(LINK_TEXT);
  // ---- Views: every region opened from the header is a tab; the active one drives the viewer ----
  const [views, setViews] = useState<ViewTab[]>([]);
  const [activeId, setActiveId] = useState<number | null>(null);
  const activeIdRef = useRef<number | null>(null);
  activeIdRef.current = activeId;
  const tabSeq = useRef(1);
  const active = views.find(v => v.id === activeId) ?? null;
  const opened: Opened | null = active?.opened ?? null;
  /** Options the next mounted viewer starts with (a reopened tab's snapshot, a session); cleared once that viewer reports. */
  const pendingSettingsRef = useRef<Partial<ViewerSettings> | undefined>();
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // ---- Variants of interest (header): drawn on every view through the viewer's known-variant layer ----
  const [knownVars, setKnownVars] = useState<KnownVariant[]>(() => EMBEDDED?.knownVariants ?? LINK?.variants ?? []);
  const knownSeq = useRef(0);
  const nextKnownId = useRef(1);
  const [knownLocus, setKnownLocus] = useState('');
  const [knownLabel, setKnownLabel] = useState('');
  const [knownError, setKnownError] = useState<string | null>(null);
  const [exportDialog, setExportDialog] = useState<ExportOptions | null>(null);
  /** A running export (HTML or SVG): what it is doing and how far it is, shown in a blocking overlay. */
  const [exportProgress, setExportProgress] = useState<{ title: string; step: string; done: number; total: number } | null>(null);
  // ---- Sessions: the viewer's latest state (options + navigation), a loaded session waiting for its files ----
  const viewerStateRef = useRef<ViewerState | null>(null);
  const [sessionName, setSessionName] = useState(() => defaultSessionName());
  const [sessionNameEdited, setSessionNameEdited] = useState(false);
  const [pendingSession, setPendingSession] = useState<{ file: SessionFile; name: string } | null>(null);
  const [sessionSeq, setSessionSeq] = useState(0);   // bumps the viewer key so a loaded session always remounts the viewer
  /** The run folder the sample paths are relative to; the handle (Chromium) lets a session reopen it later. */
  const [runFolder, setRunFolder] = useState<{ name: string; handle: FSDirHandle | null } | null>(null);
  const fileHandlesRef = useRef(new Map<string, FSHandle>());   // Chromium bookmarks of individually added files, by name
  const [reopenState, setReopenState] = useState<{ folder: string; ready: boolean; note?: string } | null>(null);
  const nextId = useRef(EMBEDDED ? Math.max(0, ...EMBEDDED.samples.map(s => s.id)) + 1 : 1);
  const dsRef = useRef<LocalDataSource>();
  if (!dsRef.current) { dsRef.current = EMBEDDED ? new EmbeddedDataSource(EMBEDDED, { build }) : new LocalDataSource({ build }); if (LINK && !EMBEDDED) dsRef.current.knownVariants = LINK.variants; }
  const ds = dsRef.current;

  /** The tabs with the active one's state refreshed from the viewer (label and region follow where the user went). */
  const snapshot = useCallback((tabs: ViewTab[]): ViewTab[] => {
    const st = viewerStateRef.current, id = activeIdRef.current;
    if (!st || id == null) return tabs;
    return tabs.map(v => (v.id === id && st.gene.chrom ? { ...v, state: st, opened: openedOfState(st, v.opened), label: labelOf(openedOfState(st, v.opened)) } : v));
  }, []);
  /** Opens a region in a new tab; the current view folds into its tab. The new viewer keeps the current options (a chosen transcript only for the same gene). */
  const openNew = useCallback((o: Opened) => {
    const id = tabSeq.current++;
    const prev = viewerStateRef.current;
    pendingSettingsRef.current = prev ? { ...prev, transcriptId: prev.gene.name === o.geneName ? prev.transcriptId : undefined } : undefined;
    setViews(tabs => [...snapshot(tabs), { id, label: labelOf(o), opened: o, state: null }]);
    setActiveId(id);
    setGene(o.mark ? fmtLocus(o.chrom, o.mark.start, o.mark.end) : o.geneName);
  }, [snapshot]);
  /** Reopens a tab with everything as it was left (full snapshot of the options). */
  const activateTab = useCallback((id: number) => {
    if (id === activeIdRef.current) return;
    setViews(tabs => {
      const snap = snapshot(tabs);
      const tab = snap.find(v => v.id === id);
      if (tab) { pendingSettingsRef.current = tab.state ?? tab.settings; setGene(tab.opened.mark ? fmtLocus(tab.opened.chrom, tab.opened.mark.start, tab.opened.mark.end) : tab.opened.geneName); }
      return snap;
    });
    viewerStateRef.current = null;
    setActiveId(id);
  }, [snapshot]);
  const closeTab = useCallback((id: number) => {
    setViews(tabs => {
      const idx = tabs.findIndex(v => v.id === id);
      const rest = tabs.filter(v => v.id !== id);
      if (id === activeIdRef.current) {
        const next = rest[Math.max(0, idx - 1)] ?? null;
        if (next) { pendingSettingsRef.current = next.state ?? next.settings; setGene(next.opened.geneName); }
        viewerStateRef.current = null;
        setActiveId(next?.id ?? null);
      }
      return rest;
    });
  }, []);
  (window as any).__sashimiDs = ds; // for debugging from the console

  const addFiles = useCallback((list: FileList | File[] | PathedFile[]) => {
    const files = Array.from(list as ArrayLike<File | PathedFile>);
    const { samples: found, unmatched, fasta: fa, fastaMissing } = pairFiles(files);
    const msgs: string[] = [];
    const added: LocalSample[] = found.map(s => ({ id: nextId.current++, ...s, pending: true }));
    for (const s of added) ds.addSample(s);
    if (added.length) setSamples(prev => [...prev, ...added]);
    // library type from the header (the aligner); the reads of the first gene opened confirm or correct it. The chip pulses until the file is open.
    for (const s of added) {
      const settle = (ev: LibraryEvidence | null) => setSamples(prev => prev.map(x => (x.id === s.id ? { ...x, pending: false, lib: ev && x.lib?.source !== 'user' && (!x.lib || x.lib.source === 'none') ? ev : x.lib } : x)));
      if (ds.getLibraryType) ds.getLibraryType(s.id).then(settle).catch(() => settle(null));
      else settle(null);
    }
    if (unmatched.length) msgs.push(`No index found for ${unmatched.join(', ')} (add the .bai / .crai file together with it)`);
    if (fa) { setFasta(fa); ds.setReference({ build, fasta: fa }); msgs.push(`Reference FASTA: ${fa.fa.name}`); }
    if (fastaMissing) msgs.push(`${fastaMissing} needs its .fai index${fastaMissing.toLowerCase().endsWith('.gz') ? ' and .gzi' : ''}`);
    setNotes(msgs);
  }, [ds, build]);

  /** A folder given through the Chromium picker or the folder input, or dropped: its files carry relative paths. */
  const addFolder = useCallback((name: string | null, handle: FSDirHandle | null, files: PathedFile[]) => {
    if (name) { setRunFolder({ name, handle }); if (handle) rememberFolder(name, handle); }
    addFiles(files);
  }, [addFiles]);
  /** The browser's folder dialog speaks of "uploading": say what really happens, at the moment the dialog opens. */
  const UPLOAD_NOTE = 'If the browser asks to "upload" the folder, that is the browser\'s own wording: nothing is sent anywhere. The page only lists the files and reads them on this computer.';
  const openFolderInput = useCallback(() => { setNotes([UPLOAD_NOTE]); folderInputRef.current?.click(); }, []);
  const chooseFolder = useCallback(async () => {
    if (!hasFileSystemAccess()) { openFolderInput(); return; }
    let h: FSDirHandle;
    try { h = await pickFolder(); } catch (e: any) {
      if (e?.name === 'AbortError') return;
      // the picker is refused in some contexts (a page opened from disk, an iframe): the plain folder input still works
      openFolderInput();
      return;
    }
    setIntake(`Listing the folder ${h.name}…`);
    setNotes([`Reading ${h.name}…`]);
    try { addFolder(h.name, h, await filesInFolder(h)); }
    catch (e: any) { setError(`Could not list the folder ${h.name}: ${e.message}`); }
    finally { setIntake(null); }
  }, [addFolder, openFolderInput]);
  const chooseFiles = useCallback(async () => {
    if (typeof (window as any).showOpenFilePicker !== 'function') { fileInputRef.current?.click(); return; }
    try {
      const picked = await pickFiles();
      for (const p of picked) fileHandlesRef.current.set(p.file.name, p.handle);
      addFiles(picked.map(p => p.file));
    } catch (e: any) { if (e?.name !== 'AbortError') setError(`Could not open the files: ${e.message}`); }
  }, [addFiles]);
  const folderInputRef = useRef<HTMLInputElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const onFolderInput = useCallback((list: FileList) => {
    setIntake(`Listing ${list.length.toLocaleString()} file${list.length === 1 ? '' : 's'}…`);
    try { const { folder, files } = filesFromFolderInput(list); addFolder(folder, null, files); } finally { setIntake(null); }
  }, [addFolder]);
  const onDropFiles = useCallback(async (dt: DataTransfer) => {
    setIntake('Reading what was dropped…');
    try {
      const { folder, folderHandle, files, fileHandles, others } = await filesFromDrop(dt);
      for (const [n, h] of fileHandles) fileHandlesRef.current.set(n, h);
      if (folder) addFolder(folder, folderHandle, files); else if (files.length) addFiles(files);
      const session = others.find(f => /\.json$/i.test(f.name));
      if (session) void loadSessionRef.current?.(session);
      else if (!files.length && !folder && others.length) setNotes([`Nothing usable in the drop (${others.map(f => f.name).join(', ')}): expected BAM/CRAM files with their index, a FASTA, a folder, or a session .json.`]);
    } catch (e: any) { setError(`Could not read the drop: ${e.message}`); }
    finally { setIntake(null); }
  }, [addFolder, addFiles]);
  /** loadSession is defined below (it depends on the samples); the drop handler reaches it through a ref. */
  const loadSessionRef = useRef<((file: File) => Promise<void>) | null>(null);

  const renameSample = useCallback((id: number, raw: string) => {
    const name = raw.trim();
    setRenaming(null);
    if (!name) return;
    ds.renameSample(id, name);
    setSamples(prev => prev.map(s => (s.id === id ? { ...s, name } : s)));
  }, [ds]);
  const sampleNames = useMemo(() => Object.fromEntries(samples.map(s => [s.id, s.name])) as Record<number, string>, [samples]);
  const sampleTypes = useMemo(() => Object.fromEntries(samples.map(s => [s.id, s.lib?.type ?? 'unknown'])) as Record<number, LibraryType>, [samples]);
  /** Evidence from the reads of a window: decisive fractions settle the type unless the user chose it; the header only breaks ties. */
  const onLibraryEvidence = useCallback((sid: number, ev: { reads: number; fraction: number; multiExon: boolean }) => {
    if (ev.reads < 200 || !ev.multiExon) return;
    const type: LibraryType | null = ev.fraction >= 0.02 ? 'rna' : ev.fraction < 0.002 ? 'dna' : null;
    if (!type) return;
    const note = `${(ev.fraction * 100).toFixed(ev.fraction < 0.01 ? 2 : 1)} % of ${ev.reads.toLocaleString()} reads spliced`;
    setSamples(prev => prev.map(s => (s.id === sid && s.lib?.source !== 'user' && (s.lib?.type !== type || s.lib.source !== 'reads') ? { ...s, lib: { type, source: 'reads', note } } : s)));
  }, []);
  /** The user decides the type of a sample (a click on its badge): RNA ↔ DNA, an undetermined sample becomes RNA. */
  const cycleLibrary = useCallback((sid: number) => setSamples(prev => prev.map(s => {
    if (s.id !== sid) return s;
    const type: LibraryType = s.lib?.type === 'rna' ? 'dna' : 'rna';
    return { ...s, lib: { type, source: 'user', note: 'chosen by the user' } };
  })), []);
  const libBadge = (lib: LibraryEvidence | undefined) => (lib?.type === 'rna' ? 'RNA' : lib?.type === 'dna' ? 'DNA' : '?');

  const removeSample = useCallback((id: number) => { ds.removeSample(id); setSamples(prev => prev.filter(s => s.id !== id)); }, [ds]);

  /** Sets the variants of interest: the data source serves them for every sample and the viewer fetches them again. */
  const setKnown = useCallback((list: KnownVariant[]) => { ds.knownVariants = list; knownSeq.current++; setKnownVars(list); }, [ds]);
  const addKnown = useCallback(() => {
    const text = knownLocus.trim();
    if (!text) return;
    const v = variantOfInterest(text, knownLabel, `user${nextKnownId.current++}`);
    if (!v) { setKnownError(`"${text}" is not a locus (chr17:43,094,464 or chr17:43,094,464-43,094,470), an HGVS genomic notation or a VCF-like line`); return; }
    setKnownError(null);
    setKnown([...knownVars, v]);
    setKnownLocus(''); setKnownLabel('');
  }, [knownLocus, knownLabel, knownVars, setKnown]);
  const removeKnown = useCallback((id: string) => setKnown(knownVars.filter(v => v.id !== id)), [knownVars, setKnown]);

  const changeBuild = (b: GenomeBuild) => { setBuild(b); ds.setReference({ build: b, fasta }); };

  /** Open the window a deep link asked for: the gene at the locus (coding first), else the `gene` parameter, with the variant pinned. */
  const openLink = useCallback(async () => {
    if (!LINK) return;
    setBusy(true); setError(null);
    try {
      const { view, mark } = LINK;
      const genes = await ds.getRegionGenes(view.chrom, view.start, view.end);
      const ov = (g: { start: number; end: number }) => Math.min(g.end, view.end) - Math.max(g.start, view.start);
      let best = [...genes].sort((a, b) => Number(b.biotype === 'protein_coding') - Number(a.biotype === 'protein_coding') || ov(b) - ov(a))[0];
      if (!best && LINK.gene) {
        const t = await ds.getTranscript(LINK.gene, LINK.gene.toUpperCase().startsWith('ENSG') ? LINK.gene : undefined);
        best = { gene_name: t.gene_name, start: t.start, end: t.end } as any;
      }
      if (!best) throw new Error(`no RefSeq gene at ${view.chrom}:${view.start.toLocaleString()}-${view.end.toLocaleString()} (add gene=SYMBOL to the link)`);
      openNew({ geneName: best.gene_name, chrom: view.chrom, start: best.start, end: best.end, view: { start: view.start, end: view.end }, mark: { start: mark.start, end: mark.end }, reads: LINK.reads });
    } catch (e: any) {
      setError(`Could not open the linked locus: ${e.message}. Check the genome build (${LINK.build}) and that api.genome.ucsc.edu (or rest.ensembl.org) is reachable.`);
    }
    setBusy(false);
  }, [ds, openNew]);

  // A link opens the browser right away; files added afterwards become tracks
  const autoOpened = useRef(false);
  useEffect(() => {
    if (LINK && !autoOpened.current) { autoOpened.current = true; openLink(); }
  }, [openLink]);

  const open = useCallback(async () => {
    const q = gene.trim();
    if (!q) return;
    if (LINK && q === LINK_TEXT) return openLink();
    setBusy(true); setError(null);
    try {
      const locus = parseLocus(q);
      if (locus) {
        // coordinates: open the gene at the locus (coding first, largest overlap) and show the locus itself
        const genes = await ds.getRegionGenes(locus.chrom, locus.start, locus.end);
        const ov = (g: { start: number; end: number }) => Math.min(g.end, locus.end) - Math.max(g.start, locus.start);
        const best = [...genes].sort((a, b) => Number(b.biotype === 'protein_coding') - Number(a.biotype === 'protein_coding') || ov(b) - ov(a))[0];
        if (!best) throw new Error(`no RefSeq gene at ${q}`);
        const point = locus.start === locus.end;
        openNew({ geneName: best.gene_name, chrom: locus.chrom, start: best.start, end: best.end, view: { start: point ? Math.max(1, locus.start - 500) : locus.start, end: point ? locus.start + 500 : locus.end }, mark: { start: locus.start, end: locus.end } });
      } else {
        const t = await ds.getTranscript(q, q.toUpperCase().startsWith('ENSG') ? q : undefined);
        openNew({ geneName: t.gene_name, chrom: t.chrom, start: t.start, end: t.end });
      }
    } catch (e: any) {
      setError(`Gene lookup failed: ${e.message}. Check the symbol, the genome build and that api.genome.ucsc.edu (or rest.ensembl.org) is reachable.`);
    }
    setBusy(false);
  }, [gene, ds, openLink, openNew]);

  // Options the viewer starts with: a reopened tab's snapshot or a loaded session, else the live options of the
  // viewer being remounted (a sample added, the build changed)
  const viewerInit: Partial<ViewerSettings> | undefined = pendingSettingsRef.current ?? (viewerStateRef.current ? { ...viewerStateRef.current } : undefined);
  // order-independent: promoting another sample to primary keeps the viewer (and its view) mounted
  const viewerKey = useMemo(() => `${activeId}|${opened?.geneName}|${opened?.view ? `${opened.view.start}-${opened.view.end}` : ''}|${opened?.mark ? `${opened.mark.start}-${opened.mark.end}` : ''}|${[...samples.map(s => s.id)].sort((a, b) => a - b).join(',')}|${build}|${fasta?.fa.name || ''}|${sessionSeq}`, [activeId, opened, samples, build, fasta, sessionSeq]);
  /** Every tab with a full state (the active one refreshed from the viewer, never-opened ones from their region and options), and the index of the active one. */
  const currentViews = useCallback(() => {
    const tabs = snapshot(views).map(v => ({ ...v, state: stateOfTab(v) }));
    return { tabs, activeIndex: Math.max(0, tabs.findIndex(v => v.id === activeIdRef.current)) };
  }, [views, snapshot]);
  const makePrimary = useCallback((id: number) => setSamples(prev => [...prev.filter(s => s.id === id), ...prev.filter(s => s.id !== id)]), []);

  // ---- Save / load a session (JSON) ----
  const saveSession = useCallback(() => {
    const { tabs, activeIndex } = currentViews();
    const session = buildSession({ build, folder: runFolder?.name ?? null, samples, fasta, state: opened ? viewerStateRef.current : null, views: tabs.map(t => ({ label: t.label, state: t.state! })), activeView: activeIndex, knownVariants: knownVars });
    if (hasFileSystemAccess()) {
      const entries = [...samples.flatMap(s => [s.file, s.index]), ...(fasta ? [fasta.fa, fasta.fai, ...(fasta.gzi ? [fasta.gzi] : [])] : [])]
        .map(f => ({ name: f.name, size: f.size, handle: fileHandlesRef.current.get(f.name)! })).filter(e => e.handle);
      rememberFiles(entries);
    }
    const name = (sessionName.trim() || defaultSessionName(session.gene?.name)).replace(/\.json$/i, '') + '.json';
    const url = URL.createObjectURL(new Blob([JSON.stringify(session, null, 2)], { type: 'application/json' }));
    const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    setNotes([`Session saved as ${name} (${session.samples.length} sample${session.samples.length === 1 ? '' : 's'}, ${tabs.length} view${tabs.length === 1 ? '' : 's'}). Files are recorded by name: add them again when loading.`]);
  }, [build, samples, fasta, opened, sessionName, runFolder, currentViews]);

  // ---- Export: this page with the data of every view embedded, for a reader without the alignment files ----
  const exportHtml = useCallback(async (opts: ExportOptions) => {
    setBusy(true); setError(null);
    let nReads = 0, nReadSets = 0;
    const skipped: string[] = [];
    let done = 0, total = 1;
    const progress = (step: string) => setExportProgress({ title: 'Exporting HTML…', step, done, total });
    try {
      const { tabs, activeIndex } = currentViews();
      if (!tabs.length) throw new Error('open a gene first');
      // one step per view (annotation), per sample coverage, per sample reads when exported, plus the file itself
      total = 1 + tabs.reduce((n, t) => n + 1 + samples.length + (t.state!.reads && t.state!.view.end - t.state!.view.start + 1 <= READS_MAX_VIEW_BP ? samples.length : 0), 0);
      progress('Preparing…');
      const session = buildSession({ build, folder: null, samples, fasta, state: tabs[activeIndex].state, views: tabs.map(t => ({ label: t.label, state: t.state! })), activeView: activeIndex, knownVariants: knownVars });
      const evs: EmbeddedView[] = [];
      for (const t of tabs) {
        const st = t.state!;
        setNotes([`Exporting ${t.label}…`]);
        progress(`${t.label}: gene models`);
        const vs = st.view.start - 1, ve = st.view.end, span = ve - vs;
        // the window exported around the view: as shown, with a half-width margin, or what the viewer itself fetched
        const margin = opts.window === 'view' ? 0 : opts.window === 'margin' ? Math.floor(span / 2) : Math.min(span, Math.max(0, Math.floor((MAX_FETCH_BP - span) / 2)));
        const ws = Math.max(0, vs - margin), we = ve + margin;
        const hint = { chrom: st.gene.chrom, start: st.gene.start, end: st.gene.end };
        const transcript = await ds.getTranscript(st.gene.name, st.gene.id, hint);
        let allTranscripts: EmbeddedView['allTranscripts']; try { allTranscripts = await ds.getAllTranscripts(st.gene.name, st.gene.id, hint); } catch { /* optional */ }
        let regionGenes: EmbeddedView['regionGenes']; try { regionGenes = await ds.getRegionGenes(st.gene.chrom, ws + 1, we, transcript.gene_name); } catch { /* optional */ }
        // boundary-spanning reads are counted at every exon boundary of every model, so the reader may pick any reference transcript
        const intronStarts = new Set<number>(), intronEnds = new Set<number>();
        const exonSets = [toTxModel(transcript).exons.map(e => ({ start: e.start, end: e.end })), ...(allTranscripts?.transcripts ?? []).map(m => m.exons.map(e => ({ start: e.start - 1, end: e.end })))];
        for (const ex of exonSets) { const sorted = [...ex].sort((a, b) => a.start - b.start); for (let i = 0; i + 1 < sorted.length; i++) { intronStarts.add(sorted[i].end); intronEnds.add(sorted[i + 1].start); } }
        const coverage: Record<string, EncodedCoverage | EncodedCoverageV2> = {};
        done++;
        for (const smp of samples) {
          progress(`${t.label}: coverage and junctions of ${smp.name}`);
          try {
            const c = await ds.getCoverage(smp.id, st.gene.chrom, ws, we, st.uniqueOnly, { intronStarts: [...intronStarts], intronEnds: [...intronEnds] }, { core: { start: vs, end: ve }, maxReads: 250_000 });
            coverage[String(smp.id)] = await encodeCoverageV2(c, { start: ws, end: we });
          } catch (e: any) {
            coverage[String(smp.id)] = { start: ws, len: [], depth: [], junctions: [], window: { start: ws, end: we }, error: e?.message || String(e) };
          }
          done++;
        }
        // reads of every loaded sample when the view shows its reads track (window and cap from the dialog)
        let reads: Record<string, EncodedReadsV2> | undefined;
        if (st.reads) {
          if (span > READS_MAX_VIEW_BP) skipped.push(`${t.label} (window of ${(span / 1000).toFixed(0)} kb, above the ${READS_MAX_VIEW_BP / 1000} kb reads limit)`);
          else {
            const half = opts.window === 'view' ? 0 : opts.window === 'margin' ? Math.floor(span / 2) : Math.floor((READS_MAX_VIEW_BP - span) / 2);
            const rs = Math.max(0, vs - half), re = ve + half;
            reads = {};
            for (const smp of samples) {
              setNotes([`Exporting ${t.label}: reads of ${smp.name}…`]);
              progress(`${t.label}: reads of ${smp.name}`);
              try {
                const r = await ds.getReads(smp.id, st.gene.chrom, rs, re, st.uniqueOnly, READS_CAPS[opts.readsCap], 'reads', 1, 0.05);
                reads[String(smp.id)] = await encodeReadsV2({ window: { start: rs, end: re }, total: r.total, reads: r.reads, reference: r.reference, reference_source: r.reference_source });
                nReads += r.reads.length; nReadSets++;
              } catch (e: any) { skipped.push(`${t.label} / ${smp.name}: ${e?.message || e}`); }
              done++;
            }
          }
        }
        evs.push({ label: t.label, gene: { name: st.gene.name, id: st.gene.id, chrom: st.gene.chrom, start: st.gene.start, end: st.gene.end }, transcript, allTranscripts, regionGenes, window: { chrom: st.gene.chrom, start: ws, end: we }, uniqueOnly: st.uniqueOnly, coverage, reads });
      }
      const payload: EmbeddedExport = {
        app: EMBEDDED_APP, version: EMBEDDED_VERSION, saved: new Date().toISOString(), build,
        samples: samples.map(s => ({ id: s.id, name: s.name, kind: s.kind, file: s.file.name, index: s.index.name, library: s.lib })),
        session, views: evs, knownVariants: knownVars,
      };
      progress('Writing the file…');
      const html = await buildExportHtml(payload);
      done = total; progress('Starting the download…');
      const name = (sessionName.trim() || defaultSessionName(tabs[activeIndex].state!.gene.name)).replace(/\.json$/i, '').replace(/\.html$/i, '') + '.html';
      const url = URL.createObjectURL(new Blob([html], { type: 'text/html' }));
      const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
      setNotes([
        `Exported ${name} (${(html.length / 1048576).toFixed(1)} MB): ${evs.length} view${evs.length === 1 ? '' : 's'}, ${samples.length} sample${samples.length === 1 ? '' : 's'}${nReadSets ? `, ${nReads.toLocaleString()} reads in ${nReadSets} reads track${nReadSets === 1 ? '' : 's'} (with reference bases and mismatches; read names replaced by numbers)` : ''}. The file opens in any browser without the alignment files.`,
        ...(skipped.length ? [`Reads not exported for ${skipped.join('; ')}`] : []),
      ]);
    } catch (e: any) {
      setError(`Export failed: ${e?.message || String(e)}`);
    }
    setExportProgress(null);
    setBusy(false);
  }, [build, samples, fasta, sessionName, ds, currentViews]);
  const viewsWithReads = views.filter(v => (v.id === activeId ? viewerStateRef.current?.reads : stateOfTab(v).reads)).length;

  // ---- Every view on one SVG page: each tab is shown in turn, its plot serialised once loaded, then the current view comes back ----
  const exportAllSvg = useCallback(async () => {
    if (!views.length) return;
    setBusy(true); setError(null);
    const original = activeIdRef.current;
    const plotEl = () => document.querySelector<SVGSVGElement>('svg[data-sashimi-plot]');
    const sleep = (ms: number) => new Promise(r => setTimeout(r, ms));
    /** waits for a plot element other than `prev` that has finished loading, coverage and reads alike, in two consecutive polls (or gives up after 30 s) */
    const quiet = (el: SVGSVGElement) => el.dataset.loading === '0' && !/loading…|updating…|loading reads/.test(el.textContent || '');
    const settled = async (prev: SVGSVGElement | null) => {
      const t0 = Date.now();
      let calm = 0;
      for (;;) {
        const el = plotEl();
        if (el && el !== prev && quiet(el)) { if (++calm >= 2) return el; } else calm = 0;
        if (Date.now() - t0 > 30000) return el;
        await sleep(200);
      }
    };
    const plots: { title: string; svg: string }[] = [];
    try {
      let prev: SVGSVGElement | null = null;
      for (let i = 0; i < views.length; i++) {
        const v = views[i];
        setNotes([`SVG of ${v.label} (${i + 1}/${views.length})…`]);
        setExportProgress({ title: 'Exporting SVG…', step: `${v.label}: showing the view and capturing its plot`, done: i, total: views.length + 1 });
        if (v.id !== activeIdRef.current) { activateTab(v.id); await sleep(50); }
        const el = await settled(v.id === original && i === 0 ? null : prev);
        if (!el) { setError(`No plot for ${v.label}`); continue; }
        prev = el;
        const st = v.id === activeIdRef.current ? viewerStateRef.current : v.state;
        const where = st ? fmtLocus(st.view.chrom, st.view.start, st.view.end) : '';
        plots.push({ title: where && !v.label.includes(where) ? `${v.label} · ${where}` : v.label, svg: serializePlotSvg(el) });
      }
      if (original != null && original !== activeIdRef.current) activateTab(original);
      const stem = (sessionName.trim() || defaultSessionName(opened?.geneName)).replace(/\.(json|html)$/i, '');
      if (!plots.length) throw new Error('no plot captured');
      setExportProgress({ title: 'Exporting SVG…', step: 'Writing the file…', done: views.length, total: views.length + 1 });
      const blob = new Blob([stackSvgs(plots)], { type: 'image/svg+xml;charset=utf-8' });
      const name = `${safeFileName(stem)}-views.svg`;
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
      setNotes([`${name} saved: ${plots.length} view${plots.length === 1 ? '' : 's'} on one SVG page, each as drawn with its own options, under its title.`]);
    } catch (e: any) {
      setError(`SVG export failed: ${e?.message || String(e)}`);
    }
    setExportProgress(null);
    setBusy(false);
  }, [views, activateTab, sessionName, opened]);

  /** Applies a loaded session with the files present: names, order, options, then the gene and window. */
  const applySession = useCallback((session: SessionFile) => {
    if (session.knownVariants) setKnown(session.knownVariants.map(k => variantOfInterest(k.text, k.label, `user${nextKnownId.current++}`)).filter((v): v is KnownVariant => !!v).map((v, i) => ({ ...v, id: `session${i + 1}` })));
    const { matched } = matchSession(session, samples);
    for (const { entry, sample } of matched) if (sample.name !== entry.name) ds.renameSample(sample.id, entry.name);
    const renamed = samples.map(s => { const m = matched.find(x => x.sample.id === s.id); return m ? { ...s, name: m.entry.name, lib: m.entry.library ?? s.lib } : s; });
    const ordered = [...matched.map(m => renamed.find(s => s.id === m.sample.id)!), ...renamed.filter(s => !matched.some(m => m.sample.id === s.id))];
    setSamples(ordered);
    setSessionSeq(n => n + 1);
    setPendingSession(null);
    setReopenState(null);
    const list = session.views ?? (session.gene ? [{ label: session.gene.name, gene: session.gene, viewer: session.viewer }] : []);
    const tabs: ViewTab[] = list.map(v => ({
      id: tabSeq.current++, label: v.label, state: null,
      opened: { geneName: v.gene.name, geneId: v.gene.id, chrom: v.gene.chrom, start: v.gene.start, end: v.gene.end, view: v.gene.view, mark: v.gene.mark ?? null },
      settings: viewerSettingsOf({ viewer: v.viewer }, ordered),
    }));
    const act = tabs[Math.min(session.activeView ?? 0, Math.max(0, tabs.length - 1))] ?? null;
    viewerStateRef.current = null;
    pendingSettingsRef.current = act?.settings;
    setViews(tabs);
    setActiveId(act?.id ?? null);
    if (act) setGene(act.opened.geneName);
    setNotes([`Session loaded: ${matched.length}/${session.samples.length} sample${session.samples.length === 1 ? '' : 's'}, ${tabs.length} view${tabs.length === 1 ? '' : 's'}.`]);
  }, [samples, ds]);

  /**
   * Gets a session's files back by itself (Chromium): from the remembered run folder by relative path, else from
   * remembered file bookmarks. `ask` lets the browser show its permission prompt (needs a user gesture).
   */
  const reopenSession = useCallback(async (session: SessionFile, ask: boolean) => {
    if (!hasFileSystemAccess()) { setReopenState(session.folder ? { folder: session.folder, ready: false } : null); return; }
    const { missing } = matchSession(session, samples);
    if (!missing.length) return;
    const out: PathedFile[] = [];
    let dir: FSDirHandle | null = null, needPermission = false;
    if (session.folder) setIntake(`Looking for the session's files in ${session.folder}…`);
    try {
    if (session.folder) {
      dir = await recallFolder(session.folder);
      if (dir) {
        if (await permitted(dir, ask)) {
          setRunFolder({ name: session.folder, handle: dir });
          for (const m of missing) {
            if (!m.path || !m.indexPath) continue;
            const [f, i] = await Promise.all([fileInFolder(dir, m.path), fileInFolder(dir, m.indexPath)]);
            if (f && i) out.push({ file: f, path: m.path }, { file: i, path: m.indexPath });
          }
        } else needPermission = true;
      }
    }
    // individual bookmarks for what the folder did not give
    const got = new Set(out.map(x => x.file.name));
    for (const m of missing) {
      if (got.has(m.file)) continue;
      const [fh, ih] = await Promise.all([recallFile(m.file), recallFile(m.index)]);
      if (!fh || !ih) continue;
      if (await permitted(fh, ask) && await permitted(ih, ask)) {
        out.push({ file: await fh.getFile() }, { file: await ih.getFile() });
        fileHandlesRef.current.set(m.file, fh); fileHandlesRef.current.set(m.index, ih);
      } else needPermission = true;
    }
    if (out.length) addFiles(out);
    const stillMissing = missing.filter(m => !out.some(x => x.file.name === m.file));
    setReopenState(stillMissing.length ? {
      folder: session.folder ?? '', ready: needPermission,
      note: needPermission ? 'this browser remembers the files: click Reopen to allow access' : session.folder ? (dir ? 'some files were not found in the remembered folder' : 'the folder is not remembered by this browser yet') : undefined,
    } : null);
    } finally { setIntake(null); }
  }, [samples, addFiles]);

  const loadSession = useCallback(async (file: File) => {
    setNotes([`Loading session ${file.name}…`]);
    setError(null);
    try {
      const session = parseSession(await file.text());
      setSessionName(file.name); setSessionNameEdited(true);
      if (session.build !== build) { setBuild(session.build); ds.setReference({ build: session.build, fasta }); }
      setPendingSession({ file: session, name: file.name });
      void reopenSession(session, true);
    } catch (e: any) {
      setError(`Could not load the session ${file.name}: ${e.message}`);
    }
  }, [build, ds, fasta]);

  // An exported page opens its views at once (its samples are embedded)
  const embeddedOpened = useRef(false);
  useEffect(() => {
    if (EMBEDDED && !embeddedOpened.current) {
      embeddedOpened.current = true; applySession(EMBEDDED.session);
      // the active view's blocks first, the others in idle moments: the page never waits for a decode
      if (ds instanceof EmbeddedDataSource) ds.prefetchAll(EMBEDDED.session.activeView ?? 0);
    }
  }, [applySession]);

  loadSessionRef.current = loadSession;

  // A pending session applies itself as soon as every one of its files is present (or right away when it needs none)
  useEffect(() => {
    if (!pendingSession) return;
    const { missing } = matchSession(pendingSession.file, samples);
    if (!missing.length) applySession(pendingSession.file);
  }, [pendingSession, samples, applySession]);

  // Default file name follows the gene until the user types one
  useEffect(() => { if (!sessionNameEdited) setSessionName(defaultSessionName(opened?.geneName)); }, [opened?.geneName, sessionNameEdited]);

  /** Drops are accepted anywhere on the page: files, a folder, or a session .json. */
  const onDrop = (e: React.DragEvent) => { e.preventDefault(); e.stopPropagation(); void onDropFiles(e.dataTransfer); };

  return (
    <div className="min-h-screen bg-gray-100 text-gray-900" onDragOver={e => e.preventDefault()} onDrop={onDrop}>
      {DEV && (
        <div className="px-5 py-1.5 text-xs font-semibold text-white flex flex-wrap items-center gap-x-3 gap-y-1"
          style={{ background: 'repeating-linear-gradient(135deg, #b91c1c 0 14px, #dc2626 14px 28px)' }}>
          <span className="px-1.5 py-0.5 rounded bg-white text-red-700 tracking-wider">DEV MODE</span>
          <span>Unreleased build from branch <code>{DEV.branch}</code>{DEV.sha ? <> · commit <code>{DEV.sha}</code></> : null}{DEV.date ? ` · built ${DEV.date}` : ''} · not for clinical use.</span>
          <a href={DEV.stable} className="underline ml-auto">Go to the stable version →</a>
        </div>
      )}
      <header className="bg-white border-b border-gray-200 px-5 py-3 space-y-2">
        {/* Row 1: title, build, gene search. Row 2: files and sample chips, which may wrap over several lines without moving the search. */}
        <div className="flex flex-wrap items-center gap-x-6 gap-y-2">
        <div className="flex items-center gap-2.5 flex-1 min-w-[320px]">
          <Logo size={34} />
          <div className="min-w-0">
            <h1 className="text-lg font-bold leading-tight">Sashimi <span className="font-normal">viewer</span></h1>
            <p className="text-xs text-gray-500">Files are read in your browser and never uploaded. Gene models (RefSeq, UCSC API) and reference bases come from the network unless you add a FASTA.</p>
            <p className="text-xs font-medium text-amber-700" role="note">⚠ Check the HGVS nomenclature, the predicted transcript, amino-acid changes and the NMD verdict before anything from it goes into a clinical report.</p>
          </div>
        </div>
        <div className="flex items-center gap-4 shrink-0 ml-auto">
        <label className="flex items-center gap-1 text-xs text-gray-600">Build
          <select value={build} onChange={e => changeBuild(e.target.value as GenomeBuild)} className="border border-gray-300 rounded px-1 py-0.5 text-xs bg-white">
            <option value="GRCh38">GRCh38 / hg38</option>
            <option value="GRCh37">GRCh37 / hg19</option>
          </select>
        </label>
        <form onSubmit={e => { e.preventDefault(); open(); }} className="flex items-center gap-1">
          <input value={gene} onChange={e => setGene(e.target.value)} placeholder="Gene, ENSG or chr:pos…" title="A gene symbol, an ENSG id, or coordinates (chr17:43,094,464 or chr17:43,000,000-43,100,000: the gene at the locus is opened)" className="border border-gray-300 rounded px-2 py-1 text-sm w-48 bg-white" />
          <button type="submit" disabled={busy || !gene.trim()} className="px-3 py-1 text-sm rounded bg-indigo-600 text-white disabled:opacity-40 hover:bg-indigo-700">{busy ? '…' : 'Open'}</button>
        </form>
        </div>
        </div>
        <div className="flex flex-wrap items-center gap-x-3 gap-y-1.5">
        <button onClick={chooseFolder} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium"
          title="Choose the run folder: every BAM/CRAM (with its index) and FASTA inside it is listed, without being read, and sessions record the files by their path inside this folder. Chrome and Edge remember the folder so a session reopens it after one click. Nothing is uploaded: if the browser's dialog says so, that is its own wording for letting this page read the files.">
          {runFolder ? `Run folder · ${runFolder.name}` : '+ Run folder…'}
        </button>
        <input ref={folderInputRef} type="file" className="hidden" {...({ webkitdirectory: '', directory: '' } as any)} onChange={e => { if (e.target.files) onFolderInput(e.target.files); e.target.value = ''; }} />
        <button onClick={chooseFiles} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium" title="Add individual files: BAM/CRAM with their .bai/.crai, a FASTA with its .fai (and .gzi)">
          + Files…
        </button>
        <input ref={fileInputRef} type="file" multiple className="hidden" accept=".bam,.bai,.cram,.crai,.fa,.fasta,.fna,.gz,.fai,.gzi" onChange={e => { if (e.target.files) addFiles(e.target.files); e.target.value = ''; }} />
        {intake && (
          <span className="flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs border border-indigo-300 bg-indigo-50 text-indigo-800" role="status" data-intake>
            <span className="inline-block w-3 h-3 rounded-full border-2 border-indigo-500 border-t-transparent animate-spin" />
            {intake}
          </span>
        )}
        <div className="flex flex-wrap items-center gap-1.5">
          {samples.map((s, i) => (
            <span key={s.id} className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border ${i === 0 ? 'bg-indigo-50 border-indigo-300 text-indigo-800' : 'bg-gray-50 border-gray-300 text-gray-700 hover:border-indigo-300 cursor-pointer'}`}
              title={`${s.file.name} · ${s.embedded ? 'data embedded in this exported page' : `${(s.file.size / 1e9).toFixed(2)} GB`} · ${s.kind.toUpperCase()}${i === 0 ? ' · primary sample' : ' · click to make it the primary sample'} · double-click to rename`}
              onClick={() => { if (i !== 0 && renaming?.id !== s.id) makePrimary(s.id); }}
              onDoubleClick={e => { e.stopPropagation(); setRenaming({ id: s.id, value: s.name }); }}>
              {i === 0 && <span title="primary sample">★</span>}
              {s.pending && <span className="inline-block w-2.5 h-2.5 rounded-full border-2 border-indigo-500 border-t-transparent animate-spin" title="Opening the file: reading its header and index…" />}
              <button onClick={e => { e.stopPropagation(); cycleLibrary(s.id); }}
                className={`px-1 rounded text-[9px] font-bold leading-4 ${s.pending ? 'bg-indigo-100 text-indigo-700 animate-pulse' : s.lib?.type === 'dna' ? 'bg-slate-700 text-white' : s.lib?.type === 'rna' ? 'bg-emerald-600 text-white' : 'bg-gray-300 text-gray-700'}`}
                title={`${s.lib?.type === 'dna' ? 'Genomic DNA' : s.lib?.type === 'rna' ? 'RNA-seq' : 'Library type not determined yet (treated as RNA)'} · ${s.lib?.note ?? 'decided from the header and the first gene opened'} · click to switch (RNA-seq shows junction arcs and usage; DNA shows depth and reads only)`}>
                {s.pending ? '…' : libBadge(s.lib)}
              </button>
              {renaming?.id === s.id ? (
                <input autoFocus value={renaming.value} onChange={e => setRenaming({ id: s.id, value: e.target.value })}
                  onBlur={() => renameSample(s.id, renaming.value)}
                  onKeyDown={e => { if (e.key === 'Enter') renameSample(s.id, renaming.value); else if (e.key === 'Escape') setRenaming(null); }}
                  onClick={e => e.stopPropagation()} className="w-32 px-1 py-0 text-xs rounded border border-indigo-300 bg-white text-gray-900" title="Enter to confirm, Esc to cancel" />
              ) : s.name}
              <button onClick={e => { e.stopPropagation(); setRenaming({ id: s.id, value: s.name }); }} className="text-gray-400 hover:text-indigo-600" title="Rename">✎</button>
              <button onClick={e => { e.stopPropagation(); removeSample(s.id); }} className="text-gray-400 hover:text-red-500" title="Remove">×</button>
            </span>
          ))}
          {fasta && <span className="px-2 py-0.5 rounded-full text-xs border bg-emerald-50 border-emerald-300 text-emerald-800" title={fasta.fa.name}>FASTA · {fasta.fa.name}</span>}
          {samples.some(s => s.pending) && (
            <span className="flex items-center gap-1 text-xs text-indigo-700" role="status">
              <span className="inline-block w-3 h-3 rounded-full border-2 border-indigo-500 border-t-transparent animate-spin" />
              Processing {samples.filter(s => s.pending).length === 1 ? 'the file' : `${samples.filter(s => s.pending).length} files`}…
            </span>
          )}
        </div>
        <div className="flex items-center gap-1.5 ml-auto" title="A session file (JSON) records the alignment files by name, the sample names and order, the FASTA, the gene and window, and every option of the viewer. Load it later and add the same files again.">
          <span className="text-xs text-gray-600">Session</span>
          <input value={sessionName} onChange={e => { setSessionName(e.target.value); setSessionNameEdited(true); }} spellCheck={false}
            className="border border-gray-300 rounded px-2 py-0.5 text-xs w-64 bg-white font-mono" title="File name of the session to save (.json)" />
          <button onClick={saveSession} disabled={!samples.length && !opened} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium disabled:opacity-40" title="Download the session as a JSON file">Save session</button>
          <label className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 cursor-pointer font-medium" title="Load a session JSON file, then add the alignment files it names">
            Load session
            <input type="file" accept=".json,application/json" className="hidden" onChange={e => { const f = e.target.files?.[0]; if (f) loadSession(f); e.target.value = ''; }} />
          </label>
          <button onClick={() => setExportDialog({ window: 'margin', readsCap: 'shown' })} disabled={busy || !opened || pageIsUnbuilt()} className="px-3 py-1 text-xs rounded border border-emerald-300 bg-emerald-50 text-emerald-900 hover:bg-emerald-100 font-medium disabled:opacity-40"
            title={pageIsUnbuilt() ? 'The export needs the built viewer (sashimi-viewer.html), not the development page.' : 'Download a copy of this viewer with the data of every registered view embedded (coverage, junctions, retention counts of each window, gene models, options and groups, and the reads of the views whose reads track is on). Anyone can open it in a browser without the alignment files and switch between the views, zoom and pan inside them.'}>
            {busy ? '…' : 'Export HTML'}
          </button>
          <button onClick={exportAllSvg} disabled={busy || !views.length} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium disabled:opacity-40"
            title="Save every registered view on one SVG page, stacked under their titles (vector, publication-ready). Each tab is shown in turn while its plot is captured; the current view comes back at the end. The SVG button inside the plot saves the current view alone.">
            {views.length > 1 ? `SVG · ${views.length} views` : 'SVG'}
          </button>
        </div>
        </div>
        <div className="flex flex-wrap items-center gap-x-3 gap-y-1.5" title="Variants of interest: each is drawn on every view (a panel under the transcript with its label, guide lines through the tracks) and can be jumped to from the viewer's Known variants menu. Saved with the session and the exported page.">
          <span className="text-xs text-gray-600">Known variants</span>
          <form onSubmit={e => { e.preventDefault(); addKnown(); }} className="flex items-center gap-1">
            <input value={knownLocus} onChange={e => { setKnownLocus(e.target.value); setKnownError(null); }} placeholder="chr17:43,094,464" spellCheck={false}
              className={`border rounded px-2 py-0.5 text-xs w-44 bg-white font-mono ${knownError ? 'border-red-400' : 'border-gray-300'}`}
              title="Where the variant is: chr:position, chr:start-end, an HGVS genomic notation (chr17:g.43094464A>G, NC_000017.11:g.43094464A>G) or a VCF-like line (chr17 43094464 A G)" />
            <input value={knownLabel} onChange={e => setKnownLabel(e.target.value)} placeholder="label, e.g. BRCA1 p.Glu23Asp" className="border border-gray-300 rounded px-2 py-0.5 text-xs w-52 bg-white"
              title="Text drawn next to the variant (gene and change, sample, anything); the locus itself when empty" />
            <button type="submit" disabled={!knownLocus.trim()} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium disabled:opacity-40">+ Add</button>
          </form>
          {knownError && <span className="text-xs text-red-600">{knownError}</span>}
          {knownVars.map(v => (
            <span key={v.id} className="flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border bg-rose-50 border-rose-300 text-rose-900" title={`${v.text}${v.label !== v.text ? ` · ${v.label}` : ''}`}>
              <span className="font-medium">{v.label}</span>
              <span className="font-mono text-rose-700/80">{v.chrom || '?'}:{(v.start + 1).toLocaleString()}{v.end > v.start + 1 ? `-${v.end.toLocaleString()}` : ''}</span>
              <button onClick={() => removeKnown(v.id)} className="text-rose-400 hover:text-red-600" title="Remove">×</button>
            </span>
          ))}
        </div>
        {views.length > 1 && (
          <div className="flex flex-wrap items-center gap-1.5" title="Registered views: each gene or locus opened from the search box above is kept as a tab with its own options. Click one to reopen it, × to forget it. Views are saved with the session and in the HTML export.">
            <span className="text-xs text-gray-600">Views</span>
            {views.map(v => (
              <span key={v.id} onClick={() => activateTab(v.id)}
                className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border cursor-pointer ${v.id === activeId ? 'bg-indigo-600 border-indigo-600 text-white' : 'bg-white border-gray-300 text-gray-700 hover:border-indigo-400 hover:bg-indigo-50'}`}
                title={`${v.opened.geneName} · ${fmtLocus(v.opened.chrom, v.opened.view?.start ?? v.opened.start, v.opened.view?.end ?? v.opened.end)}${v.id === activeId ? ' · shown' : ' · click to reopen with its options'}`}>
                {v.label}
                <button onClick={e => { e.stopPropagation(); closeTab(v.id); }} className={v.id === activeId ? 'text-indigo-200 hover:text-white' : 'text-gray-400 hover:text-red-500'} title="Forget this view">×</button>
              </span>
            ))}
          </div>
        )}
      </header>
      {pendingSession && (() => {
        const { matched, missing } = matchSession(pendingSession.file, samples);
        return (
          <div className="mx-5 mt-2 px-3 py-2 text-xs rounded-lg border border-indigo-200 bg-indigo-50 text-indigo-900 flex flex-wrap items-center gap-x-3 gap-y-1">
            <span className="font-semibold">Session {pendingSession.name}</span>
            <span>{matched.length}/{pendingSession.file.samples.length} files present{missing.length ? ` · missing: ${missing.map(m => m.file).join(', ')}` : ''}</span>
            {missing.length > 0 && reopenState?.ready && (
              <button onClick={() => reopenSession(pendingSession.file, true)} className="px-2 py-0.5 rounded bg-indigo-600 text-white hover:bg-indigo-700 font-medium">Reopen the files</button>
            )}
            {missing.length > 0 && !reopenState?.ready && (
              <span className="flex items-center gap-1.5">
                {pendingSession.file.folder ? <>drop the folder <b>{pendingSession.file.folder}</b> on the page, or</> : 'drop the files on the page, or'}
                <button onClick={chooseFolder} className="px-2 py-0.5 rounded border border-indigo-300 bg-white hover:bg-indigo-100 font-medium">{pendingSession.file.folder ? 'Choose the folder' : 'Choose a folder'}</button>
                <button onClick={chooseFiles} className="px-2 py-0.5 rounded border border-indigo-300 bg-white hover:bg-indigo-100 font-medium">Choose files</button>
                <span className="text-indigo-500">(the files stay on this computer; a browser dialog mentioning an upload only grants this page read access)</span>
              </span>
            )}
            {reopenState?.note && <span className="text-indigo-600">{reopenState.note}</span>}
            {matched.length > 0 && missing.length > 0 && <button onClick={() => applySession(pendingSession.file)} className="px-2 py-0.5 rounded border border-indigo-300 bg-white hover:bg-indigo-100 font-medium">Open with the files present</button>}
            <button onClick={() => setPendingSession(null)} className="text-indigo-500 hover:text-indigo-800" title="Forget this session">×</button>
          </div>
        );
      })()}
      {EMBEDDED && (
        <div className="mx-5 mt-2 px-3 py-2 text-xs rounded-lg border border-emerald-200 bg-emerald-50 text-emerald-900">
          <b>Exported viewer</b> · {EMBEDDED.views.length} view{EMBEDDED.views.length === 1 ? '' : 's'} and {EMBEDDED.samples.length} sample{EMBEDDED.samples.length === 1 ? '' : 's'} embedded{EMBEDDED.saved ? ` on ${EMBEDDED.saved.slice(0, 10)}` : ''}: coverage, junctions and retention counts of each view's window are in this file, the alignment files are not needed for them.
          Other genes, the reads track and exon depths need the original BAM/CRAM files (add them with the buttons above); gene lookups, common SNPs and GTEx use the network when it is available.
        </div>
      )}
      {(notes.length > 0 || error || (opened && !samples.length)) && (
        <div className="px-5 py-2 text-xs space-y-0.5">
          {opened && !samples.length && (
            <div className="text-indigo-800">
              {LINK ? `${describeLink(LINK)} · ` : ''}No alignment yet: add the BAM or CRAM files (with their index) with the button above, or drop them on the page. They stay on this computer.
            </div>
          )}
          {notes.map((n, i) => <div key={i} className="text-gray-600">{n}</div>)}
          {error && <div className="text-red-600">{error}</div>}
        </div>
      )}
      {!opened ? (
        <div className="m-6 p-10 border-2 border-dashed border-indigo-300 rounded-2xl bg-white text-center">
          <div className="text-xl font-semibold text-indigo-700">Drop a run folder, BAM / CRAM files, or a session file anywhere on this page</div>
          <div className="text-sm text-gray-600 mt-2 max-w-2xl mx-auto">
            Add each alignment with its index (<code>.bam</code> + <code>.bai</code>, or <code>.cram</code> + <code>.crai</code>). The first file is the primary sample, the others are comparison samples; click a sample chip (or "make primary" on its track) to switch.
            CRAM needs the reference: add an indexed FASTA (<code>.fa</code> + <code>.fai</code>, bgzipped with <code>.gzi</code>) or let the page fetch it from the UCSC API.
            Then type a gene and press Open.
          </div>
          <div className="text-xs text-gray-500 mt-3 max-w-2xl mx-auto">
            Everything stays on this computer: the page lists the files and reads the parts it draws. A browser dialog that speaks of "uploading" a folder is only asking whether this page may read it.
          </div>
          <div className="text-xs text-gray-500 mt-4">
            {samples.length ? `${samples.length} sample${samples.length > 1 ? 's' : ''} ready` : 'No files yet'} · reads are decoded locally with the GMOD BAM/CRAM libraries
          </div>
        </div>
      ) : (
        <div className="p-3">
          <SashimiViewer key={viewerKey} geneName={opened.geneName} geneId={opened.geneId} chrom={opened.chrom} geneStart={opened.start} geneEnd={opened.end}
            sampleId={samples[0]?.id ?? 0} sampleName={samples[0]?.name ?? ''} runId={0} darkMode={false} onClose={() => { if (activeId != null) closeTab(activeId); }} embedded dataSource={ds} allowPrimarySwitch onPrimaryChange={makePrimary} initialView={opened.view} initialMark={opened.mark} initialReads={opened.reads} sampleNames={sampleNames} knownVariantsVersion={knownSeq.current} sampleTypes={sampleTypes} onLibraryEvidence={onLibraryEvidence}
            initialSettings={viewerInit} onStateChange={s => { viewerStateRef.current = s; pendingSettingsRef.current = undefined; }} />
        </div>
      )}
      {exportProgress && (
        <div className="fixed inset-0 z-[60] bg-black/40 flex items-start justify-center p-6" role="status" aria-live="polite" data-export-progress>
          <div className="bg-white rounded-xl shadow-2xl border border-gray-200 w-full max-w-md text-gray-900 text-xs mt-16" onMouseDown={e => e.stopPropagation()}>
            <div className="flex items-center gap-3 px-4 py-3 border-b border-gray-200">
              <span className="inline-block w-5 h-5 rounded-full border-[3px] border-indigo-500 border-t-transparent animate-spin shrink-0" />
              <div>
                <div className="font-bold text-sm">{exportProgress.title}</div>
                <div className="text-[11px] text-gray-500 mt-0.5">The data of every registered view is read from the files and written into one page. The download starts when it is finished: keep this tab open.</div>
              </div>
            </div>
            <div className="px-4 py-3 space-y-2">
              <div className="text-gray-700 truncate" title={exportProgress.step}>{exportProgress.step}</div>
              <div className="h-2 rounded-full bg-gray-200 overflow-hidden">
                <div className="h-full bg-indigo-600 transition-[width] duration-300" style={{ width: `${Math.round(100 * Math.min(1, exportProgress.done / Math.max(1, exportProgress.total)))}%` }} />
              </div>
              <div className="text-[11px] text-gray-500">{Math.min(exportProgress.done, exportProgress.total)} of {exportProgress.total} step{exportProgress.total === 1 ? '' : 's'}</div>
            </div>
          </div>
        </div>
      )}
      {exportDialog && (
        <div className="fixed inset-0 z-50 bg-black/40 flex items-start justify-center p-6 overflow-y-auto" onMouseDown={() => setExportDialog(null)}>
          <div className="bg-white rounded-xl shadow-2xl border border-gray-200 w-full max-w-xl text-gray-900 text-xs" onMouseDown={e => e.stopPropagation()}>
            <div className="flex items-start justify-between gap-4 px-4 py-3 border-b border-gray-200">
              <div>
                <div className="font-bold text-sm">Export HTML</div>
                <div className="text-[11px] text-gray-500 mt-0.5">One file, this viewer with the data of the {views.length || 1} registered view{views.length === 1 ? '' : 's'} embedded: gene models, and for every sample the coverage, junctions and retention counts of each view's window with its margins. The recipient opens it in a browser without the alignment files.</div>
              </div>
              <button onClick={() => setExportDialog(null)} className="text-gray-400 hover:text-gray-700 text-lg leading-none px-1" title="Close">×</button>
            </div>
            <div className="px-4 py-3 space-y-3">
              <fieldset className="space-y-1">
                <legend className="font-medium mb-1">Window exported around each view · coverage, junctions, retention counts and reads</legend>
                {([['view', 'the view as shown: the recipient cannot pan outside it'], ['margin', 'the view with a margin of half its width on each side'], ['max', `the widest window: what the viewer itself loads around the view (up to ${MAX_FETCH_BP / 1e6} Mb for coverage and junctions, ${READS_MAX_VIEW_BP / 1000} kb for reads)`]] as const).map(([v, label]) => (
                  <label key={v} className="flex items-center gap-2"><input type="radio" name="exportWindow" checked={exportDialog.window === v} onChange={() => setExportDialog({ ...exportDialog, window: v })} />{label}</label>
                ))}
              </fieldset>
              <div className="font-semibold">Reads track · {viewsWithReads ? `${viewsWithReads} view${viewsWithReads === 1 ? '' : 's'} with the reads track on` : 'no view has the reads track on'}</div>
              <div className="text-[11px] text-gray-500">For those views the reads of every loaded sample are embedded over the same window, with the reference bases and the mismatches (the recipient can switch reads / collapsed and change Min VAF); read names are replaced by numbers. Views above {READS_MAX_VIEW_BP / 1000} kb have no reads track and are skipped.</div>
              <fieldset className="space-y-1">
                <legend className="font-medium mb-1">Reads per sample in that window</legend>
                {([['shown', 'as displayed: up to 20,000 reads (sampled evenly when the window holds more) · about 150 kB per sample and view'], ['dense', 'dense: up to 100,000 reads, for deep windows · under 1 MB per sample and view'], ['all', 'every read of the window · exact at any zoom, a few MB per sample and view on a deep window']] as const).map(([v, label]) => (
                  <label key={v} className="flex items-center gap-2"><input type="radio" name="readsCap" checked={exportDialog.readsCap === v} onChange={() => setExportDialog({ ...exportDialog, readsCap: v })} disabled={!viewsWithReads} />{label}</label>
                ))}
              </fieldset>
            </div>
            <div className="flex items-center gap-2 px-4 py-3 border-t border-gray-200">
              <span className="text-[11px] text-gray-500">File name: <code>{(sessionName.trim() || defaultSessionName(opened?.geneName)).replace(/\.json$/i, '').replace(/\.html$/i, '')}.html</code></span>
              <span className="ml-auto flex gap-2">
                <button onClick={() => setExportDialog(null)} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-gray-50">Cancel</button>
                <button onClick={() => { const o = exportDialog; setExportDialog(null); void exportHtml(o); }} className="px-3 py-1 text-xs rounded bg-emerald-600 text-white hover:bg-emerald-700 font-medium">Export</button>
              </span>
            </div>
          </div>
        </div>
      )}
      <footer className="px-5 py-3 text-[11px] text-gray-500 flex flex-wrap gap-x-3 gap-y-1">
        <span>Benjamin Cogné (CHU Nantes, 2026)</span>
        <span>·</span>
        <span>made with Claude Opus 5 and Fable 5.1</span>
        <span>·</span>
        <a href="https://github.com/benjamin-cogne/Sashimi-viewer" target="_blank" rel="noopener noreferrer" className="hover:text-indigo-700 underline decoration-dotted">github.com/benjamin-cogne/Sashimi-viewer</a>
        <span>·</span>
        <span>CC BY-NC 4.0</span>
        {DEV && <span className="text-red-700 font-semibold">DEV build {DEV.branch}{DEV.sha ? ` @ ${DEV.sha}` : ''}</span>}
      </footer>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById('root')!).render(<App />);
