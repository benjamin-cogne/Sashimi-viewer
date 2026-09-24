/**
 * Session files: everything needed to come back to a view later, as JSON.
 *
 * A browser cannot read or store file paths, so the alignments are recorded by file name (with
 * their index name and size) and the user adds the files again when loading; the session then
 * matches them by name, restores the sample names, the order (first = primary), the gene and
 * window, and every viewer option (groups and reads track refer to samples by name).
 */
import type { ViewerSettings, ViewerState } from '../components/SashimiViewer';
import type { GenomeBuild } from './ensembl';
import type { LibraryEvidence } from '../components/sashimi/types';

export const SESSION_APP = 'sashimi-viewer';
export const SESSION_VERSION = 1;

export interface SessionSample { name: string; file: string; index: string; size: number; kind: 'bam' | 'cram'; /** paths relative to the run folder, when known */ path?: string; indexPath?: string; /** RNA-seq or DNA, as decided or chosen in the saving page */ library?: LibraryEvidence }
export interface SessionGene {
  name: string; id?: string; chrom: string;
  /** gene bounds, 1-based inclusive */
  start: number; end: number;
  /** window shown, 1-based inclusive */
  view: { start: number; end: number };
  mark: { start: number; end: number } | null;
}
/** Viewer options with samples referred to by name (ids are per page load). */
export type SessionViewer = Omit<ViewerSettings, 'groups' | 'readsSample'> & { groups: { name: string; samples: string[]; color?: string }[]; readsSample: string | null };
export interface SessionFile {
  app: typeof SESSION_APP;
  version: number;
  saved: string;
  build: GenomeBuild;
  /** the run folder the sample paths are relative to (its name only; Chromium remembers the folder itself) */
  folder: string | null;
  /** in display order: the first sample is the primary one */
  samples: SessionSample[];
  fasta: { file: string; index: string; gzi?: string } | null;
  gene: SessionGene | null;
  viewer: SessionViewer | null;
  /** every registered view (tab), in order; `gene` / `viewer` above describe the active one (older readers use them) */
  views?: SessionView[];
  activeView?: number;
  /** variants of interest entered in the header (notation as typed, label shown), drawn on every view */
  knownVariants?: { text: string; label: string }[];
}
/** One registered view: the region and the full set of options it was left with. */
export interface SessionView { label: string; gene: SessionGene; viewer: SessionViewer | null }

/** The minimum a loaded sample must expose (LocalSample without the File objects' behaviour). */
export interface SampleLike { id: number; name: string; kind: 'bam' | 'cram'; file: { name: string; size: number }; index: { name: string }; path?: string; indexPath?: string; lib?: LibraryEvidence }

export function defaultSessionName(geneName?: string | null): string {
  const day = new Date().toISOString().slice(0, 10);
  return `sashimi-session${geneName ? `-${geneName.replace(/[^\w.-]+/g, '_')}` : ''}-${day}.json`;
}

export function buildSession(args: {
  build: GenomeBuild;
  folder?: string | null;
  samples: SampleLike[];
  fasta?: { fa: { name: string }; fai: { name: string }; gzi?: { name: string } };
  /** the viewer's latest state (gene, window, options); null when no gene is open */
  state: ViewerState | null;
  /** every registered view with its own state (the active one included) */
  views?: { label: string; state: ViewerState }[];
  activeView?: number;
  knownVariants?: { text: string; label: string }[];
}): SessionFile {
  const { build, samples, fasta, state } = args;
  const folder = args.folder ?? null;
  const nameOf = (id: number | null | undefined) => samples.find(s => s.id === id)?.name ?? null;
  const viewerOf = (state: ViewerState): SessionViewer => ({
    equalIntrons: state.equalIntrons, intronWidth: state.intronWidth ?? null, allTranscripts: state.allTranscripts, commonSnps: state.commonSnps, snpMinAf: state.snpMinAf,
    depthAxis: state.depthAxis, uniqueOnly: state.uniqueOnly,
    reads: state.reads, readsAll: state.readsAll, readsSample: nameOf(state.readsSample), collapseReads: state.collapseReads, minVafPct: state.minVafPct,
    minJunctionReads: state.minJunctionReads, minJunctionReadsSet: state.minJunctionReadsSet, minUsagePct: state.minUsagePct, arcLabels: state.arcLabels, intronRetention: state.intronRetention,
    viewMode: state.viewMode,
    groups: state.groups.map(g => ({ name: g.name, samples: g.sampleIds.map(nameOf).filter((n): n is string => !!n), color: g.color })),
    knownVariants: state.knownVariants, hiddenJunctions: state.hiddenJunctions ?? [], labelScales: state.labelScales && Object.keys(state.labelScales).length ? state.labelScales : undefined, hiddenTranscripts: state.hiddenTranscripts?.length ? state.hiddenTranscripts : undefined, transcriptId: state.transcriptId,
    consensusMode: state.consensusMode, minIndelBp: state.minIndelBp, longReadMinVafPct: state.longReadMinVafPct, coverageVariants: state.coverageVariants, methylation: state.methylation, methylIslands: state.methylIslands, pairs: state.pairs, haplotypes: state.haplotypes, phaseSource: state.phaseSource, readsGroup: state.readsGroup, clippedBases: state.clippedBases, insertedBases: state.insertedBases,
  });
  const geneOf = (state: ViewerState): SessionGene => ({ name: state.gene.name, id: state.gene.id, chrom: state.gene.chrom, start: state.gene.start, end: state.gene.end, view: { ...state.view }, mark: state.mark });
  const views = args.views?.map(v => ({ label: v.label, gene: geneOf(v.state), viewer: viewerOf(v.state) }));
  return {
    app: SESSION_APP, version: SESSION_VERSION, saved: new Date().toISOString(), build, folder,
    samples: samples.map(s => ({ name: s.name, file: s.file.name, index: s.index.name, size: s.file.size, kind: s.kind, path: s.path, indexPath: s.indexPath, library: s.lib && s.lib.type !== 'unknown' ? s.lib : undefined })),
    fasta: fasta ? { file: fasta.fa.name, index: fasta.fai.name, gzi: fasta.gzi?.name } : null,
    gene: state ? geneOf(state) : null,
    viewer: state ? viewerOf(state) : null,
    views: views && views.length ? views : undefined,
    activeView: views && views.length ? Math.min(Math.max(0, args.activeView ?? 0), views.length - 1) : undefined,
    knownVariants: args.knownVariants?.length ? args.knownVariants.map(k => ({ text: k.text, label: k.label })) : undefined,
  };
}

/** Parses and checks a session file; throws a readable error for anything else. */
export function parseSession(text: string): SessionFile {
  let raw: any;
  try { raw = JSON.parse(text); } catch { throw new Error('not a JSON file'); }
  if (!raw || raw.app !== SESSION_APP) throw new Error('not a Sashimi viewer session file (missing "app": "sashimi-viewer")');
  if (typeof raw.version !== 'number' || raw.version > SESSION_VERSION) throw new Error(`session version ${raw.version} is newer than this viewer (${SESSION_VERSION})`);
  if (!Array.isArray(raw.samples)) throw new Error('session has no "samples" list');
  const build: GenomeBuild = raw.build === 'GRCh37' ? 'GRCh37' : 'GRCh38';
  const samples: SessionSample[] = raw.samples
    .filter((s: any) => s && typeof s.file === 'string')
    .map((s: any) => ({ name: String(s.name || s.file), file: String(s.file), index: String(s.index || ''), size: Number(s.size) || 0, kind: s.kind === 'cram' ? 'cram' : 'bam',
      path: typeof s.path === 'string' ? s.path : undefined, indexPath: typeof s.indexPath === 'string' ? s.indexPath : undefined,
      library: s.library && (s.library.type === 'rna' || s.library.type === 'dna') ? { type: s.library.type, source: ['header', 'reads', 'user'].includes(s.library.source) ? s.library.source : 'user', note: String(s.library.note ?? '') } : undefined }));
  const parseGene = (g: any): SessionGene | null => g && typeof g.name === 'string' && typeof g.chrom === 'string' && Number.isFinite(g.start) && Number.isFinite(g.end)
    ? { name: g.name, id: typeof g.id === 'string' ? g.id : undefined, chrom: g.chrom, start: g.start, end: g.end,
        view: g.view && Number.isFinite(g.view.start) && Number.isFinite(g.view.end) ? { start: g.view.start, end: g.view.end } : { start: g.start, end: g.end },
        mark: g.mark && Number.isFinite(g.mark.start) && Number.isFinite(g.mark.end) ? { start: g.mark.start, end: g.mark.end } : null }
    : null;
  const parseViewer = (v: any): SessionViewer | null => v && typeof v === 'object' ? {
    ...v,
    groups: Array.isArray(v.groups) ? v.groups.filter((x: any) => x && typeof x.name === 'string').map((x: any) => ({ name: x.name, samples: Array.isArray(x.samples) ? x.samples.map(String) : [], color: typeof x.color === 'string' && /^#[0-9a-f]{6}$/i.test(x.color) ? x.color : undefined })) : [],
    readsSample: typeof v.readsSample === 'string' ? v.readsSample : null,
    hiddenJunctions: Array.isArray(v.hiddenJunctions) ? v.hiddenJunctions.filter((x: any) => typeof x === 'string') : [],
    labelScales: v.labelScales && typeof v.labelScales === 'object' ? Object.fromEntries(Object.entries(v.labelScales).filter(([, n]) => typeof n === 'number' && Number.isFinite(n) && n > 0)) as Record<string, number> : undefined,
    hiddenTranscripts: Array.isArray(v.hiddenTranscripts) ? v.hiddenTranscripts.filter((x: any) => typeof x === 'string') : undefined,
  } : null;
  const gene = parseGene(raw.gene);
  const viewer = parseViewer(raw.viewer);
  const views: SessionView[] | undefined = Array.isArray(raw.views)
    ? raw.views.map((x: any) => { const g = parseGene(x?.gene); return g ? { label: typeof x.label === 'string' && x.label ? x.label : g.name, gene: g, viewer: parseViewer(x.viewer) } : null; }).filter((x: any): x is SessionView => !!x)
    : undefined;
  const activeView = views?.length ? Math.min(Math.max(0, Number(raw.activeView) || 0), views.length - 1) : undefined;
  const knownVariants = Array.isArray(raw.knownVariants)
    ? raw.knownVariants.filter((k: any) => k && typeof k.text === 'string' && k.text.trim()).map((k: any) => ({ text: String(k.text), label: typeof k.label === 'string' ? k.label : '' }))
    : undefined;
  return { app: SESSION_APP, version: raw.version, saved: String(raw.saved || ''), build, folder: typeof raw.folder === 'string' ? raw.folder : null, samples, fasta: raw.fasta && typeof raw.fasta.file === 'string' ? raw.fasta : null, gene, viewer, views: views?.length ? views : undefined, activeView, knownVariants: knownVariants?.length ? knownVariants : undefined };
}

/** Which session samples are present among the loaded files (by relative path, else by file name, else by index name), and which are still missing. */
export function matchSession(session: SessionFile, samples: SampleLike[]): { matched: { entry: SessionSample; sample: SampleLike }[]; missing: SessionSample[] } {
  const used = new Set<number>();
  const matched: { entry: SessionSample; sample: SampleLike }[] = [];
  const missing: SessionSample[] = [];
  for (const entry of session.samples) {
    const hit = (entry.path ? samples.find(s => !used.has(s.id) && s.path === entry.path) : undefined)
      ?? samples.find(s => !used.has(s.id) && (s.file.name === entry.file || (!!entry.index && s.index.name === entry.index)));
    if (hit) { used.add(hit.id); matched.push({ entry, sample: hit }); } else missing.push(entry);
  }
  return { matched, missing };
}

/** The viewer options of a session with sample names turned into the ids of the loaded samples. */
export function viewerSettingsOf(session: SessionFile | { viewer: SessionViewer | null }, samples: SampleLike[]): Partial<ViewerSettings> | undefined {
  const v = session.viewer;
  if (!v) return undefined;
  const idOf = (name: string) => samples.find(s => s.name === name)?.id;
  const { groups, readsSample, ...rest } = v;
  return {
    ...rest,
    readsSample: readsSample ? idOf(readsSample) ?? null : null,
    groups: groups.map(g => ({ name: g.name, sampleIds: g.samples.map(idOf).filter((id): id is number => id != null), color: g.color })),
  };
}
