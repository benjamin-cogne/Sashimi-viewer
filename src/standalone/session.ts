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

export const SESSION_APP = 'sashimi-viewer';
export const SESSION_VERSION = 1;

export interface SessionSample { name: string; file: string; index: string; size: number; kind: 'bam' | 'cram'; /** paths relative to the run folder, when known */ path?: string; indexPath?: string }
export interface SessionGene {
  name: string; chrom: string;
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
}

/** The minimum a loaded sample must expose (LocalSample without the File objects' behaviour). */
export interface SampleLike { id: number; name: string; kind: 'bam' | 'cram'; file: { name: string; size: number }; index: { name: string }; path?: string; indexPath?: string }

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
}): SessionFile {
  const { build, samples, fasta, state } = args;
  const folder = args.folder ?? null;
  const nameOf = (id: number | null | undefined) => samples.find(s => s.id === id)?.name ?? null;
  const viewer: SessionViewer | null = state ? {
    equalIntrons: state.equalIntrons, intronWidth: state.intronWidth ?? null, allTranscripts: state.allTranscripts, commonSnps: state.commonSnps, snpMinAf: state.snpMinAf,
    depthAxis: state.depthAxis, uniqueOnly: state.uniqueOnly,
    reads: state.reads, readsAll: state.readsAll, readsSample: nameOf(state.readsSample), collapseReads: state.collapseReads, minVafPct: state.minVafPct,
    minJunctionReads: state.minJunctionReads, minUsagePct: state.minUsagePct, arcLabels: state.arcLabels, intronRetention: state.intronRetention,
    viewMode: state.viewMode,
    groups: state.groups.map(g => ({ name: g.name, samples: g.sampleIds.map(nameOf).filter((n): n is string => !!n), color: g.color })),
    knownVariants: state.knownVariants, hiddenJunctions: state.hiddenJunctions ?? [], transcriptId: state.transcriptId,
  } : null;
  return {
    app: SESSION_APP, version: SESSION_VERSION, saved: new Date().toISOString(), build, folder,
    samples: samples.map(s => ({ name: s.name, file: s.file.name, index: s.index.name, size: s.file.size, kind: s.kind, path: s.path, indexPath: s.indexPath })),
    fasta: fasta ? { file: fasta.fa.name, index: fasta.fai.name, gzi: fasta.gzi?.name } : null,
    gene: state ? { name: state.gene.name, chrom: state.gene.chrom, start: state.gene.start, end: state.gene.end, view: { ...state.view }, mark: state.mark } : null,
    viewer,
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
      path: typeof s.path === 'string' ? s.path : undefined, indexPath: typeof s.indexPath === 'string' ? s.indexPath : undefined }));
  const g = raw.gene;
  const gene: SessionGene | null = g && typeof g.name === 'string' && typeof g.chrom === 'string' && Number.isFinite(g.start) && Number.isFinite(g.end)
    ? { name: g.name, chrom: g.chrom, start: g.start, end: g.end,
        view: g.view && Number.isFinite(g.view.start) && Number.isFinite(g.view.end) ? { start: g.view.start, end: g.view.end } : { start: g.start, end: g.end },
        mark: g.mark && Number.isFinite(g.mark.start) && Number.isFinite(g.mark.end) ? { start: g.mark.start, end: g.mark.end } : null }
    : null;
  const v = raw.viewer && typeof raw.viewer === 'object' ? raw.viewer : null;
  const viewer: SessionViewer | null = v ? {
    ...v,
    groups: Array.isArray(v.groups) ? v.groups.filter((x: any) => x && typeof x.name === 'string').map((x: any) => ({ name: x.name, samples: Array.isArray(x.samples) ? x.samples.map(String) : [], color: typeof x.color === 'string' && /^#[0-9a-f]{6}$/i.test(x.color) ? x.color : undefined })) : [],
    readsSample: typeof v.readsSample === 'string' ? v.readsSample : null,
    hiddenJunctions: Array.isArray(v.hiddenJunctions) ? v.hiddenJunctions.filter((x: any) => typeof x === 'string') : [],
  } : null;
  return { app: SESSION_APP, version: raw.version, saved: String(raw.saved || ''), build, folder: typeof raw.folder === 'string' ? raw.folder : null, samples, fasta: raw.fasta && typeof raw.fasta.file === 'string' ? raw.fasta : null, gene, viewer };
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
export function viewerSettingsOf(session: SessionFile, samples: SampleLike[]): Partial<ViewerSettings> | undefined {
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
