import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { SashimiDataSource } from './sashimi/datasource';
import type { TranscriptData, CoverageRun, JunctionArc, BoundarySpanning, BoundaryHint, ReadsResponse, AlignedRead, ReadGroup, VariantSite, AllTranscripts, TranscriptModel, GeneModel, ExonUsageResponse, CommonSnp, GtexTissue, KnownVariant, RegionHint } from './sashimi/types';
import {
  LINEAR_AXIS, equalIntronAxis, defaultIntronV, makeScale, toTxModel, intronsOf,
  buildCoveragePaths, depthAt, maxDepthIn,
  classifyJunction, layerJunctions, junctionKey, arcGeom, arcYAtX,
  niceTicks, niceMax, formatBp, packReads, cdnaPosition, junctionHgvs, exonPsi, junctionAlternative, junctionFrame, codonsInWindow, parseLocus,
  usageIntervals, referenceExons, exonUsage, usageCohort, usageZ, exonSiteUsage,
  type TxModel, type Scale, type VirtualAxis, type JunctionInfo, type FrameInfo,
} from './sashimi/geometry';
import SpliceCartoon from './sashimi/SpliceCartoon';
import { spliceEvent, spliceStory, storyWindows, type SpliceStory } from './sashimi/spliceModel';
import { SNP_MAX_WINDOW, snpSourceLabel } from '../standalone/snps';
import { SPAN_EXON_ANCHOR, SPAN_INTRON_ANCHOR } from '../standalone/alignments';
import { KNOWN_VARIANT_COLORS, KNOWN_VARIANT_KIND_NAMES, isPointVariant, knownVariantTitle } from './sashimi/knownVariants';
import { GTEX_DEFAULT_FAVOURITES } from '../standalone/gtex';
import { sumCoverage, poolJunctions, poolSpanning, aggregateJunctions, pctLabel, AGG_CLASS_LABEL, PSEUDO_EXON_MAX_BP, type AggEvent, type AggResult } from './sashimi/aggregate';

// ======================== Types ========================

interface SashimiViewerProps {
  geneName: string;
  geneId?: string;
  chrom: string;
  /** 1-based inclusive gene bounds, as shown in the OUTRIDER / FRASER tables. */
  geneStart: number;
  geneEnd: number;
  /** Primary sample; 0 (or negative) means no alignment yet: the annotation, reference, known variants and GTEx still show, coverage waits for samples. */
  sampleId: number;
  sampleName: string;
  runId: number;
  darkMode: boolean;
  onClose: () => void;
  embedded?: boolean;
  /** Called with a PNG of the current plot and the viewer context when the user clicks the basket camera. */
  onSnapshot?: (png: Blob, context: SashimiSnapshotContext) => Promise<void>;
  /** Where the data comes from: the application backend (apiDataSource) or local files (standalone). */
  dataSource: SashimiDataSource;
  /** Hide the "Add sample" picker (standalone viewer manages its own files). */
  hideSamplePicker?: boolean;
  /** Let the user promote any loaded sample to primary (first track, reference of the "unique junction" comparison). */
  allowPrimarySwitch?: boolean;
  /** Called when the user promotes a sample from inside the viewer, so the host can mirror it. */
  onPrimaryChange?: (sampleId: number) => void;
  /** Window to show at first instead of the whole gene (1-based inclusive), e.g. the locus the host was asked for. */
  initialView?: { start: number; end: number };
  /** Locus label to pin (1-based inclusive) when it is narrower than the initial window, e.g. the variant a deep link asked for. */
  initialMark?: { start: number; end: number };
  /** Open with the reads track on (deep links at base resolution). */
  initialReads?: boolean;
  /** Display names chosen by the host (renamed samples), by sample id; tracks follow without remounting. */
  sampleNames?: Record<number, string>;
  /** Options to start with (a saved session, or the previous viewer's options when the host remounts it). */
  initialSettings?: Partial<ViewerSettings>;
  /** Called whenever an option or the navigation changes, with everything a session file needs. */
  onStateChange?: (state: ViewerState) => void;
}

/** Every user option of the viewer, as stored in a session file. Samples are referred to by id (the host maps names ↔ ids). */
export interface ViewerSettings {
  equalIntrons: boolean; allTranscripts: boolean; commonSnps: boolean; snpMinAf: number;
  /** width of every intron in equal-introns mode, bp-equivalents; null = default (median exon length, 80–300) */
  intronWidth: number | null;
  depthAxis: DepthAxis; uniqueOnly: boolean;
  reads: boolean; readsAll: boolean; readsSample: number | null; collapseReads: boolean; minVafPct: number;
  minJunctionReads: number; minUsagePct: number; arcLabels: 'reads' | 'usage'; intronRetention: boolean;
  viewMode: 'samples' | 'groups'; groups: { name: string; sampleIds: number[] }[];
  knownVariants: boolean;
  /** reference transcript chosen in the transcript list; absent = the default model of the gene */
  transcriptId?: string;
}
/** The options plus where the viewer is: gene, window and pinned locus, 1-based inclusive. */
export interface ViewerState extends ViewerSettings {
  gene: { name: string; id?: string; chrom: string; start: number; end: number };
  view: { chrom: string; start: number; end: number };
  mark: { start: number; end: number } | null;
}

/** What a basket screenshot documents: the region, the samples and every option in effect. */
export type DepthAxis = 'shared' | 'own' | 'relative';

export interface SashimiSnapshotContext {
  viewer: 'sashimi';
  gene: string;
  transcript?: string;
  region: { chrom: string; start: number; end: number };
  samples: string[];
  primarySample: string;
  options: {
    equalIntrons: boolean; allTranscripts: boolean;
    /** shared: one depth axis for every sample; own: each sample scaled to its own maximum; relative: each sample as % of its own maximum */
    depthAxis: DepthAxis;
    /** kept for consumers of older snapshots: depthAxis === 'shared' */
    sharedY: boolean; uniqueOnly: boolean;
    reads: boolean; readsSample?: string; collapsed: boolean; minJunctionReads: number; minVafPct: number;
    /** minimum allele frequency of the common-SNP track, null when the track is off */
    commonSnpsMinAf?: number | null;
    /** aggregate view: one pooled track per sample group instead of one track per sample */
    aggregate?: boolean;
    /** arc labels of the sample tracks: spliced reads or % usage */
    arcLabels?: 'reads' | 'usage';
    /** intron retention counted in the usage percentages */
    intronRetention?: boolean;
  };
  /** sample groups defined for the aggregate view */
  groups?: { name: string; samples: string[] }[];
  variantSites?: { pos: number; ref: string; alt: string; vaf: number; depth: number }[];
  /** variants previously identified in the primary sample that lie in the region (1-based start) */
  knownVariants?: { kind: string; chrom: string; start: number; end: number; label: string; text: string }[];
  timestamp: string;
}

/** Rasterise an SVG element to a PNG blob (white background, `scale`× the CSS size). */
async function svgToPng(svg: SVGSVGElement, scale = 2): Promise<Blob> {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.querySelectorAll('[data-export="skip"]').forEach(n => n.remove());
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  const w = svg.width.baseVal.value || svg.clientWidth, h = svg.height.baseVal.value || svg.clientHeight;
  clone.setAttribute('width', String(w));
  clone.setAttribute('height', String(h));
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' }));
  try {
    const img = new Image();
    await new Promise<void>((resolve, reject) => { img.onload = () => resolve(); img.onerror = () => reject(new Error('SVG rasterisation failed')); img.src = url; });
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(w * scale); canvas.height = Math.round(h * scale);
    const ctx = canvas.getContext('2d')!;
    ctx.scale(scale, scale);
    ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);
    return await new Promise<Blob>((resolve, reject) => canvas.toBlob(b => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png'));
  } finally {
    URL.revokeObjectURL(url);
  }
}

interface FetchWindow { chrom: string; start: number; end: number; uniqueOnly: boolean }

interface TrackData {
  sampleId: number;
  sampleName: string;
  coverage: CoverageRun[];
  junctions: JunctionArc[];
  loading: boolean;
  error?: string;
  /** Window the current coverage/junctions were fetched for (with margin). */
  fetched?: FetchWindow;
  /** Unspliced reads through the exon–intron boundaries (intron retention), when the source counts them. */
  spanning?: BoundarySpanning;
  /** The source decoded one read in `rate` of this window: depths and counts are scaled estimates. */
  sampled?: { rate: number; total: number; decoded: number };
  /** GTEx tissue track (median junction reads + reads-per-base exon profile); sampleId is negative */
  gtex?: { tissue: GtexTissue; dataset: string; unit: string; warning?: string; tpm: number | null; lowCoverage: boolean };
  /** Pooled track of a sample group (aggregate view); sampleId is negative */
  group?: { id: number; n: number; loaded: number; agg: AggResult; samplesWith: Map<string, number> };
}

/** A named set of samples pooled into one track in the aggregate ("Groups") view. */
interface SampleGroup { id: number; name: string; sampleIds: number[] }
const GROUP_ID_BASE = -100000;   // group tracks use sampleId = GROUP_ID_BASE - group id (negative, like GTEx tracks)
const PSEUDO_EXON_COLOR = '#7c3aed';
const RETENTION_COLOR = '#0d9488';
/** Exon–intron boundaries of a model, for the unspliced-read counts of the coverage request. */
const boundariesOf = (t: TxModel | null): BoundaryHint | undefined =>
  t && t.exons.length > 1 ? { intronStarts: t.exons.slice(0, -1).map(e => e.end), intronEnds: t.exons.slice(1).map(e => e.start) } : undefined;
/** Same rule as the decoder's boundary counts: one aligned block through an annotated boundary with the exon and intron anchors. */
const readSpansBoundary = (r: AlignedRead, b: BoundaryHint | undefined): boolean =>
  !!b && r.b.some(([bs, be]) =>
    b.intronStarts.some(p => bs <= p - SPAN_EXON_ANCHOR && be >= p + SPAN_INTRON_ANCHOR) ||
    b.intronEnds.some(q => bs <= q - SPAN_INTRON_ANCHOR && be >= q + SPAN_EXON_ANCHOR));

const GTEX_FAV_KEY = 'sashimi.gtex.favourites';
function loadGtexFavourites(): string[] {
  try { const v = JSON.parse(localStorage.getItem(GTEX_FAV_KEY) || 'null'); if (Array.isArray(v) && v.length) return v.map(String); } catch { /* ignore */ }
  return GTEX_DEFAULT_FAVOURITES;
}
function saveGtexFavourites(ids: string[]) { try { localStorage.setItem(GTEX_FAV_KEY, JSON.stringify(ids.slice(0, 8))); } catch { /* ignore */ } }

// ======================== Layout constants ========================

const PLOT_LEFT = 64;       // room for the depth axis
const PLOT_RIGHT_PAD = 36;  // room for the per-track remove button
const RULER_H = 42;
const COVERAGE_H = 130;
const TRACK_LABEL_H = 20;   // band at the top of each track reserved for the sample label (arcs and coverage stay below it)
const JUNC_LEVEL_STEP = 17; // extra apex height per arc nesting level
const JUNC_MIN_H = 6;       // smallest junction area; it otherwise grows to what the arcs and pills really occupy
const JUNC_PAD = 5;         // clearance between the highest arc or pill and the sample-label band
const LABEL_H = 15;         // read-count pill height
const TRACK_GAP = 10;       // gap between panels (transcript, variants, reads)
const SASHIMI_GAP = 4;      // gap between consecutive sample tracks, kept small so the samples read as one group
const TRANSCRIPT_H = 78;
const NEIGHBOUR_ROW_H = 22;    // one row per neighbouring gene drawn under the queried gene
const ALT_TX_ROW_H = 18;       // one row per transcript model in the "All transcripts" panel
const ALT_TX_HEADER_H = 22;
const ALT_TX_MAX_ROWS = 40;
const GTEX_MIN_TPM = 1;      // below this median TPM a GTEx tissue track only says "low coverage"
const SNP_PANEL_H = 46;        // common-SNP track: header + lollipops
const KNOWN_ROW_H = 15;        // known-variant panel: one row per stacked variant; the first row shares the line with the panel title
const KNOWN_PAD = 3;           // padding above the first and below the last row of the known-variant panel
const SNP_MAX_MARKS = 4000;    // beyond this many variants in view, ticks only
const LEGEND_ROW_H = 22;
const MIN_V_SPAN = 40;            // smallest zoom window, in virtual (bp-equivalent) units
const MAX_VIEW_BP = 4_000_000;    // largest zoom-out window
const MAX_FETCH_BP = 2_000_000;   // largest window fetched at once (view + margins)
const MAX_READS_PER_TRACK = 250_000; // reads decoded per coverage request; the source shrinks the margins and then samples 1 in 2, 4, 8… past it

// Reads track (IGV-like alignment view)
const READS_MAX_VIEW_BP = 100_000; // reads load only below this window size (IGV's "visibility window")
const READS_MAX_ROWS = 120;
const READS_HEADER_H = 22;
const READS_SEQ_ROW_H = 18;
const READS_MAX = 2500;
const SITES_STRIP_H = 18;      // strip above the sashimi holding the variant-site stars
const GROUP_ROW_H = 16;        // consensus row height in collapsed mode
const AA_ROW_H = 16;           // amino-acid row above the reference bases

// ======================== Colours ========================
// Categorical slots validated for colour-vision deficiency (adjacent-pair ΔE ≥ 8);
// red is reserved for "junction unique to the primary sample" and never used for a track.

const TRACK_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7'];
const UNIQUE_COLOR = '#e34948';
const INK = {
  bg: '#ffffff', text: '#1f2937', muted: '#6b7280', faint: '#9ca3af',
  grid: '#e5e7eb', gridStrong: '#d1d5db', exon: '#334155', utr: '#94a3b8', intron: '#94a3b8',
  geneBand: 'rgba(99, 102, 241, 0.035)', select: '#6366f1',
};
const FONT = 'Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
// Nucleotide colours follow the IGV convention (A green, C blue, G orange, T red)
const BASE_COLORS: Record<string, string> = { A: '#1a9e37', C: '#2452d6', G: '#d9861c', T: '#d6332b', N: '#6b7280' };
const COMPLEMENT: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const READ_FILL = '#c8cdd6';
const INSERTION_COLOR = '#7c3aed';
const STAR_COLOR = '#f59e0b';
const SNP_SNV_COLOR = '#2563eb';
const SNP_INDEL_COLOR = '#b45309';
const SNP_KNOWN_RING = '#2563eb';
// Neighbouring genes: same transcription direction as the queried gene vs antisense
const SAME_SENSE_COLOR = '#475569';
const ANTISENSE_COLOR = '#7c3aed';
// Reading-frame glyphs next to the read counts of non-canonical junctions
const FRAME_IN_COLOR = '#16a34a';
const FRAME_OUT_COLOR = '#dc2626';
const FRAME_UTR_COLOR = '#9ca3af';
const FRAME_GLYPH_R = 6.5;
// Amino acids: alternating codon fills, start and stop highlighted
const AA_FILLS = ['#cbd5e1', '#94a3b8'];
const AA_NAMES: Record<string, string> = {
  A: 'Alanine', R: 'Arginine', N: 'Asparagine', D: 'Aspartate', C: 'Cysteine', Q: 'Glutamine', E: 'Glutamate', G: 'Glycine',
  H: 'Histidine', I: 'Isoleucine', L: 'Leucine', K: 'Lysine', M: 'Methionine', F: 'Phenylalanine', P: 'Proline', S: 'Serine',
  T: 'Threonine', W: 'Tryptophan', Y: 'Tyrosine', V: 'Valine', X: 'unknown',
};
const AA_START_COLOR = '#16a34a';
const AA_STOP_COLOR = '#dc2626';

/** Five-point star path centred on (cx, cy). */
const starPath = (cx: number, cy: number, R: number): string => {
  const r = R * 0.45;
  let d = '';
  for (let i = 0; i < 10; i++) {
    const rad = i % 2 ? r : R, a = -Math.PI / 2 + (i * Math.PI) / 5;
    d += `${i ? 'L' : 'M'}${(cx + rad * Math.cos(a)).toFixed(1)},${(cy + rad * Math.sin(a)).toFixed(1)}`;
  }
  return d + 'Z';
};

const withAlpha = (hex: string, alpha: number) => {
  const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
};

const isEnsemblId = (id: string) => /^ENST/i.test(id);
/** How the displayed model was chosen, worded for its source (RefSeq via UCSC, or Ensembl as fallback). */
function modelKindLabel(m: { transcriptId: string; modelKind: string }): string {
  if (m.modelKind === 'mane') return 'MANE Select';
  if (m.modelKind === 'chosen') return 'chosen in the transcript list';
  if (m.modelKind === 'canonical') return isEnsemblId(m.transcriptId) ? 'Ensembl canonical (no MANE Select)' : 'RefSeq Select (no MANE Select)';
  return isEnsemblId(m.transcriptId) ? 'longest CDS (no MANE Select, no canonical flag)' : 'longest CDS (no MANE Select, no RefSeq Select)';
}

/** A neighbouring gene (1-based, from the data source) as a 0-based transcript model; exons ranked in transcription order. */
function toNeighbourModel(g: GeneModel, chrom: string): TxModel & { geneId: string; biotype: string; isCanonical: boolean } {
  const sorted = [...g.exons].map(e => ({ start: e.start - 1, end: e.end })).sort((a, b) => a.start - b.start);
  const n = sorted.length;
  const exons = sorted.map((e, i) => ({ ...e, rank: g.strand >= 0 ? i + 1 : n - i }));
  const hasCds = g.cds_start != null && g.cds_end != null;
  return {
    geneName: g.gene_name, transcriptId: g.transcript_id, isMane: false, modelKind: g.is_canonical ? 'canonical' : 'longest', chrom,
    strand: g.strand < 0 ? -1 : 1, start: g.start - 1, end: g.end, exons,
    cdsStart: hasCds ? (g.cds_start as number) - 1 : null, cdsEnd: hasCds ? (g.cds_end as number) : null,
    geneId: g.gene_id, biotype: g.biotype, isCanonical: g.is_canonical,
  };
}

const frameLabel = (f: FrameInfo) => f.frame === 'in' ? 'in frame' : f.frame === 'out' ? 'out of frame' : f.frame === 'utr' ? 'UTR only' : 'unknown';

/**
 * Reading-frame glyph: a red "no entry" disc (white bar) for a frameshift, a green disc with "=" when the
 * frame is kept, a grey "=" when the change stays in the UTR. Nothing is drawn when the frame cannot be judged.
 */
function renderFrameGlyph(cx: number, cy: number, f: FrameInfo, key: string): JSX.Element | null {
  const r = FRAME_GLYPH_R;
  if (f.frame === 'out') {
    return (
      <g key={key}>
        <circle cx={cx} cy={cy} r={r} fill={FRAME_OUT_COLOR} stroke="#ffffff" strokeWidth={1} />
        <rect x={cx - r * 0.62} y={cy - 1.3} width={r * 1.24} height={2.6} rx={1} fill="#ffffff" />
      </g>
    );
  }
  if (f.frame === 'in' || f.frame === 'utr') {
    const color = f.frame === 'in' ? FRAME_IN_COLOR : FRAME_UTR_COLOR;
    return (
      <g key={key}>
        <circle cx={cx} cy={cy} r={r} fill={color} stroke="#ffffff" strokeWidth={1} />
        <rect x={cx - r * 0.55} y={cy - 2.9} width={r * 1.1} height={1.8} fill="#ffffff" />
        <rect x={cx - r * 0.55} y={cy + 1.1} width={r * 1.1} height={1.8} fill="#ffffff" />
      </g>
    );
  }
  return null;
}

// ======================== Component ========================

export default function SashimiViewer({
  geneName, geneId, chrom, geneStart, geneEnd, sampleId, sampleName, runId, onClose, embedded, onSnapshot,
  dataSource, hideSamplePicker, allowPrimarySwitch, onPrimaryChange, initialView, initialMark, initialReads, sampleNames, initialSettings, onStateChange,
}: SashimiViewerProps) {
  const init = initialSettings ?? {};
  const ds = dataSource;
  const [snapshotState, setSnapshotState] = useState<'idle' | 'busy' | 'done' | 'error'>('idle');
  const svgRef = useRef<SVGSVGElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  // ---- Gene / region state (0-based half-open internally) ----
  const [currentGeneName, setCurrentGeneName] = useState(geneName);
  const [currentGeneId, setCurrentGeneId] = useState(geneId);
  const [currentChrom, setCurrentChrom] = useState(chrom);
  const [currentGeneStart, setCurrentGeneStart] = useState(geneStart - 1);
  const [currentGeneEnd, setCurrentGeneEnd] = useState(geneEnd);

  const linearDefault = (gs: number, ge: number): [number, number] => {
    const p = Math.max(Math.round((ge - gs) * 0.1), 500);
    return [Math.max(0, gs - p), ge + p];
  };
  const [viewStart, setViewStart] = useState(() => (initialView ? Math.max(0, initialView.start - 1) : linearDefault(geneStart - 1, geneEnd)[0]));
  const [viewEnd, setViewEnd] = useState(() => (initialView ? Math.max(initialView.end, initialView.start) : linearDefault(geneStart - 1, geneEnd)[1]));
  /** Locus the user asked for by coordinates, drawn as a band / line across the plot until the next gene search. */
  const [locusMark, setLocusMark] = useState<{ chrom: string; start: number; end: number } | null>(() => (initialMark ? { chrom, start: initialMark.start - 1, end: initialMark.end } : initialView ? { chrom, start: initialView.start - 1, end: initialView.end } : null));
  const [searchError, setSearchError] = useState<string | null>(null);
  const [svgWidth, setSvgWidth] = useState(1200);

  // ---- Options ----
  const [equalIntrons, setEqualIntrons] = useState(init.equalIntrons ?? false);
  const [intronWidth, setIntronWidth] = useState<number | null>(init.intronWidth ?? null);
  const [depthAxis, setDepthAxis] = useState<DepthAxis>(init.depthAxis ?? 'shared');
  // ---- Sample groups (aggregate view): one pooled track per group ----
  const [groups, setGroups] = useState<SampleGroup[]>(() => (init.groups ?? []).map((g, i) => ({ id: i + 1, name: g.name, sampleIds: [...g.sampleIds] })));
  const groupIdSeq = useRef((init.groups?.length ?? 0) + 1);
  const [viewMode, setViewMode] = useState<'samples' | 'groups'>(init.viewMode === 'groups' && (init.groups ?? []).some(g => g.sampleIds.length) ? 'groups' : 'samples');
  const [showGroupsDialog, setShowGroupsDialog] = useState(false);
  /** Arc labels of the sample tracks: spliced reads, or the usage of each event against its canonical junction (the Groups view always shows usage). */
  const [arcLabel, setArcLabel] = useState<'reads' | 'usage'>(init.arcLabels ?? 'reads');
  const showUsage = viewMode === 'groups' || arcLabel === 'usage';
  const [uniqueOnly, setUniqueOnly] = useState(init.uniqueOnly ?? false);
  const [minJunctionCount, setMinJunctionCount] = useState(init.minJunctionReads ?? 3);
  /** In % usage mode, events below this usage are hidden (junctions without a share fall back to Min reads). */
  const [minUsagePct, setMinUsagePct] = useState(init.minUsagePct ?? 1);
  /** Count intron retention in the usage percentages (IR pills, and retention in the canonical arc's denominator). */
  const [includeRetention, setIncludeRetention] = useState(init.intronRetention ?? true);
  const [showReads, setShowReads] = useState(init.reads ?? !!initialReads);
  const [readsSampleId, setReadsSampleId] = useState<number | null>(init.readsSample ?? null);
  const [readsAll, setReadsAll] = useState(init.readsAll ?? false); // one reads track under every sample (primary only by default)
  const [collapseReads, setCollapseReads] = useState(init.collapseReads ?? false);
  const [minVafPct, setMinVafPct] = useState(init.minVafPct ?? 10); // variant sites need at least this alternate-allele fraction
  const [showAllTx, setShowAllTx] = useState(init.allTranscripts ?? false);
  const [showSnps, setShowSnps] = useState(init.commonSnps ?? false);
  const [snpMinAf, setSnpMinAf] = useState(init.snpMinAf ?? 0.01);
  const [snps, setSnps] = useState<{ chrom: string; start: number; end: number; list: CommonSnp[] } | null>(null);
  const [snpStatus, setSnpStatus] = useState<{ loading: boolean; error?: string }>({ loading: false });
  const snpSeq = useRef(0);
  const [altTx, setAltTx] = useState<{ geneName: string; data: AllTranscripts } | null>(null);
  const [altTxError, setAltTxError] = useState<string | undefined>();

  // ---- Data ----
  const [transcript, setTranscript] = useState<TranscriptData | null>(null);
  const [transcriptMissing, setTranscriptMissing] = useState<string | false>(false);
  const [tracks, setTracks] = useState<TrackData[]>([]);
  const [runSamples, setRunSamples] = useState<{ id: number; name: string }[]>([]);
  type ReadsEntry = { sampleId: number; fetched: FetchWindow; mode: 'reads' | 'collapsed'; minSupport: number; minVaf: number; data: ReadsResponse };
  // Per sample, so that "all samples" keeps one reads track under each coverage track
  const [readsData, setReadsData] = useState<Record<number, ReadsEntry>>({});
  const [readsLoading, setReadsLoading] = useState<Record<number, boolean>>({});
  const [readsError, setReadsError] = useState<Record<number, string | undefined>>({});
  const readsSeq = useRef(new Map<number, number>());

  // ---- UI state ----
  const [geneSearch, setGeneSearch] = useState('');
  const [geneSearchLoading, setGeneSearchLoading] = useState(false);
  const [showPicker, setShowPicker] = useState(false);
  /** Dropdowns open to the right of their button unless that would leave the window (the toolbar wraps on narrow screens). */
  const [pickerSide, setPickerSide] = useState<'left' | 'right'>('left');
  const openDropdown = (e: React.MouseEvent, width: number, toggle: (v: (p: boolean) => boolean) => void) => {
    const r = (e.currentTarget as HTMLElement).getBoundingClientRect();
    setPickerSide(r.left + width + 16 > window.innerWidth ? 'right' : 'left');
    toggle(p => !p);
  };
  // ---- GTEx tissue tracks (median junction reads + exon reads-per-base), chosen from the GTEx tissue list ----
  const [gtexTracks, setGtexTracks] = useState<TrackData[]>([]);
  const [gtexTissues, setGtexTissues] = useState<GtexTissue[] | null>(null);
  const [gtexError, setGtexError] = useState<string | undefined>();
  const [showGtexPicker, setShowGtexPicker] = useState(false);
  const [gtexSearch, setGtexSearch] = useState('');
  const [gtexFavourites, setGtexFavourites] = useState<string[]>(loadGtexFavourites);
  const gtexIdSeq = useRef(-1);
  const [pickerSearch, setPickerSearch] = useState('');
  const [dragging, setDragging] = useState(false);
  const [regionSelect, setRegionSelect] = useState<{ startX: number; currentX: number } | null>(null);
  const [hover, setHover] = useState<{ px: number; py: number } | null>(null);
  /** Detail popover opened by clicking a junction arc or an exon (HTML, never exported). */
  const [popover, setPopover] = useState<
    | { kind: 'junction'; key: string; j: JunctionArc; x: number; y: number }
    | { kind: 'exon'; exon: { start: number; end: number; rank: number }; x: number; y: number }
    | null>(null);
  const dragMoved = useRef(false);
  /** Splicing cartoon opened from a junction popover (experimental). */
  const [cartoon, setCartoon] = useState<{ j: JunctionArc; model: TxModel; label: string; color: string; sample: string } | null>(null);
  const [cartoonState, setCartoonState] = useState<{ story: SpliceStory | null; loading: string | null; error?: string }>({ story: null, loading: null });
  /** Per-exon depth of every sample of the run (the cohort of the exon-usage statistics), one fetch per gene. */
  const [exonDepths, setExonDepths] = useState<{ key: string; data: ExonUsageResponse | null; error?: string } | null>(null);
  const [junctionOffsets, setJunctionOffsets] = useState<Record<string, number>>({});
  const junctionDrag = useRef<{ key: string; startY: number; startOffset: number } | null>(null);
  const dragStart = useRef<{ x: number; vStart: number; vEnd: number } | null>(null);
  const rafPending = useRef<number | null>(null);
  const latestMouse = useRef<{ x: number; y: number }>({ x: 0, y: 0 });
  const [reloadTrigger, setReloadTrigger] = useState(0);

  // ---- Responsive width ----
  useEffect(() => {
    const obs = new ResizeObserver(entries => {
      for (const e of entries) setSvgWidth(Math.max(480, e.contentRect.width));
    });
    if (containerRef.current) obs.observe(containerRef.current);
    return () => obs.disconnect();
  }, []);

  // ---- Transcript model, axis and scale ----
  const tx: TxModel | null = useMemo(() => (transcript ? toTxModel(transcript) : null), [transcript]);
  const txRef = useRef(tx);
  txRef.current = tx;
  const reverse = tx?.strand === -1;
  const axis: VirtualAxis = useMemo(() => (equalIntrons && tx ? equalIntronAxis(tx, intronWidth) : LINEAR_AXIS), [equalIntrons, tx, intronWidth]);
  const plotWidth = svgWidth - PLOT_LEFT - PLOT_RIGHT_PAD;
  const scale: Scale = useMemo(
    () => makeScale(axis, viewStart, viewEnd, PLOT_LEFT, plotWidth, reverse),
    [axis, viewStart, viewEnd, plotWidth, reverse],
  );

  /** Default window for a gene in a given axis mode. */
  const defaultView = useCallback((gs: number, ge: number, ax: VirtualAxis): [number, number] => {
    if (ax.kind === 'equal-intron') {
      const pad = ax.intronV * 0.5;
      return [Math.max(0, Math.round(ax.fromV(ax.toV(gs) - pad))), Math.round(ax.fromV(ax.toV(ge) + pad))];
    }
    return linearDefault(gs, ge);
  }, []);

  /** Set the view window from virtual coordinates, enforcing zoom limits. */
  const setViewFromV = useCallback((vs: number, ve: number, ax: VirtualAxis = axis) => {
    if (ve - vs < MIN_V_SPAN) {
      const c = (vs + ve) / 2;
      vs = c - MIN_V_SPAN / 2; ve = c + MIN_V_SPAN / 2;
    }
    let gs = Math.round(ax.fromV(vs));
    let ge = Math.round(ax.fromV(ve));
    if (ge - gs > MAX_VIEW_BP) return;
    if (gs < 0) { ge -= gs; gs = 0; }
    if (ge <= gs) ge = gs + 1;
    setViewStart(gs);
    setViewEnd(ge);
  }, [axis]);

  const regionStr = useMemo(() => {
    const fmt = (n: number) => n.toLocaleString();
    return `${currentChrom}:${fmt(viewStart + 1)}-${fmt(viewEnd)} · ${formatBp(viewEnd - viewStart)}`;
  }, [currentChrom, viewStart, viewEnd]);

  // ---- Load transcript ----
  // The gene span handed by the host (outlier tables) lets the data source skip the symbol lookup.
  const hintRef = useRef<RegionHint | undefined>(chrom && geneEnd > geneStart ? { chrom, start: geneStart, end: geneEnd } : undefined);
  useEffect(() => {
    let cancelled = false;
    setTranscriptMissing(false);
    if (parseLocus(currentGeneName)) { setTranscript(null); setTranscriptMissing('no RefSeq gene at this locus'); return; }
    ds.getTranscript(currentGeneName, currentGeneId, hintRef.current)
      .then(t => { if (!cancelled) setTranscript(t); })
      .catch(e => { if (!cancelled) { setTranscript(null); setTranscriptMissing(String(e?.message || e || 'lookup failed').replace(/^.*?(No transcript model|not found)/, '$1')); } });
    return () => { cancelled = true; };
  }, [currentGeneName, currentGeneId]);

  // ---- All transcript models of the gene (on demand) ----
  useEffect(() => {
    if (!showAllTx || (altTx && altTx.geneName === currentGeneName)) return;
    let cancelled = false;
    setAltTxError(undefined);
    ds.getAllTranscripts(currentGeneName, currentGeneId, hintRef.current)
      .then(d => { if (!cancelled) setAltTx({ geneName: currentGeneName, data: d }); })
      .catch(e => { if (!cancelled) setAltTxError(e.message); });
    return () => { cancelled = true; };
  }, [showAllTx, currentGeneName, currentGeneId, altTx]);

  /** Transcript models in 0-based half-open coordinates, for the current gene. */
  const altModels = useMemo(() => {
    if (!showAllTx || !altTx || altTx.geneName !== currentGeneName) return null;
    return altTx.data.transcripts.slice(0, ALT_TX_MAX_ROWS).map(m => ({
      ...m,
      exons: [...m.exons].map(e => ({ start: e.start - 1, end: e.end })).sort((a, b) => a.start - b.start),
      start: m.start - 1,
      cdsStart: m.cds_start != null && m.cds_end != null ? m.cds_start - 1 : null,
      cdsEnd: m.cds_start != null && m.cds_end != null ? m.cds_end : null,
    }));
  }, [showAllTx, altTx, currentGeneName]);
  /** Make one of the listed transcripts the displayed reference model (exon numbering, junction classes, HGVS, usage percentages). */
  const applyModel = useCallback((t: TranscriptModel) => {
    const sorted = [...t.exons].sort((a, b) => a.start - b.start);
    setTranscript(prev => ({
      gene_name: prev?.gene_name ?? currentGeneName, transcript_id: t.id, translation_id: null, is_mane_select: t.is_mane,
      model_kind: t.is_mane ? 'mane' : 'chosen', biotype: t.biotype, source: t.source, chrom: currentChrom, strand: t.strand,
      start: t.start, end: t.end,
      exons: sorted.map((e, i) => ({ start: e.start, end: e.end, rank: t.strand > 0 ? i + 1 : sorted.length - i })),
      cds_start: t.cds_start ?? null, cds_end: t.cds_end ?? null,
    }));
  }, [currentGeneName, currentChrom]);
  const chooseModel = useCallback((id: string) => {
    const t = altTx?.data.transcripts.find(x => x.id === id);
    if (t) applyModel(t);
  }, [altTx, applyModel]);
  // A session that chose a reference transcript: look it up in the gene's transcript list once the default model is in
  const wantedTxRef = useRef(init.transcriptId);
  useEffect(() => {
    const want = wantedTxRef.current;
    if (!want || !tx || tx.transcriptId === want) return;
    let cancelled = false;
    const list = altTx && altTx.geneName === currentGeneName ? Promise.resolve(altTx.data) : ds.getAllTranscripts(currentGeneName, currentGeneId, hintRef.current);
    list.then(d => {
      if (cancelled) return;
      wantedTxRef.current = undefined;
      const t = d.transcripts.find(x => x.id === want);
      if (t) applyModel(t);
    }).catch(() => { wantedTxRef.current = undefined; });
    return () => { cancelled = true; };
  }, [tx, altTx, currentGeneName, currentGeneId, applyModel]);
  /** junction key → transcript ids whose consecutive exons form that intron (for arc tooltips). */
  const altJunctionIndex = useMemo(() => {
    const idx = new Map<string, string[]>();
    for (const m of altModels || []) {
      for (let i = 0; i + 1 < m.exons.length; i++) {
        const k = `${m.exons[i].end}-${m.exons[i + 1].start}`;
        idx.set(k, [...(idx.get(k) || []), m.id]);
      }
    }
    return idx;
  }, [altModels]);
  const altPanelH = showAllTx ? ALT_TX_HEADER_H + Math.max(1, altModels?.length ?? 1) * ALT_TX_ROW_H + 8 : 0;

  // ---- Common SNPs of the window (browser-side, cached per 100 kb chunk) ----
  useEffect(() => {
    if (!showSnps) return;
    if (viewEnd - viewStart > SNP_MAX_WINDOW) { setSnpStatus({ loading: false, error: `zoom in below ${formatBp(SNP_MAX_WINDOW)} to load common SNPs` }); return; }
    if (snps && snps.chrom === currentChrom && snps.start <= viewStart && snps.end >= viewEnd) return;
    const seq = ++snpSeq.current;
    const margin = Math.min(200_000, viewEnd - viewStart);
    const s0 = Math.max(0, viewStart - margin), e0 = viewEnd + margin;
    const timer = window.setTimeout(() => {
      setSnpStatus({ loading: true });
      ds.getCommonSnps(currentChrom, s0, e0)
        .then(list => { if (seq === snpSeq.current) { setSnps({ chrom: currentChrom, start: s0, end: e0, list }); setSnpStatus({ loading: false }); } })
        .catch(e => { if (seq === snpSeq.current) setSnpStatus({ loading: false, error: e?.message || String(e) }); });
    }, 250);
    return () => window.clearTimeout(timer);
  }, [showSnps, currentChrom, viewStart, viewEnd, snps]); // eslint-disable-line react-hooks/exhaustive-deps

  /** Common SNPs in the window above the frequency threshold, and a position index for the read variant sites. */
  const visibleSnps = useMemo(() => {
    if (!showSnps || !snps || snps.chrom !== currentChrom) return [];
    return snps.list.filter(v => v.maxAf >= snpMinAf && v.end > viewStart && v.start < viewEnd);
  }, [showSnps, snps, currentChrom, snpMinAf, viewStart, viewEnd]);
  const snpByPos = useMemo(() => {
    const m = new Map<number, CommonSnp[]>();
    for (const v of visibleSnps) m.set(v.start, [...(m.get(v.start) || []), v]);
    return m;
  }, [visibleSnps]);
  /** The common SNP matching a variant site called from the reads (same position, same alternate allele for SNVs). */
  const knownSnp = useCallback((st: VariantSite): CommonSnp | null => {
    const at = snpByPos.get(st.pos) || [];
    if (st.kind === 'snv') return at.find(v => v.cls === 'snv' && v.alts.includes(st.alt)) ?? null;
    return at.find(v => v.cls !== 'snv') ?? at[0] ?? null;
  }, [snpByPos]);
  const snpText = (v: CommonSnp) => {
    const af = (x: number) => `${(x * 100).toFixed(x < 0.01 ? 2 : 1)}%`;
    const main = v.afs.filter(a => ['1000Genomes', 'TOPMED', 'GnomAD', 'GnomAD_exomes'].includes(a.source));
    const shown = (main.length ? main : v.afs.slice(0, 3)).map(a => `${a.source} ${af(a.af)}`).join(', ');
    return `${v.id} · ${v.ref || '-'}>${v.alts.join('/') || '-'} · ${v.cls}${v.impact ? ` · ${v.impact}` : ''}\nmax AF ${af(v.maxAf)}${shown ? ` (${shown})` : ''}`;
  };

  // ---- Known variants of the loaded samples (clinical indication, diagnostic, chromosome map) ----
  const [knownVariants, setKnownVariants] = useState<Map<number, KnownVariant[]>>(() => new Map());
  const [showKnown, setShowKnown] = useState(init.knownVariants ?? true);
  const knownRequested = useRef<Set<number>>(new Set());
  useEffect(() => {
    if (!ds.getKnownVariants) return;
    // Without any alignment (sampleId 0) the host may still hand over variants, e.g. from a deep link
    const ids = tracks.length ? tracks.map(t => t.sampleId) : [sampleId];
    for (const sid of ids) {
      if (sid < 0 || knownRequested.current.has(sid)) continue;
      knownRequested.current.add(sid);
      ds.getKnownVariants(sid)
        .then(list => setKnownVariants(prev => new Map(prev).set(sid, list)))
        .catch(e => { knownRequested.current.delete(sid); console.warn('[sashimi] known variants unavailable:', e?.message || e); });
    }
  }, [tracks, ds, sampleId]);
  const chromKey = (c: string) => (c.startsWith('chr') ? c : `chr${c}`).replace(/^chrMT$/, 'chrM');
  /** Variants of a sample placed on the current chromosome (a bare g. notation of the queried gene counts as here). */
  const knownOnChrom = useCallback((sid: number): KnownVariant[] => {
    const here = chromKey(currentChrom);
    return (knownVariants.get(sid) ?? []).filter(v => (v.chrom ? chromKey(v.chrom) === here : v.gene?.toUpperCase() === currentGeneName.toUpperCase()));
  }, [knownVariants, currentChrom, currentGeneName]);
  const primaryId = tracks[0]?.sampleId ?? sampleId;
  const primaryKnown = knownVariants.get(primaryId) ?? [];
  const primaryKnownHere = useMemo(() => knownOnChrom(primaryId), [knownOnChrom, primaryId]);
  const primaryKnownElsewhere = primaryKnown.length - primaryKnownHere.length;
  /** Title line of the known-variant panel: the sample and how many of its variants fall on this chromosome. */
  const knownStatus = `${tracks[0]?.sampleName ?? sampleName} · ${primaryKnownHere.length} on ${currentChrom}${primaryKnownElsewhere ? ` · ${primaryKnownElsewhere} elsewhere (${[...new Set(primaryKnown.filter(v => !primaryKnownHere.includes(v)).map(v => v.chrom || '?'))].join(', ')})` : ''}`;
  /** Rows of the known-variant panel: variants stacked so their marks and labels do not overlap on screen.
   *  The first row is the title line: variants that would sit under the title text move to the next row. */
  const knownRows = useMemo((): KnownVariant[][] => {
    if (!showKnown || !primaryKnown.length) return [];
    const titleW = 8 + 14 * 6.2 + knownStatus.length * 5.3 + 8;
    const rows: { items: KnownVariant[]; spans: [number, number][] }[] = [{ items: [], spans: [[PLOT_LEFT, PLOT_LEFT + titleW]] }];
    const sorted = [...primaryKnownHere].sort((a, b) => a.start - b.start);
    for (const v of sorted) {
      const xa = scale.x(v.start), xb = scale.x(v.end);
      const left = Math.max(PLOT_LEFT, Math.min(xa, xb)), right = Math.min(PLOT_LEFT + plotWidth, Math.max(xa, xb));
      const lo = Math.min(left, right) - 6, hi = Math.max(left + 8, right) + v.label.length * 5.4 + 12;
      let row = rows.find(r => r.spans.every(([a, b]) => hi < a || lo > b));
      if (!row) { row = { items: [], spans: [] }; rows.push(row); }
      row.items.push(v); row.spans.push([lo, hi]);
    }
    return rows.map(r => r.items);
  }, [showKnown, primaryKnown.length, primaryKnownHere, scale, plotWidth, knownStatus]);
  const knownPanelH = knownRows.length ? KNOWN_PAD * 2 + knownRows.length * KNOWN_ROW_H : 0;
  /** Vertical centre of row `i` of the known-variant panel, relative to the panel top. */
  const knownRowMid = (i: number) => KNOWN_PAD + i * KNOWN_ROW_H + KNOWN_ROW_H / 2;
  /** Centre the window on a known variant (bands get a 10 % margin, points a 1 kb window at most). */
  const jumpToVariant = useCallback((v: KnownVariant) => {
    const span = v.end - v.start;
    let s: number, e: number;
    if (isPointVariant(v)) { const half = Math.min(500, Math.max(MIN_V_SPAN, (viewEnd - viewStart) / 2)); s = Math.floor((v.start + v.end) / 2 - half); e = Math.ceil((v.start + v.end) / 2 + half); }
    else { const m = Math.max(50, Math.round(span * 0.1)); s = v.start - m; e = v.end + m; }
    if (e - s > MAX_VIEW_BP) { const c = (s + e) / 2; s = Math.round(c - MAX_VIEW_BP / 2); e = Math.round(c + MAX_VIEW_BP / 2); }
    setViewFromV(axis.toV(Math.max(0, s)), axis.toV(e));
  }, [axis, setViewFromV, viewStart, viewEnd]);

  // ---- Neighbouring genes: canonical transcript of every other gene overlapping the window ----
  // The data sources cache 500 kb chunks, so panning only costs a request when a new chunk is entered.
  interface NeighbourModel extends TxModel { geneId: string; biotype: string; isCanonical: boolean }
  const [neighbours, setNeighbours] = useState<{ chrom: string; start: number; end: number; models: NeighbourModel[] } | null>(null);
  const [neighbourError, setNeighbourError] = useState<string | undefined>();
  const neighbourSeq = useRef(0);
  useEffect(() => {
    if (!tx) return;
    // Nothing to look for while the window sits inside the queried gene body (no other gene can show its exons there
    // unless it overlaps the gene, which the chunk cache makes cheap to check anyway on the first pan).
    const seq = ++neighbourSeq.current;
    const timer = window.setTimeout(() => {
      ds.getRegionGenes(currentChrom, viewStart + 1, viewEnd, currentGeneName)
        .then(list => {
          if (seq !== neighbourSeq.current) return;
          const models: NeighbourModel[] = list.map(g => toNeighbourModel(g, currentChrom));
          setNeighbours({ chrom: currentChrom, start: viewStart, end: viewEnd, models });
          setNeighbourError(undefined);
        })
        .catch(e => { if (seq === neighbourSeq.current) setNeighbourError(e.message); });
    }, 250);
    return () => window.clearTimeout(timer);
  }, [tx, currentChrom, viewStart, viewEnd, currentGeneName]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- Exon usage cohort: fetched when a popover first needs it, cached per gene / read filter ----
  const usageKey = tx ? `${runId}|${tx.transcriptId}|${uniqueOnly ? 'u' : 'a'}` : null;
  const usageRequested = useRef<string | null>(null);
  const popoverOpen = popover != null;
  useEffect(() => {
    if (!popoverOpen || !tx || !usageKey || usageRequested.current === usageKey) return;
    usageRequested.current = usageKey;
    const key = usageKey;
    setExonDepths({ key, data: null });
    const ivs = usageIntervals(tx);
    ds.getExonUsage(runId, currentChrom, tx.strand, ivs.map(iv => [iv.start, iv.end] as [number, number]), uniqueOnly)
      .then(d => setExonDepths(prev => (prev?.key === key ? { key, data: d } : prev)))
      .catch(e => {
        if (usageRequested.current === key) usageRequested.current = null;   // allow a retry on the next click
        setExonDepths(prev => (prev?.key === key ? { key, data: null, error: e?.message || String(e) } : prev));
      });
  }, [popoverOpen, usageKey]); // eslint-disable-line react-hooks/exhaustive-deps

  /**
   * Depth-based usage of one exon (by index in the MANE model) for every sample of the run:
   * usage = median depth of the exon / median of the median depths of the other coding exons.
   * Controls of a sample = all other samples of the run; relative = usage / median(controls).
   */
  const usageOf = useCallback((idx: number) => {
    const data = exonDepths?.data;
    if (!tx || !data || exonDepths?.key !== usageKey) return null;
    const ivs = usageIntervals(tx);
    const refIdx = referenceExons(ivs, idx);
    const per = data.samples.filter(smp => !smp.error && smp.exons.length === ivs.length)
      .map(smp => ({ sample: smp, u: exonUsage(smp.exons, idx, refIdx) }));
    const values = per.map(p => p.u.usage).filter((v): v is number => v != null);
    const forSample = (sampleId: number) => {
      const me = per.find(p => p.sample.sample_id === sampleId);
      const controls = usageCohort(per.filter(p => p.sample.sample_id !== sampleId).map(p => p.u.usage).filter((v): v is number => v != null));
      const rel = me?.u.usage != null && controls.median ? me.u.usage / controls.median : null;
      return { me, controls, rel, z: usageZ(me?.u.usage ?? null, controls) };
    };
    return { interval: ivs[idx], refIdx, per, values, forSample, failed: data.samples.filter(smp => smp.error).length };
  }, [exonDepths, tx, usageKey]);

  // ---- Cartoon data: reference sequence of the blocks of both isoforms + protein domains, then the story ----
  useEffect(() => {
    if (!cartoon) { setCartoonState({ story: null, loading: null }); return; }
    let cancelled = false;
    const { j, model } = cartoon;
    const others = tracks.flatMap(t => t.junctions);
    const event = spliceEvent(j, model, others);
    setCartoonState({ story: null, loading: 'fetching the reference sequence…' });
    (async () => {
      const windows = storyWindows(event);
      const seqs = await Promise.all(windows.map(([a, b]) => ds.getReference(currentChrom, a, b)));
      if (cancelled) return;
      if (seqs.some(x => x == null)) { setCartoonState({ story: null, loading: null, error: 'No reference sequence available for this gene (set REFERENCE_FASTA or allow access to the UCSC / Ensembl APIs).' }); return; }
      const seqAt = (a: number, b: number): string => {
        for (let i = 0; i < windows.length; i++) if (a >= windows[i][0] && b <= windows[i][1]) return seqs[i]!.slice(a - windows[i][0], b - windows[i][0]);
        return 'N'.repeat(b - a);
      };
      setCartoonState({ story: null, loading: 'fetching protein domains…' });
      let features: import('./sashimi/types').ProteinDomain[] = [];
      if (model.cdsStart != null && model.cdsEnd != null) {
        try {
          features = await ds.getProteinDomains({ transcriptId: model.transcriptId, translationId: model.translationId, chrom: currentChrom, strand: model.strand, exons: model.exons, cdsStart: model.cdsStart, cdsEnd: model.cdsEnd });
        } catch { features = []; }
      }
      if (cancelled) return;
      try {
        setCartoonState({ story: spliceStory(j, model, others, seqAt, features), loading: null });
      } catch (e: any) {
        setCartoonState({ story: null, loading: null, error: e?.message || String(e) });
      }
    })();
    return () => { cancelled = true; };
  }, [cartoon]); // eslint-disable-line react-hooks/exhaustive-deps

  /** Neighbouring genes visible in the window, packed into rows (a row holds genes that do not overlap on screen). */
  const neighbourRows = useMemo((): NeighbourModel[][] => {
    if (!neighbours || neighbours.chrom !== currentChrom || !tx) return [];
    const plotRight = PLOT_LEFT + plotWidth;
    const visible = neighbours.models
      .filter(m => m.end > viewStart && m.start < viewEnd && m.geneName !== tx.geneName)
      .sort((a, b) => a.start - b.start);
    const rows: { items: NeighbourModel[]; spans: [number, number][] }[] = [];
    for (const m of visible) {
      const a = scale.x(m.start), b = scale.x(m.end);
      const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
      const lw = (m.geneName.length + 2) * 5.6 + 10;
      const inside = right - left > lw + 60;
      const span: [number, number] = inside ? [left - 2, right + 2] : right + 4 + lw <= plotRight ? [left - 2, right + lw + 6] : [left - lw - 6, right + 2]; // room for the name pill
      let row = rows.find(r => r.spans.every(([s0, e0]) => span[1] < s0 || span[0] > e0));
      if (!row) { row = { items: [], spans: [] }; rows.push(row); }
      row.items.push(m); row.spans.push(span);
    }
    return rows.map(r => r.items);
  }, [neighbours, currentChrom, tx, viewStart, viewEnd, scale, plotWidth]);
  const transcriptPanelH = TRANSCRIPT_H + neighbourRows.length * NEIGHBOUR_ROW_H;

  /**
   * Which gene model a junction belongs to. The queried gene is tried first; when it does not
   * annotate both ends, the neighbouring genes are tried and the best match wins (canonical >
   * exon skipping > one annotated site > novel). Junctions of a neighbour are then labelled and
   * framed with that gene's model, and their canonical ones stop being reported as "novel".
   */
  const junctionContext = useCallback((j: JunctionArc): { model: TxModel | null; info: JunctionInfo; foreign: NeighbourModel | null } => {
    const score = (i: JunctionInfo) => i.cls === 'canonical' ? 3 : i.cls === 'exon_skipping' ? 2 : i.cls === 'novel' ? 0 : 1;
    const main = classifyJunction(j, tx);
    let best = { model: tx, info: main, foreign: null as NeighbourModel | null, score: score(main) };
    if (best.score < 3) {
      for (const m of neighbours?.models || []) {
        if (m.end <= j.start || m.start >= j.end) continue;
        const info = classifyJunction(j, m);
        const sc = score(info);
        if (sc > best.score) best = { model: m, info: { ...info, label: `${m.geneName} · ${info.label}` }, foreign: m, score: sc };
      }
    }
    return { model: best.model, info: best.info, foreign: best.foreign };
  }, [tx, neighbours]);

  // ---- Coverage loading (with margin, stale-response protection) ----
  const viewRef = useRef({ chrom: currentChrom, start: viewStart, end: viewEnd, uniqueOnly });
  viewRef.current = { chrom: currentChrom, start: viewStart, end: viewEnd, uniqueOnly };
  const reqSeq = useRef<Map<number, number>>(new Map());

  const fetchWindowFor = (v: { chrom: string; start: number; end: number; uniqueOnly: boolean }): FetchWindow => {
    const span = v.end - v.start;
    const margin = Math.min(span, Math.max(0, Math.floor((MAX_FETCH_BP - span) / 2)));
    return { chrom: v.chrom, start: Math.max(0, v.start - margin), end: v.end + margin, uniqueOnly: v.uniqueOnly };
  };
  const covers = (f: FetchWindow | undefined, v: { chrom: string; start: number; end: number; uniqueOnly: boolean }) =>
    !!f && f.chrom === v.chrom && f.uniqueOnly === v.uniqueOnly && f.start <= v.start && f.end >= v.end;

  const loadCoverage = useCallback(async (sid: number, sname: string) => {
    const view = viewRef.current;
    const win = fetchWindowFor(view);
    const seq = (reqSeq.current.get(sid) || 0) + 1;
    reqSeq.current.set(sid, seq);
    setTracks(prev => {
      const existing = prev.find(t => t.sampleId === sid);
      if (existing) return prev.map(t => t.sampleId === sid ? { ...t, loading: true, error: undefined } : t);
      return [...prev, { sampleId: sid, sampleName: sname, coverage: [], junctions: [], loading: true }];
    });
    try {
      const data = await ds.getCoverage(sid, win.chrom, win.start, win.end, win.uniqueOnly, boundariesOf(txRef.current),
        { core: { start: view.start, end: view.end }, maxReads: MAX_READS_PER_TRACK });
      if (reqSeq.current.get(sid) !== seq) return; // a newer request superseded this one
      // the source may have read less margin than asked for (deep library): remember what it really covered
      const fetched: FetchWindow = data.window ? { ...win, start: data.window.start, end: data.window.end } : win;
      setTracks(prev => prev.map(t => t.sampleId === sid ? {
        ...t, coverage: data.coverage, junctions: data.junctions, spanning: data.spanning, sampled: data.sampled, loading: false, error: data.error, fetched,
      } : t));
    } catch (err: any) {
      if (reqSeq.current.get(sid) !== seq) return;
      setTracks(prev => prev.map(t => t.sampleId === sid ? { ...t, loading: false, error: err.message } : t));
    }
  }, []);

  // ---- Initial load: primary sample + one random comparison ----
  const didInit = useRef(false);
  useEffect(() => {
    if (didInit.current) return;
    didInit.current = true;
    if (sampleId > 0) {
      loadCoverage(sampleId, sampleName);
      ds.getRandomSample(runId, sampleId)
        .then(s => loadCoverage(s.id, s.name))
        .catch(() => {}); // no other samples is fine
    }
    ds.getRunSamples(runId).then(setRunSamples).catch(() => {});
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  // Renamed samples: relabel the tracks in place
  useEffect(() => {
    if (!sampleNames) return;
    setTracks(prev => prev.some(t => sampleNames[t.sampleId] && sampleNames[t.sampleId] !== t.sampleName)
      ? prev.map(t => (sampleNames[t.sampleId] && sampleNames[t.sampleId] !== t.sampleName ? { ...t, sampleName: sampleNames[t.sampleId] } : t))
      : prev);
  }, [sampleNames]);

  /** Make a sample the primary track (first, reference of the comparisons); loads it first when not yet a track. */
  const primaryRef = useRef(sampleId);
  const setPrimary = useCallback((sid: number, name?: string) => {
    primaryRef.current = sid;
    if (tracksRef.current.some(t => t.sampleId === sid)) {
      setTracks(prev => [...prev.filter(t => t.sampleId === sid), ...prev.filter(t => t.sampleId !== sid)]);
    } else if (name) {
      loadCoverage(sid, name);
      setTracks(prev => [...prev.filter(t => t.sampleId === sid), ...prev.filter(t => t.sampleId !== sid)]);
    }
    setJunctionOffsets({});
  }, [loadCoverage]);
  // the host changed the primary sample (standalone chips): promote it without remounting
  useEffect(() => {
    if (primaryRef.current === sampleId) return;
    primaryRef.current = sampleId;
    if (sampleId > 0) setPrimary(sampleId, sampleName);
  }, [sampleId, sampleName, setPrimary]);

  // ---- Reload tracks whose fetched window no longer covers the view (debounced) ----
  const tracksRef = useRef(tracks);
  tracksRef.current = tracks;
  useEffect(() => {
    if (tracks.length === 0) return;
    const timer = setTimeout(() => {
      for (const t of tracksRef.current) {
        if (!covers(t.fetched, viewRef.current)) loadCoverage(t.sampleId, t.sampleName);
      }
    }, 250);
    return () => clearTimeout(timer);
  }, [viewStart, viewEnd, currentChrom, uniqueOnly, reloadTrigger]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- Reads track loading (primary sample by default, or every sample; only below the visibility window) ----
  const effectiveReadsSampleId = tracks.some(t => t.sampleId === readsSampleId) ? readsSampleId : (tracks[0]?.sampleId ?? null);
  // ---- Sample groups ----
  const addGroup = useCallback(() => setGroups(prev => [...prev, { id: groupIdSeq.current++, name: `Group ${prev.length + 1}`, sampleIds: [] }]), []);
  const renameGroup = useCallback((id: number, name: string) => setGroups(prev => prev.map(g => g.id === id ? { ...g, name } : g)), []);
  const deleteGroup = useCallback((id: number) => setGroups(prev => prev.filter(g => g.id !== id)), []);
  const removeFromGroup = useCallback((id: number, sid: number) => setGroups(prev => prev.map(g => g.id === id ? { ...g, sampleIds: g.sampleIds.filter(x => x !== sid) } : g)), []);
  /** A sample belongs to one group: adding it moves it out of any other. Its coverage loads at once so the group track can be drawn. */
  const addToGroup = useCallback((id: number, sid: number) => {
    setGroups(prev => prev.map(g => g.id === id
      ? { ...g, sampleIds: g.sampleIds.includes(sid) ? g.sampleIds : [...g.sampleIds, sid] }
      : { ...g, sampleIds: g.sampleIds.filter(x => x !== sid) }));
    if (!tracksRef.current.some(t => t.sampleId === sid)) { const s = runSamples.find(x => x.id === sid); if (s) loadCoverage(s.id, s.name); }
  }, [runSamples, loadCoverage]);
  // Aggregate view: every member needs its coverage (loaded as a sample track, hidden while the groups are shown)
  useEffect(() => {
    if (viewMode !== 'groups') return;
    if (!groups.length) { setViewMode('samples'); return; }
    for (const g of groups) for (const sid of g.sampleIds) {
      if (tracksRef.current.some(t => t.sampleId === sid)) continue;
      const s = runSamples.find(x => x.id === sid);
      if (s) loadCoverage(s.id, s.name);
    }
  }, [viewMode, groups, runSamples, loadCoverage]);

  // A new gene model: tracks whose unspliced-read counts miss one of its boundaries reload (sources that never count them are left alone)
  useEffect(() => {
    const b = boundariesOf(tx);
    if (!b) return;
    for (const t of tracksRef.current) {
      if (t.gtex || t.loading || !t.fetched || !t.spanning) continue;
      const f = t.fetched;
      const missing = b.intronStarts.some(p => p >= f.start && p < f.end && t.spanning!.intronStart[p] == null)
        || b.intronEnds.some(p => p > f.start && p <= f.end && t.spanning!.intronEnd[p] == null);
      if (missing) loadCoverage(t.sampleId, t.sampleName);
    }
  }, [tx, loadCoverage]);

  /** Samples whose reads are shown, in track order. */
  const readsSampleIds = useMemo(() => !showReads ? [] : readsAll ? tracks.map(t => t.sampleId) : effectiveReadsSampleId == null ? [] : [effectiveReadsSampleId],
    [showReads, readsAll, tracks, effectiveReadsSampleId]);
  useEffect(() => {
    if (!readsSampleIds.length) return;
    const v = viewRef.current;
    const span = v.end - v.start;
    if (span > READS_MAX_VIEW_BP) return;
    const mode = collapseReads ? 'collapsed' : 'reads';
    const minVaf = Math.min(1, Math.max(0, minVafPct / 100));
    const stale = readsSampleIds.filter(sid => {
      const cur = readsData[sid];
      return !(cur && cur.mode === mode && cur.minVaf === minVaf && (mode === 'reads' || cur.minSupport === minJunctionCount) && covers(cur.fetched, v));
    });
    if (!stale.length) return;
    // Collapsed groups are computed for the exact window (counts are per window); raw reads get a pan margin
    const margin = collapseReads ? 0 : Math.floor(span * 0.25);
    const want: FetchWindow = { chrom: v.chrom, start: Math.max(0, v.start - margin), end: v.end + margin, uniqueOnly: v.uniqueOnly };
    const timer = setTimeout(() => {
      for (const sid of stale) {
        const seq = (readsSeq.current.get(sid) ?? 0) + 1;
        readsSeq.current.set(sid, seq);
        setReadsLoading(p => ({ ...p, [sid]: true }));
        setReadsError(p => ({ ...p, [sid]: undefined }));
        ds.getReads(sid, want.chrom, want.start, want.end, want.uniqueOnly, READS_MAX, mode, minJunctionCount, minVaf)
          .then(data => { if (readsSeq.current.get(sid) === seq) setReadsData(p => ({ ...p, [sid]: { sampleId: sid, fetched: want, mode, minSupport: minJunctionCount, minVaf, data } })); })
          .catch((err: any) => { if (readsSeq.current.get(sid) === seq) setReadsError(p => ({ ...p, [sid]: err.message })); })
          .finally(() => { if (readsSeq.current.get(sid) === seq) setReadsLoading(p => ({ ...p, [sid]: false })); });
      }
    }, 250);
    return () => clearTimeout(timer);
  }, [readsSampleIds, viewStart, viewEnd, currentChrom, uniqueOnly, readsData, collapseReads, minJunctionCount, minVafPct]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- Mouse interaction ----
  const svgPoint = (e: { clientX: number; clientY: number }) => {
    const rect = svgRef.current?.getBoundingClientRect();
    return { x: rect ? e.clientX - rect.left : 0, y: rect ? e.clientY - rect.top : 0 };
  };

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const { x } = svgPoint(e);
    setPopover(null);
    if (e.ctrlKey || e.metaKey) {
      setRegionSelect({ startX: x, currentX: x });
    } else {
      setDragging(true);
      dragStart.current = { x: e.clientX, vStart: scale.vStart, vEnd: scale.vEnd };
    }
    setHover(null);
  }, [scale]);

  const handleMouseMove = useCallback((e: React.MouseEvent) => {
    latestMouse.current = { x: e.clientX, y: e.clientY };
    if (rafPending.current !== null) return;
    rafPending.current = requestAnimationFrame(() => {
      rafPending.current = null;
      const m = latestMouse.current;
      if (junctionDrag.current) {
        const dy = m.y - junctionDrag.current.startY;
        if (Math.abs(dy) > 3) dragMoved.current = true;
        const key = junctionDrag.current.key, base = junctionDrag.current.startOffset;
        setJunctionOffsets(prev => ({ ...prev, [key]: base + dy }));
        return;
      }
      const p = svgPoint({ clientX: m.x, clientY: m.y });
      if (regionSelect) {
        setRegionSelect(prev => prev ? { ...prev, currentX: p.x } : null);
        return;
      }
      if (dragging && dragStart.current) {
        const dx = m.x - dragStart.current.x;
        const vSpan = dragStart.current.vEnd - dragStart.current.vStart;
        const dv = (reverse ? 1 : -1) * (dx / plotWidth) * vSpan;
        setViewFromV(dragStart.current.vStart + dv, dragStart.current.vEnd + dv);
        return;
      }
      setHover(p.x >= PLOT_LEFT && p.x <= PLOT_LEFT + plotWidth ? { px: p.x, py: p.y } : null);
    });
  }, [dragging, regionSelect, reverse, plotWidth, setViewFromV]);

  const handleMouseUp = useCallback(() => {
    if (junctionDrag.current) { junctionDrag.current = null; return; }
    if (regionSelect) {
      const x1 = Math.min(regionSelect.startX, regionSelect.currentX);
      const x2 = Math.max(regionSelect.startX, regionSelect.currentX);
      if (x2 - x1 > 5) {
        const va = scale.pxToV(x1), vb = scale.pxToV(x2);
        setViewFromV(Math.min(va, vb), Math.max(va, vb));
      }
      setRegionSelect(null);
      return;
    }
    setDragging(false);
    dragStart.current = null;
  }, [regionSelect, scale, setViewFromV]);

  const handleMouseLeave = useCallback(() => { handleMouseUp(); setHover(null); }, [handleMouseUp]);

  // Ctrl+wheel zoom around the cursor. Registered natively because React's wheel
  // listener is passive and cannot prevent the page from scrolling.
  const scaleRef = useRef(scale);
  scaleRef.current = scale;
  const setViewFromVRef = useRef(setViewFromV);
  setViewFromVRef.current = setViewFromV;
  useEffect(() => {
    const el = svgRef.current;
    if (!el) return;
    const onWheel = (e: WheelEvent) => {
      if (!(e.ctrlKey || e.metaKey)) return;
      e.preventDefault();
      const sc = scaleRef.current;
      const rect = el.getBoundingClientRect();
      const px = Math.min(sc.plotRight, Math.max(sc.plotLeft, e.clientX - rect.left));
      const vAt = sc.pxToV(px);
      const factor = e.deltaY > 0 ? 1.18 : 1 / 1.18;
      setViewFromVRef.current(vAt - (vAt - sc.vStart) * factor, vAt + (sc.vEnd - vAt) * factor);
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, []);

  const zoomBy = useCallback((factor: number) => {
    const c = (scale.vStart + scale.vEnd) / 2, half = (scale.vEnd - scale.vStart) / 2 * factor;
    setViewFromV(c - half, c + half);
  }, [scale, setViewFromV]);

  const resetZoom = useCallback(() => {
    const [s, e] = defaultView(currentGeneStart, currentGeneEnd, axis);
    setViewStart(s); setViewEnd(e);
  }, [currentGeneStart, currentGeneEnd, axis, defaultView]);

  // Keeps the current genomic window: the same bases stay in view, only the axis changes.
  const toggleEqualIntrons = useCallback((on: boolean) => {
    setEqualIntrons(on);
    setJunctionOffsets({});
  }, []);

  // ---- Track management ----
  const removeTrack = useCallback((sid: number) => {
    if (sid < 0) { setGtexTracks(prev => prev.filter(t => t.sampleId !== sid)); return; }
    setTracks(prev => prev.filter(t => t.sampleId !== sid));
    reqSeq.current.set(sid, (reqSeq.current.get(sid) || 0) + 1); // drop in-flight responses
  }, []);

  /** Load (or reload for a new gene) the GTEx profile of one tissue track. */
  const loadGtex = useCallback((sid: number, tissue: GtexTissue, geneName: string, geneId?: string) => {
    setGtexTracks(prev => prev.map(t => (t.sampleId === sid ? { ...t, loading: true, error: undefined } : t)));
    ds.getGtexProfile(geneName, geneId, tissue)
      .then(p => setGtexTracks(prev => prev.map(t => (t.sampleId === sid ? {
        ...t, loading: false, coverage: [], junctions: p.junctions,
        gtex: { tissue, dataset: p.dataset, unit: p.unit, warning: p.warning, tpm: p.tpm, lowCoverage: p.tpm != null && p.tpm < GTEX_MIN_TPM },
      } : t))))
      .catch(e => setGtexTracks(prev => prev.map(t => (t.sampleId === sid ? { ...t, loading: false, coverage: [], junctions: [], error: e?.message || String(e) } : t))));
  }, [ds]);
  const addGtexTissue = useCallback((tissue: GtexTissue) => {
    setShowGtexPicker(false); setGtexSearch('');
    if (gtexTracks.some(t => t.gtex?.tissue.id === tissue.id)) return;
    const sid = gtexIdSeq.current--;
    setGtexTracks(prev => [...prev, { sampleId: sid, sampleName: tissue.name, coverage: [], junctions: [], loading: true, gtex: { tissue, dataset: 'gtex_v10', unit: '', tpm: null, lowCoverage: false } }]);
    setGtexFavourites(prev => { const next = [tissue.id, ...prev.filter(x => x !== tissue.id)].slice(0, 8); saveGtexFavourites(next); return next; });
    loadGtex(sid, tissue, currentGeneName, currentGeneId);
  }, [gtexTracks, loadGtex, currentGeneName, currentGeneId]);
  useEffect(() => {
    if (!showGtexPicker || gtexTissues) return;
    setGtexError(undefined);
    ds.getGtexTissues().then(setGtexTissues).catch(e => setGtexError(e?.message || String(e)));
  }, [showGtexPicker, gtexTissues]); // eslint-disable-line react-hooks/exhaustive-deps
  // a new gene: every tissue track is reloaded for it
  const gtexGeneRef = useRef(currentGeneName);
  useEffect(() => {
    if (gtexGeneRef.current === currentGeneName) return;
    gtexGeneRef.current = currentGeneName;
    for (const t of gtexTracks) if (t.gtex) loadGtex(t.sampleId, t.gtex.tissue, currentGeneName, currentGeneId);
  }, [currentGeneName, currentGeneId, gtexTracks, loadGtex]);
  /**
   * Sample tracks followed by the GTEx tissue tracks, as drawn. A tissue track's profile is the usage of
   * each MANE exon (0–1) computed from the tissue's junction medians at both exon boundaries; tissues
   * where the gene's median TPM is below GTEX_MIN_TPM show nothing but "low coverage".
   */
  /** One pooled track per sample group: summed coverage, summed junction reads, per-intron shares of every splicing event. */
  const groupTracks = useMemo((): TrackData[] => groups.map(g => {
    const members = g.sampleIds.map(sid => tracks.find(t => t.sampleId === sid)).filter((t): t is TrackData => !!t);
    const { junctions, samplesWith } = poolJunctions(members);
    const spanning = poolSpanning(members);
    const failed = members.filter(m => m.error && !m.coverage.length);
    const pending = g.sampleIds.filter(sid => !members.some(m => m.sampleId === sid) && runSamples.some(x => x.id === sid));
    const sampledMembers = members.filter(m => m.sampled);
    const sampled = sampledMembers.length
      ? { rate: Math.max(...sampledMembers.map(m => m.sampled!.rate)), total: sampledMembers.reduce((a, m) => a + m.sampled!.total, 0), decoded: sampledMembers.reduce((a, m) => a + m.sampled!.decoded, 0) }
      : undefined;
    return {
      sampleId: GROUP_ID_BASE - g.id, sampleName: g.name || `Group ${g.id}`,
      coverage: sumCoverage(members.map(m => m.coverage)), junctions, spanning, sampled,
      loading: members.some(m => m.loading) || pending.length > 0,
      error: failed.length ? `${failed.map(m => m.sampleName).join(', ')}: ${failed[0].error}` : undefined,
      group: { id: g.id, n: g.sampleIds.length, loaded: members.length, agg: aggregateJunctions(junctions, tx, includeRetention ? spanning : undefined), samplesWith },
    };
  }), [groups, tracks, tx, runSamples, includeRetention]);

  const displayTracks = useMemo(() => {
    const withProfile = gtexTracks.map(t => {
      if (!t.gtex || t.gtex.lowCoverage) return { ...t, coverage: [], junctions: t.gtex?.lowCoverage ? [] : t.junctions };
      if (!tx) return t;
      const coverage: CoverageRun[] = [];
      let cursor = tx.exons.length ? tx.exons[0].start : 0;
      for (const ex of tx.exons) {
        const u = exonSiteUsage(ex, t.junctions).usage;
        if (ex.start > cursor) coverage.push({ start: cursor, end: ex.start, depth: 0 });
        coverage.push({ start: ex.start, end: ex.end, depth: u ?? 0 });
        cursor = ex.end;
      }
      return { ...t, coverage };
    });
    return [...(viewMode === 'groups' ? groupTracks : tracks), ...withProfile];
  }, [tracks, gtexTracks, tx, viewMode, groupTracks]);

  // ---- Gene navigation ----
  const navigateToGene = useCallback(async () => {
    const query = geneSearch.trim();
    if (!query) return;
    setGeneSearchLoading(true);
    setSearchError(null);
    try {
      const locus = parseLocus(query);
      if (locus) {
        // Coordinates: a position gets a 1 kb window, a range is shown as typed (capped to the zoom limit).
        const point = locus.start === locus.end;
        let s = point ? Math.max(0, locus.start - 1 - 500) : locus.start - 1, e = point ? locus.start + 500 : locus.end;
        if (e - s > MAX_VIEW_BP) { const c = (s + e) / 2; s = Math.max(0, Math.round(c - MAX_VIEW_BP / 2)); e = s + MAX_VIEW_BP; }
        let ax = axis;
        if (chromKey(locus.chrom) !== chromKey(currentChrom)) {
          // Another chromosome: the gene at the locus (coding first, then the largest overlap) becomes the queried gene
          const genes = await ds.getRegionGenes(locus.chrom, locus.start, locus.end);
          const ov = (g: GeneModel) => Math.min(g.end, locus.end) - Math.max(g.start, locus.start);
          const best = [...genes].sort((a, b) => Number(b.biotype === 'protein_coding') - Number(a.biotype === 'protein_coding') || ov(b) - ov(a))[0];
          setCurrentChrom(locus.chrom);
          setCurrentGeneId(undefined);
          if (best) {
            hintRef.current = { chrom: locus.chrom, start: best.start, end: best.end };
            const txData = await ds.getTranscript(best.gene_name, undefined, hintRef.current);
            const model = toTxModel(txData);
            setCurrentGeneName(txData.gene_name);
            setCurrentGeneStart(model.start); setCurrentGeneEnd(model.end);
            setTranscript(txData);
            ax = equalIntrons ? equalIntronAxis(model, intronWidth) : LINEAR_AXIS;
          } else {
            hintRef.current = undefined;
            setCurrentGeneName(`${locus.chrom}:${locus.start.toLocaleString()}`);
            setCurrentGeneStart(s); setCurrentGeneEnd(e);
            setTranscript(null); setTranscriptMissing('no RefSeq gene at this locus');
            ax = LINEAR_AXIS;
          }
        }
        setLocusMark({ chrom: locus.chrom, start: locus.start - 1, end: locus.end });
        setViewFromV(ax.toV(s), ax.toV(e), ax);
        setJunctionOffsets({});
        setReloadTrigger(n => n + 1);
        setGeneSearch('');
        setGeneSearchLoading(false);
        return;
      }
      const isEnsg = query.toUpperCase().startsWith('ENSG');
      const txData = await ds.getTranscript(query, isEnsg ? query : undefined);
      const model = toTxModel(txData);
      hintRef.current = { chrom: txData.chrom, start: txData.start, end: txData.end };
      setLocusMark(null);
      setCurrentGeneName(txData.gene_name);
      setCurrentGeneId(undefined);
      setCurrentChrom(txData.chrom);
      setCurrentGeneStart(model.start);
      setCurrentGeneEnd(model.end);
      setTranscript(txData);
      const ax = equalIntrons ? equalIntronAxis(model, intronWidth) : LINEAR_AXIS;
      const [s, e] = defaultView(model.start, model.end, ax);
      setViewStart(s); setViewEnd(e);
      setJunctionOffsets({});
      setReloadTrigger(n => n + 1);
      setGeneSearch('');
    } catch (e: any) {
      // not found: keep the current view, say why
      setSearchError(String(e?.message || e || 'not found').replace(/^Error:\s*/, ''));
    }
    setGeneSearchLoading(false);
  }, [geneSearch, equalIntrons, intronWidth, defaultView, axis, currentChrom, setViewFromV]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- SVG export (white background, full plot) ----
  const exportSvg = useCallback(() => {
    if (!svgRef.current) return;
    const svgEl = svgRef.current.cloneNode(true) as SVGSVGElement;
    svgEl.querySelectorAll('[data-export="skip"]').forEach(n => n.remove());
    svgEl.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const svgStr = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svgEl);
    const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `sashimi_${currentGeneName}_${currentChrom}_${viewStart + 1}-${viewEnd}${equalIntrons ? '_equal-introns' : ''}.svg`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }, [currentGeneName, currentChrom, viewStart, viewEnd, equalIntrons]);

  // ======================== Derived rendering data ========================

  const globalMaxDepth = useMemo(
    () => Math.max(1, ...tracks.map(t => maxDepthIn(t.coverage, viewStart, viewEnd))),
    [tracks, viewStart, viewEnd],
  );

  /** Junction keys seen (above threshold) in the comparison tracks, for "unique to primary" highlighting. */
  /** Tracks compared for the "unique to the first track" highlight: the samples, or the groups in aggregate view. */
  const comparedTracks = viewMode === 'groups' ? groupTracks : tracks;
  const otherTrackJunctionKeys = useMemo(() => {
    const keys = new Set<string>();
    for (let i = 1; i < comparedTracks.length; i++) {
      for (const j of comparedTracks[i].junctions) if (j.count >= minJunctionCount) keys.add(junctionKey(j));
    }
    return keys;
  }, [comparedTracks, minJunctionCount]);

  const plotRight = PLOT_LEFT + plotWidth;

  // ---- Reads track (memoised: thousands of nodes must not re-render on hover) ----
  /** Tooltip text for a variant site. */
  const siteLabel = useCallback((st: VariantSite) => {
    const what = st.kind === 'snv' ? `${st.ref}>${st.alt}` : st.kind === 'ins' ? `insertion ${st.alt} bp` : `deletion ${st.alt.slice(1)} bp`;
    const k = knownSnp(st);
    return `${currentChrom}:${(st.pos + 1).toLocaleString()} · ${what} · ${st.alt_count}/${st.depth} reads (${(st.vaf * 100).toFixed(0)}%)` +
      (k ? `\nknown common variant ${k.id} · max AF ${(k.maxAf * 100).toFixed(1)}%` : showSnps && visibleSnps.length ? '\nnot a common variant (dbSNP 155 common)' : '');
  }, [currentChrom, knownSnp, showSnps, visibleSnps.length]);

  type ReadsTrack = { height: number; el: JSX.Element; sites: VariantSite[] };
  /** One reads track per shown sample (Map in track order); each is drawn under its sample's coverage track. */
  const readsTracks = useMemo((): Map<number, ReadsTrack> => {
    const out = new Map<number, ReadsTrack>();
    if (!showReads) return out;
    const yOff = 0; // drawn in local coordinates; placed under its sample's coverage track with a translate
    const span = viewEnd - viewStart;
    const build = (sid: number): ReadsTrack => {
    const name = tracks.find(t => t.sampleId === sid)?.sampleName ?? '';
    const clipId = `sashimi-clip-reads-${sid}`;
    const frame = (h: number) => <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={h} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />;
    const header = (text: string, color = INK.muted) => (
      <text x={PLOT_LEFT + 8} y={yOff + 14} fontSize={10}>
        <tspan fill={INK.text} fontWeight={700}>{collapseReads ? 'Consensus reads' : 'Reads'}</tspan>
        <tspan fill={color}>{'  '}{name}{name ? ' · ' : ''}{text}</tspan>
      </text>
    );
    const message = (text: string, color?: string) => ({ height: 36, el: <g key={`reads${sid}`} fontFamily={FONT}>{frame(36)}{header(text, color)}</g>, sites: [] as VariantSite[] });

    if (span > READS_MAX_VIEW_BP) return message(`zoom in below ${formatBp(READS_MAX_VIEW_BP)} to load reads (window is ${formatBp(span)})`);
    const mode = collapseReads ? 'collapsed' : 'reads';
    const entry = readsData[sid];
    const current = entry && entry.fetched.chrom === currentChrom && entry.mode === mode ? entry.data : null;
    if (readsError[sid] && !current) return message(readsError[sid]!, UNIQUE_COLOR);
    if (!current) return message(collapseReads ? 'collapsing reads…' : 'loading reads…');

    const ref = current.reference;
    const sites: VariantSite[] = (current.sites || []).filter(st => st.pos >= viewStart && st.pos < viewEnd);
    const sitesRowH = 0; // stars are drawn in a strip above the sample's sashimi track, not here
    const basePx = (pos: number) => { const a = scale.x(pos), b = scale.x(pos + 1); return { left: Math.min(a, b), w: Math.max(1, Math.abs(b - a)) }; };

    // Reference segments with sequence: the whole window on a linear axis; only the exons in equal-intron
    // mode (compressed introns and flanks have no linear scale). Letters appear from 7 px per base.
    const segments: [number, number][] = [];
    if (ref) {
      const refEnd = ref.start + ref.seq.length;
      if (axis.kind === 'equal-intron' && tx) {
        for (const ex of tx.exons) {
          const a = Math.max(viewStart, ref.start, ex.start), b = Math.min(viewEnd, refEnd, ex.end);
          if (b > a) segments.push([a, b]);
        }
      } else {
        segments.push([Math.max(viewStart, ref.start), Math.min(viewEnd, refEnd)]);
      }
    }
    const pxPerBaseOf = (from: number) => Math.abs(scale.x(from + 1) - scale.x(from));
    const lettersShown = segments.some(([from, to]) => to - from <= 8000 && pxPerBaseOf(from) >= 7);
    // ---- Amino-acid row (MANE CDS translated codon by codon) above the reference bases, when letters are legible ----
    const aaRowH = ref && tx && tx.cdsStart != null && lettersShown ? AA_ROW_H : 0;
    // Minus-strand gene: the transcript-strand bases (complement, read 5′→3′ on the flipped axis) sit above the genomic + strand
    const revRowH = ref && reverse ? READS_SEQ_ROW_H : 0;
    const seqRowH = ref ? READS_SEQ_ROW_H : 0;
    const bodyTop = yOff + READS_HEADER_H + sitesRowH + aaRowH + revRowH + seqRowH;
    const aaEls: JSX.Element[] = [];
    if (aaRowH && ref && tx) {
      const rowY = yOff + READS_HEADER_H + sitesRowH;
      const seen = new Set<number>();
      for (const [from, to] of segments) {
        if (to - from > 8000 || pxPerBaseOf(from) < 7) continue;
        for (const c of codonsInWindow(tx, ref.seq, ref.start, from, to)) {
          if (seen.has(c.index)) continue;
          seen.add(c.index);
          const fill = c.aa === '*' ? AA_STOP_COLOR : c.index === 1 && c.aa === 'M' ? AA_START_COLOR : AA_FILLS[c.index % 2];
          const title = c.aa === '*' ? `stop codon (codon ${c.index})` : `${AA_NAMES[c.aa] ?? c.aa}${c.aa !== 'X' ? ` (${c.aa})` : ''} · codon ${c.index}${c.segments.length > 1 ? ' · spans a splice junction' : ''}`;
          const spans = c.segments.map(([gs, ge]) => { const a = scale.x(gs), b = scale.x(ge); return { left: Math.min(a, b), w: Math.abs(b - a) }; });
          const widest = spans.reduce((best, sp) => (sp.w > best.w ? sp : best), spans[0]);
          aaEls.push(
            <g key={`aa${c.index}`}>
              <title>{title}</title>
              {spans.map((sp, k) => <rect key={k} x={sp.left + 0.5} y={rowY + 2} width={Math.max(1, sp.w - 1)} height={aaRowH - 4} rx={2.5} fill={fill} opacity={c.aa === '*' || (c.index === 1 && c.aa === 'M') ? 0.9 : 0.75} />)}
              {widest.w >= 9 && <text x={widest.left + widest.w / 2} y={rowY + aaRowH - 4.5} textAnchor="middle" fill={c.aa === '*' || (c.index === 1 && c.aa === 'M') ? '#ffffff' : INK.text} fontSize={10} fontWeight={700}>{c.aa}</text>}
            </g>,
          );
        }
      }
    }

    // ---- Reference row: bars when ≥ 1 px per base, letters when ≥ 7 px. ----
    const refEls: JSX.Element[] = [];
    let refNote = '';
    if (ref) {
      const revY = yOff + READS_HEADER_H + sitesRowH + aaRowH;
      const rowY = revY + revRowH;
      let drawn = false;
      for (const [from, to] of segments) {
        if (to - from > 8000) continue;
        const pxPerBase = pxPerBaseOf(from);
        if (pxPerBase < 1) continue;
        drawn = true;
        for (let pos = from; pos < to; pos++) {
          const base = ref.seq[pos - ref.start];
          const { left, w } = basePx(pos);
          refEls.push(<rect key={`rb${pos}`} x={left} y={rowY + 2} width={Math.max(0.5, w - (w > 3 ? 0.5 : 0))} height={seqRowH - 4} fill={BASE_COLORS[base] || BASE_COLORS.N} opacity={pxPerBase >= 7 ? 0.22 : 0.85} />);
          if (pxPerBase >= 7) refEls.push(<text key={`rt${pos}`} x={left + w / 2} y={rowY + seqRowH - 5} textAnchor="middle" fill={BASE_COLORS[base] || BASE_COLORS.N} fontSize={Math.min(11, pxPerBase * 0.9)} fontWeight={700}>{base}</text>);
          if (revRowH) {
            const comp = COMPLEMENT[base] || 'N';
            refEls.push(<rect key={`vb${pos}`} x={left} y={revY + 2} width={Math.max(0.5, w - (w > 3 ? 0.5 : 0))} height={revRowH - 4} fill={BASE_COLORS[comp] || BASE_COLORS.N} opacity={pxPerBase >= 7 ? 0.22 : 0.85} />);
            if (pxPerBase >= 7) refEls.push(<text key={`vt${pos}`} x={left + w / 2} y={revY + revRowH - 5} textAnchor="middle" fill={BASE_COLORS[comp] || BASE_COLORS.N} fontSize={Math.min(11, pxPerBase * 0.9)} fontWeight={700}>{comp}</text>);
          }
        }
      }
      if (!drawn) refNote = axis.kind === 'equal-intron' ? 'reference · zoom in to see exon bases (introns are compressed)' : 'reference · zoom in to see bases';
    }

    const refSourceLabel: Record<string, string> = { fasta: 'REFERENCE_FASTA', ensembl: 'Ensembl (server)', browser: 'UCSC API (browser)' };
    const commonInfo = (ref ? ` · reference: ${refSourceLabel[current.reference_source ?? ''] ?? current.reference_source}` : ' · no reference genome (no REFERENCE_FASTA on the server, and the browser could not fetch bases from the UCSC / Ensembl APIs); mismatches only from MD tags') +
      (sites.length ? ` · ${sites.length} variant site${sites.length > 1 ? 's' : ''} ★` : '') + (readsLoading[sid] ? ' · updating…' : '');

    const wrap = (height: number, info: string, body: JSX.Element[], _bodyHeight: number) => ({
      height,
      sites,
      el: (
        <g key={`reads${sid}`} fontFamily={FONT}>
          <defs><clipPath id={clipId}><rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={height} /></clipPath></defs>
          {frame(height)}
          {header(info)}
          {aaRowH > 0 && <text transform={`translate(12, ${yOff + READS_HEADER_H + sitesRowH + aaRowH / 2}) rotate(-90)`} textAnchor="middle" fill={INK.faint} fontSize={8}>aa</text>}
          {revRowH > 0 && <text transform={`translate(12, ${yOff + READS_HEADER_H + sitesRowH + aaRowH + revRowH / 2}) rotate(-90)`} textAnchor="middle" fill={INK.faint} fontSize={8}><title>transcript strand (−): complement of the genomic bases, 5′→3′ left to right</title>ref −</text>}
          {ref && <text transform={`translate(12, ${yOff + READS_HEADER_H + sitesRowH + aaRowH + revRowH + seqRowH / 2}) rotate(-90)`} textAnchor="middle" fill={INK.faint} fontSize={8}>{revRowH > 0 ? 'ref +' : 'ref'}</text>}
          {refNote && <text x={PLOT_LEFT + 8} y={yOff + READS_HEADER_H + sitesRowH + aaRowH + revRowH + seqRowH - 5} fill={INK.faint} fontSize={9}>{refNote}</text>}
          <g clipPath={`url(#${clipId})`}>{aaEls}{refEls}{body}</g>
        </g>
      ),
    });

    // ======================= Collapsed mode: consensus rows =======================
    if (collapseReads) {
      // Consensus rows in display order (left to right on the plot) so that the groups of one phase
      // block sit next to each other; ambiguous buckets follow, the minor bucket last.
      const leftPx = (g: ReadGroup) => Math.min(...g.blocks.map(([bs, be]) => Math.min(scale.x(bs), scale.x(be))));
      const rank = { consensus: 0, ambiguous: 1, minor: 2 } as const;
      const groups: ReadGroup[] = [...(current.groups || [])].sort((a, b) => rank[a.kind] - rank[b.kind] || (a.kind === 'consensus' ? leftPx(a) - leftPx(b) : b.n - a.n));
      const allSites: VariantSite[] = current.sites || [];
      const hNum = (id: string) => parseInt(id.slice(1)) || 0;
      const ambLabel = (g: ReadGroup) => {
        const ids = [...(g.compatible || [])].sort((a, b) => hNum(a) - hNum(b));
        return ids.length > 3 ? `? · ${ids.length} groups` : `${ids.join('|')} ?`;
      };
      const rowStep = GROUP_ROW_H + 4;
      const bodyHeight = Math.max(1, groups.length) * rowStep + 6;
      const height = READS_HEADER_H + sitesRowH + aaRowH + revRowH + seqRowH + bodyHeight + 4;
      const rows = groups.map((g, gi) => {
        const top = bodyTop + 4 + gi * rowStep, mid = top + GROUP_ROW_H / 2;
        const isCons = g.kind === 'consensus';
        const fill = isCons ? READ_FILL : '#e2e5ea';
        const parts: JSX.Element[] = [];
        for (const [bs, be] of g.blocks) {
          const a = scale.x(bs), b = scale.x(be);
          parts.push(<rect key={`u${bs}`} x={Math.min(a, b)} y={top + 3} width={Math.max(1, Math.abs(b - a))} height={GROUP_ROW_H - 6} fill={fill} opacity={0.4} rx={1} />);
        }
        for (const [bs, be] of g.dense) {
          const a = scale.x(bs), b = scale.x(be);
          parts.push(<rect key={`d${bs}`} x={Math.min(a, b)} y={top} width={Math.max(1, Math.abs(b - a))} height={GROUP_ROW_H} fill={fill} rx={1.5}
            strokeDasharray={isCons ? undefined : '3 2'} stroke={isCons ? undefined : INK.faint} strokeWidth={isCons ? 0 : 0.8} />);
        }
        for (const [js, je] of g.chain) {
          parts.push(<line key={`j${js}`} x1={scale.x(js)} y1={mid} x2={scale.x(je)} y2={mid} stroke="#6b7280" strokeWidth={1.4} />);
        }
        if (!g.chain.length) {
          for (let k = 0; k + 1 < g.blocks.length; k++) {
            parts.push(<line key={`c${k}`} x1={scale.x(g.blocks[k][1])} y1={mid} x2={scale.x(g.blocks[k + 1][0])} y2={mid} stroke={INK.faint} strokeWidth={1} strokeDasharray="2 3" />);
          }
        }
        // Alleles at the variable sites
        const alleleTxt: string[] = [];
        allSites.forEach((st, si) => {
          const al = g.alleles[si];
          if (!al || al === '.') return;
          const { left, w } = basePx(st.pos);
          const ww = Math.max(w, 9), cx = left + w / 2;
          const isAlt = al !== 'ref';
          const letter = st.kind === 'snv' ? (al === 'alt' ? st.alt : al === 'ref' ? st.ref : al) : (isAlt ? st.alt : '=');
          const color = st.kind === 'snv' ? (BASE_COLORS[letter] || BASE_COLORS.N) : st.kind === 'ins' ? INSERTION_COLOR : '#111827';
          parts.push(
            <g key={`al${si}`}>
              <rect x={cx - ww / 2} y={top} width={ww} height={GROUP_ROW_H} fill={isAlt ? color : INK.bg} stroke={isAlt ? color : INK.faint} strokeWidth={isAlt ? 0 : 0.8} rx={2} />
              <text x={cx} y={top + GROUP_ROW_H - 4.5} textAnchor="middle" fill={isAlt ? '#fff' : INK.muted} fontSize={letter.length > 1 ? 7 : 9} fontWeight={700}>{letter}</text>
            </g>,
          );
          alleleTxt.push(`${currentChrom}:${(st.pos + 1).toLocaleString()} ${st.kind === 'snv' ? `${st.ref}>${letter}` : `${isAlt ? st.alt : 'ref'}`}${isAlt ? '' : ' (ref)'}`);
        });
        // Support badge just after the last block (or pinned inside the right edge)
        const xs = g.blocks.flatMap(([bs, be]) => [scale.x(bs), scale.x(be)]);
        const rightEnd = Math.min(plotRight - 4, Math.max(...xs) + 6);
        const badge = `${g.n.toLocaleString()} reads · ${(g.frac * 100).toFixed(g.frac < 0.1 ? 1 : 0)}%`;
        const bw = badge.length * 5.6 + 10;
        const bx = rightEnd + bw > plotRight - 2 ? plotRight - 2 - bw : rightEnd;
        parts.push(
          <g key="badge">
            <rect x={bx} y={mid - 7} width={bw} height={14} rx={7} fill={INK.bg} stroke={isCons ? '#6b7280' : INK.faint} strokeWidth={0.8} />
            <text x={bx + bw / 2} y={mid + 3.5} textAnchor="middle" fill={INK.text} fontSize={9} fontWeight={700}>{badge}</text>
          </g>,
        );
        // Row label at the left edge
        const label = g.kind === 'minor' ? `minor · ${g.patterns} pattern${g.patterns === 1 ? '' : 's'}` : g.kind === 'ambiguous' ? ambLabel(g) : g.id;
        const lw = label.length * 5.6 + 8;
        parts.push(
          <g key="label">
            <rect x={PLOT_LEFT + 3} y={mid - 7} width={lw} height={14} rx={3} fill={INK.bg} opacity={0.9} />
            <text x={PLOT_LEFT + 7} y={mid + 3.5} fill={isCons ? INK.text : INK.muted} fontSize={9} fontWeight={700}>{label}</text>
          </g>,
        );
        const chainTxt = g.chain.map(([js, je]) => junctionContext({ start: js, end: je, count: g.n }).info.label).join('; ');
        const title = (g.kind === 'consensus' ? `${g.id}: consensus of ${g.n.toLocaleString()} reads (${(g.frac * 100).toFixed(1)}%)${g.absorbed ? `, ${g.absorbed.toLocaleString()} absorbed as compatible` : ''}`
          : g.kind === 'ambiguous' ? `${g.n.toLocaleString()} reads compatible with ${[...(g.compatible || [])].sort((a, b) => hNum(a) - hNum(b)).join(', ')} (not informative between them)`
          : `${g.n.toLocaleString()} reads in ${g.patterns} pattern${g.patterns === 1 ? '' : 's'} below ${minJunctionCount} reads`) +
          (chainTxt ? `\nsplicing: ${chainTxt}` : '\nsplicing: none (unspliced)') +
          (alleleTxt.length ? `\nalleles: ${alleleTxt.join(', ')}` : '');
        return <g key={g.id}><title>{title}</title>{parts}</g>;
      });
      if (!groups.length) rows.push(<text key="none" x={PLOT_LEFT + 8} y={bodyTop + 14} fill={INK.muted} fontSize={10}>no reads in this window</text>);
      const nCons = groups.filter(g => g.kind === 'consensus').length;
      const info = `${current.total.toLocaleString()} reads → ${nCons} consensus group${nCons === 1 ? '' : 's'}` +
        (current.shown < current.total ? ` (from ${current.shown.toLocaleString()} sampled reads)` : '') + ` · min ${minJunctionCount} reads` + commonInfo;
      return wrap(height, info, rows, bodyHeight);
    }

    // ======================= Raw mode: packed alignments =======================
    const visible = current.reads.filter(r => r.e > viewStart && r.s < viewEnd);
    const { rows, nRows, hidden } = packReads(visible, READS_MAX_ROWS);
    const rowH = nRows > 60 ? 5 : 9; // squished rows beyond 60, like IGV's squished mode
    const bodyHeight = Math.max(1, nRows) * (rowH + 1) + 6;
    const height = READS_HEADER_H + sitesRowH + aaRowH + revRowH + seqRowH + bodyHeight + 4;
    const showLetters = rowH >= 9;
    const readsTop = bodyTop + 2;

    const modelBoundaries = boundariesOf(tx);
    const nSpan = visible.filter((r, i) => rows[i] >= 0 && readSpansBoundary(r, modelBoundaries)).length;
    const readEls = visible.map((r: AlignedRead, idx: number) => {
      const row = rows[idx];
      if (row < 0) return null;
      const top = readsTop + row * (rowH + 1);
      const mid = top + rowH / 2;
      const forward = r.r === 0;
      const tipX = scale.x(forward ? r.e : r.s);
      const parts: JSX.Element[] = [];
      const lowMapq = r.q === 0;
      const spans = readSpansBoundary(r, modelBoundaries);
      r.b.forEach(([bs, be], k) => {
        const xa = scale.x(bs), xb = scale.x(be);
        const left = Math.min(xa, xb), right = Math.max(xa, xb), w = Math.max(1, right - left);
        const isTipBlock = forward ? k === r.b.length - 1 : k === 0;
        if (isTipBlock && w > 8) {
          const tipRight = tipX >= right - 0.5;
          const d = tipRight
            ? `M${left},${top} H${right - 4} L${right},${mid} L${right - 4},${top + rowH} H${left} Z`
            : `M${right},${top} H${left + 4} L${left},${mid} L${left + 4},${top + rowH} H${right} Z`;
          parts.push(<path key={`b${k}`} d={d} fill={READ_FILL} opacity={lowMapq ? 0.35 : 1} />);
        } else {
          parts.push(<rect key={`b${k}`} x={left} y={top} width={w} height={rowH} fill={READ_FILL} opacity={lowMapq ? 0.35 : 1} />);
        }
        if (spans) parts.push(<rect key={`s${k}`} x={left} y={top} width={w} height={rowH} fill="none" stroke={RETENTION_COLOR} strokeWidth={1.2} />);
        if (k + 1 < r.b.length) {
          const gs = be, ge = r.b[k + 1][0];
          const isDel = r.d.some(dd => dd[0] === gs);
          const g1 = scale.x(gs), g2 = scale.x(ge);
          if (Math.abs(g2 - g1) > 0.5) parts.push(<line key={`g${k}`} x1={g1} y1={mid} x2={g2} y2={mid} stroke={isDel ? '#111827' : '#9ca3af'} strokeWidth={isDel ? 2 : 1} />);
        }
      });
      for (const [pos, base, qual] of r.m) {
        const { left, w } = basePx(pos);
        const color = BASE_COLORS[base] || BASE_COLORS.N;
        const alpha = qual < 10 ? 0.3 : qual < 20 ? 0.6 : 1;
        parts.push(<rect key={`m${pos}`} x={left} y={top} width={w} height={rowH} fill={color} opacity={alpha} />);
        if (showLetters && w >= 7) parts.push(<text key={`mt${pos}`} x={left + w / 2} y={top + rowH - 1.5} textAnchor="middle" fill="#fff" fontSize={Math.min(9, w)} fontWeight={700}>{base}</text>);
      }
      for (const [pos, len] of r.i) {
        const x = scale.x(pos);
        parts.push(
          <g key={`i${pos}`}>
            <title>{`insertion of ${len} bp at ${currentChrom}:${pos.toLocaleString()}`}</title>
            <rect x={x - 1} y={top} width={2} height={rowH} fill={INSERTION_COLOR} />
            <rect x={x - 2.5} y={top} width={5} height={1.5} fill={INSERTION_COLOR} />
            <rect x={x - 2.5} y={top + rowH - 1.5} width={5} height={1.5} fill={INSERTION_COLOR} />
          </g>,
        );
      }
      const title = `${r.n}\n${currentChrom}:${(r.s + 1).toLocaleString()}-${r.e.toLocaleString()} · ${forward ? '+' : '−'} strand · MAPQ ${r.q}${r.nh != null ? ` · NH ${r.nh}` : ''}\n` +
        `${r.b.length - 1 - r.d.length} splice gap${r.b.length - 1 - r.d.length === 1 ? '' : 's'} · ${r.m.length} mismatch${r.m.length === 1 ? '' : 'es'} · ${r.i.length} ins · ${r.d.length} del` +
        `${r.c[0] || r.c[1] ? ` · soft clips ${r.c[0]}/${r.c[1]}` : ''}` +
        (spans ? `\nruns unspliced through an exon–intron boundary of the model (≥ ${SPAN_EXON_ANCHOR} exonic and ≥ ${SPAN_INTRON_ANCHOR} intronic bases): counted for intron retention` : '');
      return <g key={r.n + r.s + r.f}><title>{title}</title>{parts}</g>;
    });

    const info = `${current.shown.toLocaleString()} of ${current.total.toLocaleString()} reads` +
      (current.shown < current.total ? ' (downsampled, zoom in for all)' : '') +
      (hidden ? ` · ${hidden.toLocaleString()} more not drawn (${READS_MAX_ROWS} rows max)` : '') +
      (modelBoundaries ? ` · ${nSpan.toLocaleString()} drawn read${nSpan === 1 ? '' : 's'} through an exon–intron boundary (teal outline)` : '') + commonInfo;
    return wrap(height, info, readEls.filter((e): e is JSX.Element => e !== null), bodyHeight);
    };
    for (const sid of readsSampleIds) out.set(sid, build(sid));
    return out;
  }, [showReads, readsSampleIds, collapseReads, minJunctionCount, viewStart, viewEnd, tracks, readsData, readsError, readsLoading, scale, plotWidth, currentChrom, tx, axis, junctionContext, reverse]);

  // ---- Screenshot to the basket (PNG + the viewer state and effect it documents) ----
  const takeSnapshot = useCallback(async () => {
    if (!svgRef.current || !onSnapshot || snapshotState === 'busy') return;
    setSnapshotState('busy');
    try {
      const png = await svgToPng(svgRef.current, 2);
      const readsName = readsAll ? 'all samples' : tracks.find(t => t.sampleId === effectiveReadsSampleId)?.sampleName;
      const context: SashimiSnapshotContext = {
        viewer: 'sashimi',
        gene: currentGeneName,
        transcript: tx?.transcriptId,
        region: { chrom: currentChrom, start: viewStart + 1, end: viewEnd },
        samples: displayTracks.map(t => t.sampleName),
        primarySample: tracks[0]?.sampleName ?? sampleName,
        options: {
          equalIntrons, allTranscripts: showAllTx, depthAxis, sharedY: depthAxis === 'shared', uniqueOnly,
          reads: showReads, readsSample: showReads ? readsName : undefined, collapsed: showReads && collapseReads,
          minJunctionReads: minJunctionCount, minVafPct, commonSnpsMinAf: showSnps ? snpMinAf : null,
          aggregate: viewMode === 'groups', arcLabels: showUsage ? 'usage' : 'reads', intronRetention: showUsage ? includeRetention : undefined,
        },
        groups: groups.length ? groups.map(g => ({ name: g.name, samples: g.sampleIds.map(id => runSamples.find(x => x.id === id)?.name ?? String(id)) })) : undefined,
        variantSites: readsTracks.get(readsSampleIds[0])?.sites.map(st => ({ pos: st.pos + 1, ref: st.ref, alt: st.alt, vaf: st.vaf, depth: st.depth })),
        knownVariants: showKnown ? primaryKnownHere.filter(v => v.end > viewStart && v.start < viewEnd).map(v => ({ kind: v.kind, chrom: v.chrom || currentChrom, start: v.start + 1, end: v.end, label: v.label, text: v.text })) : undefined,
        timestamp: new Date().toISOString(),
      };
      await onSnapshot(png, context);
      setSnapshotState('done');
    } catch (err) {
      console.error('Snapshot:', err);
      setSnapshotState('error');
    }
    setTimeout(() => setSnapshotState('idle'), 2500);
  }, [onSnapshot, snapshotState, tracks, effectiveReadsSampleId, readsAll, readsSampleIds, currentGeneName, tx, currentChrom, viewStart, viewEnd, sampleName,
    equalIntrons, showAllTx, depthAxis, uniqueOnly, showReads, collapseReads, minJunctionCount, minVafPct, readsTracks, showSnps, snpMinAf, displayTracks, showKnown, primaryKnownHere, viewMode, groups, runSamples]);

  interface ArcRender {
    j: JunctionArc; key: string; dragKey: string; level: number; color: string; dashed: boolean; unique: boolean;
    title: string; strokeW: number; geom: ReturnType<typeof arcGeom>;
    label: { x: number; y: number } | null;
    edge: { side: 'left' | 'right'; y: number; title: string } | null;
    offset: number;
    /** apex height above the higher of the two arc ends, in px */
    apexH: number;
    /** pill text: spliced reads, or the share at the intron for a group track */
    text: string;
    /** aggregate-view event of this junction (group tracks only) */
    agg?: AggEvent;
    /** Reading-frame consequence, for non-canonical junctions of a coding model. */
    frame: FrameInfo | null;
  }
  interface TrackLayout {
    track: TrackData; idx: number; color: string; yOff: number; juncH: number; yMax: number;
    paths: { fill: string; stroke: string }; arcs: ArcRender[]; height: number;
    /** Height of the variant-site strip at the top of the track (0 when none). */
    strip: number;
    /** Intron-retention pills on the baseline (usage mode). */
    retention: { x: number; y: number; text: string; title: string }[];
  }

  const transcriptY = RULER_H;
  const knownY = transcriptY + transcriptPanelH + TRACK_GAP;
  const snpY = knownY + (knownPanelH ? knownPanelH + TRACK_GAP : 0);
  const snpPanelH = showSnps ? SNP_PANEL_H : 0;
  const altY = snpY + (snpPanelH ? snpPanelH + TRACK_GAP : 0);
  const tracksTop = altY + (altPanelH ? altPanelH + TRACK_GAP : 0);

  /** Per-sample usage events (arc labels in %), same computation as the group tracks on the sample's own junctions. */
  const usageEvents = useMemo(() => {
    const m = new Map<number, AggResult>();
    if (arcLabel === 'usage') for (const t of tracks) m.set(t.sampleId, aggregateJunctions(t.junctions, tx, includeRetention ? t.spanning : undefined));
    return m;
  }, [arcLabel, tracks, tx, includeRetention]);

  const layouts: TrackLayout[] = useMemo(() => {
    let y = tracksTop;
    const out: TrackLayout[] = [];
    const plotRight = PLOT_LEFT + plotWidth;
    displayTracks.forEach((track, idx) => {
      const color = track.gtex ? track.gtex.tissue.color : TRACK_COLORS[idx % TRACK_COLORS.length];
      const trackAgg = track.group ? track.group.agg : !track.gtex ? usageEvents.get(track.sampleId) : undefined;
      const trackEvents = trackAgg?.events;
      const passes = (j: JunctionArc) => {
        if (track.gtex) return j.count >= 1;
        const ev = trackEvents?.get(junctionKey(j));
        if (ev && ev.shares.length) return Math.max(...ev.shares.map(sh => sh.pct)) * 100 >= minUsagePct;
        return j.count >= minJunctionCount;
      };
      const visible = track.junctions.filter(j => passes(j) && j.end > viewStart && j.start < viewEnd);
      const levels = layerJunctions(visible);
      const maxLevel = Math.max(1, ...levels.values());
      const readsBelow = readsTracks.get(track.sampleId);
      const strip = readsBelow && readsBelow.sites.length ? SITES_STRIP_H : 0;
      void maxLevel;
      // GTEx profiles are median reads per base and can sit well below 10: their axis floors at 0.1.
      // Relative mode draws each sample as a fraction of its own maximum in view (the axis reads 0–100 %).
      const ownMax = maxDepthIn(track.coverage, viewStart, viewEnd);
      const yMax = track.gtex ? 1 : (depthAxis === 'relative' || track.group) ? Math.max(1, ownMax) : niceMax(depthAxis === 'shared' ? globalMaxDepth : ownMax);
      // Arcs are first laid out against a baseline at y = 0 and measured; the junction area is then sized to
      // what they and their pills really occupy (rather than a fixed height per nesting level), and everything
      // is shifted down once the real baseline is known. This keeps the samples close together without any
      // arc ever reaching the sample-label band.
      const depthToY = (d: number) => -(Math.min(d, yMax) / yMax) * (COVERAGE_H - 12);

      const arcs: ArcRender[] = visible.map(j => {
        const key = junctionKey(j);
        const dragKey = `${track.sampleId}:${key}`;
        const { model, info, foreign } = junctionContext(j);
        const frame = model && info.cls !== 'canonical' ? junctionFrame(j, model, track.junctions) : null;
        const unique = idx === 0 && comparedTracks.length > 1 && !otherTrackJunctionKeys.has(key);
        const agg = trackEvents?.get(key);
        const share = agg?.shares[0];
        const approx = track.sampled ? '≈' : '';
        const text = agg ? (share ? pctLabel(share.pct) : `n=${approx}${j.count.toLocaleString()}`) : approx + j.count.toLocaleString();
        const x1 = scale.x(j.start), x2 = scale.x(j.end);
        const y1 = depthToY(depthAt(track.coverage, j.start - 1));
        const y2 = depthToY(depthAt(track.coverage, j.end));
        const level = levels.get(key) || 1;
        const apexH = 18 + (level - 1) * JUNC_LEVEL_STEP;
        const geom = arcGeom(x1, y1, x2, y2, apexH);
        const lo = Math.min(x1, x2), hi = Math.max(x1, x2);
        const visLo = Math.max(lo, PLOT_LEFT), visHi = Math.min(hi, plotRight);
        const offset = junctionOffsets[dragKey] || 0;
        const labelX = (visLo + visHi) / 2;
        const label = !track.gtex && visHi - visLo > 26 ? { x: labelX, y: arcYAtX(geom, labelX) } : null;   // GTEx arcs carry no number
        // Which genomic end is off-screen? Pixel-left is the genomic start unless the axis is flipped.
        let edge: ArcRender['edge'] = null;
        const partner = (side: 'left' | 'right') => {
          const atStart = (side === 'left') !== reverse;
          const pos = atStart ? j.start : j.end;
          const exon = atStart ? info.leftExon : info.rightExon;
          return `continues to ${currentChrom}:${(atStart ? pos + 1 : pos).toLocaleString()}${exon != null ? ` (exon ${exon})` : ''}`;
        };
        if (lo < PLOT_LEFT && hi > PLOT_LEFT) edge = { side: 'left', y: arcYAtX(geom, PLOT_LEFT), title: partner('left') };
        else if (hi > plotRight && lo < plotRight) edge = { side: 'right', y: arcYAtX(geom, plotRight), title: partner('right') };
        const inAlt = altJunctionIndex.get(key);
        const aggText = agg
          ? (agg.shares.length
            ? agg.shares.map(sh => `${pctLabel(sh.pct)} ${sh.note}`).join('\n')
            : 'touches no annotated splice site: no share') +
            `\n${AGG_CLASS_LABEL[agg.cls]}${agg.partner ? ' (two arcs paired)' : ''} · ${j.count.toLocaleString()} ${track.group ? `pooled reads in ${track.group.samplesWith.get(key) ?? 0}/${track.group.loaded} samples` : `spliced read${j.count > 1 ? 's' : ''}`}\n`
          : null;
        const title = (track.gtex ? `median ${j.count.toLocaleString()} junction reads per sample (${track.sampleName})\n` : aggText ?? `${j.count.toLocaleString()} spliced read${j.count > 1 ? 's' : ''}\n`) +
          `${currentChrom}:${(j.start + 1).toLocaleString()}-${j.end.toLocaleString()} · intron ${formatBp(j.end - j.start)}\n` +
          info.label + (foreign ? ` (${foreign.strand === tx?.strand ? 'same strand as' : 'antisense to'} ${tx?.geneName ?? 'the queried gene'})` : '') +
          (inAlt ? `\nannotated in ${inAlt.slice(0, 4).join(', ')}${inAlt.length > 4 ? ` +${inAlt.length - 4}` : ''}` : '') +
          (frame ? `\nreading frame: ${frameLabel(frame)} · ${frame.text}` : '') +
          (unique ? `\nnot seen in the comparison ${track.group ? 'group' : 'sample'}${comparedTracks.length > 2 ? 's' : ''}` : '') +
          (track.sampled ? `\n≈ deep window: 1 read in ${track.sampled.rate} decoded${track.group ? ' in at least one sample' : ''}, counts scaled back (estimates)` : '');
        return {
          j, key, dragKey, level, color: unique ? UNIQUE_COLOR : agg?.cls === 'pseudo_exon' ? PSEUDO_EXON_COLOR : color, dashed: info.cls !== 'canonical', unique, title,
          strokeW: agg ? 1 + 3.5 * (share?.pct ?? 0) : Math.min(4.5, 1 + Math.log2(j.count) * 0.55), geom, label, edge, offset, frame, apexH, text, agg,
        };
      });
      // Push colliding read-count pills upward (lower arcs keep their place) so every count stays legible,
      // notably when several arcs leave the window and share the same visible midpoint.
      const placed: { x: number; y: number; w: number }[] = [];
      for (const a of [...arcs].sort((p, q) => p.level - q.level)) {
        if (!a.label) continue;
        const w = a.text.length * 6 + 10 + (a.frame && a.frame.frame !== 'unknown' ? FRAME_GLYPH_R * 2 + 6 : 0);
        for (let iter = 0; iter < 24; iter++) {
          const hit = placed.some(p => Math.abs(p.x - a.label!.x) < (p.w + w) / 2 + 4 && Math.abs(p.y - a.label!.y) < LABEL_H + 2);
          if (!hit) break;
          a.label.y -= LABEL_H + 3;
        }
        placed.push({ x: a.label.x, y: a.label.y, w });
      }
      // Highest point of the track content relative to the baseline (negative): the coverage area itself,
      // every arc apex, read-count pill and edge chevron, including the offset of arcs the user dragged.
      let top = -COVERAGE_H;
      for (const a of arcs) {
        top = Math.min(top, Math.min(a.geom.y1, a.geom.y2) - a.apexH + a.offset);
        if (a.label) top = Math.min(top, a.label.y + a.offset - LABEL_H / 2);
        if (a.edge) top = Math.min(top, a.edge.y + a.offset - 4);
      }
      const juncH = TRACK_LABEL_H + strip + Math.max(JUNC_MIN_H, -top - COVERAGE_H + JUNC_PAD);
      const baseline = y + juncH + COVERAGE_H;
      for (const a of arcs) {
        a.geom = arcGeom(a.geom.x1, a.geom.y1 + baseline, a.geom.x2, a.geom.y2 + baseline, a.apexH);
        if (a.label) a.label.y += baseline;
        if (a.edge) a.edge.y += baseline;
      }
      const paths = buildCoveragePaths(track.coverage, scale, viewStart, viewEnd, baseline, d => baseline + depthToY(d));
      // intron retention pills: on the baseline at the middle of the visible part of each intron, above the Min % threshold
      const retention = (trackAgg?.retention ?? [])
        .filter(r => r.pct * 100 >= minUsagePct && r.pct > 0 && r.end > viewStart && r.start < viewEnd)
        .map(r => {
          const xa = scale.x(Math.max(r.start, viewStart)), xb = scale.x(Math.min(r.end, viewEnd));
          return {
            x: (xa + xb) / 2, y: baseline - LABEL_H / 2 - 3, text: `IR ${pctLabel(r.pct)}`,
            title: `${pctLabel(r.pct)} ${r.note}\n${track.group ? 'reads pooled over the group' : track.sampleName} · unspliced through both boundaries, ≥ 6 aligned bases on the exon side and ≥ 10 on the intron side`,
          };
        });
      const height = juncH + COVERAGE_H;
      out.push({ track, idx, color, yOff: y, juncH, yMax, paths, arcs, height, strip, retention });
      y += height + SASHIMI_GAP;
      if (readsBelow) y += readsBelow.height + TRACK_GAP;
    });
    return out;
  }, [comparedTracks, displayTracks, minJunctionCount, viewStart, viewEnd, scale, depthAxis, globalMaxDepth, tx, otherTrackJunctionKeys, junctionOffsets, plotWidth, reverse, currentChrom, readsTracks, tracksTop, altJunctionIndex, junctionContext, usageEvents, minUsagePct]);

  const lastTrackBottom = layouts.length ? layouts[layouts.length - 1].yOff + layouts[layouts.length - 1].height + TRACK_GAP : tracksTop;
  /** Each reads track sits right under the coverage track of its sample. */
  const readsPlacements = [...readsTracks].map(([sid, rt]) => {
    const L = layouts.find(l => l.track.sampleId === sid);
    return { sid, rt, y: L ? L.yOff + L.height + TRACK_GAP : lastTrackBottom };
  });
  // A reads track may sit under the last coverage track, so it can extend the stack
  const tracksBottom = readsPlacements.reduce((m, p) => Math.max(m, p.y + p.rt.height + TRACK_GAP), lastTrackBottom);


  /** Legend items (local y), laid out left to right and wrapped into rows that fit the plot width. */
  const legendItems = (() => {
    const y = 14;
    const items: { w: number; el: JSX.Element }[] = [];
    const primaryColor = TRACK_COLORS[0];
    const line = (color: string, dashed: boolean, label: string, key: string) => items.push({
      w: label.length * 5.6 + 40,
      el: (
        <g key={key}>
          <line x1={0} y1={y} x2={22} y2={y} stroke={color} strokeWidth={2} strokeDasharray={dashed ? '5 3.5' : undefined} strokeLinecap="round" />
          <text x={28} y={y + 3.5} fill={INK.muted} fontSize={9.5}>{label}</text>
        </g>
      ),
    });
    line(primaryColor, false, 'canonical junction (consecutive exons)', 'l1');
    line(primaryColor, true, 'non-canonical (exon skipping, novel site)', 'l2');
    if (showUsage) line(PSEUDO_EXON_COLOR, true, `pseudo-exon (alt 3′ in + alt 5′ out, ≤ ${PSEUDO_EXON_MAX_BP} bp, paired)`, 'l2b');
    if (showUsage && includeRetention) items.push({
      w: 300, el: (
        <g key="lir">
          <rect x={0} y={y - 6.5} width={32} height={13} rx={6.5} fill={INK.bg} stroke={RETENTION_COLOR} strokeWidth={1} />
          <text x={16} y={y + 3} textAnchor="middle" fill={RETENTION_COLOR} fontSize={8.5} fontWeight={700}>IR %</text>
          <text x={38} y={y + 3.5} fill={INK.muted} fontSize={9.5}>intron retention: unspliced reads through both boundaries, share of the intron's reads</text>
        </g>
      ),
    });
    items.push({
      w: 150, el: (
        <g key="lf">
          {renderFrameGlyph(FRAME_GLYPH_R, y, { frame: 'in', delta: 0, cdsBases: 0, text: '' }, 'lf1')}
          <text x={FRAME_GLYPH_R * 2 + 4} y={y + 3.5} fill={INK.muted} fontSize={9.5}>in frame</text>
          {renderFrameGlyph(FRAME_GLYPH_R * 2 + 58, y, { frame: 'out', delta: 0, cdsBases: 0, text: '' }, 'lf2')}
          <text x={FRAME_GLYPH_R * 3 + 62} y={y + 3.5} fill={INK.muted} fontSize={9.5}>frameshift</text>
        </g>
      ),
    });
    if (neighbourRows.length) {
      const g = (color: string, label: string, key: string) => items.push({
        w: label.length * 5.6 + 40,
        el: (
          <g key={key}>
            <line x1={0} y1={y} x2={22} y2={y} stroke={color} strokeWidth={1.2} />
            <rect x={4} y={y - 5} width={6} height={10} fill={color} rx={1} />
            <rect x={14} y={y - 5} width={5} height={10} fill={color} rx={1} />
            <text x={28} y={y + 3.5} fill={INK.muted} fontSize={9.5}>{label}</text>
          </g>
        ),
      });
      g(SAME_SENSE_COLOR, 'neighbouring gene, same strand', 'ln1');
      g(ANTISENSE_COLOR, 'neighbouring gene, antisense', 'ln2');
    }
    if (comparedTracks.length > 1) line(UNIQUE_COLOR, false, `only in ${comparedTracks[0].sampleName}`, 'l3');
    items.push({
      w: 60, el: (
        <g key="l4">
          <rect x={0} y={y - 6} width={14} height={12} fill={INK.exon} rx={1.5} />
          <text x={19} y={y + 3.5} fill={INK.muted} fontSize={9.5}>CDS</text>
        </g>
      ),
    });
    items.push({
      w: 60, el: (
        <g key="l5">
          <rect x={0} y={y - 3} width={14} height={6} fill={INK.utr} rx={1.5} />
          <text x={19} y={y + 3.5} fill={INK.muted} fontSize={9.5}>UTR</text>
        </g>
      ),
    });
    if (showReads) {
      items.push({
        w: 118, el: (
          <g key="l7">
            {(['A', 'C', 'G', 'T'] as const).map((b, i) => (
              <g key={b}>
                <rect x={i * 14} y={y - 6} width={12} height={12} fill={BASE_COLORS[b]} rx={1.5} />
                <text x={i * 14 + 6} y={y + 3.5} textAnchor="middle" fill="#fff" fontSize={8.5} fontWeight={700}>{b}</text>
              </g>
            ))}
            <text x={60} y={y + 3.5} fill={INK.muted} fontSize={9.5}>mismatch</text>
          </g>
        ),
      });
      items.push({
        w: 76, el: (
          <g key="l8">
            <rect x={0} y={y - 5} width={2} height={10} fill={INSERTION_COLOR} />
            <rect x={-1.5} y={y - 5} width={5} height={1.5} fill={INSERTION_COLOR} />
            <rect x={-1.5} y={y + 3.5} width={5} height={1.5} fill={INSERTION_COLOR} />
            <text x={8} y={y + 3.5} fill={INK.muted} fontSize={9.5}>insertion</text>
          </g>
        ),
      });
      items.push({
        w: 72, el: (
          <g key="l9">
            <line x1={0} y1={y} x2={14} y2={y} stroke="#111827" strokeWidth={2} />
            <text x={19} y={y + 3.5} fill={INK.muted} fontSize={9.5}>deletion</text>
          </g>
        ),
      });
      items.push({
        w: 90, el: (
          <g key="l10">
            <path d={starPath(6, y, 6.5)} fill={STAR_COLOR} stroke="#92400e" strokeWidth={0.8} />
            <text x={16} y={y + 3.5} fill={INK.muted} fontSize={9.5}>variant site</text>
          </g>
        ),
      });
    }
    if (showSnps) {
      items.push({
        w: 176, el: (
          <g key="l11">
            <line x1={6} y1={y + 6} x2={6} y2={y - 4} stroke={SNP_SNV_COLOR} strokeWidth={1} /><circle cx={6} cy={y - 4} r={2.6} fill={SNP_SNV_COLOR} />
            <line x1={22} y1={y + 6} x2={22} y2={y - 4} stroke={SNP_INDEL_COLOR} strokeWidth={1} /><rect x={19} y={y - 6.5} width={6} height={5} rx={1} fill={SNP_INDEL_COLOR} />
            <text x={31} y={y + 3.5} fill={INK.muted} fontSize={9.5}>common SNV / indel (height ∝ AF)</text>
          </g>
        ),
      });
      if (showReads) items.push({
        w: 118, el: (
          <g key="l12">
            <circle cx={7} cy={y} r={8} fill="none" stroke={SNP_KNOWN_RING} strokeWidth={1.4} />
            <path d={starPath(7, y, 5.5)} fill={STAR_COLOR} stroke="#92400e" strokeWidth={0.7} />
            <text x={20} y={y + 3.5} fill={INK.muted} fontSize={9.5}>known common SNP</text>
          </g>
        ),
      });
    }
    if (knownPanelH > 0) {
      items.push({
        w: 'known variant of the sample: SNV / indel, CNV / SV band'.length * 5.6 + 48, el: (
          <g key="l13">
            <path d={`M6,${y - 5} L11,${y} L6,${y + 5} L1,${y} Z`} fill={KNOWN_VARIANT_COLORS.snv} />
            <rect x={17} y={y - 4} width={16} height={8} rx={1.5} fill={KNOWN_VARIANT_COLORS.del} opacity={0.55} stroke={KNOWN_VARIANT_COLORS.del} strokeWidth={0.8} />
            <text x={39} y={y + 3.5} fill={INK.muted} fontSize={9.5}>known variant of the sample: SNV / indel, CNV / SV band</text>
          </g>
        ),
      });
    }
    items.push({
      w: 0, el: <text key="l6" x={0} y={y + 3.5} fill={INK.faint} fontSize={9}>{viewMode === 'groups' ? 'arc width ∝ usage · label = % of the reads competing at the intron, 100 % per intron (reads pooled over the group)' : showUsage ? 'arc width ∝ usage · label = % of the reads competing at the intron, 100 % per intron (sample reads)' : 'arc width ∝ log₂ reads · label = spliced reads'}</text>,
    });
    return items;
  })();
  const legendPlaced = (() => {
    let x = PLOT_LEFT, row = 0;
    const out: { it: { w: number; el: JSX.Element }; x: number; row: number }[] = [];
    for (const it of legendItems) {
      if (x + it.w > PLOT_LEFT + plotWidth && x > PLOT_LEFT) { x = PLOT_LEFT; row++; }
      out.push({ it, x, row });
      x += it.w;
    }
    return { placed: out, rows: row + 1 };
  })();
  const legendH = legendPlaced.rows * LEGEND_ROW_H + 8;
  const renderLegend = (yOff: number) => (
    <g fontFamily={FONT}>
      {legendPlaced.placed.map(({ it, x, row }) => <g key={it.el.key} transform={`translate(${x}, ${yOff + row * LEGEND_ROW_H})`}>{it.el}</g>)}
    </g>
  );

  const legendY = tracksBottom;
  const totalHeight = legendY + legendH;

  // ======================== Render helpers ========================

  const renderRuler = () => {
    const baseY = RULER_H - 1;
    const items: JSX.Element[] = [];
    const labelY = baseY - 10;
    const tick = (x: number, label: string | null, key: string, anchor: 'middle' | 'start' | 'end' = 'middle') => (
      <g key={key}>
        <line x1={x} y1={baseY - 6} x2={x} y2={baseY} stroke={INK.gridStrong} strokeWidth={1} />
        {label && <text x={x} y={labelY} textAnchor={anchor} fill={INK.muted} fontSize={9.5}>{label}</text>}
      </g>
    );
    if (axis.kind === 'linear') {
      for (const t of niceTicks(viewStart + 1, viewEnd, Math.max(3, Math.floor(plotWidth / 130)))) {
        items.push(tick(scale.x(t - 0.5), t.toLocaleString(), `t${t}`));
      }
    } else if (tx) {
      // Exon boundaries carry the coordinates; introns get a broken-axis mark.
      let lastLabelX = -Infinity;
      const bounds: { pos: number; label: number }[] = [];
      for (const ex of tx.exons) bounds.push({ pos: ex.start, label: ex.start + 1 }, { pos: ex.end, label: ex.end });
      bounds.sort((a, b) => scale.x(a.pos) - scale.x(b.pos));
      for (const b of bounds) {
        const x = scale.x(b.pos);
        if (x < PLOT_LEFT || x > plotRight) continue;
        const show = Math.abs(x - lastLabelX) >= 68;
        if (show) lastLabelX = x;
        items.push(tick(x, show ? b.label.toLocaleString() : null, `b${b.pos}`));
      }
      for (const intron of intronsOf(tx)) {
        const mx = scale.x((intron.start + intron.end) / 2);
        if (mx < PLOT_LEFT + 6 || mx > plotRight - 6) continue;
        items.push(
          <g key={`brk${intron.start}`} stroke={INK.gridStrong} strokeWidth={1.2} fill="none">
            <path d={`M${mx - 4},${baseY + 3} l3,-6 M${mx + 1},${baseY + 3} l3,-6`} />
          </g>,
        );
      }
    }
    return (
      <g fontFamily={FONT}>
        <line x1={PLOT_LEFT} y1={baseY} x2={plotRight} y2={baseY} stroke={INK.gridStrong} strokeWidth={1} />
        {items}
        <text x={PLOT_LEFT} y={11} fill={INK.faint} fontSize={9} fontWeight={600} letterSpacing={0.4}>
          {currentChrom.toUpperCase()}{reverse ? '  ·  5′ → 3′ (reverse strand, axis flipped)' : ''}{axis.kind === 'equal-intron' ? '  ·  INTRONS DRAWN AT EQUAL WIDTH' : ''}
        </text>
      </g>
    );
  };

  const renderTrack = (L: TrackLayout) => {
    const { track, idx, color, yOff, juncH, yMax, paths, arcs, height } = L;
    const clipId = `sashimi-clip-${track.sampleId}`;
    void L.strip;
    const baseline = yOff + juncH + COVERAGE_H;
    const depthToY = (d: number) => baseline - (Math.min(d, yMax) / yMax) * (COVERAGE_H - 12);
    const isPrimary = idx === 0 && tracks.length > 1;
    const gtexNote = track.gtex ? `  GTEx ${track.gtex.dataset.replace('gtex_', '')} · n=${track.gtex.tissue.samples}${track.gtex.tpm != null ? ` · median ${track.gtex.tpm < 10 ? track.gtex.tpm.toFixed(2) : track.gtex.tpm.toFixed(0)} TPM` : ''}${track.gtex.lowCoverage ? ' · LOW COVERAGE (TPM < 1)' : ' · exon usage from junction medians'}` : '';
    const gtexWarn = track.gtex ? (track.error || track.gtex.warning || '') : '';
    const relative = (depthAxis === 'relative' || !!track.group) && !track.gtex;
    const groupNote = track.group ? `  ·  ${track.group.loaded}/${track.group.n} sample${track.group.n === 1 ? '' : 's'} pooled` : '';
    const axisNote = relative ? `  ·  max ${yMax.toLocaleString()}×` : '';
    const sampledNote = track.sampled ? `  ·  ≈ 1 read in ${track.sampled.rate}` : '';
    const sampledTitle = track.sampled
      ? `Deep window: ${track.sampled.decoded.toLocaleString()} of ${track.sampled.total.toLocaleString()} reads decoded (every ${track.sampled.rate === 2 ? 'other' : `${track.sampled.rate}th`} read${track.group ? ', in the deepest sample' : ''}); depths and counts are scaled back by ${track.sampled.rate} and are estimates. Zoom in for exact counts.`
      : '';
    const labelW = track.sampleName.length * 6.4 + 24 + (isPrimary ? 44 : 0) + gtexNote.length * 5.2 + groupNote.length * 5.2 + axisNote.length * 5.2 + sampledNote.length * 5.2;
    const status = track.error && track.coverage.length === 0
      ? { text: track.error, color: UNIQUE_COLOR }
      : track.loading ? { text: track.coverage.length ? 'updating…' : 'loading…', color: INK.faint } : null;

    return (
      <g key={track.sampleId} fontFamily={FONT}>
        <defs>
          <clipPath id={clipId}><rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={height} /></clipPath>
        </defs>
        {/* Track frame + subtle junction/coverage separator */}
        <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={height} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />
        {L.strip > 0 && <line x1={PLOT_LEFT} y1={yOff + TRACK_LABEL_H + L.strip} x2={plotRight} y2={yOff + TRACK_LABEL_H + L.strip} stroke={INK.grid} strokeWidth={0.8} strokeDasharray="2 3" />}
        <line x1={PLOT_LEFT} y1={baseline} x2={plotRight} y2={baseline} stroke={INK.gridStrong} strokeWidth={1} />

        {/* Depth axis */}
        {[0, yMax / 2, yMax].map(v => {
          const y = depthToY(v);
          return (
            <g key={v}>
              <line x1={PLOT_LEFT - 4} y1={y} x2={PLOT_LEFT} y2={y} stroke={INK.gridStrong} strokeWidth={1} />
              {v > 0 && <line x1={PLOT_LEFT} y1={y} x2={plotRight} y2={y} stroke={INK.grid} strokeWidth={0.6} strokeDasharray="2 4" />}
              <text x={PLOT_LEFT - 7} y={y + 3} textAnchor="end" fill={INK.muted} fontSize={9}>{relative ? `${Math.round((v / yMax) * 100)} %` : v < 10 && v % 1 ? v.toFixed(v < 1 ? 2 : 1) : v.toLocaleString()}</text>
            </g>
          );
        })}
        <text transform={`translate(${12}, ${baseline - (COVERAGE_H - 12) / 2}) rotate(-90)`} textAnchor="middle" fill={INK.faint} fontSize={8.5} letterSpacing={0.3}>{track.gtex ? 'exon usage' : relative ? 'depth · % of max' : 'depth'}</text>
        {track.gtex?.lowCoverage && <text x={PLOT_LEFT + plotWidth / 2} y={baseline - COVERAGE_H / 2 + 4} textAnchor="middle" fill={INK.faint} fontSize={13} fontWeight={600}>low coverage · median {track.gtex.tpm?.toFixed(2)} TPM in {track.sampleName}</text>}

        <g clipPath={`url(#${clipId})`}>
          {paths.fill && <path d={paths.fill} fill={withAlpha(color, 0.26)} stroke="none" />}
          {paths.stroke && <path d={paths.stroke} fill="none" stroke={color} strokeWidth={1.3} strokeLinejoin="round" />}

          {/* Junction arcs */}
          {arcs.map(a => (
            <g key={a.key} transform={a.offset ? `translate(0, ${a.offset})` : undefined}
              style={{ cursor: junctionDrag.current?.key === a.dragKey ? 'grabbing' : 'grab' }}
              onMouseDown={e => {
                e.stopPropagation();
                dragMoved.current = false;
                junctionDrag.current = { key: a.dragKey, startY: e.clientY, startOffset: a.offset };
              }}
              onClick={e => {
                e.stopPropagation();
                if (dragMoved.current) return;
                const p = svgPoint(e);
                setPopover(prev => (prev?.kind === 'junction' && prev.key === a.key ? null : { kind: 'junction', key: a.key, j: a.j, x: p.x, y: p.y }));
              }}>
              <title>{a.title}</title>
              <path d={a.geom.d} fill="none" stroke="transparent" strokeWidth={Math.max(a.strokeW + 8, 12)} />
              <path d={a.geom.d} fill="none" stroke={INK.bg} strokeWidth={a.strokeW + 2} opacity={0.9} />
              <path d={a.geom.d} fill="none" stroke={a.color} strokeWidth={a.strokeW} strokeLinecap="round"
                strokeDasharray={a.dashed ? '5 3.5' : undefined} opacity={0.92} />
            </g>
          ))}
          {/* Allele-fraction bars at the variant sites: alt allele in its base colour over the reference share */}
          {L.strip > 0 && readsTracks.get(track.sampleId)?.sites.map(st => {
            const xa = scale.x(st.pos), xb = scale.x(st.pos + 1);
            let left = Math.min(xa, xb), w = Math.abs(xb - xa);
            if (w < 3) { left += w / 2 - 1.5; w = 3; }
            const d = Math.max(depthAt(track.coverage, st.pos), st.depth);
            const top = depthToY(d), full = baseline - top, altH = full * st.vaf;
            const colr = st.kind === 'snv' ? (BASE_COLORS[st.alt] || BASE_COLORS.N) : st.kind === 'ins' ? INSERTION_COLOR : '#111827';
            const pct = `${(st.vaf * 100).toFixed(st.vaf < 0.1 ? 1 : 0)}%`;
            return (
              <g key={`vaf${st.pos}${st.kind}`}>
                <title>{siteLabel(st)}</title>
                <rect x={left - 1} y={top - 1} width={w + 2} height={full + 2} fill={INK.bg} opacity={0.9} />
                <rect x={left} y={top} width={w} height={Math.max(0, full - altH)} fill="#9ca3af" />
                <rect x={left} y={baseline - altH} width={w} height={altH} fill={colr} />
                <text x={left + w / 2} y={top - 5} textAnchor="middle" fontSize={9} fontWeight={700} fill={colr} stroke={INK.bg} strokeWidth={3} paintOrder="stroke">{pct}</text>
              </g>
            );
          })}
          {/* Read-count pills, drawn after every arc so no stroke paints over a number */}
          {arcs.filter(a => a.label).map(a => {
            const txt = a.text;
            const w = txt.length * 6 + 10;
            const lx = a.label!.x, ly = a.label!.y + a.offset;
            const glyph = a.frame && a.frame.frame !== 'unknown' ? a.frame : null;
            return (
              <g key={`l-${a.key}`} pointerEvents="none">
                <rect x={lx - w / 2} y={ly - LABEL_H / 2} width={w} height={LABEL_H} rx={LABEL_H / 2} fill={INK.bg} stroke={a.color} strokeWidth={1} />
                <text x={lx} y={ly + 3.5} textAnchor="middle" fill={INK.text} fontSize={9.5} fontWeight={700}>{txt}</text>
                {glyph && renderFrameGlyph(lx + w / 2 + FRAME_GLYPH_R + 3, ly, glyph, `fg-${a.key}`)}
              </g>
            );
          })}
          {/* Intron-retention pills on the baseline (usage mode) */}
          {L.retention.map((r, i) => {
            const w = r.text.length * 6 + 10;
            return (
              <g key={`ir-${i}`}>
                <title>{r.title}</title>
                <rect x={r.x - w / 2} y={r.y - LABEL_H / 2} width={w} height={LABEL_H} rx={LABEL_H / 2} fill={INK.bg} stroke={RETENTION_COLOR} strokeWidth={1} />
                <text x={r.x} y={r.y + 3.5} textAnchor="middle" fill={RETENTION_COLOR} fontSize={9.5} fontWeight={700}>{r.text}</text>
              </g>
            );
          })}
        </g>

        {/* Edge chevrons for arcs continuing beyond the window (drawn outside the clip) */}
        {arcs.filter(a => a.edge).map(a => {
          const e = a.edge!;
          const y = e.y + a.offset;
          const d = e.side === 'left'
            ? `M${PLOT_LEFT - 1},${y} l-6,-4 v8 z`
            : `M${plotRight + 1},${y} l6,-4 v8 z`;
          return (
            <g key={`e-${a.key}`}>
              <title>{`${a.j.count.toLocaleString()} reads · ${e.title}`}</title>
              <path d={d} fill={a.color} />
            </g>
          );
        })}

        {/* Sample label, in the band reserved above the junctions so it never sits on arcs or coverage */}
        <g transform={`translate(${PLOT_LEFT + 8}, ${yOff + 2})`}>
          <rect x={-4} y={-2} width={labelW} height={18} rx={4} fill={INK.bg} opacity={0.96} />
          <rect x={0} y={2} width={10} height={10} rx={2} fill={color} />
          <text x={15} y={11.5} fontSize={10.5}>
            <tspan fill={INK.text} fontWeight={600}>{track.sampleName}</tspan>
            {isPrimary && <tspan fill={INK.faint} fontSize={9}>{'  primary'}</tspan>}
            {gtexNote && <tspan fill={INK.faint} fontSize={9}>{gtexNote}</tspan>}
            {groupNote && <tspan fill={INK.faint} fontSize={9}>{groupNote}</tspan>}
            {axisNote && <tspan fill={INK.faint} fontSize={9}>{axisNote}</tspan>}
            {sampledNote && <tspan fill={SNP_INDEL_COLOR} fontSize={9} fontWeight={600}>{sampledNote}<title>{sampledTitle}</title></tspan>}

          </text>
          {/* Make primary chip (standalone): promote this sample to the first track */}
          {allowPrimarySwitch && !track.gtex && !track.group && idx > 0 && (
            <g data-export="skip" transform={`translate(${labelW + 52}, 0)`} style={{ cursor: 'pointer' }}
              onClick={e => { e.stopPropagation(); setPrimary(track.sampleId); onPrimaryChange?.(track.sampleId); }} onMouseDown={e => e.stopPropagation()}>
              <title>{`Make ${track.sampleName} the primary sample (first track; "unique" junctions are judged against the others)`}</title>
              <rect x={0} y={0} width={76} height={14} rx={7} fill={INK.bg} stroke={INK.faint} strokeWidth={0.8} />
              <text x={38} y={10} textAnchor="middle" fill={INK.muted} fontSize={8.5} fontWeight={600}>★ make primary</text>
            </g>
          )}
          {/* Reads chip: show this sample's alignments right under its coverage */}
          {!track.gtex && !track.group && (() => {
            const active = readsSampleIds.includes(track.sampleId);
            return (
              <g data-export="skip" transform={`translate(${labelW + 4}, 0)`} style={{ cursor: 'pointer' }}
                onClick={e => { e.stopPropagation(); if (active && !readsAll) setShowReads(false); else { setReadsAll(false); setReadsSampleId(track.sampleId); setShowReads(true); } }}
                onMouseDown={e => e.stopPropagation()}>
                <title>{active ? (readsAll ? `Show only ${track.sampleName} reads` : 'Hide the reads track') : `Show ${track.sampleName} reads under this track`}</title>
                <rect x={0} y={0} width={44} height={14} rx={7} fill={active ? color : INK.bg} stroke={active ? color : INK.faint} strokeWidth={0.8} />
                <text x={22} y={10} textAnchor="middle" fill={active ? '#fff' : INK.muted} fontSize={8.5} fontWeight={600}>{active ? 'reads ✓' : 'reads'}</text>
              </g>
            );
          })()}
        </g>
        {status && !gtexWarn && (
          <text x={plotRight - 6} y={yOff + 14} textAnchor="end" fill={status.color} fontSize={9.5}>{status.text}</text>
        )}
        {gtexWarn && (
          <g transform={`translate(${PLOT_LEFT + 8}, ${baseline - COVERAGE_H + 30})`}>
            <title>{gtexWarn}</title>
            <rect x={-4} y={-11} width={Math.min(plotWidth - 12, gtexWarn.length * 5.3 + 24)} height={15} rx={3} fill={INK.bg} opacity={0.92} />
            <text x={0} y={0} fill={UNIQUE_COLOR} fontSize={9}>⚠ {gtexWarn.length > Math.floor((plotWidth - 40) / 5.3) ? gtexWarn.slice(0, Math.floor((plotWidth - 40) / 5.3) - 1) + '…' : gtexWarn}</text>
          </g>
        )}

        {/* Remove (group tracks are managed in the groups dialog) */}
        {!track.group && (
          <g data-export="skip" style={{ cursor: 'pointer' }} onClick={() => removeTrack(track.sampleId)}>
            <title>Remove {track.sampleName}</title>
            <circle cx={plotRight + 16} cy={yOff + 12} r={8} fill={INK.bg} stroke={INK.grid} />
            <text x={plotRight + 16} y={yOff + 15.5} textAnchor="middle" fill={INK.muted} fontSize={12}>×</text>
          </g>
        )}
      </g>
    );
  };

  const renderTranscript = (yOff: number) => {
    const midY = yOff + TRANSCRIPT_H / 2 + 8;
    const exonH = 18, utrH = 9;
    if (!tx) {
      return (
        <text x={PLOT_LEFT + 4} y={yOff + 20} fill={INK.muted} fontSize={11} fontFamily={FONT}>
          {transcriptMissing ? `No transcript model for ${currentGeneName}: ${transcriptMissing}` : 'Loading transcript…'}
        </text>
      );
    }
    const gA = scale.x(tx.start), gB = scale.x(tx.end);
    const gLeft = Math.max(PLOT_LEFT, Math.min(gA, gB)), gRight = Math.min(plotRight, Math.max(gA, gB));

    // Direction chevrons every 56 px along the intron line (5' → 3' after the axis flip)
    const arrows: JSX.Element[] = [];
    for (let ax = gLeft + 28; ax < gRight - 8; ax += 56) {
      arrows.push(<path key={`arr${ax}`} d={`M${ax - 2.5},${midY - 3.5} L${ax + 2.5},${midY} L${ax - 2.5},${midY + 3.5}`} fill="none" stroke={INK.intron} strokeWidth={1.2} />);
    }

    // Exon boxes: CDS at full height, UTR portions at half height
    const boxes: JSX.Element[] = [];
    const openExon = (ex: { start: number; end: number; rank: number }, e: React.MouseEvent) => {
      e.stopPropagation();
      const p = svgPoint(e);
      setPopover(prev => (prev?.kind === 'exon' && prev.exon.rank === ex.rank ? null : { kind: 'exon', exon: ex, x: p.x, y: p.y }));
    };
    const box = (ex: { start: number; end: number; rank: number }, s: number, e: number, h: number, fill: string, key: string) => {
      const a = scale.x(s), b = scale.x(e);
      const left = Math.min(a, b), w = Math.max(1, Math.abs(b - a));
      if (left > plotRight || left + w < PLOT_LEFT) return;
      boxes.push(<rect key={key} x={left} y={midY - h / 2} width={w} height={h} fill={fill} rx={1.5} style={{ cursor: 'pointer' }}
        onMouseDown={ev => ev.stopPropagation()} onClick={ev => openExon(ex, ev)} />);
    };
    for (const ex of tx.exons) {
      if (tx.cdsStart == null || tx.cdsEnd == null) { box(ex, ex.start, ex.end, exonH, INK.exon, `x${ex.rank}`); continue; }
      const cs = Math.max(ex.start, tx.cdsStart), ce = Math.min(ex.end, tx.cdsEnd);
      if (ce > cs) box(ex, cs, ce, exonH, INK.exon, `c${ex.rank}`);
      if (ex.start < Math.min(ex.end, tx.cdsStart)) box(ex, ex.start, Math.min(ex.end, tx.cdsStart), utrH, INK.utr, `u5${ex.rank}`);
      if (Math.max(ex.start, tx.cdsEnd) < ex.end) box(ex, Math.max(ex.start, tx.cdsEnd), ex.end, utrH, INK.utr, `u3${ex.rank}`);
    }

    const introns = intronsOf(tx);
    return (
      <g fontFamily={FONT}>
        <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={transcriptPanelH} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />
        <text x={PLOT_LEFT + 8} y={yOff + 14} fontSize={10}>
          <tspan fill={INK.text} fontWeight={700}>{tx.geneName}</tspan>
          <tspan fill={INK.muted}>{'  '}{tx.transcriptId} · {modelKindLabel(tx)} · {tx.strand > 0 ? '+' : '−'} strand · {tx.exons.length} exons{tx.cdsStart == null ? ' · non-coding' : ''}</tspan>
        </text>
        {tx.strand < 0 && (
          <text x={plotRight - 6} y={yOff + 14} textAnchor="end" fill={UNIQUE_COLOR} fontSize={9.5} fontWeight={700}>
            <title>{`${tx.geneName} is transcribed from the minus strand. The axis is reversed so that the transcript reads 5′→3′ from left to right: genomic positions decrease towards the right, which is the opposite of IGV and of the UCSC browser. Reference bases in the reads track are shown on both strands.`}</title>
            ⚠ antisense gene (− strand): axis reversed, 5′→3′ left to right, genomic positions decrease to the right
          </text>
        )}
        {neighbourError && <text x={plotRight - 6} y={yOff + (tx.strand < 0 ? 26 : 14)} textAnchor="end" fill={UNIQUE_COLOR} fontSize={9}>neighbouring genes unavailable: {neighbourError}</text>}
        <defs><clipPath id="sashimi-clip-tx"><rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={transcriptPanelH} /></clipPath></defs>
        <g clipPath="url(#sashimi-clip-tx)">
        {gRight > gLeft && <line x1={gLeft} y1={midY} x2={gRight} y2={midY} stroke={INK.intron} strokeWidth={1.5} />}
        {arrows}
        {/* Intron lengths (equal-intron mode only, where the axis no longer tells) */}
        {axis.kind === 'equal-intron' && introns.map(i => {
          const a = scale.x(i.start), b = scale.x(i.end);
          const w = Math.abs(b - a), mx = (a + b) / 2;
          if (w < 30 || mx < PLOT_LEFT || mx > plotRight) return null;
          return <text key={`il${i.start}`} x={mx} y={midY - 12} textAnchor="middle" fill={INK.faint} fontSize={8.5}>{formatBp(i.end - i.start)}</text>;
        })}
        {boxes}
        {tx.exons.map(ex => {
          const inCds = tx.cdsStart != null && tx.cdsEnd != null && ex.end > tx.cdsStart && ex.start < tx.cdsEnd;
          // Centre the number on the coding part, so it never straddles the CDS/UTR step
          const s0 = inCds ? Math.max(ex.start, tx.cdsStart!) : ex.start;
          const e0 = inCds ? Math.min(ex.end, tx.cdsEnd!) : ex.end;
          const a = scale.x(s0), b = scale.x(e0);
          const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
          const w = right - left;
          const cx = left + w / 2;
          if (w < 14) return null;
          return (
            <g key={`n${ex.rank}`} style={{ cursor: 'pointer' }} onMouseDown={ev => ev.stopPropagation()} onClick={ev => openExon(ex, ev)}>
              <title>{`Exon ${ex.rank} · ${currentChrom}:${(ex.start + 1).toLocaleString()}-${ex.end.toLocaleString()} · ${ex.end - ex.start} bp · click for ψ`}</title>
              <text x={cx} y={inCds ? midY + 3.5 : midY + 20} textAnchor="middle" fill={inCds ? '#ffffff' : INK.muted} fontSize={8.5} fontWeight={600}>{ex.rank}</text>
            </g>
          );
        })}
        {neighbourRows.map((row, ri) => renderNeighbourRow(row, yOff + TRANSCRIPT_H + ri * NEIGHBOUR_ROW_H))}
        </g>
      </g>
    );
  };

  /** One row of neighbouring genes: exons (CDS tall, UTR short), a line with chevrons in the gene's own direction, and a name pill. */
  const renderNeighbourRow = (row: NeighbourModel[], top: number) => {
    const mid = top + NEIGHBOUR_ROW_H / 2 + 1;
    return row.map(m => {
      const same = tx ? m.strand === tx.strand : true;
      const color = same ? SAME_SENSE_COLOR : ANTISENSE_COLOR;
      const a = scale.x(m.start), b = scale.x(m.end);
      const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
      const parts: JSX.Element[] = [];
      if (right > left) parts.push(<line key="l" x1={left} y1={mid} x2={right} y2={mid} stroke={color} strokeWidth={1.2} opacity={0.8} />);
      // Chevrons point where the gene's 3′ end lies on screen: its own strand, after the axis flip of the queried gene
      const toRight = (m.strand > 0) !== reverse;
      for (let ax = left + 22; ax < right - 6; ax += 48) {
        const d = toRight ? `M${ax - 2.5},${mid - 3} L${ax + 2.5},${mid} L${ax - 2.5},${mid + 3}` : `M${ax + 2.5},${mid - 3} L${ax - 2.5},${mid} L${ax + 2.5},${mid + 3}`;
        parts.push(<path key={`a${ax}`} d={d} fill="none" stroke={color} strokeWidth={1.2} />);
      }
      const box = (s0: number, e0: number, hh: number, key: string) => {
        const xa = scale.x(s0), xb = scale.x(e0);
        const l = Math.min(xa, xb), w = Math.max(1, Math.abs(xb - xa));
        if (l > plotRight || l + w < PLOT_LEFT) return;
        parts.push(<rect key={key} x={l} y={mid - hh / 2} width={w} height={hh} fill={color} rx={1} />);
      };
      m.exons.forEach((ex, k) => {
        if (m.cdsStart == null || m.cdsEnd == null) { box(ex.start, ex.end, 12, `x${k}`); return; }
        const cs = Math.max(ex.start, m.cdsStart), ce = Math.min(ex.end, m.cdsEnd);
        if (ce > cs) box(cs, ce, 12, `c${k}`);
        if (ex.start < Math.min(ex.end, m.cdsStart)) box(ex.start, Math.min(ex.end, m.cdsStart), 6, `u5${k}`);
        if (Math.max(ex.start, m.cdsEnd) < ex.end) box(Math.max(ex.start, m.cdsEnd), ex.end, 6, `u3${k}`);
      });
      const label = `${m.geneName} ${toRight ? '→' : '←'}`;
      const lw = label.length * 5.6 + 10;
      // Name pill: beside a short gene (so its exons stay visible), inside the extent of a long one
      const inside = right - left > lw + 60;
      const lx = inside ? Math.max(PLOT_LEFT + 3, Math.min(left, plotRight - lw - 3))
        : right + 4 + lw <= plotRight ? right + 4 : Math.max(PLOT_LEFT + 3, left - lw - 4);
      const relation = tx ? (same ? `same strand as ${tx.geneName}` : `antisense to ${tx.geneName}`) : '';
      const title = `${m.geneName} · ${m.transcriptId}${m.isCanonical ? (isEnsemblId(m.transcriptId) ? ' (Ensembl canonical)' : ' (RefSeq Select)') : ''} · ${m.biotype}\n` +
        `${currentChrom}:${(m.start + 1).toLocaleString()}-${m.end.toLocaleString()} · ${m.strand > 0 ? '+' : '−'} strand · ${relation}\n` +
        `${m.exons.length} exon${m.exons.length > 1 ? 's' : ''}${m.cdsStart == null ? ' · non-coding' : ''}`;
      return (
        <g key={m.geneId}>
          <title>{title}</title>
          {parts}
          <rect x={lx} y={mid - 7} width={lw} height={14} rx={3} fill={INK.bg} opacity={0.92} stroke={color} strokeWidth={0.6} />
          <text x={lx + 5} y={mid + 3.5} fill={color} fontSize={9} fontWeight={700}>{label}</text>
        </g>
      );
    });
  };

  /** Common-SNP track: one lollipop per variant, stem height ∝ log allele frequency, SNVs blue, indels amber. */
  const renderSnps = (yOff: number) => {
    const base = yOff + SNP_PANEL_H - 6;
    const list = visibleSnps;
    const dense = list.length > SNP_MAX_MARKS;
    const stem = (af: number) => 6 + 18 * Math.min(1, Math.log10(Math.max(af, 0.001) / 0.001) / Math.log10(0.5 / 0.001));
    const status = snpStatus.error ? snpStatus.error : snpStatus.loading ? 'loading…' : !snps ? '' : `${list.length.toLocaleString()} variant${list.length === 1 ? '' : 's'} with AF ≥ ${(snpMinAf * 100).toFixed(snpMinAf < 0.01 ? 1 : 0)}%${dense ? ' (ticks only, zoom in)' : ''}`;
    const marks: JSX.Element[] = [];
    for (const v of list) {
      const xa = scale.x(v.start), xb = scale.x(v.end);
      const left = Math.min(xa, xb), w = Math.abs(xb - xa);
      if (left + w < PLOT_LEFT || left > plotRight) continue;
      const cx = left + w / 2;
      const isSnv = v.cls === 'snv';
      const color = isSnv ? SNP_SNV_COLOR : SNP_INDEL_COLOR;
      if (dense) { marks.push(<line key={v.id + v.start} x1={cx} y1={base - 10} x2={cx} y2={base} stroke={color} strokeWidth={1} opacity={0.6} />); continue; }
      const h = stem(v.maxAf);
      marks.push(
        <g key={v.id + v.start}>
          <title>{snpText(v)}</title>
          <line x1={cx} y1={base} x2={cx} y2={base - h} stroke={color} strokeWidth={1} opacity={0.7} />
          {isSnv
            ? <circle cx={cx} cy={base - h} r={2.6} fill={color} stroke={INK.bg} strokeWidth={0.6} />
            : <rect x={cx - Math.max(2, w / 2)} y={base - h - 2.5} width={Math.max(4, w)} height={5} rx={1} fill={color} stroke={INK.bg} strokeWidth={0.6} />}
        </g>,
      );
    }
    return (
      <g fontFamily={FONT}>
        <defs><clipPath id="sashimi-clip-snp"><rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={SNP_PANEL_H} /></clipPath></defs>
        <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={SNP_PANEL_H} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />
        <text x={PLOT_LEFT + 8} y={yOff + 13} fontSize={10}>
          <tspan fill={INK.text} fontWeight={700}>Common SNPs</tspan>
          <tspan fill={snpStatus.error ? UNIQUE_COLOR : INK.muted}>{'  '}{snpSourceLabel()} · {status}</tspan>
        </text>
        <line x1={PLOT_LEFT} y1={base} x2={plotRight} y2={base} stroke={INK.grid} strokeWidth={1} />
        <g clipPath="url(#sashimi-clip-snp)">{marks}</g>
      </g>
    );
  };

  /** A known variant as a mark: a diamond for point variants, a band for CNVs / SVs; clamped to the plot edges. */
  const knownMark = (v: KnownVariant, mid: number, small: boolean, withLabel: boolean): JSX.Element | null => {
    const color = KNOWN_VARIANT_COLORS[v.kind];
    const xa = scale.x(v.start), xb = scale.x(v.end);
    const left = Math.min(xa, xb), right = Math.max(xa, xb);
    if (right < PLOT_LEFT - 1 && left < PLOT_LEFT - 1 && !withLabel) return null;
    if (left > plotRight + 1 && right > plotRight + 1 && !withLabel) return null;
    const point = isPointVariant(v);
    const cx = (left + right) / 2;
    const off = cx < PLOT_LEFT ? 'left' : cx > plotRight ? 'right' : null;
    const r = small ? 3.5 : 5;
    const stop = (e: React.MouseEvent) => e.stopPropagation();
    const onClick = (e: React.MouseEvent) => { e.stopPropagation(); jumpToVariant(v); };
    if (off) {
      if (!withLabel) return null;
      const x = off === 'left' ? PLOT_LEFT + 6 : plotRight - 6;
      const dir = off === 'left' ? -1 : 1;
      // distance from the window edge, in genomic bp (reverse-strand genes flip which end is which)
      const edgePos = scale.invert(off === 'left' ? PLOT_LEFT : plotRight);
      const dist = off === 'left' ? Math.abs(Math.round(edgePos - (reverse ? v.start : v.end))) : Math.abs(Math.round((reverse ? v.end : v.start) - edgePos));
      const lab = `${v.label} · ${formatBp(dist)} ${(off === 'left') !== reverse ? 'upstream' : 'downstream'}`;
      return (
        <g key={v.id} style={{ cursor: 'pointer' }} onMouseDown={stop} onClick={onClick}>
          <path d={`M${x},${mid - 5} L${x + dir * 7},${mid} L${x},${mid + 5} Z`} fill={color} />
          <text x={x + (off === 'left' ? 11 : -11)} y={mid + 3.5} textAnchor={off === 'left' ? 'start' : 'end'} fill={color} fontSize={9} fontWeight={600}>{lab}</text>
        </g>
      );
    }
    const l = Math.max(PLOT_LEFT, left), rr = Math.min(plotRight, right);
    const label = withLabel ? (
      <text x={point ? cx + r + 3 : Math.min(l + 4, plotRight - v.label.length * 5.4 - 4)} y={mid + 3.5} fill={color} fontSize={9} fontWeight={600}>{v.label}</text>
    ) : null;
    return (
      <g key={v.id} style={{ cursor: 'pointer' }} onMouseDown={stop} onClick={onClick}>
        {point
          ? <path d={`M${cx},${mid - r} L${cx + r},${mid} L${cx},${mid + r} L${cx - r},${mid} Z`} fill={color} stroke={INK.bg} strokeWidth={0.8} />
          : <rect x={l} y={mid - (small ? 3 : 5)} width={Math.max(2, rr - l)} height={small ? 6 : 10} rx={1.5} fill={color} opacity={0.55} stroke={color} strokeWidth={0.8} />}
        {label}
      </g>
    );
  };

  /** Known-variant panel: the primary sample's identified variants stacked in rows, with the ones off-screen as edge arrows. */
  const renderKnown = (yOff: number) => {
    return (
      <g fontFamily={FONT}>
        <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={knownPanelH} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />
        <text x={PLOT_LEFT + 8} y={yOff + knownRowMid(0) + 3.5} fontSize={10}>
          <tspan fill={INK.text} fontWeight={700}>Known variants</tspan>
          <tspan fill={INK.muted}>{'  '}{knownStatus}</tspan>
        </text>
        {knownRows.map((row, i) => row.map(v => knownMark(v, yOff + knownRowMid(i), false, true)))}
      </g>
    );
  };

  /** Guide lines / bands of the primary sample's known variants through every track, and small marks on the other samples' tracks. */
  const renderKnownOverlay = () => {
    if (!showKnown) return null;
    const top = knownY + knownPanelH, bottom = tracksBottom - TRACK_GAP;
    const guides = primaryKnownHere.map(v => {
      const xa = scale.x(v.start), xb = scale.x(v.end);
      const left = Math.max(PLOT_LEFT, Math.min(xa, xb)), right = Math.min(plotRight, Math.max(xa, xb));
      if (right < PLOT_LEFT || left > plotRight) return null;
      const color = KNOWN_VARIANT_COLORS[v.kind];
      if (isPointVariant(v)) {
        const cx = (Math.min(xa, xb) + Math.max(xa, xb)) / 2;
        if (cx < PLOT_LEFT || cx > plotRight) return null;
        return <line key={v.id} x1={cx} y1={top} x2={cx} y2={bottom} stroke={color} strokeWidth={1} strokeDasharray="3 3" opacity={0.7} />;
      }
      return (
        <g key={v.id}>
          <rect x={left} y={top} width={Math.max(1, right - left)} height={Math.max(0, bottom - top)} fill={color} opacity={0.07} />
          {Math.min(xa, xb) >= PLOT_LEFT && <line x1={left} y1={top} x2={left} y2={bottom} stroke={color} strokeWidth={1} strokeDasharray="4 3" opacity={0.6} />}
          {Math.max(xa, xb) <= plotRight && <line x1={right} y1={top} x2={right} y2={bottom} stroke={color} strokeWidth={1} strokeDasharray="4 3" opacity={0.6} />}
        </g>
      );
    });
    const marks = layouts.filter(L => L.track.sampleId >= 0 && !L.track.gtex).flatMap(L => knownOnChrom(L.track.sampleId).map(v => knownMark(v, L.yOff + L.juncH - 8, true, false)));
    return <g fontFamily={FONT}><g pointerEvents="none">{guides}</g>{marks}</g>;
  };

  const renderAltTranscripts = (yOff: number) => {
    const h = altPanelH;
    const src = altTx?.data.source === 'refseq' ? 'RefSeq (UCSC)' : 'Ensembl';
    const models = altModels || [];
    const manes = new Set<string>();
    for (const ex of tx?.exons || []) { manes.add(`s${ex.start}`); manes.add(`e${ex.end}`); }
    const rows = models.map((m, i) => {
      const top = yOff + ALT_TX_HEADER_H + i * ALT_TX_ROW_H, mid = top + ALT_TX_ROW_H / 2;
      const a = scale.x(m.start), b = scale.x(m.end);
      const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
      const parts: JSX.Element[] = [];
      if (right > left) parts.push(<line key="l" x1={left} y1={mid} x2={right} y2={mid} stroke={INK.intron} strokeWidth={1} />);
      m.exons.forEach((ex, k) => {
        const novel = !manes.has(`s${ex.start}`) && !manes.has(`e${ex.end}`);
        const fill = novel ? '#b45309' : '#64748b';
        const box = (s0: number, e0: number, hh: number, key: string) => {
          const xa = scale.x(s0), xb = scale.x(e0);
          const l = Math.min(xa, xb), w = Math.max(1, Math.abs(xb - xa));
          if (l > plotRight || l + w < PLOT_LEFT) return;
          parts.push(<rect key={key} x={l} y={mid - hh / 2} width={w} height={hh} fill={fill} rx={1} />);
        };
        if (m.cdsStart == null || m.cdsEnd == null) { box(ex.start, ex.end, 12, `x${k}`); return; }
        const cs = Math.max(ex.start, m.cdsStart), ce = Math.min(ex.end, m.cdsEnd);
        if (ce > cs) box(cs, ce, 12, `c${k}`);
        if (ex.start < Math.min(ex.end, m.cdsStart)) box(ex.start, Math.min(ex.end, m.cdsStart), 6, `u5${k}`);
        if (Math.max(ex.start, m.cdsEnd) < ex.end) box(Math.max(ex.start, m.cdsEnd), ex.end, 6, `u3${k}`);
      });
      const shown = !!tx && m.id === tx.transcriptId;
      const label = `${m.id}${m.is_mane ? ' · MANE' : ''}${shown ? ' · shown' : ''}`;
      const lw = label.length * 5.4 + 8;
      const nNovel = m.exons.filter(ex => !manes.has(`s${ex.start}`) && !manes.has(`e${ex.end}`)).length;
      const title = `${m.id}${m.name && m.name !== m.id ? ` · ${m.name}` : ''} · ${m.biotype}\n${currentChrom}:${(m.start + 1).toLocaleString()}-${m.end.toLocaleString()} · ${m.exons.length} exons` +
        (m.cdsStart == null ? ' · non-coding' : '') + (m.is_mane ? '\nsame exon structure as the MANE Select transcript' : shown ? '\nthe model displayed on the top track' : nNovel ? `\n${nNovel} exon${nNovel > 1 ? 's' : ''} absent from the displayed model (amber)` : '') +
        (shown ? '' : '\nclick to display this model as the reference (exon numbering, junction classes, HGVS, usage)');
      return (
        <g key={m.id} style={{ cursor: shown ? 'default' : 'pointer' }} onMouseDown={e => e.stopPropagation()}
          onClick={e => { e.stopPropagation(); if (!shown) chooseModel(m.id); }}>
          <title>{title}</title>
          <rect x={PLOT_LEFT} y={top} width={plotWidth} height={ALT_TX_ROW_H} fill={shown ? '#4f46e5' : 'transparent'} opacity={shown ? 0.08 : 0} />
          {parts}
          <rect x={PLOT_LEFT + 3} y={mid - 7} width={lw} height={14} rx={3} fill={INK.bg} opacity={0.9} stroke={shown ? '#4f46e5' : 'none'} strokeWidth={0.8} />
          <text x={PLOT_LEFT + 7} y={mid + 3.5} fill={shown ? '#4338ca' : m.is_mane ? INK.text : INK.muted} fontSize={9} fontWeight={m.is_mane || shown ? 700 : 500}>{label}</text>
        </g>
      );
    });
    const status = altTxError ? altTxError : !altModels ? 'loading…' : models.length === 0 ? 'no other transcript models' : '';
    return (
      <g fontFamily={FONT}>
        <defs><clipPath id="sashimi-clip-alt"><rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={h} /></clipPath></defs>
        <rect x={PLOT_LEFT} y={yOff} width={plotWidth} height={h} fill="none" stroke={INK.grid} strokeWidth={1} rx={4} />
        <text x={PLOT_LEFT + 8} y={yOff + 14} fontSize={10}>
          <tspan fill={INK.text} fontWeight={700}>All transcripts</tspan>
          <tspan fill={altTxError ? UNIQUE_COLOR : INK.muted}>{'  '}{altModels ? `${src} · ${altTx?.data.transcripts.length ?? 0} model${(altTx?.data.transcripts.length ?? 0) === 1 ? '' : 's'}${(altTx?.data.transcripts.length ?? 0) > ALT_TX_MAX_ROWS ? ` (first ${ALT_TX_MAX_ROWS})` : ''} · amber exon = absent from ${tx?.modelKind === 'mane' ? 'MANE Select' : 'the displayed model'} · click a model to make it the reference` : status}</tspan>
        </text>
        <g clipPath="url(#sashimi-clip-alt)">{rows}</g>
      </g>
    );
  };

  // ---- Hover crosshair + tooltip ----
  const hoverInfo = useMemo(() => {
    if (!hover || dragging || regionSelect) return null;
    const pos = Math.floor(scale.invert(hover.px));
    // Over a row of the "All transcripts" panel: c./n. position on that transcript, not on the MANE model
    let alt: { id: string; label: string; kind: string; biotype: string } | null = null;
    if (showAllTx && altModels && hover.py >= altY + ALT_TX_HEADER_H) {
      const i = Math.floor((hover.py - altY - ALT_TX_HEADER_H) / ALT_TX_ROW_H);
      const m = altModels[i];
      if (m && pos >= m.start - 1 && pos < m.end) {
        const n = m.exons.length;
        const model: TxModel = {
          geneName: currentGeneName, transcriptId: m.id, isMane: m.is_mane, modelKind: m.is_mane ? 'mane' : m.is_canonical ? 'canonical' : 'longest', biotype: m.biotype,
          chrom: currentChrom, strand: m.strand < 0 ? -1 : 1, start: m.start, end: m.end,
          exons: m.exons.map((e, k) => ({ ...e, rank: m.strand < 0 ? n - k : k + 1 })), cdsStart: m.cdsStart, cdsEnd: m.cdsEnd,
        };
        const c = cdnaPosition(pos, model);
        alt = { id: m.id, label: c.label, kind: c.kind, biotype: m.biotype };
      }
    }
    const here = showSnps ? visibleSnps.filter(v => pos >= v.start && pos < v.end) : [];
    const known: { v: KnownVariant; sample: string | null; off: boolean }[] = [];
    if (showKnown) {
      const px = hover.px;
      const extent = (v: KnownVariant) => { const xa = scale.x(v.start), xb = scale.x(v.end); return { l: Math.min(xa, xb), r: Math.max(xa, xb) }; };
      const hit = (v: KnownVariant) => { const { l, r } = extent(v); return isPointVariant(v) ? Math.abs((l + r) / 2 - px) <= 6 : px >= l && px <= r; };
      for (const v of primaryKnownHere) if (hit(v)) known.push({ v, sample: null, off: false });
      // edge arrows of the off-screen variants, in their panel row
      knownRows.forEach((row, i) => {
        const y0 = knownY + knownRowMid(i) - KNOWN_ROW_H / 2;
        if (hover.py < y0 || hover.py > y0 + KNOWN_ROW_H) return;
        for (const v of row) {
          const { l, r } = extent(v); const cx = (l + r) / 2;
          if (cx >= PLOT_LEFT && cx <= plotRight) continue;
          const left = cx < PLOT_LEFT, x = left ? PLOT_LEFT + 6 : plotRight - 6, w = (v.label.length + 24) * 5.4;
          if ((left ? px >= x && px <= x + w : px <= x && px >= x - w) && !known.some(k => k.v === v)) known.push({ v, sample: null, off: true });
        }
      });
      for (const L of layouts) {
        if (L.idx === 0 || L.track.sampleId < 0 || L.track.gtex || hover.py < L.yOff || hover.py > L.yOff + L.height) continue;
        for (const v of knownOnChrom(L.track.sampleId)) if (hit(v)) known.push({ v, sample: L.track.sampleName, off: false });
      }
    }
    return {
      x: hover.px, pos, alt, snps: here, known,
      cdna: tx ? cdnaPosition(pos, tx) : null,
      rows: layouts.map(L => ({ name: L.track.sampleName, color: L.color, depth: depthAt(L.track.coverage, pos) })),
    };
  }, [hover, dragging, regionSelect, scale, layouts, tx, showAllTx, altModels, altY, currentGeneName, currentChrom, showSnps, visibleSnps, showKnown, primaryKnownHere, knownRows, knownY, knownOnChrom, plotRight]);

  // ---- Popover content: junction (HGVS, frame, share vs canonical, usage of skipped exons) or exon (depth usage + junction ψ) ----
  interface PopTable { caption?: string; head: string[]; rows: string[][] }
  interface PopStrip { values: { value: number; name: string; color: string | null }[]; median: number | null }
  const popoverContent = useMemo((): { title: string; subtitle: string; hgvs: string[]; tables: PopTable[]; strip: PopStrip | null; note: string; status?: string; cartoon?: { j: JunctionArc; model: TxModel; label: string } | null } | null => {
    if (!popover) return null;
    const pct = (v: number | null | undefined, d = 0) => (v == null || !Number.isFinite(v) ? '—' : `${(v * 100).toFixed(d)}%`);
    const num = (v: number | null | undefined, d = 1) => (v == null || !Number.isFinite(v) ? '—' : v.toFixed(d));
    const status = exonDepths?.key === usageKey ? (exonDepths?.error ? `exon depth unavailable: ${exonDepths.error}` : exonDepths?.data ? undefined : 'computing exon depth over the run…') : 'computing exon depth over the run…';
    const strandNote = (() => {
      const d = exonDepths?.data; if (!d || exonDepths?.key !== usageKey) return '';
      const calls = new Map<string, number>();
      for (const smp of d.samples) if (!smp.error) calls.set(smp.strandness, (calls.get(smp.strandness) || 0) + 1);
      const parts = [...calls.entries()].map(([k, n]) => `${n} ${k === 'firststrand' ? 'stranded (dUTP), sense reads only' : k === 'secondstrand' ? 'stranded (forward), sense reads only' : k === 'unstranded' ? 'unstranded, all reads' : 'strand unknown, all reads'}`);
      return parts.length ? ` Library: ${parts.join('; ')}.` : '';
    })();
    /** usage table rows for the exons `idxs`, loaded tracks in detail + cohort summary */
    const usageTable = (idx: number, caption?: string): PopTable => {
      const u = usageOf(idx);
      const rows: string[][] = [];
      if (!u) return { caption, head: ['sample', 'exon depth', 'gene depth', 'usage', 'vs controls', 'z'], rows: [[status || '—', '', '', '', '', '']] };
      for (const t of tracks) {
        const r = u.forSample(t.sampleId);
        if (!r.me) { rows.push([t.sampleName, 'no data', '', '', '', '']); continue; }
        const sd = r.me.u.sd != null ? ` ± ${pct(r.me.u.sd)}` : '';
        rows.push([t.sampleName, num(r.me.u.depth, 0), num(r.me.u.refDepth, 0), `${pct(r.me.u.usage)}${sd}`, pct(r.rel), r.z == null ? '—' : r.z.toFixed(1)]);
      }
      const c = usageCohort(u.values);
      rows.push([`run (n=${c.n}${u.failed ? `, ${u.failed} without BAM` : ''})`, '', '', `median ${pct(c.median)} · MAD ${pct(c.mad)}`, '', '']);
      return { caption, head: ['sample', 'exon depth', 'gene depth', 'usage', 'vs controls', 'z'], rows };
    };
    const stripFor = (idx: number): PopStrip | null => {
      const u = usageOf(idx); if (!u) return null;
      const colorOf = new Map(tracks.map((t, i) => [t.sampleId, TRACK_COLORS[i % TRACK_COLORS.length]]));
      return { values: u.per.filter(p => p.u.usage != null).map(p => ({ value: p.u.usage!, name: p.sample.sample_name, color: colorOf.get(p.sample.sample_id) ?? null })), median: usageCohort(u.values).median };
    };
    const exonIdx = (rank: number) => tx ? tx.exons.findIndex(e => e.rank === rank) : -1;

    if (popover.kind === 'junction') {
      const j = popover.j;
      const { model, info, foreign } = junctionContext(j);
      const allJunctions = tracks.flatMap(t => t.junctions);
      const hgvs = model ? junctionHgvs(j, model, allJunctions) : null;
      const frame = model && info.cls !== 'canonical' ? junctionFrame(j, model, allJunctions) : null;
      const key = junctionKey(j);
      const inAlt = altJunctionIndex.get(key);
      const hgvsLines = hgvs ? [`${model!.transcriptId}: donor ${hgvs.donor} → acceptor ${hgvs.acceptor}`, `${hgvs.effect} · ${hgvs.summary}`] : ['no MANE Select model loaded'];
      if (foreign && tx) hgvsLines.unshift(`${foreign.geneName} (${foreign.strand === tx.strand ? 'same strand as' : 'antisense to'} ${tx.geneName}) · ${foreign.transcriptId}`);
      if (frame) hgvsLines.push(`reading frame: ${frameLabel(frame)} · ${frame.text}`);
      const alts = displayTracks.map(t => {
        const mine = t.junctions.find(k => junctionKey(k) === key);
        const alt = model ? junctionAlternative({ ...j, count: mine?.count ?? 0 }, t.junctions, model) : null;
        return { name: t.sampleName, count: mine?.count ?? 0, alt };
      });
      const altLabel = alts.find(a => a.alt)?.alt?.label ?? '';
      const tables: PopTable[] = [{
        head: ['sample', 'reads', 'canonical alternative', 'share'],
        rows: alts.map(a => [a.name, a.count.toLocaleString(), a.alt?.canonical == null ? '—' : num(a.alt.canonical, a.alt.canonical % 1 ? 1 : 0), pct(a.alt?.share)]),
      }];
      // skipped exons of the queried gene: their depth usage
      const strips: PopStrip[] = [];
      if (model && model === tx && info.leftExon != null && info.rightExon != null) {
        const li = exonIdx(info.leftExon), ri = exonIdx(info.rightExon);
        for (let i = Math.min(li, ri) + 1; i < Math.max(li, ri); i++) {
          tables.push(usageTable(i, `Exon ${tx!.exons[i].rank} · depth-based usage`));
          const st = stripFor(i); if (st) strips.push(st);
        }
      }
      const canCartoon = !!model && model.cdsStart != null && ['canonical', 'skip', 'cryptic_exon', 'exonic_site', 'intronic_site'].includes(spliceEvent(j, model, allJunctions).kind);
      return {
        title: `Junction ${currentChrom}:${(j.start + 1).toLocaleString()}-${j.end.toLocaleString()}`,
        subtitle: `intron ${formatBp(j.end - j.start)} · ${info.label}${inAlt ? ` · annotated in ${inAlt.slice(0, 3).join(', ')}${inAlt.length > 3 ? '…' : ''}` : ''}`,
        cartoon: canCartoon ? { j, model: model!, label: info.label } : null,
        hgvs: hgvsLines, tables, strip: strips[0] ?? null,
        note: (altLabel ? `Share = reads of this junction / (reads + ${altLabel}), junction reads only (rMATS-style).` : 'A canonical junction has no single alternative; click an exon for its usage.') +
          (tables.length > 1 ? ' Usage of the skipped exon(s) below is depth-based (see exon panel).' : ''),
        status: tables.length > 1 ? status : undefined,
      };
    }
    const ex = popover.exon;
    const cFirst = tx ? cdnaPosition(tx.strand > 0 ? ex.start : ex.end - 1, tx).label : null;
    const cLast = tx ? cdnaPosition(tx.strand > 0 ? ex.end - 1 : ex.start, tx).label : null;
    const idx = exonIdx(ex.rank);
    const u = idx >= 0 ? usageOf(idx) : null;
    const iv = u?.interval;
    const psiRows = displayTracks.map(t => ({ name: t.sampleName, ...exonPsi(ex, t.junctions, tx?.strand ?? 1) }));
    const tables: PopTable[] = [
      usageTable(idx, `Depth-based usage${iv ? ` · ${iv.coding ? 'coding part ' : ''}${currentChrom}:${(iv.start + 1).toLocaleString()}-${iv.end.toLocaleString()}` : ''}`),
      { caption: 'Junction reads', head: ['sample', 'inclusion 5′', 'inclusion 3′', 'skipping', 'ψ (inclusion)'],
        rows: psiRows.map(r => [r.name, r.inclusionUp.toLocaleString(), r.inclusionDown.toLocaleString(), r.exclusion.toLocaleString(), pct(r.psi)]) },
    ];
    return {
      title: `Exon ${ex.rank} · ${currentChrom}:${(ex.start + 1).toLocaleString()}-${ex.end.toLocaleString()} · ${ex.end - ex.start} bp`,
      subtitle: cFirst && cLast ? `${tx!.transcriptId}: ${cFirst}_${cLast.replace(/^[cn]\./, '')}` : '',
      hgvs: [] as string[], tables, strip: idx >= 0 ? stripFor(idx) : null,
      note: `Usage = median depth of the exon / median depth of the gene's other coding exons (${u ? u.refIdx.length : '…'} exons); ± is a delta-method sd from the read counts. ` +
        'vs controls = usage / median usage of the other samples of the run; z = robust z-score (median/MAD, ≥ 5 controls).' + strandNote +
        ' ψ = mean inclusion junction reads / (mean inclusion + skipping reads).',
      status,
    };
  }, [popover, tx, tracks, displayTracks, altJunctionIndex, currentChrom, junctionContext, exonDepths, usageKey, usageOf]);

  useEffect(() => {
    if (!popover) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') setPopover(null); };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [popover]);

  // ---- Picker samples: every sample loaded in the page, the ones shown as a track highlighted ----
  const filteredSamples = useMemo(() => {
    const q = pickerSearch.toLowerCase();
    return runSamples.filter(s => !q || s.name.toLowerCase().includes(q));
  }, [runSamples, pickerSearch]);
  /** Picker row click: show the sample as a track, or remove its track when it is already shown. */
  const toggleSample = useCallback((s: { id: number; name: string }) => {
    if (tracksRef.current.some(t => t.sampleId === s.id)) removeTrack(s.id); else loadCoverage(s.id, s.name);
  }, [removeTrack, loadCoverage]);

  // ======================== Main render (always light theme for readability) ========================

  // Report every option and the navigation to the host (session files, options kept across remounts)
  const onStateChangeRef = useRef(onStateChange);
  onStateChangeRef.current = onStateChange;
  useEffect(() => {
    onStateChangeRef.current?.({
      equalIntrons, intronWidth, allTranscripts: showAllTx, commonSnps: showSnps, snpMinAf, depthAxis, uniqueOnly,
      reads: showReads, readsAll, readsSample: readsSampleId, collapseReads, minVafPct,
      minJunctionReads: minJunctionCount, minUsagePct, arcLabels: arcLabel, intronRetention: includeRetention,
      viewMode, groups: groups.map(g => ({ name: g.name, sampleIds: [...g.sampleIds] })), knownVariants: showKnown,
      transcriptId: transcript?.model_kind === 'chosen' ? transcript.transcript_id : undefined,
      gene: { name: currentGeneName, id: currentGeneId, chrom: currentChrom, start: currentGeneStart + 1, end: currentGeneEnd },
      view: { chrom: currentChrom, start: viewStart + 1, end: viewEnd },
      mark: locusMark ? { start: locusMark.start + 1, end: locusMark.end } : null,
    });
  }, [equalIntrons, intronWidth, showAllTx, showSnps, snpMinAf, depthAxis, uniqueOnly, showReads, readsAll, readsSampleId, collapseReads, minVafPct, minJunctionCount, minUsagePct, arcLabel, includeRetention, viewMode, groups, showKnown, transcript, currentGeneName, currentGeneId, currentChrom, currentGeneStart, currentGeneEnd, viewStart, viewEnd, locusMark]);

  const t = {
    bg: 'bg-white', text: 'text-gray-900', muted: 'text-gray-500', border: 'border-gray-200',
    inp: 'bg-white text-gray-800 border-gray-300',
    btn: 'px-2 py-0.5 text-xs rounded border border-gray-200 hover:bg-indigo-50 hover:border-indigo-300 transition-colors',
  };
  /** Pill-style segmented switch: the active option is a raised white chip with an indigo label. */
  const Segmented = <T extends string>({ value, onChange, options, disabled, title }: {
    value: T; onChange: (v: T) => void; disabled?: boolean; title: string;
    options: { value: T; label: string; icon: JSX.Element; hint?: string }[];
  }) => (
    <span title={title} className={`inline-flex items-center rounded-full bg-gray-100 border border-gray-200 p-0.5 text-xs select-none ${disabled ? 'opacity-60' : ''}`}>
      {options.map(o => {
        const active = o.value === value;
        return (
          <button key={o.value} type="button" disabled={disabled} onClick={() => onChange(o.value)} title={o.hint}
            className={`flex items-center gap-1.5 px-2.5 py-0.5 rounded-full transition-all ${active ? 'bg-white text-indigo-700 font-semibold shadow-sm ring-1 ring-indigo-200' : 'text-gray-500 hover:text-gray-800'} disabled:cursor-not-allowed`}>
            <span className={active ? 'text-indigo-600' : 'text-gray-400'}>{o.icon}</span>{o.label}
          </button>
        );
      })}
    </span>
  );
  const ICON = {
    reads: <svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round"><path d="M1.5 9.5c0-4 2-7 4.5-7s4.5 3 4.5 7" /><path d="M1 9.5h10" /></svg>,
    usage: <svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"><circle cx="3.2" cy="3.2" r="1.7" /><circle cx="8.8" cy="8.8" r="1.7" /><path d="M10 2 2 10" /></svg>,
    samples: <svg width="12" height="12" viewBox="0 0 12 12" fill="currentColor"><rect x="1" y="1.5" width="10" height="2" rx="1" /><rect x="1" y="5" width="10" height="2" rx="1" /><rect x="1" y="8.5" width="10" height="2" rx="1" /></svg>,
    groups: <svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" strokeWidth="1.5"><rect x="1" y="1" width="4.2" height="4.2" rx="1" /><rect x="6.8" y="1" width="4.2" height="4.2" rx="1" /><rect x="1" y="6.8" width="4.2" height="4.2" rx="1" /><rect x="6.8" y="6.8" width="4.2" height="4.2" rx="1" /></svg>,
  };
  const Toggle = ({ checked, onChange, label, title, disabled }: { checked: boolean; onChange: (v: boolean) => void; label: string; title: string; disabled?: boolean }) => (
    <label className={`flex items-center gap-1 text-xs ${disabled ? 'text-gray-300' : t.muted} select-none`} title={title}>
      <input type="checkbox" checked={checked} disabled={disabled} onChange={e => onChange(e.target.checked)} className="accent-indigo-600" />
      {label}
    </label>
  );

  const content = (
    <div className={embedded ? `${t.bg} ${t.text} w-full` : `${t.bg} ${t.text} rounded-xl shadow-2xl w-full max-w-[95vw] border ${t.border}`}>
      {/* Header */}
      <div className={`flex flex-wrap items-center justify-between gap-x-4 gap-y-2 px-5 py-3 border-b ${t.border}`}>
        <div className="flex flex-wrap items-center gap-x-4 gap-y-2">
          <div>
            <h2 className="text-lg font-bold leading-tight">{currentGeneName} <span className={`text-sm font-normal ${t.muted}`}>Sashimi plot</span></h2>
            <p className={`text-xs ${t.muted} font-mono`}>{regionStr}</p>
          </div>
          <form onSubmit={e => { e.preventDefault(); navigateToGene(); }} className="flex items-center gap-1">
            <input type="text" value={geneSearch} onChange={e => { setGeneSearch(e.target.value); setSearchError(null); }}
              placeholder="Gene, ENSG or chr:pos…" title="A gene symbol, an ENSG id, or genomic coordinates: chr17:43,094,464 (1 kb window) or chr17:43,000,000-43,100,000; on another chromosome the gene at the locus is opened"
              className={`${t.inp} w-40 px-2 py-0.5 text-xs rounded border ${searchError ? 'border-red-400' : ''}`} />
            <button type="submit" disabled={geneSearchLoading} className={t.btn}>{geneSearchLoading ? '…' : 'Go'}</button>
            {searchError && <span className="text-[10px] text-red-600 max-w-[260px] truncate" title={searchError}>{searchError}</span>}
          </form>
          <div className="flex items-center gap-1.5">
            <button onClick={() => zoomBy(1 / 1.4)} className={`${t.btn} font-bold`} title="Zoom in (Ctrl + scroll up)">+</button>
            <button onClick={() => zoomBy(1.4)} className={`${t.btn} font-bold`} title="Zoom out (Ctrl + scroll down)">&minus;</button>
            <button onClick={resetZoom} className={t.btn} title="Reset to the whole gene (or double-click the plot)">Reset</button>
          </div>
          <div className="flex flex-wrap items-center gap-3">
            <Toggle checked={equalIntrons} onChange={toggleEqualIntrons} disabled={!tx || intronsOf(tx).length === 0} label="Equal introns"
              title="Draw every intron at the same width so exons and junctions dominate the plot. Intronic signal (retention, cryptic exons) is compressed; switch off to inspect it." />
            {equalIntrons && tx && intronsOf(tx).length > 0 && (
              <label className={`flex items-center gap-1 text-xs ${t.muted}`}
                title={`Width given to every intron, in bp-equivalents (one exon base = one unit). Default ${defaultIntronV(tx)}: the median exon length of the model, kept between 80 and 300. Clear the box for the default.`}>
                Intron width
                <input type="number" min={10} max={20000} step={10} value={intronWidth ?? ''} placeholder={String(defaultIntronV(tx))}
                  onChange={e => setIntronWidth(e.target.value === '' ? null : Math.min(20000, Math.max(10, parseInt(e.target.value) || 10)))}
                  className={`${t.inp} w-20 px-1.5 py-0.5 text-xs rounded border`} />
                {intronWidth != null && <button onClick={() => setIntronWidth(null)} className="text-gray-400 hover:text-gray-700" title="Back to the default width">×</button>}
              </label>
            )}
            <Toggle checked={showAllTx} onChange={setShowAllTx} label="All transcripts"
              title="Show every transcript model of the gene under the MANE Select track: RefSeq models (NM_/NR_) from the UCSC API, Ensembl transcripts when the UCSC API is unreachable. Exons absent from the displayed model are amber." />
            <Toggle checked={showSnps} onChange={setShowSnps} label="Common SNPs"
              title="Track of common human variants under the transcript (dbSNP 155 common via the UCSC API, Ensembl variation as fallback): lollipop height follows the highest allele frequency across frequency projects, SNVs blue, indels amber. Variant sites called from the reads that match a common SNP get a blue ring." />
            {showSnps && (
              <label className={`flex items-center gap-1 ${t.muted}`} title="Minimum allele frequency (highest over the frequency projects) for a variant to be shown">AF ≥
                <select value={snpMinAf} onChange={e => setSnpMinAf(parseFloat(e.target.value))} className={`${t.inp} px-1 py-0.5 text-xs rounded border`}>
                  {[[0.001, '0.1%'], [0.01, '1%'], [0.05, '5%']].map(([v, l]) => <option key={String(v)} value={v as number}>{l as string}</option>)}
                </select>
              </label>
            )}
            {primaryKnown.length > 0 && (
              <span className="flex items-center gap-1">
                <Toggle checked={showKnown} onChange={setShowKnown} label="Known variants"
                  title="Variants previously identified in the primary sample (clinical indication, diagnostic conclusion, chromosome-map CNVs / SVs): a panel under the transcript, guide lines through every track, and small marks on the other samples' tracks for their own variants." />
                {showKnown && (
                  <select value="" onChange={e => { const v = primaryKnownHere.find(k => k.id === e.target.value); if (v) jumpToVariant(v); }}
                    className={`${t.inp} px-1 py-0.5 text-xs rounded border`} title="Centre the view on one of the sample's known variants">
                    <option value="">go to…</option>
                    {primaryKnownHere.map(v => <option key={v.id} value={v.id}>{v.label} · {KNOWN_VARIANT_KIND_NAMES[v.kind]} · {(v.start + 1).toLocaleString()}</option>)}
                    {primaryKnown.filter(v => !primaryKnownHere.includes(v)).map(v => <option key={v.id} value={v.id} disabled>{v.label} · {v.chrom || '?'} (other chromosome)</option>)}
                  </select>
                )}
              </span>
            )}
            <label className={`flex items-center gap-1 text-xs ${t.muted} select-none`}
              title="Depth axis. Shared: one axis for all samples (heights comparable). Per sample: each sample scales to its own maximum, rounded to a round number. Relative: each sample drawn as a percentage of its own maximum in view, axis 0–100 %, so profiles are comparable whatever their depth.">
              Depth axis
              <select value={depthAxis} onChange={e => setDepthAxis(e.target.value as DepthAxis)} className={`${t.inp} px-1 py-0.5 text-xs rounded border`}>
                <option value="shared">shared</option>
                <option value="own">per sample</option>
                <option value="relative">relative (% of max)</option>
              </select>
            </label>
            <Toggle checked={uniqueOnly} onChange={setUniqueOnly} label="Unique reads"
              title="Count only uniquely mapped reads (NH:1, or MAPQ ≥ 30 when NH is absent) for coverage, junctions and the reads track." />
            <span className="flex items-center gap-1">
              <Toggle checked={showReads} onChange={setShowReads} label="Reads"
                title={`Show the alignments of the primary sample (or of every sample) in a track below its coverage, IGV-style: base mismatches against the reference genome, insertions, deletions and splice gaps. Loads when the window is below ${formatBp(READS_MAX_VIEW_BP)}.`} />
              {showReads && tracks.length > 1 && (
                <select value={readsAll ? 'all' : (effectiveReadsSampleId ?? '')}
                  onChange={e => { if (e.target.value === 'all') setReadsAll(true); else { setReadsAll(false); setReadsSampleId(parseInt(e.target.value)); } }}
                  className={`${t.inp} px-1 py-0.5 text-xs rounded border`} title="Sample shown in the reads track, or all samples (one reads track under each coverage track; each sample is decoded separately, so it takes longer)">
                  {tracks.map(tr => <option key={tr.sampleId} value={tr.sampleId}>{tr.sampleName}</option>)}
                  <option value="all">All samples</option>
                </select>
              )}
              {showReads && (
                <label className={`flex items-center gap-1 text-xs ${t.muted}`} title="Minimum alternate-allele fraction for a variant site to be shown (★, allele bar on the coverage) and used to collapse reads. Sites also need at least 3 alternate reads with base quality ≥ 20.">
                  Min VAF
                  <input type="number" min={1} max={100} value={minVafPct} onChange={e => setMinVafPct(Math.min(100, Math.max(1, parseInt(e.target.value) || 1)))}
                    className={`${t.inp} w-14 px-1.5 py-0.5 text-xs rounded border`} />%
                </label>
              )}
              {showReads && (
                <Toggle checked={collapseReads} onChange={setCollapseReads} label="Collapse"
                  title={`Collapse the reads into consensus groups: one row per local haplotype × splice pattern with its number of supporting reads. Variable sites (★) need at least 3 alternate reads and the Min VAF fraction of the depth; groups below "Min reads" fold into a minor bucket. Sites never co-covered by a read stay in separate groups (no invented phase).`} />
              )}
            </span>
            <Segmented value={showUsage ? 'usage' : 'reads'} onChange={setArcLabel} disabled={viewMode === 'groups'}
              title={viewMode === 'groups' ? 'The Groups view always shows % usage.' : 'What the arc pills show.'}
              options={[
                { value: 'reads', label: 'Reads', icon: ICON.reads, hint: 'Spliced reads of each junction' },
                { value: 'usage', label: 'Usage', icon: ICON.usage, hint: 'Each arc labelled with its share of the reads competing at its intron, so the labels of one intron add up to 100 %: canonical C, alternative site n, pseudo-exon (A + B) / 2 on both arcs, exon skipping S, intron retention (R5 + R3) / 2 shown as IR pills on the baseline. A skipping arc shows 2·S over the totals of the two introns it spans (the rMATS value when nothing else competes). Tooltips also give each event against the canonical junction alone.' },
              ]} />
            {showUsage ? (
              <>
                <label className={`flex items-center gap-1 text-xs ${t.muted}`} title="Hide events whose usage is below this percentage (junctions without a usage value, touching no annotated splice site, follow Min reads instead). Hidden events still count in the denominators.">
                  Min %
                  <input type="number" min={0} max={100} step={0.5} value={minUsagePct} onChange={e => setMinUsagePct(Math.min(100, Math.max(0, parseFloat(e.target.value) || 0)))}
                    className={`${t.inp} w-16 px-1.5 py-0.5 text-xs rounded border`} />
                </label>
                <Toggle checked={includeRetention} onChange={setIncludeRetention} label="Intron retention"
                  title="Count intron retention in the usage percentages: IR pills on the intron baselines, (R5 + R3) / (R5 + R3 + 2·C) from the reads running unspliced through both boundaries, and retention in the canonical arc's denominator. Off: junction-only percentages, no IR pill." />
              </>
            ) : (
              <label className={`flex items-center gap-1 text-xs ${t.muted}`} title="Hide junctions supported by fewer spliced reads">
                Min reads
                <input type="number" min={1} value={minJunctionCount} onChange={e => setMinJunctionCount(Math.max(1, parseInt(e.target.value) || 1))}
                  className={`${t.inp} w-14 px-1.5 py-0.5 text-xs rounded border`} />
              </label>
            )}
          </div>
        </div>
        <div className="flex items-center gap-2">
          {!hideSamplePicker && (
            <Segmented value={viewMode} title="One track per sample, or one pooled track per sample group"
              onChange={v => { if (v === 'groups' && !groups.some(g => g.sampleIds.length)) setShowGroupsDialog(true); else setViewMode(v); }}
              options={[
                { value: 'samples', label: 'Samples', icon: ICON.samples, hint: 'One track per sample' },
                { value: 'groups', label: groups.length ? `Groups · ${groups.length}` : 'Groups', icon: ICON.groups, hint: 'One pooled track per sample group, arcs labelled with % usage' },
              ]} />
          )}
          {!hideSamplePicker && (
            <button onClick={() => setShowGroupsDialog(true)} className={`${t.btn} px-3 py-1 font-medium`} title="Create and edit the sample groups of the aggregate view">Groups…</button>
          )}
          {!hideSamplePicker && <div className="relative">
            <button onClick={e => openDropdown(e, 256, setShowPicker)} className={`${t.btn} px-3 py-1 font-medium ${showPicker ? 'bg-indigo-50 border-indigo-300' : ''}`}
              title="Samples loaded in the page: click one to show it as a track, click it again to remove the track">
              {runSamples.length ? `Samples · ${tracks.length}/${runSamples.length} shown` : '+ Add sample'}
            </button>
            {showPicker && (
              <div className={`absolute top-full ${pickerSide === 'right' ? 'right-0' : 'left-0'} mt-1 bg-white border-gray-200 border rounded-lg shadow-xl z-20 w-64 overflow-hidden`}>
                <input type="text" value={pickerSearch} onChange={e => setPickerSearch(e.target.value)}
                  placeholder="Search samples…" autoFocus className={`${t.inp} border-b w-full px-3 py-2 text-xs`} />
                <div className={`px-3 py-1 text-[10px] ${t.muted} border-b border-gray-100`}>click to show as a track · click again to remove</div>
                <div className="max-h-48 overflow-y-auto">
                  {filteredSamples.length === 0 ? (
                    <div className={`px-3 py-2 text-xs ${t.muted}`}>{runSamples.length ? 'No sample matches' : 'No samples loaded yet'}</div>
                  ) : filteredSamples.slice(0, 50).map(s => {
                    const idx = tracks.findIndex(x => x.sampleId === s.id);
                    const shown = idx >= 0;
                    return (
                      <button key={s.id} onClick={() => toggleSample(s)}
                        title={shown ? `Shown as track ${idx + 1}${idx === 0 ? ' (primary)' : ''} · click to remove it from the plot` : 'Click to show this sample as a track'}
                        className={`w-full text-left px-3 py-1.5 text-xs flex items-center gap-2 ${shown ? 'bg-indigo-50 text-indigo-900 font-semibold hover:bg-indigo-100' : `${t.text} hover:bg-gray-50`}`}>
                        <span className="inline-block w-2.5 h-2.5 rounded-sm shrink-0" style={shown ? { background: TRACK_COLORS[idx % TRACK_COLORS.length] } : { border: '1px solid #cbd5e1' }} />
                        <span className="flex-1 truncate">{s.name}</span>
                        {shown && <span className="text-indigo-600 text-[10px] font-medium">✓ shown</span>}
                      </button>
                    );
                  })}
                </div>
              </div>
            )}
          </div>}
          <div className="relative">
            <button onClick={e => openDropdown(e, 288, setShowGtexPicker)} className={`${t.btn} px-3 py-1 font-medium`}
              title="Add a GTEx tissue as a track: median junction read counts (arcs) and median exon reads per base (profile) over all samples of the tissue (GTEx v10, hg38)">+ GTEx tissue</button>
            {showGtexPicker && (() => {
              const q = gtexSearch.trim().toLowerCase();
              const all = gtexTissues || [];
              const match = (x: GtexTissue) => !q || x.name.toLowerCase().includes(q) || x.site.toLowerCase().includes(q) || x.id.toLowerCase().includes(q);
              const favs = gtexFavourites.map(id => all.find(x => x.id === id)).filter((x): x is GtexTissue => !!x && match(x));
              const rest = all.filter(x => match(x) && !gtexFavourites.includes(x.id));
              const loaded = new Set(gtexTracks.map(g => g.gtex?.tissue.id));
              const row = (x: GtexTissue) => (
                <button key={x.id} onClick={() => addGtexTissue(x)} disabled={loaded.has(x.id)}
                  className={`w-full text-left px-3 py-1.5 text-xs hover:bg-indigo-50 flex items-center gap-2 ${t.text} disabled:opacity-40`}>
                  <span className="inline-block w-2.5 h-2.5 rounded-sm" style={{ background: x.color }} />
                  <span className="flex-1 truncate">{x.name}</span>
                  <span className={`${t.muted} text-[10px]`}>n={x.samples}</span>
                </button>
              );
              return (
                <div className={`absolute top-full ${pickerSide === 'right' ? 'right-0' : 'left-0'} mt-1 bg-white border-gray-200 border rounded-lg shadow-xl z-20 w-72 overflow-hidden`}>
                  <input type="text" value={gtexSearch} onChange={e => setGtexSearch(e.target.value)} placeholder="Search tissues…" autoFocus className={`${t.inp} border-b w-full px-3 py-2 text-xs`} />
                  <div className="max-h-64 overflow-y-auto">
                    {gtexError ? <div className="px-3 py-2 text-xs text-red-600">{gtexError}</div>
                      : !gtexTissues ? <div className={`px-3 py-2 text-xs ${t.muted}`}>loading GTEx tissues…</div>
                      : <>
                        {favs.length > 0 && <div className={`px-3 pt-1.5 text-[10px] uppercase tracking-wide ${t.muted}`}>Favourites</div>}
                        {favs.map(row)}
                        {rest.length > 0 && <div className={`px-3 pt-1.5 text-[10px] uppercase tracking-wide ${t.muted}`}>All tissues</div>}
                        {rest.map(row)}
                        {favs.length + rest.length === 0 && <div className={`px-3 py-2 text-xs ${t.muted}`}>No tissue matches</div>}
                      </>}
                  </div>
                </div>
              );
            })()}
          </div>
          {onSnapshot && (
            <button onClick={takeSnapshot} disabled={snapshotState === 'busy'}
              className={`${t.btn} px-3 py-1 font-medium ${snapshotState === 'done' ? 'bg-green-50 border-green-300 text-green-700' : snapshotState === 'error' ? 'bg-red-50 border-red-300 text-red-700' : ''}`}
              title="Add a screenshot of this plot to the basket, together with the region, the options in effect and the outlier effect it was opened from">
              {snapshotState === 'busy' ? '📷 …' : snapshotState === 'done' ? '✓ Added to basket' : snapshotState === 'error' ? 'Screenshot failed' : '📷 Basket'}
            </button>
          )}
          <button onClick={exportSvg} className={`${t.btn} px-3 py-1 font-medium`} title="Export the plot as SVG (vector, publication-ready)">SVG</button>
          {!embedded && <button onClick={onClose} className={`${t.muted} text-2xl leading-none hover:text-red-400 px-2`}>&times;</button>}
        </div>
      </div>

      {/* Plot */}
      <div ref={containerRef} className="relative overflow-x-auto px-2 pb-2 pt-1"
        onMouseUp={handleMouseUp} onMouseLeave={handleMouseLeave}>
        <svg
          ref={svgRef}
          width={svgWidth}
          height={totalHeight}
          viewBox={`0 0 ${svgWidth} ${totalHeight}`}
          style={{ cursor: regionSelect ? 'col-resize' : dragging ? 'grabbing' : 'crosshair', userSelect: 'none', display: 'block' }}
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onDoubleClick={resetZoom}
          xmlns="http://www.w3.org/2000/svg"
          fontFamily={FONT}
        >
          <rect width={svgWidth} height={totalHeight} fill={INK.bg} />

          {/* Gene band */}
          {(() => {
            const a = scale.x(currentGeneStart), b = scale.x(currentGeneEnd);
            const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
            return right > left ? <rect x={left} y={RULER_H} width={right - left} height={legendY - RULER_H} fill={INK.geneBand} /> : null;
          })()}

          {/* Locus asked for by coordinates */}
          {locusMark && chromKey(locusMark.chrom) === chromKey(currentChrom) && (() => {
            const a = scale.x(locusMark.start), b = scale.x(locusMark.end);
            const left = Math.max(PLOT_LEFT, Math.min(a, b)), right = Math.min(plotRight, Math.max(a, b));
            if (right < PLOT_LEFT || left > plotRight) return null;
            const point = locusMark.end - locusMark.start <= 1;
            const label = point ? `${locusMark.chrom}:${(locusMark.start + 1).toLocaleString()}` : `${locusMark.chrom}:${(locusMark.start + 1).toLocaleString()}-${locusMark.end.toLocaleString()}`;
            const w = label.length * 5.6 + 10;
            const lx = Math.min(plotRight - w - 2, Math.max(PLOT_LEFT + 2, (left + right) / 2 - w / 2));
            return (
              <g pointerEvents="none" fontFamily={FONT}>
                {point
                  ? <line x1={(a + b) / 2} y1={RULER_H} x2={(a + b) / 2} y2={legendY} stroke={INK.select} strokeWidth={1.2} strokeDasharray="5 3" opacity={0.8} />
                  : <rect x={left} y={RULER_H} width={Math.max(1, right - left)} height={legendY - RULER_H} fill={withAlpha(INK.select, 0.08)} stroke={INK.select} strokeWidth={0.8} strokeDasharray="5 3" opacity={0.9} />}
                <rect x={lx} y={RULER_H - 40} width={w} height={14} rx={3} fill={INK.select} opacity={0.9} />
                <text x={lx + w / 2} y={RULER_H - 29.5} textAnchor="middle" fill="#fff" fontSize={9} fontWeight={600}>{label}</text>
              </g>
            );
          })()}

          {renderRuler()}
          {renderTranscript(transcriptY)}
          {knownPanelH > 0 && renderKnown(knownY)}
          {showSnps && renderSnps(snpY)}
          {showAllTx && renderAltTranscripts(altY)}
          {layouts.map(renderTrack)}
          {readsPlacements.map(p => <g key={`reads-${p.sid}`} transform={`translate(0, ${p.y})`}>{p.rt.el}</g>)}

          {/* Variant sites: stars in the strip above each sample's sashimi, guide lines through coverage and reads */}
          {readsPlacements.map(p => {
            const L = layouts.find(l => l.strip > 0 && l.track.sampleId === p.sid);
            if (!L) return null;
            const cy = L.yOff + TRACK_LABEL_H + L.strip / 2;
            const bottom = p.y + p.rt.height;
            return (
              <g key={`sites-${p.sid}`} fontFamily={FONT}>
                <text x={PLOT_LEFT + 8} y={cy + 3.5} fill={INK.faint} fontSize={8.5} letterSpacing={0.3}>VARIANT SITES</text>
                {p.rt.sites.map(st => {
                  const cx = scale.x(st.pos + 0.5);
                  if (cx < PLOT_LEFT || cx > plotRight) return null;
                  return (
                    <g key={`site${st.pos}${st.kind}`}>
                      <title>{siteLabel(st)}</title>
                      <line x1={cx} y1={cy + 7} x2={cx} y2={bottom} stroke={STAR_COLOR} strokeWidth={1} strokeDasharray="2 3" opacity={0.75} />
                      {knownSnp(st) && <circle cx={cx} cy={cy} r={9.5} fill="none" stroke={SNP_KNOWN_RING} strokeWidth={1.6} />}
                      <path d={starPath(cx, cy, 7)} fill={STAR_COLOR} stroke="#92400e" strokeWidth={0.8} />
                    </g>
                  );
                })}
              </g>
            );
          })}
          {renderKnownOverlay()}
          {renderLegend(legendY)}

          {tracks.length === 0 && (
            <text x={PLOT_LEFT + 8} y={RULER_H + 24} fill={INK.muted} fontSize={11}>{sampleId > 0 ? 'Loading coverage…' : 'No alignment yet: add BAM or CRAM files (with their index) to see coverage, junctions and reads here.'}</text>
          )}

          {/* Hover crosshair */}
          {hoverInfo && (
            <g data-export="skip" pointerEvents="none">
              <line x1={hoverInfo.x} y1={RULER_H} x2={hoverInfo.x} y2={tracksBottom - TRACK_GAP} stroke={INK.select} strokeWidth={1} strokeDasharray="3 3" opacity={0.7} />
              {(() => {
                const txt = `${currentChrom}:${(hoverInfo.pos + 1).toLocaleString()}`;
                const w = txt.length * 6.2 + 12;
                const x = Math.min(plotRight - w / 2, Math.max(PLOT_LEFT + w / 2, hoverInfo.x));
                return (
                  <g>
                    <rect x={x - w / 2} y={2} width={w} height={16} rx={3} fill={INK.select} />
                    <text x={x} y={13.5} textAnchor="middle" fill="#fff" fontSize={9.5} fontWeight={600}>{txt}</text>
                  </g>
                );
              })()}
            </g>
          )}

          {/* Region selection overlay (Ctrl+drag) */}
          {regionSelect && (() => {
            const x1 = Math.min(regionSelect.startX, regionSelect.currentX);
            const x2 = Math.max(regionSelect.startX, regionSelect.currentX);
            return <rect data-export="skip" x={x1} y={0} width={x2 - x1} height={totalHeight} fill={withAlpha(INK.select, 0.12)} stroke={INK.select} strokeWidth={1} strokeDasharray="4 2" />;
          })()}
        </svg>

        {/* HTML tooltip (crisper text than SVG, never exported) */}
        {hoverInfo && !popover && (
          <div className="pointer-events-none absolute z-20 min-w-[190px] rounded-md border border-gray-200 bg-white/95 px-2.5 py-1.5 text-xs shadow-lg"
            style={{ left: Math.min(hover!.px + 14, svgWidth - 200), top: Math.max(RULER_H + 4, hover!.py - 10) }}>
            <div className="font-mono text-[11px] text-gray-900 leading-5">
              {currentChrom}:{(hoverInfo.pos + 1).toLocaleString()}
              {hoverInfo.alt ? <span className="ml-2 font-semibold text-amber-700">{hoverInfo.alt.label}</span> : hoverInfo.cdna && <span className="ml-2 font-semibold text-indigo-700">{hoverInfo.cdna.label}</span>}
            </div>
            {hoverInfo.known.map(k => (
              <div key={k.v.id + (k.sample || '')} className="text-[10px] leading-4 mb-1 whitespace-pre-line border-l-2 pl-1.5" style={{ borderColor: KNOWN_VARIANT_COLORS[k.v.kind] }}>
                <span className="font-semibold text-gray-900">{k.sample ? `${k.sample} · ` : ''}{k.v.label}</span>
                {'\n'}{knownVariantTitle(k.v)}
                <span className="text-gray-400">{'\n'}{k.off ? 'click the arrow to go there' : 'click the mark to centre the view on it'}</span>
              </div>
            ))}
            {hoverInfo.snps.map(v => <div key={v.id + v.start} className="text-[10px] text-blue-700 leading-4 mb-0.5 whitespace-pre-line">{snpText(v)}</div>)}
            {hoverInfo.alt && <div className="text-[10px] text-amber-700 leading-4 mb-0.5">{hoverInfo.alt.id} ({hoverInfo.alt.biotype}) · {hoverInfo.alt.kind === 'cds' ? 'coding' : hoverInfo.alt.kind === 'intron' ? 'intronic' : hoverInfo.alt.kind.replace('utr', "UTR ")}</div>}
            {hoverInfo.cdna && tx && <div className="text-[10px] text-gray-500 leading-4 mb-0.5">{hoverInfo.alt ? `${hoverInfo.cdna.label} on ` : ''}{tx.transcriptId} ({modelKindLabel(tx)}) · {hoverInfo.cdna.kind === 'cds' ? 'coding' : hoverInfo.cdna.kind === 'intron' ? 'intronic' : hoverInfo.cdna.kind.replace('utr', "UTR ")}</div>}
            {hoverInfo.rows.map(r => (
              <div key={r.name} className="flex items-center gap-3 leading-5">
                <span className="inline-block h-2.5 w-2.5 rounded-sm" style={{ background: r.color }} />
                <span className="text-gray-700">{r.name}</span>
                <span className="ml-auto font-mono font-semibold text-gray-900">{r.depth.toLocaleString()}</span>
              </div>
            ))}
          </div>
        )}
        {cartoon && (
          <SpliceCartoon story={cartoonState.story} loading={cartoonState.loading} error={cartoonState.error} tx={cartoon.model}
            sampleName={cartoon.sample} sampleColor={cartoon.color} junctionLabel={cartoon.label} onClose={() => setCartoon(null)} />
        )}
        {popover && popoverContent && (
          <div className="absolute z-30 w-[540px] max-w-[95%] rounded-lg border border-gray-300 bg-white shadow-2xl text-xs"
            style={{ left: Math.min(popover.x + 12, Math.max(8, svgWidth - 552)), top: Math.max(RULER_H, popover.y + 12) }}
            onMouseDown={e => e.stopPropagation()}>
            <div className="flex items-start justify-between gap-2 px-3 pt-2">
              <div>
                <div className="font-semibold text-gray-900">{popoverContent.title}</div>
                {popoverContent.subtitle && <div className="text-gray-500">{popoverContent.subtitle}</div>}
              </div>
              <div className="flex items-center gap-2">
                {popoverContent.cartoon && (
                  <button onClick={() => { const c = popoverContent.cartoon!; const idx = tracks.findIndex(t => t.junctions.some(k => junctionKey(k) === junctionKey(c.j))); setCartoon({ j: c.j, model: c.model, label: c.label, color: TRACK_COLORS[Math.max(0, idx) % TRACK_COLORS.length], sample: tracks[Math.max(0, idx)]?.sampleName ?? '' }); }}
                    className="px-2 py-0.5 rounded-full bg-indigo-600 text-white text-[11px] font-semibold hover:bg-indigo-700" title="Animated cartoon: splicing, translation, NMD or protein consequence (experimental)">🎬 Cartoon</button>
                )}
                <button onClick={() => setPopover(null)} className="text-gray-400 hover:text-gray-700 text-base leading-none" title="Close (Esc)">×</button>
              </div>
            </div>
            {popoverContent.hgvs.length > 0 && (
              <div className="mx-3 mt-2 rounded bg-indigo-50 px-2 py-1.5 font-mono text-[11px] text-indigo-900">
                {popoverContent.hgvs.map((l, i) => <div key={i}>{l}</div>)}
              </div>
            )}
            {popoverContent.tables.map((tb, ti) => (
              <table key={ti} className="mx-3 my-2 w-[calc(100%-1.5rem)] border-collapse">
                {tb.caption && <caption className="text-left text-[10.5px] font-semibold text-gray-700 pb-0.5">{tb.caption}</caption>}
                <thead><tr>{tb.head.map(h => <th key={h} className="text-left text-[10px] font-semibold uppercase tracking-wide text-gray-500 py-1 pr-2">{h}</th>)}</tr></thead>
                <tbody>{tb.rows.map((r, i) => (
                  <tr key={i} className="border-t border-gray-100">{r.map((c, k) => <td key={k} className={`py-1 pr-2 ${k === 0 ? 'text-gray-800' : 'font-mono text-gray-900'}`}>{c}</td>)}</tr>
                ))}</tbody>
              </table>
            ))}
            {popoverContent.strip && popoverContent.strip.values.length > 1 && (() => {
              const st = popoverContent.strip!;
              const W = 400, H = 30, L = 8, R = 8;
              const max = Math.max(1.2, ...st.values.map(v => v.value));
              const x = (v: number) => L + (Math.min(v, max) / max) * (W - L - R);
              return (
                <div className="mx-3 mb-1">
                  <div className="text-[10px] text-gray-500">usage of every sample of the run (grey = controls, coloured = loaded tracks, bar = median)</div>
                  <svg width="100%" viewBox={`0 0 ${W} ${H}`} className="block">
                    <line x1={x(0)} y1={16} x2={x(max)} y2={16} stroke="#d1d5db" strokeWidth={1} />
                    {[0, 0.5, 1].map(t => <g key={t}><line x1={x(t)} y1={13} x2={x(t)} y2={19} stroke="#9ca3af" strokeWidth={1} /><text x={x(t)} y={28} textAnchor="middle" fontSize={8} fill="#6b7280">{t * 100}%</text></g>)}
                    {st.median != null && <line x1={x(st.median)} y1={6} x2={x(st.median)} y2={26} stroke="#374151" strokeWidth={1.5} />}
                    {st.values.filter(v => !v.color).map((v, i) => <circle key={`c${i}`} cx={x(v.value)} cy={16} r={3} fill="#9ca3af" opacity={0.55}><title>{`${v.name}: ${(v.value * 100).toFixed(0)}%`}</title></circle>)}
                    {st.values.filter(v => v.color).map((v, i) => <circle key={`t${i}`} cx={x(v.value)} cy={16} r={4.5} fill={v.color!} stroke="#fff" strokeWidth={1}><title>{`${v.name}: ${(v.value * 100).toFixed(0)}%`}</title></circle>)}
                  </svg>
                </div>
              );
            })()}
            {popoverContent.status && <div className="px-3 pb-1 text-[10px] text-indigo-600">{popoverContent.status}</div>}
            <div className="px-3 pb-2 text-[10px] text-gray-500">{popoverContent.note}</div>
          </div>
        )}
      </div>
      {showGroupsDialog && (
        <div className="fixed inset-0 z-50 bg-black/40 flex items-start justify-center p-6 overflow-y-auto" onMouseDown={() => setShowGroupsDialog(false)}>
          <div className="bg-white rounded-xl shadow-2xl border border-gray-200 w-full max-w-2xl text-gray-900" onMouseDown={e => e.stopPropagation()}>
            <div className="flex items-start justify-between gap-4 px-4 py-3 border-b border-gray-200">
              <div>
                <div className="font-bold text-sm">Sample groups</div>
                <div className="text-[11px] text-gray-500">Each group becomes one pooled track in the Groups view: coverage and junction reads summed over its samples, every arc labelled with a percentage instead of a read count. At each intron of the reference model the events using its donor or acceptor compete, and every arc shows its share of their reads, so the labels of one intron add up to 100 %: canonical C, alternative site n, pseudo-exon (A + B) / 2 on both arcs, exon skipping S, intron retention (R5 + R3) / 2 from the unspliced reads through both boundaries. A skipping arc shows 2·S over the totals of the two introns it spans, the rMATS value 2·S / (I₁ + I₂ + 2·S) when nothing else competes there. Tooltips also give each event against the canonical junction alone.</div>
              </div>
              <button onClick={() => setShowGroupsDialog(false)} className="text-gray-400 hover:text-gray-700 text-lg leading-none px-1" title="Close">×</button>
            </div>
            <div className="px-4 py-3 space-y-3 max-h-[60vh] overflow-y-auto">
              {groups.length === 0 && <div className="text-xs text-gray-500">No group yet. Create one and add samples to it; a sample belongs to one group at a time.</div>}
              {groups.map((g, gi) => {
                const color = TRACK_COLORS[gi % TRACK_COLORS.length];
                const nameOf = (sid: number) => runSamples.find(x => x.id === sid)?.name ?? tracks.find(x => x.sampleId === sid)?.sampleName ?? `#${sid}`;
                const free = runSamples.filter(x => !g.sampleIds.includes(x.id));
                return (
                  <div key={g.id} className="border border-gray-200 rounded-lg p-3">
                    <div className="flex items-center gap-2">
                      <span className="inline-block w-3 h-3 rounded-sm shrink-0" style={{ background: color }} />
                      <input value={g.name} onChange={e => renameGroup(g.id, e.target.value)} placeholder="Group name"
                        className={`${t.inp} border rounded px-2 py-1 text-sm font-semibold flex-1 min-w-0`} />
                      <span className="text-[11px] text-gray-500 whitespace-nowrap">{g.sampleIds.length} sample{g.sampleIds.length === 1 ? '' : 's'}</span>
                      <button onClick={() => deleteGroup(g.id)} className="text-xs text-red-600 hover:underline">delete</button>
                    </div>
                    <div className="flex flex-wrap items-center gap-1.5 mt-2">
                      {g.sampleIds.map(sid => (
                        <span key={sid} className="flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border bg-gray-50 border-gray-300">
                          {nameOf(sid)}
                          <button onClick={() => removeFromGroup(g.id, sid)} className="text-gray-400 hover:text-red-500" title="Remove from the group">×</button>
                        </span>
                      ))}
                      <select value="" onChange={e => { const id = parseInt(e.target.value); if (id) addToGroup(g.id, id); }}
                        className={`${t.inp} border rounded px-1 py-0.5 text-xs`} title="Add a sample loaded in the page to this group">
                        <option value="">+ add sample…</option>
                        {free.map(x => { const other = groups.find(o => o.id !== g.id && o.sampleIds.includes(x.id)); return <option key={x.id} value={x.id}>{x.name}{other ? ` (moves from ${other.name})` : ''}</option>; })}
                      </select>
                    </div>
                  </div>
                );
              })}
            </div>
            <div className="flex flex-wrap items-center gap-2 px-4 py-3 border-t border-gray-200">
              <button onClick={addGroup} className={`${t.btn} px-3 py-1 font-medium`}>+ New group</button>
              <span className="text-[11px] text-gray-500">{runSamples.length} sample{runSamples.length === 1 ? '' : 's'} loaded in the page</span>
              <span className="ml-auto flex gap-2">
                <button onClick={() => setShowGroupsDialog(false)} className={`${t.btn} px-3 py-1`}>Close</button>
                <button onClick={() => { setShowGroupsDialog(false); setViewMode('groups'); }} disabled={!groups.some(g => g.sampleIds.length)}
                  className="px-3 py-1 text-xs rounded bg-indigo-600 text-white disabled:opacity-40 hover:bg-indigo-700 font-medium">Show groups</button>
              </span>
            </div>
          </div>
        </div>
      )}
      <div className={`px-5 pb-2 text-[10.5px] ${t.muted}`}>
        Drag to pan · Ctrl+drag to zoom into a region · Ctrl+scroll to zoom around the cursor · double-click to reset · drag an arc vertically to untangle it · hover for c. positions · click an arc (HGVS, frame, share vs canonical) or an exon (depth-based usage) for details
      </div>
    </div>
  );

  if (embedded) return content;
  return (
    <div className="fixed inset-0 bg-black/70 backdrop-blur-sm flex items-start justify-center p-4 z-50 overflow-y-auto">
      {content}
    </div>
  );
}
