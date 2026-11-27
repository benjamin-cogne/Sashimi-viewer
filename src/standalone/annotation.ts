/**
 * Gene models from the user's annotation files (GFF3, GTF, BED 3–12, plain or gzipped), for a custom genome (any FASTA)
 * or ahead of RefSeq on a human build. Files are read in the page, line by line (a gzipped file is inflated as a stream,
 * never held whole as text), and indexed by contig: gene name, gene id, transcript id and ID search, overlap queries.
 *
 * - GFF3: features linked by ID / Parent (a feature may have several parents). A transcript is an mRNA, a transcript or
 *   any *RNA, or any feature with exon or CDS children; its exons are its `exon` children, else the union of its CDS,
 *   UTR and codon children; its CDS the span of its CDS children. A gene with exon or CDS children of its own (a
 *   prokaryote's gene → CDS) is its own transcript. A feature with neither parent nor children (a region of interest, a
 *   gene without transcripts) is a model of one exon. Attributes are URL-decoded; the name is Name, gene_name, gene,
 *   locus_tag, else ID. Lines after ##FASTA are not read.
 * - GTF: rows grouped by transcript_id under gene_id; exons from `exon` rows, else from CDS / UTR / codon rows; the CDS
 *   from the CDS rows and the stop codon; the name from gene_name, else gene_id.
 * - BED: one model per line, named by column 4 (else chrom:start-end); thickStart / thickEnd the CDS; BED12 blocks the
 *   exons.
 * Coordinates are kept 0-based half-open here; the outputs follow the shapes of ucsc.ts (1-based inclusive, except the
 * SvTranscript of a window, 0-based).
 *
 * A large file (a whole human GTF is ~1 GB of text) is better given bgzipped with its tabix index (.tbi or .csi;
 * `sort -k1,1 -k4,4n x.gtf | bgzip > x.gtf.gz && tabix -p gff x.gtf.gz`, `-p bed` for a BED): it is then never held
 * in memory. The genes of a window are read from the index when the window is drawn (the rows reaching past the
 * window widen the read to the genes' whole span), and the names to find genes by come from one pass over the file in
 * a worker (annotationScan.ts), while the windows already work.
 */
import type { TranscriptData, AllTranscripts, TranscriptModel, GeneModel, RegionHint } from '../components/sashimi/types';
import type { SvTranscript } from '../components/sashimi/svGenes';
import { resolveContig, isGeneticCode } from '../components/sashimi/geometry';
import { TabixIndexedFile } from '@gmod/tabix';
import { BlobFile } from 'generic-filehandle2';
import { readLines, isGzipped, formatOf, scanNames, own, SKIP_LONE, type AnnotationFormat, type GeneLoc, type NameIndex } from './annotationScan';
export { resolveContig };
export type { AnnotationFormat };

export const ANNOTATION_EXT = /\.(gff3?|gtf|bed)(\.gz)?$/i;
/** the tabix index of a bgzipped annotation file */
export const ANNOTATION_INDEX_EXT = /\.(gff3?|gtf|bed)\.gz\.(tbi|csi)$/i;
/** a gene's rows read from an index reach at most this far past the window asked */
const MAX_GENE_SPAN = 10_000_000;

interface ATx {
  id: string; name: string; chrom: string; strand: 1 | -1;
  /** 0-based half-open, sorted, disjoint */
  exons: [number, number][];
  cdsStart: number | null; cdsEnd: number | null;
  biotype: string; canonical: boolean;
  /** the genetic code of its CDS (transl_table), when one other than the standard is given */
  translTable?: number;
}
/** What a track shows of an entity besides its model: its type, colour, score, attributes, and no strand ('.'). */
interface Extra { type?: string; color?: string; score?: number; attrs?: string; unstranded?: boolean;
  /** drawn in the file's track only, never a gene to find or a neighbour (a TF binding site under a gene) */
  trackOnly?: boolean }
interface AGene extends Extra { key: string; id: string; name: string;
  /** other names it is found by (locus_tag, gene, Alias, gene_synonym, old_locus_tag) */
  aliases?: string[]; chrom: string; strand: 1 | -1; start: number; end: number; biotype: string; txs: ATx[] }

/**
 * An entity of an annotation file as its track draws it (0-based half-open, on the viewer's chromosome): a gene with
 * its transcripts, or any other feature (a GFF3 repeat_region, enhancer, region of interest; a BED line) as one
 * transcript of its blocks.
 */
export interface TrackItem {
  key: string; name: string; type: string;
  /** 0: no strand given ('.') */
  strand: 1 | -1 | 0;
  start: number; end: number;
  color?: string; score?: number; attrs?: string;
  /** the transcript drawn when the track is folded (pickTx: canonical tag, else longest CDS, else longest) */
  shown: number;
  txs: { id: string; name: string; biotype: string; strand: 1 | -1; exons: [number, number][]; cdsStart: number | null; cdsEnd: number | null }[];
}
/** One file's entities over a window. */
export interface TrackData { file: string; format: AnnotationFormat; items: TrackItem[] }

/** What a file held, for its chip. */
export interface AnnotationFileInfo {
  file: string; format: AnnotationFormat; genes: number; transcripts: number; lines: number; skipped: number;
  /** what reading the file had to make up for (exons without a Parent placed by position…), for its chip */
  notes?: string[];
  /** other names of its sequences (a GFF3's region rows: Name=, chromosome=, Alias=), to find them in the reference */
  aliases?: Record<string, string[]>;
  /** the sequences of its ##FASTA section, by name and length, when it carries them (bases kept by the source) */
  embedded?: { name: string; length: number }[];
  /** read through its tabix index (a bgzipped file with its .tbi / .csi): windows from the index, names from a pass */
  indexed?: { index: string; namesRead: number; namesReady: boolean; error?: string };
  /** contigs of the file, and those the reference does not hold (set once the reference's names are known) */
  contigs: string[]; unmatched: string[];
}

// ======================== parsers ========================

const decode = (s: string) => { try { return decodeURIComponent(s); } catch { return s; } };
const gffAttrs = (col: string): Map<string, string> => {
  const m = new Map<string, string>();
  for (const kv of col.split(';')) {
    const i = kv.indexOf('=');
    if (i > 0) m.set(kv.slice(0, i).trim(), decode(kv.slice(i + 1).trim()));
  }
  return m;
};
const gtfAttrs = (col: string): Map<string, string> => {
  const m = new Map<string, string>();
  for (const x of col.matchAll(/(\w+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?/g)) {
    const k = x[1], v = x[2] ?? x[3] ?? '';
    // repeated keys (tag "basic"; tag "Ensembl_canonical";): kept joined
    m.set(k, m.has(k) ? `${m.get(k)},${v}` : v);
  }
  return m;
};
const strandOf = (s: string): 1 | -1 => (s === '-' ? -1 : 1);
const CANONICAL = /Ensembl_canonical|MANE_Select|RefSeq Select|canonical/i;
const isCanonical = (a: Map<string, string>) => CANONICAL.test(a.get('tag') ?? '') || /^(1|true)$/i.test(a.get('is_canonical') ?? '') || /MANE/i.test(a.get('mane') ?? '');
const biotypeAttr = (a: Map<string, string>) => a.get('biotype') ?? a.get('gene_biotype') ?? a.get('transcript_biotype') ?? a.get('gene_type') ?? a.get('transcript_type');

/** Intervals merged into a sorted list of disjoint ones (0-based half-open). */
function merge(list: [number, number][]): [number, number][] {
  const out: [number, number][] = [];
  for (const e of [...list].sort((a, b) => a[0] - b[0])) {
    const last = out[out.length - 1];
    if (last && e[0] <= last[1]) last[1] = Math.max(last[1], e[1]);
    else out.push([e[0], e[1]]);
  }
  return out;
}
/**
 * The blocks of one transcript, sorted: abutting, duplicated or nested ones joined, but two that overlap in part kept
 * apart, which is how a programmed ribosomal frameshift is written (SARS-CoV-2 ORF1ab: one base read twice).
 */
function blocksOf(list: [number, number][]): [number, number][] {
  const out: [number, number][] = [];
  for (const e of [...list].sort((a, b) => a[0] - b[0] || a[1] - b[1])) {
    const last = out[out.length - 1];
    if (last && (e[0] === last[1] || e[1] <= last[1] || e[0] === last[0])) last[1] = Math.max(last[1], e[1]);
    else out.push([e[0], e[1]]);
  }
  return out;
}

interface Feat { id: string | null; parents: string[]; type: string; chrom: string; start: number; end: number; strand: 1 | -1; attrs: Map<string, string>; score?: number; raw?: string; unstranded?: boolean;
  /** the CDS phase (GFF3 column 8, GTF frame): bases to skip before the first whole codon */
  phase?: number }
const TX_TYPE = /^(mRNA|transcript|primary_transcript|pseudogenic_transcript|[a-zA-Z_]*RNA|[a-zA-Z_]*_transcript)$/;
const EXON_TYPE = /^(exon|pseudogenic_exon|noncoding_exon)$/;
const PART_TYPE = /^(CDS|five_prime_UTR|three_prime_UTR|UTR|5'UTR|3'UTR|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/;
/** parts of a transcript that are not drawn: its introns, and its polyA signals */
const NOT_ENTITY = /^(intron|polyA_site|polyA_signal_sequence)$/;
/** the other names a GFF3 gene is known by */
const aliasesOf = (a: Map<string, string> | undefined): string[] => (a ? [...new Set(['locus_tag', 'gene', 'gene_name', 'Alias', 'gene_synonym', 'old_locus_tag', 'gene_id'].flatMap(k => (a.get(k) ?? '').split(',')).map(x => x.trim()).filter(Boolean))] : []);
const geneNameOf = (a: Map<string, string>, fallback: string) => a.get('Name') ?? a.get('gene_name') ?? a.get('gene') ?? a.get('locus_tag') ?? a.get('ID') ?? fallback;

/** One transcript from its exon and part children (0-based half-open). */
function txFrom(id: string, name: string, chrom: string, strand: 1 | -1, exonKids: Feat[], partKids: Feat[], span: [number, number] | null, biotype: string | undefined, canonical: boolean): ATx {
  const cds = partKids.filter(k => k.type === 'CDS' || (k.type === 'stop_codon'));
  const cdsOnly = partKids.filter(k => k.type === 'CDS');
  let exons = blocksOf(exonKids.map(k => [k.start, k.end] as [number, number]));
  if (!exons.length) exons = blocksOf(partKids.map(k => [k.start, k.end] as [number, number]));
  if (!exons.length && span) exons = [span];
  let cdsStart = cdsOnly.length ? Math.min(...cds.map(k => k.start)) : null, cdsEnd = cdsOnly.length ? Math.max(...cds.map(k => k.end)) : null;
  // a CDS whose first piece (in the transcript's direction) starts with a phase of 1 or 2 is incomplete at its 5' end:
  // the bases before its first whole codon are left out, so that the codons are read in their frame
  if (cdsStart != null && cdsEnd != null) {
    const first = strand > 0 ? cdsOnly.reduce((a, b) => (b.start < a.start ? b : a)) : cdsOnly.reduce((a, b) => (b.end > a.end ? b : a));
    const ph = first.phase ?? 0;
    if (ph > 0 && cdsEnd - cdsStart > ph) { if (strand > 0) cdsStart += ph; else cdsEnd -= ph; }
  }
  const tt = Number(cdsOnly.map(k => k.attrs.get('transl_table')).find(Boolean) ?? NaN);
  return { id, name, chrom, strand, exons, cdsStart, cdsEnd, biotype: biotype ?? (cdsOnly.length ? 'protein_coding' : 'non_coding'), canonical, ...(isGeneticCode(tt) && tt !== 1 ? { translTable: tt } : {}) };
}

function geneOf(key: string, id: string, name: string, txs: ATx[], biotype?: string, extra: Extra = {}): AGene {
  const t0 = txs[0];
  return {
    ...extra, key, id, name, chrom: t0.chrom, strand: t0.strand, txs,
    start: Math.min(...txs.map(t => t.exons[0][0])), end: Math.max(...txs.map(t => t.exons[t.exons.length - 1][1])),
    biotype: biotype ?? (txs.some(t => t.cdsStart != null) ? 'protein_coding' : t0.biotype),
  };
}

type Lines = AsyncIterable<string> | Iterable<string>;
/** the attributes shown in a track's tooltip: at most this many characters */
const ATTRS_MAX = 400;
const scoreOf = (x: string | undefined) => { const v = x == null || x === '.' ? NaN : Number(x); return Number.isFinite(v) ? v : undefined; };
/** a colour as CSS: an itemRgb / GFF3 `color` "r,g,b", or a #hex / CSS name kept as given */
function cssColor(x: string | undefined): string | undefined {
  if (!x || x === '0' || x === '.') return undefined;
  const m = /^\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*$/.exec(x);
  if (m) return `rgb(${+m[1]},${+m[2]},${+m[3]})`;
  return /^(#[0-9a-f]{3,8}|[a-z]{3,20})$/i.test(x.trim()) ? x.trim() : undefined;
}
/** the type, colour, score and attributes of a GFF3 / GTF feature (its root's, for a gene) */
const extraOf = (f: Feat | null | undefined): Extra => (f ? {
  type: f.type, color: cssColor(f.attrs.get('color') ?? f.attrs.get('colour')), score: f.score,
  attrs: f.raw ? f.raw.slice(0, ATTRS_MAX) : undefined, unstranded: f.unstranded || undefined,
} : {});
const yieldToUi = (): Promise<void> => {
  const s = (globalThis as { scheduler?: { yield?: () => Promise<void> } }).scheduler;
  return s?.yield ? s.yield() : new Promise<void>(r => { setTimeout(r, 0); });
};

/** the sequences of a GFF3's ##FASTA section are kept up to this many bases in all (a bacterium, a plasmid, a virus) */
const EMBEDDED_FASTA_MAX = 200_000_000;

async function parseGff3(lines: Lines, info: AnnotationFileInfo, sequences?: Map<string, string>): Promise<AGene[]> {
  const feats: Feat[] = [];
  let inFasta = false, seqName = '', seqParts: string[] = [], seqTotal = 0;
  const endSeq = () => { if (sequences && seqName) sequences.set(seqName, seqParts.join('').toUpperCase()); seqName = ''; seqParts = []; };
  for await (const line of lines) {
    info.lines++;
    if (info.lines % 50000 === 0) await yieldToUi();
    if (inFasta) {
      // the sequences the file carries after ##FASTA: a reference when no FASTA is given
      if (!sequences) break;
      if (line.startsWith('>')) { endSeq(); seqName = line.slice(1).trim().split(/\s+/)[0]; continue; }
      seqTotal += line.length;
      if (seqTotal > EMBEDDED_FASTA_MAX) { sequences.clear(); seqName = ''; info.notes = [...(info.notes ?? []), 'its ##FASTA section is too large to keep: add the FASTA with its .fai instead']; break; }
      if (seqName) seqParts.push(line.trim());
      continue;
    }
    if (line.startsWith('##FASTA')) { inFasta = true; continue; }
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (c.length < 9) { info.skipped++; continue; }
    const s = parseInt(c[3]), e = parseInt(c[4]);
    if (!(s > 0) || !(e >= s)) { info.skipped++; continue; }
    const attrs = gffAttrs(c[8]);
    // a row spanning a whole sequence names it otherwise (NCBI: ID=NC_000001.11, Name=1, chromosome=1): aliases
    if (/^(region|chromosome|contig|supercontig|scaffold|source)$/i.test(c[2])) {
      const al = [attrs.get('Name'), attrs.get('chromosome'), attrs.get('Alias')?.split(',')].flat().filter((x): x is string => !!x && x !== decode(c[0]));
      if (al.length) { const k = decode(c[0]); info.aliases = { ...(info.aliases ?? {}), [k]: [...new Set([...(info.aliases?.[k] ?? []), ...al])] }; }
    }
    feats.push({ id: attrs.get('ID') ?? null, parents: (attrs.get('Parent') ?? '').split(',').filter(Boolean), type: c[2], chrom: decode(c[0]), start: s - 1, end: e, strand: strandOf(c[6]), attrs,
      score: scoreOf(c[5]), raw: decode(c[8]), unstranded: c[6] !== '+' && c[6] !== '-', phase: /^[012]$/.test(c[7]) ? +c[7] : undefined });
  }
  endSeq();
  const byId = new Map<string, Feat>();
  for (const f of feats) if (f.id && !byId.has(f.id)) byId.set(f.id, f);
  const kids = new Map<string, Feat[]>();
  for (const f of feats) for (const p of f.parents) { const l = kids.get(p); if (l) l.push(f); else kids.set(p, [f]); }
  const isLeaf = (f: Feat) => EXON_TYPE.test(f.type) || PART_TYPE.test(f.type);
  // Exons, CDS, UTRs and codons given without a Parent (or with one the file does not hold), as some tools write them:
  // each goes to the smallest gene or transcript holding it on its contig and strand, which then has them as its own
  // exons; those no feature holds make a model of their own, the lines sharing an ID (one discontinuous feature) together.
  const orphans = feats.filter(f => isLeaf(f) && !f.parents.some(p => byId.has(p)));
  if (orphans.length) {
    const holders = byChrom(feats.filter(f => !isLeaf(f) && !SKIP_LONE.test(f.type) && !NOT_ENTITY.test(f.type)));
    let placed = 0, auto = 0, guessed = 0;
    const alone = new Map<string, Feat[]>();
    const sameStrand = (x: Feat, o: Feat) => x.unstranded || o.unstranded || x.strand === o.strand;
    /** the genes or transcripts holding [start, end) of a contig on a strand, the smallest first */
    const holding = (o: Feat, start: number, end: number) => (holders.get(o.chrom)?.overlapping(start, end) ?? [])
      .filter(x => x.start <= start && x.end >= end && sameStrand(x, o))
      .sort((a, b) => (a.end - a.start) - (b.end - b.start) || Number(TX_TYPE.test(b.type)) - Number(TX_TYPE.test(a.type)));
    const place = (h: Feat, lines: Feat[]) => {
      if (!h.id) { h.id = `${h.type}@${h.chrom}:${h.start + 1}`; byId.set(h.id, h); }
      for (const o of lines) o.parents = [h.id];
      const l = kids.get(h.id); if (l) l.push(...lines); else kids.set(h.id, [...lines]);
      placed += lines.length;
    };
    const leaveAlone = (o: Feat) => { const k = `${o.id ?? `${o.type}@${o.chrom}:${o.start + 1}`}\t${o.chrom}\t${o.strand}`; const l = alone.get(k); if (l) l.push(o); else alone.set(k, [o]); };
    // 1. lines sharing an ID are one feature (a CDS in pieces): placed together, in the smallest holder of them all
    const groups = new Map<string, Feat[]>(), singles: Feat[] = [];
    for (const o of orphans) {
      if (!o.id) { singles.push(o); continue; }
      const k = `${o.id}\t${o.chrom}\t${o.strand}`;
      const l = groups.get(k); if (l) l.push(o); else groups.set(k, [o]);
    }
    for (const g of groups.values()) {
      if (g.length === 1) { singles.push(g[0]); continue; }
      const hs = holding(g[0], Math.min(...g.map(x => x.start)), Math.max(...g.map(x => x.end)));
      if (hs.length) { place(hs[0], g); if (hs.length > 1) guessed += g.length; } else for (const o of g) leaveAlone(o);
    }
    // 2. the other lines (exons with an ID each): the holder whose CDS already placed overlaps them, else the smallest
    // evidence: a CDS piece of the holder sharing the exon's bounds (2), one of them (1), overlapping it (0.5)
    const evidence = (h: Feat, o: Feat) => {
      let best = 0;
      for (const k of h.id ? kids.get(h.id) ?? [] : []) {
        if (k.type !== 'CDS' || k.start >= o.end || k.end <= o.start) continue;
        best = Math.max(best, (k.start === o.start ? 1 : 0) + (k.end === o.end ? 1 : 0) || 0.5);
      }
      return best;
    };
    for (const o of singles) {
      const hs = holding(o, o.start, o.end);
      if (!hs.length) { leaveAlone(o); continue; }
      const scored = hs.map(h => ({ h, e: evidence(h, o) })).sort((x, y) => y.e - x.e);   // stable: the smallest first among equals
      place(scored[0].h, [o]);
      if (scored.length > 1 && scored[1].e === scored[0].e) guessed++;
    }
    // a model of its own for each group no feature holds: a feature standing for the group, its lines its children
    for (const group of alone.values()) {
      const g0 = group[0], id = `${g0.id ?? g0.type}@${g0.chrom}:${g0.start + 1}`;
      const f: Feat = { ...g0, id, parents: [], type: g0.type === 'CDS' ? 'CDS' : 'transcript', start: Math.min(...group.map(x => x.start)), end: Math.max(...group.map(x => x.end)) };
      for (const x of group) x.parents = [id];
      feats.push(f); byId.set(id, f); kids.set(id, group); auto++;
    }
    info.notes = [...(info.notes ?? []), placed ? `${placed.toLocaleString('en-US')} exon / CDS line${placed === 1 ? '' : 's'} without a Parent placed in the gene holding ${placed === 1 ? 'it' : 'them'}${guessed ? ` (${guessed.toLocaleString('en-US')} where genes overlap on one strand, placed by best guess: a Parent= on each line makes it certain)` : ''}` : '',
      auto ? `${auto.toLocaleString('en-US')} group${auto === 1 ? '' : 's'} of exon / CDS lines without a Parent and outside any gene drawn as model${auto === 1 ? '' : 's'} of their own` : ''].filter(Boolean);
  }
  const childrenOf = (f: Feat) => (f.id ? kids.get(f.id) ?? [] : []);
  const isTx = (f: Feat) => TX_TYPE.test(f.type) || childrenOf(f).some(k => EXON_TYPE.test(k.type) || k.type === 'CDS');
  // genes: features holding transcripts, or transcripts without a parent
  const genes = new Map<string, { gene: Feat | null; txs: ATx[] }>();
  /** lone features without an ID, by their key */
  const lone = new Map<string, Feat>();
  /** features under a gene that are neither transcripts nor their parts (a TF binding site, a promoter): track only */
  const trackOnly = new Set<string>();
  const addTx = (geneKey: string, gene: Feat | null, tx: ATx) => { const g = genes.get(geneKey); if (g) g.txs.push(tx); else genes.set(geneKey, { gene, txs: [tx] }); };
  for (const f of feats) {
    if (!isTx(f) || isLeaf(f)) continue;
    const ks = childrenOf(f);
    const exonKids = ks.filter(k => EXON_TYPE.test(k.type)), partKids = ks.filter(k => PART_TYPE.test(k.type));
    const parentGene = f.parents.map(p => byId.get(p)).find((p): p is Feat => !!p);
    // a gene whose children are its exons / CDS (no transcript level) is its own transcript
    const gene = parentGene && !isTx(parentGene) ? parentGene : parentGene && parentGene !== f ? parentGene : null;
    const gf = gene ?? f;
    const name = geneNameOf(gf.attrs, gf.id ?? `${f.chrom}:${f.start + 1}-${f.end}`);
    const tx = txFrom(f.id ?? `${name}.${f.start}`, f.attrs.get('Name') ?? f.attrs.get('transcript_id') ?? f.id ?? name, f.chrom, f.strand, exonKids, partKids, [f.start, f.end], biotypeAttr(f.attrs) ?? (gene ? biotypeAttr(gene.attrs) : undefined), isCanonical(f.attrs));
    if (!tx.exons.length) continue;
    addTx(gf.id ?? `${name}@${f.chrom}:${f.start}`, gene, tx);
  }
  // features that gave no transcript: those without a parent (a region of interest, a gene whose children are not
  // transcripts) are a model of one block; those under a gene (not its transcripts) show in the track only
  for (const f of feats) {
    if (SKIP_LONE.test(f.type) || NOT_ENTITY.test(f.type) || isLeaf(f) || isTx(f)) continue;
    const name = geneNameOf(f.attrs, `${f.chrom}:${f.start + 1}-${f.end}`);
    const key = f.id ?? `${name}@${f.chrom}:${f.start}`;
    if (genes.has(key)) continue;
    const under = f.parents.some(p => byId.has(p));
    if (under) trackOnly.add(key);
    if (!f.id) lone.set(key, f);
    addTx(key, f, { id: f.id ?? name, name, chrom: f.chrom, strand: f.strand, exons: [[f.start, f.end]], cdsStart: null, cdsEnd: null, biotype: biotypeAttr(f.attrs) ?? f.type, canonical: false });
  }
  const out: AGene[] = [];
  for (const [key, { gene, txs }] of genes) {
    const name = gene ? geneNameOf(gene.attrs, key) : txs[0].name;
    // a gene, or a transcript with no gene above it, or a lone feature: the track shows its type and attributes
    const root = gene ?? byId.get(key) ?? null;
    const g = geneOf(key, gene?.id ?? key, name, txs, gene ? biotypeAttr(gene.attrs) : undefined, { ...extraOf(root ?? lone.get(key)), ...(trackOnly.has(key) ? { trackOnly: true } : {}) });
    const al = aliasesOf((gene ?? root ?? lone.get(key))?.attrs).filter(x => x !== name && x !== g.id);
    if (al.length) g.aliases = al;
    out.push(g);
  }
  return out;
}

async function parseGtf(lines: Lines, info: AnnotationFileInfo): Promise<AGene[]> {
  // transcripts and genes are told apart by their contig too: one ID placed on several (chrX / chrY PAR, alt contigs)
  const txRows = new Map<string, { id: string; raw?: string; gene: string; chrom: string; strand: 1 | -1; exon: Feat[]; part: Feat[]; attrs: Map<string, string>; span: [number, number] | null }>();
  const geneAttrs = new Map<string, { attrs: Map<string, string>; f: Feat }>();
  for await (const line of lines) {
    info.lines++;
    if (info.lines % 50000 === 0) await yieldToUi();
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (c.length < 9) { info.skipped++; continue; }
    const s = parseInt(c[3]), e = parseInt(c[4]);
    if (!(s > 0) || !(e >= s)) { info.skipped++; continue; }
    const attrs = gtfAttrs(c[8]);
    const at = `\t${c[0]}`;
    const f: Feat = { id: null, parents: [], type: c[2], chrom: c[0], start: s - 1, end: e, strand: strandOf(c[6]), attrs, score: scoreOf(c[5]), raw: c[8], unstranded: c[6] !== '+' && c[6] !== '-', phase: /^[012]$/.test(c[7]) ? +c[7] : undefined };
    const gid = attrs.get('gene_id') ?? attrs.get('gene_name') ?? '';
    const tid = attrs.get('transcript_id');
    if (c[2] === 'gene') { if (gid && !geneAttrs.has(gid + at)) geneAttrs.set(gid + at, { attrs, f }); continue; }
    if (!tid) { if (gid && !geneAttrs.has(gid + at)) geneAttrs.set(gid + at, { attrs, f }); continue; }
    let t = txRows.get(tid + at);
    if (!t) txRows.set(tid + at, t = { id: tid, raw: c[8], gene: (gid || tid) + at, chrom: c[0], strand: f.strand, exon: [], part: [], attrs, span: null });
    if (c[2] === 'transcript') { t.attrs = new Map([...t.attrs, ...attrs]); t.span = [f.start, f.end]; }
    else if (EXON_TYPE.test(c[2])) t.exon.push(f);
    else if (PART_TYPE.test(c[2]) || /^(UTR|five_prime_utr|three_prime_utr)$/i.test(c[2])) t.part.push(f);
  }
  const byGene = new Map<string, ATx[]>();
  for (const t of txRows.values()) {
    const tid = t.id;
    const tx = txFrom(tid, t.attrs.get('transcript_name') ?? tid, t.chrom, t.strand, t.exon, t.part, t.span, biotypeAttr(t.attrs), isCanonical(t.attrs));
    if (!tx.exons.length) continue;
    const l = byGene.get(t.gene); if (l) l.push(tx); else byGene.set(t.gene, [tx]);
  }
  const out: AGene[] = [];
  for (const [gk, txs] of byGene) {
    const gid = gk.slice(0, gk.lastIndexOf('\t')), key = `${gid}@${txs[0].chrom}`;
    const ga = geneAttrs.get(gk)?.attrs, ta = txRows.get(`${txs[0].id}\t${txs[0].chrom}`)?.attrs;
    const gf = geneAttrs.get(gk)?.f;
    out.push(geneOf(key, gid, ga?.get('gene_name') ?? ta?.get('gene_name') ?? gid, txs, ga ? biotypeAttr(ga) : undefined,
      { ...extraOf(gf), type: 'gene', attrs: (gf?.raw ?? txRows.get(`${txs[0].id}\t${txs[0].chrom}`)?.raw)?.slice(0, ATTRS_MAX) }));
  }
  // genes without transcript rows: one exon over the gene
  for (const [gk, { attrs, f }] of geneAttrs) if (!byGene.has(gk)) {
    const gid = gk.slice(0, gk.lastIndexOf('\t'));
    const name = attrs.get('gene_name') ?? gid;
    out.push(geneOf(`${gid}@${f.chrom}`, gid, name, [{ id: gid, name, chrom: f.chrom, strand: f.strand, exons: [[f.start, f.end]], cdsStart: null, cdsEnd: null, biotype: biotypeAttr(attrs) ?? 'gene', canonical: false }], undefined, extraOf(f)));
  }
  return out;
}

async function parseBed(lines: Lines, info: AnnotationFileInfo): Promise<AGene[]> {
  const out: AGene[] = [];
  for await (const line of lines) {
    info.lines++;
    if (info.lines % 50000 === 0) await yieldToUi();
    if (!line.trim() || /^(track|browser)\b/.test(line) || line.startsWith('#')) continue;
    const c = line.split(/\t/).length >= 3 ? line.split('\t') : line.trim().split(/\s+/);
    const s = parseInt(c[1]), e = parseInt(c[2]);
    if (!(s >= 0) || !(e > s)) { info.skipped++; continue; }
    const chrom = c[0], name = c[3] && c[3] !== '.' ? c[3] : `${chrom}:${s + 1}-${e}`, strand = strandOf(c[5] ?? '+');
    let exons: [number, number][] = [[s, e]];
    const n = parseInt(c[9]);
    if (n > 0 && c[10] && c[11]) {
      const sizes = c[10].split(',').filter(Boolean).map(Number), starts = c[11].split(',').filter(Boolean).map(Number);
      if (sizes.length === n && starts.length === n) exons = merge(starts.map((st, i) => [s + st, s + st + sizes[i]] as [number, number]));
    }
    const ts = parseInt(c[6]), te = parseInt(c[7]);
    const coding = ts >= s && te > ts && te <= e;
    const key = `${name}@${chrom}:${s}`;
    const extra: Extra = { type: 'bed', color: cssColor(c[8]), score: scoreOf(c[4]), unstranded: c[5] !== '+' && c[5] !== '-' || undefined,
      attrs: c.length > 12 ? c.slice(12).join('\t').slice(0, ATTRS_MAX) : undefined };
    out.push(geneOf(key, name, name, [{ id: name, name, chrom, strand, exons, cdsStart: coding ? ts : null, cdsEnd: coding ? te : null, biotype: coding ? 'protein_coding' : 'feature', canonical: false }], undefined, extra));
  }
  return out;
}

// ======================== sources: a file in memory, or read through its index ========================

/** the genes of some lines, their strings copied out of the text read (`own`) so that it can go */
async function parse(format: AnnotationFormat, lines: Lines, info: AnnotationFileInfo, sequences?: Map<string, string>): Promise<AGene[]> {
  const genes = await (format === 'gtf' ? parseGtf(lines, info) : format === 'bed' ? parseBed(lines, info) : parseGff3(lines, info, sequences));
  const chroms = new Map<string, string>();
  const chrom = (c: string) => { let o = chroms.get(c); if (o == null) chroms.set(c, o = own(c)); return o; };
  for (const g of genes) {
    g.key = own(g.key); g.id = own(g.id); g.name = own(g.name); g.chrom = chrom(g.chrom); g.biotype = own(g.biotype);
    if (g.type) g.type = own(g.type);
    if (g.color) g.color = own(g.color);
    if (g.attrs) g.attrs = own(g.attrs);
    if (g.aliases) g.aliases = g.aliases.map(own);
    for (const t of g.txs) { t.id = own(t.id); t.name = own(t.name); t.chrom = chrom(t.chrom); t.biotype = own(t.biotype); }
  }
  return genes;
}
/** the names a gene is found by: its name, ID, its transcripts' IDs and names (upper case) */
const keysOf = (g: AGene) => new Set([g.name, g.id, ...(g.aliases ?? []), ...g.txs.flatMap(t => [t.id, t.name])].map(k => k.toUpperCase()));

/** Intervals sorted by start with the running maximum of their ends: what overlaps a window, by bisection. */
class Intervals<T extends { start: number; end: number }> {
  private items: T[]; private maxEnd: number[] = [];
  constructor(items: T[]) {
    this.items = items.sort((a, b) => a.start - b.start);
    let m = 0;
    for (const g of this.items) { m = Math.max(m, g.end); this.maxEnd.push(m); }
  }
  overlapping(start: number, end: number): T[] {
    let lo = 0, hi = this.items.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (this.items[mid].start < end) lo = mid + 1; else hi = mid; }
    const out: T[] = [];
    for (let i = lo - 1; i >= 0 && this.maxEnd[i] > start; i--) if (this.items[i].end > start) out.push(this.items[i]);
    return out.reverse();
  }
}
const byChrom = <T extends { chrom: string; start: number; end: number }>(items: Iterable<T>) => {
  const per = new Map<string, T[]>();
  for (const g of items) { const l = per.get(g.chrom); if (l) l.push(g); else per.set(g.chrom, [g]); }
  return new Map([...per].map(([c, l]) => [c, new Intervals(l)]));
};

/** One annotation file, as the page asks it. Contigs are the file's spelling. */
interface GeneSource {
  readonly info: AnnotationFileInfo;
  /** genes known by this name (upper case); several when the name is on several contigs */
  byName(q: string): Promise<AGene[]>;
  /** genes of a contig overlapping [start, end) 0-based, whole */
  overlapping(chrom: string, start: number, end: number): Promise<AGene[]>;
  /** [name searched (upper case), gene name shown] of every gene, for the suggestions; none while not read */
  nameEntries(): Iterable<[string, string]>;
  /** the bases of a GFF3's ##FASTA section, by sequence name */
  readonly sequences?: Map<string, string>;
}

/** A file read whole: every gene in memory. */
class FileIndex implements GeneSource {
  private chroms: Map<string, Intervals<AGene>>;
  private names = new Map<string, AGene[]>();
  constructor(readonly info: AnnotationFileInfo, genes: AGene[], readonly sequences?: Map<string, string>) {
    if (sequences?.size) info.embedded = [...sequences].map(([name, seq]) => ({ name, length: seq.length }));
    for (const g of genes) if (!g.trackOnly) for (const k of keysOf(g)) { const l = this.names.get(k); if (l) l.push(g); else this.names.set(k, [g]); }
    this.chroms = byChrom(genes);
    const found = genes.filter(g => !g.trackOnly);
    info.genes = found.length; info.transcripts = found.reduce((n, g) => n + g.txs.length, 0); info.contigs = [...this.chroms.keys()];
  }
  async byName(q: string) { return this.names.get(q) ?? []; }
  async overlapping(chrom: string, start: number, end: number) { return this.chroms.get(chrom)?.overlapping(start, end) ?? []; }
  *nameEntries(): Iterable<[string, string]> { for (const [k, gs] of this.names) yield [k, gs[0].name]; }
}

/** How a name index is built: in a worker (annotationWorkerClient.ts), or here. */
export type NameScanner = (file: Blob, gz: boolean, format: AnnotationFormat, onProgress?: (fraction: number) => void) => Promise<NameIndex>;

/** A bgzipped file with its tabix index: the genes of a window read when asked, names from one pass over the file. */
class TabixSource implements GeneSource {
  private tbx: TabixIndexedFile;
  private names: Map<string, number[]> | null = null;
  private locs: GeneLoc[] = [];
  private locsByChrom = new Map<string, Intervals<GeneLoc & { i: number }>>();
  /** the genes of the last windows read (key chrom:start-end) */
  private cache = new Map<string, Promise<AGene[]>>();
  constructor(readonly info: AnnotationFileInfo, file: File, index: File) {
    const fh = new BlobFile(index);
    this.tbx = new TabixIndexedFile({ filehandle: new BlobFile(file), ...(/\.csi$/i.test(index.name) ? { csiFilehandle: fh } : { tbiFilehandle: fh }) });
  }
  /** the index's sequence names as the file writes them (GFF3 escapes some characters), by their decoded name */
  private raw = new Map<string, string>();
  async open(): Promise<void> {
    const meta = await this.tbx.getMetadata();
    const dec = (n: string) => (this.info.format === 'gff3' ? (() => { try { return decodeURIComponent(n); } catch { return n; } })() : n);
    for (const n of meta.refIdToName) this.raw.set(dec(n), n);
    this.info.contigs = [...this.raw.keys()];
    // the other names of each sequence (a GFF3's region row, at its start: Name=1, chromosome=1), read from the index
    if (this.info.format === 'gff3' && this.raw.size <= 2000) {
      const aliases: Record<string, string[]> = {};
      await Promise.all([...this.raw].map(async ([name, raw]) => {
        const rows: string[] = [];
        await this.tbx.getLines(raw, 0, 1, { lineCallback: l => { rows.push(l); } }).catch(() => undefined);
        for (const r of rows) {
          const c = r.split('\t');
          if (!/^(region|chromosome|contig|supercontig|scaffold|source)$/i.test(c[2] ?? '')) continue;
          const a = new Map((c[8] ?? '').split(';').map(kv => { const i = kv.indexOf('='); return [kv.slice(0, i).trim(), dec(kv.slice(i + 1).trim())] as [string, string]; }));
          const al = [a.get('Name'), a.get('chromosome'), ...(a.get('Alias')?.split(',') ?? [])].filter((x): x is string => !!x && x !== name);
          if (al.length) aliases[name] = [...new Set(al)];
        }
      }));
      if (Object.keys(aliases).length) this.info.aliases = aliases;
    }
  }
  /** builds the name index; `onChange` is told as it goes */
  async readNames(file: File, scan: NameScanner, onChange: () => void): Promise<void> {
    const ix = this.info.indexed!;
    try {
      const idx = await scan(file, true, this.info.format, f => { ix.namesRead = f; onChange(); });
      this.locs = idx.locs;
      const names = new Map<string, number[]>();
      for (const [k, i] of idx.keys) { const l = names.get(k); if (l) l.push(i); else names.set(k, [i]); }
      this.names = names;
      this.locsByChrom = byChrom(idx.locs.map((l, i) => ({ ...l, i })));
      Object.assign(this.info, { genes: idx.genes, transcripts: idx.transcripts, lines: idx.lines });
      ix.namesRead = 1; ix.namesReady = true;
    } catch (e) {
      ix.error = (e as Error)?.message || String(e);
    }
    onChange();
  }
  private async lines(chrom: string, start: number, end: number): Promise<string[]> {
    const out: string[] = [];
    await this.tbx.getLines(this.raw.get(chrom) ?? chrom, Math.max(0, start), end, { lineCallback: l => { out.push(l); } });
    return out;
  }
  /** the span a gene's rows cover: those of the rows read (gene and transcript rows give it whole) and the name index's */
  private reach(chrom: string, start: number, end: number, rows: string[]): [number, number] {
    let lo = start, hi = end;
    for (const g of this.locsByChrom.get(chrom)?.overlapping(start, end) ?? []) { lo = Math.min(lo, g.start); hi = Math.max(hi, g.end); }
    if (this.info.format !== 'bed') for (const r of rows) {
      const c = r.split('\t', 5);
      if (SKIP_LONE.test(c[2] ?? '')) continue;
      const s = +c[3] - 1, e = +c[4];
      if (s >= 0 && e > s) { lo = Math.min(lo, s); hi = Math.max(hi, e); }
    }
    return [Math.max(lo, start - MAX_GENE_SPAN), Math.min(hi, end + MAX_GENE_SPAN)];
  }
  overlapping(chrom: string, start: number, end: number): Promise<AGene[]> {
    if (!this.info.contigs.includes(chrom)) return Promise.resolve([]);
    const key = `${chrom}:${start}-${end}`;
    let p = this.cache.get(key);
    if (!p) {
      p = (async () => {
        let rows = await this.lines(chrom, start, end);
        let [lo, hi] = this.reach(chrom, start, end, rows);
        // the rows of the genes reaching past the window, until no row read reaches further (twice at most)
        for (let k = 0; k < 2 && (lo < start || hi > end); k++) {
          rows = await this.lines(chrom, lo, hi);
          const [lo2, hi2] = this.reach(chrom, lo, hi, rows);
          if (lo2 >= lo && hi2 <= hi) break;
          [lo, hi] = [lo2, hi2];
        }
        const scratch: AnnotationFileInfo = { ...this.info, lines: 0, skipped: 0 };
        const genes = await parse(this.info.format, rows, scratch);
        return genes.filter(g => g.chrom === chrom && g.end > start && g.start < end);
      })();
      this.cache.set(key, p);
      p.catch(() => this.cache.delete(key));
      if (this.cache.size > 64) this.cache.delete(this.cache.keys().next().value!);
    }
    return p;
  }
  async byName(q: string): Promise<AGene[]> {
    const hits = this.names?.get(q);
    if (!hits) return [];
    const out: AGene[] = [];
    for (const i of hits) {
      const l = this.locs[i];
      const genes = await this.overlapping(l.chrom, l.start, l.end);
      // the gene the name stands for: known by it, else the one named as the index says, else the one covering most of it
      const ov = (g: AGene) => Math.min(g.end, l.end) - Math.max(g.start, l.start);
      const real = genes.filter(x => !x.trackOnly);
      const g = real.find(x => keysOf(x).has(q)) ?? real.find(x => x.name === l.name) ?? [...real].sort((a, b) => ov(b) - ov(a))[0];
      if (g && !out.includes(g)) out.push(g);
    }
    return out;
  }
  *nameEntries(): Iterable<[string, string]> {
    if (!this.names) return;
    for (const [k, is] of this.names) yield [k, this.locs[is[0]].name];
  }
}

// ======================== the annotation of a page: its files, in order ========================

const cdsLen = (t: ATx) => (t.cdsStart == null || t.cdsEnd == null ? 0 : t.exons.reduce((n, [a, b]) => n + Math.max(0, Math.min(b, t.cdsEnd!) - Math.max(a, t.cdsStart!)), 0));
const txLen = (t: ATx) => t.exons.reduce((n, [a, b]) => n + b - a, 0);
/** The model shown: a canonical / MANE tag, else the longest CDS, else the longest transcript. */
function pickTx(g: AGene): { tx: ATx; kind: 'canonical' | 'longest' } {
  const tagged = g.txs.filter(t => t.canonical).sort((a, b) => cdsLen(b) - cdsLen(a) || txLen(b) - txLen(a))[0];
  if (tagged) return { tx: tagged, kind: 'canonical' };
  const coding = g.txs.filter(t => cdsLen(t) > 0);
  if (coding.length) return { tx: coding.reduce((b, t) => (cdsLen(t) > cdsLen(b) ? t : b)), kind: 'longest' };
  return { tx: g.txs.reduce((b, t) => (txLen(t) > txLen(b) ? t : b)), kind: 'longest' };
}

export class Annotation {
  readonly files: GeneSource[];
  private toRef = new Map<string, string>();
  /** told when a file's state changes (an index's names read) */
  onChange: (() => void) | null = null;
  constructor(files: GeneSource[]) { this.files = files; }

  get infos(): AnnotationFileInfo[] { return this.files.map(f => f.info); }
  get label(): string { return this.files.map(f => f.info.file).join(', '); }

  /**
   * The names the file's contigs are shown under: matched to the reference's (a custom genome), or spelled by `spell`
   * (a human build: chr-prefixed); with neither, as the file has them.
   */
  setReferenceNames(names: string[] | null, spell?: (c: string) => string): void {
    this.toRef.clear();
    for (const f of this.files) {
      f.info.unmatched = [];
      for (const c of f.info.contigs) {
        const r = names?.length ? resolveContig(names, c) ?? (f.info.aliases?.[c] ?? []).map(a => resolveContig(names, a)).find(Boolean) ?? null : spell ? spell(c) : c;
        if (r) this.toRef.set(c, r); else f.info.unmatched.push(c);
      }
    }
  }
  private refChrom(c: string): string { return this.toRef.get(c) ?? c; }
  /** the file contig names a viewer chromosome stands for */
  private fileChroms(f: GeneSource, chrom: string): string[] {
    return f.info.contigs.filter(c => this.refChrom(c) === chrom || resolveContig([c], chrom));
  }

  /** a gene by gene name, gene id, transcript id or ID (case-insensitive); the copy nearest the hint when several */
  async find(query: string, hint?: RegionHint): Promise<{ gene: AGene; file: GeneSource; tx?: ATx } | null> {
    const q = query.trim().toUpperCase();
    for (const f of this.files) {
      let list = await f.byName(q);
      // where the gene is known (a view reopened, a neighbour clicked): its window, whose genes are read whatever the
      // name index (an indexed file's may still be being read)
      if (!list.length && hint) list = (await Promise.all(this.fileChroms(f, hint.chrom).map(c => f.overlapping(c, hint.start - 1, hint.end)))).flat().filter(g => !g.trackOnly && keysOf(g).has(q));
      if (!list.length) continue;
      let gene = list[0];
      if (hint && list.length > 1) {
        const d = (g: AGene) => (this.refChrom(g.chrom) === hint.chrom || resolveContig([g.chrom], hint.chrom) ? Math.abs((g.start + g.end) / 2 - (hint.start + hint.end) / 2) : Infinity);
        gene = list.reduce((b, g) => (d(g) < d(b) ? g : b));
      }
      // a transcript asked by its own id or name (not the gene's, which a one-transcript gene or a BED line shares)
      const tx = gene.txs.length > 1 && q !== gene.name.toUpperCase() && q !== gene.id.toUpperCase()
        ? gene.txs.find(t => t.id.toUpperCase() === q || t.name.toUpperCase() === q) : undefined;
      return { gene, file: f, tx };
    }
    return null;
  }

  async transcript(query: string, hint?: RegionHint): Promise<TranscriptData | null> {
    const hit = await this.find(query, hint);
    if (!hit) return null;
    const { gene: g, file } = hit;
    const { tx: t, kind } = hit.tx ? { tx: hit.tx, kind: 'chosen' as const } : pickTx(g);
    const n = t.exons.length;
    return {
      gene_name: g.name, transcript_id: t.id, translation_id: null, is_mane_select: false, model_kind: kind, biotype: t.biotype,
      source: 'file', source_label: file.info.file, chrom: this.refChrom(t.chrom), strand: t.strand,
      start: t.exons[0][0] + 1, end: t.exons[n - 1][1],
      exons: t.exons.map(([a, b], i) => ({ start: a + 1, end: b, rank: t.strand > 0 ? i + 1 : n - i })),
      cds_start: t.cdsStart != null ? t.cdsStart + 1 : null, cds_end: t.cdsEnd, ...(t.translTable ? { transl_table: t.translTable } : {}),
    };
  }

  async allTranscripts(query: string, hint?: RegionHint): Promise<AllTranscripts | null> {
    const hit = await this.find(query, hint);
    if (!hit) return null;
    const { gene: g, file } = hit, picked = pickTx(g).tx;
    const models = g.txs.map((t): TranscriptModel => ({
      id: t.id, name: t.name, source: 'file', biotype: t.biotype, start: t.exons[0][0] + 1, end: t.exons[t.exons.length - 1][1], strand: t.strand,
      exons: t.exons.map(([a, b]) => ({ start: a + 1, end: b })), cds_start: t.cdsStart != null ? t.cdsStart + 1 : null, cds_end: t.cdsEnd,
      is_mane: false, is_canonical: t === picked,
    })).sort((a, b) => Number(b.is_canonical) - Number(a.is_canonical) || a.id.localeCompare(b.id, undefined, { numeric: true }));
    return { gene_name: g.name, chrom: this.refChrom(g.chrom), strand: g.strand, source: 'file', source_label: file.info.file, transcripts: models };
  }

  /** every file's genes overlapping [start, end) 0-based of a viewer chromosome, file by file (`all`: the track's entities too) */
  private async genesIn(chrom: string, start: number, end: number, all = false): Promise<AGene[][]> {
    return Promise.all(this.files.map(async f => (await Promise.all(this.fileChroms(f, chrom).map(c => f.overlapping(c, start, end)))).flat().filter(g => all || !g.trackOnly)));
  }

  /** genes overlapping [start, end] 1-based inclusive, with their shown model and the other exons merged */
  async regionGenes(chrom: string, start: number, end: number, exclude?: string): Promise<GeneModel[]> {
    // a gene of a later file with the name of an earlier one's is left out: the first file wins
    const out: GeneModel[] = [], seen = new Set<string>(), names = new Set<string>();
    (await this.genesIn(chrom, start - 1, end)).forEach((genes, fi) => {
      const mine = new Set<string>();
      for (const g of genes) {
        const nm = g.name.toUpperCase();
        if (seen.has(`${fi}|${g.key}`) || names.has(nm) || (exclude && (nm === exclude.toUpperCase() || g.id === exclude))) continue;
        seen.add(`${fi}|${g.key}`); mine.add(nm);
        const { tx: t, kind } = pickTx(g);
        out.push({
          gene_id: g.id, gene_name: g.name, biotype: g.biotype, strand: g.strand, start: g.start + 1, end: g.end,
          transcript_id: t.id, is_canonical: kind === 'canonical', exons: t.exons.map(([a, b]) => ({ start: a + 1, end: b })),
          cds_start: t.cdsStart != null ? t.cdsStart + 1 : null, cds_end: t.cdsEnd,
          other_exons: merge(g.txs.filter(x => x !== t).flatMap(x => x.exons)).map(([a, b]) => ({ start: a + 1, end: b })),
        });
      }
      for (const n of mine) names.add(n);
    });
    return out.sort((a, b) => a.start - b.start);
  }

  /** coding transcripts overlapping [start, end) 0-based (the SV browser's models) */
  async regionTranscripts(chrom: string, start: number, end: number): Promise<SvTranscript[]> {
    const out: SvTranscript[] = [];
    for (const genes of await this.genesIn(chrom, start, end)) for (const g of genes) {
      const picked = pickTx(g);
      for (const t of g.txs) {
        if (t.cdsStart == null) continue;
        const n = t.exons.length;
        out.push({ id: t.id, gene: g.name, strand: t.strand, start: t.exons[0][0], end: t.exons[n - 1][1], exons: t.exons.map(([a, b], i) => ({ start: a, end: b, rank: t.strand > 0 ? i + 1 : n - i })), cdsStart: t.cdsStart, cdsEnd: t.cdsEnd, status: t === picked.tx && picked.kind === 'canonical' ? 'select' : null });
      }
    }
    return out;
  }

  /** the sequences carried by the files (##FASTA sections), by name, when any: a reference when no FASTA is given */
  get embeddedSequences(): Map<string, string> | null {
    if (this.seqs === undefined) {
      const m = new Map<string, string>();
      for (const f of this.files) for (const [k, v] of f.sequences ?? []) if (!m.has(k)) m.set(k, v);
      this.seqs = m.size ? m : null;
    }
    return this.seqs;
  }
  private seqs: Map<string, string> | null | undefined;
  /** the file names, in order: one track each */
  get trackFiles(): string[] { return this.files.map(f => f.info.file); }
  /** every entity of each file overlapping [start, end) 0-based of a viewer chromosome, file by file (no name merged) */
  async trackItems(chrom: string, start: number, end: number): Promise<TrackData[]> {
    const per = await this.genesIn(chrom, start, end, true);
    return this.files.map((f, fi) => ({
      file: f.info.file, format: f.info.format,
      items: per[fi].map((g): TrackItem => {
        const picked = pickTx(g).tx;
        return {
          key: g.key, name: g.name, type: g.type ?? (f.info.format === 'bed' ? 'bed' : 'gene'), strand: g.unstranded ? 0 : g.strand,
          start: g.start, end: g.end, color: g.color, score: g.score, attrs: g.attrs, shown: Math.max(0, g.txs.indexOf(picked)),
          txs: g.txs.map(t => ({ id: t.id, name: t.name, biotype: t.biotype, strand: t.strand, exons: t.exons.map(e => [e[0], e[1]] as [number, number]), cdsStart: t.cdsStart, cdsEnd: t.cdsEnd })),
        };
      }).sort((a, b) => a.start - b.start),
    }));
  }

  /** names starting with the text (then containing it), for the gene box */
  suggest(text: string, n = 12): string[] {
    const q = text.trim().toUpperCase();
    if (!q) return [];
    const starts: string[] = [], contains: string[] = [], seen = new Set<string>();
    for (const f of this.files) for (const [k, name] of f.nameEntries()) {
      if (seen.has(name)) continue;
      const nm = name.toUpperCase();
      if (nm.startsWith(q)) { seen.add(name); starts.push(name); }
      else if (contains.length < 4 * n && (k.startsWith(q) || nm.includes(q))) { seen.add(name); contains.push(name); }
      if (starts.length >= 4 * n) break;
    }
    const cmp = (a: string, b: string) => a.length - b.length || a.localeCompare(b, undefined, { numeric: true });
    return [...starts.sort(cmp), ...contains.sort(cmp)].slice(0, n);
  }
}

/** An annotation file, and its tabix index when it has one. */
export interface AnnotationInput { file: File; index?: File }

/** Reads one file: whole into memory, or opened through its index (names then read in the background). */
async function loadOne(a: AnnotationInput, scan: NameScanner, onChange: () => void): Promise<GeneSource> {
  const { file, index } = a;
  const gz = await isGzipped(file, file.name);
  const format = await formatOf(file, file.name, gz);
  const info: AnnotationFileInfo = { file: file.name, format, genes: 0, transcripts: 0, lines: 0, skipped: 0, contigs: [], unmatched: [] };
  if (index) {
    if (!gz) throw new Error(`${index.name} indexes a bgzipped file: ${file.name} is not compressed`);
    info.indexed = { index: index.name, namesRead: 0, namesReady: false };
    const src = new TabixSource(info, file, index);
    await src.open();
    void src.readNames(file, scan, onChange);
    return src;
  }
  const sequences = new Map<string, string>();
  return new FileIndex(info, await parse(format, readLines(file, gz), info, format === 'gff3' ? sequences : undefined), sequences);
}

/**
 * Reads every file; a file that fails is reported, the others kept. `scan` builds the name index of an indexed file
 * (a worker's, by default here on the page); `onChange` is told as those indexes progress.
 */
export async function loadAnnotation(inputs: (AnnotationInput | File)[], opts: { scan?: NameScanner } = {}): Promise<{ annotation: Annotation; errors: string[] }> {
  const ok: GeneSource[] = [], errors: string[] = [];
  const annotation = new Annotation(ok);
  const changed = () => annotation.onChange?.();
  for (const x of inputs) {
    const a = x instanceof Blob ? { file: x } : x;
    try { ok.push(await loadOne(a, opts.scan ?? scanNames, changed)); } catch (e) { errors.push(`${a.file.name}: ${(e as Error).message || e}`); }
  }
  return { annotation, errors };
}
