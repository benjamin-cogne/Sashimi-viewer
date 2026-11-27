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
 */
import type { TranscriptData, AllTranscripts, TranscriptModel, GeneModel, RegionHint } from '../components/sashimi/types';
import type { SvTranscript } from '../components/sashimi/svGenes';
import { resolveContig } from '../components/sashimi/geometry';
export { resolveContig };

export type AnnotationFormat = 'gff3' | 'gtf' | 'bed';
export const ANNOTATION_EXT = /\.(gff3?|gtf|bed)(\.gz)?$/i;

interface ATx {
  id: string; name: string; chrom: string; strand: 1 | -1;
  /** 0-based half-open, sorted, disjoint */
  exons: [number, number][];
  cdsStart: number | null; cdsEnd: number | null;
  biotype: string; canonical: boolean;
}
interface AGene { key: string; id: string; name: string; chrom: string; strand: 1 | -1; start: number; end: number; biotype: string; txs: ATx[] }

/** What a file held, for its chip. */
export interface AnnotationFileInfo {
  file: string; format: AnnotationFormat; genes: number; transcripts: number; lines: number; skipped: number;
  /** contigs of the file, and those the reference does not hold (set once the reference's names are known) */
  contigs: string[]; unmatched: string[];
}

// ======================== reading ========================

const yieldToUi = (): Promise<void> => {
  const s = (globalThis as { scheduler?: { yield?: () => Promise<void> } }).scheduler;
  return s?.yield ? s.yield() : new Promise<void>(r => { setTimeout(r, 0); });
};

/** The file's lines, inflated as a stream when gzipped (a bgzipped file's blocks too), decoded as UTF-8. */
async function* readLines(file: Blob, gz: boolean): AsyncGenerator<string> {
  let stream: ReadableStream<Uint8Array> = file.stream();
  if (gz) stream = stream.pipeThrough(new DecompressionStream('gzip') as unknown as ReadableWritablePair<Uint8Array, Uint8Array>);
  const reader = stream.pipeThrough(new TextDecoderStream() as unknown as ReadableWritablePair<string, Uint8Array>).getReader();
  let rest = '';
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    const parts = (rest + value).split('\n');
    rest = parts.pop() ?? '';
    for (const p of parts) yield p.endsWith('\r') ? p.slice(0, -1) : p;
  }
  if (rest) yield rest;
}

async function isGzipped(file: Blob, name: string): Promise<boolean> {
  if (/\.gz$/i.test(name)) return true;
  const b = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return b[0] === 0x1f && b[1] === 0x8b;
}

/** The format from the name; a `.gff` is sniffed (key=value attributes: GFF3; key "value"; : GTF). */
async function formatOf(file: Blob, name: string, gz: boolean): Promise<AnnotationFormat> {
  const n = name.replace(/\.gz$/i, '').toLowerCase();
  if (n.endsWith('.gtf')) return 'gtf';
  if (n.endsWith('.bed')) return 'bed';
  if (n.endsWith('.gff3')) return 'gff3';
  let k = 0;
  for await (const line of readLines(file, gz)) {
    if (!line || line.startsWith('#')) continue;
    const col = line.split('\t')[8] ?? '';
    if (/^\s*\w+\s+"/.test(col)) return 'gtf';
    if (/=/.test(col)) return 'gff3';
    if (++k > 50) break;
  }
  return 'gff3';
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

interface Feat { id: string | null; parents: string[]; type: string; chrom: string; start: number; end: number; strand: 1 | -1; attrs: Map<string, string> }
const TX_TYPE = /^(mRNA|transcript|primary_transcript|pseudogenic_transcript|[a-zA-Z_]*RNA|[a-zA-Z_]*_transcript)$/;
const EXON_TYPE = /^(exon|pseudogenic_exon|noncoding_exon)$/;
const PART_TYPE = /^(CDS|five_prime_UTR|three_prime_UTR|UTR|5'UTR|3'UTR|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/;
const SKIP_LONE = /^(region|chromosome|contig|supercontig|scaffold|source|databank_entry|sequence_assembly|match|match_part|cDNA_match|biological_region)$/i;
const geneNameOf = (a: Map<string, string>, fallback: string) => a.get('Name') ?? a.get('gene_name') ?? a.get('gene') ?? a.get('locus_tag') ?? a.get('ID') ?? fallback;

/** One transcript from its exon and part children (0-based half-open). */
function txFrom(id: string, name: string, chrom: string, strand: 1 | -1, exonKids: Feat[], partKids: Feat[], span: [number, number] | null, biotype: string | undefined, canonical: boolean): ATx {
  const cds = partKids.filter(k => k.type === 'CDS' || (k.type === 'stop_codon'));
  const cdsOnly = partKids.filter(k => k.type === 'CDS');
  let exons = merge(exonKids.map(k => [k.start, k.end] as [number, number]));
  if (!exons.length) exons = merge(partKids.map(k => [k.start, k.end] as [number, number]));
  if (!exons.length && span) exons = [span];
  const cdsStart = cdsOnly.length ? Math.min(...cds.map(k => k.start)) : null, cdsEnd = cdsOnly.length ? Math.max(...cds.map(k => k.end)) : null;
  return { id, name, chrom, strand, exons, cdsStart, cdsEnd, biotype: biotype ?? (cdsOnly.length ? 'protein_coding' : 'non_coding'), canonical };
}

function geneOf(key: string, id: string, name: string, txs: ATx[], biotype?: string): AGene {
  const t0 = txs[0];
  return {
    key, id, name, chrom: t0.chrom, strand: t0.strand, txs,
    start: Math.min(...txs.map(t => t.exons[0][0])), end: Math.max(...txs.map(t => t.exons[t.exons.length - 1][1])),
    biotype: biotype ?? (txs.some(t => t.cdsStart != null) ? 'protein_coding' : t0.biotype),
  };
}

async function parseGff3(file: Blob, gz: boolean, info: AnnotationFileInfo): Promise<AGene[]> {
  const feats: Feat[] = [];
  for await (const line of readLines(file, gz)) {
    info.lines++;
    if (info.lines % 50000 === 0) await yieldToUi();
    if (line.startsWith('##FASTA')) break;
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (c.length < 9) { info.skipped++; continue; }
    const s = parseInt(c[3]), e = parseInt(c[4]);
    if (!(s > 0) || !(e >= s)) { info.skipped++; continue; }
    const attrs = gffAttrs(c[8]);
    feats.push({ id: attrs.get('ID') ?? null, parents: (attrs.get('Parent') ?? '').split(',').filter(Boolean), type: c[2], chrom: decode(c[0]), start: s - 1, end: e, strand: strandOf(c[6]), attrs });
  }
  const byId = new Map<string, Feat>();
  for (const f of feats) if (f.id && !byId.has(f.id)) byId.set(f.id, f);
  const kids = new Map<string, Feat[]>();
  for (const f of feats) for (const p of f.parents) { const l = kids.get(p); if (l) l.push(f); else kids.set(p, [f]); }
  const childrenOf = (f: Feat) => (f.id ? kids.get(f.id) ?? [] : []);
  const isTx = (f: Feat) => TX_TYPE.test(f.type) || childrenOf(f).some(k => EXON_TYPE.test(k.type) || k.type === 'CDS');
  // genes: features holding transcripts, or transcripts without a parent
  const genes = new Map<string, { gene: Feat | null; txs: ATx[] }>();
  const addTx = (geneKey: string, gene: Feat | null, tx: ATx) => { const g = genes.get(geneKey); if (g) g.txs.push(tx); else genes.set(geneKey, { gene, txs: [tx] }); };
  for (const f of feats) {
    if (!isTx(f) || EXON_TYPE.test(f.type) || PART_TYPE.test(f.type)) continue;
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
  // lone features: neither parent nor children (a region of interest, a gene with no transcript)
  for (const f of feats) {
    if (f.parents.length || childrenOf(f).length || SKIP_LONE.test(f.type) || EXON_TYPE.test(f.type) || PART_TYPE.test(f.type)) continue;
    if (f.id && genes.has(f.id)) continue;
    const name = geneNameOf(f.attrs, `${f.chrom}:${f.start + 1}-${f.end}`);
    addTx(f.id ?? `${name}@${f.chrom}:${f.start}`, f, { id: f.id ?? name, name, chrom: f.chrom, strand: f.strand, exons: [[f.start, f.end]], cdsStart: null, cdsEnd: null, biotype: biotypeAttr(f.attrs) ?? f.type, canonical: false });
  }
  const out: AGene[] = [];
  for (const [key, { gene, txs }] of genes) {
    const name = gene ? geneNameOf(gene.attrs, key) : txs[0].name;
    out.push(geneOf(key, gene?.id ?? key, name, txs, gene ? biotypeAttr(gene.attrs) : undefined));
  }
  return out;
}

async function parseGtf(file: Blob, gz: boolean, info: AnnotationFileInfo): Promise<AGene[]> {
  const txRows = new Map<string, { gene: string; chrom: string; strand: 1 | -1; exon: Feat[]; part: Feat[]; attrs: Map<string, string>; span: [number, number] | null }>();
  const geneAttrs = new Map<string, { attrs: Map<string, string>; f: Feat }>();
  for await (const line of readLines(file, gz)) {
    info.lines++;
    if (info.lines % 50000 === 0) await yieldToUi();
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (c.length < 9) { info.skipped++; continue; }
    const s = parseInt(c[3]), e = parseInt(c[4]);
    if (!(s > 0) || !(e >= s)) { info.skipped++; continue; }
    const attrs = gtfAttrs(c[8]);
    const f: Feat = { id: null, parents: [], type: c[2], chrom: c[0], start: s - 1, end: e, strand: strandOf(c[6]), attrs };
    const gid = attrs.get('gene_id') ?? attrs.get('gene_name') ?? '';
    const tid = attrs.get('transcript_id');
    if (c[2] === 'gene') { if (gid && !geneAttrs.has(gid)) geneAttrs.set(gid, { attrs, f }); continue; }
    if (!tid) { if (gid && !geneAttrs.has(gid)) geneAttrs.set(gid, { attrs, f }); continue; }
    let t = txRows.get(tid);
    if (!t) txRows.set(tid, t = { gene: gid || tid, chrom: c[0], strand: f.strand, exon: [], part: [], attrs, span: null });
    if (c[2] === 'transcript') { t.attrs = new Map([...t.attrs, ...attrs]); t.span = [f.start, f.end]; }
    else if (EXON_TYPE.test(c[2])) t.exon.push(f);
    else if (PART_TYPE.test(c[2]) || /^(UTR|five_prime_utr|three_prime_utr)$/i.test(c[2])) t.part.push(f);
  }
  const byGene = new Map<string, ATx[]>();
  for (const [tid, t] of txRows) {
    const tx = txFrom(tid, t.attrs.get('transcript_name') ?? tid, t.chrom, t.strand, t.exon, t.part, t.span, biotypeAttr(t.attrs), isCanonical(t.attrs));
    if (!tx.exons.length) continue;
    const l = byGene.get(t.gene); if (l) l.push(tx); else byGene.set(t.gene, [tx]);
  }
  const out: AGene[] = [];
  for (const [gid, txs] of byGene) {
    const ga = geneAttrs.get(gid)?.attrs, ta = txRows.get(txs[0].id)?.attrs;
    out.push(geneOf(gid, gid, ga?.get('gene_name') ?? ta?.get('gene_name') ?? gid, txs, ga ? biotypeAttr(ga) : undefined));
  }
  // genes without transcript rows: one exon over the gene
  for (const [gid, { attrs, f }] of geneAttrs) if (!byGene.has(gid)) {
    const name = attrs.get('gene_name') ?? gid;
    out.push(geneOf(gid, gid, name, [{ id: gid, name, chrom: f.chrom, strand: f.strand, exons: [[f.start, f.end]], cdsStart: null, cdsEnd: null, biotype: biotypeAttr(attrs) ?? 'gene', canonical: false }]));
  }
  return out;
}

async function parseBed(file: Blob, gz: boolean, info: AnnotationFileInfo): Promise<AGene[]> {
  const out: AGene[] = [];
  for await (const line of readLines(file, gz)) {
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
    out.push(geneOf(key, name, name, [{ id: name, name, chrom, strand, exons, cdsStart: coding ? ts : null, cdsEnd: coding ? te : null, biotype: coding ? 'protein_coding' : 'feature', canonical: false }]));
  }
  return out;
}

// ======================== index ========================

class FileIndex {
  readonly byChrom = new Map<string, { genes: AGene[]; maxEnd: number[] }>();
  readonly byName = new Map<string, AGene[]>();
  constructor(readonly info: AnnotationFileInfo, genes: AGene[]) {
    const per = new Map<string, AGene[]>();
    for (const g of genes) {
      const l = per.get(g.chrom); if (l) l.push(g); else per.set(g.chrom, [g]);
      const keys = new Set([g.name, g.id, ...g.txs.flatMap(t => [t.id, t.name])].map(k => k.toUpperCase()));
      for (const k of keys) { const l2 = this.byName.get(k); if (l2) l2.push(g); else this.byName.set(k, [g]); }
    }
    for (const [c, gs] of per) {
      gs.sort((a, b) => a.start - b.start);
      const maxEnd: number[] = []; let m = 0;
      for (const g of gs) { m = Math.max(m, g.end); maxEnd.push(m); }
      this.byChrom.set(c, { genes: gs, maxEnd });
    }
    info.genes = genes.length; info.transcripts = genes.reduce((n, g) => n + g.txs.length, 0); info.contigs = [...per.keys()];
  }
  /** genes of a contig (the file's spelling) overlapping [start, end) 0-based */
  overlapping(chrom: string, start: number, end: number): AGene[] {
    const c = this.byChrom.get(chrom);
    if (!c) return [];
    let lo = 0, hi = c.genes.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (c.genes[mid].start < end) lo = mid + 1; else hi = mid; }
    const out: AGene[] = [];
    for (let i = lo - 1; i >= 0 && c.maxEnd[i] > start; i--) if (c.genes[i].end > start) out.push(c.genes[i]);
    return out.reverse();
  }
}

/** Reads one annotation file into an index. */
export async function loadAnnotationFile(file: File): Promise<FileIndex> {
  const gz = await isGzipped(file, file.name);
  const format = await formatOf(file, file.name, gz);
  const info: AnnotationFileInfo = { file: file.name, format, genes: 0, transcripts: 0, lines: 0, skipped: 0, contigs: [], unmatched: [] };
  const genes = format === 'gtf' ? await parseGtf(file, gz, info) : format === 'bed' ? await parseBed(file, gz, info) : await parseGff3(file, gz, info);
  return new FileIndex(info, genes);
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
  readonly files: FileIndex[];
  private toRef = new Map<string, string>();
  constructor(files: FileIndex[]) { this.files = files; }

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
        const r = names?.length ? resolveContig(names, c) : spell ? spell(c) : c;
        if (r) this.toRef.set(c, r); else f.info.unmatched.push(c);
      }
    }
  }
  private refChrom(c: string): string { return this.toRef.get(c) ?? c; }
  /** the file contig names a viewer chromosome stands for */
  private fileChroms(f: FileIndex, chrom: string): string[] {
    const out: string[] = [];
    for (const c of f.byChrom.keys()) if (this.refChrom(c) === chrom || resolveContig([c], chrom)) out.push(c);
    return out;
  }

  /** a gene by gene name, gene id, transcript id or ID (case-insensitive); the copy nearest the hint when several */
  find(query: string, hint?: RegionHint): { gene: AGene; file: FileIndex; tx?: ATx } | null {
    const q = query.trim().toUpperCase();
    for (const f of this.files) {
      const list = f.byName.get(q);
      if (!list?.length) continue;
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

  transcript(query: string, hint?: RegionHint): TranscriptData | null {
    const hit = this.find(query, hint);
    if (!hit) return null;
    const { gene: g, file } = hit;
    const { tx: t, kind } = hit.tx ? { tx: hit.tx, kind: 'chosen' as const } : pickTx(g);
    const n = t.exons.length;
    return {
      gene_name: g.name, transcript_id: t.id, translation_id: null, is_mane_select: false, model_kind: kind, biotype: t.biotype,
      source: 'file', source_label: file.info.file, chrom: this.refChrom(t.chrom), strand: t.strand,
      start: t.exons[0][0] + 1, end: t.exons[n - 1][1],
      exons: t.exons.map(([a, b], i) => ({ start: a + 1, end: b, rank: t.strand > 0 ? i + 1 : n - i })),
      cds_start: t.cdsStart != null ? t.cdsStart + 1 : null, cds_end: t.cdsEnd,
    };
  }

  allTranscripts(query: string, hint?: RegionHint): AllTranscripts | null {
    const hit = this.find(query, hint);
    if (!hit) return null;
    const { gene: g, file } = hit, picked = pickTx(g).tx;
    const models = g.txs.map((t): TranscriptModel => ({
      id: t.id, name: t.name, source: 'file', biotype: t.biotype, start: t.exons[0][0] + 1, end: t.exons[t.exons.length - 1][1], strand: t.strand,
      exons: t.exons.map(([a, b]) => ({ start: a + 1, end: b })), cds_start: t.cdsStart != null ? t.cdsStart + 1 : null, cds_end: t.cdsEnd,
      is_mane: false, is_canonical: t === picked,
    })).sort((a, b) => Number(b.is_canonical) - Number(a.is_canonical) || a.id.localeCompare(b.id, undefined, { numeric: true }));
    return { gene_name: g.name, chrom: this.refChrom(g.chrom), strand: g.strand, source: 'file', source_label: file.info.file, transcripts: models };
  }

  /** genes overlapping [start, end] 1-based inclusive, with their shown model and the other exons merged */
  regionGenes(chrom: string, start: number, end: number, exclude?: string): GeneModel[] {
    // a gene of a later file with the name of an earlier one's is left out: the first file wins
    const out: GeneModel[] = [], seen = new Set<string>(), names = new Set<string>();
    this.files.forEach((f, fi) => {
      const mine = new Set<string>();
      for (const fc of this.fileChroms(f, chrom)) for (const g of f.overlapping(fc, start - 1, end)) {
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
  regionTranscripts(chrom: string, start: number, end: number): SvTranscript[] {
    const out: SvTranscript[] = [];
    for (const f of this.files) for (const fc of this.fileChroms(f, chrom)) for (const g of f.overlapping(fc, start, end)) {
      const picked = pickTx(g);
      for (const t of g.txs) {
        if (t.cdsStart == null) continue;
        const n = t.exons.length;
        out.push({ id: t.id, gene: g.name, strand: t.strand, start: t.exons[0][0], end: t.exons[n - 1][1], exons: t.exons.map(([a, b], i) => ({ start: a, end: b, rank: t.strand > 0 ? i + 1 : n - i })), cdsStart: t.cdsStart, cdsEnd: t.cdsEnd, status: t === picked.tx && picked.kind === 'canonical' ? 'select' : null });
      }
    }
    return out;
  }

  /** names starting with the text (then containing it), for the gene box */
  suggest(text: string, n = 12): string[] {
    const q = text.trim().toUpperCase();
    if (!q) return [];
    const starts: string[] = [], contains: string[] = [], seen = new Set<string>();
    for (const f of this.files) for (const [k, gs] of f.byName) {
      const name = gs[0].name;
      if (seen.has(name)) continue;
      const nm = name.toUpperCase();
      if (nm.startsWith(q)) { seen.add(name); starts.push(name); }
      else if (k.startsWith(q) || nm.includes(q)) { seen.add(name); contains.push(name); }
      if (starts.length >= n) break;
    }
    const cmp = (a: string, b: string) => a.length - b.length || a.localeCompare(b, undefined, { numeric: true });
    return [...starts.sort(cmp), ...contains.sort(cmp)].slice(0, n);
  }
}

/** Reads every file; a file that fails is reported, the others kept. */
export async function loadAnnotation(files: File[]): Promise<{ annotation: Annotation; errors: string[] }> {
  const ok: FileIndex[] = [], errors: string[] = [];
  for (const f of files) {
    try { ok.push(await loadAnnotationFile(f)); } catch (e) { errors.push(`${f.name}: ${(e as Error).message || e}`); }
  }
  return { annotation: new Annotation(ok), errors };
}
