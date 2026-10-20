/**
 * Splicing consequence model for the cartoon: the aberrant transcript implied by a junction is
 * assembled from the reference sequence, translated, and judged for nonsense-mediated decay.
 * Pure functions, shared by the application and the standalone viewer.
 *
 * NMD rules applied (see Nagy & Maquat 1998 TIBS; Lindeboom et al. 2016 Nat Genet; Kurosaki et
 * al. 2019 Nat Rev Mol Cell Biol): a premature termination codon (PTC) more than 50–55 nt upstream
 * of the last exon–exon junction triggers EJC-dependent NMD; a PTC in the last exon or within
 * ~50 nt of the last junction escapes; a PTC within ~150 nt of the start codon often escapes
 * through translation re-initiation; a long 3′ UTR (> ~1 kb) can trigger EJC-independent NMD.
 */
import { pseudoExonOf, translateCodon, type TxModel, type Exon0 } from './geometry';
import type { ProteinDomain } from './types';
import { fusionFrame, layoutOf, type FusionEnd, type FusionFrame } from './fusionFrame';

// ======================== Splice event → blocks ========================

export type BlockKind = 'exon' | 'cryptic' | 'extension';

export interface Block {
  /** genomic, 0-based half-open */
  start: number; end: number;
  kind: BlockKind;
  /** exon rank when the block is (part of) an annotated exon */
  rank?: number;
  /** true when the block is genomically contiguous with the previous one in transcription order (no splice junction between them) */
  joinsPrev?: boolean;
  /** the gene the block comes from, in a fusion transcript (exon ranks are the gene's own) */
  gene?: string;
  /** its coding part [start, end), when it is not the one of the model drawn (the 3′ partner's exons of a fusion) */
  cds?: [number, number] | null;
  /** its chromosome and strand, when they are not the model's (a fusion's 3′ partner on another chromosome or strand) */
  chrom?: string;
  strand?: 1 | -1;
}

export type EventKind = 'canonical' | 'skip' | 'cryptic_exon' | 'exonic_site' | 'intronic_site' | 'fusion' | 'unsupported';

export interface SpliceEvent {
  kind: EventKind;
  text: string;
  /** ranks of the exons removed (skip) or changed (cryptic sites) */
  affected: number[];
  /** aberrant blocks, genomic order */
  blocks: Block[];
  /** canonical blocks, genomic order */
  canonical: Block[];
}

const exonBlocks = (tx: TxModel): Block[] => tx.exons.map(e => ({ start: e.start, end: e.end, kind: 'exon' as const, rank: e.rank }));

/** Which exons and which parts of them the aberrant transcript keeps, from one junction (genomic reasoning, strand-agnostic). */
export function spliceEvent(j: { start: number; end: number }, tx: TxModel, others: { start: number; end: number }[] = []): SpliceEvent {
  const ex = tx.exons;
  const canonical = exonBlocks(tx);
  const li = ex.findIndex(e => e.end === j.start), ri = ex.findIndex(e => e.start === j.end);
  const plus = tx.strand > 0;
  const unsupported = (text: string): SpliceEvent => ({ kind: 'unsupported', text, affected: [], blocks: canonical, canonical });
  if (li >= 0 && ri >= 0) {
    if (ri === li + 1) return { kind: 'canonical', text: 'canonical junction: normal splicing', affected: [], blocks: canonical, canonical };
    if (ri > li + 1) {
      const removed = ex.slice(li + 1, ri).map(e => e.rank);
      const blocks = canonical.filter(b => !removed.includes(b.rank!));
      return { kind: 'skip', text: `exon${removed.length > 1 ? 's' : ''} ${removed.join(', ')} skipped`, affected: removed, blocks, canonical };
    }
    return unsupported('unexpected exon order (back-splice?)');
  }
  const pe = pseudoExonOf(j, tx, others);
  if (pe) {
    const blocks = [...canonical, { start: pe.start, end: pe.end, kind: 'cryptic' as const }].sort((a, b) => a.start - b.start);
    return { kind: 'cryptic_exon', text: `${pe.end - pe.start}-nt cryptic exon inserted between exons ${Math.min(pe.upExon.rank, pe.downExon.rank)} and ${Math.max(pe.upExon.rank, pe.downExon.rank)}`, affected: [pe.upExon.rank, pe.downExon.rank], blocks, canonical };
  }
  if (li >= 0 && ri < 0) {
    // novel site at j.end, downstream (genomically) of exon li
    const down = ex[li + 1];
    if (!down) return unsupported('novel site beyond the last exon');
    if (j.end > down.start && j.end < down.end) {
      const blocks = canonical.map(b => (b.rank === down.rank ? { ...b, start: j.end } : b));
      const side = plus ? 'acceptor' : 'donor';
      return { kind: 'exonic_site', text: `cryptic ${side} inside exon ${down.rank}: ${j.end - down.start} nt of the exon lost`, affected: [down.rank], blocks, canonical };
    }
    if (j.end < down.start) {
      const blocks = [...canonical, { start: j.end, end: down.start, kind: 'extension' as const, rank: down.rank }].sort((a, b) => a.start - b.start);
      const side = plus ? 'acceptor' : 'donor';
      return { kind: 'intronic_site', text: `cryptic ${side} in the intron: ${down.start - j.end} nt of intron added to exon ${down.rank}`, affected: [down.rank], blocks, canonical };
    }
    return unsupported('novel site past the next exon');
  }
  if (ri >= 0 && li < 0) {
    const up = ex[ri - 1];
    if (!up) return unsupported('novel site before the first exon');
    if (j.start > up.start && j.start < up.end) {
      const blocks = canonical.map(b => (b.rank === up.rank ? { ...b, end: j.start } : b));
      const side = plus ? 'donor' : 'acceptor';
      return { kind: 'exonic_site', text: `cryptic ${side} inside exon ${up.rank}: ${up.end - j.start} nt of the exon lost`, affected: [up.rank], blocks, canonical };
    }
    if (j.start > up.end) {
      const blocks = [...canonical, { start: up.end, end: j.start, kind: 'extension' as const, rank: up.rank }].sort((a, b) => a.start - b.start);
      const side = plus ? 'donor' : 'acceptor';
      return { kind: 'intronic_site', text: `cryptic ${side} in the intron: ${j.start - up.end} nt of intron added to exon ${up.rank}`, affected: [up.rank], blocks, canonical };
    }
    return unsupported('novel site before the previous exon');
  }
  return unsupported('neither end of the junction matches an annotated exon boundary');
}

// ======================== Isoform assembly and translation ========================

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
export const revComp = (s: string) => s.split('').reverse().map(c => COMP[c] ?? 'N').join('');

export interface Isoform {
  /** blocks in transcription order, `joinsPrev` set */
  blocks: Block[];
  /** mRNA sense sequence (T for U) */
  mrna: string;
  /** mRNA index of the first base of each block (transcription order) */
  blockStart: number[];
  /** mRNA indices where a splice junction sits (first base after the junction) */
  junctions: number[];
  /** mRNA index of the A of the start codon, null when the start codon is not in the transcript */
  cdsStart: number | null;
  /** translated protein (no stop symbol) */
  protein: string;
  /** mRNA index of the first base of the stop codon; null when translation runs off the transcript */
  stopIndex: number | null;
  /** genomic coordinate (0-based) of the first base of the stop codon, for comparison with the canonical stop */
  stopGenomic: number | null;
}

/**
 * `seqAt(start, end, chrom)` returns the + strand reference of a genomic window (`chrom` is given for the blocks that
 * carry their own). The blocks are put in transcription order, unless `inOrder` says they already are (a fusion
 * transcript, whose 3′ partner may lie upstream of the 5′ one, on the other strand or on another chromosome).
 */
export function buildIsoform(blocksGenomic: Block[], tx: Pick<TxModel, 'strand' | 'cdsStart' | 'cdsEnd'> & { chrom?: string }, seqAt: (start: number, end: number, chrom?: string) => string, inOrder = false): Isoform {
  const plus = tx.strand > 0;
  const fwd = (b: Block) => (b.strand ?? tx.strand) > 0;
  const ordered = inOrder ? blocksGenomic : plus ? [...blocksGenomic].sort((a, b) => a.start - b.start) : [...blocksGenomic].sort((a, b) => b.start - a.start);
  const blocks: Block[] = ordered.map((b, i) => {
    const prev = ordered[i - 1];
    const joins = prev && prev.chrom === b.chrom && fwd(prev) === fwd(b) ? (fwd(b) ? prev.end === b.start : prev.start === b.end) : false;
    return { ...b, joinsPrev: joins };
  });
  const parts = blocks.map(b => { const s = seqAt(b.start, b.end, b.chrom).toUpperCase(); return fwd(b) ? s : revComp(s); });
  const blockStart: number[] = [];
  let acc = 0;
  for (const p of parts) { blockStart.push(acc); acc += p.length; }
  const mrna = parts.join('');
  const junctions = blocks.map((b, i) => (i > 0 && !b.joinsPrev ? blockStart[i] : -1)).filter(x => x >= 0);
  // transcript index of a genomic base (on `tx.chrom` when the blocks carry chromosomes)
  const toMrna = (g: number): number | null => {
    for (let i = 0; i < blocks.length; i++) {
      const b = blocks[i];
      if (b.chrom != null && tx.chrom != null && b.chrom !== tx.chrom) continue;
      if (g >= b.start && g < b.end) return blockStart[i] + (fwd(b) ? g - b.start : b.end - 1 - g);
    }
    return null;
  };
  const toGenomic = (m: number): number | null => {
    for (let i = 0; i < blocks.length; i++) {
      const len = blocks[i].end - blocks[i].start;
      if (m >= blockStart[i] && m < blockStart[i] + len) return fwd(blocks[i]) ? blocks[i].start + (m - blockStart[i]) : blocks[i].end - 1 - (m - blockStart[i]);
    }
    return null;
  };
  const startG = tx.cdsStart != null && tx.cdsEnd != null ? (plus ? tx.cdsStart : tx.cdsEnd - 1) : null;
  const cdsStart = startG != null ? toMrna(startG) : null;
  let protein = '', stopIndex: number | null = null;
  if (cdsStart != null) {
    for (let i = cdsStart; i + 3 <= mrna.length; i += 3) {
      const aa = translateCodon(mrna.slice(i, i + 3));
      if (aa === '*') { stopIndex = i; break; }
      protein += aa;
    }
  }
  return { blocks, mrna, blockStart, junctions, cdsStart, protein, stopIndex, stopGenomic: stopIndex != null ? toGenomic(stopIndex) : null };
}

// ======================== NMD verdict ========================

export type Verdict = 'nmd' | 'escape_last_exon' | 'escape_near_junction' | 'escape_start_proximal' | 'no_ptc' | 'no_stop' | 'start_lost' | 'unknown';

export interface NmdVerdict {
  verdict: Verdict;
  /** the aberrant transcript is expected to be degraded */
  degraded: boolean;
  /** mRNA index of the PTC (first base) */
  ptc: number | null;
  /** codon number of the PTC in the aberrant protein (1-based) */
  ptcCodon: number | null;
  /** nt from the PTC to the last exon–exon junction (positive = upstream of it) */
  distanceToLastJunction: number | null;
  /** nt from the start codon to the PTC */
  distanceFromStart: number | null;
  utr3Length: number | null;
  longUtr: boolean;
  headline: string;
  text: string;
}

export const NMD_JUNCTION_NT = 55;
export const NMD_START_PROXIMAL_NT = 150;
export const NMD_LONG_UTR_NT = 1000;

export function nmdVerdict(iso: Isoform, canonical: Isoform): NmdVerdict {
  const base = { ptc: null, ptcCodon: null, distanceToLastJunction: null, distanceFromStart: null, utr3Length: null, longUtr: false };
  if (iso.cdsStart == null) return { ...base, verdict: 'start_lost', degraded: false, headline: 'Start codon lost', text: 'The annotated AUG is not part of the transcript: no protein from this start; translation could re-initiate at a downstream AUG.' };
  if (iso.stopIndex == null) return { ...base, verdict: 'no_stop', degraded: true, headline: 'No stop codon', text: 'Translation runs into the poly(A) tail: substrate of non-stop decay.' };
  const utr3 = iso.mrna.length - (iso.stopIndex + 3);
  const longUtr = utr3 > NMD_LONG_UTR_NT;
  const fromStart = iso.stopIndex - iso.cdsStart;
  if (canonical.stopGenomic != null && iso.stopGenomic === canonical.stopGenomic) {
    return { ...base, verdict: 'no_ptc', degraded: false, ptc: null, utr3Length: utr3, longUtr: false, distanceFromStart: fromStart,
      headline: 'Canonical stop codon reached', text: 'No premature termination codon: the transcript is translated to its normal stop.' };
  }
  const ptcCodon = Math.floor(fromStart / 3) + 1;
  const lastJ = iso.junctions.length ? iso.junctions[iso.junctions.length - 1] : null;
  const d = lastJ != null ? lastJ - iso.stopIndex : null;
  const common = { ptc: iso.stopIndex, ptcCodon, distanceToLastJunction: d, distanceFromStart: fromStart, utr3Length: utr3, longUtr };
  const utrNote = longUtr ? ` The ${utr3.toLocaleString('en-US')}-nt 3′ UTR left downstream is long, which can trigger EJC-independent NMD.` : '';
  if (lastJ == null || iso.stopIndex >= lastJ) {
    return { ...common, verdict: 'escape_last_exon', degraded: false, headline: 'PTC in the last exon: escapes NMD',
      text: `The premature stop (codon ${ptcCodon}) lies in the last exon, downstream of every exon junction complex: the truncated protein is made.${utrNote}` };
  }
  if (d! <= NMD_JUNCTION_NT) {
    return { ...common, verdict: 'escape_near_junction', degraded: false, headline: `PTC ${d} nt before the last junction: escapes NMD`,
      text: `A stop within ${NMD_JUNCTION_NT} nt of the last exon–exon junction is not recognised as premature (the ribosome displaces the last EJC): the truncated protein is made.${utrNote}` };
  }
  if (fromStart < NMD_START_PROXIMAL_NT) {
    return { ...common, verdict: 'escape_start_proximal', degraded: false, headline: `PTC ${fromStart} nt after the AUG: probable escape by re-initiation`,
      text: `Start-proximal stops (< ${NMD_START_PROXIMAL_NT} nt) often escape NMD because translation re-initiates at a downstream AUG; the outcome is an N-truncated protein rather than decay. Treat as uncertain.` };
  }
  return { ...common, verdict: 'nmd', degraded: true, headline: `PTC ${d} nt upstream of the last junction: degraded by NMD`,
    text: `The stop at codon ${ptcCodon} leaves ${iso.junctions.filter(x => x > iso.stopIndex!).length} exon junction complex${iso.junctions.filter(x => x > iso.stopIndex!).length > 1 ? 'es' : ''} downstream; UPF1 is recruited and the transcript is degraded (50–55-nt rule).` };
}

// ======================== Protein comparison and domains ========================

export type DiffKind = 'identical' | 'in_frame_deletion' | 'in_frame_insertion' | 'in_frame_indel' | 'frameshift' | 'truncation' | 'no_protein'
  /** a fusion protein: the 5′ partner's residues up to `prefix`, then the 3′ partner's, in frame to its own stop codon */
  | 'fusion';

export interface ProteinDiff {
  kind: DiffKind;
  /** identical residues from the N-terminus */
  prefix: number;
  /** identical residues from the C-terminus (in-frame changes only) */
  suffix: number;
  /** canonical residues lost, 0-based half-open */
  lost: [number, number];
  /** residues inserted (in frame) or the frameshifted neopeptide */
  inserted: string;
  canonicalLength: number;
  aberrantLength: number;
  text: string;
}

export function proteinDiff(canon: string, ab: string, canonicalStopReached: boolean): ProteinDiff {
  if (!ab.length) return { kind: 'no_protein', prefix: 0, suffix: 0, lost: [0, canon.length], inserted: '', canonicalLength: canon.length, aberrantLength: 0, text: 'no protein' };
  let prefix = 0;
  while (prefix < canon.length && prefix < ab.length && canon[prefix] === ab[prefix]) prefix++;
  if (canon === ab) return { kind: 'identical', prefix, suffix: 0, lost: [prefix, prefix], inserted: '', canonicalLength: canon.length, aberrantLength: ab.length, text: 'identical protein' };
  if (canonicalStopReached) {
    let suffix = 0;
    while (suffix < canon.length - prefix && suffix < ab.length - prefix && canon[canon.length - 1 - suffix] === ab[ab.length - 1 - suffix]) suffix++;
    const lost: [number, number] = [prefix, canon.length - suffix];
    const inserted = ab.slice(prefix, ab.length - suffix);
    const nLost = lost[1] - lost[0];
    const kind: DiffKind = nLost && inserted ? 'in_frame_indel' : nLost ? 'in_frame_deletion' : 'in_frame_insertion';
    const text = kind === 'in_frame_deletion' ? `in-frame deletion of ${nLost} aa (p.${canon[lost[0]]}${lost[0] + 1}_${canon[lost[1] - 1]}${lost[1]}del)`
      : kind === 'in_frame_insertion' ? `in-frame insertion of ${inserted.length} aa after residue ${prefix}`
      : `${nLost} aa replaced by ${inserted.length} aa (residues ${lost[0] + 1}–${lost[1]})`;
    return { kind, prefix, suffix, lost, inserted, canonicalLength: canon.length, aberrantLength: ab.length, text };
  }
  const neo = ab.slice(prefix);
  const kind: DiffKind = neo.length ? 'frameshift' : 'truncation';
  const text = kind === 'frameshift'
    ? `frameshift after residue ${prefix} (p.${canon[prefix] ?? '?'}${prefix + 1}${ab[prefix]}fs*${neo.length + 1}): ${neo.length} aberrant aa then a premature stop`
    : `truncation after residue ${prefix} (p.${canon[prefix] ?? '?'}${prefix + 1}*)`;
  return { kind, prefix, suffix: 0, lost: [prefix, canon.length], inserted: neo, canonicalLength: canon.length, aberrantLength: ab.length, text };
}

export type DomainState = 'intact' | 'lost' | 'disrupted';
export interface DomainStatus {
  domain: ProteinDomain; state: DomainState; note: string;
  /** the gene whose protein carries the domain, in a fusion */
  gene?: string;
  /** residue ranges drawn on each row (0-based half-open), when they do not follow from the diff (a fusion) */
  rows?: { canonical: [number, number][]; observed: [number, number][] };
}

/** UniProt domains first, then Pfam, SMART, PANTHER / Superfamily; overlapping hits of the same id merged. */
export function pickDomains(features: ProteinDomain[]): ProteinDomain[] {
  const order = ['UniProt', 'Pfam', 'Smart', 'SMART', 'PANTHER', 'Superfamily', 'Prosite_profiles', 'Gene3D', 'CDD'];
  for (const t of order) {
    const hits = features.filter(f => f.type === t).sort((a, b) => a.start - b.start);
    if (!hits.length) continue;
    const out: ProteinDomain[] = [];
    for (const h of hits) {
      const last = out[out.length - 1];
      if (last && last.id === h.id && h.start <= last.end + 5) last.end = Math.max(last.end, h.end);
      else out.push({ ...h });
    }
    return out;
  }
  return [];
}

export function domainStatus(domains: ProteinDomain[], diff: ProteinDiff): DomainStatus[] {
  return domains.map(dom => {
    const s = dom.start - 1, e = dom.end;   // 0-based half-open
    const [a, b] = diff.lost;
    if (diff.kind === 'identical') return { domain: dom, state: 'intact', note: 'unchanged' };
    if (diff.kind === 'no_protein') return { domain: dom, state: 'lost', note: 'no protein' };
    if (diff.kind === 'frameshift' || diff.kind === 'truncation') {
      if (e <= a) return { domain: dom, state: 'intact', note: 'upstream of the change' };
      if (s >= a) return { domain: dom, state: 'lost', note: 'downstream of the premature stop' };
      return { domain: dom, state: 'disrupted', note: `truncated after residue ${a} (${a - s} of ${e - s} aa kept)` };
    }
    // in-frame changes
    if (diff.kind === 'in_frame_insertion') {
      return s < a && e > a ? { domain: dom, state: 'disrupted', note: `${diff.inserted.length} aa inserted inside` } : { domain: dom, state: 'intact', note: e <= a ? 'unchanged' : 'shifted, sequence unchanged' };
    }
    if (e <= a || s >= b) return { domain: dom, state: 'intact', note: s >= b ? 'shifted, sequence unchanged' : 'unchanged' };
    if (s >= a && e <= b) return { domain: dom, state: 'lost', note: 'entirely deleted' };
    const kept = (e - s) - (Math.min(e, b) - Math.max(s, a));
    return { domain: dom, state: 'disrupted', note: `${e - s - kept} of ${e - s} aa deleted${diff.inserted ? `, ${diff.inserted.length} aa inserted` : ''}` };
  });
}

// ======================== Everything the cartoon needs ========================

export interface SpliceStory {
  event: SpliceEvent;
  canonical: Isoform;
  aberrant: Isoform;
  verdict: NmdVerdict;
  diff: ProteinDiff;
  domains: DomainStatus[];
  /** a fusion transcript: the canonical row is the 5′ partner's transcript, the observed row the fused one */
  fusion?: FusionStoryInfo;
}

export interface FusionStoryInfo {
  five: TxModel; three: TxModel;
  frame: FusionFrame;
  /** fused-protein residues [from, to) that come from the 3′ partner in its own frame (null: none, out of frame) */
  partnerResidues: [number, number] | null;
  /** the 3′ partner's protein length, and its first residue carried (1-based) */
  partnerLength: number;
  partnerFrom: number | null;
  /** the start codon used: the 5′ partner's, or the 3′ partner's when the junction is in the 5′ partner's 5′ UTR */
  startFrom: 'five' | 'three' | null;
}

export function spliceStory(j: { start: number; end: number }, tx: TxModel, others: { start: number; end: number }[], seqAt: (s: number, e: number) => string, features: ProteinDomain[]): SpliceStory {
  const event = spliceEvent(j, tx, others);
  const canonical = buildIsoform(event.canonical, tx, seqAt);
  const aberrant = buildIsoform(event.blocks, tx, seqAt);
  const verdict = nmdVerdict(aberrant, canonical);
  const diff = proteinDiff(canonical.protein, aberrant.protein, verdict.verdict === 'no_ptc');
  if (verdict.verdict === 'no_stop' && diff.kind === 'frameshift') diff.text = `frameshift after residue ${diff.prefix}: no stop codon before the poly(A) tail`;
  return { event, canonical, aberrant, verdict, diff, domains: domainStatus(pickDomains(features), diff) };
}

/** Genomic windows whose sequence the story needs (blocks of both isoforms), merged when close. */
export function storyWindows(event: SpliceEvent, mergeGap = 5000): [number, number][] {
  return blockWindows([...event.canonical, ...event.blocks], mergeGap);
}

function blockWindows(blocks: { start: number; end: number }[], mergeGap: number): [number, number][] {
  const all = blocks.map(b => [b.start, b.end] as [number, number]).sort((a, b) => a[0] - b[0]);
  const out: [number, number][] = [];
  for (const w of all) {
    const last = out[out.length - 1];
    if (last && w[0] <= last[1] + mergeGap) last[1] = Math.max(last[1], w[1]); else out.push([w[0], w[1]]);
  }
  return out;
}

// ======================== Fusion transcripts ========================

const cdsPart = (start: number, end: number, m: TxModel): [number, number] | null => {
  if (m.cdsStart == null || m.cdsEnd == null) return null;
  const a = Math.max(start, m.cdsStart), b = Math.min(end, m.cdsEnd);
  return b > a ? [a, b] : null;
};

/** The exons of a model in transcription order, tagged with its gene and their coding part. */
function txBlocks(m: TxModel): Block[] {
  const plus = m.strand > 0;
  return [...m.exons].sort((a, b) => (plus ? a.start - b.start : b.start - a.start))
    .map(e => ({ start: e.start, end: e.end, kind: 'exon' as const, rank: e.rank, gene: m.geneName, cds: cdsPart(e.start, e.end, m), chrom: m.chrom, strand: m.strand }));
}

/**
 * The fused transcript in transcription order: the 5′ partner's exons up to its donor, then the 3′ partner's from its
 * acceptor (the two partners on one strand; the ends are exon boundaries, as `fusionFrame` takes them).
 */
export function fusionBlocks(five: FusionEnd, three: FusionEnd): Block[] {
  // each partner on its own strand (and chromosome): an inversion or a translocation joins them all the same
  const pa = five.model.strand > 0, pb = three.model.strand > 0;
  const clip = (b: Block, start: number, end: number, m: TxModel): Block => ({ ...b, start, end, cds: cdsPart(start, end, m) });
  const a = txBlocks(five.model).filter(b => (pa ? b.start < five.pos : b.end > five.pos))
    .map(b => (pa ? clip(b, b.start, Math.min(b.end, five.pos), five.model) : clip(b, Math.max(b.start, five.pos), b.end, five.model)));
  const c = txBlocks(three.model).filter(b => (pb ? b.end > three.pos : b.start < three.pos))
    .map(b => (pb ? clip(b, Math.max(b.start, three.pos), b.end, three.model) : clip(b, b.start, Math.min(b.end, three.pos), three.model)));
  return [...a, ...c];
}

/** Genomic windows whose sequence a fusion story needs: every exon of the two partners, on each one's chromosome. */
export function fusionWindows(five: FusionEnd, three: FusionEnd, mergeGap = 5000): { chrom: string; start: number; end: number }[] {
  const A = five.model, B = three.model;
  if (A.chrom !== B.chrom) return [...blockWindows(A.exons, mergeGap).map(([start, end]) => ({ chrom: A.chrom, start, end })), ...blockWindows(B.exons, mergeGap).map(([start, end]) => ({ chrom: B.chrom, start, end }))];
  return blockWindows([...A.exons, ...B.exons], mergeGap).map(([start, end]) => ({ chrom: A.chrom, start, end }));
}

/** A fusion story can be told when both ends are on their partner's exons and a start codon is kept. */
export function fusionStoryPossible(five: FusionEnd, three: FusionEnd): boolean {
  if (five.model.geneName === three.model.geneName) return false;
  const f = fusionFrame(five, three);
  return f.aBases != null || f.verdict === 'five-utr';
}

const mod3 = (n: number) => ((n % 3) + 3) % 3;

/**
 * The fused transcript of `five` (donor end) joined to `three` (acceptor end), translated and judged like an aberrant
 * splice product: the canonical row is the 5′ partner's transcript, the observed row the fusion. Translation starts at
 * the 5′ partner's start codon, or at the 3′ partner's when the junction lies in the 5′ partner's 5′ UTR (a promoter
 * swap). The stop codon is premature unless it is the 3′ partner's own, so the NMD rules apply as to any transcript:
 * an in-frame fusion read to the 3′ partner's stop escapes; a stop more than 55 nt before the fused transcript's last
 * exon junction is degraded. The 5′ partner's domains are kept up to the junction; the 3′ partner's are carried when
 * they lie after it and the fusion is read in their frame.
 */
export function fusionStory(five: FusionEnd, three: FusionEnd, seqAt: (s: number, e: number, chrom?: string) => string, featuresFive: ProteinDomain[], featuresThree: ProteinDomain[]): SpliceStory {
  const A = five.model, B = three.model, gA = A.geneName, gB = B.geneName;
  const plusOf = (m: TxModel) => m.strand > 0;
  const frame = fusionFrame(five, three);
  const canonBlocks = txBlocks(A), blocks = fusionBlocks(five, three);
  const nA = blocks.filter(b => b.gene === gA).length;
  const canonical = buildIsoform(canonBlocks, A, seqAt, true);
  const partner = buildIsoform(txBlocks(B), B, seqAt, true);
  const startOf = (m: TxModel) => (m.cdsStart != null && m.cdsEnd != null ? (plusOf(m) ? m.cdsStart : m.cdsEnd - 1) : null);
  const holds = (g: number | null, from: number, to: number) => g != null && blocks.slice(from, to).some(b => g >= b.start && g < b.end);
  const sA = startOf(A), sB = startOf(B);
  const startFrom = holds(sA, 0, nA) ? 'five' as const : holds(sB, nA, blocks.length) ? 'three' as const : null;
  const startM = startFrom === 'five' ? A : B, startG = startFrom === 'five' ? sA : startFrom === 'three' ? sB : null;
  const aberrant = buildIsoform(blocks, { strand: startM.strand, chrom: startM.chrom, cdsStart: startG, cdsEnd: startG != null ? startG + 1 : null }, seqAt, true);
  // the stop codon is premature unless it is the 3′ partner's own (its first base, as Isoform.stopGenomic), reached in
  // the 3′ partner's part of the transcript (the 5′ one may lie on another chromosome, at the same coordinates)
  const ownStop = B.cdsStart != null && B.cdsEnd != null ? (plusOf(B) ? B.cdsEnd - 3 : B.cdsStart + 2) : null;
  const stopInB = aberrant.stopIndex != null && nA < blocks.length && aberrant.stopIndex >= aberrant.blockStart[nA];
  const v = nmdVerdict({ ...aberrant, stopGenomic: stopInB ? aberrant.stopGenomic : null }, { stopGenomic: ownStop } as Isoform);
  const pa = canonical.protein, pf = aberrant.protein;
  let prefix = 0;
  if (startFrom === 'five') while (prefix < pa.length && prefix < pf.length && pa[prefix] === pf[prefix]) prefix++;

  // the 3′ partner's residues the fusion carries in their own frame: its coding base `bOff` sits at the junction
  let partnerResidues: [number, number] | null = null, partnerFrom: number | null = null, shift: number | null = null;
  const LB = layoutOf(B), tb = LB.tOf(plusOf(B) ? three.pos : three.pos - 1);
  if (aberrant.cdsStart != null && tb != null && LB.cdsStart != null && nA < blocks.length) {
    const bOff = tb - LB.cdsStart;
    const d = aberrant.blockStart[nA] - bOff - aberrant.cdsStart;   // fused coding base of the 3′ partner's first one
    if (mod3(d) === 0) {
      shift = d / 3;
      const first = Math.max(0, Math.ceil(bOff / 3));   // its first codon wholly after the junction (0-based)
      const from = first + shift, to = Math.min(pf.length, partner.protein.length + shift);
      if (to > from && pf.slice(from, to) === partner.protein.slice(first, to - shift)) { partnerResidues = [from, to]; partnerFrom = first + 1; }
    }
  }
  const carried = v.verdict === 'no_ptc' && partnerResidues != null && partnerResidues[1] === partner.protein.length + shift!;

  let diff: ProteinDiff;
  const base = { prefix, suffix: 0, lost: [prefix, pa.length] as [number, number], canonicalLength: pa.length, aberrantLength: pf.length };
  if (!pf.length) diff = { ...base, kind: 'no_protein', inserted: '', text: `no protein: no start codon of ${gA} or ${gB} in the fused transcript` };
  else if (carried) {
    diff = { ...base, kind: 'fusion', inserted: pf.slice(prefix),
      text: startFrom === 'three'
        ? `${gB}'s whole protein (${partner.protein.length} aa) from its own start codon, under ${gA}'s promoter and 5′ UTR`
        : `${gA} p.1–${prefix} fused in frame to ${gB} p.${partnerFrom}–${partner.protein.length}: a ${pf.length}-aa fusion protein` };
  } else if (startFrom === 'five' && pf === pa) {
    diff = { ...base, kind: 'identical', inserted: '', text: `${gA}'s whole protein: the junction lies after its stop codon, ${gB}'s exons are untranslated` };
  } else {
    const neo = pf.slice(prefix);
    diff = { ...base, kind: neo.length ? 'frameshift' : 'truncation', inserted: neo,
      text: frame.verdict === 'in-frame'
        ? `in frame by the lengths, but a stop codon ${neo.length + 1} codons after ${gA} p.${prefix} ends the protein before ${gB}'s own sequence`
        : neo.length ? `${gA} p.1–${prefix}, then ${neo.length} aa of ${gB} read out of its frame, then a stop` : `${gA} p.1–${prefix}, then a stop codon at the junction` };
  }

  // the verdict, in the words of a fusion
  const verdict: NmdVerdict = carried
    ? { ...v, headline: `In-frame fusion: ${gB}'s own stop codon reached`,
        text: `Translated from ${startFrom === 'three' ? `${gB}'s` : `${gA}'s`} start codon through the junction to ${gB}'s own stop codon: no premature stop, the fused transcript escapes NMD and a ${pf.length}-aa ${startFrom === 'three' ? `${gB} protein` : 'fusion protein'} is expected (its stability is that of the protein, which this does not predict).` }
    : { ...v, text: `${frame.verdict === 'in-frame' ? 'In frame by the lengths, but the stop codon reached is not ' + gB + '\'s own. ' : frame.verdict === 'frameshift' ? `${gB} is read out of its frame after the junction. ` : frame.verdict === 'after-stop' ? `${gA}'s own stop codon lies before the junction; ${gB}'s exons become its 3′ UTR. ` : ''}${v.text}` };

  // domains: the 5′ partner's up to the junction, the 3′ partner's where they are carried in frame
  const domsA = domainStatus(pickDomains(featuresFive), diff.kind === 'identical' ? diff : { ...diff, kind: 'truncation' })
    .map(ds => {
      const s = ds.domain.start - 1, e = ds.domain.end;
      return { ...ds, gene: gA, note: ds.state === 'lost' ? `not in the fusion (after ${gA} p.${prefix})` : ds.note, rows: { canonical: [[s, e]] as [number, number][], observed: (e <= prefix ? [[s, e]] : s >= prefix ? [] : [[s, prefix]]) as [number, number][] } };
    });
  const domsB = pickDomains(featuresThree).map((dom): DomainStatus => {
    const s = dom.start - 1, e = dom.end;
    if (!partnerResidues || shift == null) return { domain: dom, state: 'lost', note: frame.verdict === 'frameshift' ? `${gB} read out of its frame` : `${gB}'s protein not carried`, gene: gB, rows: { canonical: [], observed: [] } };
    const a = Math.max(s + shift, partnerResidues[0]), b = Math.min(e + shift, partnerResidues[1]);
    const state: DomainState = b <= a ? 'lost' : b - a === e - s ? 'intact' : 'disrupted';
    const note = state === 'intact' ? `carried in frame (fusion residues ${a + 1}–${b})` : state === 'lost' ? (e + shift <= partnerResidues[0] ? `before the junction in ${gB}` : 'after the premature stop') : `${b - a} of ${e - s} aa carried`;
    return { domain: dom, state, note, gene: gB, rows: { canonical: [], observed: b > a ? [[a, b]] : [] } };
  });

  const rankA = blocks[nA - 1]?.rank, rankB = blocks[nA]?.rank;
  const event: SpliceEvent = { kind: 'fusion', text: `${gA}${rankA != null ? ` exon ${rankA}` : ''} → ${gB}${rankB != null ? ` exon ${rankB}` : ''} fusion transcript`, affected: [], blocks, canonical: canonBlocks };
  return { event, canonical, aberrant, verdict, diff, domains: [...domsA, ...domsB],
    fusion: { five: A, three: B, frame, partnerResidues, partnerLength: partner.protein.length, partnerFrom, startFrom } };
}

export { type Exon0 };
