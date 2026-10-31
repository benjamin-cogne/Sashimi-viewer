/**
 * What a structural variant does to the genes: every RefSeq NM_ transcript of the window is laid onto the rebuilt
 * chromosome (svModel.ts) and read back in its own direction of transcription.
 *
 * Each piece of the rebuilt chromosome holds a fragment of every transcript overlapping its segment. Fragments that
 * meet at a junction (or at a reference adjacency) on the same strand of the rebuilt chromosome are chained: a chain is
 * what RNA polymerase would read through. A chain is then compared with the reference transcript:
 *   - one transcript, every exon once, from its start to its end: an intact copy (perhaps inverted or moved);
 *   - one transcript with exons missing or repeated: an intragenic deletion or tandem duplication, in frame or not;
 *   - one transcript without its start or its end: the gene is cut, its parts separated;
 *   - two transcripts: a fusion, 5′ partner :: 3′ partner, in frame or not, or a promoter swap.
 * Gene by gene, the copies that stay intact give the dosage; ClinGen and gnomAD tell whether that dosage matters.
 * Chromatin domains (a BED of TADs the user loads) are laid onto the rebuilt chromosome the same way: two domain parts
 * joined at a junction make a new domain, and an intact gene whose domain gains foreign sequence or loses part of its
 * own may be misregulated (position effect).
 */
import type { SvReconstruction, SvJunction, SvPiece, SvSegment } from './svModel';

// ======================== Inputs ========================

/** One NM_ transcript, 0-based half-open, exons sorted by start with their rank in transcription order. */
export interface SvTranscript {
  id: string; gene: string; strand: 1 | -1; start: number; end: number;
  exons: { start: number; end: number; rank: number }[];
  cdsStart: number | null; cdsEnd: number | null;
  /** MANE Select / MANE Plus Clinical (GRCh38), RefSeq Select (the MANE set on GRCh38 for coding genes; the reference on GRCh37) */
  status: 'mane' | 'mane-plus' | 'select' | null;
}
export interface SvGene { symbol: string; strand: 1 | -1; start: number; end: number; transcripts: SvTranscript[]; principal: SvTranscript }

/** ClinGen dosage sensitivity: 0 no evidence, 1 little, 2 some, 3 sufficient, 30 autosomal recessive, 40 unlikely. */
export interface DosageScore { hi: number | null; ts: number | null }
export interface Constraint { pli: number | null; loeuf: number | null; version: string }
export interface PanelGene { symbol: string; confidence: 1 | 2 | 3; moi?: string; phenotypes?: string[] }
export interface GenePanel { source: string; id: string; name: string; version?: string; genes: Map<string, PanelGene> }
export interface Domain { start: number; end: number; name?: string }

export interface GeneAnnotations {
  dosage: Map<string, DosageScore>;
  constraint: Map<string, Constraint>;
  panel: GenePanel | null;
}

// ======================== Thresholds ========================
// Every one of them is printed in the panel's method note.

/**
 * gnomAD LOEUF below which a gene is called LoF-intolerant, when a hard threshold is needed. gnomAD v4: < 0.45, the
 * browser's help "for the interpretation of Mendelian disease cases" (the 85th LOEUF percentile); v2.1.1: < 0.35
 * (Karczewski et al. 2020). gnomAD advises reading LOEUF as a continuous value.
 */
export const LOEUF_INTOLERANT: Record<string, number> = { v4: 0.45, 'v2.1.1': 0.35 };
/**
 * gnomAD's general cut-off for a constrained gene: < 0.6 for v4.0 ("due to expected shifts in the LOEUF distribution
 * between gnomAD v2.1.1 and v4.0", gnomAD v4.0 gene constraint, March 2024), < 0.35 for v2.1.1.
 */
export const LOEUF_CONSTRAINED: Record<string, number> = { v4: 0.6, 'v2.1.1': 0.35 };
/** A gene's LoF constraint class: below the Mendelian threshold, below the general one only, above both, or not known. */
export type LoeufClass = 'intolerant' | 'constrained' | 'tolerant' | 'unknown';
export function loeufClass(c: Constraint | undefined): LoeufClass {
  if (!c || c.loeuf == null) return 'unknown';
  if (c.loeuf < (LOEUF_INTOLERANT[c.version] ?? LOEUF_INTOLERANT.v4)) return 'intolerant';
  if (c.loeuf < (LOEUF_CONSTRAINED[c.version] ?? LOEUF_CONSTRAINED.v4)) return 'constrained';
  return 'tolerant';
}
export const LOEUF_CLASS: Record<LoeufClass, { color: string; label: (v: string) => string }> = {
  intolerant: { color: '#b91c1c', label: v => `LOEUF < ${LOEUF_INTOLERANT[v] ?? LOEUF_INTOLERANT.v4}: LoF-intolerant` },
  constrained: { color: '#ea580c', label: v => `LOEUF ${LOEUF_INTOLERANT[v] ?? LOEUF_INTOLERANT.v4}–${LOEUF_CONSTRAINED[v] ?? LOEUF_CONSTRAINED.v4}: constrained` },
  tolerant: { color: '#475569', label: v => `LOEUF ≥ ${LOEUF_CONSTRAINED[v] ?? LOEUF_CONSTRAINED.v4}: not constrained` },
  unknown: { color: '#9aa6b5', label: () => 'no gnomAD constraint' },
};
/** pLI at or above which a gene is called LoF-intolerant (the convention gnomAD and UCSC cite). */
export const PLI_INTOLERANT = 0.9;

// ======================== Rebuilt chromosome in its own coordinates ========================

export interface DerPiece { piece: SvPiece; seg: SvSegment; derStart: number; derEnd: number; index: number }
export interface DerLayout { name: string; chrom: string; pieces: DerPiece[]; length: number; rearranged: boolean }

export function derLayouts(rec: SvReconstruction): DerLayout[] {
  return rec.ders.map(d => {
    let off = 0;
    const pieces = d.pieces.map((p, index) => {
      const len = p.seg.end - p.seg.start;
      const dp: DerPiece = { piece: p, seg: p.seg, derStart: off, derEnd: off + len, index };
      off += len;
      return dp;
    });
    return { name: d.name, chrom: d.chrom, pieces, length: off, rearranged: d.rearranged };
  });
}
/** A reference position of a piece's segment, on the rebuilt chromosome. */
export const toDer = (dp: DerPiece, refPos: number): number => (dp.piece.reversed ? dp.derStart + (dp.seg.end - refPos) : dp.derStart + (refPos - dp.seg.start));

// ======================== Fragments and chains ========================

export interface Fragment {
  tx: SvTranscript; der: number; dp: DerPiece;
  /** the part of the transcript inside the segment (reference coordinates) */
  lo: number; hi: number;
  /** strand on the rebuilt chromosome */
  strand: 1 | -1;
  /** the transcript goes on past the segment's left / right end (reference orientation), at a breakpoint or at the window's edge */
  cutLo: boolean; cutHi: boolean;
  /** the cut is the window's edge, not a breakpoint */
  edgeLo: boolean; edgeHi: boolean;
  /** exons inside (ranks), and the exon a breakpoint cuts, if any */
  ranks: number[]; cutExon: number | null;
}

/** Fragments of a transcript in transcription order: is the 5′ (or 3′) side of this fragment cut? */
const cut5 = (f: Fragment) => (f.tx.strand > 0 ? f.cutLo : f.cutHi);
const cut3 = (f: Fragment) => (f.tx.strand > 0 ? f.cutHi : f.cutLo);
const edge5 = (f: Fragment) => (f.tx.strand > 0 ? f.edgeLo : f.edgeHi);
const edge3 = (f: Fragment) => (f.tx.strand > 0 ? f.edgeHi : f.edgeLo);
/** On the rebuilt chromosome, from left to right: is the fragment cut at its left / right end, by a breakpoint? */
const derCutLeft = (f: Fragment) => (f.dp.piece.reversed ? f.cutHi && !f.edgeHi : f.cutLo && !f.edgeLo);
const derCutRight = (f: Fragment) => (f.dp.piece.reversed ? f.cutLo && !f.edgeLo : f.cutHi && !f.edgeHi);

function fragmentsOf(tx: SvTranscript, ders: DerLayout[], window: { start: number; end: number }): Fragment[] {
  const out: Fragment[] = [];
  ders.forEach((d, di) => {
    for (const dp of d.pieces) {
      const s = dp.seg;
      const lo = Math.max(tx.start, s.start), hi = Math.min(tx.end, s.end);
      if (hi <= lo) continue;
      const cutLo = tx.start < s.start, cutHi = tx.end > s.end;
      const ranks = tx.exons.filter(e => e.start < hi && e.end > lo).map(e => e.rank).sort((a, b) => a - b);
      const cutAt = [cutLo && s.start, cutHi && s.end].filter((x): x is number => typeof x === 'number');
      const ce = tx.exons.find(e => cutAt.some(p => e.start < p && e.end > p));
      out.push({ tx, der: di, dp, lo, hi, strand: (dp.piece.reversed ? -tx.strand : tx.strand) as 1 | -1, cutLo, cutHi,
        edgeLo: cutLo && s.first && s.start <= window.start, edgeHi: cutHi && s.last && s.end >= window.end, ranks, cutExon: ce ? ce.rank : null });
    }
  });
  return out;
}

export interface Chain {
  der: number;
  /** fragments in transcription order */
  parts: Fragment[];
  /** junctions crossed, in transcription order (null: a reference adjacency) */
  joins: (SvJunction | null)[];
}

/**
 * Chains of fragments on each rebuilt chromosome. At every junction, a fragment cut there is joined to the fragment
 * cut on the other side if both lie on the same strand of the rebuilt chromosome: the same transcript first (the gene
 * rejoined across an intragenic event), else the same gene, else any (a fusion).
 */
function chainsOf(frags: Fragment[], ders: DerLayout[]): Chain[] {
  const next = new Map<Fragment, Fragment>(), prev = new Map<Fragment, Fragment>();
  ders.forEach((d, di) => {
    for (let k = 0; k + 1 < d.pieces.length; k++) {
      const left = frags.filter(f => f.der === di && f.dp.index === k && derCutRight(f));
      const right = frags.filter(f => f.der === di && f.dp.index === k + 1 && derCutLeft(f));
      const pairs: [Fragment, Fragment, number][] = [];
      for (const a of left) for (const b of right) {
        if (a.strand !== b.strand) continue;
        pairs.push([a, b, a.tx.id === b.tx.id ? 0 : a.tx.gene === b.tx.gene ? 1 : 2]);
      }
      pairs.sort((x, y) => x[2] - y[2]);
      for (const [a, b] of pairs) if (!next.has(a) && !prev.has(b)) { next.set(a, b); prev.set(b, a); }
    }
  });
  const chains: Chain[] = [];
  for (const f of frags) {
    if (prev.has(f)) continue;
    const run: Fragment[] = [f];
    for (let x = next.get(f); x; x = next.get(x)) run.push(x);
    // left to right on the rebuilt chromosome; transcription runs right to left on its − strand
    const parts = run[0].strand > 0 ? run : [...run].reverse();
    const d = ders[f.der];
    const joins = parts.slice(1).map((p, i) => {
      const a = parts[i], b = p;
      const k = Math.max(a.dp.index, b.dp.index);
      const via = d.pieces[k].piece.via;
      return via && via !== 'ref' ? via : null;
    });
    chains.push({ der: f.der, parts, joins });
  }
  return chains;
}

// ======================== Reading a chain ========================

/** Coding bases of a transcript in exons of rank < `rank` (transcription order). */
function codingBefore(tx: SvTranscript, rank: number): number {
  if (tx.cdsStart == null || tx.cdsEnd == null) return 0;
  let n = 0;
  for (const e of tx.exons) if (e.rank < rank) n += Math.max(0, Math.min(e.end, tx.cdsEnd) - Math.max(e.start, tx.cdsStart));
  return n;
}
const codingIn = (tx: SvTranscript, ranks: number[]): number => {
  if (tx.cdsStart == null || tx.cdsEnd == null) return 0;
  let n = 0;
  for (const e of tx.exons) if (ranks.includes(e.rank)) n += Math.max(0, Math.min(e.end, tx.cdsEnd) - Math.max(e.start, tx.cdsStart));
  return n;
};
const cdsTotal = (tx: SvTranscript) => codingIn(tx, tx.exons.map(e => e.rank));
/** Ranks run: "exons 4–6", "exon 7", "exons 2, 5". */
export const ranksText = (r: number[]): string => {
  const s = [...new Set(r)].sort((a, b) => a - b);
  if (!s.length) return '';
  if (s.length === 1) return `exon ${s[0]}`;
  return s[s.length - 1] - s[0] === s.length - 1 ? `exons ${s[0]}–${s[s.length - 1]}` : `exons ${s.join(', ')}`;
};

export type ImpactKind = 'fusion' | 'disrupted' | 'intragenic-del' | 'deleted' | 'intragenic-dup' | 'partial-copy' | 'duplicated' | 'domain' | 'moved' | 'nearby' | 'unaffected';
export type FrameCall = 'kept' | 'shifted' | 'unknown' | 'none';

export interface FusionCall {
  five: SvTranscript; three: SvTranscript; junction: SvJunction | null;
  /** last exon of the 5′ partner and first exon of the 3′ partner kept */
  fiveLast: number; threeFirst: number;
  frame: FrameCall;
  /** what the product could be */
  kind: 'in-frame' | 'out-of-frame' | 'promoter-swap' | 'five-cds-whole' | 'uncertain';
  text: string;
}
export interface TxImpact {
  tx: SvTranscript;
  kind: ImpactKind;
  /** intact copies of the transcript on the rearranged homologue(s) */
  intact: number;
  /** functional copies in the sample (intact copies + the intact homologues), and in the reference */
  copies: number; refCopies: number;
  lost: number[]; repeated: number[];
  frame: FrameCall;
  fusions: FusionCall[];
  /** a breakpoint cuts it: where, and what goes where */
  cuts: { junction: SvJunction | null; rank: number | null; kept: number[]; text: string }[];
  /** the intact copy sits in an inverted or moved piece */
  moved: { reversed: boolean; segment: string }[];
  /** exons on no rebuilt chromosome (lost with a segment) */
  gone: number[];
  /** the transcript reaches past the window analysed: its far part is assumed unchanged */
  beyondWindow: boolean;
  text: string;
}

function fusionOf(chain: Chain): FusionCall | null {
  const p = chain.parts;
  for (let i = 0; i + 1 < p.length; i++) {
    const a = p[i], b = p[i + 1];
    if (a.tx.gene === b.tx.gene) continue;
    const five = a.tx, three = b.tx;
    // exons kept on each side of the junction (the 5′ partner's whole exons up to it, the 3′ partner's from it)
    const aRanks = p.slice(0, i + 1).filter(x => x.tx === five).flatMap(x => x.ranks);
    const bRanks = p.slice(i + 1).filter(x => x.tx === three).flatMap(x => x.ranks);
    const fiveLast = Math.max(...aRanks), threeFirst = Math.min(...bRanks);
    const exonCut = a.cutExon != null || b.cutExon != null;
    const fiveCoding = codingBefore(five, fiveLast + 1), fiveCds = cdsTotal(five);
    const threeBefore = codingBefore(three, threeFirst), threeCds = cdsTotal(three);
    let kind: FusionCall['kind'], frame: FrameCall, text: string;
    const name = `${five.gene}::${three.gene}`;
    const ex = `${five.gene} exons 1–${fiveLast} joined to ${three.gene} exons ${threeFirst}–${Math.max(...three.exons.map(e => e.rank))}`;
    if (exonCut) { kind = 'uncertain'; frame = 'unknown'; text = `${name}: a breakpoint lies inside exon ${a.cutExon ?? b.cutExon}; the splice sites used across the junction cannot be told from the DNA`; }
    else if (fiveCoding === 0) { kind = 'promoter-swap'; frame = 'none'; text = `${name}: the junction is in the 5′ UTR of ${five.gene}, so its promoter drives ${three.gene}${threeBefore === 0 && threeCds > 0 ? ` (whole coding sequence of ${three.gene}, from its own start codon)` : ` (from exon ${threeFirst}: ${threeCds - threeBefore} of its ${threeCds} coding bases)`}`; }
    else if (fiveCoding >= fiveCds) { kind = 'five-cds-whole'; frame = 'none'; text = `${name}: ${five.gene} keeps its whole coding sequence, the junction lies in its 3′ UTR: no fusion protein, its 3′ end changes`; }
    else if (threeBefore >= threeCds) { kind = 'out-of-frame'; frame = 'shifted'; text = `${name}: ${ex}, but past the stop codon of ${three.gene}: a truncated ${five.gene} with a new tail`; }
    else if (fiveCoding % 3 === threeBefore % 3) { kind = 'in-frame'; frame = 'kept'; text = `${name}: ${ex}, frame kept (${five.gene} ends ${fiveCoding} coding bases in, phase ${fiveCoding % 3}; ${three.gene} resumes ${threeBefore} bases in, phase ${threeBefore % 3}): an in-frame fusion protein is possible`; }
    else { kind = 'out-of-frame'; frame = 'shifted'; text = `${name}: ${ex}, frame shifted (phase ${fiveCoding % 3} into phase ${threeBefore % 3}): a truncated ${five.gene} with a short out-of-frame tail, NMD likely`; }
    return { five, three, junction: chain.joins[i], fiveLast, threeFirst, frame, kind, text };
  }
  return null;
}

/** Everything that happens to one transcript. */
function impactOf(tx: SvTranscript, chains: Chain[], ders: DerLayout[], copies: { flank: number; zygosity: 'het' | 'hom' }): TxImpact {
  const mine = chains.filter(c => c.parts.some(p => p.tx === tx));
  const all = tx.exons.map(e => e.rank).sort((a, b) => a - b);
  const n = all.length;
  let intact = 0;
  const lost = new Set<number>(), repeated = new Set<number>();
  const cuts: TxImpact['cuts'] = [], fusions: FusionCall[] = [], moved: TxImpact['moved'] = [];
  let frame: FrameCall = 'none';
  for (const c of mine) {
    const own = c.parts.filter(p => p.tx === tx);
    const fz = fusionOf(c);
    if (fz && (fz.five === tx || fz.three === tx)) { fusions.push(fz); continue; }
    const ranks = own.flatMap(p => p.ranks);
    const first = own[0], last = own[own.length - 1];
    const hasStart = !cut5(first) || edge5(first), hasEnd = !cut3(last) || edge3(last);
    const counts = new Map<number, number>(); for (const r of ranks) counts.set(r, (counts.get(r) ?? 0) + 1);
    const missing = all.filter(r => !counts.has(r) && r > Math.min(...ranks) && r < Math.max(...ranks));
    const twice = [...counts].filter(([, k]) => k > 1).map(([r]) => r);
    const whole = hasStart && hasEnd && !missing.length && !twice.length && own.every(p => p.cutExon == null);
    if (whole) {
      intact++;
      const dps = own.map(p => p.dp);
      const reversedAll = dps.every(dp => dp.piece.reversed);
      const movedPiece = dps.some(dp => !dp.seg.first && !dp.seg.last) || reversedAll;
      if (movedPiece) moved.push({ reversed: reversedAll, segment: dps.map(dp => `${dp.seg.letter}${dp.piece.reversed ? '′' : ''}`).join('') });
      continue;
    }
    if (hasStart && hasEnd && (missing.length || twice.length)) {
      missing.forEach(r => lost.add(r)); twice.forEach(r => repeated.add(r));
      const changed = codingIn(tx, missing) + codingIn(tx, twice);
      const exonCut = own.some(p => p.cutExon != null);
      frame = exonCut ? 'unknown' : changed === 0 ? 'none' : changed % 3 === 0 ? 'kept' : 'shifted';
      if (!twice.length || !missing.length) continue;
    }
    // the transcript stops at a breakpoint: what this part keeps
    const kept = [...new Set(ranks)].sort((a, b) => a - b);
    const brokenEnd = !hasEnd ? last : first;
    const rank = brokenEnd.cutExon;
    const j = c.joins.length ? c.joins[c.joins.length - 1] : null;
    const side = !hasStart && !hasEnd ? 'an inner part' : hasStart ? `the 5′ part (with its promoter)` : 'the 3′ part (no promoter)';
    cuts.push({ junction: j, rank, kept, text: `${side}: ${ranksText(kept) || 'no exon'}${rank != null ? ` (exon ${rank} cut)` : ''}` });
  }
  // gene never on a rebuilt chromosome, or only partly: the copies the sample keeps
  const homs = copies.zygosity === 'het' ? 1 : copies.flank, normals = copies.zygosity === 'het' ? copies.flank - 1 : 0;
  const total = normals + intact * homs;
  const onDer = mine.length > 0;
  const present = new Set(mine.flatMap(c => c.parts.filter(p => p.tx === tx).flatMap(p => p.ranks)));
  const gone = onDer ? all.filter(r => !present.has(r)) : [];
  const beyondWindow = mine.some(c => c.parts.some(p => p.tx === tx && (p.edgeLo || p.edgeHi)));
  let kind: ImpactKind;
  if (fusions.length) kind = 'fusion';
  else if (cuts.length && intact === 0) kind = 'disrupted';
  else if (cuts.length) kind = 'partial-copy';
  else if (lost.size && !repeated.size) kind = 'intragenic-del';
  else if (repeated.size && !lost.size) kind = 'intragenic-dup';
  else if (repeated.size && lost.size) kind = 'disrupted';
  else if (!onDer || (intact === 0 && !cuts.length)) kind = 'deleted';
  else if (intact > 1) kind = 'duplicated';
  else if (moved.length) kind = 'moved';
  else kind = 'unaffected';
  const t = describe(tx, kind, { intact, total, flank: copies.flank, lost: [...lost], repeated: [...repeated], frame, cuts, fusions, moved, gone, ders })
    + (beyondWindow ? ' It reaches past the window analysed: its far part is taken as unchanged.' : '');
  return { tx, kind, intact, copies: total, refCopies: copies.flank, lost: [...lost].sort((a, b) => a - b), repeated: [...repeated].sort((a, b) => a - b), frame, fusions, cuts, moved, gone, beyondWindow, text: t };
}

function describe(tx: SvTranscript, kind: ImpactKind, x: { intact: number; total: number; flank: number; lost: number[]; repeated: number[]; frame: FrameCall; cuts: TxImpact['cuts']; fusions: FusionCall[]; moved: TxImpact['moved']; gone: number[]; ders: DerLayout[] }): string {
  const cp = `${x.total} functional cop${x.total === 1 ? 'y' : 'ies'} of ${x.flank}`;
  const frameText = (bases: number) => x.frame === 'kept' ? `${bases} coding bases, frame kept (${bases / 3} codons)` : x.frame === 'shifted' ? `${bases} coding bases, not a multiple of 3: the frame shifts, premature stop and NMD expected` : x.frame === 'unknown' ? 'a breakpoint lies inside an exon: the frame cannot be told from the DNA' : 'no coding base';
  switch (kind) {
    case 'fusion': return x.fusions.map(f => f.text).join(' · ') + `. ${cp} of ${tx.gene} itself.`;
    case 'disrupted': return `Cut by the rearrangement: ${x.cuts.map(c => c.text).join('; ')}${x.gone.length ? `; ${ranksText(x.gone)} lost with the deleted segment` : ''}. No intact copy on the rearranged homologue: ${cp}.`;
    case 'partial-copy': return `An intact copy remains, plus a partial copy: ${x.cuts.map(c => c.text).join('; ')}${x.cuts.some(c => c.text.startsWith('the 5′')) ? ' (a truncated product may be made from it)' : ' (no promoter: not expected to be expressed)'}. ${cp}.`;
    case 'intragenic-del': return `${ranksText(x.lost)} deleted: ${frameText(codingIn(tx, x.lost))}.`;
    case 'intragenic-dup': return `${ranksText(x.repeated)} repeated in tandem: ${frameText(codingIn(tx, x.repeated))}.`;
    case 'deleted': return `Whole gene on a lost segment: ${cp}.`;
    case 'duplicated': return `Whole gene present ${x.intact} times on the rebuilt chromosome: ${cp}.`;
    case 'moved': return `Coding sequence intact, ${x.moved.map(m => m.reversed ? `inverted (segment ${m.segment})` : `moved (segment ${m.segment})`).join(', ')}: a position effect cannot be judged without chromatin domains.`;
    case 'domain': return 'Coding sequence intact; its chromatin domain changes.';
    default: return 'Not changed by the rearrangement.';
  }
}

// ======================== Chromatin domains ========================

export interface DomainPart { domain: number; lo: number; hi: number; der: number; dp: DerPiece; cutLeft: boolean; cutRight: boolean }
export interface DerDomain {
  der: number; derStart: number; derEnd: number;
  parts: DomainPart[];
  /** made of parts of more than one reference domain: a new domain */
  neo: boolean;
}
export interface DomainChange {
  gene: string;
  /** reference domain of the gene's start (transcription start) */
  domain: number | null;
  /** bases of other reference domains the gene now shares a domain with, and the genes there */
  gained: number; partners: string[];
  /** bases of its own domain it no longer shares a domain with, and the genes there */
  lostBases: number; separated: string[];
  text: string;
}

/** Domains of the window laid onto the rebuilt chromosomes; parts joined at a junction become one domain. */
export function derDomains(domains: Domain[], ders: DerLayout[]): DerDomain[] {
  const out: DerDomain[] = [];
  ders.forEach((d, di) => {
    let cur: DerDomain | null = null;
    for (const dp of d.pieces) {
      const parts = domains.map((dom, k) => ({ dom, k })).filter(({ dom }) => dom.end > dp.seg.start && dom.start < dp.seg.end).map(({ dom, k }) => {
        const lo = Math.max(dom.start, dp.seg.start), hi = Math.min(dom.end, dp.seg.end);
        const cutL = dp.piece.reversed ? dom.end > dp.seg.end : dom.start < dp.seg.start, cutR = dp.piece.reversed ? dom.start < dp.seg.start : dom.end > dp.seg.end;
        return { domain: k, lo, hi, der: di, dp, cutLeft: cutL, cutRight: cutR };
      });
      // left to right on the rebuilt chromosome
      parts.sort((a, b) => toDer(dp, dp.piece.reversed ? a.hi : a.lo) - toDer(dp, dp.piece.reversed ? b.hi : b.lo));
      for (const p of parts) {
        const s = toDer(dp, dp.piece.reversed ? p.hi : p.lo), e = toDer(dp, dp.piece.reversed ? p.lo : p.hi);
        // a part cut at its left end continues the domain open on the left, when that one is cut at its right end, right there
        if (cur && p.cutLeft && cur.parts[cur.parts.length - 1].cutRight && Math.abs(cur.derEnd - s) <= 1) {
          cur.parts.push(p); cur.derEnd = e; cur.neo = cur.neo || cur.parts.some(q => q.domain !== p.domain);
        } else { cur = { der: di, derStart: s, derEnd: e, parts: [p], neo: false }; out.push(cur); }
      }
    }
  });
  return out;
}

function domainChanges(genes: SvGene[], domains: Domain[], dd: DerDomain[], impacts: Map<string, TxImpact>): Map<string, DomainChange> {
  const out = new Map<string, DomainChange>();
  const tss = (g: SvGene) => (g.principal.strand > 0 ? g.principal.start : g.principal.end - 1);
  const domainOf = (pos: number) => domains.findIndex(d => pos >= d.start && pos < d.end);
  const genesIn = (k: number, lo: number, hi: number, self: string) => genes.filter(o => o.symbol !== self && tss(o) >= lo && tss(o) < hi && domainOf(tss(o)) === k).map(o => o.symbol);
  for (const g of genes) {
    const imp = impacts.get(g.symbol);
    if (!imp || (imp.kind !== 'moved' && imp.kind !== 'unaffected' && imp.kind !== 'duplicated')) continue;
    const pos = tss(g), k = domainOf(pos);
    if (k < 0) continue;
    // the rebuilt domain(s) holding the gene's start
    const homes = dd.filter(x => x.parts.some(p => p.domain === k && pos >= p.lo && pos < p.hi));
    if (!homes.length) continue;
    let gained = 0, kept = 0;
    const partners = new Set<string>();
    for (const h of homes) for (const p of h.parts) {
      if (p.domain !== k) { gained += p.hi - p.lo; genesIn(p.domain, p.lo, p.hi, g.symbol).forEach(s => partners.add(s)); }
      else kept = Math.max(kept, p.hi - p.lo);
    }
    const own = domains[k].end - domains[k].start;
    const lostBases = Math.max(0, own - kept);
    const separated = genes.filter(o => o.symbol !== g.symbol && domainOf(tss(o)) === k && !homes.some(h => h.parts.some(p => p.domain === k && tss(o) >= p.lo && tss(o) < p.hi))).map(o => o.symbol);
    if (!gained && !lostBases) continue;
    const kb = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} Mb` : `${Math.round(n / 1000)} kb`);
    const text = [
      gained ? `now shares a domain with ${kb(gained)} from ${gained && partners.size ? `another domain (genes ${[...partners].join(', ')})` : 'another domain'}: enhancer adoption possible` : '',
      lostBases ? `its domain loses ${kb(lostBases)}${separated.length ? ` (with ${separated.join(', ')})` : ''}: its own enhancers may be gone` : '',
    ].filter(Boolean).join('; ');
    out.set(g.symbol, { gene: g.symbol, domain: k, gained, partners: [...partners], lostBases, separated, text });
  }
  return out;
}

// ======================== Genes ========================

export interface GeneImpact {
  gene: SvGene;
  /** the nearest breakpoint of the rearrangement, for an intact gene (bp from the gene's span; 0 when it lies inside) */
  nearest: { distance: number; pos: number; junctions: string[] } | null;
  principal: TxImpact;
  /** every NM_ transcript, when they differ from the principal one */
  others: TxImpact[];
  kind: ImpactKind;
  domain: DomainChange | null;
  /** severity for sorting */
  rank: number;
  fragments: Fragment[];
}
export interface GeneStory {
  ders: DerLayout[];
  genes: GeneImpact[];
  chains: Chain[];
  fragments: Fragment[];
  domains: Domain[];
  derDomains: DerDomain[];
}

const SEVERITY: Record<ImpactKind, number> = { fusion: 8, disrupted: 7, deleted: 6, 'intragenic-del': 6, duplicated: 5, 'intragenic-dup': 5, 'partial-copy': 4, domain: 3, moved: 2, nearby: 1, unaffected: 0 };

/**
 * Without chromatin domains, an intact gene within this distance of a breakpoint is reported as a possible position
 * effect. 1 Mb by default: TADs average about 0.85–1.1 Mb in Hi-C (Dixon et al. 2012), and position effects have been
 * documented from breakpoints 1.3 Mb away (SOX9, campomelic dysplasia) and further.
 */
export const NEAR_BREAKPOINT_DEFAULT = 1_000_000;

export function geneStory(rec: SvReconstruction, genes: SvGene[], window: { start: number; end: number }, opts: { flankCopies: number; zygosity: 'het' | 'hom'; domains: Domain[]; near?: number }): GeneStory {
  const ders = derLayouts(rec);
  const txs = genes.flatMap(g => g.transcripts);
  const fragments = txs.flatMap(tx => fragmentsOf(tx, ders, window));
  const chains = chainsOf(fragments, ders);
  const out: GeneImpact[] = [];
  const principalImpacts = new Map<string, TxImpact>();
  for (const g of genes) {
    const imps = g.transcripts.map(tx => impactOf(tx, chains, ders, { flank: opts.flankCopies, zygosity: opts.zygosity }));
    const p = imps.find(x => x.tx === g.principal)!;
    principalImpacts.set(g.symbol, p);
    const others = imps.filter(x => x !== p && (x.kind !== p.kind || x.frame !== p.frame || x.fusions.map(f => f.kind).join() !== p.fusions.map(f => f.kind).join()));
    out.push({ gene: g, principal: p, others, kind: p.kind, domain: null, nearest: null, rank: SEVERITY[p.kind], fragments: fragments.filter(f => f.tx.gene === g.symbol) });
  }
  // the nearest breakpoint of every gene, with the junctions meeting there
  const bps = new Map<number, string[]>();
  for (const j of rec.junctions) for (const e of [j.a, j.b]) if (e.chrom === rec.chrs[0]) bps.set(e.pos, [...new Set([...(bps.get(e.pos) ?? []), j.id])]);
  for (const gi of out) {
    let best: GeneImpact['nearest'] = null;
    for (const [pos, js] of bps) {
      const d = pos < gi.gene.start ? gi.gene.start - pos : pos > gi.gene.end ? pos - gi.gene.end : 0;
      if (!best || d < best.distance) best = { distance: d, pos, junctions: js };
    }
    gi.nearest = best;
  }
  const domains = opts.domains.filter(d => d.end > window.start && d.start < window.end).sort((a, b) => a.start - b.start);
  const dd = domains.length ? derDomains(domains, ders) : [];
  if (domains.length) {
    const changes = domainChanges(genes, domains, dd, principalImpacts);
    for (const gi of out) {
      const c = changes.get(gi.gene.symbol);
      if (!c) continue;
      gi.domain = c;
      if (gi.kind === 'unaffected' || gi.kind === 'moved') {
        gi.kind = 'domain'; gi.rank = SEVERITY.domain;
        gi.principal = { ...gi.principal, text: `Coding sequence intact${gi.principal.moved.length ? `, ${gi.principal.moved.map(m => (m.reversed ? `inverted (segment ${m.segment})` : `moved (segment ${m.segment})`)).join(', ')}` : ''}; its chromatin domain changes: ${c.text}.` };
      }
    }
  }
  // no domains: an intact, unmoved gene near a breakpoint may still be misregulated
  if (!domains.length && opts.near != null && opts.near > 0) {
    for (const gi of out) {
      if (gi.kind !== 'unaffected' || !gi.nearest || gi.nearest.distance > opts.near) continue;
      const side = gi.nearest.pos < gi.gene.start ? (gi.gene.strand > 0 ? 'upstream' : 'downstream') : gi.nearest.pos > gi.gene.end ? (gi.gene.strand > 0 ? 'downstream' : 'upstream') : 'inside its span';
      const dist = gi.nearest.distance >= 1e6 ? `${(gi.nearest.distance / 1e6).toFixed(2)} Mb` : `${Math.round(gi.nearest.distance / 1000)} kb`;
      gi.kind = 'nearby'; gi.rank = SEVERITY.nearby;
      gi.principal = { ...gi.principal, text: `Coding sequence intact, ${dist} ${side} of the breakpoint at ${(gi.nearest.pos + 1).toLocaleString('en-US')} (${gi.nearest.junctions.join(', ')}): a position effect is possible (its enhancers or their insulation may have changed); load chromatin domains to judge it.` };
    }
  }
  out.sort((a, b) => b.rank - a.rank || a.gene.start - b.gene.start);
  return { ders, genes: out, chains, fragments, domains, derDomains: dd };
}

// ======================== Dosage reading ========================

const HI_TEXT: Record<number, string> = { 0: 'no evidence', 1: 'little evidence', 2: 'some evidence', 3: 'sufficient evidence', 30: 'autosomal recessive', 40: 'dosage sensitivity unlikely' };
export const dosageScoreText = (s: number | null | undefined) => (s == null ? 'not curated' : HI_TEXT[s] ?? String(s));

/** Does the change of dosage matter for this gene: one line from ClinGen and gnomAD, or null when nothing is known. */
export function dosageReading(gi: GeneImpact, a: GeneAnnotations): { level: 'high' | 'moderate' | 'low' | 'unknown'; text: string } {
  const sym = gi.gene.symbol, d = a.dosage.get(sym), c = a.constraint.get(sym);
  const p = gi.principal;
  const loss = p.copies < p.refCopies || gi.kind === 'disrupted' || (gi.kind === 'intragenic-del' && p.frame !== 'kept') || (gi.kind === 'intragenic-dup' && p.frame === 'shifted') || p.fusions.some(f => f.kind === 'out-of-frame');
  const gain = gi.kind === 'duplicated';
  const thr = c ? LOEUF_INTOLERANT[c.version] ?? LOEUF_INTOLERANT.v4 : LOEUF_INTOLERANT.v4;
  const intolerant = c?.loeuf != null && c.loeuf < thr;
  const parts: string[] = [];
  let level: 'high' | 'moderate' | 'low' | 'unknown' = 'unknown';
  if (loss) {
    if (d?.hi === 3) { level = 'high'; parts.push('ClinGen: haploinsufficient (3, sufficient evidence): a loss of one copy is a known disease mechanism'); }
    else if (d?.hi === 2) { level = 'moderate'; parts.push('ClinGen haploinsufficiency 2 (some evidence)'); }
    else if (d?.hi === 30) { level = 'low'; parts.push('ClinGen: autosomal recessive gene: one lost copy matters with a variant on the other allele'); }
    else if (d?.hi === 40) { level = 'low'; parts.push('ClinGen: dosage sensitivity unlikely'); }
    else if (d?.hi != null) parts.push(`ClinGen haploinsufficiency ${d.hi} (${dosageScoreText(d.hi)})`);
    if (c?.loeuf != null) { parts.push(`gnomAD LOEUF ${c.loeuf.toFixed(2)}${intolerant ? ` (< ${thr}: LoF-intolerant)` : ''}`); if (intolerant && level === 'unknown') level = 'moderate'; if (!intolerant && level === 'unknown') level = 'low'; }
  } else if (gain) {
    if (d?.ts === 3) { level = 'high'; parts.push('ClinGen: triplosensitive (3, sufficient evidence)'); }
    else if (d?.ts != null) { parts.push(`ClinGen triplosensitivity ${d.ts} (${dosageScoreText(d.ts)})`); level = d.ts === 2 ? 'moderate' : 'low'; }
    else parts.push('ClinGen triplosensitivity not curated');
  } else if (gi.kind === 'nearby' || gi.kind === 'domain' || gi.kind === 'moved') {
    // intact but repositioned: what a change of its expression would mean
    if (d?.hi === 3) { level = 'moderate'; parts.push('haploinsufficient gene (ClinGen 3): a position effect lowering its expression could matter'); }
    if (d?.ts === 3) { level = 'moderate'; parts.push('triplosensitive gene (ClinGen 3): a position effect raising its expression could matter'); }
    if (!parts.length && intolerant) { level = 'low'; parts.push(`LoF-intolerant (gnomAD LOEUF ${c!.loeuf!.toFixed(2)}): sensitive to a loss of expression`); }
    if (!parts.length) return { level: 'unknown', text: '' };
  } else return { level: 'unknown', text: '' };
  return { level, text: parts.join(' · ') || 'no dosage curation, no constraint data' };
}
