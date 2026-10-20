/**
 * Reading frame across a fusion junction: the 5′ partner's coding sequence up to the junction (its donor), then the
 * 3′ partner's transcript from the junction on (its acceptor).
 *
 * - From the models alone (fusionFrame): the 5′ partner keeps `aBases` coding bases from its start codon; the
 *   junction lies `bOffset` bases into the 3′ partner's CDS (negative in its 5′ UTR). The 3′ partner's protein follows
 *   in frame when the two agree modulo 3, the rule fusion callers apply (STAR-Fusion, Arriba, FusionCatcher).
 * - From the bases (fusionPeptide): the fused sequence translated in the 5′ partner's frame from a few codons before
 *   the junction, to its first stop; the peptide after the junction is compared with the 3′ partner's own translation.
 *   This sees what the lengths cannot: a stop codon in the 3′ partner's 5′ UTR before its start codon, or one formed
 *   by the junction codon itself.
 */
import { cdnaPosition, translateCodon, type TxModel } from './geometry';

/** One end of a fusion junction: the model, and the exon boundary the junction is on (genomic, 0-based). */
export interface FusionEnd { model: TxModel; pos: number }

export type FusionFrameVerdict =
  | 'in-frame' | 'frameshift'
  /** the 5′ partner's junction lies before its start codon: no protein of it in the fusion */
  | 'five-utr'
  /** the 5′ partner's junction lies after its stop codon: its whole protein, then untranslated sequence */
  | 'after-stop'
  /** the 3′ partner's junction lies after its stop codon: the 5′ protein runs into its 3′ UTR */
  | 'three-utr'
  /** a partner without a coding sequence, or a boundary off its exons */
  | 'noncoding';

export interface FusionFrame {
  verdict: FusionFrameVerdict;
  /** coding bases of the 5′ partner kept, from its start codon (null: none) */
  aBases: number | null;
  /** position of the junction in the 3′ partner's CDS, 0-based (negative: in its 5′ UTR; null: no CDS) */
  bOffset: number | null;
  /** codons of the 5′ partner's protein kept (whole codons, the junction codon left out) */
  aCodons: number | null;
  /** the 3′ partner's codon the fusion continues with when in frame (1-based, as p. numbering) */
  bCodon: number | null;
  text: string;
}

interface Layout { tOf: (pos: number) => number | null; posOf: (t: number) => number | null; cdsStart: number | null; cdsEnd: number | null; length: number; plus: boolean }

const layouts = new WeakMap<TxModel, Layout>();
/** Transcript coordinates (1-based, 5′→3′) of a model: exonic positions only (kept with the model, which does not change). */
export function layoutOf(m: TxModel): Layout {
  let L = layouts.get(m);
  if (!L) layouts.set(m, L = buildLayout(m));
  return L;
}
function buildLayout(m: TxModel): Layout {
  const plus = m.strand > 0;
  const order = (plus ? [...m.exons] : [...m.exons].reverse()).map(e => ({ start: e.start, end: e.end, t: 0 }));
  let t = 1;
  for (const e of order) { e.t = t; t += e.end - e.start; }
  const tOf = (pos: number): number | null => {
    for (const e of order) if (pos >= e.start && pos < e.end) return e.t + (plus ? pos - e.start : e.end - 1 - pos);
    return null;
  };
  const posOf = (tt: number): number | null => {
    for (const e of order) { const len = e.end - e.start; if (tt >= e.t && tt < e.t + len) return plus ? e.start + (tt - e.t) : e.end - 1 - (tt - e.t); }
    return null;
  };
  const coding = m.cdsStart != null && m.cdsEnd != null && m.cdsEnd > m.cdsStart;
  return {
    tOf, posOf, plus, length: t - 1,
    cdsStart: coding ? tOf(plus ? m.cdsStart! : m.cdsEnd! - 1) : null,
    cdsEnd: coding ? tOf(plus ? m.cdsEnd! - 1 : m.cdsStart!) : null,
  };
}

/** The last transcript base of the 5′ partner before the junction, the first of the 3′ partner after it. */
function lastBefore(a: FusionEnd, L: Layout): number | null { return L.tOf(L.plus ? a.pos - 1 : a.pos); }
function firstAfter(b: FusionEnd, L: Layout): number | null { return L.tOf(L.plus ? b.pos : b.pos - 1); }

const mod3 = (n: number) => ((n % 3) + 3) % 3;

/** The frame of a fusion of `five` (donor end) to `three` (acceptor end), from the models. */
export function fusionFrame(five: FusionEnd, three: FusionEnd): FusionFrame {
  const A = layoutOf(five.model), B = layoutOf(three.model);
  const ta = lastBefore(five, A), tb = firstAfter(three, B);
  const g5 = five.model.geneName, g3 = three.model.geneName;
  const none = (verdict: FusionFrameVerdict, text: string, extra: Partial<FusionFrame> = {}): FusionFrame => ({ verdict, aBases: null, bOffset: null, aCodons: null, bCodon: null, text, ...extra });
  if (ta == null || tb == null) return none('noncoding', 'a junction end is off its partner\'s exons: no reading frame');
  if (A.cdsStart == null || A.cdsEnd == null) return none('noncoding', `${g5} has no coding sequence: no fusion protein from its start codon`);
  if (ta < A.cdsStart) return none('five-utr', `the junction is in ${g5}'s 5′ UTR, before its start codon: no ${g5} protein in the fusion (${g3}'s own start codon may be used)`);
  if (ta >= A.cdsEnd) return none('after-stop', `the junction is after ${g5}'s stop codon: ${g5}'s whole protein, then untranslated sequence`);
  const aBases = ta - A.cdsStart + 1, aCodons = Math.floor(aBases / 3);
  if (B.cdsStart == null || B.cdsEnd == null) return none('noncoding', `${g3} has no coding sequence: ${g5}'s first ${aCodons} codons, then ${g3}'s transcript read in ${g5}'s frame`, { aBases, aCodons });
  const bOffset = tb - B.cdsStart;
  if (tb > B.cdsEnd) return none('three-utr', `the junction is after ${g3}'s stop codon: ${g5}'s first ${aCodons} codons, then ${g3}'s 3′ UTR read in ${g5}'s frame`, { aBases, bOffset, aCodons });
  const inFrame = mod3(aBases - bOffset) === 0;
  // the 3′ partner's codon the fused protein continues with: the first whole codon of it after the junction
  const bCodon = inFrame && bOffset >= 0 ? Math.floor((bOffset + mod3(-aBases)) / 3) + 1 : null;
  if (!inFrame) {
    return { verdict: 'frameshift', aBases, bOffset, aCodons, bCodon: null,
      text: `frameshift: ${g5} keeps ${aBases.toLocaleString('en-US')} coding bases (${aBases % 3 ? `codon ${aCodons + 1} split, ` : ''}phase ${aBases % 3}), ${g3} continues from ${bOffset < 0 ? `${-bOffset} bases before its start codon` : `coding base ${bOffset + 1}`} (phase ${mod3(bOffset)}): ${g3}'s protein is read in another frame` };
  }
  return { verdict: 'in-frame', aBases, bOffset, aCodons, bCodon,
    text: bOffset < 0
      ? `in frame with ${g3}'s start codon: ${g5} p.1–${aCodons}, then ${-bOffset} bases of ${g3}'s 5′ UTR and its whole protein (unless a stop codon lies in that UTR)`
      : `in frame: ${g5} p.1–${aCodons}${aBases % 3 ? ' (its last codon joined across the junction)' : ''} fused to ${g3} from p.${bCodon}` };
}

/**
 * The fusion transcript in HGVS r. notation: the 5′ partner from its first base to the junction, `::`, the 3′ partner
 * from the junction to its last base (e.g. NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924).
 */
export function fusionHgvs(five: FusionEnd, three: FusionEnd): string {
  const ends = (m: TxModel) => { const ex = [...m.exons].sort((x, y) => x.start - y.start); return m.strand > 0 ? [ex[0].start, ex[ex.length - 1].end - 1] : [ex[ex.length - 1].end - 1, ex[0].start]; };
  const r = (g: number, m: TxModel) => cdnaPosition(g, m).bare;
  const A = five.model, B = three.model, pa = A.strand > 0, pb = B.strand > 0;
  return `${A.transcriptId}:r.${r(ends(A)[0], A)}_${r(pa ? five.pos - 1 : five.pos, A)}::${B.transcriptId}:r.${r(pb ? three.pos : three.pos - 1, B)}_${r(ends(B)[1], B)}`;
}

/** The transcript bases [from, to] (1-based, inclusive) of a model, from its genomic bases. */
export async function transcriptBases(m: TxModel, from: number, to: number, getRef: (start: number, end: number, chrom: string) => Promise<string | null>): Promise<string | null> {
  const L = layoutOf(m);
  from = Math.max(1, from); to = Math.min(L.length, to);
  if (to < from) return '';
  const p1 = L.posOf(from)!, p2 = L.posOf(to)!;
  const lo = Math.min(p1, p2), hi = Math.max(p1, p2) + 1;
  // the exons over [lo, hi), each read on its own: an intron of 100 kb is not fetched
  let out = '';
  const parts = (L.plus ? [...m.exons] : [...m.exons].reverse());
  for (const e of parts) {
    const a = Math.max(e.start, lo), b = Math.min(e.end, hi);
    if (b <= a) continue;
    const s = await getRef(a, b, m.chrom);
    if (s == null || s.length !== b - a) return null;
    out += L.plus ? s.toUpperCase() : revcomp(s.toUpperCase());
  }
  return out;
}

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const revcomp = (s: string) => { let r = ''; for (let i = s.length - 1; i >= 0; i--) r += COMP[s[i]] ?? 'N'; return r; };
const translate = (s: string) => { let p = ''; for (let i = 0; i + 3 <= s.length; i += 3) p += translateCodon(s.slice(i, i + 3)); return p; };

export interface FusionPeptide {
  /** the 5′ partner's last codons before the junction */
  before: string;
  /** the codon across the junction (bases of both partners), when the junction splits a codon */
  junction: string | null;
  /** the fused protein after the junction (the junction codon left out), to its first stop (`*` included) or the end of what was read */
  after: string;
  /** codons after the 5′ partner's last whole codon to the first stop, the junction codon included (null: none in what was read) */
  stopAfter: number | null;
  /** the stop codon is the 3′ partner's own (in frame: its protein is carried to its end; read too short: null) */
  carries: boolean | null;
  text: string;
}

/**
 * The fused sequence translated across the junction: `beforeCodons` codons of the 5′ partner, then the 3′ partner read
 * in that frame for up to `afterBases` bases, to the first stop codon; for an in-frame fusion, whether that stop is the
 * 3′ partner's own (its protein carried to its end) or a premature one (in its 5′ UTR, or the junction codon).
 */
export async function fusionPeptide(five: FusionEnd, three: FusionEnd, frame: FusionFrame, getRef: (start: number, end: number, chrom: string) => Promise<string | null>,
  beforeCodons = 10, afterBases = 6000): Promise<FusionPeptide | null> {
  if (frame.aBases == null) return null;
  const A = layoutOf(five.model), B = layoutOf(three.model);
  const ta = lastBefore(five, A)!, tb = firstAfter(three, B)!;
  const rest = frame.aBases % 3, keep = Math.min(frame.aCodons ?? 0, beforeCodons);
  const aSeq = await transcriptBases(five.model, ta - rest - 3 * keep + 1, ta, getRef);
  // read on to the 3′ partner's own stop codon when it lies further (a long CDS), so an in-frame fusion can be judged
  if (B.cdsEnd != null) afterBases = Math.max(afterBases, Math.min(150000, B.cdsEnd - tb + 4));
  const bSeq = await transcriptBases(three.model, tb, tb + afterBases - 1, getRef);
  if (aSeq == null || bSeq == null) return null;
  const before = translate(aSeq.slice(0, 3 * keep));
  // from the 5′ partner's last whole codon on: the junction codon (when split), then the 3′ partner in that frame
  const fused = translate(aSeq.slice(3 * keep) + bSeq);
  const junction = rest ? fused[0] ?? null : null;
  const stop = fused.indexOf('*');
  const stopAfter = stop >= 0 ? stop + 1 : null;
  const after = (stop >= 0 ? fused.slice(0, stop + 1) : fused).slice(rest ? 1 : 0);
  const readAll = tb + afterBases - 1 >= B.length;
  let carries: boolean | null = null;
  if (frame.verdict === 'in-frame' && B.cdsStart != null && B.cdsEnd != null && frame.bOffset != null) {
    // the transcript base of the 3′ partner where the fused reading reaches its stop codon: its own stop when carried
    const stopBase = stop >= 0 ? tb + 3 * stop - rest : null;   // first base of the stop codon, in the 3′ partner's coordinates
    const ownStop = B.cdsEnd - 2;
    carries = stopBase != null ? stopBase === ownStop : readAll ? false : null;
  }
  const g3 = three.model.geneName, g5 = five.model.geneName;
  const nOwn = B.cdsStart != null && B.cdsEnd != null ? Math.round((B.cdsEnd - B.cdsStart + 1) / 3) : 0;
  const text = frame.verdict === 'in-frame'
    ? carries === true ? `translated across the junction, the fused protein reads on into ${g3}'s and ends at ${g3}'s own stop codon (p.${nOwn}*)`
      : carries === false ? (junction === '*' ? 'in frame by the lengths, but the codon across the junction is a stop: the protein ends there'
        : `in frame by the lengths, but a stop codon ${stopAfter} codons after ${g5}'s part ends the protein before ${g3}'s own sequence${(frame.bOffset ?? 0) < 0 ? ` (in ${g3}'s 5′ UTR)` : ''}`)
      : `in frame; no stop codon within the ${afterBases.toLocaleString('en-US')} bases of ${g3} read`
    : stopAfter != null ? `translated across the junction, the fused protein ends ${stopAfter} codon${stopAfter === 1 ? '' : 's'} after ${g5}'s part (a stop in ${g3}'s sequence read out of its frame)`
      : `no stop codon within the ${afterBases.toLocaleString('en-US')} bases of ${g3} read`;
  return { before, junction, after: after.length > 30 ? `${after.slice(0, 30)}…` : after, stopAfter, carries, text };
}
