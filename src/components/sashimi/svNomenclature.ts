/**
 * Nomenclature of a rebuilt structural variant, and the reverse: the rearrangement a typed description implies.
 *
 * HGVS (the variant, on the reference): following the HGVS recommendations' complex (HGVS/ISCN) page, a rearranged
 * chromosome is described in the allele form, its changes listed from pter to qter on reference positions
 * (`NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]` for an inversion with a deletion at one
 * breakpoint); a tandem copy is a `dup`, any other copy or moved piece an insertion (`g.234_235ins123_234inv`), and a
 * chromosome end replaced by another chromosome a `delins[NC_…:g.…_qter]`. Positions are the best ones the reads give,
 * not shifted by the 3′ rule: their uncertainty goes along as "±N bp", outside the description.
 *
 * ISCN (the resulting chromosome): the bands of the breakpoints, `seq[GRCh38] del(20)(p12.2p12.2)` for a single
 * deletion, duplication or inversion, otherwise the detailed form of the derivative from pter to qter
 * (`der(20)(pter->p12.2::p12.2->p12.2::p12.2->qter)`), and its sequence form with the positions
 * (`chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::chr20:g.10,808,806_cen_qter`).
 *
 * Every description is a draft: the viewer says so next to it.
 */
import type { SvJunction, SvReconstruction, SvSegment, SvPiece } from './svModel';
import { ncAccession, chromOfAccession } from './knownVariants';

export interface Cytoband { start: number; end: number; name: string; stain: string }
export interface SvNomenclature {
  /** one HGVS description per rearranged chromosome */
  hgvs: string[];
  /** ISCN with bands (when the cytobands are known), and the sequence form */
  iscn: string | null; iscnSeq: string[];
  /** "±5–45 bp" over the breakpoints described */
  uncertainty: string;
}
export const NOMENCLATURE_WARNING = 'Automated draft: check before any clinical report.';

const num = (c: string) => c.replace(/^chr/i, '');
const comma = (n: number) => n.toLocaleString('en-US');

/** A piece of sequence on a derivative: 0-based half-open, `pter`/`qter` when it reaches the chromosome end. */
interface Span { chrom: string; start: number; end: number; reversed: boolean; pter: boolean; qter: boolean }
const spanOf = (p: SvPiece): Span => ({ chrom: p.seg.chrom, start: p.seg.start, end: p.seg.end, reversed: p.reversed, pter: p.seg.first && p.seg.i === 0, qter: p.seg.last });

/** `10489766_10783340` (1-based) with pter / qter at the chromosome ends. */
const range = (s: { start: number; end: number; pter?: boolean; qter?: boolean }, sep = (n: number) => String(n)) =>
  `${s.pter ? 'pter' : sep(s.start + 1)}_${s.qter ? 'qter' : sep(s.end)}`;

/**
 * The HGVS description of one rebuilt chromosome. The pieces kept forward in the reference order form the backbone
 * (longest such run); what lies between two backbone pieces on the reference and on the derivative gives the change:
 * nothing there → del; the same segments, some reversed in place → del / inv; a forward repeat of what precedes →
 * dup; anything else → ins (between two backbone pieces) or delins (replacing reference segments).
 */
function hgvsOf(der: SvReconstruction['ders'][number], segs: SvSegment[], build: string): string | null {
  const own = der.chrom;
  const acc = ncAccession(build, own) ?? own;
  const chromSegs = segs.filter(s => s.chrom === own);
  const pieces = der.pieces;
  // backbone: longest increasing run (by segment index) of forward, first-visit pieces of this chromosome
  const cand = pieces.map((p, k) => ({ p, k })).filter(x => x.p.seg.chrom === own && !x.p.reversed && x.p.copy === 0);
  const best: number[] = cand.map(() => 1), prev: number[] = cand.map(() => -1);
  for (let i = 0; i < cand.length; i++) for (let j = 0; j < i; j++) if (cand[j].p.seg.i < cand[i].p.seg.i && best[j] + 1 > best[i]) { best[i] = best[j] + 1; prev[i] = j; }
  let end = best.indexOf(Math.max(...best));
  const bb: { p: SvPiece; k: number }[] = [];
  while (end >= 0) { bb.unshift(cand[end]); end = prev[end]; }
  if (!bb.length || bb[0].p.seg.i !== 0) return null;
  const ops: string[] = [];
  const content = (run: SvPiece[]): string => {
    const items = run.map(spanOf).map(s => `${s.chrom === own ? '' : `${ncAccession(build, s.chrom) ?? s.chrom}:g.`}${range(s)}${s.reversed ? 'inv' : ''}`);
    return items.length === 1 && run[0].seg.chrom === own ? items[0] : `[${items.join(';')}]`;
  };
  const gapOps = (left: SvPiece, right: SvPiece | null, run: SvPiece[]) => {
    const gap = chromSegs.filter(s => s.i > left.seg.i && (right ? s.i < right.seg.i : true));
    if (!run.length) {
      if (gap.length) ops.push(`${range({ start: gap[0].start, end: gap[gap.length - 1].end, qter: !right })}del`);
      return;
    }
    const inGap = new Set(gap.map(s => s.id));
    // the same segments as the gap, some reversed in place
    if (right && run.every(p => p.seg.chrom === own && inGap.has(p.seg.id)) && new Set(run.map(p => p.seg.id)).size === run.length) {
      const blocks: SvPiece[][] = [];
      for (const p of run) {
        const b = blocks[blocks.length - 1], last = b?.[b.length - 1];
        if (last && last.reversed === p.reversed && p.seg.i === last.seg.i + (p.reversed ? -1 : 1)) b.push(p); else blocks.push([p]);
      }
      const lo = (b: SvPiece[]) => Math.min(...b.map(p => p.seg.i)), hi = (b: SvPiece[]) => Math.max(...b.map(p => p.seg.i));
      if (blocks.every((b, k) => !k || lo(b) > hi(blocks[k - 1]))) {
        const kept = new Map<number, 'inv' | 'keep'>();
        for (const b of blocks) for (const p of b) kept.set(p.seg.id, p.reversed && b.length >= 1 ? 'inv' : 'keep');
        let k = 0;
        while (k < gap.length) {
          const what = kept.get(gap[k].id) ?? 'del';
          let e = k;
          if (what === 'inv') { const b = blocks.find(bl => bl.some(p => p.seg.id === gap[k].id))!; e = k + b.length - 1; }
          else while (e + 1 < gap.length && (kept.get(gap[e + 1].id) ?? 'del') === what) e++;
          if (what !== 'keep') ops.push(`${range({ start: gap[k].start, end: gap[e].end })}${what}`);
          k = e + 1;
        }
        return;
      }
    }
    // a forward repeat of the segments just before: tandem duplication
    if (!gap.length && run.every(p => p.seg.chrom === own && !p.reversed) && run.every((p, k) => !k || p.seg.i === run[k - 1].seg.i + 1) && run[run.length - 1].seg.i === left.seg.i) {
      ops.push(`${range({ start: run[0].seg.start, end: run[run.length - 1].seg.end })}dup`);
      return;
    }
    if (!gap.length) { ops.push(`${left.seg.end}_${left.seg.end + 1}ins${content(run)}`); return; }
    ops.push(`${range({ start: gap[0].start, end: gap[gap.length - 1].end, qter: !right })}delins${content(run)}`);
  };
  for (let k = 0; k < bb.length; k++) {
    const left = bb[k].p, right = bb[k + 1]?.p ?? null;
    const run = pieces.slice(bb[k].k + 1, right ? bb[k + 1].k : pieces.length);
    gapOps(left, right, run);
  }
  if (!ops.length) return null;
  return `${acc}:g.${ops.length === 1 ? ops[0] : `[${ops.join(';')}]`}`;
}

const bandAt = (bands: Cytoband[] | undefined, pos: number): string | null => bands?.find(b => pos >= b.start && pos < b.end)?.name ?? null;
const centromere = (bands: Cytoband[] | undefined): number | null => {
  const acen = bands?.filter(b => b.stain === 'acen') ?? [];
  return acen.length ? acen[0].end : null;
};

/** ISCN of the rebuilt chromosomes: a short form for a single del / dup / inv, else the detailed form. */
function iscnOf(rec: SvReconstruction, bandsOf: (chrom: string) => Cytoband[] | undefined, build: string): { band: string | null; seq: string[] } {
  const ders = rec.ders.filter(d => d.rearranged);
  const chroms = [...new Set(ders.flatMap(d => d.pieces.map(p => p.seg.chrom)))];
  const many = chroms.length > 1;
  const seq: string[] = [];
  const bandForms: string[] = [];
  let allBands = true;
  for (const d of ders) {
    const spans = d.pieces.map(spanOf);
    // merge pieces contiguous on the reference
    const merged: Span[] = [];
    for (const s of spans) {
      const m = merged[merged.length - 1];
      if (m && m.chrom === s.chrom && !m.reversed && !s.reversed && m.end === s.start) { m.end = s.end; m.qter = s.qter; }
      else if (m && m.chrom === s.chrom && m.reversed && s.reversed && s.end === m.start) { m.start = s.start; m.pter = s.pter; }
      else merged.push({ ...s });
    }
    seq.push(merged.map(s => {
      const cen = centromere(bandsOf(s.chrom));
      const a = s.pter ? 'pter' : comma(s.start + 1), b = s.qter ? 'qter' : comma(s.end);
      const mid = cen != null && cen > (s.pter ? 0 : s.start) && (s.qter || cen < s.end) ? '_cen' : '';
      return `chr${num(s.chrom)}:g.${a}${mid}_${b}${s.reversed ? 'inv' : ''}`;
    }).join('::'));
    const bandName = (s: Span, which: 'from' | 'to'): string => {
      const atStart = (which === 'from') !== s.reversed;
      if (atStart && s.pter) return `${many ? num(s.chrom) : ''}pter`;
      if (!atStart && s.qter) return `${many ? num(s.chrom) : ''}qter`;
      const b = bandAt(bandsOf(s.chrom), atStart ? s.start : Math.max(s.start, s.end - 1));
      if (!b) allBands = false;
      return `${many ? num(s.chrom) : ''}${b ?? '?'}`;
    };
    bandForms.push(`der(${num(d.chrom)})(${merged.map(s => `${bandName(s, 'from')}->${bandName(s, 'to')}`).join('::')})`);
  }
  // a single deletion, duplication or inversion in one chromosome: the short form
  const hg = rec.ders.filter(d => d.rearranged);
  if (hg.length === 1 && !many) {
    const h = hgvsOf(hg[0], rec.segs, build);
    const m = h?.match(/:g\.(\d+|pter)_(\d+|qter)(del|dup|inv)$/);
    const bands = bandsOf(hg[0].chrom);
    if (m && bands) {
      const a = m[1] === 'pter' ? null : bandAt(bands, +m[1] - 1), b = m[2] === 'qter' ? null : bandAt(bands, +m[2] - 1);
      if ((a || m[1] === 'pter') && (b || m[2] === 'qter')) {
        const terminal = m[1] === 'pter' || m[2] === 'qter';
        return { band: `seq[${build}] ${m[3]}(${num(hg[0].chrom)})(${terminal ? (a ?? b) : `${a}${b}`})`, seq };
      }
    }
  }
  return { band: allBands && bandForms.length ? `seq[${build}] ${bandForms.join(';')}` : null, seq };
}

export function nomenclatureOf(rec: SvReconstruction, build: string, bandsOf: (chrom: string) => Cytoband[] | undefined = () => undefined): SvNomenclature {
  const hgvs = rec.ders.filter(d => d.rearranged).map(d => hgvsOf(d, rec.segs, build)).filter((x): x is string => !!x);
  const { band, seq } = iscnOf(rec, bandsOf, build);
  const errs = Object.values(rec.bpError ?? {});
  const lo = errs.length ? Math.min(...errs) : 0, hi = errs.length ? Math.max(...errs) : 0;
  return { hgvs, iscn: band, iscnSeq: seq, uncertainty: !errs.length ? '' : lo === hi ? `±${hi} bp` : `±${lo}–${hi} bp` };
}

// ======================== Reading a description ========================

export interface ParsedSv {
  chrom: string; junctions: SvJunction[]; window: { start: number; end: number }; notes: string[]; rec: SvReconstruction;
  /** the same derivative over another window (the gene browser widens it) */
  rebuild: (w: { start: number; end: number }) => SvReconstruction;
}

const CHROM_RE = /^(?:(NC_\d+)(?:\.\d+)?|chr([0-9XYM]+)|([0-9XYM]{1,2}))$/i;
function chromOf(ref: string): string | null {
  const m = ref.trim().match(CHROM_RE);
  if (!m) return null;
  if (m[1]) { const c = chromOfAccession(m[1].toUpperCase()); return c ? `chr${c}` : null; }
  return `chr${(m[2] ?? m[3]).toUpperCase()}`;
}

/**
 * The rearrangement a description implies, as junctions on one chromosome, rebuilt as the cartoon rebuilds picked
 * arcs. Reads HGVS genomic descriptions (`NC_…:g.` or `chrN:g.`: del, dup, inv, ins, delins, the allele form
 * `[…;…]`, pter / qter, uncertain positions `(a_b)` taken at their middle) and the sequence form of ISCN
 * (`chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::chr20:g.10,808,806_qter`, after an optional
 * `seq[GRCh38]`).
 */
export function parseSv(text: string, build: string, chromLen: (chrom: string) => number | null): ParsedSv {
  const notes: string[] = [];
  const raw = text.trim().replace(/^seq\[[^\]]*\]\s*/i, '').replace(/\s+/g, '').replace(/,/g, '');
  if (!raw) throw new Error('Type an HGVS or ISCN description.');
  if (/^[a-z]+\(\w+(;\w+)?\)\(/i.test(raw) || /->|→/.test(raw)) throw new Error('Band-level ISCN gives no positions: type the HGVS description, or the ISCN sequence form (chr20:g.pter_10,489,765::chr20:g.10,783,341_10,807,815inv::…).');
  const lenOf = (c: string) => { const l = chromLen(c); if (!l) throw new Error(`Unknown chromosome ${c} on ${build}.`); return l; };
  // a position: number, pter, qter, cen (dropped), or (a_b) taken at its middle
  const pos = (s: string, c: string): number => {
    if (/^pter$/i.test(s)) return 0;
    if (/^qter$/i.test(s)) return lenOf(c);
    const u = s.match(/^\((\d+|\?)_(\d+|\?)\)$/);
    if (u) {
      const a = u[1] === '?' ? NaN : +u[1], b = u[2] === '?' ? NaN : +u[2];
      const v = Number.isFinite(a) && Number.isFinite(b) ? Math.round((a + b) / 2) : Number.isFinite(a) ? a : b;
      if (!Number.isFinite(v)) throw new Error(`Position ${s} is unknown on both sides.`);
      notes.push(`${s} taken at ${comma(v)}`);
      return v;
    }
    if (!/^\d+$/.test(s)) throw new Error(`Cannot read the position "${s}".`);
    return +s;
  };
  /** `a_b` (1-based inclusive) → 0-based half-open, `cen` ignored. */
  const rng = (s: string, c: string): [number, number] => {
    const parts = s.split('_').filter(x => !/^cen$/i.test(x));
    // uncertain positions contain "_" inside parentheses: rejoin them
    const toks: string[] = [];
    for (const p of parts) { const last = toks[toks.length - 1]; if (last && last.startsWith('(') && !last.endsWith(')')) toks[toks.length - 1] = `${last}_${p}`; else toks.push(p); }
    if (toks.length === 1) { const x = pos(toks[0], c); return [/^pter$/i.test(toks[0]) ? 0 : x - 1, /^qter$/i.test(toks[0]) ? x : x]; }
    if (toks.length !== 2) throw new Error(`Cannot read the range "${s}".`);
    const a = pos(toks[0], c), b = pos(toks[1], c);
    return [/^pter$/i.test(toks[0]) ? 0 : a - 1, b];
  };
  // inserted content: x_y, x_yinv, NC_…:g.x_y(inv), a literal sequence (not drawn), N[12]
  const contentOf = (s: string, c: string): Span[] => {
    const items = s.startsWith('[') && s.endsWith(']') ? splitTop(s.slice(1, -1)) : [s];
    const out: Span[] = [];
    for (const it of items) {
      const m = it.match(/^(?:([^:]+):g\.)?(.+?)(inv)?$/);
      const ic = m?.[1] ? chromOf(m[1]) : c;
      if (!m || !ic) throw new Error(`Cannot read "${it}".`);
      if (/^[ACGTN]+$|^N\[\d+\]$/i.test(m[2]) && !m[3]) { notes.push(`inserted bases ${it} not drawn`); continue; }
      const [a, b] = rng(m[2], ic);
      out.push({ chrom: ic, start: a, end: b, reversed: !!m[3], pter: a === 0, qter: b === lenOf(ic) });
    }
    return out;
  };

  let pieces: Span[];
  let chrom: string;
  if (raw.includes('::')) {
    // ISCN sequence form: the derivative from pter to qter
    const parts = raw.split('::');
    pieces = parts.map(p => {
      const m = p.match(/^(?:der\(\w+\))?([^:]+):g\.(.+?)(inv)?$/);
      const c = m ? chromOf(m[1]) : null;
      if (!m || !c) throw new Error(`Cannot read "${p}".`);
      const [a, b] = rng(m[2], c);
      return { chrom: c, start: a, end: b, reversed: !!m[3], pter: a === 0, qter: b === lenOf(c) };
    });
    chrom = pieces[0].chrom;
  } else {
    const m = raw.match(/^([^:]+):g\.(.+)$/);
    if (!m) throw new Error('Expected NC_…:g.… or chrN:g.… (HGVS), or the ISCN sequence form chrN:g.…::chrN:g.…');
    const c = chromOf(m[1]);
    if (!c) throw new Error(`Unknown reference "${m[1]}".`);
    chrom = c;
    const body = m[2];
    const ops = body.startsWith('[') && body.endsWith(']') ? splitTop(body.slice(1, -1)) : [body];
    type Op = { a: number; b: number; kind: 'del' | 'dup' | 'inv' | 'ins' | 'delins'; content: Span[] };
    const parsed: Op[] = ops.map(o => {
      const mm = o.match(/^(.+?)(delins|del|dup|inv|ins)(.*)$/);
      if (!mm) throw new Error(`Cannot read "${o}": expected del, dup, inv, ins or delins (substitutions are not structural).`);
      const kind = mm[2] as Op['kind'];
      const [a, b] = rng(mm[1], c);
      if ((kind === 'ins' || kind === 'delins') && !mm[3]) throw new Error(`"${o}": the inserted sequence is missing.`);
      if (kind !== 'ins' && kind !== 'delins' && mm[3]) throw new Error(`Cannot read "${o}".`);
      return { a, b, kind, content: kind === 'ins' || kind === 'delins' ? contentOf(mm[3], c) : [] };
    }).sort((x, y) => x.a - y.a);
    const L = lenOf(c);
    pieces = [];
    let cur = 0;
    const ref = (s: number, e: number) => { if (e > s) pieces.push({ chrom: c, start: s, end: e, reversed: false, pter: s === 0, qter: e === L }); };
    for (const o of parsed) {
      if (o.a < cur) throw new Error('The changes overlap: describe each stretch of the reference once.');
      if (o.kind === 'ins') { ref(cur, o.a + 1); pieces.push(...o.content); cur = o.a + 1; continue; }
      ref(cur, o.a);
      if (o.kind === 'inv') pieces.push({ chrom: c, start: o.a, end: o.b, reversed: true, pter: o.a === 0, qter: o.b === L });
      else if (o.kind === 'dup') { ref(o.a, o.b); ref(o.a, o.b); }
      else if (o.kind === 'delins') pieces.push(...o.content);
      cur = o.b;
    }
    ref(cur, L);
  }
  // junctions between pieces not contiguous on the reference
  const junctions: SvJunction[] = [];
  /** the junction that joins piece k to the one before */
  const joins: (SvJunction | null)[] = pieces.map(() => null);
  for (let k = 0; k + 1 < pieces.length; k++) {
    const p = pieces[k], q = pieces[k + 1];
    if (p.chrom === q.chrom && !p.reversed && !q.reversed && p.end === q.start) continue;
    if (p.chrom === q.chrom && p.reversed && q.reversed && q.end === p.start) continue;
    const a = p.reversed ? { chrom: p.chrom, pos: p.start, side: 'head' as const } : { chrom: p.chrom, pos: p.end, side: 'tail' as const };
    const b = q.reversed ? { chrom: q.chrom, pos: q.end, side: 'tail' as const } : { chrom: q.chrom, pos: q.start, side: 'head' as const };
    const kind: SvJunction['kind'] = a.chrom !== b.chrom ? 'translocation' : a.side === b.side ? 'inversion' : (a.side === 'tail' ? b.pos > a.pos : a.pos > b.pos) ? 'deletion' : 'duplication';
    junctions.push({ id: `J${junctions.length + 1}`, a, b, reads: 1, kind, evidence: 'typed', precision: 0 });
    joins[k + 1] = junctions[junctions.length - 1];
  }
  if (!junctions.length) throw new Error('This description changes no structure.');
  if (junctions.some(j => j.kind === 'translocation')) throw new Error('Translocations are not drawn yet: the browser shows one chromosome.');
  const bps = junctions.flatMap(j => [j.a.pos, j.b.pos]);
  const L = lenOf(chrom);
  const lo = Math.min(...bps), hi = Math.max(...bps), pad = Math.max(5000, Math.round((hi - lo) * 0.1));
  const window = { start: Math.max(0, lo - pad), end: Math.min(L, hi + pad) };
  const rebuild = (w: { start: number; end: number }) => recFromPieces(pieces, joins, junctions, chrom, w);
  const rec = rebuild(window);
  return { chrom, junctions, window, notes, rec, rebuild };
}

/** The rebuilt chromosome of known pieces: segments cut wherever a piece starts or ends inside the window. */
function recFromPieces(pieces: Span[], joins: (SvJunction | null)[], junctions: SvJunction[], chrom: string, window: { start: number; end: number }): SvReconstruction {
  const cuts = [...new Set(pieces.flatMap(p => [p.start, p.end]))].filter(x => x > window.start && x < window.end).sort((x, y) => x - y);
  const bounds = [window.start, ...cuts, window.end];
  const segs: SvSegment[] = bounds.slice(0, -1).map((b, i) => ({ id: i, chrom, letter: i < 26 ? String.fromCharCode(65 + i) : `S${i + 1}`, start: b, end: bounds[i + 1], first: i === 0, last: i === bounds.length - 2, i }));
  const visits = new Map<number, number>();
  const derPieces: SvPiece[] = [];
  pieces.forEach((p, k) => {
    const a = Math.max(p.start, window.start), b = Math.min(p.end, window.end);
    const inside = segs.filter(x => x.start >= a && x.end <= b);
    (p.reversed ? [...inside].reverse() : inside).forEach((seg, n) => {
      const copy = visits.get(seg.id) ?? 0; visits.set(seg.id, copy + 1);
      derPieces.push({ seg, reversed: p.reversed, copy, via: !derPieces.length ? null : n === 0 && joins[k] ? joins[k]! : 'ref' });
    });
  });
  return {
    segs, byChr: { [chrom]: segs }, chrs: [chrom],
    ders: [{ chrom, name: `der(${num(chrom)})`, pieces: derPieces, ends: '', rearranged: true }],
    lost: segs.filter(x => !visits.has(x.id)), warnings: [], junctions, merged: [], choice: [], options: [], unplaced: [], bpError: {},
  };
}

/** Splits `a;[b;c];d` at the top-level semicolons. */
function splitTop(s: string): string[] {
  const out: string[] = [];
  let depth = 0, cur = '';
  for (const ch of s) {
    if (ch === '[') depth++;
    if (ch === ']') depth--;
    if (ch === ';' && depth === 0) { out.push(cur); cur = ''; } else cur += ch;
  }
  if (cur) out.push(cur);
  return out;
}
