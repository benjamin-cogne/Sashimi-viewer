/**
 * Aggregate ("group") view: several samples pooled into one track.
 *
 * - Coverage is the sum of the members' depth profiles (drawn relative to its own maximum).
 * - Junction reads are pooled (summed) over the members.
 * - For every annotated intron of the displayed model, the events that use its donor or its
 *   acceptor compete: the canonical junction (C), alternative 5′ / 3′ sites (n), pseudo-exons (an
 *   alternative-3′ arc into a short cryptic exon and an alternative-5′ arc out of it, paired and
 *   weighted by the mean of their two arcs), exon skipping (S) and intron retention (the unspliced
 *   reads through the two boundaries, weighted by their mean). Each event is labelled with its
 *   weight over the sum of every weight at that intron, so the labels of one intron add up to 100 %.
 * - An exon-skipping arc spans two introns and shows 2·S over the two totals, which is the
 *   rMATS value 2·S / (I₁ + I₂ + 2·S) when nothing else competes at those introns; the two inclusion
 *   junctions then show the inclusion level (I₁ + I₂) / (I₁ + I₂ + 2·S), the same on both sides and
 *   the complement of the skipping label, instead of their own per-intron shares.
 * - Tooltips also give the rMATS-style value of each event against the canonical form alone.
 */
import type { BoundarySpanning, CoverageRun, JunctionArc, StructuralEvidence } from './types';
import { classifyJunction, junctionKey, type TxModel } from './geometry';
import { mergeSvArcs, poolPairArcs } from '../../standalone/svmerge';

/** Longest cryptic exon that two facing alternative-site arcs are paired into a pseudo-exon event. */
export const PSEUDO_EXON_MAX_BP = 500;

/**
 * Sum of several run-length coverage profiles (0-based half-open runs, each sorted). The result is
 * contiguous from the first to the last covered base, zero-depth runs included: the coverage path
 * joins consecutive runs with a straight line, so a gap would draw as a ramp across the intron.
 */
export function sumCoverage(profiles: CoverageRun[][]): CoverageRun[] {
  const events: { pos: number; delta: number }[] = [];
  for (const runs of profiles) {
    for (const r of runs) {
      if (r.end <= r.start || !r.depth) continue;
      events.push({ pos: r.start, delta: r.depth });
      events.push({ pos: r.end, delta: -r.depth });
    }
  }
  if (!events.length) return [];
  events.sort((a, b) => a.pos - b.pos || a.delta - b.delta);
  const out: CoverageRun[] = [];
  let depth = 0, from = events[0].pos;
  for (let i = 0; i < events.length; i++) {
    const e = events[i];
    if (e.pos !== from) {
      const last = out[out.length - 1];
      if (last && last.end === from && last.depth === depth) last.end = e.pos; else out.push({ start: from, end: e.pos, depth });
      from = e.pos;
    }
    depth += e.delta;
  }
  return out;
}

/** Junction reads summed over samples, plus how many samples carry each junction. */
export function poolJunctions(samples: { junctions: JunctionArc[] }[]): { junctions: JunctionArc[]; samplesWith: Map<string, number> } {
  const pooled = new Map<string, JunctionArc>();
  const samplesWith = new Map<string, number>();
  for (const s of samples) {
    for (const j of s.junctions) {
      if (j.count <= 0) continue;
      const k = junctionKey(j);
      const p = pooled.get(k);
      if (p) p.count += j.count; else pooled.set(k, { start: j.start, end: j.end, count: j.count });
      samplesWith.set(k, (samplesWith.get(k) || 0) + 1);
    }
  }
  return { junctions: [...pooled.values()].sort((a, b) => a.start - b.start || a.end - b.end), samplesWith };
}

/** Unspliced reads through the exon–intron boundaries, summed over samples (positions counted by any member are kept). */
export function poolSpanning(samples: { spanning?: BoundarySpanning }[]): BoundarySpanning | undefined {
  const withData = samples.filter(s => s.spanning);
  if (!withData.length) return undefined;
  const out: BoundarySpanning = { intronStart: {}, intronEnd: {} };
  for (const s of withData) {
    for (const [k, v] of Object.entries(s.spanning!.intronStart)) out.intronStart[+k] = (out.intronStart[+k] ?? 0) + v;
    for (const [k, v] of Object.entries(s.spanning!.intronEnd)) out.intronEnd[+k] = (out.intronEnd[+k] ?? 0) + v;
  }
  return out;
}

/** Structural evidence summed over samples (DNA groups): counts added by span, clip position or target chromosome. */
export function poolStructural(samples: { structural?: StructuralEvidence }[]): StructuralEvidence | undefined {
  const withData = samples.map(s => s.structural).filter((x): x is StructuralEvidence => !!x);
  if (!withData.length) return undefined;
  const clips = new Map<string, StructuralEvidence['clips'][number]>(), elsewhere = new Map<string, StructuralEvidence['elsewhere'][number]>(), ins = new Map<number, StructuralEvidence['insertions'][number]>();
  for (const s of withData) {
    for (const c of s.clips) { const k = `${c.side}${c.pos}`; const p = clips.get(k); if (p) p.count += c.count; else clips.set(k, { ...c }); }
    for (const e of s.elsewhere) { const k = `${e.kind}${e.chrom}@${e.pos}`; const p = elsewhere.get(k); if (p) p.count += e.count; else elsewhere.set(k, { ...e }); }
    for (const x of s.insertions) { const p = ins.get(x.pos); if (p) { p.len = Math.round((p.len * p.count + x.len * x.count) / (p.count + x.count)); p.count += x.count; } else ins.set(x.pos, { ...x }); }
  }
  const medians = withData.map(s => s.insertMedian).filter((x): x is number => x != null).sort((a, b) => a - b);
  return {
    // events of the members pooled with the same breakpoint tolerance as within a sample (svmerge.ts)
    deletions: mergeSvArcs(withData.map(s => [...s.deletions, ...s.splits])), splits: [], duplications: mergeSvArcs(withData.map(s => s.duplications ?? [])), inversions: mergeSvArcs(withData.map(s => s.inversions ?? [])), discordant: poolPairArcs(withData.map(s => s.discordant)),
    insertions: [...ins.values()].sort((a, b) => a.pos - b.pos),
    clips: [...clips.values()].sort((a, b) => a.pos - b.pos), elsewhere: [...elsewhere.values()].sort((a, b) => a.pos - b.pos),
    insertMedian: medians.length ? medians[medians.length >> 1] : null, reads: withData.reduce((a, s) => a + s.reads, 0),
  };
}

export type AggClass = 'canonical' | 'exon_skipping' | 'alt5' | 'alt3' | 'pseudo_exon' | 'novel';

export const AGG_CLASS_LABEL: Record<AggClass, string> = {
  canonical: 'canonical',
  exon_skipping: 'exon skipping',
  alt5: 'alternative 5′ site',
  alt3: 'alternative 3′ site',
  pseudo_exon: 'pseudo-exon',
  novel: 'novel junction',
};

/** Share of one event of the reads competing at its intron (both introns for an exon-skipping arc). */
export interface AggShare {
  /** event weight / sum of the weights of every competing event */
  pct: number;
  /** sum of the weights (the denominator) */
  total: number;
  /** what the percentage is relative to, for the tooltip */
  note: string;
}

export interface AggEvent {
  key: string;
  /** pooled reads of this junction */
  count: number;
  cls: AggClass;
  /** the other arc of a pseudo-exon pair */
  partner?: string;
  /** weight of the event at its intron (mean of the two arcs for a pseudo-exon, own reads otherwise) */
  eventCount: number;
  /** the event's share (one entry; an exon-skipping arc's single share spans its two introns) */
  shares: AggShare[];
}

/** Intron retention of one annotated intron: unspliced reads through both boundaries against the canonical junction. */
export interface AggRetention {
  /** intron index and genomic bounds (0-based half-open) */
  intron: number; start: number; end: number;
  fromExon: number; toExon: number;
  /** unspliced reads through the genomic-left boundary (intron start) and the genomic-right one (intron end); which is the donor depends on the strand */
  rStart: number; rEnd: number;
  canonical: number;
  /** (rStart + rEnd) / 2 over every weight competing at the intron */
  pct: number;
  note: string;
}

export interface AggResult {
  events: Map<string, AggEvent>;
  /** one entry per annotated intron with spanning data (introns without any unspliced read are listed with pct 0) */
  retention: AggRetention[];
}

/** Percent text for an arc label: "90 %", "0.4 %". */
export function pctLabel(p: number): string {
  const v = p * 100;
  return `${v >= 9.95 ? v.toFixed(0) : v >= 0.95 ? v.toFixed(1) : v.toFixed(v > 0 ? 1 : 0)} %`;
}

/**
 * Shares of every junction against the displayed model: at each annotated intron, every competing
 * event's weight over the sum of the weights (see the module comment). Weights: canonical C, alternative
 * site n, pseudo-exon (A + B) / 2 on both arcs, exon skipping S (at each of the two introns it spans),
 * intron retention (R5 + R3) / 2. A skipping arc is labelled 2·S / (total₁ + total₂). Junctions
 * touching no annotated splice site get no share (drawn with their pooled reads).
 */
export function aggregateJunctions(junctions: JunctionArc[], tx: TxModel | null, spanning?: BoundarySpanning, maxPseudoExon = PSEUDO_EXON_MAX_BP): AggResult {
  const out = new Map<string, AggEvent>();
  const retention: AggRetention[] = [];
  for (const j of junctions) {
    const info = classifyJunction(j, tx);
    const cls: AggClass = info.cls === 'canonical' ? 'canonical' : info.cls === 'exon_skipping' ? 'exon_skipping'
      : info.cls === 'novel_donor' ? 'alt5' : info.cls === 'novel_acceptor' ? 'alt3' : 'novel';
    out.set(junctionKey(j), { key: junctionKey(j), count: j.count, cls, eventCount: j.count, shares: [] });
  }
  if (!tx || tx.exons.length < 2) return { events: out, retention };
  const nIntrons = tx.exons.length - 1;
  const intronStart = (k: number) => tx.exons[k].end, intronEnd = (k: number) => tx.exons[k + 1].start;
  const exonLabel = (k: number, k2: number) => { const a = tx.exons[k].rank, b = tx.exons[k2].rank; return `${Math.min(a, b)}→${Math.max(a, b)}`; };
  /** canonical junction reads of intron k (0 when absent) */
  const canonical = Array.from({ length: nIntrons }, (_, k) => junctions.find(j => j.start === intronStart(k) && j.end === intronEnd(k))?.count ?? 0);

  // ---- exon skipping junctions: both ends annotated with at least one exon between them ----
  const exonIndexAtEnd = new Map<number, number>(), exonIndexAtStart = new Map<number, number>();
  tx.exons.forEach((e, i) => { exonIndexAtEnd.set(e.end, i); exonIndexAtStart.set(e.start, i); });
  const skips: { j: JunctionArc; k1: number; k2: number }[] = [];
  for (const j of junctions) {
    const li = exonIndexAtEnd.get(j.start), ri = exonIndexAtStart.get(j.end);
    if (li == null || ri == null || ri <= li + 1) continue;
    skips.push({ j, k1: li, k2: ri - 1 });               // flanking introns: exon li → li+1 and exon ri−1 → ri
  }
  const isSkip = (j: JunctionArc) => skips.some(sk => sk.j === j);

  // ---- pass 1, per intron: the competing events and their weights ----
  interface Pool {
    k: number; iStart: number; iEnd: number; label: string;
    C: number;
    /** unspliced reads through the two boundaries (when the source counted both) */
    hasSpan: boolean; rStart: number; rEnd: number; R: number; sideNote: string;
    pairs: { a: JunctionArc; b: JunctionArc; n: number }[];
    /** competing junctions that are neither canonical, nor paired, nor skips: alternative sites and other single junctions */
    alts: JunctionArc[];
    skips: { j: JunctionArc; k1: number; k2: number }[];
    /** every weight summed: C + Σ alt + Σ (A + B)/2 + Σ S + R/2 */
    total: number;
  }
  const pools: (Pool | null)[] = [];
  for (let k = 0; k < nIntrons; k++) {
    const iStart = intronStart(k), iEnd = intronEnd(k);
    if (iEnd <= iStart) { pools.push(null); continue; }
    const rStart = spanning?.intronStart[iStart], rEnd = spanning?.intronEnd[iEnd];
    const hasSpan = rStart != null && rEnd != null;
    const R = hasSpan ? rStart + rEnd : 0;
    // the donor is the genomic-left boundary on the + strand and the genomic-right one on the − strand
    const plus = tx.strand > 0;
    // last exonic base of each site, 1-based: intron start → iStart, intron end → iEnd + 1
    const donorPos = plus ? iStart : iEnd + 1, acceptorPos = plus ? iEnd + 1 : iStart;
    const sideNote = hasSpan
      ? `${(plus ? rStart : rEnd).toLocaleString('en-US')} unspliced reads through the donor at ${donorPos.toLocaleString('en-US')} + ${(plus ? rEnd : rStart).toLocaleString('en-US')} through the acceptor at ${acceptorPos.toLocaleString('en-US')}`
      : '';
    const competing = junctions.filter(j => j.start === iStart || j.end === iEnd);
    // pseudo-exon pairing: left-anchored arc ending inside the intron + right-anchored arc starting
    // inside it, in that order, with a short cryptic exon between them
    const lefts = competing.filter(j => j.start === iStart && j.end < iEnd).sort((a, b) => b.count - a.count);
    const rights = competing.filter(j => j.end === iEnd && j.start > iStart);
    const paired = new Set<string>();
    const pairs: Pool['pairs'] = [];
    for (const a of lefts) {
      const kA = junctionKey(a);
      if (paired.has(kA)) continue;
      const candidates = rights.filter(b => !paired.has(junctionKey(b)) && b.start > a.end && b.start - a.end <= maxPseudoExon)
        .sort((x, y) => (x.start - a.end) - (y.start - a.end));
      if (!candidates.length) continue;
      const b = candidates[0];
      paired.add(kA); paired.add(junctionKey(b));
      pairs.push({ a, b, n: (a.count + b.count) / 2 });
    }
    const alts = competing.filter(j => !paired.has(junctionKey(j)) && !isSkip(j) && !(j.start === iStart && j.end === iEnd));
    const here = skips.filter(sk => sk.k1 === k || sk.k2 === k);
    const total = canonical[k] + alts.reduce((a, j) => a + j.count, 0) + pairs.reduce((a, p) => a + p.n, 0) + here.reduce((a, sk) => a + sk.j.count, 0) + R / 2;
    pools.push({ k, iStart, iEnd, label: exonLabel(k, k + 1), C: canonical[k], hasSpan, rStart: rStart ?? 0, rEnd: rEnd ?? 0, R, sideNote, pairs, alts, skips: here, total });
  }

  // ---- pass 2: every event's share of the reads competing at its intron ----
  const pctTxt = (num: number, den: number) => pctLabel(den > 0 ? num / den : 0);
  const num = (x: number) => (Number.isInteger(x) ? x.toLocaleString('en-US') : x.toFixed(1));
  for (const pool of pools) {
    if (!pool) continue;
    const { k, label, C, R, total } = pool;
    const totalTxt = `${num(total)} reads competing at intron ${label}`;
    if (pool.hasSpan) {
      retention.push({
        intron: k, start: pool.iStart, end: pool.iEnd, fromExon: tx.exons[k].rank, toExon: tx.exons[k + 1].rank, rStart: pool.rStart, rEnd: pool.rEnd, canonical: C,
        pct: total > 0 ? (R / 2) / total : 0,
        note: `intron retention ${label} = (${pool.sideNote}; positions = last exonic base, 1-based) / 2 / ${totalTxt}` +
          `\nvs the canonical junction alone (rMATS-style): (${R.toLocaleString('en-US')}) / (${R.toLocaleString('en-US')} + 2 × ${C.toLocaleString('en-US')}) = ${pctTxt(R, R + 2 * C)}`,
      });
    }
    if (!total) continue;
    // pseudo-exon: one value on both arcs, the mean of the two junctions over the intron's total
    for (const p of pool.pairs) {
      const sum = p.a.count + p.b.count;
      for (const [me, other] of [[p.a, p.b], [p.b, p.a]] as const) {
        const ev = out.get(junctionKey(me))!;
        ev.cls = 'pseudo_exon'; ev.partner = junctionKey(other); ev.eventCount = p.n;
        ev.shares.push({
          pct: p.n / total, total,
          note: `pseudo-exon = (${p.a.count.toLocaleString('en-US')} + ${p.b.count.toLocaleString('en-US')} reads of its two junctions) / 2 / ${totalTxt}` +
            `\nvs the canonical junction alone (rMATS-style): (${sum.toLocaleString('en-US')}) / (${sum.toLocaleString('en-US')} + 2 × ${C.toLocaleString('en-US')}) = ${pctTxt(sum, sum + 2 * C)}`,
        });
      }
    }
    // alternative 5′ / 3′ site, or another single junction using one site of the intron
    for (const j of pool.alts) {
      out.get(junctionKey(j))!.shares.push({
        pct: j.count / total, total,
        note: `= ${j.count.toLocaleString('en-US')} reads / ${totalTxt}` +
          `\nvs the canonical junction alone (rMATS-style): ${j.count.toLocaleString('en-US')} / (${j.count.toLocaleString('en-US')} + ${C.toLocaleString('en-US')}) = ${pctTxt(j.count, j.count + C)}`,
      });
    }
    // canonical: its own reads over the total; the competitors are listed so a value below 100 % is explained even when they are hidden or off-screen
    const cj = junctions.find(j => j.start === pool.iStart && j.end === pool.iEnd);
    if (cj) {
      const ev = out.get(junctionKey(cj))!;
      ev.eventCount = C;
      // an inclusion junction whose only competitor is one skipping arc, on both sides of the skipped exon, carries the
      // inclusion level of that event pooled over the two flanking introns (rMATS ψ): the same value on both inclusion
      // junctions, the complement of the skipping arc's label
      const pure = (q: Pool | null, sk: { j: JunctionArc }) => !!q && q.skips.length === 1 && q.skips[0].j === sk.j && !q.alts.length && !q.pairs.length && q.R === 0;
      const sk = pool.skips.length === 1 ? pool.skips[0] : null;
      const other = sk ? pools[sk.k1 === k ? sk.k2 : sk.k1] : null;
      if (sk && other && pure(pool, sk) && pure(other, sk)) {
        const S = sk.j.count, Cs = C + other.C, den = Cs + 2 * S;
        const skipped = sk.k2 - sk.k1 > 1 ? `exons ${exonLabel(sk.k1 + 1, sk.k2)}` : `exon ${tx.exons[sk.k1 + 1].rank}`;
        ev.shares.push({
          pct: den > 0 ? Cs / den : 0, total: den,
          note: `inclusion of ${skipped} = (${C.toLocaleString('en-US')} + ${other.C.toLocaleString('en-US')} reads of the two inclusion junctions) / (${C.toLocaleString('en-US')} + ${other.C.toLocaleString('en-US')} + 2 × ${S.toLocaleString('en-US')} skipping reads) = ${pctTxt(Cs, den)}` +
            `\nthe same value on both inclusion junctions (rMATS ψ), the complement of the skipping arc's ${pctTxt(2 * S, den)}` +
            `\nthis junction alone at intron ${label}: ${C.toLocaleString('en-US')} / ${totalTxt} = ${pctTxt(C, total)}`,
        });
        continue;
      }
      const inclusion = pool.skips.map(x => { const o = pools[x.k1 === k ? x.k2 : x.k1]; if (!o) return ''; const S = x.j.count, Cs = C + o.C; return `\ninclusion level of the skipping ${(x.j.start + 1).toLocaleString('en-US')}-${x.j.end.toLocaleString('en-US')}, pooled over its two introns (rMATS ψ): (${C.toLocaleString('en-US')} + ${o.C.toLocaleString('en-US')}) / (${C.toLocaleString('en-US')} + ${o.C.toLocaleString('en-US')} + 2 × ${S.toLocaleString('en-US')}) = ${pctTxt(Cs, Cs + 2 * S)}`; }).join('');
      const others = [
        ...(R > 0 ? [`intron retention (${pool.sideNote}; counted as their mean)`] : []),
        ...pool.pairs.map(p => `pseudo-exon ${p.a.end.toLocaleString('en-US')}-${p.b.start.toLocaleString('en-US')} (${p.a.count.toLocaleString('en-US')} + ${p.b.count.toLocaleString('en-US')} reads, counted as their mean)`),
        ...pool.alts.map(x => `${AGG_CLASS_LABEL[out.get(junctionKey(x))!.cls]} ${(x.start + 1).toLocaleString('en-US')}-${x.end.toLocaleString('en-US')} (${x.count.toLocaleString('en-US')} reads)`),
        ...pool.skips.map(sk => `exon skipping ${(sk.j.start + 1).toLocaleString('en-US')}-${sk.j.end.toLocaleString('en-US')} (${sk.j.count.toLocaleString('en-US')} reads)`),
      ];
      ev.shares.push({
        pct: C / total, total,
        note: `= ${C.toLocaleString('en-US')} canonical reads / ${totalTxt}` +
          (others.length ? `\nother events at this intron, shown or not: ${others.slice(0, 5).join('; ')}${others.length > 5 ? `; +${others.length - 5} more` : ''}` : '') + inclusion,
      });
    }
  }
  // exon skipping: one value for the arc, 2·S over the totals of the two introns it spans (the rMATS value when nothing else competes there)
  for (const sk of skips) {
    const p1 = pools[sk.k1], p2 = pools[sk.k2];
    if (!p1 || !p2) continue;
    const S = sk.j.count, den = p1.total + p2.total, I1 = canonical[sk.k1], I2 = canonical[sk.k2];
    const li = sk.k1, ri = sk.k2 + 1;
    if (den <= 0) continue;
    out.get(junctionKey(sk.j))!.shares.push({
      pct: (2 * S) / den, total: den,
      note: `skipping of exon${ri - li > 2 ? `s ${exonLabel(li + 1, ri - 1)}` : ` ${tx.exons[li + 1].rank}`} = 2 × ${S.toLocaleString('en-US')} skipping reads / (${num(p1.total)} reads competing at intron ${p1.label} + ${num(p2.total)} at intron ${p2.label})` +
        `\nshare at intron ${p1.label}: ${pctTxt(S, p1.total)} · at intron ${p2.label}: ${pctTxt(S, p2.total)}` +
        `\nvs the inclusion junctions alone (rMATS-style): 2 × ${S.toLocaleString('en-US')} / (${I1.toLocaleString('en-US')} + ${I2.toLocaleString('en-US')} + 2 × ${S.toLocaleString('en-US')}) = ${pctTxt(2 * S, I1 + I2 + 2 * S)}`,
    });
  }
  return { events: out, retention };
}
