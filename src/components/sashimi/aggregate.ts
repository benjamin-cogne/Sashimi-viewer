/**
 * Aggregate ("group") view: several samples pooled into one track.
 *
 * - Coverage is the sum of the members' depth profiles (drawn relative to its own maximum).
 * - Junction reads are pooled (summed) over the members.
 * - For every annotated intron of the displayed model, the junctions that use its donor or its
 *   acceptor compete: the canonical intron, alternative 5′ / 3′ sites, exon skipping and
 *   pseudo-exons (an alternative-3′ arc into a short cryptic exon followed by an alternative-5′
 *   arc out of it, paired into one event). Each event's share of that intron's reads is the
 *   percentage shown on the arc instead of a read count.
 */
import type { CoverageRun, JunctionArc } from './types';
import { classifyJunction, junctionKey, type TxModel } from './geometry';

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

export type AggClass = 'canonical' | 'exon_skipping' | 'alt5' | 'alt3' | 'pseudo_exon' | 'novel';

export const AGG_CLASS_LABEL: Record<AggClass, string> = {
  canonical: 'canonical',
  exon_skipping: 'exon skipping',
  alt5: 'alternative 5′ site',
  alt3: 'alternative 3′ site',
  pseudo_exon: 'pseudo-exon',
  novel: 'novel junction',
};

/** Share of one event at one annotated intron. */
export interface AggShare {
  /** index of the intron in the model (between exons[k] and exons[k+1], genomic order) */
  intron: number;
  /** exon ranks flanking the intron, in transcript order */
  fromExon: number; toExon: number;
  /** event reads / reads of every event competing at this intron */
  pct: number;
  /** reads of every competing event at this intron */
  total: number;
}

export interface AggEvent {
  key: string;
  /** pooled reads of this junction */
  count: number;
  cls: AggClass;
  /** the other arc of a pseudo-exon pair */
  partner?: string;
  /** reads counted for the event (mean of the two arcs for a pseudo-exon) */
  eventCount: number;
  /** one share per intron the event competes at (an exon-skipping arc competes at two) */
  shares: AggShare[];
}

/** Percent text for an arc label: "90 %", "0.4 %". */
export function pctLabel(p: number): string {
  const v = p * 100;
  return `${v >= 9.95 ? v.toFixed(0) : v >= 0.95 ? v.toFixed(1) : v.toFixed(v > 0 ? 1 : 0)} %`;
}

/**
 * Per-intron shares of every junction against the displayed model. Junctions that touch no
 * annotated splice site get no share (they are drawn with their pooled read count).
 */
export function aggregateJunctions(junctions: JunctionArc[], tx: TxModel | null, maxPseudoExon = PSEUDO_EXON_MAX_BP): Map<string, AggEvent> {
  const out = new Map<string, AggEvent>();
  for (const j of junctions) {
    const info = classifyJunction(j, tx);
    const cls: AggClass = info.cls === 'canonical' ? 'canonical' : info.cls === 'exon_skipping' ? 'exon_skipping'
      : info.cls === 'novel_donor' ? 'alt5' : info.cls === 'novel_acceptor' ? 'alt3' : 'novel';
    out.set(junctionKey(j), { key: junctionKey(j), count: j.count, cls, eventCount: j.count, shares: [] });
  }
  if (!tx || tx.exons.length < 2) return out;

  for (let k = 0; k + 1 < tx.exons.length; k++) {
    const iStart = tx.exons[k].end, iEnd = tx.exons[k + 1].start;
    if (iEnd <= iStart) continue;
    const rL = tx.exons[k].rank, rR = tx.exons[k + 1].rank;
    const fromExon = Math.min(rL, rR), toExon = Math.max(rL, rR);
    // competitors: every junction using the intron's left or right boundary
    const competing = junctions.filter(j => j.start === iStart || j.end === iEnd);
    if (!competing.length) continue;
    // pseudo-exon pairing: left-anchored arc ending inside the intron + right-anchored arc starting
    // inside it, in that order, with a short cryptic exon between them
    const lefts = competing.filter(j => j.start === iStart && j.end < iEnd).sort((a, b) => b.count - a.count);
    const rights = competing.filter(j => j.end === iEnd && j.start > iStart);
    const paired = new Set<string>();
    const pairs: { a: JunctionArc; b: JunctionArc; n: number }[] = [];
    for (const a of lefts) {
      const kA = junctionKey(a);
      if (paired.has(kA)) continue;
      const candidates = rights.filter(b => !paired.has(junctionKey(b)) && b.start > a.end && b.start - a.end <= maxPseudoExon)
        .sort((x, y) => (x.start - a.end) - (y.start - a.end));
      if (!candidates.length) continue;
      const b = candidates[0];
      paired.add(kA); paired.add(junctionKey(b));
      pairs.push({ a, b, n: Math.round((a.count + b.count) / 2) });
    }
    const events: { keys: string[]; n: number }[] = [
      ...pairs.map(p => ({ keys: [junctionKey(p.a), junctionKey(p.b)], n: p.n })),
      ...competing.filter(j => !paired.has(junctionKey(j))).map(j => ({ keys: [junctionKey(j)], n: j.count })),
    ];
    const total = events.reduce((s, e) => s + e.n, 0);
    if (!total) continue;
    for (const p of pairs) {
      for (const [me, other] of [[p.a, p.b], [p.b, p.a]] as const) {
        const ev = out.get(junctionKey(me))!;
        ev.cls = 'pseudo_exon'; ev.partner = junctionKey(other); ev.eventCount = p.n;
      }
    }
    for (const e of events) {
      for (const key of e.keys) {
        const ev = out.get(key)!;
        ev.shares.push({ intron: k, fromExon, toExon, pct: e.n / total, total });
      }
    }
  }
  return out;
}
