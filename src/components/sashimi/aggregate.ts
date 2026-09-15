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
 * - Exon skipping is measured rMATS-style against its two inclusion junctions pooled:
 *   2·S / (I₁ + I₂ + 2·S), one value for the skip and one for both inclusion arcs.
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

/** Share of one event: at one annotated intron, or of one exon-skipping event against its two inclusion junctions. */
export interface AggShare {
  /** event reads / reads of every competing event */
  pct: number;
  /** reads of every competing event (the denominator) */
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
  /** reads counted for the event (mean of the two arcs for a pseudo-exon, mean of the two inclusion junctions for a canonical arc flanking a skip) */
  eventCount: number;
  /** one share per intron the event competes at; an exon-skipping arc carries a single share */
  shares: AggShare[];
}

/** Percent text for an arc label: "90 %", "0.4 %". */
export function pctLabel(p: number): string {
  const v = p * 100;
  return `${v >= 9.95 ? v.toFixed(0) : v >= 0.95 ? v.toFixed(1) : v.toFixed(v > 0 ? 1 : 0)} %`;
}

/**
 * Shares of every junction against the displayed model.
 *
 * - Exon skipping (exon a → exon b, b > a + 1): one value, rMATS-style, from the two inclusion
 *   junctions pooled: skipping = 2·S / (I₁ + I₂ + 2·S), where I₁ and I₂ are the canonical
 *   junctions a → a+1 and b−1 → b and S the skipping reads. Both inclusion arcs then show the
 *   same inclusion level, (I₁ + I₂) / (I₁ + I₂ + 2·S), since a canonical arc flanking a skip is
 *   weighted by the mean of the two inclusion junctions, I/2, wherever it competes.
 * - Pseudo-exon (alternative-3′ arc A into a cryptic exon, alternative-5′ arc B out of it): the same
 *   rule against the canonical junction C that excludes it, inclusion = (A + B) / (A + B + 2·C),
 *   one value on both arcs.
 * - Alternative 5′ / 3′ site (one junction n): n / (n + C) against the canonical junction.
 * - The canonical arc shows a multi-way share at its intron: its weight (I/2 when it flanks a
 *   skip, its own reads otherwise) over every competitor (pseudo-exons as the mean of their two
 *   arcs, alternative sites, skips using one of its sites), so it reflects all of them at once.
 * - Junctions touching no annotated splice site get no share (drawn with their pooled reads).
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
  const nIntrons = tx.exons.length - 1;
  const intronStart = (k: number) => tx.exons[k].end, intronEnd = (k: number) => tx.exons[k + 1].start;
  const exonLabel = (k: number, k2: number) => { const a = tx.exons[k].rank, b = tx.exons[k2].rank; return `${Math.min(a, b)}→${Math.max(a, b)}`; };
  /** canonical junction reads of intron k (0 when absent) */
  const canonical = Array.from({ length: nIntrons }, (_, k) => junctions.find(j => j.start === intronStart(k) && j.end === intronEnd(k))?.count ?? 0);

  // ---- exon skipping: one share per skip, from the two inclusion junctions pooled ----
  const exonIndexAtEnd = new Map<number, number>(), exonIndexAtStart = new Map<number, number>();
  tx.exons.forEach((e, i) => { exonIndexAtEnd.set(e.end, i); exonIndexAtStart.set(e.start, i); });
  const skips: { j: JunctionArc; k1: number; k2: number; inclusion: number }[] = [];
  for (const j of junctions) {
    const li = exonIndexAtEnd.get(j.start), ri = exonIndexAtStart.get(j.end);
    if (li == null || ri == null || ri <= li + 1) continue;
    const k1 = li, k2 = ri - 1;                       // flanking introns: exon li → li+1 and exon ri−1 → ri
    const inclusion = canonical[k1] + canonical[k2];
    skips.push({ j, k1, k2, inclusion });
    const total = inclusion + 2 * j.count;
    if (total > 0) out.get(junctionKey(j))!.shares.push({
      pct: (2 * j.count) / total, total,
      note: `skipping of exon${ri - li > 2 ? `s ${exonLabel(li + 1, ri - 1)}` : ` ${tx.exons[li + 1].rank}`} = 2 × ${j.count.toLocaleString()} skipping reads / (${canonical[k1].toLocaleString()} + ${canonical[k2].toLocaleString()} inclusion reads at introns ${exonLabel(k1, k1 + 1)} and ${exonLabel(k2, k2 + 1)} + 2 × ${j.count.toLocaleString()})`,
    });
  }

  // ---- per intron: canonical vs alternative sites, pseudo-exons and skips using one of its sites ----
  for (let k = 0; k < nIntrons; k++) {
    const iStart = intronStart(k), iEnd = intronEnd(k);
    if (iEnd <= iStart) continue;
    const label = exonLabel(k, k + 1);
    const competing = junctions.filter(j => j.start === iStart || j.end === iEnd);
    if (!competing.length) continue;
    const touching = skips.filter(sk => sk.k1 === k || sk.k2 === k);
    // a canonical arc flanking one or more skips weighs the mean of the two inclusion junctions
    const canonicalWeight = touching.length ? touching.reduce((a, sk) => a + sk.inclusion / 2, 0) / touching.length : canonical[k];
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
    const isSkip = (j: JunctionArc) => skips.some(sk => sk.j === j);
    const isCanonical = (j: JunctionArc) => j.start === iStart && j.end === iEnd;
    const singles = competing.filter(j => !paired.has(junctionKey(j)));
    // multi-way total at the intron: canonical (weighted), pseudo-exons (mean of both arcs), alternative sites, skips
    const total = pairs.reduce((a, p) => a + p.n, 0) + singles.reduce((a, j) => a + (isCanonical(j) ? canonicalWeight : j.count), 0);
    if (!total) continue;
    const C = canonical[k];
    const cTxt = `${C.toLocaleString()} canonical reads at intron ${label}`;
    // pseudo-exon: the two inclusion arcs against the canonical junction that excludes the cryptic exon
    for (const p of pairs) {
      const sum = p.a.count + p.b.count, denom = sum + 2 * C;
      for (const [me, other] of [[p.a, p.b], [p.b, p.a]] as const) {
        const ev = out.get(junctionKey(me))!;
        ev.cls = 'pseudo_exon'; ev.partner = junctionKey(other); ev.eventCount = p.n;
        if (denom > 0) ev.shares.push({ pct: sum / denom, total: denom, note: `pseudo-exon inclusion = (${p.a.count.toLocaleString()} + ${p.b.count.toLocaleString()} reads of its two junctions) / (${sum.toLocaleString()} + 2 × ${cTxt})` });
      }
    }
    for (const j of singles) {
      if (isSkip(j)) continue;                          // a skip carries its own single share
      const ev = out.get(junctionKey(j))!;
      if (isCanonical(j)) {
        ev.eventCount = Math.round(canonicalWeight);
        // the competitors, so a canonical arc below 100 % is explained even when they are hidden (below the threshold) or off-screen
        const others = [
          ...pairs.map(p => `pseudo-exon ${p.a.end.toLocaleString()}-${p.b.start.toLocaleString()} (${p.a.count.toLocaleString()} + ${p.b.count.toLocaleString()} reads)`),
          ...singles.filter(x => x !== j).map(x => {
            const cls = out.get(junctionKey(x))!.cls;
            return `${AGG_CLASS_LABEL[cls]} ${(x.start + 1).toLocaleString()}-${x.end.toLocaleString()} (${x.count.toLocaleString()} reads)`;
          }),
        ];
        ev.shares.push({
          pct: canonicalWeight / total, total: Math.round(total),
          note: `of the ${Math.round(total).toLocaleString()} reads competing at intron ${label}${touching.length ? ' (canonical weighted by the mean of the two inclusion junctions of the skip)' : ''}` +
            (others.length ? `\nother events at this intron, shown or not: ${others.slice(0, 5).join('; ')}${others.length > 5 ? `; +${others.length - 5} more` : ''}` : ''),
        });
      } else {
        // alternative 5′ / 3′ site (or another single junction using one site of the intron) against the canonical junction
        const denom = j.count + C;
        if (denom > 0) ev.shares.push({ pct: j.count / denom, total: denom, note: `= ${j.count.toLocaleString()} reads / (${j.count.toLocaleString()} + ${cTxt})` });
      }
    }
  }
  return out;
}
