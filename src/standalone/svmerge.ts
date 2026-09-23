/**
 * Structural arcs with nearby breakpoints merged into one event.
 *
 * The reads of one structural variant rarely agree to the base on its breakpoints: long reads place them with a few
 * to a few tens of bases of scatter (more in repeats and at microhomologies), and one event shows through several
 * kinds of evidence at once, a deletion inside the CIGAR in some reads, a split alignment (hard- or soft-clipped
 * supplementary part, SA tag) in others, clipped reads placed by realignment in others still; an inversion shows
 * two junctions (+ then −, − then +) a few bases apart. Drawn as they come, one event made several arcs of a few
 * reads each, none of which reached the support threshold.
 *
 * Arcs of one family (deletion: CIGAR D and deletion-type splits; duplication; inversion) are merged when both
 * their starts and their ends lie within `svMergeTolerance` of each other: 5 % of the event's length, at least
 * SV_MERGE_MIN_BP and at most SV_MERGE_MAX_BP. The ceiling follows the long-read SV callers' ONT settings (cuteSV
 * recommends --max_cluster_bias_DEL 100 for ONT; Sniffles2 merges within 150 bp by default); the relative term
 * and the floor keep two distinct small events (two 60 bp deletions 40 bp apart) apart. The tolerance of the
 * smaller of the two events applies, and both ends must match, so nested events stay separate.
 */
import type { SvArc } from '../components/sashimi/types';

export const SV_MERGE_MIN_BP = 20, SV_MERGE_MAX_BP = 100, SV_MERGE_FRACTION = 0.05;
export const svMergeTolerance = (len: number): number => Math.min(SV_MERGE_MAX_BP, Math.max(SV_MERGE_MIN_BP, Math.round(Math.abs(len) * SV_MERGE_FRACTION)));

/** Two arcs are one event: both ends within the smaller event's tolerance. */
export function sameSvEvent(a: { start: number; end: number }, b: { start: number; end: number }): boolean {
  const tol = Math.min(svMergeTolerance(a.end - a.start), svMergeTolerance(b.end - b.start));
  return Math.abs(a.start - b.start) <= tol && Math.abs(a.end - b.end) <= tol;
}

/** The event of `list` that `j` belongs to (the closest when several could), or undefined. */
export function findSvEvent<T extends { start: number; end: number }>(list: T[], j: { start: number; end: number }): T | undefined {
  let best: T | undefined, bestD = Infinity;
  for (const x of list) {
    if (!sameSvEvent(x, j)) continue;
    const d = Math.abs(x.start - j.start) + Math.abs(x.end - j.end);
    if (d < bestD) { best = x; bestD = d; }
  }
  return best;
}

/** One arc before merging: its evidence units by source, and the reads behind them (by name; unnamed ones counted apart). */
export interface SvMember { start: number; end: number; count: number; src: Record<string, number>; names: Set<string>; anon: number; /** a member that is an event already */ spread?: [number, number, number, number] }

/** Weighted median of values (weights ≥ 0). */
function weightedMedian(v: { x: number; w: number }[]): number {
  const s = [...v].sort((a, b) => a.x - b.x), total = s.reduce((t, e) => t + e.w, 0);
  let acc = 0;
  for (const e of s) { acc += e.w; if (acc * 2 >= total) return e.x; }
  return s[s.length - 1].x;
}

/**
 * Clusters the members of one family into events, in start order, each member joining the nearest event it is within
 * tolerance of on both ends. An event's breakpoints are the read-weighted medians of its members';
 * its count the distinct reads behind them (a read showing both junctions of an inversion counts once).
 */
export function clusterSv(members: SvMember[]): { event: SvArc; members: SvMember[] }[] {
  const sorted = [...members].sort((a, b) => a.start - b.start || a.end - b.end);
  // complete linkage: a member joins a group only when it lies within tolerance of every member already in it (of the
  // group's extreme starts and ends), so a group never drifts wider than the tolerance by chaining
  const groups: { members: SvMember[]; s: number; e: number; w: number; s0: number; s1: number; e0: number; e1: number }[] = [];
  for (const m of sorted) {
    let best: (typeof groups)[number] | null = null, bestD = Infinity;
    for (let k = groups.length - 1; k >= 0; k--) {
      const g = groups[k];
      if (m.start - g.s0 > SV_MERGE_MAX_BP) break;   // sorted by start: older groups are out of reach
      const tol = Math.min(svMergeTolerance(m.end - m.start), svMergeTolerance(g.e - g.s));
      if (Math.max(Math.abs(m.start - g.s0), Math.abs(m.start - g.s1)) > tol || Math.max(Math.abs(m.end - g.e0), Math.abs(m.end - g.e1)) > tol) continue;
      const d = Math.abs(g.s - m.start) + Math.abs(g.e - m.end);
      if (d < bestD) { best = g; bestD = d; }
    }
    const w = Math.max(1e-9, m.count);
    if (best) {
      best.members.push(m); best.s = (best.s * best.w + m.start * w) / (best.w + w); best.e = (best.e * best.w + m.end * w) / (best.w + w); best.w += w;
      best.s0 = Math.min(best.s0, m.start); best.s1 = Math.max(best.s1, m.start); best.e0 = Math.min(best.e0, m.end); best.e1 = Math.max(best.e1, m.end);
    } else groups.push({ members: [m], s: m.start, e: m.end, w, s0: m.start, s1: m.start, e0: m.end, e1: m.end });
  }
  return groups.map(g => {
    const names = new Set<string>(); let anon = 0;
    const sources: Record<string, number> = {};
    for (const m of g.members) { for (const n of m.names) names.add(n); anon += m.anon; for (const [k, v] of Object.entries(m.src)) sources[k] = (sources[k] ?? 0) + v; }
    const start = weightedMedian(g.members.map(m => ({ x: m.start, w: m.count }))), end = weightedMedian(g.members.map(m => ({ x: m.end, w: m.count })));
    const event: SvArc = { start, end: Math.max(end, start + 1), count: names.size + anon, sources };
    const merged = g.members.reduce((n, m) => n + (m.spread ? 2 : 1), 0) > 1;
    if (merged) {
      event.merged = g.members.length;
      const sp = g.members.map(m => m.spread ?? [m.start, m.start, m.end, m.end]);
      event.spread = [Math.min(...sp.map(x => x[0])), Math.max(...sp.map(x => x[1])), Math.min(...sp.map(x => x[2])), Math.max(...sp.map(x => x[3]))];
    }
    return { event, members: g.members };
  });
}

/** Events of several lists (samples pooled in a group) merged the same way; counts add up. */
export function mergeSvArcs(lists: SvArc[][]): SvArc[] {
  const members: SvMember[] = lists.flat().map(a => ({ start: a.start, end: a.end, count: a.count, src: { ...(a.sources ?? {}) }, names: new Set<string>(), anon: a.count, ...(a.spread ? { spread: a.spread } : {}) }));
  return clusterSv(members).map(x => x.event).sort((a, b) => a.start - b.start || a.end - b.end);
}
