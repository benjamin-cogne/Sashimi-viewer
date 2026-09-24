/**
 * Long reads place a junction a few bases off where the bases next to it carry errors (ONT especially, where minimap2
 * aligns through a noisy exon end): one junction then shows as a fan of arcs, the true one and neighbours a few bases
 * away with a handful of reads each, and no two reads of one isoform share a splice pattern.
 *
 * A junction within JUNCTION_SNAP_BP at both ends of one seen at least JUNCTION_SNAP_RATIO times as often is taken
 * for that one. IsoQuant corrects to the annotation within 6 bp on ONT data (its `delta`, 4 on PacBio) and FLAIR
 * within 15; here there is no annotation, only the reads. The ratio is strict on purpose: a variant creating a
 * cryptic donor or acceptor a few bases from the canonical one gives a junction at one exact place, used by a real
 * share of the reads (a whole haplotype's, often), while alignment jitter spreads a few per cent over many offsets.
 * A cryptic site used by at least 1 / JUNCTION_SNAP_RATIO as many reads as its canonical neighbour stays its own
 * junction; one used less is merged into it (a known limit, said in the arc's tooltip).
 */
export const JUNCTION_SNAP_BP = 6, JUNCTION_SNAP_RATIO = 20;
/** Junction key: start and length packed in one safe integer (introns up to 8.4 Mb, starts up to 2^30). */
export const JUNCTION_KEY = 8_388_608;
export const junctionKey = (s: number, e: number) => s * JUNCTION_KEY + Math.min(e - s, JUNCTION_KEY - 1);

/** For each junction to be taken for a more common one (by key), that one's [start, end). */
export function snapJunctions(list: { start: number; end: number; count: number }[]): Map<number, [number, number]> {
  const byCount = [...list].sort((a, b) => b.count - a.count || a.start - b.start || a.end - b.end);
  const byStart = [...list].sort((a, b) => a.start - b.start);
  const snap = new Map<number, [number, number]>(), done = new Set<number>();
  for (const c of byCount) {
    const ck = junctionKey(c.start, c.end);
    if (done.has(ck)) continue;
    done.add(ck);
    let lo = 0, hi = byStart.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (byStart[mid].start < c.start - JUNCTION_SNAP_BP) lo = mid + 1; else hi = mid; }
    for (let i = lo; i < byStart.length && byStart[i].start <= c.start + JUNCTION_SNAP_BP; i++) {
      const o = byStart[i], ok = junctionKey(o.start, o.end);
      if (done.has(ok) || Math.abs(o.end - c.end) > JUNCTION_SNAP_BP || o.count * JUNCTION_SNAP_RATIO > c.count) continue;
      done.add(ok); snap.set(ok, [c.start, c.end]);
    }
  }
  return snap;
}
