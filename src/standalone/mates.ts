/**
 * Which read of a set is whose mate: the one rule behind the pile-up's joined pairs and the fragments of phasing and
 * consensus groups.
 *
 * A read that carries a mate key (`mk`, set by the columnar decoder on two reads its stream links to each other) is
 * joined to the read carrying the same key and to no other. When that read is not in the set (outside the window,
 * filtered, sampled out), the keyed read stays alone: its mate is known, and it is not here.
 *
 * Every other read is joined by position: a read whose mate starts where this one says, pointing back at this one's
 * start, with the other pair bit (64/128) and no mate on another chromosome. Several reads can meet all of that when
 * fragments share their starts, and position alone then joins the first it finds — about one pair in 80,000 on a
 * paired RNA-seq file, and 37 % on a fixture built of look-alike pairs, measured. So the candidates are taken in
 * three passes over the whole set:
 *  1. the one with the same name, which makes a BAM's or CRAM's pairing exact;
 *  2. then the one whose template length is this read's negated: the two ends of one fragment carry opposite TLENs
 *     whatever the aligner's convention, and look-alikes of another length do not (the reads of an exported page
 *     are numbered, not named, and rely on this);
 *  3. then position alone.
 * Passes rather than one greedy walk ranking each read's candidates: a read whose own mate is missing would
 * otherwise take, by TLEN or position, a read whose same-named mate comes later in the set, leaving that one alone.
 *
 * Returns, per read, the index of its mate in `reads` or -1.
 */
import type { AlignedRead } from '../components/sashimi/types';

export function pairMates(reads: AlignedRead[]): Int32Array {
  const mate = new Int32Array(reads.length).fill(-1);
  const byKey = new Map<string, number>();
  const byStart = new Map<number, number[]>();
  reads.forEach((r, i) => {
    if (r.mk != null) {
      const j = byKey.get(r.mk);
      if (j == null) byKey.set(r.mk, i);
      else if (mate[j] < 0) { mate[i] = j; mate[j] = i; }
      return;
    }
    const l = byStart.get(r.s);
    if (l) l.push(i); else byStart.set(r.s, [i]);
  });
  /** 0 same name, 1 opposite template length, 2 position only */
  const rank = (r: AlignedRead, m: AlignedRead) => (m.n === r.n ? 0 : (m.tl ?? 0) === -(r.tl ?? 0) ? 1 : 2);
  for (let pass = 0; pass < 3; pass++) {
    reads.forEach((r, i) => {
      if (mate[i] >= 0 || r.mk != null || r.mp == null || r.mc) return;
      for (const j of byStart.get(r.mp) ?? []) {
        const m = reads[j];
        if (j === i || mate[j] >= 0 || m.mp !== r.s || m.mc || (m.f & 192) === (r.f & 192) || rank(r, m) > pass) continue;
        mate[i] = j; mate[j] = i;
        break;
      }
    });
  }
  return mate;
}
