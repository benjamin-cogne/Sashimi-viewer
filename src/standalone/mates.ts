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
 *  1. the one with the same name, which makes a BAM's or CRAM's pairing exact. It is looked up among the reads of
 *     that name, as IGV does (a map keyed by read name), not among the reads starting at the mate's position: the
 *     same read is found, but a pile-up of D reads on one start no longer costs D² comparisons;
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
  const fits = (r: AlignedRead, i: number, m: AlignedRead, j: number) =>
    j !== i && mate[j] < 0 && m.s === r.mp && m.mp === r.s && !m.mc && (m.f & 192) !== (r.f & 192);
  const alone = (r: AlignedRead, i: number) => mate[i] < 0 && r.mk == null && r.mp != null && !r.mc;
  // Every pass takes, for a read, the first read in the set's order that fits, among those starting at its mate's
  // position. Where more than PILE start there (a pile-up: rRNA, amplicons, duplicates), going through them all for
  // each read cost D² for D reads; they are then indexed once, on the key the pass asks for, and each key's list is
  // walked past its reads taken already. The indices are those of that start alone, built when first asked for.
  const index = new Map<string, Map<string, { l: number[]; p: number }>>();
  const indexed = (pass: number, s: number, key: (m: AlignedRead) => string | null) => {
    const id = `${pass}:${s}`;
    let x = index.get(id);
    if (!x) {
      index.set(id, x = new Map());
      for (const j of byStart.get(s)!) {
        const k = mate[j] < 0 ? key(reads[j]) : null;
        if (k == null) continue;
        const e = x.get(k);
        if (e) e.l.push(j); else x.set(k, { l: [j], p: 0 });
      }
    }
    return x;
  };
  const firstFree = (x: Map<string, { l: number[]; p: number }>, k: string) => {
    const e = x.get(k);
    if (!e) return -1;
    while (e.p < e.l.length && mate[e.l[e.p]] >= 0) e.p++;
    return e.p < e.l.length ? e.l[e.p] : -1;
  };
  // 1. the read of the same name
  reads.forEach((r, i) => {
    if (!alone(r, i)) return;
    const at = byStart.get(r.mp!);
    if (!at) return;
    if (at.length > PILE) {
      // same name, and the fit tested on the read found (the start is the list's; the rest is not in the key)
      const e = indexed(0, r.mp!, m => m.n).get(r.n);
      if (!e) return;
      for (let q = e.p; q < e.l.length; q++) {
        const j = e.l[q];
        if (!fits(r, i, reads[j], j)) continue;
        mate[i] = j; mate[j] = i;
        break;
      }
      return;
    }
    for (const j of at) {
      const m = reads[j];
      if (m.n !== r.n || !fits(r, i, m, j)) continue;
      mate[i] = j; mate[j] = i;
      break;
    }
  });
  // 2 and 3. A read still alone has no free read of its name that fits (none was free when its turn came, and none
  // has been freed since); it takes the first free read that fits with the opposite template length, then the first
  // that fits. In a pile the reads that fit are those of this start whose mate starts at the read's own start, of
  // either other pair bit and no mate on another chromosome: the key.
  for (let pass = 2; pass <= 3; pass++) {
    const key = (m: AlignedRead) => (m.mp == null || m.mc ? null
      : pass === 2 ? `${m.mp},${m.f & 192},${m.tl ?? 0}` : `${m.mp},${m.f & 192}`);
    reads.forEach((r, i) => {
      if (!alone(r, i)) return;
      const at = byStart.get(r.mp!);
      if (!at) return;
      if (at.length > PILE) {
        const x = indexed(pass, r.mp!, key);
        let best = -1;
        for (const bits of [0, 64, 128, 192]) {
          if (bits === (r.f & 192)) continue;
          const j = firstFree(x, pass === 2 ? `${r.s},${bits},${-(r.tl ?? 0)}` : `${r.s},${bits}`);
          if (j >= 0 && (best < 0 || j < best)) best = j;
        }
        if (best >= 0) { mate[i] = best; mate[best] = i; }
        return;
      }
      for (const j of at) {
        const m = reads[j];
        if (!fits(r, i, m, j) || (pass === 2 && (m.tl ?? 0) !== -(r.tl ?? 0))) continue;
        mate[i] = j; mate[j] = i;
        break;
      }
    });
  }
  return mate;
}

/** Reads on one start past which a pass indexes them rather than going through them for each read. */
const PILE = 16;

/**
 * A 32-bit hash of a read name (FNV-1a, then MurmurHash3's finalizer, so that its low bits are as mixed as the rest):
 * the two mates of a pair share their name, so a sample keeping the reads whose hash is a multiple of a power of two
 * keeps or drops both, as `samtools view --subsample` does and as IGV keeps a sampled read's mate.
 */
export function nameHash(name: string): number {
  let h = 0x811c9dc5;
  for (let i = 0; i < name.length; i++) h = Math.imul(h ^ name.charCodeAt(i), 0x01000193);
  h ^= h >>> 16; h = Math.imul(h, 0x85ebca6b);
  h ^= h >>> 13; h = Math.imul(h, 0xc2b2ae35);
  return (h ^ (h >>> 16)) >>> 0;
}
