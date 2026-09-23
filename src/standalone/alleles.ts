/**
 * Variant sites from allele counts streamed out of the records, and kept.
 *
 * The full variant scan used to decode every read of a window into an object (name, sequence,
 * qualities, blocks, mismatches found by comparing strings), then call sites from those objects,
 * one 100 kb tile at a time without handing the thread back. Measured in the browser on a 1 000×
 * whole-gene capture, one zoom-out that widened the window froze the page for 11.5 s: encodeRead
 * 11.2 s, garbage collection 7.7 s, callSites 1.9 s.
 *
 * Here each record's bases go straight into per-position counts: the depth (block starts and ends,
 * as for coverage), each mismatching base with a quality of MIN_BQ or more, insertions and deletions
 * by position and length, and where reads start and end (the reads a range holds). Sites are then
 * read back with exactly the rules of callSites (collapse.ts): at least MIN_ALT supporting reads, an
 * alternate fraction of at least Min VAF over the block depth (plus the deleting reads themselves for
 * a deletion). The counts stay with the sample, so moving or zooming back reads nothing again.
 *
 * Nothing here decodes: the caller hands over each record's start, packed CIGAR, base codes and
 * qualities, and the reference over the record's span.
 */
import type { VariantSite } from '../components/sashimi/types';
import { Hist } from './coverage';

/** callSites' fixed thresholds: supporting reads, and base quality of a counted mismatch. */
export const MIN_ALT = 3, MIN_BQ = 20;
const KEY = 8_388_608;          // position × KEY + length, for insertions and deletions
const N_CODE = 78;              // 'N'
const BASE_INDEX: Record<number, number> = { 65: 0, 67: 1, 71: 2, 84: 3, 78: 4 };   // A C G T N
const BASES = ['A', 'C', 'G', 'T', 'N'];

/** Reference bases (upper case character codes) of [start, start + codes.length). */
export interface RefWindow { start: number; codes: Uint8Array }
export function refWindow(start: number, seq: string | null): RefWindow | null {
  if (seq == null) return null;
  const codes = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) codes[i] = seq.charCodeAt(i);
  return { start, codes };
}

class AlleleCounts {
  /** aligned blocks: depth at x = starts ≤ x − ends ≤ x */
  S = new Hist();
  E = new Hist();
  /** reads: where they start, where their span ends — the reads over [a, b) are starts < b − ends ≤ a */
  RS = new Hist();
  RE = new Hist();
  /** mismatching bases A, C, G, T, N with quality ≥ MIN_BQ; any other code in `other` (position × 256 + code) */
  base = [new Hist(), new Hist(), new Hist(), new Hist(), new Hist()];
  other = new Map<number, number>();
  ins = new Map<number, number>();
  del = new Map<number, number>();

  get bytes(): number {
    let b = this.S.bytes + this.E.bytes + this.RS.bytes + this.RE.bytes;
    for (const h of this.base) b += h.bytes;
    return b + (this.other.size + this.ins.size + this.del.size) * 32;
  }

  /**
   * One read. `codes[0..seqLen)` are its bases as character codes (seqLen 0 when the record has no sequence:
   * no mismatch is looked for, as encodeRead did), `quals` its base qualities (null: all taken as 30).
   */
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null): void {
    let pos = start, q = 0;
    for (let k = 0; k < ops.length; k++) {
      const v = ops[k], len = v >>> 4, op = v & 15;
      if (op === 0 || op === 7 || op === 8) {
        this.S.add(pos); this.E.add(pos + len);
        if (seqLen && ref) {
          const rc = ref.codes, r0 = ref.start;
          for (let i = 0; i < len; i++) {
            const qi = q + i;
            if (qi >= seqLen) break;
            const ri = pos + i - r0;
            if (ri < 0 || ri >= rc.length) continue;
            const rb = rc[ri], b = codes[qi];
            if (b === rb || rb === N_CODE) continue;
            const qual = quals && quals.length > qi ? quals[qi] : 30;
            if (qual < MIN_BQ) continue;
            const bi = BASE_INDEX[b];
            if (bi !== undefined) this.base[bi].add(pos + i);
            else { const key = (pos + i) * 256 + b; this.other.set(key, (this.other.get(key) ?? 0) + 1); }
          }
        }
        pos += len; q += len;
      } else if (op === 1) { const key = pos * KEY + Math.min(len, KEY - 1); this.ins.set(key, (this.ins.get(key) ?? 0) + 1); q += len; }
      else if (op === 2) { const key = pos * KEY + Math.min(len, KEY - 1); this.del.set(key, (this.del.get(key) ?? 0) + 1); pos += len; }
      else if (op === 3) pos += len;
      else if (op === 4) q += len;
    }
    this.RS.add(start); this.RE.add(pos);
  }
}

/** The reads of one part of a scan: every read goes to `all`, a read that is not uniquely mapped to `multi` as well. */
export class AlleleLayer {
  all = new AlleleCounts();
  multi = new AlleleCounts();
  get bytes(): number { return this.all.bytes + this.multi.bytes; }
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null, unique: boolean): void {
    this.all.add(start, ops, codes, seqLen, quals, ref);
    if (!unique) this.multi.add(start, ops, codes, seqLen, quals, ref);
  }
}

/** Allele counts of one sample over one stretch of a chromosome: owned reads start in [ps, pe), the spill reaches in from the left (see CoverageState). */
export class AlleleState {
  owned = new AlleleLayer();
  spill = new AlleleLayer();
  ps = 0;
  pe = 0;
  lastUsed = 0;
  /** long reads (median aligned span above 1 kb), decided on the first reads counted */
  longReads: boolean | null = null;
  /** aligned spans of the first reads, until the decision */
  spans: number[] = [];
  constructor(readonly chrom: string) {}
  get empty(): boolean { return this.pe <= this.ps; }
  get bytes(): number { return this.owned.bytes + this.spill.bytes; }
  covers(start: number, end: number): boolean { return !this.empty && this.ps <= start && this.pe >= end; }
}

/**
 * The sites of [start, end), from the layers holding every read over it once, by callSites' rules:
 * SNV alt count n and block depth d at its position, n ≥ MIN_ALT and n / d ≥ minVaf; an insertion the
 * same with the depth at its position; a deletion with d + n (the deleting reads are not in the blocks).
 * `ref` gives the reference base of an SNV ('?' without one). Also the reads over [start, end).
 */
export function sitesFromCounts(layers: AlleleLayer[], uniqueOnly: boolean, start: number, end: number, ref: RefWindow | null,
  minVaf: number, minIndel: number): { sites: VariantSite[]; reads: number } {
  const parts: [AlleleCounts, number][] = [];
  for (const l of layers) { parts.push([l.all, 1]); if (uniqueOnly) parts.push([l.multi, -1]); }
  const w = Math.max(0, end - start);
  const depth = new Int32Array(w);
  let base = 0, reads = 0;
  for (const [c, sign] of parts) {
    c.S.addInto(depth, start, sign); c.E.addInto(depth, start, -sign);
    base += sign * (c.S.below(start) - c.E.below(start));
    reads += sign * (c.RS.below(end) - c.RE.below(start + 1));
  }
  for (let i = 0, d = base; i < w; i++) { d += depth[i]; depth[i] = d; }
  const sites: VariantSite[] = [];
  const refBase = (pos: number) => (ref && pos - ref.start >= 0 && pos - ref.start < ref.codes.length ? String.fromCharCode(ref.codes[pos - ref.start]) : '?');
  const snv = (pos: number, alt: string, n: number) => {
    const d = depth[pos - start];
    if (n >= MIN_ALT && d && n / d >= minVaf) sites.push({ pos, kind: 'snv', ref: refBase(pos), alt, length: 0, alt_count: n, depth: d, vaf: n / d });
  };
  const counts = new Int32Array(w);
  for (let b = 0; b < BASES.length; b++) {
    counts.fill(0);
    for (const [c, sign] of parts) c.base[b].addInto(counts, start, sign);
    for (let i = 0; i < w; i++) if (counts[i] >= MIN_ALT) snv(start + i, BASES[b], counts[i]);
  }
  const merged = (pick: (c: AlleleCounts) => Map<number, number>) => {
    const m = new Map<number, number>();
    for (const [c, sign] of parts) for (const [k, v] of pick(c)) m.set(k, (m.get(k) ?? 0) + sign * v);
    return m;
  };
  for (const [k, n] of merged(c => c.other)) {
    const pos = Math.floor(k / 256);
    if (pos >= start && pos < end) snv(pos, String.fromCharCode(k % 256), n);
  }
  for (const [k, n] of merged(c => c.ins)) {
    const pos = Math.floor(k / KEY), len = k % KEY;
    if (len < minIndel || pos < start || pos >= end) continue;
    const d = depth[pos - start];
    if (n >= MIN_ALT && d && n / d >= minVaf) sites.push({ pos, kind: 'ins', ref: '', alt: `+${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  for (const [k, n] of merged(c => c.del)) {
    const ds = Math.floor(k / KEY), len = k % KEY;
    if (len < minIndel || ds < start || ds >= end) continue;
    const d = depth[ds - start] + n;
    if (n >= MIN_ALT && d && n / d >= minVaf) sites.push({ pos: ds, kind: 'del', ref: '', alt: `-${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  sites.sort((x, y) => x.pos - y.pos || x.kind.localeCompare(y.kind) || (x.alt < y.alt ? -1 : x.alt > y.alt ? 1 : 0));
  return { sites, reads };
}
