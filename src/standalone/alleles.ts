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
 * as for coverage), each mismatching base whatever its quality (those under MIN_BQ flagged for the base-quality
 * check), insertions and deletions
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

/**
 * Supporting reads a site needs; the base quality under which an alternate base is flagged. The scan counts every
 * alternate base in the allele fraction (a low-quality one is flagged, the BQ check says how much of the allele they
 * are); callSites, on the reads of the reads track, still leaves them out.
 */
export const MIN_ALT = 3, MIN_BQ = 20;
/**
 * Quality evidence of each call, for the variants track: of the reads carrying an alternate allele, how many are on the
 * + strand, have a mapping quality under LOW_MQ, place it within END_BP of their alignment's ends, and (SNVs) how many
 * of them carried it with a base quality under MIN_BQ (counted in the allele fraction all the same). Kept per (position, allele)
 * as one number, four 13-bit counters (saturating), in a Map: only positions where an alternate base was seen take
 * room. The reference side of each comparison is counted over the blocks, like the depth: + strand, low mapping
 * quality and near-end depth.
 */
export const LOW_MQ = 20, END_BP = 10;
const QF = 8192, Q_FWD = 1, Q_LMQ = QF, Q_END = QF * QF, Q_LBQ = QF * QF * QF;
const qField = (v: number, unit: number) => Math.floor(v / unit) % QF;
/**
 * The flags of the first two alternate calls seen at a (position, base), 16 bits per position and base, in chunks:
 * most mismatches of a noisy long-read library are isolated errors, and a Map entry for each would take far more room
 * than the counts themselves. A site needs MIN_ALT (3) calls anyway: a (position, base) gets its entry at its third
 * call, the first two calls' flags folded in. Bits 0–3: the first call's + strand, low MAPQ, near-end, low base
 * quality; 4–7: the second's; 8: one call seen; 9: two.
 */
const FIRST_SET = 256, SECOND_SET = 512, F_FWD = 1, F_LMQ = 2, F_END = 4, F_LBQ = 8;
const flagUnits = (fl: number): number[] => [...(fl & F_FWD ? [Q_FWD] : []), ...(fl & F_LMQ ? [Q_LMQ] : []), ...(fl & F_END ? [Q_END] : []), ...(fl & F_LBQ ? [Q_LBQ] : [])];
class FirstFlags {
  private chunks = new Map<number, Uint16Array>();
  get bytes(): number { return this.chunks.size * 8192; }
  get(pos: number): number { return this.chunks.get(pos >> 12)?.[pos & 4095] ?? 0; }
  set(pos: number, v: number): void { let a = this.chunks.get(pos >> 12); if (!a) this.chunks.set(pos >> 12, a = new Uint16Array(4096)); a[pos & 4095] = v; }
}
const qBump = (m: Map<number, number>, key: number, units: number[]) => {
  if (!units.length && !m.has(key)) return;
  let v = m.get(key) ?? 0;
  for (const u of units) if (qField(v, u) < QF - 1) v += u;
  m.set(key, v);
};
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
  /** mismatching bases A, C, G, T, N, whatever their quality; any other code in `other` (position × 256 + code) */
  base = [new Hist(), new Hist(), new Hist(), new Hist(), new Hist()];
  other = new Map<number, number>();
  ins = new Map<number, number>();
  del = new Map<number, number>();
  /** quality evidence of the alternate calls: SNVs by position × 8 + base index (5: other codes), from their second call; indels by their key */
  snvQ = new Map<number, number>();
  first = [new FirstFlags(), new FirstFlags(), new FirstFlags(), new FirstFlags(), new FirstFlags()];
  insQ = new Map<number, number>();
  delQ = new Map<number, number>();
  /** blocks of + strand reads, of reads with a mapping quality under LOW_MQ, and the END_BP at each end of every alignment */
  FS = new Hist(); FE = new Hist();
  LS = new Hist(); LE = new Hist();
  NS = new Hist(); NE = new Hist();

  get bytes(): number {
    let b = this.S.bytes + this.E.bytes + this.RS.bytes + this.RE.bytes + this.FS.bytes + this.FE.bytes + this.LS.bytes + this.LE.bytes + this.NS.bytes + this.NE.bytes;
    for (const h of this.base) b += h.bytes;
    for (const f of this.first) b += f.bytes;
    return b + (this.other.size + this.ins.size + this.del.size) * 32 + (this.snvQ.size + this.insQ.size + this.delQ.size) * 48;
  }

  /**
   * One read. `codes[0..seqLen)` are its bases as character codes (seqLen 0 when the record has no sequence:
   * no mismatch is looked for, as encodeRead did), `quals` its base qualities (null: all taken as 30).
   */
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null, reverse = false, mapq = 60): void {
    // the alignment's end first: calls within END_BP of either end are flagged
    let aEnd = start;
    for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) aEnd += ops[k] >>> 4; }
    const fwd = !reverse, lowMq = mapq < LOW_MQ;
    const units = (p: number): number[] => { const u: number[] = []; if (fwd) u.push(Q_FWD); if (lowMq) u.push(Q_LMQ); if (p - start < END_BP || aEnd - p <= END_BP) u.push(Q_END); return u; };
    if (aEnd - start <= 2 * END_BP) { this.NS.add(start); this.NE.add(aEnd); }
    else { this.NS.add(start); this.NE.add(start + END_BP); this.NS.add(aEnd - END_BP); this.NE.add(aEnd); }
    let pos = start, q = 0;
    for (let k = 0; k < ops.length; k++) {
      const v = ops[k], len = v >>> 4, op = v & 15;
      if (op === 0 || op === 7 || op === 8) {
        this.S.add(pos); this.E.add(pos + len);
        if (fwd) { this.FS.add(pos); this.FE.add(pos + len); }
        if (lowMq) { this.LS.add(pos); this.LE.add(pos + len); }
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
            // every alternate base counts in the allele fraction, whatever its quality; a low one is flagged (the BQ check)
            const bi = BASE_INDEX[b], p = pos + i, qkey = p * 8 + (bi ?? 5);
            const u = units(p);
            if (qual < MIN_BQ) u.push(Q_LBQ);
            if (bi === undefined) { const key = p * 256 + b; this.other.set(key, (this.other.get(key) ?? 0) + 1); qBump(this.snvQ, qkey, u); continue; }
            this.base[bi].add(p);
            if (this.snvQ.has(qkey)) { qBump(this.snvQ, qkey, u); continue; }
            const f = this.first[bi].get(p), bits = (u.includes(Q_FWD) ? F_FWD : 0) | (u.includes(Q_LMQ) ? F_LMQ : 0) | (u.includes(Q_END) ? F_END : 0) | (u.includes(Q_LBQ) ? F_LBQ : 0);
            if (!(f & FIRST_SET)) { this.first[bi].set(p, FIRST_SET | bits); continue; }
            if (!(f & SECOND_SET)) { this.first[bi].set(p, f | SECOND_SET | (bits << 4)); continue; }
            // the third call: the entry starts with the first two calls' flags
            this.snvQ.set(qkey, 0);
            for (const fl of [f & 15, (f >> 4) & 15]) qBump(this.snvQ, qkey, flagUnits(fl));
            qBump(this.snvQ, qkey, u);
          }
        }
        pos += len; q += len;
      } else if (op === 1) { const key = pos * KEY + Math.min(len, KEY - 1); this.ins.set(key, (this.ins.get(key) ?? 0) + 1); qBump(this.insQ, key, units(pos)); q += len; }
      else if (op === 2) { const key = pos * KEY + Math.min(len, KEY - 1); this.del.set(key, (this.del.get(key) ?? 0) + 1); qBump(this.delQ, key, units(pos)); pos += len; }
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
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null, unique: boolean, reverse = false, mapq = 60): void {
    this.all.add(start, ops, codes, seqLen, quals, ref, reverse, mapq);
    if (!unique) this.multi.add(start, ops, codes, seqLen, quals, ref, reverse, mapq);
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
  /** a block quantity (depth, + strand depth, …) at every position of the window */
  const dense = (s: (c: AlleleCounts) => Hist, e: (c: AlleleCounts) => Hist): Int32Array => {
    const out = new Int32Array(w);
    let base = 0;
    for (const [c, sign] of parts) { s(c).addInto(out, start, sign); e(c).addInto(out, start, -sign); base += sign * (s(c).below(start) - e(c).below(start)); }
    for (let i = 0, d = base; i < w; i++) { d += out[i]; out[i] = d; }
    return out;
  };
  const depth = dense(c => c.S, c => c.E);
  let reads = 0;
  for (const [c, sign] of parts) reads += sign * (c.RS.below(end) - c.RE.below(start + 1));
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
  if (sites.length) addSiteQuality(sites, parts, dense, depth, start, ref);
  sites.sort((x, y) => x.pos - y.pos || x.kind.localeCompare(y.kind) || (x.alt < y.alt ? -1 : x.alt > y.alt ? 1 : 0));
  return { sites, reads };
}

/**
 * The quality evidence of each site (VariantSite.q): shares of the alternate reads on the + strand, with a low
 * mapping quality, with the call near an alignment end, and (SNVs) of the alternate bases under MIN_BQ; each
 * against the same share among the other reads over the position (the reference side), so an amplicon library
 * whose reads all run one way, or a repeat where every read maps poorly, is not taken for an artefact. For indels,
 * the length of the reference homopolymer at the site.
 */
function addSiteQuality(sites: VariantSite[], parts: [AlleleCounts, number][], dense: (s: (c: AlleleCounts) => Hist, e: (c: AlleleCounts) => Hist) => Int32Array,
  depth: Int32Array, start: number, ref: RefWindow | null): void {
  const fwdD = dense(c => c.FS, c => c.FE), lmqD = dense(c => c.LS, c => c.LE), endD = dense(c => c.NS, c => c.NE);
  const qOf = (pick: (c: AlleleCounts) => Map<number, number>, key: number) => { let v = 0; for (const [c, sign] of parts) v += sign * (pick(c).get(key) ?? 0); return v; };
  /** an SNV's evidence: each layer's entry, or the one or two calls it still holds as flags (a site's calls can be split between layers) */
  const snvOf = (pos: number, bi: number | undefined) => {
    let v = 0;
    for (const [c, sign] of parts) {
      const e = c.snvQ.get(pos * 8 + (bi ?? 5));
      if (e != null) { v += sign * e; continue; }
      if (bi === undefined) continue;
      const f = c.first[bi].get(pos);
      for (const [set, fl] of [[FIRST_SET, f & 15], [SECOND_SET, (f >> 4) & 15]]) if (f & set) v += sign * flagUnits(fl).reduce((a, x) => a + x, 0);
    }
    return v;
  };
  const share = (a: number, b: number) => (b > 0 ? Math.max(0, Math.min(1, a / b)) : null);
  for (const s of sites) {
    const i = s.pos - start, n = s.alt_count;
    let v: number;
    if (s.kind === 'snv') v = snvOf(s.pos, BASE_INDEX[s.alt.charCodeAt(0)]);
    else { const key = s.pos * KEY + Math.min(s.length, KEY - 1); v = s.kind === 'ins' ? qOf(c => c.insQ, key) : qOf(c => c.delQ, key); }
    const f = qField(v, Q_FWD), l = qField(v, Q_LMQ), e = qField(v, Q_END), lb = qField(v, Q_LBQ);
    // the other reads over the position: its blocks, less the alternate reads when they are among them (SNV, insertion)
    const inBlocks = s.kind !== 'del', bd = depth[i] - (inBlocks ? n : 0);
    s.q = {
      fwd: n ? f / n : 0, fwdRef: share(fwdD[i] - (inBlocks ? f : 0), bd),
      lowMq: n ? l / n : 0, lowMqRef: share(lmqD[i] - (inBlocks ? l : 0), bd) ?? 0,
      end: n ? e / n : 0, endRef: share(endD[i] - (inBlocks ? e : 0), bd), nRef: Math.max(0, bd),
      ...(s.kind === 'snv' ? { lowBq: n ? Math.min(1, lb / n) : 0 } : { hp: homopolymer(ref, s.pos, s.kind === 'ins') }),
    };
  }
}

/** Length of the reference homopolymer at an indel: the run through `pos` (deletion), or the longer of the runs on either side of it (insertion). */
function homopolymer(ref: RefWindow | null, pos: number, between: boolean): number {
  if (!ref) return 0;
  const at = (p: number) => (p - ref.start >= 0 && p - ref.start < ref.codes.length ? ref.codes[p - ref.start] : -1);
  const through = (p: number) => {
    const b = at(p);
    if (b < 0) return 0;
    let l = p, r = p;
    while (at(l - 1) === b) l--;
    while (at(r + 1) === b) r++;
    return r - l + 1;
  };
  return between ? Math.max(through(pos - 1), through(pos)) : through(pos);
}
