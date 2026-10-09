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
import type { SiteHaplotypes, VariantSite } from '../components/sashimi/types';
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

/**
 * Haplotags (HP, PS) of a phased file: each (phase set, haplotype) met in a stretch is numbered from 1 (0: untagged), at
 * most HAP_GROUPS of them; the allele counts below are kept per number.
 */
export const HAP_GROUPS = 32767;
const HAP_KEY = 32768, HAP_LEN = 1024;
export class HapGroups {
  private index = new Map<string, number>();
  /** PS (null: tagged without one) and HP of each number */
  ps: (number | null)[] = [null];
  hp: number[] = [0];
  of(hp: number, ps: number | null): number {
    const key = `${ps ?? ''}:${hp}`;
    let g = this.index.get(key);
    if (g == null) {
      if (this.hp.length > HAP_GROUPS) return 0;
      g = this.hp.length; this.index.set(key, g); this.ps.push(ps); this.hp.push(hp);
    }
    return g;
  }
}
/** A growable list of 32-bit integers. */
class IntList {
  a = new Int32Array(64);
  n = 0;
  push(v: number): void { if (this.n === this.a.length) { const b = new Int32Array(this.a.length * 2); b.set(this.a); this.a = b; } this.a[this.n++] = v; }
  /** the values, sorted (cached until the next push) */
  private s: Int32Array | null = null;
  private sn = -1;
  sorted(): Int32Array { if (this.sn !== this.n) { this.s = this.a.slice(0, this.n).sort(); this.sn = this.n; } return this.s!; }
  /** values ≤ x */
  upTo(x: number): number { const v = this.sorted(); let l = 0, h = v.length; while (l < h) { const m = (l + h) >> 1; if (v[m] <= x) l = m + 1; else h = m; } return l; }
}
/** The groups of the first two calls at a (position, base): 16 bits each (see FirstFlags). */
class FirstGroups {
  private chunks = new Map<number, Uint32Array>();
  get bytes(): number { return this.chunks.size * 16384; }
  get(pos: number): number { return this.chunks.get(pos >> 12)?.[pos & 4095] ?? 0; }
  set(pos: number, v: number): void { let a = this.chunks.get(pos >> 12); if (!a) this.chunks.set(pos >> 12, a = new Uint32Array(4096)); a[pos & 4095] = v; }
}
/**
 * The counts of the tagged reads, by group: where each read starts and ends (a read spans a site when it starts at or
 * before it and ends after it), and the alternate calls of each group. SNV calls are kept the way their quality evidence
 * is: from a (position, base)'s third call, the first two calls' groups held in `first`; indels from their first call.
 */
class HapCounts {
  starts = new Map<number, IntList>();
  ends = new Map<number, IntList>();
  /** first and last base covered by each group's reads */
  span = new Map<number, [number, number]>();
  snv = new Map<number, number>();
  ins = new Map<number, number>();
  del = new Map<number, number>();
  first = [new FirstGroups(), new FirstGroups(), new FirstGroups(), new FirstGroups(), new FirstGroups()];
  get bytes(): number {
    let b = (this.snv.size + this.ins.size + this.del.size) * 32;
    for (const l of this.starts.values()) b += l.a.length * 8;
    for (const f of this.first) b += f.bytes;
    return b;
  }
  read(g: number, start: number, end: number): void {
    let s = this.starts.get(g), e = this.ends.get(g);
    if (!s || !e) { this.starts.set(g, s = new IntList()); this.ends.set(g, e = new IntList()); }
    s.push(start); e.push(end);
    const sp = this.span.get(g);
    if (!sp) this.span.set(g, [start, end]); else { if (start < sp[0]) sp[0] = start; if (end > sp[1]) sp[1] = end; }
  }
  /** the group's reads spanning `pos` */
  depth(g: number, pos: number): number { const s = this.starts.get(g), e = this.ends.get(g); return s && e ? s.upTo(pos) - e.upTo(pos) : 0; }
}
const bump = (m: Map<number, number>, k: number) => m.set(k, (m.get(k) ?? 0) + 1);
const indelHapKey = (pos: number, len: number, g: number) => (pos * HAP_LEN + Math.min(len, HAP_LEN - 1)) * HAP_KEY + g;

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
  /** the reads carrying a haplotag, by group (allocated with the first one) */
  hap: HapCounts | null = null;

  get bytes(): number {
    let b = this.S.bytes + this.E.bytes + this.RS.bytes + this.RE.bytes + this.FS.bytes + this.FE.bytes + this.LS.bytes + this.LE.bytes + this.NS.bytes + this.NE.bytes;
    for (const h of this.base) b += h.bytes;
    for (const f of this.first) b += f.bytes;
    if (this.hap) b += this.hap.bytes;
    return b + (this.other.size + this.ins.size + this.del.size) * 32 + (this.snvQ.size + this.insQ.size + this.delQ.size) * 48;
  }

  /**
   * One read. `codes[0..seqLen)` are its bases as character codes (seqLen 0 when the record has no sequence:
   * no mismatch is looked for, as encodeRead did), `quals` its base qualities (null: all taken as 30).
   */
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null, reverse = false, mapq = 60, g = 0): void {
    // the alignment's end first: calls within END_BP of either end are flagged
    let aEnd = start;
    for (let k = 0; k < ops.length; k++) { const op = ops[k] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) aEnd += ops[k] >>> 4; }
    // a tagged read (g: its haplotag's group, 0 untagged)
    const h = g ? (this.hap ??= new HapCounts()) : this.hap;
    if (g) h!.read(g, start, aEnd);
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
            if (this.snvQ.has(qkey)) { qBump(this.snvQ, qkey, u); if (g) bump(h!.snv, qkey * HAP_KEY + g); continue; }
            const f = this.first[bi].get(p), bits = (u.includes(Q_FWD) ? F_FWD : 0) | (u.includes(Q_LMQ) ? F_LMQ : 0) | (u.includes(Q_END) ? F_END : 0) | (u.includes(Q_LBQ) ? F_LBQ : 0);
            if (!(f & FIRST_SET)) { this.first[bi].set(p, FIRST_SET | bits); if (g) h!.first[bi].set(p, g); continue; }
            if (!(f & SECOND_SET)) { this.first[bi].set(p, f | SECOND_SET | (bits << 4)); if (g) h!.first[bi].set(p, h!.first[bi].get(p) | (g << 16)); continue; }
            // the third call: the entry starts with the first two calls' flags (and groups)
            this.snvQ.set(qkey, 0);
            for (const fl of [f & 15, (f >> 4) & 15]) qBump(this.snvQ, qkey, flagUnits(fl));
            qBump(this.snvQ, qkey, u);
            if (h) {
              const fg = h.first[bi].get(p);
              for (const x of [fg & 0xffff, fg >>> 16, g]) if (x) bump(h.snv, qkey * HAP_KEY + x);
            }
          }
        }
        pos += len; q += len;
      } else if (op === 1) { const key = pos * KEY + Math.min(len, KEY - 1); this.ins.set(key, (this.ins.get(key) ?? 0) + 1); qBump(this.insQ, key, units(pos)); if (g) bump(h!.ins, indelHapKey(pos, len, g)); q += len; }
      else if (op === 2) { const key = pos * KEY + Math.min(len, KEY - 1); this.del.set(key, (this.del.get(key) ?? 0) + 1); qBump(this.delQ, key, units(pos)); if (g) bump(h!.del, indelHapKey(pos, len, g)); pos += len; }
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
  add(start: number, ops: ArrayLike<number>, codes: Uint8Array, seqLen: number, quals: ArrayLike<number> | null, ref: RefWindow | null, unique: boolean, reverse = false, mapq = 60, g = 0): void {
    this.all.add(start, ops, codes, seqLen, quals, ref, reverse, mapq, g);
    if (!unique) this.multi.add(start, ops, codes, seqLen, quals, ref, reverse, mapq, g);
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
  /** the haplotags met (a phased file) */
  groups = new HapGroups();
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
  minVaf: number, minIndel: number, groups?: HapGroups): { sites: VariantSite[]; reads: number; haplotagged: number } {
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
  const haplotagged = groups ? addSiteHaplotypes(sites, parts, groups, start, end) : 0;
  return { sites, reads, haplotagged };
}

/**
 * The site's reads by haplotag (VariantSite.hap), when reads over the window carry one: per group, the reads spanning
 * the site and those with its allele; the phase set holding the most tagged reads over the site is the site's, its HP 1
 * and HP 2 are given apart, any other tagged read is pooled. Returns the tagged reads over [start, end).
 */
function addSiteHaplotypes(sites: VariantSite[], parts: [AlleleCounts, number][], groups: HapGroups, start: number, end: number): number {
  const hparts = parts.filter(([c]) => c.hap) as [AlleleCounts & { hap: HapCounts }, number][];
  if (!hparts.length) return 0;
  // the groups and their spans, over every part
  const spans = new Map<number, [number, number]>();
  for (const [c] of hparts) for (const [g, [a, b]] of c.hap.span) { const sp = spans.get(g); if (!sp) spans.set(g, [a, b]); else { sp[0] = Math.min(sp[0], a); sp[1] = Math.max(sp[1], b); } }
  let tagged = 0;
  for (const g of spans.keys()) for (const [c, sign] of hparts) { const s = c.hap.starts.get(g), e = c.hap.ends.get(g); if (s && e) tagged += sign * (s.upTo(end - 1) - e.upTo(start)); }
  if (tagged <= 0) return 0;
  const glist = [...spans.entries()].sort((x, y) => x[1][0] - y[1][0]);
  for (const s of sites) {
    const over = glist.filter(([, [a, b]]) => a <= s.pos && b > s.pos).map(([g]) => g);
    if (!over.length) continue;
    const depthOf = (g: number) => { let d = 0; for (const [c, sign] of hparts) d += sign * c.hap.depth(g, s.pos); return d; };
    const altOf = (g: number) => {
      let n = 0;
      for (const [c, sign] of hparts) {
        if (s.kind === 'snv') {
          const bi = BASE_INDEX[s.alt.charCodeAt(0)];
          if (bi === undefined) continue;
          const qkey = s.pos * 8 + bi;
          if (c.snvQ.has(qkey)) { n += sign * (c.hap.snv.get(qkey * HAP_KEY + g) ?? 0); continue; }
          // one or two calls at the position in this part: their groups
          const f = c.first[bi].get(s.pos), fg = c.hap.first[bi].get(s.pos);
          if (f & FIRST_SET && (fg & 0xffff) === g) n += sign;
          if (f & SECOND_SET && fg >>> 16 === g) n += sign;
        } else n += sign * ((s.kind === 'ins' ? c.hap.ins : c.hap.del).get(indelHapKey(s.pos, s.length, g)) ?? 0);
      }
      return n;
    };
    const hap = siteHaplotypesOf(over.map(g => ({ ps: groups.ps[g], hp: groups.hp[g], depth: depthOf(g), alt: altOf(g) })));
    if (hap) s.hap = hap;
  }
  return tagged;
}

/**
 * A site's counts by haplotag from its counts per (phase set, haplotype): the phase set with the most HP 1 and HP 2 reads
 * over the site is the site's (ties: the lower PS), its two haplotypes apart, every other tagged read pooled.
 */
export function siteHaplotypesOf(per: { ps: number | null; hp: number; depth: number; alt: number }[]): SiteHaplotypes | null {
  const over = per.filter(x => x.depth > 0);
  if (!over.length) return null;
  const bySet = new Map<number | null, number>();
  for (const x of over) if (x.hp === 1 || x.hp === 2) bySet.set(x.ps, (bySet.get(x.ps) ?? 0) + x.depth);
  const rank = (ps: number | null) => ps ?? Infinity;
  const best = [...bySet.entries()].sort((p, q) => q[1] - p[1] || rank(p[0]) - rank(q[0]))[0];
  const hap: SiteHaplotypes = { ps: best ? best[0] : null, alt: [0, 0], depth: [0, 0], otherAlt: 0, otherDepth: 0 };
  for (const x of over) {
    if (best && x.ps === best[0] && (x.hp === 1 || x.hp === 2)) { hap.alt[x.hp - 1] += x.alt; hap.depth[x.hp - 1] += x.depth; }
    else { hap.otherAlt += x.alt; hap.otherDepth += x.depth; }
  }
  return hap;
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
