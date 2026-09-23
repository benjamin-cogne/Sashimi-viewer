/**
 * CpG methylation from the base-modification tags of long reads (ONT: dorado; PacBio: primrose / jasmine), counted the
 * way ONT's modkit builds a CpG table (`modkit pileup --cpg --ref … --combine-strands`, its "traditional" preset):
 *
 * - **Sites from the reference, not from the reads**: every CpG of the reference; a read counts at the C of a CpG when
 *   it aligns there (+ strand reads) or at its G (− strand reads: the C of the other strand), both strands summed on the
 *   C, and only when the read carries the canonical base there (a SNP destroying the CpG is not a call).
 * - **Calls from MM / ML** (SAM tags specification): MM lists, per modification, how many bases of its canonical kind
 *   to skip between two listed ones, counted on the read as sequenced (for a read aligned on the − strand, from the
 *   end of SEQ, on its complement); ML gives each listed base a probability (N ↦ (N + 0.5) / 256). A `?` entry says
 *   nothing of the bases it skips; a `.` entry (or none) says they are unmodified. The legacy Mm / Ml names are read too.
 *   Records whose SEQ no longer matches the read the tags describe (hard clips without an MN tag, MN ≠ SEQ length) are
 *   left out, as the specification asks.
 * - **5mC against unmodified C, 5hmC set aside**: the 5hmC probability is shared out between the two others (their
 *   probabilities renormalised), then the call is the likelier state and its probability the call's confidence.
 * - **Low-confidence calls filtered**: calls below a threshold estimated from the data, the 10th percentile of the
 *   sample's call confidences (modkit's default filter percentile), are counted apart and left out of the fractions.
 *
 * Counts are kept per CpG, per haplotag (untagged, HP 1, HP 2) and per confidence bin, so the threshold can be applied
 * when a window is read back, from whatever has been counted.
 */

/** Confidence bins over [0.5, 1]: 32nds. */
export const CONF_BINS = 16;
const HAPS = 3;                         // 0 untagged (or HP > 2), 1, 2
/** Share of the sample's calls the threshold filters: its lowest-confidence 10 %, as modkit does by default. */
export const FILTER_PERCENTILE = 0.1;
const C = 67, G = 71;

/** Counts of one CpG: [hap][state 0 unmodified / 1 modified][confidence bin]. */
export type CpgCounts = Uint16Array;
const idx = (hap: number, state: number, bin: number) => (hap * 2 + state) * CONF_BINS + bin;
export const confBin = (conf: number) => Math.max(0, Math.min(CONF_BINS - 1, Math.floor((conf - 0.5) * 2 * CONF_BINS)));

/** CpG sites of a reference stretch: the positions of their C (0-based). */
export function cpgSites(start: number, seq: string): Int32Array<ArrayBuffer> {
  const out: number[] = [];
  for (let i = 0; i + 1 < seq.length; i++) {
    const a = seq.charCodeAt(i) & 0xdf, b = seq.charCodeAt(i + 1) & 0xdf;   // upper case
    if (a === C && b === G) out.push(start + i);
  }
  return Int32Array.from(out);
}

/** Scratch buffers reused from read to read. */
export interface ModScratch { pm: Float32Array; ph: Float32Array; listed: Uint8Array; cIdx: Int32Array }
export const newModScratch = (): ModScratch => ({ pm: new Float32Array(1024), ph: new Float32Array(1024), listed: new Uint8Array(1024), cIdx: new Int32Array(1024) });

/**
 * Parses the MM / ML tags of a read for C modifications on the sequenced strand into `sc`, indexed by SEQ position:
 * pm / ph the 5mC / 5hmC probabilities, listed[i] 1 when the base was listed. Returns whether unlisted C's are
 * implicitly unmodified (a `.` or bare 5mC entry), or null when the tags do not describe this SEQ.
 */
export function parseModTags(mm: string, ml: ArrayLike<number> | null, codes: Uint8Array, n: number, reverse: boolean, sc: ModScratch): boolean | null {
  if (sc.pm.length < n) { const m = Math.max(n, sc.pm.length * 2); sc.pm = new Float32Array(m); sc.ph = new Float32Array(m); sc.listed = new Uint8Array(m); sc.cIdx = new Int32Array(m); }
  sc.listed.fill(0, 0, n);
  // SEQ positions of the C's of the read as sequenced, in sequencing order: C's of SEQ forward, or G's of SEQ from its end
  let nc = 0;
  if (!reverse) { for (let i = 0; i < n; i++) if (codes[i] === C) sc.cIdx[nc++] = i; }
  else { for (let i = n - 1; i >= 0; i--) if (codes[i] === G) sc.cIdx[nc++] = i; }
  let mlAt = 0, implicit = false, any = false;
  for (const entry of mm.split(';')) {
    if (!entry) continue;
    const parts = entry.split(',');
    const head = parts[0];
    // head: base, strand, modification codes (letters, or one ChEBI number), optional '.' / '?'
    const m = /^([ACGTUN])([+-])([a-z]+|\d+)([.?]?)$/.exec(head);
    const nSkips = parts.length - 1;
    if (!m) { return null; }
    const codesInEntry = /^\d+$/.test(m[3]) ? [m[3]] : m[3].split('');
    const per = codesInEntry.length;
    if (m[1] !== 'C' || m[2] !== '+') { mlAt += nSkips * per; continue; }
    const im = codesInEntry.indexOf('m'), ih = codesInEntry.indexOf('h');
    if (im < 0 && ih < 0) { mlAt += nSkips * per; continue; }
    if (im >= 0 && m[4] !== '?') implicit = true;
    let p = -1;
    for (let k = 1; k <= nSkips; k++) {
      p += parseInt(parts[k], 10) + 1;
      if (!(p < nc)) return null;    // more listed C's than the SEQ has: the tags describe another sequence
      const q = sc.cIdx[p];
      if (!sc.listed[q]) { sc.listed[q] = 1; sc.pm[q] = 0; sc.ph[q] = 0; }
      if (ml) {
        if (im >= 0) sc.pm[q] = (ml[mlAt + (k - 1) * per + im] + 0.5) / 256;
        if (ih >= 0) sc.ph[q] = (ml[mlAt + (k - 1) * per + ih] + 0.5) / 256;
      }
      any = true;
    }
    mlAt += nSkips * per;
  }
  return any || implicit ? implicit : null;
}

/** The 5mC call of one base: whether it is modified, and the call's confidence (5hmC set aside, the two others renormalised). */
export function callOf(pm: number, ph: number): { mod: boolean; conf: number } | null {
  const pc = Math.max(0, 1 - pm - ph), s = pm + pc;
  if (s <= 1e-6) return null;
  const m = pm / s;
  return m > 0.5 ? { mod: true, conf: m } : { mod: false, conf: 1 - m };
}

/** CpG counts of one set of reads, with the histogram of their call confidences (for the filter threshold). */
export class MethylCounts {
  sites = new Map<number, CpgCounts>();
  hist = new Float64Array(CONF_BINS);
  calls = 0;
  get bytes(): number { return this.sites.size * (HAPS * 2 * CONF_BINS * 2 + 48); }
  add(pos: number, hap: number, mod: boolean, conf: number): void {
    let c = this.sites.get(pos);
    if (!c) this.sites.set(pos, c = new Uint16Array(HAPS * 2 * CONF_BINS));
    const b = confBin(conf), i = idx(hap, mod ? 1 : 0, b);
    if (c[i] < 65535) c[i]++;
    this.hist[b]++; this.calls++;
  }
}

/**
 * Counts one read's CpG calls. `ops` the packed CIGAR, `codes[0..n)` its SEQ as character codes, `cpgs` the sorted CpG
 * positions of the reference over it, `hap` its haplotag (0 untagged).
 */
export function countRead(out: MethylCounts, start: number, ops: ArrayLike<number>, codes: Uint8Array, n: number, reverse: boolean,
  mm: string, ml: ArrayLike<number> | null, cpgs: Int32Array, hap: number, sc: ModScratch): number {
  const h = hap > 0 && hap <= 2 ? hap : 0;
  return visitReadCalls(start, ops, codes, n, reverse, mm, ml, cpgs, sc, (cpg, mod, conf) => out.add(cpg, h, mod, conf));
}

/**
 * Each CpG call of one read, in reference order: `visit(cpg, modified, confidence)` with `cpg` the C of the CpG on the
 * + strand. Returns the number of calls.
 */
export function visitReadCalls(start: number, ops: ArrayLike<number>, codes: Uint8Array, n: number, reverse: boolean,
  mm: string, ml: ArrayLike<number> | null, cpgs: Int32Array, sc: ModScratch, visit: (cpg: number, mod: boolean, conf: number) => void): number {
  const implicit = parseModTags(mm, ml, codes, n, reverse, sc);
  if (implicit === null) return 0;
  let pos = start, q = 0, counted = 0;
  // first CpG at or after the read start
  let lo = 0, hi = cpgs.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (cpgs[mid] < start - 1) lo = mid + 1; else hi = mid; }
  let k = lo;
  for (let o = 0; o < ops.length; o++) {
    const v = ops[o], len = v >>> 4, op = v & 15;
    if (op === 0 || op === 7 || op === 8) {
      const bEnd = pos + len;
      // + strand: the C at c must be in the block; − strand: the G at c + 1
      while (k < cpgs.length && cpgs[k] + (reverse ? 1 : 0) < pos) k++;
      for (let j = k; j < cpgs.length; j++) {
        const at = cpgs[j] + (reverse ? 1 : 0);
        if (at >= bEnd) break;
        const qi = q + (at - pos);
        if (qi >= n || codes[qi] !== (reverse ? G : C)) continue;
        let pm = 0, ph = 0;
        if (sc.listed[qi]) { pm = sc.pm[qi]; ph = sc.ph[qi]; }
        else if (!implicit) continue;
        const call = callOf(pm, ph);
        if (!call) continue;
        visit(cpgs[j], call.mod, call.conf);
        counted++;
      }
      pos = bEnd; q += len;
    } else if (op === 1 || op === 4) q += len;
    else if (op === 2 || op === 3) pos += len;
  }
  return counted;
}

/** One window read back: per CpG, per haplotag, the modified and total calls passing the threshold. */
export interface MethylWindow {
  start: number; end: number;
  /** CpG positions (C of the CpG, 0-based) with at least one call */
  pos: Int32Array;
  /** per haplotag 0 untagged, 1, 2: modified calls and all passing calls at each CpG */
  mod: [Uint16Array, Uint16Array, Uint16Array];
  total: [Uint16Array, Uint16Array, Uint16Array];
  /** confidence threshold applied (the lower edge of its bin), and the calls it left out of the window */
  threshold: number; filtered: number; calls: number;
  /** CpG islands of the reference over the window (Gardiner-Garden & Frommer 1987) */
  islands: [number, number][];
  /** CpG sites of the reference in the window, called or not */
  cpgs: number;
}

/** The threshold bin: the 10th percentile of the call confidences of `hist`. */
export function thresholdBin(hist: Float64Array): number {
  let total = 0;
  for (const x of hist) total += x;
  if (!total) return 0;
  let acc = 0;
  for (let b = 0; b < hist.length; b++) { acc += hist[b]; if (acc >= total * FILTER_PERCENTILE) return b; }
  return 0;
}

/** Reads a window back from counts (layers each holding reads once), applying the threshold of their combined histogram. */
export function methylWindow(layers: MethylCounts[], start: number, end: number, refSeq: string | null, refStart: number): MethylWindow {
  const hist = new Float64Array(CONF_BINS);
  for (const l of layers) for (let b = 0; b < CONF_BINS; b++) hist[b] += l.hist[b];
  const tb = thresholdBin(hist);
  const merged = new Map<number, Uint16Array>();
  for (const l of layers) for (const [p, c] of l.sites) {
    if (p < start || p >= end) continue;
    const m = merged.get(p);
    if (!m) merged.set(p, c.slice()); else for (let i = 0; i < c.length; i++) m[i] = Math.min(65535, m[i] + c[i]);
  }
  const pos = Int32Array.from([...merged.keys()].sort((a, b) => a - b));
  const mod = [new Uint16Array(pos.length), new Uint16Array(pos.length), new Uint16Array(pos.length)] as MethylWindow['mod'];
  const total = [new Uint16Array(pos.length), new Uint16Array(pos.length), new Uint16Array(pos.length)] as MethylWindow['total'];
  let filtered = 0, calls = 0;
  pos.forEach((p, i) => {
    const c = merged.get(p)!;
    for (let h = 0; h < HAPS; h++) for (let st = 0; st < 2; st++) for (let b = 0; b < CONF_BINS; b++) {
      const x = c[idx(h, st, b)];
      if (!x) continue;
      calls += x;
      if (b < tb) { filtered += x; continue; }
      total[h][i] += x; if (st === 1) mod[h][i] += x;
    }
  });
  let cpgs = 0;
  if (refSeq) for (let i = Math.max(0, start - refStart); i + 1 < refSeq.length && refStart + i < end; i++) if ((refSeq.charCodeAt(i) & 0xdf) === C && (refSeq.charCodeAt(i + 1) & 0xdf) === G) cpgs++;
  return { start, end, pos, mod, total, threshold: 0.5 + tb / (2 * CONF_BINS), filtered, calls, islands: refSeq ? cpgIslands(refSeq, refStart, start, end) : [], cpgs };
}

/**
 * CpG islands of a reference stretch by the criteria of Gardiner-Garden & Frommer (1987): 200 bp windows with a GC
 * content of at least 50 % and an observed / expected CpG ratio of at least 0.6, merged where they overlap.
 */
export function cpgIslands(seq: string, seqStart: number, start: number, end: number): [number, number][] {
  const W = 200, a0 = Math.max(0, start - seqStart), a1 = Math.min(seq.length, end - seqStart);
  if (a1 - a0 < W) return [];
  const up = (i: number) => seq.charCodeAt(i) & 0xdf;
  let c = 0, g = 0, cg = 0;
  const addAt = (i: number, s: number) => {
    const x = up(i);
    if (x === C) c += s; else if (x === G) g += s;
    if (i + 1 < a1 && x === C && up(i + 1) === G) cg += s;
  };
  for (let i = a0; i < a0 + W; i++) addAt(i, 1);
  const out: [number, number][] = [];
  let open = -1, lastEnd = -1;
  for (let i = a0; i + W <= a1; i++) {
    if (i > a0) { addAt(i - 1, -1); addAt(i + W - 1, 1); }
    // the CpG straddling the right edge was counted with its C only once the G came in: close enough for a 200 bp scan
    const gc = (c + g) / W, oe = c && g ? (cg * W) / (c * g) : 0;
    const ok = gc >= 0.5 && oe >= 0.6;
    if (ok) { if (open < 0) open = i; lastEnd = i + W; }
    else if (open >= 0 && i >= lastEnd) { out.push([seqStart + open, seqStart + lastEnd]); open = -1; }
  }
  if (open >= 0) out.push([seqStart + open, seqStart + lastEnd]);
  return out;
}

/** Methylation counts of one sample over one stretch of a chromosome: owned reads start in [ps, pe), the spill reaches in (see CoverageState). */
export class MethylState {
  owned = new MethylCounts();
  spill = new MethylCounts();
  ps = 0;
  pe = 0;
  lastUsed = 0;
  constructor(readonly chrom: string) {}
  get empty(): boolean { return this.pe <= this.ps; }
  get bytes(): number { return this.owned.bytes + this.spill.bytes; }
}
