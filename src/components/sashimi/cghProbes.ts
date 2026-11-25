/**
 * The probes' log ratios of an array (CytoGenomics' CGH_Log_Ratio.txt: `AMADID <design> hg19`, then Probe Name, ChrName,
 * Start, Stop, FeatureNum, Description, Name of Gene, Normalized log ratio(<sample>); NaN for a probe left out). Kept
 * as the file's bytes and read on first use into typed arrays per chromosome (about 12 bytes a probe), for a probe
 * the probe panel above the browser (CghProbePanel.tsx) and the DLRS.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { decodeText } from './cghFormats';
import { cghBuild, cghChrom, leafOf } from './cghIntervals';
import { norm, splitDelimited } from './cghMapping';

export interface CghProbeChrom { start: Uint32Array; end: Uint32Array; lr: Float32Array; name: string[] }
export interface CghProbes {
  sample: string; build: GenomeBuild | null;
  /** probes read (NaN included) */
  n: number;
  /** chr1…chrY, each sorted by start */
  byChrom: Map<string, CghProbeChrom>;
  /** derivative log ratio spread over all chromosomes (null: under 10 probes) */
  dlrs: number | null;
}

const LR_HEAD = /normalized\s+log\s*ratio\s*\(([^)]*)\)/i;
/** the sample a probe log-ratio file is of (its `Normalized log ratio(<sample>)` header), or null when the bytes are not one */
export function probeSampleOf(bytes: Uint8Array): string | null {
  const head = decodeText(bytes.subarray(0, 8192));
  for (const l of head.split(/\r\n|\n|\r/).slice(0, 30)) if (/^probe\s*name\t/i.test(l)) return l.match(LR_HEAD)?.[1].trim() ?? null;
  return null;
}

/** the file read: probes per chromosome, sorted by start (`spec`: any other layout, by its columns) */
export function parseProbes(bytes: Uint8Array, spec?: ProbeSpec): CghProbes {
  if (spec) return parseBySpec(bytes, spec);
  const lines = decodeText(bytes).split(/\r\n|\n|\r/);
  const hi = lines.findIndex(l => /^probe\s*name\t/i.test(l));
  if (hi < 0) throw new Error('not a probe log-ratio file');
  const h = lines[hi].split('\t').map(x => x.trim().toLowerCase());
  const cChr = h.findIndex(x => x === 'chrname' || x === 'chr' || x === 'chromosome'), cS = h.indexOf('start'), cE = h.indexOf('stop'), cL = h.findIndex(x => LR_HEAD.test(x));
  const sample = lines[hi].match(LR_HEAD)?.[1].trim() ?? '';
  const build = cghBuild(lines.slice(0, hi).join(' '));
  const acc = new Map<string, { s: number[]; e: number[]; v: number[]; p: string[] }>();
  const cName = Math.max(0, h.findIndex(x => x === 'probe name' || x === 'probename'));
  let n = 0;
  for (let i = hi + 1; i < lines.length; i++) {
    const f = lines[i].split('\t');
    const chrom = cghChrom(f[cChr] ?? ''), s = Number(f[cS]), e = Number(f[cE]);
    if (!chrom || chrom === 'chrM' || !Number.isFinite(s) || !Number.isFinite(e)) continue;
    let a = acc.get(chrom); if (!a) acc.set(chrom, a = { s: [], e: [], v: [], p: [] });
    const v = Number(f[cL]);
    a.s.push(s); a.e.push(e); a.v.push(Number.isFinite(v) ? v : NaN); a.p.push(f[cName] ?? ''); n++;
  }
  return sorted(acc, sample, build, n);
}
type Acc = Map<string, { s: number[]; e: number[]; v: number[]; p: string[] }>;
const CHR_RANK = (c: string) => { const x = c.slice(3); return x === 'X' ? 23 : x === 'Y' ? 24 : Number(x); };
/** each chromosome's probes sorted by start; `ordered`: the chromosomes too (chr1…chr22, X, Y), for a file in any order */
function sorted(acc: Acc, sample: string, build: GenomeBuild | null, n: number, ordered = false): CghProbes {
  const byChrom = new Map<string, CghProbeChrom>();
  for (const [chrom, a] of ordered ? [...acc].sort((x, y) => CHR_RANK(x[0]) - CHR_RANK(y[0])) : acc) {
    const order = a.s.map((_, i) => i).sort((x, y) => a.s[x] - a.s[y]);
    byChrom.set(chrom, { start: Uint32Array.from(order, i => a.s[i]), end: Uint32Array.from(order, i => a.e[i]), lr: Float32Array.from(order, i => a.v[i]), name: order.map(i => a.p[i]) });
  }
  return { sample, build, n, byChrom, dlrs: dlrsOf(byChrom) };
}

/**
 * Any other probe-level layout (Agilent Feature Extraction-like, Illumina Log R Ratio, Affymetrix ChAS-like, CNVkit .cnr,
 * DNAcopy / CBS inputs, a column per sample): the columns of one sample. `end` -1: one position column (start = end);
 * `zeroBased`: BED-like starts (CNVkit), moved to 1-based.
 */
export interface ProbeSpec { delim: string; headerRow: number; chrom: number; start: number; end: number; lr: number; name: number; sample: string; zeroBased: boolean; build: GenomeBuild | null; layout: string }
/** a log ratio as a cell writes it: empty, NA, NaN, null, '.' are NaN (a probe left out); a decimal comma when the delimiter is not one */
const lrOf = (t: string, dec: boolean) => { const v = t.trim().replace(/^"|"$/g, ''); if (!v || /^(na|nan|null|none|\.|-|inf|-inf)$/i.test(v)) return NaN; const n = Number(dec ? v.replace(',', '.') : v); return Number.isFinite(n) ? n : NaN; };
const posOf = (t: string | undefined) => { const v = (t ?? '').trim().replace(/^"|"$/g, ''); return /^\d+(\.0+)?$/.test(v) ? Number(v) : NaN; };
function parseBySpec(bytes: Uint8Array, sp: ProbeSpec): CghProbes {
  const lines = decodeText(bytes).replace(/^\uFEFF/, '').split(/\r\n|\n|\r/), acc: Acc = new Map(), dec = sp.delim !== ',', shift = sp.zeroBased ? 1 : 0;
  let n = 0;
  for (let i = sp.headerRow + 1; i < lines.length; i++) {
    const l = lines[i];
    if (!l || l[0] === '#') continue;
    const f = splitDelimited(l, sp.delim);
    const chrom = cghChrom((f[sp.chrom] ?? '').replace(/^"|"$/g, '')), s = posOf(f[sp.start]), e = sp.end >= 0 ? posOf(f[sp.end]) : s;
    if (!chrom || chrom === 'chrM' || !Number.isFinite(s) || !Number.isFinite(e)) continue;
    let a = acc.get(chrom); if (!a) acc.set(chrom, a = { s: [], e: [], v: [], p: [] });
    a.s.push(s + shift); a.e.push(Math.max(s + shift, e + (sp.end >= 0 ? 0 : shift))); a.v.push(lrOf(f[sp.lr] ?? '', dec)); a.p.push(sp.name >= 0 ? (f[sp.name] ?? '').replace(/^"|"$/g, '') : ''); n++;
  }
  return sorted(acc, sp.sample, sp.build, n, true);
}

/** a header cell naming a log-ratio column: the sample it names ('' for none), else null */
export function lrColumnSample(h: string): string | null {
  const t = h.trim().replace(/^"|"$/g, '');
  const pre = t.match(/^(?:normali[sz]ed\s*)?(?:log\s*r\s*ratio|log\s*2?\s*[-_]?\s*ratio|log2|lrr)\s*[([]\s*(.+?)\s*[)\]]$/i);
  if (pre) return pre[1];
  if (/^(lr|value|log2\s*copy\s*ratio)$/i.test(t)) return '';
  const suf = t.match(/^(.*?)[\s._:|-]*(log\s*r\s*ratio|weighted\s*log\s*2\s*ratio|log\s*2?\s*[-_]?\s*ratio|log\s*2\s*r|log2|logr|lrr|l2r)$/i);
  if (!suf) return null;
  return /^(normali[sz]ed|mean|smoothed|weighted|raw|cn|probe|segmented|seg)?$/i.test(suf[1].trim()) ? '' : suf[1].trim();
}
const SYN = {
  chrom: ['chromosome', 'chrom', 'chr', 'chrname', 'chromosomename', 'seqnames', 'seqname', 'contig'],
  start: ['start', 'position', 'pos', 'startposition', 'chromstart', 'chrstart', 'physicalposition', 'maploc', 'locstart', 'coordinate'],
  end: ['stop', 'end', 'endposition', 'chromend', 'chrend', 'locend'],
  name: ['probename', 'probeid', 'probesetname', 'probesetid', 'probe', 'probeset', 'name', 'snpname', 'snp', 'marker', 'markername', 'id'],
};
/** columns of a table of calls, never of probes: a type, a size, a number of probes / markers, a segment mean */
const CALL_COLS = new Set(['type', 'cnvtype', 'aberrationtype', 'call', 'calls', 'size', 'sizekb', 'sizebp', 'length', 'probes', 'nprobes', 'numprobes', 'numberofprobes', 'nummark', 'nummarkers', 'nbprobes', 'segmean', 'aberrationno']);
/** fewer data rows than this: a table of calls (a probe file has tens of thousands) */
export const PROBE_MIN_ROWS = 2000;
const countRows = (b: Uint8Array) => { let n = 0; for (let i = 0; i < b.length; i++) if (b[i] === 10) n++; return n; };
/** the file name without folders and extensions: a probe file's sample when nothing else names one */
const baseName = (fileName: string) => leafOf(fileName).replace(/(\.gz)?$/i, '').replace(/\.(txt|tsv|csv|cnr|cnn|tab|dat)$/i, '');

/**
 * The samples of a probe-level file other than CytoGenomics' CGH_Log_Ratio.txt: the columns chosen once in the mapping
 * dialog for this header (remembered by its column names), else recognised by its header (any delimiter; `#` comments or
 * `key: value` lines before it): a chromosome, a position (or start / end) and one or several log-ratio columns, at least
 * PROBE_MIN_ROWS rows and no column of a calls table; one spec per sample. Null: not one.
 */
export function probeLayoutOf(bytes: Uint8Array, fileName: string): ProbeSpec[] | null {
  const lines = decodeText(bytes.subarray(0, 65536)).replace(/^\uFEFF/, '').split(/\r\n|\n|\r/).slice(0, 80);
  const kept = findProbeMap(lines, fileName);
  if (kept) return kept;
  if (countRows(bytes) < PROBE_MIN_ROWS) return null;
  for (let hi = 0; hi < lines.length - 1; hi++) {
    const l = lines[hi];
    if (!l.trim() || l[0] === '#') continue;
    for (const delim of ['\t', ';', ',']) {
      if (!l.includes(delim)) continue;
      const raw = splitDelimited(l, delim), h = raw.map(norm);
      if (h.length < 3 || h.some(x => CALL_COLS.has(x))) continue;
      const col = (names: string[], not: number[] = []) => { for (const s of names) { const i = h.findIndex((x, j) => x === s && !not.includes(j)); if (i >= 0) return i; } return -1; };
      const chrom = col(SYN.chrom), start = col(SYN.start, [chrom]), end = col(SYN.end, [chrom, start]), name = col(SYN.name, [chrom, start, end]);
      if (chrom < 0 || start < 0) continue;
      const lrs = raw.map((x, j) => ({ j, s: [chrom, start, end, name].includes(j) ? null : lrColumnSample(x) })).filter((x): x is { j: number; s: string } => x.s != null);
      if (!lrs.length) continue;
      const cnvkit = h.includes('log2') && (h.includes('depth') || h.includes('weight')) || /\.(cnr|cnn)(\.gz)?$/i.test(fileName);
      const zeroBased = cnvkit || h[start] === 'chromstart';
      const layout = cnvkit ? 'CNVkit (.cnr)' : raw.some(x => /log\s*r\s*ratio/i.test(x)) ? 'Illumina-like (Log R Ratio)' : raw.some(x => /log2\s*ratio/i.test(x)) && h.some(x => x.startsWith('probeset')) ? 'Affymetrix ChAS-like' : h.includes('probename') ? 'Agilent-like' : 'generic (chromosome, position, log ratio)';
      return specsOf(lines, hi, delim, raw, { chrom, start, end, name, zeroBased, layout }, lrs, fileName);
    }
  }
  return null;
}
/** one spec per sample: columns naming no sample (Log2Ratio and WeightedLog2Ratio of one array) give one, the plain one first */
function specsOf(lines: string[], hi: number, delim: string, raw: string[], c: Pick<ProbeSpec, 'chrom' | 'start' | 'end' | 'name' | 'zeroBased' | 'layout'> & { build?: ProbeSpec['build'] }, lrs: { j: number; s: string }[], fileName: string): ProbeSpec[] {
  const pre = lines.slice(0, hi), build = cghBuild(pre.join(' ')) ?? c.build ?? null;
  // a `Sample: …` line before the header, else the file's name, for a column naming none
  const named = pre.map(l => l.replace(/^#+\s*/, '').match(/^(?:sample|sample\s*name|sample\s*id|patient|array\s*name|barcode)\s*[:=\t]\s*(.+?)\s*$/i)?.[1]).find(Boolean);
  const anon = lrs.filter(x => !x.s).sort((a, b) => Number(/weight|smooth/i.test(raw[a.j])) - Number(/weight|smooth/i.test(raw[b.j])))[0];
  const seen = new Set<string>(), out: ProbeSpec[] = [];
  for (const x of [...(anon ? [anon] : []), ...lrs.filter(y => y.s)]) {
    const sample = x.s || named || baseName(fileName);
    if (seen.has(sample)) continue;
    seen.add(sample);
    out.push({ delim, headerRow: hi, chrom: c.chrom, start: c.start, end: c.end, lr: x.j, name: c.name, sample, zeroBased: c.zeroBased, build, layout: c.layout });
  }
  return out;
}

// ---- remembered probe mappings (the dialog's "probe log ratios" mode): by the header's column names only ----
const PROBE_MAPS_KEY = 'sashimi.cghProbeMaps';
/** what is kept: the columns' indexes and the options; never a row, never a sample name typed */
export interface StoredProbeMap { sig: string; delim: string; chrom: number; start: number; end: number; lr: number[]; name: number; zeroBased: boolean; build: GenomeBuild | null }
const loadProbeMaps = (): StoredProbeMap[] => { try { const v = JSON.parse(localStorage.getItem(PROBE_MAPS_KEY) ?? '[]'); return Array.isArray(v) ? v : []; } catch { return []; } };
const sigOf = (cells: string[]) => cells.map(norm).join('|');
export function rememberProbeMap(header: string[], m: Omit<StoredProbeMap, 'sig'>) {
  const sig = sigOf(header);
  try { localStorage.setItem(PROBE_MAPS_KEY, JSON.stringify([{ ...m, sig }, ...loadProbeMaps().filter(x => x.sig !== sig)].slice(0, 20))); } catch { /* storage blocked */ }
}
/** the sample a log-ratio column is named for: the name in its header, else the header itself when there are several, else the file's */
export const sampleOfColumn = (h: string, several: boolean, fileName: string) => lrColumnSample(h) || (several ? h.trim() : '') || baseName(fileName);
function findProbeMap(lines: string[], fileName: string): ProbeSpec[] | null {
  for (const m of loadProbeMaps()) {
    if (!m?.delim || !Array.isArray(m.lr)) continue;
    const hi = lines.findIndex(l => l.includes(m.delim) && sigOf(splitDelimited(l, m.delim)) === m.sig);
    if (hi < 0) continue;
    const raw = splitDelimited(lines[hi], m.delim);
    return specsOf(lines, hi, m.delim, raw, { ...m, layout: 'columns chosen before', build: m.build ?? undefined }, m.lr.map(j => ({ j, s: sampleOfColumn(raw[j] ?? '', m.lr.length > 1, fileName) })), fileName);
  }
  return null;
}

/**
 * DLRS (derivative log ratio spread, Agilent's QC metric): the standard deviation of the differences between consecutive
 * probes (by position, within a chromosome, NaN skipped), divided by √2 (the difference of two noises doubles the variance).
 */
export function dlrsOf(byChrom: Map<string, CghProbeChrom>): number | null {
  let n = 0, sum = 0, sq = 0;
  for (const { lr } of byChrom.values()) {
    let prev = NaN;
    for (let i = 0; i < lr.length; i++) {
      const v = lr[i]; if (Number.isNaN(v)) continue;
      if (!Number.isNaN(prev)) { const d = v - prev; n++; sum += d; sq += d * d; }
      prev = v;
    }
  }
  return n < 10 ? null : Math.sqrt(Math.max(0, (sq - sum * sum / n) / (n - 1))) / Math.SQRT2;
}
/** Agilent guidance on the DLRS */
export const dlrsWord = (d: number) => (d < 0.2 ? 'excellent' : d <= 0.3 ? 'good' : 'poor');
export const DLRS_TIP = 'DLRS (derivative log ratio spread): the standard deviation of the log-ratio differences between consecutive probes, divided by √2; the noise of the array. Agilent guidance: below 0.2 excellent, 0.2 to 0.3 good, above 0.3 poor.';

/** the first probe whose start is ≥ pos (binary search) */
export function lowerBound(a: Uint32Array, pos: number): number {
  let lo = 0, hi = a.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < pos) lo = m + 1; else hi = m; }
  return lo;
}
/**
 * the mean log ratio of the probes in [from, to] (NaN skipped), and their number: a probe counts when its midpoint lies in
 * [from, to], ends included, so a probe straddling an end counts on the side holding most of it, never twice nor never
 * (CytoGenomics' calls run from their first probe's start to their last probe's stop: its #Probes, on every real call)
 */
export function meanIn(p: CghProbeChrom, from: number, to: number): { mean: number | null; n: number } {
  let n = 0, sum = 0, i = lowerBound(p.start, from);
  while (i > 0 && p.end[i - 1] >= from) i--;
  for (; i < p.start.length && p.start[i] <= to; i++) {
    const v = p.lr[i], m = (p.start[i] + p.end[i]) / 2;
    if (m >= from && m <= to && !Number.isNaN(v)) { n++; sum += v; }
  }
  return { mean: n ? sum / n : null, n };
}
/** the non-NaN probes of [from, to] (those starting before `from` but reaching into it too), with k/2 more on each side for the median's edges */
export function visible(p: CghProbeChrom, from: number, to: number, k: number) {
  let i0 = lowerBound(p.start, from), extra = 0;
  while (i0 > 0 && p.end[i0 - 1] >= from) i0--;
  while (i0 > 0 && extra < k >> 1) { i0--; if (!Number.isNaN(p.lr[i0])) extra++; }
  const idx: number[] = [];
  let after = 0;
  for (let i = i0; i < p.start.length; i++) {
    if (Number.isNaN(p.lr[i])) continue;
    idx.push(i);
    if (p.start[i] > to && ++after > k >> 1) break;
  }
  return idx;
}
/** the running median of k probes (k odd), centred; at the ends of the list over the probes there are (fewer) */
export function runningMedian(v: number[], k: number): number[] {
  const h = k >> 1, out = new Array<number>(v.length), w: number[] = [];
  for (let i = 0; i < v.length; i++) {
    w.length = 0;
    for (let j = Math.max(0, i - h); j <= Math.min(v.length - 1, i + h); j++) w.push(v[j]);
    w.sort((a, b) => a - b);
    out[i] = w.length % 2 ? w[w.length >> 1] : (w[w.length / 2 - 1] + w[w.length / 2]) / 2;
  }
  return out;
}

/** the probe columns of a table by names, then by content (the mapping dialog's probe mode): chromosome, position, end, probe name, log ratio(s) */
export function guessProbeCols(header: string[], rows: string[][]): { chrom: number; start: number; end: number; name: number; lr: number[] } {
  const h = header.map(norm), sample = rows.slice(0, 50), used = new Set<number>();
  const col = (names: string[]) => { for (const s of names) { const i = h.findIndex((x, j) => x === s && !used.has(j)); if (i >= 0) { used.add(i); return i; } } return -1; };
  const cells = (j: number) => sample.map(r => (r[j] ?? '').trim().replace(/^"|"$/g, '')).filter(Boolean);
  const share = (j: number, ok: (v: string) => boolean) => { const v = cells(j); return v.length ? v.filter(ok).length / v.length : 0; };
  const free = (ok: (j: number) => boolean) => { for (let j = 0; j < header.length; j++) if (!used.has(j) && ok(j)) { used.add(j); return j; } return -1; };
  const isInt = (v: string) => /^\d{3,}$/.test(v), isLr = (v: string) => /^(na|nan|null)$/i.test(v) || (/^[-+]?\d*[.,]?\d+(e[-+]?\d+)?$/i.test(v) && Math.abs(Number(v.replace(',', '.'))) <= 10);
  let chrom = col(SYN.chrom); if (chrom < 0) chrom = free(j => share(j, v => /^(chr)?(\d{1,2}|X|Y)$/i.test(v)) >= 0.8);
  let start = col(SYN.start); if (start < 0) start = free(j => share(j, isInt) >= 0.95);
  let end = col(SYN.end); if (end < 0 && start >= 0 && start + 1 < header.length && !used.has(start + 1) && share(start + 1, isInt) >= 0.95 && sample.every(r => !isInt(r[start + 1] ?? '') || Number(r[start + 1]) >= Number(r[start]))) { end = start + 1; used.add(end); }
  const name = col(SYN.name);
  let lr = header.flatMap((x, j) => (!used.has(j) && lrColumnSample(x) != null ? [j] : []));
  if (!lr.length) lr = header.flatMap((_, j) => (!used.has(j) && share(j, isLr) >= 0.9 && cells(j).some(v => /[.,]/.test(v)) ? [j] : []));
  return { chrom, start, end, name, lr };
}
