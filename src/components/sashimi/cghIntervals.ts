/**
 * Array-CGH calls for the SV browser: the interval tables exported per sample by Agilent CytoGenomics (a few header
 * lines, `Genome hg19`, `Record Name …`, then one row per call: Chromosome, Start, Stop, Cytoband, Gene Name, Size(kb),
 * Type, #Probes, Mean Log Ratio, ISCN, Suppress, Classification), any other table once its columns are mapped, or a BED
 * file (chrom, start, end, name). The SV VCF is read in cghVcf.ts, the bytes (encodings, .xlsx, .gz) in cghFormats.ts.
 *
 * Columns are found by their header, not by their place. Positions of the interval table are read as the ISCN of the
 * same row writes them (1-based, the first and last probe of the call); a BED's start is 0-based and moved to 1-based.
 * The copy number comes from the ISCN (x0, x1, x3, x4…; `amp` for a high-level gain), the type from the Type column
 * (else from the copy number, else from the log ratio's sign). A file is never sent anywhere: the calls stay in this
 * tab, with the classification typed in the table, and can be written back as the same table.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { norm, parseNum, parsePos, parseRegion, copiesOf, splitDelimited, joinDelimited, splitRecords, eolOf, detectTable, autoCols, AUTO_OPTS, SIZE_MULT, NeedsMapping, type MapSpec } from './cghMapping';

export type CghKind = 'loss' | 'gain';
export type CghSex = 'male' | 'female';
/** a sex as a header line may write it (Gender / Sex: Male, F, Femme…) */
const sexOf = (v: string | null | undefined): CghSex | null => (!v ? null : /^\s*(m|male|homme|masculin)\b/i.test(v) ? 'male' : /^\s*(f|female|femme|f[ée]minin)\b/i.test(v) ? 'female' : null);
export interface CghCnv {
  /** unique within the tab: file, row */
  id: string;
  sample: string;
  /** the sample's key (its file and name) */
  sampleKey: string;
  /** chr1…chr22, chrX, chrY */
  chrom: string;
  /** 1-based inclusive, as in the file's ISCN */
  start: number; end: number;
  /** as written: p31.1, p11.2 - p11.1 */
  cytoband: string;
  genes: string[];
  /** the file's Size(kb), else (end - start) / 1000 */
  sizeKb: number;
  kind: CghKind;
  /** the Type column as written (Deletion, Amplification), or the BED name */
  typeText: string;
  /** null when the file gives no probe count (a BED) */
  probes: number | null;
  logRatio: number | null;
  iscn: string | null;
  /** copies from the ISCN (x1 → 1); null for `amp` or no ISCN; `copiesHi` for a range (x1~2, mosaic) */
  copies: number | null; copiesHi: number | null; amp: boolean;
  suppressed: boolean;
  /** the classification as given (the text before `;;`), and what followed it, kept for the export */
  classification: string; classRest: string;
  /** the file's comment on the call (a Comment / Note column), else '' */
  note: string;
  /** the file's line of this call (0-based), for the export */
  line: number;
  /** the call's p-value, when the file gives one (CytoGenomics' interval-based report) */
  pval?: number | null;
  /** detected from the probes (cghSegment.ts), not in the file: drawn, never exported */
  fromProbes?: boolean;
}
export interface CghSample { key: string; name: string; fileId: string; declaredCalls: number | null; cnvs: CghCnv[];
  /** the sex a header line gives (Gender, Sex), else null */
  sex: CghSex | null }
export interface CghFile {
  /** the file's name (for a file of a zip: zip › path inside), and the name its export is saved under */
  id: string; name: string; saveAs: string; format: 'interval' | 'bed' | 'vcf';
  /** the build the file states (Genome hg19 → GRCh37), null when it states none */
  build: GenomeBuild | null; buildText: string | null;
  /** the header lines as key / value (Array Design Info, Creation Date…) */
  meta: [string, string][];
  lines: string[]; eol: string;
  /** the Classification column (interval table), or null (a BED: the class is added as a last column; a VCF: in INFO) */
  classCol: number | null;
  /** whether the file's classifications end with `;;` (written the same way back) */
  classSuffix: string;
  /** the header line (interval table) and its Comment / Note column, if any (else one is added at the export) */
  headerLine: number | null; commentCol: number | null;
  samples: CghSample[];
  warnings: string[];
  /** the delimiter (tab for the software's own export), and the columns this file was read with (null for a BED and an interval-based report) */
  delim: string; spec: MapSpec | null;
  /** the zip's patient name, to read the file again with other columns */
  fallback?: string;
  /** rows not read (unreadable, or a Type value set to ignore) and the first reason; whether a remembered mapping was used */
  skipped: number; firstSkip: string | null; remembered?: boolean;
}

/** chr1, 1, CHR1, Chr1, chrX, X, 23, 24 → chr1, chrX, chrY; chrM for M / MT (skipped by the readers); null for anything else */
export function cghChrom(s: string): string | null {
  const m = s.trim().match(/^(?:chr)?(\d{1,2}|X|Y|M|MT)$/i);
  if (!m) return null;
  const c = m[1].toUpperCase();
  if (c === '23') return 'chrX';
  if (c === '24') return 'chrY';
  if (/^\d+$/.test(c) && (+c < 1 || +c > 22)) return null;
  return `chr${c === 'MT' ? 'M' : c}`;
}
export function cghBuild(s: string | null | undefined): GenomeBuild | null {
  if (!s) return null;
  if (/hg19|grch37|b37/i.test(s)) return 'GRCh37';
  if (/hg38|grch38/i.test(s)) return 'GRCh38';
  return null;
}
/** copies from an ISCN: …(72523846_72667555)x1, x1~2, amp */
function iscnCopies(iscn: string | null): { copies: number | null; copiesHi: number | null; amp: boolean } {
  const m = iscn?.match(/\)\s*x\s*(\d+)(?:\s*~\s*(\d+))?/i);
  if (m) return { copies: +m[1], copiesHi: m[2] ? +m[2] : null, amp: false };
  return { copies: null, copiesHi: null, amp: !!iscn && /\)\s*amp\b/i.test(iscn) };
}
/** a loss or a gain as a Type cell says it (deletion, CN loss, hemizygous, perte, x1, CN=1, 1 copy; dup, amp, trip, CN gain, x3…) */
export const kindOfText = (raw: string): CghKind | null => {
  const t = raw.normalize('NFD').replace(/[\u0300-\u036f]/g, '');
  return /del|loss|perte|hemi|\bx[01]\b|\bcn\s*[=:]?\s*[01]\b|\b[01]\s*cop/i.test(t) ? 'loss'
    : /amp|gain|dup|trip|\bx[3-9]\b|\bcn\s*[=:]?\s*([3-9]|\d\d)\b|\b[3-9]\s*cop/i.test(t) ? 'gain' : null;
};
/** a copy-neutral call (LOH, ROH, AOH, cnLOH, UPD): skipped, counted */
export const isLohText = (t: string) => /\b(?:cn-?)?loh\b|\b[ar]oh\b|heterozygosity|homozygosity|copy[- ]?neutral|uniparental|\bupd\b/i.test(t);
/** a count said in a warning: 3 rows, 1 row */
const count = (n: number, one: string, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;
/** the classification cell: the class before `;;`, and the rest */
function splitClass(cell: string): { classification: string; classRest: string } {
  const i = cell.indexOf(';;');
  return i < 0 ? { classification: cell.trim(), classRest: '' } : { classification: cell.slice(0, i).trim(), classRest: cell.slice(i) };
}

let fileSeq = 0;
/** a file's id, unique within the tab */
export const nextFileId = (prefix = 'cgh') => `${prefix}${++fileSeq}`;
/** a BED of three columns: no type to read */
class UntypedBed extends Error {}
/** the file's name without its folders (a zip's path) */
export const leafOf = (fileName: string) => fileName.split(/[›/\\]/).pop()!.trim();
/**
 * One file's text: an interval table (CytoGenomics, or mapped by `spec`) or a BED. Throws when no call can be read,
 * NeedsMapping when it is a delimited table whose columns are not recognised. The sample is the table's Record Name,
 * else `fallbackSample` (a zip's name: one patient per zip), else the file's name. An .xlsx (read as tab-separated
 * text) is exported as .classified.txt.
 */
export function parseCghFile(text: string, fileName: string, fallbackSample?: string, spec?: MapSpec): CghFile {
  const eol = eolOf(text), lines = splitRecords(text);
  const id = nextFileId();
  const leaf = leafOf(fileName);
  const base = fallbackSample ?? leaf.replace(/\.(txt|tsv|csv|bed|xlsx?)$/i, '');
  let f: Omit<CghFile, 'saveAs'>;
  const ab = spec ? -1 : aberrationHeader(lines);
  if (ab >= 0) f = parseAberrations(lines, id, fileName, base, eol, ab);
  else if (spec) f = parseInterval(lines, id, fileName, base, eol, spec);
  else {
    const hi = lines.findIndex(l => { const c = l.split('\t').map(norm); return c.includes('chromosome') && c.includes('start') && (c.includes('stop') || c.includes('end')); });
    try {
      f = hi >= 0 ? parseInterval(lines, id, fileName, base, eol, { cols: autoCols(lines[hi].split('\t').map(norm)), opts: AUTO_OPTS(hi) }) : parseBed(lines, id, fileName, base, eol);
    } catch (e) {
      // a delimited table whose columns are not recognised (or with no header row): the caller asks which column is which
      const t = !(e instanceof UntypedBed) && detectTable(lines);
      if (t) throw new NeedsMapping((e as Error).message);
      throw e;
    }
  }
  return { ...f, fallback: fallbackSample, saveAs: leaf.replace(/\.xlsx$/i, '.txt').replace(/(\.[^.]+)?$/, m => `.classified${m || '.txt'}`) };
}

function parseInterval(lines: string[], id: string, fileName: string, base: string, eol: string, spec: MapSpec): Omit<CghFile, 'saveAs'> {
  const { cols: c, opts } = spec, hi = opts.headerRow, d = opts.delim;
  const sp = (l: string) => splitDelimited(l, d);
  // a decimal comma only where the delimiter is not the comma (a log ratio, a size, a copy number)
  const dec = d !== ',', nn = (s: string | undefined) => parseNum(s, dec);
  const meta: [string, string][] = [];
  for (const l of lines.slice(0, Math.max(0, hi))) {
    const [k, ...v] = sp(l);
    // `Key: value` in one cell (a CSV's header lines) as well as Key, value
    const kv = !v.join('').trim() && k?.match(/^\s*([^:]+?)\s*:\s*(.+)$/);
    if (kv) meta.push([kv[1], kv[2].trim()]);
    else if (k?.trim()) meta.push([k.trim().replace(/:$/, ''), v.join(' ').trim()]);
  }
  const metaOf = (re: RegExp) => meta.find(([k]) => re.test(k))?.[1] || null;
  const buildText = metaOf(/^genome$/i);
  const record = metaOf(/^record name$/i);
  const warnings: string[] = [];
  let skipped = 0, firstSkip: string | null = null, ignored = 0, loh = 0, mito = 0, swapped = 0, untyped = 0, normal = 0, iscnBuild: GenomeBuild | null = null;
  // the log-ratio columns (one, or gain / loss ones): each with the type its name says (Gain, Loss…), else the sign's
  const head = hi >= 0 ? sp(lines[hi]) : [];
  const lrCols = [c.lr, ...(spec.lrMore ?? [])].filter((j, k, a) => j >= 0 && a.indexOf(j) === k).map(j => { const h = head[j] ?? ''; return { j, kind: /del|loss|perte/i.test(h) && /dup|gain|amp/i.test(h) ? null : kindOfText(h), name: h.trim() }; });
  const sex = sexOf(meta.find(([k]) => /^(gender|sex|sexe)$/i.test(k))?.[1]);
  const skip = (why: string) => { skipped++; firstSkip ??= why; };
  if (c.probes < 0) warnings.push('no #Probes column: the minimum number of probes cannot be applied');
  // stored mappings from before the one-cell position have no `pos`
  const posCol = c.pos ?? -1;
  const samples = new Map<string, CghSample>();
  const declared = new Map<string, number>();
  let suffixYes = 0, suffixNo = 0;
  for (let i = hi + 1; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim()) continue;
    // the footer's count per sample (26A0657=36): a line of one cell, not a call's CN=1
    const foot = sp(l).filter(x => x.trim()).length === 1 ? l.match(/^\s*([^\t=]+?)\s*=\s*(\d+)\s*$/) : null;
    if (foot) { declared.set(foot[1], +foot[2]); continue; }
    if (/^number of calls/i.test(l)) continue;
    const f = sp(l);
    const iscnRaw = c.iscn >= 0 ? (f[c.iscn] ?? '').trim() : '';
    const iscn = iscnRaw && !/^na$/i.test(iscnRaw) ? iscnRaw : null;
    // the position: Chromosome / Start / End, else one cell (chr1:1-2), else the ISCN's (1-based as written)
    let chromText = c.chrom >= 0 ? f[c.chrom] ?? '' : '', s0 = c.start >= 0 ? parsePos(f[c.start]) : null, end = c.stop >= 0 ? parsePos(f[c.stop]) : null, shift = opts.zeroBased;
    if (!cghChrom(chromText) || s0 == null || end == null) {
      const r = (posCol >= 0 ? parseRegion(f[posCol]) : null) ?? (iscn ? parseRegion(iscn) : null);
      if (r) { chromText = r.chrom; s0 = r.start; end = r.end; shift = opts.zeroBased && !r.iscn; }
    }
    const chrom = cghChrom(chromText);
    if (chrom === 'chrM') { mito++; skip(`line ${i + 1}: chrM`); continue; }
    if (s0 != null && end != null && end < s0) { [s0, end] = [end, s0]; swapped++; }
    const start = s0 != null && shift ? s0 + 1 : s0;
    if (!chrom || start == null || end == null || end < start) { skip(`line ${i + 1}: ${!chrom ? 'chromosome' : 'start / end'} not a position (${l.slice(0, 40)})`); warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    if (iscn) iscnBuild ??= cghBuild(iscn.match(/arr\s*\[([^\]]+)\]/i)?.[1]);
    const typeText = c.type >= 0 ? (f[c.type] ?? '').trim() : '';
    let cn = iscnCopies(iscn);
    if (cn.copies == null && !cn.amp) { const n = c.copies >= 0 ? copiesOf(f[c.copies]) : null; if (n != null) cn = { copies: n, copiesHi: null, amp: false }; }
    // the non-zero one (the largest in absolute value if several), as in the interval-based report
    let lr: number | null = null, lrBy: (typeof lrCols)[number] | null = null;
    for (const x of lrCols) { const v = nn(f[x.j]); if (v != null && (lr == null || Math.abs(v) > Math.abs(lr))) { lr = v; lrBy = x; } }
    const choice = typeText ? opts.typeMap[typeText.toLowerCase()] ?? (opts.ignoreOtherTypes && !kindOfText(typeText) && !isLohText(typeText) ? 'ignore' : undefined) : undefined;
    if (choice === 'ignore') { ignored++; skip(`line ${i + 1}: type "${typeText}" ignored`); continue; }
    if (!choice && isLohText(typeText)) { loh++; skip(`line ${i + 1}: copy-neutral (${typeText})`); continue; }
    // no type: the copy number against the normal one (a male's X, and a Y, have one), else the log ratio's sign (0: none)
    const cn2 = chrom === 'chrY' || (chrom === 'chrX' && sex === 'male') ? 1 : 2;
    const byCn = cn.amp ? 'gain' : cn.copies != null ? (cn.copies < cn2 ? 'loss' : cn.copies > cn2 ? 'gain' : 'normal') : null;
    const kindRaw = choice ?? kindOfText(typeText) ?? byCn ?? (lr ? lrBy?.kind ?? (lr < 0 ? 'loss' : 'gain') : null);
    if (kindRaw === 'normal') { normal++; skip(`line ${i + 1}: copy number ${cn.copies}, the normal one`); continue; }
    if (!kindRaw) { untyped++; skip(`line ${i + 1}: no type, copy number or non-zero log ratio`); continue; }
    const kind: CghKind = kindRaw;
    const cell = c.cls >= 0 ? (f[c.cls] ?? '') : '';
    if (cell.trim()) { if (cell.includes(';;')) suffixYes++; else suffixNo++; }
    const name = (c.sample >= 0 && f[c.sample]?.trim()) || record || base;
    let s = samples.get(name);
    if (!s) { s = { key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs: [], sex }; samples.set(name, s); }
    const size = c.size >= 0 ? nn(f[c.size]) : null;
    s.cnvs.push({
      id: `${id}:${i}`, sample: name, sampleKey: s.key, chrom, start, end,
      cytoband: c.band >= 0 ? (f[c.band] ?? '').trim() : '',
      genes: c.genes >= 0 ? (f[c.genes] ?? '').split(/[,;]/).map(g => g.trim()).filter(Boolean) : [],
      sizeKb: size != null ? size * SIZE_MULT[opts.sizeUnit] : (end - start) / 1000, kind, typeText: typeText || (lrBy?.kind && lr ? lrBy.name : kind === 'loss' ? 'Deletion' : 'Gain'),
      probes: c.probes >= 0 ? nn(f[c.probes]) : null, logRatio: lr, iscn, ...cn,
      suppressed: c.sup >= 0 && /^(true|yes|1)$/i.test((f[c.sup] ?? '').trim()),
      ...splitClass(cell), note: c.note >= 0 ? (f[c.note] ?? '').trim() : '', line: i,
    });
  }
  for (const s of samples.values()) {
    s.declaredCalls = declared.get(s.name) ?? (samples.size === 1 && declared.size === 1 ? [...declared.values()][0] : null);
    // the footer counts the calls the software made; a different count is said, not corrected
    if (s.declaredCalls != null && s.declaredCalls !== s.cnvs.length) warnings.push(`${s.name}: the file announces ${s.declaredCalls} calls, ${s.cnvs.length} were read`);
  }
  if (untyped) warnings.push(`${count(untyped, 'row')} skipped: no type could be deduced (no type, no copy number, log ratio 0 or missing)`);
  if (normal) warnings.push(`${count(normal, 'row')} skipped: copy number normal (2; 1 on a male's X or a Y)`);
  if (ignored) warnings.push(`${count(ignored, 'row')} ignored by their type`);
  if (loh) warnings.push(`${count(loh, 'LOH / ROH row')} skipped (copy-neutral)`);
  if (mito) warnings.push(`${count(mito, 'call')} on chrM skipped`);
  if (swapped) warnings.push(`${count(swapped, 'row')} with start after end: the two swapped`);
  if (!samples.size) throw new Error(`${fileName}: no call could be read (an interval table needs Chromosome, Start, Stop, or a position in one column, and a Type, a copy number, an ISCN or a log ratio).`);
  return { id, name: fileName, format: 'interval', build: cghBuild(buildText) ?? iscnBuild ?? opts.build, buildText, meta, lines, eol, classCol: c.cls >= 0 ? c.cls : null, classSuffix: suffixYes > suffixNo ? ';;' : '', headerLine: hi >= 0 ? hi : null, commentCol: c.note >= 0 ? c.note : null, samples: [...samples.values()], warnings, delim: d, spec, skipped, firstSkip };
}

/** the four columns of CytoGenomics' interval-based report holding a call's log ratio, and the type each one says */
const ABERR_COLS: [string, CghKind][] = [['amplification', 'gain'], ['gain', 'gain'], ['loss', 'loss'], ['deletion', 'loss']];
/**
 * The header line of CytoGenomics' automated interval-based report (IntervalBasedReport, often named .xls although it is
 * tab-separated text): AberrationNo, Chr, Start, Stop, #Probes and Amplification / Gain / Loss / Deletion; -1 when none.
 */
export function aberrationHeader(lines: string[]): number {
  for (let i = 0; i < Math.min(lines.length, 300); i++) {
    if (!/aberrationno/i.test(lines[i])) continue;
    const c = lines[i].split('\t').map(norm);
    if (['aberrationno', 'chr', 'start', 'stop', 'probes'].every(k => c.includes(k)) && ABERR_COLS.some(([k]) => c.includes(k))) return i;
  }
  return -1;
}
/**
 * The interval-based report: `Key: value` header lines (Genome: hg19, Array Design Info: sample[design]…), the table's
 * header, then for each sample a line holding only its name and its calls, each with its log ratio in ONE of
 * Amplification / Gain / Loss / Deletion (the others 0): the non-zero one (the largest in absolute value if several),
 * the type from its column (Amplification / Gain: gain; Loss / Deletion: loss), else from its sign. No copy number.
 */
function parseAberrations(lines: string[], id: string, fileName: string, base: string, eol: string, hi: number): Omit<CghFile, 'saveAs'> {
  const meta: [string, string][] = [];
  for (const l of lines.slice(0, hi)) {
    const m = l.match(/^\s*([^:\t]+?)\s*(?::|\t)\s*(.*)$/);
    if (m) meta.push([m[1].trim(), m[2].replace(/\t+/g, ' ').trim()]);
  }
  const metaOf = (re: RegExp) => meta.find(([k]) => re.test(k))?.[1] || null;
  const buildText = metaOf(/^genome$/i);
  // Array Design Info: 26X0001[012345_20200101]: the sample before its design
  const design = metaOf(/^array design info$/i)?.match(/^\s*([^[\],;]+?)\s*\[/)?.[1] ?? null;
  const head = lines[hi].split('\t').map(norm), col = (k: string) => head.indexOf(k);
  const cChr = col('chr'), cBand = col('cytoband'), cStart = col('start'), cStop = col('stop'), cProbes = col('probes'), cP = col('pval'), cGenes = head.findIndex(h => h === 'genenames' || h === 'genename' || h === 'genes');
  const lrCols = ABERR_COLS.map(([k, kind]) => ({ i: col(k), kind, name: k[0].toUpperCase() + k.slice(1) })).filter(x => x.i >= 0);
  const warnings: string[] = [], samples = new Map<string, CghSample>(), declared = new Map<string, number>();
  let skipped = 0, firstSkip: string | null = null, mito = 0, current: string | null = null;
  const skip = (why: string) => { skipped++; firstSkip ??= why; };
  for (let i = hi + 1; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim() || /^number of calls/i.test(l)) continue;
    const f = l.split('\t'), filled = f.filter(x => x.trim());
    if (filled.length === 1) {
      // the footer's count per sample (26X0001 =12), else the name of the sample whose calls follow
      const foot = l.match(/^\s*([^\t=]+?)\s*=\s*(\d+)\s*$/);
      if (foot) declared.set(foot[1], +foot[2]); else current = filled[0].trim();
      continue;
    }
    const chrom = cghChrom(f[cChr] ?? ''), start = parsePos(f[cStart]), end = parsePos(f[cStop]);
    if (chrom === 'chrM') { mito++; skip(`line ${i + 1}: chrM`); continue; }
    if (!chrom || start == null || end == null || end < start) { skip(`line ${i + 1}: not a position (${l.slice(0, 40)})`); warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    let best: (typeof lrCols)[number] | null = null, lr: number | null = null;
    for (const x of lrCols) { const v = parseNum(f[x.i]); if (v != null && v !== 0 && (lr == null || Math.abs(v) > Math.abs(lr))) { lr = v; best = x; } }
    const kind = best?.kind ?? (lr != null ? (lr < 0 ? 'loss' : 'gain') : null);
    if (!kind) { skip(`line ${i + 1}: no log ratio in Amplification / Gain / Loss / Deletion`); warnings.push(`line ${i + 1}: no log ratio, not read`); continue; }
    const name = current ?? design ?? base;
    let s = samples.get(name);
    if (!s) { s = { key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs: [], sex: sexOf(metaOf(/^(gender|sex|sexe)$/i)) }; samples.set(name, s); }
    s.cnvs.push({
      id: `${id}:${i}`, sample: name, sampleKey: s.key, chrom, start, end,
      cytoband: cBand >= 0 ? (f[cBand] ?? '').trim() : '',
      genes: cGenes >= 0 ? (f[cGenes] ?? '').split(/[,;]/).map(g => g.trim()).filter(Boolean) : [],
      sizeKb: (end - start) / 1000, kind, typeText: best?.name ?? (kind === 'loss' ? 'Loss' : 'Gain'),
      probes: cProbes >= 0 ? parseNum(f[cProbes]) : null, logRatio: lr, iscn: null, copies: null, copiesHi: null, amp: false,
      suppressed: false, classification: '', classRest: '', note: '', line: i, pval: cP >= 0 ? parseNum(f[cP]) : null,
    });
  }
  for (const s of samples.values()) {
    s.declaredCalls = declared.get(s.name) ?? null;
    if (s.declaredCalls != null && s.declaredCalls !== s.cnvs.length) warnings.push(`${s.name}: the file announces ${s.declaredCalls} calls, ${s.cnvs.length} were read`);
  }
  if (mito) warnings.push(`${count(mito, 'call')} on chrM skipped`);
  if (!samples.size) throw new Error(`${fileName}: an interval-based report with no call that could be read.`);
  return { id, name: fileName, format: 'interval', build: cghBuild(buildText), buildText, meta, lines, eol, classCol: null, classSuffix: '', headerLine: hi, commentCol: null, samples: [...samples.values()], warnings, delim: '\t', spec: null, skipped, firstSkip };
}

/** A BED: chrom, start (0-based), end, then an optional name telling the type (del, dup, loss, gain, x1, x3…). */
function parseBed(lines: string[], id: string, fileName: string, base: string, eol: string): Omit<CghFile, 'saveAs'> {
  const s: CghSample = { key: `${id}:${base}`, name: base, fileId: id, declaredCalls: null, cnvs: [], sex: null };
  const warnings: string[] = [];
  let noType = 0, loh = 0, mito = 0, swapped = 0, wide = false;
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim() || /^(track|browser|#)/i.test(l)) continue;
    const f = l.split(/\t/).length >= 3 ? l.split('\t') : l.trim().split(/\s+/);
    const chrom = cghChrom(f[0] ?? '');
    let a = parsePos(f[1]), b = parsePos(f[2]);
    if (chrom === 'chrM') { mito++; continue; }
    if (a != null && b != null && b < a) { [a, b] = [b, a]; swapped++; }
    if (!chrom || a == null || b == null || b <= a) { warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    const name = (f[3] ?? '').trim();
    if (f.slice(3).some(x => x.trim())) wide = true;
    if (isLohText(name)) { loh++; continue; }
    // an ISCN as name, else x1 / x3 / x1~2 alone, else amp
    const bare = name.match(/(?:^|[^a-z])x(\d+)(?:~(\d+))?(?![a-z0-9])/i);
    const cn = name.includes('(') ? iscnCopies(name) : bare ? { copies: +bare[1], copiesHi: bare[2] ? +bare[2] : null, amp: false } : { copies: null, copiesHi: null, amp: /\bamp/i.test(name) };
    const kind = kindOfText(`${name} ${f.slice(4).join(' ')}`) ?? (cn.copies != null ? (cn.copies < 2 ? 'loss' : 'gain') : null);
    if (!kind) { noType++; continue; }
    s.cnvs.push({ id: `${id}:${i}`, sample: base, sampleKey: s.key, chrom, start: a + 1, end: b, cytoband: '', genes: [], sizeKb: (b - a) / 1000, kind, typeText: name || kind, probes: null, logRatio: null, iscn: null, ...cn, suppressed: false, classification: '', classRest: '', note: '', line: i });
  }
  if (noType) warnings.push(`${count(noType, 'interval')} without a type (the 4th column should say del / dup, loss / gain or x1 / x3), not read`);
  if (loh) warnings.push(`${count(loh, 'LOH / ROH interval')} skipped (copy-neutral)`);
  if (mito) warnings.push(`${count(mito, 'interval')} on chrM skipped`);
  if (swapped) warnings.push(`${count(swapped, 'interval')} with start after end: the two swapped`);
  if (!s.cnvs.length && noType && !wide) throw new UntypedBed(`${fileName}: a BED of three columns (chrom, start, end) says no type: add a 4th column saying del / dup, loss / gain or x1 / x3.`);
  if (!s.cnvs.length) throw new Error(`${fileName}: neither an interval table (no Chromosome / Start / Stop header) nor a BED with typed intervals.`);
  warnings.unshift('BED: no probe count, the minimum number of probes does not apply');
  return { id, name: fileName, format: 'bed', build: null, buildText: null, meta: [], lines, eol, classCol: null, classSuffix: '', headerLine: null, commentCol: null, samples: [s], warnings, delim: '\t', spec: null, skipped: noType, firstSkip: noType ? 'no type in the 4th column' : null };
}

/**
 * The file written back, each call with the class and the note given in the browser, the other lines as they were. An
 * interval table gets its Classification column filled and a Comment column (added after the last one when it has none);
 * a BED gets the class and the note as two last columns; a VCF gets them in INFO (CGH_CLASS, CGH_NOTE, percent-encoded).
 */
export function exportCghFile(file: CghFile, classOf: (cnv: CghCnv) => string, noteOf: (cnv: CghCnv) => string = c => c.note): string {
  const out = file.lines.slice(), d = file.delim, sp = (l: string) => splitDelimited(l, d);
  const flat = (t: string) => t.replace(/[\t\r\n]+/g, ' ').trim();
  const byLine = new Map(file.samples.flatMap(s => s.cnvs).map(c => [c.line, c]));
  if (file.format === 'vcf') return exportVcf(file, out, byLine, classOf, c => flat(noteOf(c)));
  // an interval table: its Classification and Comment columns, added after the last one when missing
  let classCol = file.classCol, noteCol = file.commentCol;
  if (file.headerLine != null && (classCol == null || noteCol == null)) {
    const h = sp(out[file.headerLine]);
    while (h.length && h[h.length - 1] === '') h.pop();
    if (classCol == null) { classCol = h.length; h.push('Classification'); }
    if (noteCol == null) { noteCol = h.length; h.push('Comment'); }
    out[file.headerLine] = joinDelimited(h, d);
  } else if (file.format === 'interval' && (classCol == null || noteCol == null)) {
    // a table without a header row: the missing columns after the widest call row
    const w = Math.max(...[...byLine.keys()].map(i => { const f = sp(out[i]); let n = f.length; while (n > 0 && !f[n - 1].trim()) n--; return n; }));
    if (classCol == null) classCol = noteCol === w ? w + 1 : w;
    if (noteCol == null) noteCol = Math.max(w, classCol + 1);
  }
  for (const [i, cnv] of byLine) {
    const cls = classOf(cnv).trim(), note = flat(noteOf(cnv));
    const cell = cls ? `${cls}${cnv.classRest || file.classSuffix}` : '';
    if (classCol != null && noteCol != null) {
      const f = sp(out[i]);
      while (f.length <= Math.max(classCol, noteCol)) f.push('');
      f[classCol] = cell; f[noteCol] = note; out[i] = joinDelimited(f, d);
    } else out[i] = `${out[i].replace(/\t+$/, '')}\t${cell}\t${note}`;
  }
  return out.join(file.eol);
}
/** INFO values: no blank, `;`, `=`, `,` or `%` as such */
export const vcfEncode = (t: string) => t.replace(/[%;=,\s]/g, ch => `%${ch.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);
function exportVcf(file: CghFile, out: string[], byLine: Map<number, CghCnv>, classOf: (cnv: CghCnv) => string, noteOf: (cnv: CghCnv) => string): string {
  for (const [i, cnv] of byLine) {
    const f = out[i].split('\t'), cls = classOf(cnv).trim(), note = noteOf(cnv);
    const info = (f[7] ?? '').split(';').filter(x => x && x !== '.' && !/^CGH_(CLASS|NOTE)=/.test(x));
    if (cls) info.push(`CGH_CLASS=${vcfEncode(cls)}`);
    if (note) info.push(`CGH_NOTE=${vcfEncode(note)}`);
    while (f.length < 8) f.push('.');
    f[7] = info.join(';') || '.';
    out[i] = f.join('\t');
  }
  const h = out.findIndex(l => l.startsWith('#CHROM'));
  const add = [['CGH_CLASS', 'Classification given in the array-CGH page'], ['CGH_NOTE', 'Note given in the array-CGH page']]
    .filter(([id]) => !out.some(l => l.startsWith(`##INFO=<ID=${id},`)))
    .map(([id, d]) => `##INFO=<ID=${id},Number=1,Type=String,Description="${d} (percent-encoded)">`);
  if (h >= 0) out.splice(h, 0, ...add);
  return out.join(file.eol);
}

/** The CNV-Hub result page of a CNV (Dijon University Hospital's platform: AnnotSV, ClassifyCNV, ISV, X-CNV and the
 * AChro-Puce criteria). The address is the one of its result pages, cnvhub.net/CNV/Result/hg19-chr16-24700000-36600000-loss;
 * only hg19 addresses have been seen, so GRCh38 gets none. */
export function cnvHubUrl(cnv: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>, build: GenomeBuild): string | null {
  if (build !== 'GRCh37') return null;
  return `https://cnvhub.net/CNV/Result/hg19-${cnv.chrom}-${cnv.start}-${cnv.end}-${cnv.kind}`;
}

/** The copy number as the ISCN writes it, for the table: x1, x0, x3, x1~2, amp, else the type. */
export const copyText = (c: CghCnv) => (c.copies != null ? `x${c.copies}${c.copiesHi != null ? `~${c.copiesHi}` : ''}` : c.amp ? 'amp' : c.kind === 'loss' ? 'loss' : 'gain');

// ---- classes ----
/** The classes, as the laboratory writes them (and as they are exported). */
export const CGH_CLASSES: { value: string; color: string; hint?: string }[] = [
  { value: 'Benign', color: '#15803d' },
  { value: 'Likely benign', color: '#22a35a' },
  { value: 'Uncertain significance', color: '#64748b' },
  { value: 'Likely pathogenic', color: '#ea580c' },
  { value: 'Pathogenic', color: '#b91c1c' },
  { value: 'PIEV', color: '#7c3aed', hint: 'susceptibility: incomplete penetrance, variable expressivity' },
  { value: 'Male Reference', color: '#0284c7', hint: 'a call of the male reference' },
  { value: 'Female Reference', color: '#db2777', hint: 'a call of the female reference' },
  { value: 'NS', color: '#94a3b8', hint: 'non significant: recurring technical artifact / atypical log ratio' },
];
/** a class a file writes in other words (French or English, with or without its number), as one of the classes; null if none */
export function toClass(v: string): string | null {
  const t = norm(v);
  if (!t) return null;
  const known = CGH_CLASSES.find(c => norm(c.value) === t);
  if (known) return known.value;
  if (/probablementpathog|likelypathog|classe?4/.test(t)) return 'Likely pathogenic';
  if (/pathog|classe?5/.test(t)) return 'Pathogenic';
  if (/probablementbenin|likelybenign|classe?2/.test(t)) return 'Likely benign';
  if (/benin|benign|classe?1/.test(t)) return 'Benign';
  if (/vsi|vus|incertain|uncertain|inconnu|classe?3/.test(t)) return 'Uncertain significance';
  if (/piev|susceptib/.test(t)) return 'PIEV';
  if (/temoinfeminin|referencefeminine?|femalereference|femalecontrol/.test(t)) return 'Female Reference';
  if (/temoinmasculin|referencemasculine?|malereference|malecontrol/.test(t)) return 'Male Reference';
  if (/^ns$|nonsignificati|nonsignificant|artefact|artifact/.test(t)) return 'NS';
  return null;
}
export const classColor = (v: string): string | null => (v.trim() ? CGH_CLASSES.find(c => c.value === v)?.color ?? '#94a3b8' : null);
