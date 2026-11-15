/**
 * Array-CGH calls for the SV browser: the interval tables exported per sample by Agilent CytoGenomics (a few header
 * lines, `Genome hg19`, `Record Name …`, then one row per call: Chromosome, Start, Stop, Cytoband, Gene Name, Size(kb),
 * Type, #Probes, Mean Log Ratio, ISCN, Suppress, Classification), or a BED file (chrom, start, end, name).
 *
 * Columns are found by their header, not by their place. Positions of the interval table are read as the ISCN of the
 * same row writes them (1-based, the first and last probe of the call); a BED's start is 0-based and moved to 1-based.
 * The copy number comes from the ISCN (x0, x1, x3, x4…; `amp` for a high-level gain), the type from the Type column
 * (else from the copy number, else from the log ratio's sign). A file is never sent anywhere: the calls stay in this
 * tab, with the classification typed in the table, and can be written back as the same table.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { norm, parseNum, splitDelimited, joinDelimited, detectTable, autoCols, AUTO_OPTS, SIZE_MULT, NeedsMapping, type MapSpec } from './cghMapping';

export type CghKind = 'loss' | 'gain';
export type CghSex = 'male' | 'female';
/** a sex as a header line may write it (Gender / Sex: Male, F, Femme…) */
export const sexOf = (v: string | null | undefined): CghSex | null => (!v ? null : /^\s*(m|male|homme|masculin)\b/i.test(v) ? 'male' : /^\s*(f|female|femme|f[ée]minin)\b/i.test(v) ? 'female' : null);
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
}
export interface CghSample { key: string; name: string; fileId: string; declaredCalls: number | null; cnvs: CghCnv[];
  /** the sex a header line gives (Gender, Sex), else null */
  sex: CghSex | null }
export interface CghFile {
  /** the file's name (for a file of a zip: zip › path inside), and the name its export is saved under */
  id: string; name: string; saveAs: string; format: 'interval' | 'bed';
  /** the build the file states (Genome hg19 → GRCh37), null when it states none */
  build: GenomeBuild | null; buildText: string | null;
  /** the header lines as key / value (Array Design Info, Creation Date…) */
  meta: [string, string][];
  lines: string[]; eol: string;
  /** the Classification column (interval table), or null (a BED: the class is added as a last column) */
  classCol: number | null;
  /** whether the file's classifications end with `;;` (written the same way back) */
  classSuffix: string;
  /** the header line (interval table) and its Comment / Note column, if any (else one is added at the export) */
  headerLine: number | null; commentCol: number | null;
  samples: CghSample[];
  warnings: string[];
  /** the delimiter (tab for the software's own export), and the columns this file was read with (null for a BED) */
  delim: string; spec: MapSpec | null;
  /** the zip's patient name, to read the file again with other columns */
  fallback?: string;
  /** rows not read (unreadable, or a Type value set to ignore) and the first reason; whether a remembered mapping was used */
  skipped: number; firstSkip: string | null; remembered?: boolean;
}

const num = (s: string | undefined): number | null => parseNum(s);
/** chr1, 1, chrX, X → chr1, chrX; null for anything else */
export function cghChrom(s: string): string | null {
  const m = s.trim().match(/^(?:chr)?(\d{1,2}|X|Y|M|MT)$/i);
  if (!m) return null;
  const c = m[1].toUpperCase();
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
export function iscnCopies(iscn: string | null): { copies: number | null; copiesHi: number | null; amp: boolean } {
  const m = iscn?.match(/\)\s*x\s*(\d+)(?:\s*~\s*(\d+))?/i);
  if (m) return { copies: +m[1], copiesHi: m[2] ? +m[2] : null, amp: false };
  return { copies: null, copiesHi: null, amp: !!iscn && /\)\s*amp\b/i.test(iscn) };
}
export const kindOfText = (raw: string): CghKind | null => {
  const t = raw.normalize('NFD').replace(/[̀-ͯ]/g, '');
  return /del|loss|perte|\bx[01]\b/i.test(t) ? 'loss' : /amp|gain|dup|\bx[3-9]\b/i.test(t) ? 'gain' : null;
};
/** the classification cell: the class before `;;`, and the rest */
function splitClass(cell: string): { classification: string; classRest: string } {
  const i = cell.indexOf(';;');
  return i < 0 ? { classification: cell.trim(), classRest: '' } : { classification: cell.slice(0, i).trim(), classRest: cell.slice(i) };
}

let fileSeq = 0;
/**
 * One file: an interval table (CytoGenomics) or a BED. Throws when no call can be read from it. The sample is the
 * table's Record Name, else `fallbackSample` (a zip's name: one patient per zip), else the file's name.
 */
export function parseCghFile(text: string, fileName: string, fallbackSample?: string, spec?: MapSpec): CghFile {
  // the commoner line ending (a file may mix them), written back the same
  const crlf = (text.match(/\r\n/g) ?? []).length, eol = crlf * 2 > (text.match(/\n/g) ?? []).length ? '\r\n' : '\n';
  const lines = text.replace(/^﻿/, '').split(/\r?\n/);
  const id = `cgh${++fileSeq}`;
  const leaf = fileName.split(/[›/\\]/).pop()!.trim();
  const base = fallbackSample ?? leaf.replace(/\.(txt|tsv|bed|xls)$/i, '');
  let f: Omit<CghFile, 'saveAs'>;
  if (spec) f = parseInterval(lines, id, fileName, base, eol, spec);
  else {
    const hi = lines.findIndex(l => { const c = l.split('\t').map(norm); return c.includes('chromosome') && c.includes('start') && (c.includes('stop') || c.includes('end')); });
    try {
      f = hi >= 0 ? parseInterval(lines, id, fileName, base, eol, { cols: autoCols(lines[hi].split('\t').map(norm)), opts: AUTO_OPTS(hi) }) : parseBed(lines, id, fileName, base, eol);
    } catch (e) {
      // a delimited table whose columns are not recognised: the caller asks which column is which
      const t = detectTable(lines);
      if (t && t.ncols >= 3) throw new NeedsMapping((e as Error).message);
      throw e;
    }
  }
  return { ...f, fallback: fallbackSample, saveAs: leaf.replace(/(\.[^.]+)?$/, m => `.classified${m || '.txt'}`) };
}

function parseInterval(lines: string[], id: string, fileName: string, base: string, eol: string, spec: MapSpec): Omit<CghFile, 'saveAs'> {
  const { cols: c, opts } = spec, hi = opts.headerRow, d = opts.delim;
  const sp = (l: string) => (d === '\t' ? l.split('\t') : splitDelimited(l, d));
  // a decimal comma only where the delimiter is not the comma (a log ratio, a size, a copy number)
  const dec = d !== ',', nn = (s: string | undefined) => parseNum(s, dec);
  const meta: [string, string][] = [];
  for (const l of lines.slice(0, hi)) {
    const [k, ...v] = sp(l);
    if (k?.trim()) meta.push([k.trim().replace(/:$/, ''), v.join(' ').trim()]);
  }
  const metaOf = (re: RegExp) => meta.find(([k]) => re.test(k))?.[1] || null;
  const buildText = metaOf(/^genome$/i);
  const record = metaOf(/^record name$/i);
  const warnings: string[] = [];
  let skipped = 0, firstSkip: string | null = null, ignored = 0;
  const skip = (why: string) => { skipped++; firstSkip ??= why; };
  if (c.probes < 0) warnings.push('no #Probes column: the minimum number of probes cannot be applied');
  const samples = new Map<string, CghSample>();
  const declared = new Map<string, number>();
  let suffixYes = 0, suffixNo = 0;
  for (let i = hi + 1; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim()) continue;
    const foot = l.match(/^\s*([^\t=]+?)\s*=\s*(\d+)\s*$/);
    if (foot) { declared.set(foot[1], +foot[2]); continue; }
    if (/^number of calls/i.test(l)) continue;
    const f = sp(l);
    const chrom = c.chrom >= 0 ? cghChrom(f[c.chrom] ?? '') : null;
    const s0 = num(f[c.start]), end = num(f[c.stop]);
    const start = s0 != null && opts.zeroBased ? s0 + 1 : s0;
    if (!chrom || start == null || end == null || end < start) { skip(`line ${i + 1}: ${!chrom ? 'chromosome' : 'start / end'} not a position (${l.slice(0, 40)})`); warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    const iscnRaw = c.iscn >= 0 ? (f[c.iscn] ?? '').trim() : '';
    const iscn = iscnRaw && !/^na$/i.test(iscnRaw) ? iscnRaw : null;
    let cn = iscnCopies(iscn);
    if (cn.copies == null && !cn.amp && c.copies >= 0) { const n = nn(f[c.copies]); if (n != null) cn = { copies: Math.round(n), copiesHi: null, amp: false }; }
    const typeText = c.type >= 0 ? (f[c.type] ?? '').trim() : '';
    const lr = c.lr >= 0 ? nn(f[c.lr]) : null;
    const choice = typeText ? opts.typeMap[typeText.toLowerCase()] : undefined;
    if (choice === 'ignore') { ignored++; skip(`line ${i + 1}: type "${typeText}" ignored`); continue; }
    const kind = choice ?? kindOfText(typeText) ?? (cn.amp ? 'gain' : cn.copies != null ? (cn.copies < 2 ? 'loss' : 'gain') : lr != null ? (lr < 0 ? 'loss' : 'gain') : null);
    if (!kind) { skip(`line ${i + 1}: no type, copy number or log ratio`); warnings.push(`line ${i + 1}: no type, copy number or log ratio, not read`); continue; }
    const cell = c.cls >= 0 ? (f[c.cls] ?? '') : '';
    if (cell.trim()) { if (cell.includes(';;')) suffixYes++; else suffixNo++; }
    const name = (c.sample >= 0 && f[c.sample]?.trim()) || record || base;
    let s = samples.get(name);
    if (!s) { s = { key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs: [], sex: sexOf(metaOf(/^(gender|sex|sexe)$/i)) }; samples.set(name, s); }
    const size = c.size >= 0 ? nn(f[c.size]) : null;
    s.cnvs.push({
      id: `${id}:${i}`, sample: name, sampleKey: s.key, chrom, start, end,
      cytoband: c.band >= 0 ? (f[c.band] ?? '').trim() : '',
      genes: c.genes >= 0 ? (f[c.genes] ?? '').split(/[,;]/).map(g => g.trim()).filter(Boolean) : [],
      sizeKb: size != null ? size * SIZE_MULT[opts.sizeUnit] : (end - start) / 1000, kind, typeText: typeText || (kind === 'loss' ? 'Deletion' : 'Gain'),
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
  if (ignored) warnings.push(`${ignored} row${ignored > 1 ? 's' : ''} ignored by their type`);
  if (!samples.size) throw new Error(`${fileName}: no call could be read (an interval table needs Chromosome, Start, Stop and a Type, an ISCN or a log ratio).`);
  return { id, name: fileName, format: 'interval', build: cghBuild(buildText) ?? opts.build, buildText, meta, lines, eol, classCol: c.cls >= 0 ? c.cls : null, classSuffix: suffixYes > suffixNo ? ';;' : '', headerLine: hi, commentCol: c.note >= 0 ? c.note : null, samples: [...samples.values()], warnings, delim: d, spec, skipped, firstSkip };
}

/** A BED: chrom, start (0-based), end, then an optional name telling the type (del, dup, loss, gain, x1, x3…). */
function parseBed(lines: string[], id: string, fileName: string, base: string, eol: string): Omit<CghFile, 'saveAs'> {
  const s: CghSample = { key: `${id}:${base}`, name: base, fileId: id, declaredCalls: null, cnvs: [], sex: null };
  const warnings: string[] = [];
  let noType = 0;
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim() || /^(track|browser|#)/i.test(l)) continue;
    const f = l.split(/\t/).length >= 3 ? l.split('\t') : l.trim().split(/\s+/);
    const chrom = cghChrom(f[0] ?? ''), a = num(f[1]), b = num(f[2]);
    if (!chrom || a == null || b == null || b <= a) { warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    const name = (f[3] ?? '').trim();
    // an ISCN as name, else x1 / x3 / x1~2 alone, else amp
    const bare = name.match(/(?:^|[^a-z])x(\d+)(?:~(\d+))?(?![a-z0-9])/i);
    const cn = name.includes('(') ? iscnCopies(name) : bare ? { copies: +bare[1], copiesHi: bare[2] ? +bare[2] : null, amp: false } : { copies: null, copiesHi: null, amp: /\bamp/i.test(name) };
    const kind = kindOfText(`${name} ${f.slice(4).join(' ')}`) ?? (cn.copies != null ? (cn.copies < 2 ? 'loss' : 'gain') : null);
    if (!kind) { noType++; continue; }
    s.cnvs.push({ id: `${id}:${i}`, sample: base, sampleKey: s.key, chrom, start: a + 1, end: b, cytoband: '', genes: [], sizeKb: (b - a) / 1000, kind, typeText: name || kind, probes: null, logRatio: null, iscn: null, ...cn, suppressed: false, classification: '', classRest: '', note: '', line: i });
  }
  if (noType) warnings.push(`${noType} interval${noType > 1 ? 's' : ''} without a type (the 4th column should say del / dup, loss / gain or x1 / x3), not read`);
  if (!s.cnvs.length) throw new Error(`${fileName}: neither an interval table (no Chromosome / Start / Stop header) nor a BED with typed intervals.`);
  warnings.unshift('BED: no probe count, the minimum number of probes does not apply');
  return { id, name: fileName, format: 'bed', build: null, buildText: null, meta: [], lines, eol, classCol: null, classSuffix: '', headerLine: null, commentCol: null, samples: [s], warnings, delim: '\t', spec: null, skipped: noType, firstSkip: noType ? 'no type in the 4th column' : null };
}

/**
 * The file written back, each call with the class and the note given in the browser, the other lines as they were. An
 * interval table gets its Classification column filled and a Comment column (added after the last one when it has none);
 * a BED gets the class and the note as two last columns.
 */
export function exportCghFile(file: CghFile, classOf: (cnv: CghCnv) => string, noteOf: (cnv: CghCnv) => string = c => c.note): string {
  const out = file.lines.slice(), d = file.delim, sp = (l: string) => (d === '\t' ? l.split('\t') : splitDelimited(l, d));
  const flat = (t: string) => t.replace(/[\t\r\n]+/g, ' ').trim();
  const byLine = new Map(file.samples.flatMap(s => s.cnvs).map(c => [c.line, c]));
  // an interval table: its Classification and Comment columns, added after the last one when missing
  let classCol = file.classCol, noteCol = file.commentCol;
  if (file.headerLine != null && (classCol == null || noteCol == null)) {
    const h = sp(out[file.headerLine]);
    while (h.length && h[h.length - 1] === '') h.pop();
    if (classCol == null) { classCol = h.length; h.push('Classification'); }
    if (noteCol == null) { noteCol = h.length; h.push('Comment'); }
    out[file.headerLine] = joinDelimited(h, d);
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

/** The CNV-Hub result page of a CNV (Dijon University Hospital's platform: AnnotSV, ClassifyCNV, ISV, X-CNV and the
 * AChro-Puce criteria). The address is the one of its result pages, cnvhub.net/CNV/Result/hg19-chr16-24700000-36600000-loss;
 * only hg19 addresses have been seen, so GRCh38 gets none. */
export function cnvHubUrl(cnv: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>, build: GenomeBuild): string | null {
  if (build !== 'GRCh37') return null;
  return `https://cnvhub.net/CNV/Result/hg19-${cnv.chrom}-${cnv.start}-${cnv.end}-${cnv.kind}`;
}

/** The text the SV browser's box takes for a CNV: DEL / DUP before its coordinates. */
export const cnvSearchText = (cnv: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>) => `${cnv.kind === 'loss' ? 'DEL' : 'DUP'} ${cnv.chrom}:${cnv.start}-${cnv.end}`;

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
