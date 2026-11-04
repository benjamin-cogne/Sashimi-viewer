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

export type CghKind = 'loss' | 'gain';
export interface CghCnv {
  /** unique within the tab: file, row */
  id: string;
  sample: string;
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
  /** the file's line of this call (0-based), for the export */
  line: number;
}
export interface CghSample { key: string; name: string; fileId: string; declaredCalls: number | null; cnvs: CghCnv[] }
export interface CghFile {
  id: string; name: string; format: 'interval' | 'bed';
  /** the build the file states (Genome hg19 → GRCh37), null when it states none */
  build: GenomeBuild | null; buildText: string | null;
  /** the header lines as key / value (Array Design Info, Creation Date…) */
  meta: [string, string][];
  lines: string[]; eol: string;
  /** the Classification column (interval table), or null (a BED: the class is added as a last column) */
  classCol: number | null;
  /** whether the file's classifications end with `;;` (written the same way back) */
  classSuffix: string;
  samples: CghSample[];
  warnings: string[];
}

const norm = (s: string) => s.toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '').replace(/[^a-z0-9]/g, '');
const num = (s: string | undefined): number | null => {
  const t = (s ?? '').trim().replace(/,/g, '');
  if (!t || /^(na|nan|-)$/i.test(t)) return null;
  const n = Number(t);
  return Number.isFinite(n) ? n : null;
};
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
const kindOfText = (t: string): CghKind | null => (/del|loss|perte|\bx[01]\b/i.test(t) ? 'loss' : /amp|gain|dup|\bx[3-9]\b/i.test(t) ? 'gain' : null);
/** the classification cell: the class before `;;`, and the rest */
function splitClass(cell: string): { classification: string; classRest: string } {
  const i = cell.indexOf(';;');
  return i < 0 ? { classification: cell.trim(), classRest: '' } : { classification: cell.slice(0, i).trim(), classRest: cell.slice(i) };
}

let fileSeq = 0;
/** One file: an interval table (CytoGenomics) or a BED. Throws when no call can be read from it. */
export function parseCghFile(text: string, fileName: string): CghFile {
  // the commoner line ending (a file may mix them), written back the same
  const crlf = (text.match(/\r\n/g) ?? []).length, eol = crlf * 2 > (text.match(/\n/g) ?? []).length ? '\r\n' : '\n';
  const lines = text.replace(/^﻿/, '').split(/\r?\n/);
  const id = `cgh${++fileSeq}`;
  const base = fileName.replace(/\.(txt|tsv|bed|xls)$/i, '');
  const hi = lines.findIndex(l => { const c = l.split('\t').map(norm); return c.includes('chromosome') && c.includes('start') && (c.includes('stop') || c.includes('end')); });
  return hi >= 0 ? parseInterval(lines, hi, id, fileName, base, eol) : parseBed(lines, id, fileName, base, eol);
}

function parseInterval(lines: string[], hi: number, id: string, fileName: string, base: string, eol: string): CghFile {
  const meta: [string, string][] = [];
  for (const l of lines.slice(0, hi)) {
    const [k, ...v] = l.split('\t');
    if (k?.trim()) meta.push([k.trim().replace(/:$/, ''), v.join(' ').trim()]);
  }
  const metaOf = (re: RegExp) => meta.find(([k]) => re.test(k))?.[1] || null;
  const buildText = metaOf(/^genome$/i);
  const record = metaOf(/^record name$/i);
  const head = lines[hi].split('\t').map(norm);
  const col = (...names: string[]) => { for (const n of names) { const i = head.indexOf(n); if (i >= 0) return i; } return -1; };
  const c = {
    chrom: col('chromosome', 'chr', 'chrom'), start: col('start'), stop: col('stop', 'end'), band: col('cytoband', 'band'),
    genes: col('genename', 'genenames', 'genes', 'gene'), size: col('sizekb', 'size'), type: col('type', 'aberrationtype'),
    probes: col('probes', 'probe', 'nprobes', 'numprobes', 'numberofprobes'), lr: col('meanlogratio', 'logratio', 'log2ratio'),
    iscn: col('iscn'), sup: col('suppress', 'suppressed'), cls: col('classification', 'class'),
    sample: col('samplename', 'sample', 'recordname', 'arrayname'),
  };
  const warnings: string[] = [];
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
    const f = l.split('\t');
    const chrom = cghChrom(f[c.chrom] ?? '');
    const start = num(f[c.start]), end = num(f[c.stop]);
    if (!chrom || start == null || end == null || end < start) { warnings.push(`line ${i + 1} not read: ${l.slice(0, 60)}`); continue; }
    const iscnRaw = c.iscn >= 0 ? (f[c.iscn] ?? '').trim() : '';
    const iscn = iscnRaw && !/^na$/i.test(iscnRaw) ? iscnRaw : null;
    const cn = iscnCopies(iscn);
    const typeText = c.type >= 0 ? (f[c.type] ?? '').trim() : '';
    const lr = c.lr >= 0 ? num(f[c.lr]) : null;
    const kind = kindOfText(typeText) ?? (cn.amp ? 'gain' : cn.copies != null ? (cn.copies < 2 ? 'loss' : 'gain') : lr != null ? (lr < 0 ? 'loss' : 'gain') : null);
    if (!kind) { warnings.push(`line ${i + 1}: no type, copy number or log ratio, not read`); continue; }
    const cell = c.cls >= 0 ? (f[c.cls] ?? '') : '';
    if (cell.trim()) { if (cell.includes(';;')) suffixYes++; else suffixNo++; }
    const name = (c.sample >= 0 && f[c.sample]?.trim()) || record || base;
    let s = samples.get(name);
    if (!s) { s = { key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs: [] }; samples.set(name, s); }
    const size = c.size >= 0 ? num(f[c.size]) : null;
    s.cnvs.push({
      id: `${id}:${i}`, sample: name, chrom, start, end,
      cytoband: c.band >= 0 ? (f[c.band] ?? '').trim() : '',
      genes: c.genes >= 0 ? (f[c.genes] ?? '').split(/[,;]/).map(g => g.trim()).filter(Boolean) : [],
      sizeKb: size ?? (end - start) / 1000, kind, typeText: typeText || (kind === 'loss' ? 'Deletion' : 'Gain'),
      probes: c.probes >= 0 ? num(f[c.probes]) : null, logRatio: lr, iscn, ...cn,
      suppressed: c.sup >= 0 && /^(true|yes|1)$/i.test((f[c.sup] ?? '').trim()),
      ...splitClass(cell), line: i,
    });
  }
  for (const s of samples.values()) {
    s.declaredCalls = declared.get(s.name) ?? (samples.size === 1 && declared.size === 1 ? [...declared.values()][0] : null);
    // the footer counts the calls the software made; a different count is said, not corrected
    if (s.declaredCalls != null && s.declaredCalls !== s.cnvs.length) warnings.push(`${s.name}: the file announces ${s.declaredCalls} calls, ${s.cnvs.length} were read`);
  }
  if (!samples.size) throw new Error(`${fileName}: no call could be read (an interval table needs Chromosome, Start, Stop and a Type, an ISCN or a log ratio).`);
  return { id, name: fileName, format: 'interval', build: cghBuild(buildText), buildText, meta, lines, eol, classCol: c.cls >= 0 ? c.cls : null, classSuffix: suffixYes > suffixNo ? ';;' : '', samples: [...samples.values()], warnings };
}

/** A BED: chrom, start (0-based), end, then an optional name telling the type (del, dup, loss, gain, x1, x3…). */
function parseBed(lines: string[], id: string, fileName: string, base: string, eol: string): CghFile {
  const s: CghSample = { key: `${id}:${base}`, name: base, fileId: id, declaredCalls: null, cnvs: [] };
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
    s.cnvs.push({ id: `${id}:${i}`, sample: base, chrom, start: a + 1, end: b, cytoband: '', genes: [], sizeKb: (b - a) / 1000, kind, typeText: name || kind, probes: null, logRatio: null, iscn: null, ...cn, suppressed: false, classification: '', classRest: '', line: i });
  }
  if (noType) warnings.push(`${noType} interval${noType > 1 ? 's' : ''} without a type (the 4th column should say del / dup, loss / gain or x1 / x3), not read`);
  if (!s.cnvs.length) throw new Error(`${fileName}: neither an interval table (no Chromosome / Start / Stop header) nor a BED with typed intervals.`);
  warnings.unshift('BED: no probe count, the minimum number of probes does not apply');
  return { id, name: fileName, format: 'bed', build: null, buildText: null, meta: [], lines, eol, classCol: null, classSuffix: '', samples: [s], warnings };
}

/** The file written back, each call with the class given in the browser (the other lines as they were). */
export function exportCghFile(file: CghFile, classOf: (cnv: CghCnv) => string): string {
  const out = file.lines.slice();
  const byLine = new Map(file.samples.flatMap(s => s.cnvs).map(c => [c.line, c]));
  for (const [i, cnv] of byLine) {
    const cls = classOf(cnv).trim();
    const cell = cls ? `${cls}${cnv.classRest || file.classSuffix}` : '';
    const f = out[i].split('\t');
    if (file.classCol != null) { while (f.length <= file.classCol) f.push(''); f[file.classCol] = cell; out[i] = f.join('\t'); }
    else out[i] = `${out[i].replace(/\t+$/, '')}\t${cell}`;
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
/** The classes offered: the AChro-Puce classes (the French network's, 5 to 1 and PIEV), then the files' own words. */
export const CGH_CLASSES: { value: string; color: string }[] = [
  { value: 'Pathogenic (class 5)', color: '#b91c1c' },
  { value: 'Likely pathogenic (class 4)', color: '#ea580c' },
  { value: 'PIEV (susceptibility)', color: '#7c3aed' },
  { value: 'Uncertain significance (class 3)', color: '#64748b' },
  { value: 'Likely benign (class 2)', color: '#22a35a' },
  { value: 'Benign (class 1)', color: '#15803d' },
];
/** a class's colour: one of the offered classes, else read from its words (a file's own vocabulary, French or English) */
export function classColor(v: string): string | null {
  if (!v.trim()) return null;
  const known = CGH_CLASSES.find(c => c.value === v);
  if (known) return known.color;
  const t = norm(v);
  if (/probablementpathog|likelypathog|class4|classe4/.test(t)) return '#ea580c';
  if (/pathog|class5|classe5/.test(t)) return '#b91c1c';
  if (/piev|susceptib/.test(t)) return '#7c3aed';
  if (/probablementbenin|likelybenign|class2|classe2/.test(t)) return '#22a35a';
  if (/benin|benign|class1|classe1/.test(t)) return '#15803d';
  if (/vsi|vus|incertain|uncertain|inconnu|class3|classe3/.test(t)) return '#64748b';
  return '#94a3b8';
}
