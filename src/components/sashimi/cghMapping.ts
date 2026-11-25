/**
 * Column mapping of an array-CGH table: the delimited text (tab, `;` or `,`, quoted fields) is split, its header row
 * and delimiter guessed, its columns matched to the fields of a call by name (English and French variants), and what the
 * user chose is remembered by the header's column names (and the Type words given loss or gain), never a patient's row.
 */
import type { GenomeBuild } from '../../standalone/ensembl';

export const norm = (s: string) => s.toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '').replace(/[^a-z0-9]/g, '');

export type Field = 'chrom' | 'start' | 'stop' | 'pos' | 'band' | 'genes' | 'size' | 'type' | 'copies' | 'probes' | 'lr' | 'iscn' | 'sup' | 'cls' | 'note' | 'sample';
/** the column index of each field, -1 when the table has none */
export type ColMap = Record<Field, number>;
const FIELDS: Field[] = ['chrom', 'start', 'stop', 'pos', 'type', 'copies', 'iscn', 'lr', 'probes', 'size', 'band', 'genes', 'cls', 'note', 'sample', 'sup'];
const noCols = (): ColMap => Object.fromEntries(FIELDS.map(f => [f, -1])) as ColMap;
export type TypeChoice = 'loss' | 'gain' | 'ignore';
export interface MapOpts {
  delim: string;
  /** line of the header (0-based, in the file's lines); -1: no header row, every line is a call */
  headerRow: number;
  /** BED-like start, moved to 1-based */
  zeroBased: boolean;
  sizeUnit: 'bp' | 'kb' | 'Mb';
  /** the build to take when the file states none */
  build: GenomeBuild | null;
  /** what a Type value that does not say loss or gain means (lower case): only the values given loss or gain are kept */
  typeMap: Record<string, TypeChoice>;
  /** a Type value neither recognised nor in typeMap is ignored (the mapping dialog's default), else read from the copy number / log ratio */
  ignoreOtherTypes?: boolean;
}
/** `lrMore`: further log-ratio columns (gain / loss ones, as Amplification / Gain / Loss / Deletion), merged with `lr`: the
 * non-zero value; absent in mappings saved before it */
export interface MapSpec { cols: ColMap; opts: MapOpts; lrMore?: number[] }
/** thrown by the automatic reading when a table looks delimited but its columns are not recognised */
export class NeedsMapping extends Error {}

/**
 * The text's records: lines ending in CR LF, LF or CR alone, a quoted field ("…" opened at a field's start) holding
 * line breaks kept in its record. A quote left open to the end is not one: the plain lines are returned.
 */
export function splitRecords(text: string): string[] {
  const t = text.replace(/^\uFEFF/, '');
  if (!t.includes('"')) return t.split(/\r\n|\n|\r/);
  const out: string[] = [];
  let start = 0, q = false;
  for (let i = 0; i < t.length; i++) {
    const ch = t[i];
    if (q) { if (ch === '"') { if (t[i + 1] === '"') i++; else q = false; } continue; }
    if (ch === '"') { const p = t[i - 1]; if (i === start || p === '\t' || p === ';' || p === ',') q = true; }
    else if (ch === '\n' || ch === '\r') { out.push(t.slice(start, i)); if (ch === '\r' && t[i + 1] === '\n') i++; start = i + 1; }
  }
  if (q) return t.split(/\r\n|\n|\r/);
  out.push(t.slice(start));
  return out;
}
/** the commoner line ending of a text (a file may mix them), written back the same */
export function eolOf(text: string): string {
  const crlf = (text.match(/\r\n/g) ?? []).length, lf = (text.match(/\n/g) ?? []).length - crlf, cr = (text.match(/\r/g) ?? []).length - crlf;
  return cr > crlf && cr > lf ? '\r' : crlf > lf ? '\r\n' : '\n';
}

/** a line split at a one-character delimiter, "quoted, fields" and "" kept; a quote only opens at a field's start */
export function splitDelimited(line: string, d: string): string[] {
  if (!line.includes('"')) return line.split(d);
  const out: string[] = []; let cur = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (q) { if (ch === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += ch; }
    else if (ch === '"' && cur === '') q = true;
    else if (ch === d) { out.push(cur); cur = ''; }
    else cur += ch;
  }
  out.push(cur);
  return out;
}
export const joinDelimited = (cells: string[], d: string) => cells.map(c => (d === '\t' ? /[\t\r\n]/.test(c) : c.includes(d) || c.includes('"') || /[\r\n]/.test(c)) ? `"${c.replace(/"/g, '""')}"` : c).join(d);

/** a number, with a decimal comma read as one when `decimal` (else commas are thousands separators) */
export function parseNum(s: string | undefined, decimal = false): number | null {
  let t = (s ?? '').trim();
  if (decimal && /^[-+]?\d+,\d+$/.test(t)) t = t.replace(',', '.');
  t = t.replace(/,/g, '');
  if (!t || /^(na|nan|-)$/i.test(t)) return null;
  const n = Number(t);
  return Number.isFinite(n) ? n : null;
}
const isNum = (s: string) => /^[-+]?\d+(?:[.,]\d+)*$/.test(s.trim());
/** a position: digits, or digits grouped by threes with one separator (, . space ' ’), else a plain number */
export function parsePos(s: string | undefined): number | null {
  const t = (s ?? '').trim();
  if (/^\d+$/.test(t)) return Number(t);
  if (/^\d{1,3}([,. '’\u00a0\u202f])\d{3}(?:\1\d{3})*$/.test(t)) return Number(t.replace(/\D/g, ''));
  return parseNum(t);
}
const POS = String.raw`\d[\d,. '’\u00a0\u202f]*?`;
const REGION = new RegExp(String.raw`^(?:chr)?(\d{1,2}|X|Y|MT?)\s*:\s*(?:g\.)?(${POS})\s*(?:-|–|\.\.|_)\s*(${POS})\s*$`, 'i');
const ISCN_REGION = /(?:^|[\s\]])(\d{1,2}|X|Y)(?:[pq][\d.]*){0,2}\s*\(\s*([\d,]+)\s*[_-]\s*([\d,]+)\s*\)/i;
/** a region in one cell: chr1:72523846-72667555, 1:72,523,846-72,667,555, chr1:72523846..72667555, or an ISCN's
 * 1p31.1(72523846_72667555); the chromosome as written (1, X, MT), positions as written (1-based for an ISCN) */
export function parseRegion(s: string | undefined): { chrom: string; start: number; end: number; iscn: boolean } | null {
  const t = (s ?? '').trim();
  const m = t.match(REGION), i = m ? null : t.match(ISCN_REGION);
  const r = m ?? i;
  if (!r) return null;
  const a = parsePos(r[2]), b = parsePos(r[3]);
  return a == null || b == null ? null : { chrom: r[1], start: a, end: b, iscn: !!i };
}
/** a copy number as a cell writes it: x1, 1, 1.0, CN=1, CN:1, 1 copy, 3 copies; null for anything else */
export function copiesOf(s: string | undefined): number | null {
  const m = (s ?? '').trim().match(/^(?:x|cn\s*[=:]?\s*)?(\d+(?:[.,]\d+)?)\s*(?:cop(?:y|ies))?$/i);
  return m ? Math.round(Number(m[1].replace(',', '.'))) : null;
}
/** a row's cells without the empty ones a trailing delimiter leaves */
const width = (l: string, d: string) => { const c = splitDelimited(l, d); let n = c.length; while (n > 0 && !c[n - 1].trim()) n--; return n; };
/** how much a row looks like a call: numbers, and a region in one cell counting for two */
const dataScore = (cells: string[]) => cells.reduce((k, c) => k + (parseRegion(c) ? 2 : isNum(c) ? 1 : 0), 0);
/** a header row: most of its non-empty cells neither numbers nor positions (a data row must never be remembered) */
export function looksLikeHeader(cells: string[]): boolean {
  const t = cells.map(c => c.trim()).filter(Boolean);
  return t.length > 0 && t.filter(c => !isNum(c) && !parseRegion(c) && parsePos(c) == null).length * 2 > t.length;
}

/**
 * The delimiter and the header row of a delimited table (≥ 2 columns, a trailing delimiter not counted; 2 only with a
 * position in one cell), headerRow -1 when it has none (its rows are calls); null when the text does not look like one.
 */
export function detectTable(lines: string[]): { delim: string; headerRow: number; ncols: number } | null {
  const rows = lines.map((l, i) => ({ l, i })).filter(r => r.l.trim()).slice(0, 40);
  let best: { delim: string; mode: number; n: number } | null = null;
  for (const delim of ['\t', ';', ',']) {
    const freq = new Map<number, number>();
    for (const r of rows) { const k = width(r.l, delim); if (k >= 3 || (k === 2 && splitDelimited(r.l, delim).some(c => parseRegion(c)))) freq.set(k, (freq.get(k) ?? 0) + 1); }
    for (const [mode, n] of freq) if (n >= 2 && (!best || n > best.n || (n === best.n && mode > best.mode && best.delim !== '\t'))) best = { delim, mode, n };
  }
  if (!best) return null;
  const { delim, mode } = best;
  const same = rows.filter(r => width(r.l, delim) === mode);
  for (let k = 0; k < same.length - 1; k++) {
    const cells = splitDelimited(same[k].l, delim).map(c => c.trim()).filter(Boolean);
    // a header: mostly words, no call's numbers; the row after it a call
    if (cells.length && cells.filter(c => !isNum(c)).length * 2 > cells.length && dataScore(cells) < 2 && dataScore(splitDelimited(same[k + 1].l, delim)) >= 2) return { delim, headerRow: same[k].i, ncols: mode };
  }
  // no header: rows that are calls from the first one
  return same.filter(r => dataScore(splitDelimited(r.l, delim)) >= 2).length * 2 >= same.length ? { delim, headerRow: -1, ncols: mode } : null;
}

const SYN: Record<Field, string[]> = {
  chrom: ['chromosome', 'chrom', 'chr', 'chromosom', 'contig'],
  start: ['start', 'debut', 'begin', 'startpos', 'startposition', 'chromstart', 'from', 'pos', 'position'],
  stop: ['stop', 'end', 'fin', 'endpos', 'endposition', 'chromend', 'finish', 'to'],
  pos: ['location', 'coordinates', 'coordinate', 'region', 'locus', 'interval', 'genomicposition', 'genomiccoordinates', 'cnvposition', 'localisation'],
  type: ['type', 'aberrationtype', 'aberration', 'cnvtype', 'event', 'evenement', 'variationtype', 'vartype', 'alteration', 'anomalie', 'call', 'status'],
  copies: ['copynumber', 'copies', 'copynb', 'nbcopies', 'cnstate', 'cn', 'ploidy'],
  iscn: ['iscn', 'nomenclature'],
  lr: ['meanlogratio', 'meanl2r', 'meanlog2ratio', 'meanlog2', 'meanlog', 'logratio', 'log2ratio', 'log2', 'lr', 'l2r', 'ratio'],
  probes: ['probes', 'nprobes', 'numprobes', 'numberofprobes', 'nbprobes', 'nbsondes', 'sondes', 'markers', 'nbmarkers', 'probe'],
  size: ['sizekb', 'size', 'taille', 'length', 'len', 'sizebp', 'sizemb'],
  band: ['cytoband', 'band', 'bande', 'cytobande'],
  genes: ['genename', 'genenames', 'genes', 'gene', 'symbol', 'symbols', 'genesymbol'],
  cls: ['classification', 'class', 'classe'],
  note: ['comment', 'comments', 'commentaire', 'commentaires', 'note', 'notes', 'remarque', 'remarks'],
  sample: ['samplename', 'sample', 'recordname', 'arrayname', 'echantillon', 'patient', 'patientid', 'sampleid', 'id'],
  sup: ['suppress', 'suppressed'],
};
/**
 * The columns by their names, loosely (English and French synonyms, then a name beginning with one), then by their
 * cells (`rows`: a region in one cell, an ISCN, a chromosome followed by two positions, a type); one column, one field.
 */
/** the name of the ISCN column added by this page at export (read back, it is not the file's own ISCN) */
export const OWN_ISCN_HEADER = 'ISCN (Sashimi viewer)';
export const isOwnIscnHeader = (h: string) => /^\s*iscn\s*\(\s*sashimi viewer\s*\)\s*$/i.test(h);
export function guessCols(header: string[], rows: string[][] = []): ColMap {
  // the ISCN column this page added at an export is not the file's own: never guessed (it is rewritten at the next export)
  const h = header.map(norm), cols = noCols(), used = new Set<number>(header.flatMap((x, j) => (isOwnIscnHeader(x) ? [j] : [])));
  for (const pass of [0, 1]) for (const f of FIELDS) {
    if (cols[f] >= 0) continue;
    for (const s of SYN[f]) {
      const i = h.findIndex((x, j) => !used.has(j) && x !== '' && (pass === 0 ? x === s : s.length >= 3 && x.startsWith(s)));
      if (i >= 0) { cols[f] = i; used.add(i); break; }
    }
  }
  const sample = rows.slice(0, 30), n = Math.max(header.length, ...sample.map(r => r.length));
  /** the share of a column's filled cells that pass a test */
  const share = (j: number, ok: (v: string) => boolean) => { const v = sample.map(r => (r[j] ?? '').trim()).filter(Boolean); return v.length ? v.filter(ok).length / v.length : 0; };
  const take = (f: Field, j: number) => { for (const g of FIELDS) if (cols[g] === j) cols[g] = -1; cols[f] = j; used.add(j); };
  for (let j = 0; j < n && sample.length; j++) {
    // a column named for another field keeps it, unless named as a position (Position, Chr) and holding whole regions
    if (used.has(j) && cols.chrom !== j && cols.start !== j && cols.stop !== j) continue;
    if (cols.iscn < 0 && share(j, v => /arr\s*\[|\)\s*(x\d|amp)/i.test(v)) >= 0.8) take('iscn', j);
    else if (cols.pos < 0 && share(j, v => !!parseRegion(v)) >= 0.8) take('pos', j);
  }
  if (cols.chrom < 0 && cols.pos < 0 && sample.length) {
    const j = [...Array(n).keys()].find(k => !used.has(k) && share(k, v => /^(chr)?(\d{1,2}|X|Y)$/i.test(v)) >= 0.8);
    if (j != null && share(j + 1, isNum) >= 0.8 && share(j + 2, isNum) >= 0.8 && !used.has(j + 1) && !used.has(j + 2)) { take('chrom', j); if (cols.start < 0) take('start', j + 1); if (cols.stop < 0) take('stop', j + 2); }
  }
  if (cols.type < 0 && sample.length) {
    const j = [...Array(n).keys()].find(k => !used.has(k) && share(k, v => /del|dup|loss|gain|amp|perte|^x\d|^cn\b/i.test(v)) >= 0.8);
    if (j != null) take('type', j);
  }
  return cols;
}
/** log-ratio columns besides `cols.lr`: unused columns holding numbers, named like one (log, lr, ratio), or named gain, loss,
 * amplification, deletion with values within ±5; the first goes to `cols.lr` when it has none */
export function guessLrMore(header: string[], rows: string[][], cols: ColMap): { cols: ColMap; lrMore: number[] } {
  const used = new Set(FIELDS.map(f => cols[f]).filter(i => i >= 0)), more: number[] = [], sample = rows.slice(0, 30);
  const numeric = (j: number) => { const v = sample.map(r => (r[j] ?? '').trim()).filter(Boolean); return v.length > 0 && v.filter(isNum).length / v.length >= 0.8; };
  header.forEach((h, j) => {
    const n = norm(h);
    // a name saying log ratio is enough; one saying only gain / loss (a count, a size?) needs values a log ratio can take (|v| ≤ 5)
    const lrName = /log|^l2r|^lr|ratio/.test(n);
    if (used.has(j) || !(lrName || /^(amplification|gain|loss|deletion|dup|del)/.test(n)) || !numeric(j)) return;
    if (lrName || sample.every(r => { const v = parseNum(r[j]); return v == null || Math.abs(v) <= 5; })) more.push(j);
  });
  const c = { ...cols };
  if (c.lr < 0 && more.length) c.lr = more.shift()!;
  return { cols: c, lrMore: more };
}
/** the columns of an interval table as the software exports it: by exact name (the header already normalised) */
export function autoCols(head: string[]): ColMap {
  const col = (...names: string[]) => { for (const n of names) { const i = head.indexOf(n); if (i >= 0) return i; } return -1; };
  return {
    chrom: col('chromosome', 'chr', 'chrom'), start: col('start'), stop: col('stop', 'end'), band: col('cytoband', 'band'),
    genes: col('genename', 'genenames', 'genes', 'gene'), size: col('sizekb', 'size'), type: col('type', 'aberrationtype'), copies: -1, pos: -1,
    probes: col('probes', 'probe', 'nprobes', 'numprobes', 'numberofprobes'), lr: col('meanlogratio', 'logratio', 'log2ratio'),
    iscn: col('iscn'), sup: col('suppress', 'suppressed'), cls: col('classification', 'class'),
    note: col('comment', 'comments', 'note', 'notes', 'commentaire', 'commentaires'),
    sample: col('samplename', 'sample', 'recordname', 'arrayname'),
  };
}
export const AUTO_OPTS = (headerRow: number): MapOpts => ({ delim: '\t', headerRow, zeroBased: false, sizeUnit: 'kb', build: null, typeMap: {} });
export const SIZE_MULT = { bp: 1 / 1000, kb: 1, Mb: 1000 } as const;

/** the unit a Size column is in: the one whose values come closest to end − start */
export function guessSizeUnit(rows: string[][], c: ColMap, decimal: boolean): 'bp' | 'kb' | 'Mb' {
  const dev: Record<string, number[]> = { bp: [], kb: [], Mb: [] };
  for (const r of rows.slice(0, 60)) {
    const reg = c.pos >= 0 ? parseRegion(r[c.pos]) : null;
    const s = parseNum(r[c.size], decimal), a = reg?.start ?? parsePos(r[c.start]), b = reg?.end ?? parsePos(r[c.stop]);
    if (s == null || a == null || b == null || s <= 0 || b <= a) continue;
    for (const u of ['bp', 'kb', 'Mb'] as const) dev[u].push(Math.abs(Math.log(s * 1000 * SIZE_MULT[u] / (b - a))));
  }
  const med = (v: number[]) => (v.length ? v.sort((x, y) => x - y)[v.length >> 1] : Infinity);
  return (['kb', 'bp', 'Mb'] as const).reduce((m, u) => (med(dev[u]) < med(dev[m]) ? u : m), 'kb' as 'bp' | 'kb' | 'Mb');
}

// ---- remembered mappings: by the header's column names only ----
const MAPS_KEY = 'sashimi.cghColumnMaps';
const MAPS_MAX = 20;
interface Stored { sig: string; spec: MapSpec }
const loadMaps = (): Stored[] => { try { const v = JSON.parse(localStorage.getItem(MAPS_KEY) ?? '[]'); return Array.isArray(v) ? v : []; } catch { return []; } };
const sigOf = (cells: string[]) => cells.map(norm).join('|');
export function rememberMap(header: string[], spec: MapSpec) {
  // a data row taken for the header would store a patient's values, and match nothing anyway
  if (header.filter(h => h.trim()).length < 3 || !looksLikeHeader(header)) return;
  const sig = sigOf(header);
  try { localStorage.setItem(MAPS_KEY, JSON.stringify([{ sig, spec }, ...loadMaps().filter(m => m.sig !== sig)].slice(0, MAPS_MAX))); } catch { /* storage blocked */ }
}
/** the mapping remembered for a header found in the first lines (its line given), else null */
export function findMap(lines: string[]): { spec: MapSpec; row: number } | null {
  const maps = loadMaps();
  if (!maps.length) return null;
  for (const m of maps) {
    if (!m?.spec?.opts?.delim) continue;
    for (let i = 0; i < Math.min(lines.length, 40); i++) {
      if (lines[i].includes(m.spec.opts.delim) && sigOf(splitDelimited(lines[i], m.spec.opts.delim)) === m.sig) return { spec: m.spec, row: i };
    }
  }
  return null;
}
