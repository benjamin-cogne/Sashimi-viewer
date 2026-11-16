/**
 * Column mapping of an array-CGH table: the delimited text (tab, `;` or `,`, quoted fields) is split, its header row
 * and delimiter guessed, its columns matched to the fields of a call by name (English and French variants), and what the
 * user chose is remembered by the header's column names (and the Type words given loss or gain), never a patient's row.
 */
import type { GenomeBuild } from '../../standalone/ensembl';

export const norm = (s: string) => s.toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '').replace(/[^a-z0-9]/g, '');

export type Field = 'chrom' | 'start' | 'stop' | 'band' | 'genes' | 'size' | 'type' | 'copies' | 'probes' | 'lr' | 'iscn' | 'sup' | 'cls' | 'note' | 'sample';
/** the column index of each field, -1 when the table has none */
export type ColMap = Record<Field, number>;
const FIELDS: Field[] = ['chrom', 'start', 'stop', 'type', 'copies', 'iscn', 'lr', 'probes', 'size', 'band', 'genes', 'cls', 'note', 'sample', 'sup'];
const noCols = (): ColMap => Object.fromEntries(FIELDS.map(f => [f, -1])) as ColMap;
export type TypeChoice = 'loss' | 'gain' | 'ignore';
export interface MapOpts {
  delim: string;
  /** line of the header (0-based, in the file's lines) */
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
export interface MapSpec { cols: ColMap; opts: MapOpts }
/** thrown by the automatic reading when a table looks delimited but its columns are not recognised */
export class NeedsMapping extends Error {}

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
export const joinDelimited = (cells: string[], d: string) => (d === '\t' ? cells.join(d) : cells.map(c => (c.includes(d) || c.includes('"') || /[\r\n]/.test(c) ? `"${c.replace(/"/g, '""')}"` : c)).join(d));

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

/** the delimiter and the header row of a delimited table, null when the text does not look like one (≥ 3 columns) */
export function detectTable(lines: string[]): { delim: string; headerRow: number; ncols: number } | null {
  const rows = lines.map((l, i) => ({ l, i })).filter(r => r.l.trim()).slice(0, 40);
  let best: { delim: string; mode: number; n: number } | null = null;
  for (const delim of ['\t', ';', ',']) {
    const freq = new Map<number, number>();
    for (const r of rows) { const k = splitDelimited(r.l, delim).length; if (k >= 3) freq.set(k, (freq.get(k) ?? 0) + 1); }
    for (const [mode, n] of freq) if (n >= 2 && (!best || n > best.n || (n === best.n && mode > best.mode && best.delim !== '\t'))) best = { delim, mode, n };
  }
  if (!best) return null;
  const same = rows.filter(r => splitDelimited(r.l, best!.delim).length === best!.mode);
  for (let k = 0; k < same.length - 1; k++) {
    const cells = splitDelimited(same[k].l, best.delim).map(c => c.trim()).filter(Boolean);
    const next = splitDelimited(same[k + 1].l, best.delim);
    if (cells.length && cells.filter(c => !isNum(c)).length * 2 > cells.length && next.filter(isNum).length >= 2) return { delim: best.delim, headerRow: same[k].i, ncols: best.mode };
  }
  return null;
}

const SYN: Record<Field, string[]> = {
  chrom: ['chromosome', 'chrom', 'chr', 'chromosom', 'contig'],
  start: ['start', 'debut', 'begin', 'startpos', 'startposition', 'chromstart', 'from', 'pos', 'position'],
  stop: ['stop', 'end', 'fin', 'endpos', 'endposition', 'chromend', 'finish', 'to'],
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
/** the columns by their names, loosely (English and French synonyms, then a name beginning with one); one column serves one field */
export function guessCols(header: string[]): ColMap {
  const h = header.map(norm), cols = noCols(), used = new Set<number>();
  for (const pass of [0, 1]) for (const f of FIELDS) {
    if (cols[f] >= 0) continue;
    for (const s of SYN[f]) {
      const i = h.findIndex((x, j) => !used.has(j) && x !== '' && (pass === 0 ? x === s : s.length >= 3 && x.startsWith(s)));
      if (i >= 0) { cols[f] = i; used.add(i); break; }
    }
  }
  return cols;
}
/** the columns of an interval table as the software exports it: by exact name (the header already normalised) */
export function autoCols(head: string[]): ColMap {
  const col = (...names: string[]) => { for (const n of names) { const i = head.indexOf(n); if (i >= 0) return i; } return -1; };
  return {
    chrom: col('chromosome', 'chr', 'chrom'), start: col('start'), stop: col('stop', 'end'), band: col('cytoband', 'band'),
    genes: col('genename', 'genenames', 'genes', 'gene'), size: col('sizekb', 'size'), type: col('type', 'aberrationtype'), copies: -1,
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
    const s = parseNum(r[c.size], decimal), a = parseNum(r[c.start]), b = parseNum(r[c.stop]);
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
  if (header.filter(h => h.trim()).length < 3) return;
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
