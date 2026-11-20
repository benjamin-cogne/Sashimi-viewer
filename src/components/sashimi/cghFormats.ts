/**
 * The bytes of an array-CGH file as text: UTF-8 (with or without BOM), UTF-16 (BOM, or NUL bytes every other byte),
 * else Windows-1252; a gzip unpacked; the first worksheet of an .xlsx as tab-separated text (so that it goes through the
 * same table reading). An old binary .xls is recognised to be refused with a clear message.
 */
import { unzipSync, gunzipSync } from 'fflate';

const starts = (b: Uint8Array, ...sig: number[]) => sig.every((x, i) => b[i] === x);
export const isGzip = (b: Uint8Array) => starts(b, 0x1f, 0x8b);
export const isZip = (b: Uint8Array) => starts(b, 0x50, 0x4b, 0x03, 0x04);
/** an OLE2 compound file: an Excel 97-2003 .xls (or another Office binary) */
export const isOle = (b: Uint8Array) => starts(b, 0xd0, 0xcf, 0x11, 0xe0);
export const OLE_MESSAGE = 'an Excel 97-2003 binary file (.xls) is not read: save it as .xlsx, .csv or tab-separated text';
export const gunzip = (b: Uint8Array) => gunzipSync(b);

/** UTF-16 without a BOM: NUL bytes on one side of most pairs of the start */
function utf16Guess(b: Uint8Array): 'utf-16le' | 'utf-16be' | null {
  const n = Math.min(b.length, 4000) & ~1;
  if (n < 8) return null;
  let even = 0, odd = 0;
  for (let i = 0; i < n; i += 2) { if (b[i] === 0) even++; if (b[i + 1] === 0) odd++; }
  const pairs = n / 2;
  return odd > pairs * 0.3 && even < pairs * 0.05 ? 'utf-16le' : even > pairs * 0.3 && odd < pairs * 0.05 ? 'utf-16be' : null;
}
/** a file's text: by its BOM, else UTF-16 by its NUL bytes, else UTF-8, else Windows-1252 (an export from a Windows PC) */
export function decodeText(b: Uint8Array): string {
  if (starts(b, 0xef, 0xbb, 0xbf)) return new TextDecoder('utf-8').decode(b.subarray(3));
  if (starts(b, 0xff, 0xfe)) return new TextDecoder('utf-16le').decode(b.subarray(2));
  if (starts(b, 0xfe, 0xff)) return new TextDecoder('utf-16be').decode(b.subarray(2));
  const u16 = utf16Guess(b);
  if (u16) return new TextDecoder(u16).decode(b);
  try { return new TextDecoder('utf-8', { fatal: true }).decode(b); } catch { return new TextDecoder('windows-1252').decode(b); }
}

// ---- .xlsx: the first worksheet ----
const ENT: Record<string, string> = { lt: '<', gt: '>', amp: '&', quot: '"', apos: "'" };
/** XML text: entities, and OOXML's _x000D_ escapes */
const xmlText = (s: string) => s
  .replace(/&(#x[0-9a-f]+|#\d+|lt|gt|amp|quot|apos);/gi, (_, e: string) => (e[0] === '#' ? String.fromCodePoint(e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : +e.slice(1)) : ENT[e.toLowerCase()]))
  .replace(/_x([0-9A-F]{4})_/g, (_, h: string) => String.fromCharCode(parseInt(h, 16)));
/** the text runs of a string item (<t>, rich-text <r><t>), phonetic runs (<rPh>) left out */
const runs = (xml: string) => [...xml.replace(/<rPh\b[\s\S]*?<\/rPh>/g, '').matchAll(/<t\b[^>]*>([\s\S]*?)<\/t>/g)].map(m => xmlText(m[1])).join('');
/** A1 → 0, AB7 → 27 */
const colIndex = (ref: string) => { let n = 0; for (const ch of ref.toUpperCase().match(/^[A-Z]+/)?.[0] ?? '') n = n * 26 + ch.charCodeAt(0) - 64; return n - 1; };
/** a number cell as the sheet shows it: binary float noise rounded away (0.30000000000000004 → 0.3) */
const numText = (v: string) => { const n = Number(v); return Number.isFinite(n) && !Number.isInteger(n) ? String(Number(n.toPrecision(15))) : v; };

/** the rows of an .xlsx's first worksheet (shared strings, inline strings, numbers, booleans), trailing empty cells dropped */
export function xlsxRows(bytes: Uint8Array): string[][] {
  const z = unzipSync(bytes, { filter: f => /^xl\/(workbook\.xml|_rels\/workbook\.xml\.rels|sharedStrings\.xml|worksheets\/[^/]+\.xml)$/i.test(f.name) });
  const get = (p: string) => { const k = Object.keys(z).find(x => x.toLowerCase() === p.toLowerCase()); return k ? new TextDecoder().decode(z[k]) : ''; };
  // the first sheet of the workbook, by its relationship; else the first worksheet file
  const rid = get('xl/workbook.xml').match(/<sheet\b[^>]*?\br:id="([^"]+)"/)?.[1];
  const rel = rid ? [...get('xl/_rels/workbook.xml.rels').matchAll(/<Relationship\b[^>]*>/g)].map(m => m[0]).find(r => r.includes(`Id="${rid}"`)) : undefined;
  const target = rel?.match(/Target="([^"]+)"/)?.[1];
  const sheetPath = target ? (target.startsWith('/') ? target.slice(1) : `xl/${target.replace(/^\.\//, '')}`) : Object.keys(z).filter(k => /^xl\/worksheets\//i.test(k)).sort((a, b) => a.localeCompare(b, 'en', { numeric: true }))[0];
  const sheet = sheetPath ? get(sheetPath) : '';
  if (!sheet) throw new Error('no worksheet in this .xlsx');
  const shared = [...get('xl/sharedStrings.xml').matchAll(/<si\b[^>]*>([\s\S]*?)<\/si>|<si\b[^>]*\/>/g)].map(m => (m[1] ? runs(m[1]) : ''));
  const rows: string[][] = [];
  let next = 0;
  for (const rm of sheet.matchAll(/<row\b([^>]*?)(?:\/>|>([\s\S]*?)<\/row>)/g)) {
    const r = Number(rm[1].match(/\br="(\d+)"/)?.[1] ?? next + 1) - 1;
    next = r + 1;
    const cells: string[] = [];
    let ci = 0;
    for (const cm of (rm[2] ?? '').matchAll(/<c\b([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)) {
      const ref = cm[1].match(/\br="([A-Z]+\d*)"/i)?.[1];
      const j = ref ? colIndex(ref) : ci;
      ci = j + 1;
      const t = cm[1].match(/\bt="(\w+)"/)?.[1] ?? 'n', body = cm[2] ?? '';
      const v = body.match(/<v>([\s\S]*?)<\/v>/)?.[1];
      const text = t === 's' ? shared[Number(v)] ?? '' : t === 'inlineStr' ? runs(body.match(/<is>([\s\S]*?)<\/is>/)?.[1] ?? '') : t === 'b' ? (v === '1' ? 'TRUE' : 'FALSE') : v == null ? '' : t === 'n' ? numText(v) : xmlText(v);
      while (cells.length < j) cells.push('');
      cells[j] = text;
    }
    while (cells.length && !cells[cells.length - 1].trim()) cells.pop();
    rows[r] = cells;
  }
  return Array.from(rows, r => r ?? []);
}
/** an .xlsx's first worksheet as tab-separated text (a tab or a line break inside a cell made a space) */
export const xlsxText = (bytes: Uint8Array) => xlsxRows(bytes).map(r => r.map(c => c.replace(/[\t\r\n]+/g, ' ')).join('\t')).join('\n');
