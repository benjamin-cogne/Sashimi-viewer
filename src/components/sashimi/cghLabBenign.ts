/**
 * The laboratory's list of benign CNVs: a BED (chrom, start 0-based, end, optional name that may say del / dup /
 * loss / gain; a name saying neither makes the entry apply to both types; further columns ignored). A laboratory
 * reference, not patient data: it is kept in localStorage and reloaded at the next session. A call matching an entry
 * is only marked and suggested as Benign, never classified, and the priority is unchanged.
 */
import { cghChrom, type CghCnv, type CghFile, type CghKind } from './cghIntervals';

/** an entry matches a call when it shares this much of the CALL (overlap / call length) */
export const LAB_COVER = 0.8;
export const LAB_KEY = 'sashimi.cghLabBenign';
/** above this many characters the list is kept in memory only */
export const LAB_MAX_STORED = 2_000_000;

export interface LabEntry { chrom: string; /** 1-based inclusive */ start: number; end: number; name: string; type: CghKind | 'both' }
export interface LabBenign { fileName: string; entries: LabEntry[]; /** lines that held an interval / lines skipped as unreadable (headers and comments not counted) */ read: number; skipped: number }

const typeOf = (t: string): CghKind | 'both' => (/del|loss|perte|\bx[01]\b/i.test(t) ? 'loss' : /amp|gain|dup|\bx[3-9]\b/i.test(t) ? 'gain' : 'both');

/** tolerant: tab or blank separated, commas in numbers, chr prefix optional, start and end swapped if need be */
export function parseLabBenign(text: string, fileName: string): LabBenign {
  const entries: LabEntry[] = [];
  let skipped = 0;
  for (const raw of text.replace(/^﻿/, '').split(/\r\n|\n|\r/)) {
    const l = raw.trim();
    if (!l || /^(track|browser|#)/i.test(l)) continue;
    const f = raw.includes('\t') ? raw.split('\t').map(x => x.trim()) : l.split(/\s+/);
    const chrom = cghChrom(f[0] ?? ''), n = (x: string | undefined) => (x && /^\d[\d,]*$/.test(x.trim()) ? Number(x.replace(/,/g, '')) : NaN);
    let a = n(f[1]), b = n(f[2]);
    if (!chrom || !Number.isFinite(a) || !Number.isFinite(b)) { skipped++; continue; }
    if (b < a) [a, b] = [b, a];
    if (b === a) { skipped++; continue; }
    const name = (f[3] ?? '').trim();
    entries.push({ chrom, start: a + 1, end: b, name, type: typeOf(name) });
  }
  return { fileName, entries, read: entries.length, skipped };
}

export interface LabHit { entry: LabEntry; pct: number }
/** the best entry of a compatible type sharing LAB_COVER of the call or more */
export function labMatch(c: Pick<CghCnv, 'chrom' | 'start' | 'end' | 'kind'>, entries: LabEntry[]): LabHit | null {
  const len = c.end - c.start + 1;
  let best: LabHit | null = null;
  for (const e of entries) {
    if (e.chrom !== c.chrom || (e.type !== 'both' && e.type !== c.kind)) continue;
    const ov = Math.min(c.end, e.end) - Math.max(c.start, e.start) + 1;
    if (ov <= 0) continue;
    const pct = ov / len;
    if (pct >= LAB_COVER && (!best || pct > best.pct)) best = { entry: e, pct };
  }
  return best;
}

// memoised on the files and the list: the table redraws often
let memo: { files: CghFile[]; lab: LabBenign; hits: Map<string, LabHit> } | null = null;
export function labHits(files: CghFile[], lab: LabBenign | null): Map<string, LabHit> {
  if (!lab || !lab.entries.length) return new Map();
  if (memo && memo.files === files && memo.lab === lab) return memo.hits;
  const hits = new Map<string, LabHit>();
  for (const fl of files) for (const sm of fl.samples) for (const c of sm.cnvs) { const h = labMatch(c, lab.entries); if (h) hits.set(c.id, h); }
  memo = { files, lab, hits };
  return hits;
}

export function loadStoredLab(): LabBenign | null {
  try {
    const v = JSON.parse(localStorage.getItem(LAB_KEY) ?? 'null') as LabBenign | null;
    return v && Array.isArray(v.entries) && typeof v.fileName === 'string' ? v : null;
  } catch { return null; }
}
/** false when the list was not stored (too large, or storage blocked) */
export function storeLab(lab: LabBenign | null): boolean {
  try {
    if (!lab) { localStorage.removeItem(LAB_KEY); return true; }
    const s = JSON.stringify(lab);
    if (s.length > LAB_MAX_STORED) { localStorage.removeItem(LAB_KEY); return false; }
    localStorage.setItem(LAB_KEY, s);
    return true;
  } catch { return false; }
}
