/**
 * Gene-name suggestions while a gene is typed: the genes the page already knows first (its views, the recent searches),
 * then HGNC's (approved symbols, then aliases and previous symbols) from the NLM Clinical Table Search Service, made for
 * autocompletion (https://clinicaltables.nlm.nih.gov/apidoc/genes/v4/doc.html). Only the letters typed are sent. When the
 * service cannot be reached the page's own genes are all that is suggested, without a message.
 */

export interface GeneSuggestion {
  /** the approved symbol, what is searched when the suggestion is picked */
  symbol: string;
  /** the gene's full name */
  name?: string;
  /** cytogenetic band */
  location?: string;
  /** an alias or previous symbol the letters matched, when not the symbol itself */
  via?: { kind: 'alias' | 'previous'; text: string };
  /** the page knows it (a view, a recent search) */
  local?: boolean;
}

/** suggestions start from this many letters, and at most this many are listed */
export const SUGGEST_MIN = 3;
export const SUGGEST_MAX = 10;

const HGNC_URL = 'https://clinicaltables.nlm.nih.gov/api/genes/v4/search';
const RECENT_KEY = 'sashimi.recentGenes';

/** Text that reads as a gene symbol: not coordinates, a c./n./g. position, an exon number or an Ensembl id. */
export function looksLikeGene(text: string): boolean {
  const t = text.trim();
  if (t.length < SUGGEST_MIN || t.length > 30) return false;
  if (/^(chr|ens[a-z]*\d|nm_|nr_|xm_|xr_|nc_|[cngmp]\.|exons?\b)/i.test(t)) return false;
  return /^[A-Za-z][A-Za-z0-9.\-@]*$/.test(t);
}

/** Genes searched before in this browser, most recent first. */
export function recentGenes(): string[] {
  try { const v = JSON.parse(localStorage.getItem(RECENT_KEY) ?? '[]'); return Array.isArray(v) ? v.filter((x): x is string => typeof x === 'string') : []; } catch { return []; }
}
/** Remembers a gene opened, for the suggestions of the next searches. */
export function rememberGene(symbol: string): void {
  if (!looksLikeGene(symbol)) return;
  try { localStorage.setItem(RECENT_KEY, JSON.stringify([symbol, ...recentGenes().filter(g => g.toUpperCase() !== symbol.toUpperCase())].slice(0, 30))); } catch { /* private mode */ }
}

const cache = new Map<string, Promise<GeneSuggestion[]>>();
/** the service failed: not asked again for a minute */
let downUntil = 0;

/** HGNC genes matching the letters typed, ranked: the symbol itself, symbols starting with them, then aliases and previous symbols. */
function hgnc(q: string, signal?: AbortSignal): Promise<GeneSuggestion[]> {
  const key = q.toUpperCase();
  let p = cache.get(key);
  if (p) return p;
  if (Date.now() < downUntil) return Promise.resolve([]);
  const url = `${HGNC_URL}?terms=${encodeURIComponent(q)}&maxList=${4 * SUGGEST_MAX}&sf=symbol,alias_symbol,prev_symbol&df=symbol,name&ef=alias_symbol,prev_symbol,location`;
  p = fetch(url, { signal }).then(r => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.json(); }).then((d: unknown) => parseHgnc(d, key));
  p.catch(e => { cache.delete(key); if ((e as { name?: string })?.name !== 'AbortError') downUntil = Date.now() + 60_000; });
  cache.set(key, p);
  return p;
}

/** The service's answer: [total, codes, extra fields by name, display rows [symbol, name]]. */
function parseHgnc(d: unknown, key: string): GeneSuggestion[] {
  if (!Array.isArray(d) || !Array.isArray(d[3])) return [];
  const extra = (d[2] && typeof d[2] === 'object' ? d[2] : {}) as Record<string, unknown[]>;
  const list = (v: unknown): string[] => (Array.isArray(v) ? v : typeof v === 'string' ? v.split(/[|,]\s*/) : []).map(String).filter(Boolean);
  const out: { s: GeneSuggestion; rank: number }[] = [];
  (d[3] as unknown[]).forEach((row, i) => {
    if (!Array.isArray(row) || !row[0]) return;
    const symbol = String(row[0]), up = symbol.toUpperCase();
    const s: GeneSuggestion = { symbol, name: row[1] ? String(row[1]) : undefined, location: extra.location?.[i] ? String(extra.location[i]) : undefined };
    let rank = up === key ? 0 : up.startsWith(key) ? 1 : 3;
    if (rank === 3) {
      const alias = list(extra.alias_symbol?.[i]).find(a => a.toUpperCase().startsWith(key));
      const prev = list(extra.prev_symbol?.[i]).find(a => a.toUpperCase().startsWith(key));
      if (prev) { s.via = { kind: 'previous', text: prev }; rank = 2; } else if (alias) { s.via = { kind: 'alias', text: alias }; rank = 2; }
    }
    out.push({ s, rank });
  });
  return out.sort((a, b) => a.rank - b.rank || a.s.symbol.length - b.s.symbol.length || a.s.symbol.localeCompare(b.s.symbol)).map(x => x.s);
}

/** Suggestions for the letters typed: the page's genes starting with them, then HGNC's, SUGGEST_MAX at most. */
export async function suggestGenes(text: string, local: string[], signal?: AbortSignal): Promise<GeneSuggestion[]> {
  const q = text.trim();
  if (!looksLikeGene(q)) return [];
  const key = q.toUpperCase();
  const seen = new Set<string>();
  const out: GeneSuggestion[] = [];
  for (const g of [...local, ...recentGenes()]) {
    const up = g.toUpperCase();
    if (!up.startsWith(key) || seen.has(up)) continue;
    seen.add(up); out.push({ symbol: g, local: true });
  }
  let remote: GeneSuggestion[] = [];
  try { remote = await hgnc(q, signal); } catch { /* offline, blocked: the page's own genes only */ }
  for (const r of remote) {
    const up = r.symbol.toUpperCase();
    const mine = out.find(o => o.symbol.toUpperCase() === up);
    if (mine) { mine.name ??= r.name; mine.location ??= r.location; continue; }
    if (seen.has(up)) continue;
    seen.add(up); out.push(r);
  }
  return out.slice(0, SUGGEST_MAX);
}
