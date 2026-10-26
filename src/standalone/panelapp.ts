/**
 * Gene panels from PanelApp (Genomics England; PanelApp Australia runs the same code). REST API v1, GET only, open to
 * any origin for /api/v1/: `/panels/` lists the active panels 100 per page (`name` filters them; the list is small
 * enough to be fetched whole and searched here), `/panels/<id>/` returns a panel with its genes inline, each with
 * `entity_name` (the symbol), `confidence_level` ('3' green, '2' amber, '1' red, '0' no list), `mode_of_inheritance`
 * and `phenotypes`. Without network access, a panel downloaded from PanelApp as TSV, or a plain list of symbols, is
 * read from a file instead.
 */
export interface PanelSource { id: string; label: string; base: string }
export const PANEL_SOURCES: PanelSource[] = [
  { id: 'ge', label: 'PanelApp (Genomics England)', base: 'https://panelapp.genomicsengland.co.uk/api/v1' },
  { id: 'aus', label: 'PanelApp Australia', base: 'https://panelapp-aus.org/api/v1' },
];
export interface PanelSummary { source: string; id: string; name: string; version: string; genes: number; disorders: string[] }
export interface PanelGeneEntry { symbol: string; confidence: 1 | 2 | 3; moi?: string; phenotypes?: string[] }
export interface LoadedPanel { source: string; id: string; name: string; version: string; genes: PanelGeneEntry[] }

const listCache = new Map<string, Promise<PanelSummary[]>>();
/** Every active panel of a source (all pages), cached for the page's life. */
export function listPanels(src: PanelSource): Promise<PanelSummary[]> {
  if (!listCache.has(src.id)) {
    listCache.set(src.id, (async () => {
      const out: PanelSummary[] = [];
      let url: string | null = `${src.base}/panels/?page=1`;
      for (let guard = 0; url && guard < 30; guard++) {
        const r = await fetch(url, { headers: { Accept: 'application/json' } });
        if (!r.ok) throw new Error(`${src.label}: HTTP ${r.status}`);
        const j: any = await r.json();
        for (const p of j?.results ?? []) out.push({ source: src.id, id: String(p.id), name: String(p.name ?? ''), version: String(p.version ?? ''), genes: Number(p.stats?.number_of_genes ?? 0), disorders: Array.isArray(p.relevant_disorders) ? p.relevant_disorders.map(String) : [] });
        url = typeof j?.next === 'string' ? j.next : null;
      }
      return out.sort((a, b) => a.name.localeCompare(b.name));
    })().catch(e => { listCache.delete(src.id); throw e; }));
  }
  return listCache.get(src.id)!;
}

const conf = (v: unknown): 1 | 2 | 3 | null => { const n = Number(v); return n >= 3 ? 3 : n === 2 ? 2 : n === 1 ? 1 : null; };

/** A panel's genes (confidence 0, "no list", left out). */
export async function loadPanel(src: PanelSource, id: string): Promise<LoadedPanel> {
  const r = await fetch(`${src.base}/panels/${encodeURIComponent(id)}/`, { headers: { Accept: 'application/json' } });
  if (!r.ok) throw new Error(`${src.label}: HTTP ${r.status}`);
  const p: any = await r.json();
  const genes: PanelGeneEntry[] = [];
  for (const g of p?.genes ?? []) {
    const symbol = String(g?.gene_data?.gene_symbol ?? g?.entity_name ?? '');
    const c = conf(g?.confidence_level);
    if (symbol && c) genes.push({ symbol, confidence: c, moi: g.mode_of_inheritance || undefined, phenotypes: Array.isArray(g.phenotypes) ? g.phenotypes.map(String) : undefined });
  }
  return { source: src.label, id: String(p?.id ?? id), name: String(p?.name ?? `panel ${id}`), version: String(p?.version ?? ''), genes };
}

/**
 * A panel from a file: a PanelApp TSV download (columns "Gene Symbol" / "Entity Name" and "Model_Of_Inheritance",
 * "Phenotypes"; the rating from a "GEL_Status" / confidence column or from "Sources" mentioning Expert Review Green /
 * Amber / Red), or a plain list of symbols (one per line or separated by commas), taken as green.
 */
export function parsePanelFile(name: string, text: string): LoadedPanel {
  const lines = text.split(/\r?\n/).filter(l => l.trim());
  const head = lines[0]?.split('\t') ?? [];
  const col = (re: RegExp) => head.findIndex(h => re.test(h.trim()));
  const iSym = col(/^(gene symbol|entity name|gene_symbol|entity_name|symbol)$/i);
  const genes: PanelGeneEntry[] = [];
  if (head.length > 1 && iSym >= 0) {
    const iConf = col(/gel_status|confidence|level/i), iMoi = col(/model_of_inheritance|mode_of_inheritance|moi/i), iPh = col(/^phenotypes?$/i), iSrc = col(/^sources?/i);
    for (const l of lines.slice(1)) {
      const f = l.split('\t');
      const symbol = (f[iSym] ?? '').trim();
      if (!symbol) continue;
      let c = iConf >= 0 ? conf(f[iConf]) : null;
      if (c == null && iSrc >= 0) { const s = f[iSrc] ?? ''; c = /expert review green/i.test(s) ? 3 : /expert review amber/i.test(s) ? 2 : /expert review red/i.test(s) ? 1 : null; }
      genes.push({ symbol, confidence: c ?? 3, moi: iMoi >= 0 ? f[iMoi]?.trim() || undefined : undefined, phenotypes: iPh >= 0 && f[iPh] ? f[iPh].split(/;\s*/).filter(Boolean) : undefined });
    }
  } else {
    for (const t of text.split(/[\s,;]+/)) if (/^[A-Za-z0-9][\w.-]*$/.test(t)) genes.push({ symbol: t, confidence: 3 });
  }
  return { source: 'file', id: name, name: name.replace(/\.(tsv|txt|csv)$/i, ''), version: '', genes };
}
