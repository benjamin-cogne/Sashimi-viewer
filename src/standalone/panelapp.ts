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

// ======================== A gene's modes of inheritance across panels ========================

/**
 * PanelApp's mode_of_inheritance values, by what they require:
 * AD     MONOALLELIC, autosomal or pseudoautosomal (NOT imprinted, or imprinted status unknown)
 * AD-pat MONOALLELIC … maternally imprinted (paternal allele expressed)
 * AD-mat MONOALLELIC … paternally imprinted (maternal allele expressed)
 * AR     BIALLELIC, autosomal or pseudoautosomal
 * AD+AR  BOTH monoallelic and biallelic (including "biallelic mutations cause a more SEVERE disease form")
 * XLR    X-LINKED: hemizygous mutation in males, biallelic mutations in females
 * XLD    X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease
 * MT     MITOCHONDRIAL
 * ?      Unknown, Other, or not given
 */
export type MoiClass = 'AD' | 'AD-pat' | 'AD-mat' | 'AR' | 'AD+AR' | 'XLR' | 'XLD' | 'MT' | '?';
export function classifyMoi(text: string | undefined | null): MoiClass {
  const t = String(text ?? '').toLowerCase();
  if (!t.trim()) return '?';
  if (t.includes('x-linked') || t.includes('x linked')) return /monoallelic mutations in females|females may/.test(t) ? 'XLD' : 'XLR';
  if (t.includes('mitochondrial')) return 'MT';
  if (t.startsWith('both')) return 'AD+AR';
  if (t.startsWith('biallelic')) return 'AR';
  if (t.startsWith('monoallelic')) return t.includes('maternally imprinted') ? 'AD-pat' : t.includes('paternally imprinted') ? 'AD-mat' : 'AD';
  return '?';
}
export interface GeneMode { panel: string; panelId: string; confidence: 1 | 2 | 3; moi: string; cls: MoiClass }
export interface GeneModes { symbol: string; source: string; entries: GeneMode[]; main: MoiClass; parts: { cls: MoiClass; panels: number; green: number }[] }

/** The modes of a gene's entries, the green ones first: the main mode is the commonest among the best-rated entries. */
export function summarizeModes(symbol: string, source: string, entries: GeneMode[]): GeneModes {
  const by = new Map<MoiClass, { cls: MoiClass; panels: number; green: number }>();
  for (const e of entries) { const p = by.get(e.cls) ?? { cls: e.cls, panels: 0, green: 0 }; p.panels++; if (e.confidence === 3) p.green++; by.set(e.cls, p); }
  const parts = [...by.values()].sort((a, b) => Number(a.cls === '?') - Number(b.cls === '?') || b.green - a.green || b.panels - a.panels);
  return { symbol, source, entries, main: parts[0]?.cls ?? '?', parts };
}

const modeCache = new Map<string, Promise<GeneModes>>();
/**
 * Every panel entry of a gene (rated green, amber or red) with its mode of inheritance: `/genes/<symbol>/`, every page;
 * PanelApp notes that this list can come in varying order with duplicates, so entries are kept once per panel; when it
 * fails, `/entities/?entity_name=<symbol>`. Cached for the page's life.
 */
export function geneModes(src: PanelSource, symbol: string): Promise<GeneModes> {
  const key = `${src.id}|${symbol.toUpperCase()}`;
  if (!modeCache.has(key)) {
    modeCache.set(key, (async () => {
      const read = async (first: string): Promise<GeneMode[]> => {
        const out = new Map<string, GeneMode>();
        let url: string | null = first;
        for (let guard = 0; url && guard < 20; guard++) {
          const r = await fetch(url, { headers: { Accept: 'application/json' } });
          if (!r.ok) throw new Error(`${src.label}: HTTP ${r.status}`);
          const j: any = await r.json();
          for (const g of Array.isArray(j) ? j : j?.results ?? []) {
            const sym = String(g?.entity_name ?? g?.gene_data?.gene_symbol ?? '');
            if (sym.toUpperCase() !== symbol.toUpperCase() || (g?.entity_type && g.entity_type !== 'gene')) continue;
            const c = conf(g?.confidence_level), pid = String(g?.panel?.id ?? g?.panel ?? '');
            if (!c || !pid || out.has(pid)) continue;
            out.set(pid, { panel: String(g?.panel?.name ?? `panel ${pid}`), panelId: pid, confidence: c, moi: String(g?.mode_of_inheritance ?? ''), cls: classifyMoi(g?.mode_of_inheritance) });
          }
          url = typeof j?.next === 'string' ? j.next : null;
        }
        return [...out.values()];
      };
      let entries: GeneMode[];
      try { entries = await read(`${src.base}/genes/${encodeURIComponent(symbol)}/`); }
      catch { entries = await read(`${src.base}/entities/?entity_name=${encodeURIComponent(symbol)}`); }
      return summarizeModes(symbol, src.label, entries);
    })().catch(e => { modeCache.delete(key); throw e; }));
  }
  return modeCache.get(key)!;
}
