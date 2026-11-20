/**
 * Gene-level loss-of-function constraint from the gnomAD GraphQL API (https://gnomad.broadinstitute.org/api):
 * `gene(gene_symbol, reference_genome) { gnomad_constraint { pli oe_lof_upper } }`, gnomAD v4 on GRCh38 and v2.1.1 on
 * GRCh37. The API answers any origin (the server registers `cors()` app-wide) and limits each client to about 30
 * requests and a query cost of 100 per minute, each gene costing 1 with at most 25 per request: genes are asked 20 at a
 * time, and kept for the page's life. When the API cannot be reached, the UCSC `pliByGene` track gives gnomAD v2.1.1.
 */
import type { GenomeBuild } from './ensembl';
import { getConstraintV2 } from './ucsc';

export const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';
const BATCH = 20;
export interface ConstraintValue { pli: number | null; loeuf: number | null; version: 'v4' | 'v2.1.1' }

const cache = new Map<string, ConstraintValue | null>();
const key = (build: GenomeBuild, s: string) => `${build}|${s.toUpperCase()}`;

async function graphql(build: GenomeBuild, symbols: string[]): Promise<Map<string, ConstraintValue | null>> {
  const ref = build === 'GRCh38' ? 'GRCh38' : 'GRCh37';
  const body = `query {\n${symbols.map((s, i) => `  g${i}: gene(gene_symbol: ${JSON.stringify(s)}, reference_genome: ${ref}) { symbol gnomad_constraint { pli oe_lof_upper } }`).join('\n')}\n}`;
  const r = await fetch(GNOMAD_API, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query: body }) });
  if (!r.ok && r.status !== 400) throw new Error(`gnomAD API ${r.status}`);
  const j = await r.json();
  if (!j?.data) throw new Error(`gnomAD API: ${(j?.errors ?? []).map((e: any) => e.message).join('; ') || 'no data'}`);
  const out = new Map<string, ConstraintValue | null>();
  symbols.forEach((s, i) => {
    const c = j.data[`g${i}`]?.gnomad_constraint;
    out.set(s, c ? { pli: c.pli ?? null, loeuf: c.oe_lof_upper ?? null, version: build === 'GRCh38' ? 'v4' : 'v2.1.1' } : null);
  });
  return out;
}

/**
 * Constraint of each gene; `window` lets the UCSC fallback read the track there. Genes gnomAD does not know are left
 * out of the map. `source` says where the values came from.
 */
export async function getConstraint(build: GenomeBuild, symbols: string[], window: { chrom: string; start: number; end: number }): Promise<{ values: Map<string, ConstraintValue>; source: string }> {
  const values = new Map<string, ConstraintValue>();
  const todo = [...new Set(symbols)].filter(s => !cache.has(key(build, s)));
  let source = `gnomAD ${build === 'GRCh38' ? 'v4' : 'v2.1.1'} (gnomAD API)`;
  try {
    for (let i = 0; i < todo.length; i += BATCH) {
      const got = await graphql(build, todo.slice(i, i + BATCH));
      for (const [s, v] of got) cache.set(key(build, s), v);
    }
  } catch {
    // the API is out of reach (a proxy, a rate limit): gnomAD v2.1.1 from UCSC
    const v2 = await getConstraintV2(build, window.chrom, window.start, window.end);
    for (const s of symbols) { const x = v2.get(s); if (x) values.set(s, { ...x, version: 'v2.1.1' }); }
    return { values, source: 'gnomAD v2.1.1 (UCSC pliByGene; the gnomAD API was not reachable)' };
  }
  for (const s of symbols) { const v = cache.get(key(build, s)); if (v) values.set(s, v); }
  return { values, source };
}
