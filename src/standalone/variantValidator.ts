/**
 * A transcript HGVS (`NM_001018069.2:c.1180T>G`) moved to the genome by VariantValidator's REST API
 * (https://rest.variantvalidator.org; Freeman et al., Hum Mutat 2018). The service allows any origin (CORS) and asks
 * for at most 3 requests a second.
 *
 * Its answer (VariantValidator `valoutput.format_as_dict`): a `flag` ('gene_variant', 'warning', 'intergenic', …), one
 * key per result (the transcript description, or `validation_warning_N` / `validation_error_N`) and `metadata`. A result
 * has `validation_warnings` and `primary_assembly_loci`, keyed by build in lower case (`grch38`, `grch37`, `hg38`,
 * `hg19`), each with `hgvs_genomic_description` (on the NC_ accession) and `vcf` {chr, pos, ref, alt}.
 */
import type { GenomeBuild } from './ensembl';

const VV_HOST = 'https://rest.variantvalidator.org';
/** past this, the next service is asked */
export const VV_TIMEOUT_MS = 10_000;

export async function vvRecode(build: GenomeBuild, hgvs: string, timeoutMs = VV_TIMEOUT_MS): Promise<{ hgvsg: string; vcf: string | null; gene: string | null; warnings: string[] }> {
  const tx = hgvs.split(':')[0].replace(/\(.*\)$/, '');
  // Ensembl transcripts have their own endpoint; RefSeq and LRG the main one, limited to the transcript typed
  const endpoint = /^ENST/i.test(tx) ? 'variantvalidator_ensembl' : 'variantvalidator';
  const url = `${VV_HOST}/VariantValidator/${endpoint}/${build}/${encodeURIComponent(hgvs)}/${encodeURIComponent(tx)}?content-type=application%2Fjson`;
  const ctrl = new AbortController(), timer = setTimeout(() => ctrl.abort(), timeoutMs);
  let resp: Response;
  try { resp = await fetch(url, { headers: { Accept: 'application/json' }, signal: ctrl.signal }); }
  catch (e) { throw new Error((e as Error).name === 'AbortError' ? `VariantValidator gave no answer in ${timeoutMs / 1000} s` : 'VariantValidator unreachable'); }
  finally { clearTimeout(timer); }
  const body: any = await resp.json().catch(() => null);
  if (resp.status === 429) throw new Error('VariantValidator: too many requests, try again in a moment');
  if (!resp.ok || !body || typeof body !== 'object') throw new Error(`VariantValidator ${resp.status}${body?.error ? `: ${String(body.error).slice(0, 200)}` : ''}`);
  const results = Object.entries(body).filter(([k, v]) => k !== 'flag' && k !== 'metadata' && v && typeof v === 'object') as [string, any][];
  const warningsOf = (r: any): string[] => (Array.isArray(r?.validation_warnings) ? r.validation_warnings.map(String).filter(Boolean) : []);
  const key = build.toLowerCase();
  for (const [, r] of results) {
    const locus = r?.primary_assembly_loci?.[key];
    const g = locus?.hgvs_genomic_description ? String(locus.hgvs_genomic_description) : '';
    if (!/^NC_\d+\.\d+:g\./.test(g)) continue;
    const v = locus.vcf;
    const vcf = v && v.pos ? `${v.chr}-${v.pos}-${v.ref}-${v.alt}` : null;
    return { hgvsg: g, vcf, gene: r.gene_symbol ? String(r.gene_symbol) : null, warnings: warningsOf(r) };
  }
  const why = results.flatMap(([, r]) => warningsOf(r));
  if (body.flag === 'warning' && why.length) throw new Error(`VariantValidator: ${why.join('; ').slice(0, 400)}`);
  throw new Error(`VariantValidator gave no ${build} position${why.length ? `: ${why.join('; ').slice(0, 300)}` : ''}`);
}
