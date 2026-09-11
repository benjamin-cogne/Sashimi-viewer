/**
 * Variants previously identified in a sample, turned into genomic intervals the Sashimi viewer can
 * draw. Sources (all carried by the sample record of the API):
 *   - the clinical indication: `hgvs_genomic` (+ `extra_variants`) with the c. / p. notations
 *   - the diagnostic conclusion: `diag_hgvs_genomic` (+ `diag_extra_variants`), and for CNVs / SVs
 *     the free-text `diag_nomenclature` (ISCN `arr[GRCh38] 7q11.23(73,132,891-74,608,254)x1`,
 *     `DEL chr7:73.1-74.6 Mb`, or an HGVS `g.` description)
 *   - the chromosome map: `sv_regions` (typed intervals)
 * Notations are parsed here, on the client, so the viewer does not depend on a new API shape.
 * Coordinates out are 0-based half-open, chromosomes `chr`-prefixed.
 */
import type { KnownVariant, KnownVariantKind } from './types';

/** RefSeq chromosome accession versions of the two builds, to infer the build of an NC_ notation. */
const NC_VERSION: Record<'GRCh38' | 'GRCh37', Record<string, number>> = {
  GRCh38: { 1: 11, 2: 12, 3: 12, 4: 12, 5: 10, 6: 12, 7: 14, 8: 11, 9: 12, 10: 11, 11: 10, 12: 12, 13: 11, 14: 9, 15: 10, 16: 10, 17: 11, 18: 10, 19: 10, 20: 11, 21: 9, 22: 11, X: 11, Y: 10 },
  GRCh37: { 1: 10, 2: 11, 3: 11, 4: 11, 5: 9, 6: 11, 7: 13, 8: 10, 9: 11, 10: 10, 11: 9, 12: 11, 13: 10, 14: 8, 15: 9, 16: 9, 17: 10, 18: 9, 19: 9, 20: 10, 21: 8, 22: 10, X: 10, Y: 9 },
};

/** `NC_000017` → `17`, `NC_000023` → `X`, `NC_012920` → `M`. */
function chromOfAccession(acc: string): string | null {
  const n = Number(acc.replace(/^NC_0*/, ''));
  if (acc.startsWith('NC_012920')) return 'M';
  if (!Number.isFinite(n)) return null;
  if (n >= 1 && n <= 22) return String(n);
  if (n === 23) return 'X';
  if (n === 24) return 'Y';
  return null;
}

/** `chr17`, `17`, `chrM`, `MT` → `chr17` / `chrM`; empty when not a chromosome name. */
export function normalizeChrom(raw: string): string {
  let c = String(raw || '').trim().replace(/^chr/i, '').toUpperCase();
  if (c === 'MT') c = 'M';
  if (c === '23') c = 'X';
  if (c === '24') c = 'Y';
  return /^([1-9]|1\d|2[0-2]|X|Y|M)$/.test(c) ? `chr${c}` : '';
}

const buildOfAccession = (chrom: string, version: number): KnownVariant['build'] => {
  const c = chrom.replace(/^chr/, '');
  if (NC_VERSION.GRCh38[c] === version) return 'GRCh38';
  if (NC_VERSION.GRCh37[c] === version) return 'GRCh37';
  return null;
};

const LARGE = 50;   // bp: deletions / duplications at least this long are drawn as bands, below as marks

/** `41,500,000`, `41.5 Mb`, `41500 kb` → bp. */
function toBp(raw: string, unit?: string): number | null {
  const v = Number(raw.replace(/[,_'\s]/g, ''));
  if (!Number.isFinite(v) || v < 0) return null;
  const scale: Record<string, number> = { bp: 1, b: 1, kb: 1e3, k: 1e3, mb: 1e6, m: 1e6 };
  return Math.round(v * (unit ? scale[unit.toLowerCase()] ?? 1 : 1));
}

interface Parsed { kind: KnownVariantKind; chrom: string; start: number; end: number; build: KnownVariant['build']; label: string }

/**
 * An HGVS genomic description anywhere in the text: `NC_000017.11:g.43094464C>T`, `chr17:g.43094464C>T`,
 * `g.43094464_43094470del` (chromosome from `fallbackChrom`), `g.(?_123)_(456_?)del`, `dup`, `ins`, `delins`, `inv`.
 */
export function parseHgvsGenomic(text: string, fallbackChrom = ''): Parsed | null {
  const m = String(text || '').match(/(?:(NC_\d{6})\.(\d+)|(?:chr)?([0-9]{1,2}|X|Y|MT|M))?\s*:?\s*g\.([^\s;,]+)/i);
  if (!m) return null;
  let chrom = '', build: KnownVariant['build'] = null;
  if (m[1]) { const c = chromOfAccession(m[1]); chrom = c ? `chr${c}` : ''; if (chrom) build = buildOfAccession(chrom, Number(m[2])); }
  else if (m[3]) chrom = normalizeChrom(m[3]);
  else chrom = normalizeChrom(fallbackChrom);
  const desc = m[4];
  const posTok = desc.match(/^[\d_()?]+/)?.[0] ?? '';
  const ints = (posTok.match(/\d+/g) || []).map(Number).filter(n => n > 0);
  if (!ints.length) return null;
  const p1 = Math.min(...ints), p2 = Math.max(...ints);
  const op = desc.slice(posTok.length);
  let kind: KnownVariantKind = 'other';
  let start = p1 - 1, end = p2;
  const sub = op.match(/^([ACGTN]+)>([ACGTN]+)/i);
  if (sub) kind = sub[1].length === 1 && sub[2].length === 1 && p1 === p2 ? 'snv' : 'indel';
  else if (/^delins/i.test(op)) kind = 'indel';
  else if (/^del/i.test(op)) kind = end - start >= LARGE ? 'del' : 'indel';
  else if (/^dup/i.test(op)) kind = end - start >= LARGE ? 'dup' : 'indel';
  else if (/^ins/i.test(op)) { kind = 'ins'; end = Math.max(start + 1, end); }
  else if (/^inv/i.test(op)) kind = 'inv';
  else if (/^=/.test(op)) kind = 'other';
  else if (/^[ACGTN]+$/i.test(op) || op === '') kind = p1 === p2 ? 'snv' : 'other';
  const label = `g.${desc}`;
  return { kind, chrom, start, end: Math.max(end, start + 1), build, label };
}

/** ISCN array notation: `arr[GRCh38] 7q11.23(73,132,891-74,608,254)x1` (x0/x1 → deletion, x3+ → duplication). */
export function parseIscnArray(text: string): Parsed | null {
  const m = String(text || '').match(/(?:arr|seq)\s*\[\s*(GRCh38|GRCh37|hg38|hg19)\s*\]\s*([0-9]{1,2}|X|Y)\s*[pq][\d.]*(?:[pq][\d.]*)?\s*\(\s*([\d,.]+)\s*[-–_]\s*([\d,.]+)\s*\)\s*x\s*(\d+)/i);
  if (!m) return null;
  const chrom = normalizeChrom(m[2]);
  const a = toBp(m[3]), b = toBp(m[4]);
  if (!chrom || a == null || b == null) return null;
  const copies = Number(m[5]);
  const build = /38/.test(m[1]) ? 'GRCh38' : 'GRCh37';
  return { kind: copies <= 1 ? 'del' : copies >= 3 ? 'dup' : 'cnv', chrom, start: Math.min(a, b) - 1, end: Math.max(a, b), build, label: `${copies <= 1 ? 'DEL' : copies >= 3 ? 'DUP' : 'CNV'} ×${copies}` };
}

const TYPE_WORDS: Record<string, KnownVariantKind> = {
  DEL: 'del', DELETION: 'del', LOSS: 'del', DUP: 'dup', DUPLICATION: 'dup', GAIN: 'dup', AMP: 'dup', AMPLIFICATION: 'dup',
  INV: 'inv', INVERSION: 'inv', INS: 'ins', INSERTION: 'ins', CNV: 'cnv', BND: 'bnd', BREAKEND: 'bnd', TRA: 'bnd', TRANSLOCATION: 'bnd',
};

/** A span with a type word: `DEL chr7:73,132,891-74,608,254`, `chr7:73.1-74.6 Mb DUP`, `7 73132891 74608254 DEL`. */
export function parseSpan(text: string): Parsed | null {
  const t = String(text || '');
  const tw = t.match(/\b(DEL|DELETION|LOSS|DUP|DUPLICATION|GAIN|AMP|AMPLIFICATION|INV|INVERSION|INS|INSERTION|CNV|BND|BREAKEND|TRA|TRANSLOCATION)\b/i);
  const kind: KnownVariantKind = tw ? TYPE_WORDS[tw[1].toUpperCase()] : 'cnv';
  const num = "\\d[\\d,_']*(?:\\.\\d+)?", unit = '(?:\\s*(bp|kb|mb|b|k|m)\\b)?';
  const colon = t.match(new RegExp(`(?:chr)?([0-9]{1,2}|X|Y|MT|M)\\s*:\\s*(${num})${unit}\\s*(?:-|–|—|\\.\\.|_)\\s*(${num})${unit}`, 'i'));
  const bed = colon ? null : t.match(new RegExp(`(?:^|\\s)(?:chr)?([0-9]{1,2}|X|Y|MT|M)[ \\t]+(${num})[ \\t]+(${num})(?=\\s|$)`, 'i'));
  let chrom = '', a: number | null = null, b: number | null = null;
  if (colon) { chrom = normalizeChrom(colon[1]); const u = colon[3] || colon[5]; a = toBp(colon[2], colon[3] || u); b = toBp(colon[4], colon[5] || u); }
  else if (bed) { chrom = normalizeChrom(bed[1]); a = toBp(bed[2]); b = toBp(bed[3]); }
  if (!chrom || a == null || b == null) return null;
  const build: KnownVariant['build'] = /GRCh38|hg38/i.test(t) ? 'GRCh38' : /GRCh37|hg19/i.test(t) ? 'GRCh37' : null;
  return { kind, chrom, start: Math.max(0, Math.min(a, b) - 1), end: Math.max(a, b), build, label: tw ? tw[1].toUpperCase() : 'CNV' };
}

/** Whatever notation the text holds: HGVS g. first, then ISCN array, then a typed span. */
export function parseVariantText(text: string, fallbackChrom = ''): Parsed | null {
  return parseHgvsGenomic(text, fallbackChrom) ?? parseIscnArray(text) ?? parseSpan(text);
}

const SV_KIND: Record<string, KnownVariantKind> = { DEL: 'del', DUP: 'dup', INV: 'inv', INS: 'ins', CNV: 'cnv', BND: 'bnd' };

/** Short c. notation (`c.123A>G`) without its transcript prefix. */
const shortC = (s: string) => String(s || '').trim().replace(/^[^:]+:/, '');

/**
 * Every placeable variant of a sample record (the `/samples/{id}` payload), clinical indication first,
 * then diagnostic, then chromosome map; the same locus is listed once.
 */
export function knownVariantsOfSample(sample: any): KnownVariant[] {
  const out: KnownVariant[] = [];
  const seen = new Map<string, KnownVariant>();
  let seq = 0;
  const push = (p: Parsed | null, text: string, source: KnownVariant['source'], extra: { gene?: string; cdna?: string; protein?: string; label?: string } = {}) => {
    if (!p) return;
    const key = `${p.chrom}:${p.start}-${p.end}:${p.kind}`;
    const prev = seen.get(key);
    if (prev) { prev.cdna ||= extra.cdna; prev.protein ||= extra.protein; prev.gene ||= extra.gene; return; }
    const label = extra.cdna ? shortC(extra.cdna) : extra.protein ? shortC(extra.protein) : extra.label ?? p.label;
    const v: KnownVariant = { id: `kv${++seq}`, kind: p.kind, chrom: p.chrom, start: p.start, end: p.end, label, text, source, build: p.build,
      gene: extra.gene || undefined, cdna: extra.cdna || undefined, protein: extra.protein || undefined };
    seen.set(key, v);
    out.push(v);
  };
  const snv = (gene: string, g: string, c: string, pr: string, source: KnownVariant['source']) => {
    if (!String(g || '').trim()) return;
    push(parseVariantText(g), g, source, { gene, cdna: c, protein: pr });
  };
  if (sample) {
    snv(sample.gene_name, sample.hgvs_genomic, sample.hgvs_cdna, sample.hgvs_protein, 'indication');
    for (const v of Array.isArray(sample.extra_variants) ? sample.extra_variants : []) snv(v?.gene_name, v?.hgvs_genomic, v?.hgvs_cdna, v?.hgvs_protein, 'indication');
    const d = sample.diagnostic;
    if (d) {
      snv(d.diag_gene_name, d.diag_hgvs_genomic, d.diag_hgvs_cdna, d.diag_hgvs_protein, 'diagnostic');
      for (const v of Array.isArray(d.diag_extra_variants) ? d.diag_extra_variants : []) snv(v?.gene_name, v?.hgvs_genomic, v?.hgvs_cdna, v?.hgvs_protein, 'diagnostic');
      for (const line of String(d.diag_nomenclature || '').split(/[\n;]+/)) {
        const t = line.trim();
        if (t) push(parseVariantText(t), t, 'diagnostic', { gene: d.diag_gene_name || undefined });
      }
    }
    for (const r of Array.isArray(sample.sv_regions) ? sample.sv_regions : []) {
      const chrom = normalizeChrom(String(r?.chrom ?? ''));
      const start = Number(r?.start), end = Number(r?.end);
      if (!chrom || !Number.isFinite(start) || !Number.isFinite(end)) continue;
      const type = String(r?.type ?? 'CNV').toUpperCase();
      const kind = SV_KIND[type] ?? 'cnv';
      push({ kind, chrom, start: Math.min(start, end), end: Math.max(start, end, Math.min(start, end) + 1), build: null, label: type }, `${type} ${chrom}:${start}-${end}${r?.label ? ` ${r.label}` : ''}`, 'chromosome_map', { label: r?.label ? `${type} · ${r.label}` : type });
    }
  }
  return out;
}

/** Colours of the markers: small variants amber-red, structural types as on the chromosome map. */
export const KNOWN_VARIANT_COLORS: Record<KnownVariantKind, string> = {
  snv: '#dc2626', indel: '#dc2626', ins: '#06b6d4', del: '#a855f7', dup: '#10b981', inv: '#f59e0b', cnv: '#ec4899', bnd: '#64748b', other: '#dc2626',
};

export const KNOWN_VARIANT_KIND_NAMES: Record<KnownVariantKind, string> = {
  snv: 'SNV', indel: 'small indel', ins: 'insertion', del: 'deletion', dup: 'duplication', inv: 'inversion', cnv: 'CNV', bnd: 'breakend', other: 'variant',
};

/** True for variants drawn as a point mark (one base or a few), false for bands. */
export const isPointVariant = (v: KnownVariant) => v.end - v.start < LARGE;

const SOURCE_NAMES: Record<KnownVariant['source'], string> = { indication: 'clinical indication', diagnostic: 'diagnostic conclusion', chromosome_map: 'chromosome map' };

/** Tooltip text of a known variant. */
export function knownVariantTitle(v: KnownVariant): string {
  const span = v.end - v.start;
  const lines = [
    `${v.gene ? `${v.gene} · ` : ''}${KNOWN_VARIANT_KIND_NAMES[v.kind]} · ${SOURCE_NAMES[v.source]}`,
    `${v.chrom}:${(v.start + 1).toLocaleString()}${span > 1 ? `-${v.end.toLocaleString()} (${span >= 1e6 ? `${(span / 1e6).toFixed(2)} Mb` : span >= 1e3 ? `${Math.round(span / 1e3)} kb` : `${span} bp`})` : ''}`,
    v.text,
  ];
  if (v.cdna) lines.push(v.cdna);
  if (v.protein) lines.push(v.protein);
  if (v.build === 'GRCh37') lines.push('⚠ GRCh37 notation: the viewer is GRCh38, the position may be off');
  return lines.join('\n');
}
