/**
 * Deep links: another tool (a variant interpretation site such as MobiDetails, a report generator,
 * a database export) can open the viewer on a variant or a locus. The parameters travel in the URL
 * fragment (`#variant=…`), which never reaches a server log, or in the query string.
 *
 *   #variant=NC_000017.11:g.43094464G>A&label=BRCA1%20c.5266dupC&pad=100
 *
 *   variant   HGVS g. with an NC_ accession (its version gives the build), `chr17:g.43094464G>A`,
 *             a pseudo-VCF `17-43094464-G-A` / `chr17:43094464:G:A`, or a bare `chr17:43094464`.
 *             Several separated by commas. Deletions / duplications draw as bands.
 *   locus     A position or interval (`chr17:43094464`, `chr17:43000000-43100000`) to show, when the
 *             window should not be derived from the variants.
 *   pad       Flank in bp around a point variant or position (default 100).
 *   label     Text drawn next to the marker (c. / p. notation); one per variant, separated by commas.
 *   gene      Symbol used when no RefSeq gene covers the locus (fallback lookup).
 *   build     GRCh38 (default) or GRCh37; only needed when no variant carries an accession.
 *   reads     `1` opens with the reads track on (off by default).
 *
 * The alignments still come from the user's files: a browser page cannot fetch a BAM by itself, so
 * the page waits for the first file and then opens the view on the requested window.
 */
import type { KnownVariant } from '../components/sashimi/types';
import { parseHgvsGenomic, normalizeChrom } from '../components/sashimi/knownVariants';
import { parseLocus } from '../components/sashimi/geometry';
import type { GenomeBuild } from './ensembl';

export interface LinkRequest {
  build: GenomeBuild;
  /** Variants to draw (0-based half-open, as the viewer's known-variant layer expects). */
  variants: KnownVariant[];
  /** Window to open, 1-based inclusive. */
  view: { chrom: string; start: number; end: number };
  /** Point or interval the locus label points at, 1-based inclusive. */
  mark: { chrom: string; start: number; end: number };
  pad: number;
  gene?: string;
  reads: boolean;
}

const DEFAULT_PAD = 100;
const MAX_PAD = 5_000_000;

function readParams(hash: string, search: string): URLSearchParams {
  const h = hash.replace(/^#\/?/, '');
  const p = new URLSearchParams(h.includes('=') ? h : '');
  if (!p.has('variant') && !p.has('locus')) return new URLSearchParams(search.replace(/^\?/, ''));
  return p;
}

/** Pseudo-VCF `17-43094464-G-A`, `chr17:43094464:G:A`, `17 43094464 G A` (1-based position). */
function parsePseudoVcf(text: string): { chrom: string; pos: number; ref: string; alt: string } | null {
  const m = text.trim().match(/^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[-:\s_]+(\d+)[-:\s_]+([ACGTN]+)[-:\s_]+([ACGTN]+)$/i);
  if (!m) return null;
  const chrom = normalizeChrom(m[1]);
  return chrom ? { chrom, pos: Number(m[2]), ref: m[3].toUpperCase(), alt: m[4].toUpperCase() } : null;
}

function variantOf(text: string, index: number, label: string | undefined, fallbackChrom: string): KnownVariant | null {
  const t = text.trim();
  if (!t) return null;
  const hgvs = parseHgvsGenomic(t, fallbackChrom);
  if (hgvs && hgvs.chrom) {
    return { id: `link${index}`, kind: hgvs.kind, chrom: hgvs.chrom, start: hgvs.start, end: hgvs.end, label: label || hgvs.label, text: t, source: 'indication', build: hgvs.build };
  }
  const vcf = parsePseudoVcf(t);
  if (vcf) {
    // VCF anchors indels on the base before the event; the affected span is the reference allele
    const isSnv = vcf.ref.length === 1 && vcf.alt.length === 1;
    const start = isSnv ? vcf.pos - 1 : vcf.pos;
    const end = Math.max(start + 1, vcf.pos - 1 + vcf.ref.length);
    const kind: KnownVariant['kind'] = isSnv ? 'snv' : vcf.ref.length > vcf.alt.length && end - start >= 50 ? 'del' : vcf.alt.length > vcf.ref.length && vcf.alt.length - vcf.ref.length >= 50 ? 'dup' : 'indel';
    return { id: `link${index}`, kind, chrom: vcf.chrom, start, end, label: label || `${vcf.chrom}:${vcf.pos} ${vcf.ref}>${vcf.alt}`, text: t, source: 'indication', build: null };
  }
  const locus = parseLocus(t);
  if (locus) {
    const point = locus.start === locus.end;
    return { id: `link${index}`, kind: point ? 'snv' : 'other', chrom: locus.chrom, start: locus.start - 1, end: locus.end, label: label || `${locus.chrom}:${locus.start.toLocaleString('en-US')}${point ? '' : `-${locus.end.toLocaleString('en-US')}`}`, text: t, source: 'indication', build: null };
  }
  return null;
}

/**
 * A variant of interest typed by the user: a locus (chr17:43,094,464 or chr17:43,094,464-43,094,470), an HGVS
 * genomic notation (NC_000017.11:g.43094464A>G, chr17:g.43094464A>G) or a VCF-like line (chr17 43094464 A G).
 * `label` is the text drawn next to it (a gene and protein change, say); null when nothing parses.
 */
export function variantOfInterest(text: string, label: string, id: string): KnownVariant | null {
  const v = variantOf(text, 0, label.trim() || undefined, '');
  return v ? { ...v, id } : null;
}

/** The request carried by the page URL, or null when there is none (or nothing parses). */
export function parseLink(hash: string, search: string): LinkRequest | null {
  const p = readParams(hash, search);
  const variantText = p.get('variant') || p.get('variants') || '';
  const locusText = p.get('locus') || p.get('region') || p.get('pos') || '';
  if (!variantText && !locusText) return null;

  const padRaw = Number(p.get('pad'));
  const pad = Number.isFinite(padRaw) && padRaw > 0 ? Math.min(MAX_PAD, Math.round(padRaw)) : DEFAULT_PAD;
  const labels = (p.get('label') || '').split(',').map(s => s.trim());
  const locus = locusText ? parseLocus(locusText) : null;
  const fallbackChrom = locus?.chrom ?? '';
  const variants = variantText.split(',').map((v, i) => variantOf(v, i + 1, labels[i] || undefined, fallbackChrom)).filter((v): v is KnownVariant => !!v);
  if (!variants.length && !locus) return null;

  const buildParam = (p.get('build') || p.get('genome') || '').toUpperCase();
  const build: GenomeBuild = /37|19/.test(buildParam) ? 'GRCh37' : /38/.test(buildParam) ? 'GRCh38' : (variants.find(v => v.build)?.build ?? 'GRCh38');
  // The build warning of the tooltip is for notations that disagree with the page build
  for (const v of variants) if (v.build === build) v.build = null;

  let mark: LinkRequest['mark'];
  if (locus) mark = { chrom: locus.chrom, start: locus.start, end: locus.end };
  else {
    const chrom = variants[0].chrom;
    const same = variants.filter(v => v.chrom === chrom);
    mark = { chrom, start: Math.min(...same.map(v => v.start)) + 1, end: Math.max(...same.map(v => v.end)) };
  }
  // The window is the mark plus the flank on each side, for a single base as for a band
  const view = { chrom: mark.chrom, start: Math.max(1, mark.start - pad), end: mark.end + pad };

  const gene = (p.get('gene') || '').trim() || undefined;
  const reads = /^(1|true|yes|on)$/i.test(p.get('reads') || '');
  return { build, variants, view, mark, pad, gene, reads };
}

/** One line describing the request, for the waiting notice. */
export function describeLink(req: LinkRequest): string {
  const what = req.variants.length ? req.variants.map(v => v.label === v.text ? v.text : `${v.text} (${v.label})`).join(', ') : `${req.mark.chrom}:${req.mark.start.toLocaleString('en-US')}${req.mark.end > req.mark.start ? `-${req.mark.end.toLocaleString('en-US')}` : ''}`;
  return `${what} · ${req.build} · window ${req.view.chrom}:${req.view.start.toLocaleString('en-US')}-${req.view.end.toLocaleString('en-US')}`;
}
