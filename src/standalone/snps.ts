/**
 * Common human variants for the Sashimi "Common SNPs" track, queried by the browser.
 *
 * Primary source: dbSNP build 155 "common" subset served by the UCSC Genome Browser REST API
 * (track dbSnp155Common on hg38 / hg19: variants with a minor allele frequency ≥ 1 % in at least
 * one frequency project; the track carries one MAF per project). Fallback when that API refuses
 * the request: Ensembl variation features of the window with 1000 Genomes MAF fetched in batches.
 * Results are cached per 100 kb chunk.
 */
import type { CommonSnp } from '../components/sashimi/types';
import type { GenomeBuild } from './ensembl';

const UCSC_API = 'https://api.genome.ucsc.edu';
const ENSEMBL: Record<GenomeBuild, string> = { GRCh38: 'https://rest.ensembl.org', GRCh37: 'https://grch37.rest.ensembl.org' };
const UCSC_GENOME: Record<GenomeBuild, string> = { GRCh38: 'hg38', GRCh37: 'hg19' };
export const SNP_CHUNK = 100_000;
/** Windows larger than this are not queried (the track asks to zoom in). */
export const SNP_MAX_WINDOW = 1_000_000;
export const SNP_MAX_INDEL = 20;

/**
 * Frequency projects of dbSNP 155 in the order used by the bigDbSnp `minorAlleleFreq` array
 * (dbSNP "freqSourceOrder"). Unknown positions are reported as "source N".
 */
export const DBSNP_FREQ_SOURCES = ['1000Genomes', 'dbGaP_PopFreq', 'TOPMED', 'KOREAN', 'SGDP_PRJ', 'Qatari', 'NorthernSweden', 'Siberian', 'TWINSUK', 'TOMMO',
  'ALSPAC', 'GENOME_DK', 'GnomAD', 'GoNL', 'Estonian', 'Vietnamese', 'Korea1K', 'HapMap', 'PRJEB36033', 'HGDP_Stanford', 'Daghestan', 'PAGE_STUDY',
  'Chileans', 'MGP', 'PRJEB37584', 'GnomAD_exomes', 'FINRISK', 'PharmGKB', 'PRJEB37766'];

const cache = new Map<string, Promise<CommonSnp[]>>();
let ucscBroken: string | null = null;   // set when the UCSC API refused the track, so later chunks go straight to Ensembl

const numList = (v: unknown): number[] => Array.isArray(v) ? v.map(Number) : typeof v === 'string' ? v.split(',').filter(Boolean).map(Number) : typeof v === 'number' ? [v] : [];
const strList = (v: unknown): string[] => Array.isArray(v) ? v.map(String) : typeof v === 'string' ? v.split(',').filter(Boolean) : [];

function classOf(ref: string, alts: string[], hint?: string): CommonSnp['cls'] {
  const h = (hint || '').toLowerCase();
  if (h === 'snv' || h === 'snp') return 'snv';
  if (h === 'ins' || h === 'insertion') return 'ins';
  if (h === 'del' || h === 'deletion') return 'del';
  if (h === 'delins' || h === 'indel') return 'delins';
  if (h === 'mnv' || h === 'substitution') return 'mnv';
  const lens = alts.map(a => a.replace('-', '').length), r = ref.replace('-', '').length;
  if (r === 1 && lens.every(l => l === 1)) return 'snv';
  if (r === 0 || lens.every(l => l > r)) return 'ins';
  if (lens.every(l => l < r)) return 'del';
  return r === lens[0] ? 'mnv' : 'delins';
}

/** One UCSC dbSnp155Common item → CommonSnp (chromStart is 0-based). */
function fromUcsc(it: any): CommonSnp | null {
  const start = Number(it.chromStart), end = Number(it.chromEnd);
  if (!Number.isFinite(start)) return null;
  const alts = strList(it.alts);
  const mafs = numList(it.minorAlleleFreq);
  const afs = mafs.map((af, i) => ({ source: DBSNP_FREQ_SOURCES[i] ?? `source ${i + 1}`, af })).filter(x => x.af > 0 && Number.isFinite(x.af));
  const maxAf = afs.length ? Math.max(...afs.map(x => x.af)) : 0;
  return { id: String(it.name || ''), start, end: Math.max(end, start + 1), ref: String(it.ref ?? ''), alts, cls: classOf(String(it.ref ?? ''), alts, it.class), maxAf, afs, source: 'dbSNP155', impact: it.maxFuncImpact != null ? String(it.maxFuncImpact) : undefined };
}

async function ucscChunk(build: GenomeBuild, chrom: string, idx: number): Promise<CommonSnp[]> {
  const start = idx * SNP_CHUNK, end = start + SNP_CHUNK;
  const c = chrom.startsWith('chr') ? chrom : `chr${chrom}`;
  const url = `${UCSC_API}/getData/track?genome=${UCSC_GENOME[build]};track=dbSnp155Common;chrom=${encodeURIComponent(c)};start=${start};end=${end};maxItemsOutput=100000`;
  const resp = await fetch(url, { headers: { Accept: 'application/json' } });
  if (!resp.ok) throw new Error(`UCSC API ${resp.status}`);
  const data = await resp.json();
  const items = data?.dbSnp155Common ?? (Array.isArray(data) ? data : null);
  if (!Array.isArray(items)) throw new Error(data?.error ? String(data.error) : 'UCSC API: unexpected answer');
  return items.map(fromUcsc).filter((x): x is CommonSnp => !!x);
}

/** Ensembl fallback: dbSNP variants of the chunk, then their 1000 Genomes MAF in batches of 200. */
async function ensemblChunk(build: GenomeBuild, chrom: string, idx: number): Promise<CommonSnp[]> {
  const host = ENSEMBL[build];
  const bare = chrom.startsWith('chr') ? chrom.slice(3) : chrom;
  const start = idx * SNP_CHUNK + 1, end = (idx + 1) * SNP_CHUNK;
  const r = await fetch(`${host}/overlap/region/human/${bare}:${start}-${end}?feature=variation;content-type=application/json`, { headers: { Accept: 'application/json' } });
  if (!r.ok) throw new Error(`Ensembl ${r.status}`);
  const feats: any[] = await r.json();
  const rs = feats.filter(f => String(f.id || '').startsWith('rs')).slice(0, 4000);
  const out: CommonSnp[] = [];
  for (let i = 0; i < rs.length; i += 200) {
    const ids = rs.slice(i, i + 200).map(f => f.id);
    const p = await fetch(`${host}/variation/homo_sapiens?content-type=application/json`, { method: 'POST', headers: { 'Content-Type': 'application/json', Accept: 'application/json' }, body: JSON.stringify({ ids }) });
    if (!p.ok) throw new Error(`Ensembl ${p.status}`);
    const info: Record<string, any> = await p.json();
    for (const f of rs.slice(i, i + 200)) {
      const v = info[f.id];
      const maf = v?.MAF != null ? Number(v.MAF) : NaN;
      if (!Number.isFinite(maf) || maf <= 0) continue;
      const alleles = String(f.alleles?.join?.('/') ?? f.alleles ?? '').split('/');
      const ref = alleles[0] ?? '', alts = alleles.slice(1);
      out.push({ id: f.id, start: Number(f.start) - 1, end: Number(f.end), ref, alts, cls: classOf(ref, alts), maxAf: maf, afs: [{ source: '1000Genomes', af: maf }], source: 'ensembl', impact: f.consequence_type });
    }
  }
  return out;
}

function chunk(build: GenomeBuild, chrom: string, idx: number): Promise<CommonSnp[]> {
  const key = `${build}|${chrom}|${idx}`;
  if (!cache.has(key)) {
    const run = async () => {
      if (!ucscBroken) {
        try { return await ucscChunk(build, chrom, idx); }
        catch (e: any) { ucscBroken = e?.message || String(e); console.warn(`[sashimi] UCSC dbSnp155Common unavailable (${ucscBroken}); using Ensembl variation`); }
      }
      return ensemblChunk(build, chrom, idx);
    };
    cache.set(key, run().catch(e => { cache.delete(key); throw e; }));
  }
  return cache.get(key)!;
}

/** Common variants overlapping the 0-based half-open window (SNVs and indels ≤ SNP_MAX_INDEL bp), sorted by position. */
export async function getCommonSnps(build: GenomeBuild, chrom: string, start: number, end: number): Promise<CommonSnp[]> {
  if (end <= start) return [];
  if (end - start > SNP_MAX_WINDOW) throw new Error(`window larger than ${(SNP_MAX_WINDOW / 1000).toFixed(0)} kb`);
  const first = Math.floor(Math.max(0, start) / SNP_CHUNK), last = Math.floor((end - 1) / SNP_CHUNK);
  const parts: CommonSnp[][] = [];
  for (let i = first; i <= last; i++) parts.push(await chunk(build, chrom, i));
  const seen = new Set<string>();
  return parts.flat()
    .filter(v => v.end > start && v.start < end && v.end - v.start <= SNP_MAX_INDEL + 1 && v.alts.every(a => a.length <= SNP_MAX_INDEL + 1))
    .filter(v => { const k = `${v.id}:${v.start}`; if (seen.has(k)) return false; seen.add(k); return true; })
    .sort((a, b) => a.start - b.start);
}

/** Which source is in use, for the track header. */
export function snpSourceLabel(): string { return ucscBroken ? 'Ensembl variation · 1000 Genomes MAF' : 'dbSNP 155 common (UCSC)'; }
