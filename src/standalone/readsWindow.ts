/**
 * The part of a reads answer every source shares once it holds the window's reads: the supporting-reads selection and
 * the cap for sources that hold whole windows, then the long-read noise filters, the variant sites, and in collapsed
 * mode the consensus groups or the phased haplotypes. LocalDataSource (BAM, CRAM), exported pages and any registered
 * file kind (fileKinds.ts) answer through it, so a window reads the same whatever file it came from.
 */
import type { AlignedRead, ReadsResponse } from '../components/sashimi/types';
import type { ReadsOptions } from '../components/sashimi/datasource';
import { callSites, collapseReads } from './collapse';
import { arcReadFromAligned, supportsArc } from './arcSupport';
import { phaseReads } from './phasing';
import { haplotagCounts, windowHaplotypes } from './haplotypes';

/**
 * Long reads (ONT, PacBio): median aligned length above 1 kb. The length counts aligned bases only, not the skipped
 * introns (N) nor the deletions (D): a 2×100 RNA-seq read spliced over a 3 kb intron spans 3.2 kb of the genome but is
 * a short read (counting its intron took most RNA-seq tracks of multi-exon genes for long reads), and so is a 2×150 read
 * carrying an 8 kb deletion (the supporting reads of a deletion, shown alone, were taken for long reads).
 */
export function isLongRead(reads: { b: [number, number][]; d: [number, number][] }[]): boolean {
  if (!reads.length) return false;
  const lens = reads.map(r => { let n = 0; for (const [a, b] of r.b) n += b - a; return n; }).sort((a, b) => a - b);
  return lens[lens.length >> 1] > 1000;
}

/**
 * The reads a window answer keeps, from all the reads a source holds over it: with `support`, the reads supporting
 * that arc (every k-th past `maxReads`) and then their mates in the window; otherwise every k-th read past the cap
 * (40,000 in collapsed mode, at least 100 else). `total` counts the reads before the cap (the supporting ones with
 * `support`), `mates` the mates added.
 */
export function supportAndCap(all: AlignedRead[], chrom: string, maxReads: number, collapsed: boolean, support: ReadsOptions['support'] | undefined):
  { reads: AlignedRead[]; total: number; mates: number } {
  let reads = all;
  let total = reads.length, mates = 0;
  if (support) {
    let hits = reads.filter(r => supportsArc(arcReadFromAligned(r, chrom), chrom, support));
    total = hits.length;
    const cap = Math.max(1, maxReads);
    if (hits.length > cap) { const step = hits.length / cap; hits = Array.from({ length: cap }, (_, i) => hits[Math.floor(i * step)]); }
    const kept = new Set(hits), names = new Set(hits.filter(r => r.mp != null).map(r => r.n));
    const extra = reads.filter(r => !kept.has(r) && names.has(r.n));
    mates = extra.length;
    reads = [...hits, ...extra].sort((a, b) => a.s - b.s);
  } else {
    const cap = collapsed ? 40000 : Math.max(100, maxReads);
    if (reads.length > cap) { const step = reads.length / cap; reads = Array.from({ length: cap }, (_, i) => reads[Math.floor(i * step)]); }
  }
  return { reads, total, mates };
}

/** What the source says of the window besides its reads. */
export type ReadsAnswerBase = Pick<ReadsResponse, 'sample_id' | 'sample_name' | 'total' | 'shown' | 'supporting' | 'reference' | 'reference_source'>;

/** The answer for `reads` (already filtered and capped): sites, and in collapsed mode groups or haplotypes instead of reads. */
export function answerReads(base: ReadsAnswerBase, reads: AlignedRead[], start: number, end: number, mode: 'reads' | 'collapsed',
  minSupport: number, minVaf: number, opts?: ReadsOptions): ReadsResponse {
  const ref = base.reference?.seq ?? null, refStart = base.reference?.start ?? 0;
  const longReads = isLongRead(reads);
  const minIndel = longReads ? Math.max(1, opts?.longReadMinIndel ?? 1) : 1;
  const vaf = longReads ? Math.max(minVaf, opts?.longReadMinVaf ?? 0.2) : minVaf;
  const full = { ...base, long_reads: longReads, haplotags: haplotagCounts(reads) };
  if (mode === 'collapsed') {
    if (opts?.haplotypes !== 'any') {
      // the window's sites, called once: the phasing, the haplotypes' checks and the answer share them
      const sites = callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel);
      const phaseOf = () => phaseReads(reads, start, end, ref, refStart, 3, vaf, 20, minIndel, sites);
      const { phase, haplotypes } = windowHaplotypes(reads, start, end, ref, refStart, vaf, minIndel, opts?.phaseSource ?? 'auto', phaseOf, sites, longReads);
      return { ...full, reads: [], sites, groups: [], phase, haplotypes };
    }
    const summary = collapseReads(reads, start, end, ref, refStart, 3, vaf, 20, Math.max(1, minSupport), minIndel, longReads);
    return { ...full, reads: [], sites: summary.sites, groups: summary.groups };
  }
  return { ...full, reads, sites: callSites(reads, start, end, ref, refStart, 3, vaf, 20, minIndel), groups: [] };
}
