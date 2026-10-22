/**
 * The human build a file was aligned on, from the lengths of the chromosomes in its header (@SQ LN): each primary
 * chromosome has a different length in GRCh38 and GRCh37 (and hg38 / hg19, their UCSC names), so one is enough. An
 * @SQ AS tag naming the build is taken first. Files on another reference (T2T-CHM13, NCBI36, a subset with renamed
 * sequences) are reported as such when known, else unknown: the page then says nothing.
 *
 * Lengths: Genome Reference Consortium assembly reports (GCA_000001405.15 for GRCh38, GCA_000001405.1 for GRCh37).
 */
import type { GenomeBuild } from './ensembl';

const GRCH38 = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
  114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415];
const GRCH37 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
  115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566];
/** chr1 of references the viewer has no annotation for */
const OTHER: Record<number, string> = { 248387328: 'T2T-CHM13', 247249719: 'NCBI36 (hg18)' };
const NAMES = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];

export interface AssemblyCall {
  /** the build, when it is one the page has annotation for */
  build: GenomeBuild | null;
  /** what the header says, in words ("chr1 is 248,956,422 bp: GRCh38") */
  note: string;
  /** another reference recognised (T2T-CHM13, NCBI36) */
  other?: string;
}

/** The build of a header's sequences, from its @SQ lines (or from `lengths` when the text has none: BAM binary refs). */
export function assemblyOfHeader(text: string, lengths?: { name: string; length: number }[]): AssemblyCall | null {
  const seqs: { name: string; length: number; as?: string }[] = [];
  for (const line of text.split('\n')) {
    if (!line.startsWith('@SQ')) continue;
    const f = Object.fromEntries(line.split('\t').slice(1).map(x => [x.slice(0, 2), x.slice(3)]));
    if (f.SN && f.LN) seqs.push({ name: f.SN, length: +f.LN, as: f.AS });
  }
  if (!seqs.length && lengths) seqs.push(...lengths);
  if (!seqs.length) return null;
  const as = seqs.find(s => s.as)?.as ?? '';
  if (/grch38|hg38/i.test(as)) return { build: 'GRCh38', note: `the header names the assembly ${as}` };
  if (/grch37|hg19|b37|hs37/i.test(as)) return { build: 'GRCh37', note: `the header names the assembly ${as}` };
  for (const s of seqs) {
    const i = NAMES.indexOf(s.name.replace(/^chr/i, '').toUpperCase());
    if (i < 0) continue;
    const len = s.length.toLocaleString('en-US');
    if (GRCH38[i] === s.length) return { build: 'GRCh38', note: `${s.name} is ${len} bp: GRCh38` };
    if (GRCH37[i] === s.length) return { build: 'GRCh37', note: `${s.name} is ${len} bp: GRCh37` };
    if (i === 0 && OTHER[s.length]) return { build: null, other: OTHER[s.length], note: `${s.name} is ${len} bp: ${OTHER[s.length]}` };
  }
  return null;
}
