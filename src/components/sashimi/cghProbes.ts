/**
 * The probes' log ratios of an array (CytoGenomics' CGH_Log_Ratio.txt: `AMADID <design> hg19`, then Probe Name, ChrName,
 * Start, Stop, FeatureNum, Description, Name of Gene, Normalized log ratio(<sample>); NaN for a probe left out). Kept
 * as the file's bytes and read on first use into typed arrays per chromosome (about 12 bytes a probe), for a probe
 * view. Nothing is drawn from them yet.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { decodeText } from './cghFormats';
import { cghBuild, cghChrom } from './cghIntervals';

export interface CghProbeChrom { start: Uint32Array; end: Uint32Array; lr: Float32Array }
export interface CghProbes {
  sample: string; build: GenomeBuild | null;
  /** probes read (NaN included) */
  n: number;
  /** chr1…chrY, each sorted by start */
  byChrom: Map<string, CghProbeChrom>;
}

const LR_HEAD = /normalized\s+log\s*ratio\s*\(([^)]*)\)/i;
/** the sample a probe log-ratio file is of (its `Normalized log ratio(<sample>)` header), or null when the bytes are not one */
export function probeSampleOf(bytes: Uint8Array): string | null {
  const head = decodeText(bytes.subarray(0, 8192));
  for (const l of head.split(/\r\n|\n|\r/).slice(0, 30)) if (/^probe\s*name\t/i.test(l)) return l.match(LR_HEAD)?.[1].trim() ?? null;
  return null;
}

/** the file read: probes per chromosome, sorted by start */
export function parseProbes(bytes: Uint8Array): CghProbes {
  const lines = decodeText(bytes).split(/\r\n|\n|\r/);
  const hi = lines.findIndex(l => /^probe\s*name\t/i.test(l));
  if (hi < 0) throw new Error('not a probe log-ratio file');
  const h = lines[hi].split('\t').map(x => x.trim().toLowerCase());
  const cChr = h.findIndex(x => x === 'chrname' || x === 'chr' || x === 'chromosome'), cS = h.indexOf('start'), cE = h.indexOf('stop'), cL = h.findIndex(x => LR_HEAD.test(x));
  const sample = lines[hi].match(LR_HEAD)?.[1].trim() ?? '';
  const build = cghBuild(lines.slice(0, hi).join(' '));
  const acc = new Map<string, { s: number[]; e: number[]; v: number[] }>();
  let n = 0;
  for (let i = hi + 1; i < lines.length; i++) {
    const f = lines[i].split('\t');
    const chrom = cghChrom(f[cChr] ?? ''), s = Number(f[cS]), e = Number(f[cE]);
    if (!chrom || chrom === 'chrM' || !Number.isFinite(s) || !Number.isFinite(e)) continue;
    let a = acc.get(chrom); if (!a) acc.set(chrom, a = { s: [], e: [], v: [] });
    const v = Number(f[cL]);
    a.s.push(s); a.e.push(e); a.v.push(Number.isFinite(v) ? v : NaN); n++;
  }
  const byChrom = new Map<string, CghProbeChrom>();
  for (const [chrom, a] of acc) {
    const order = a.s.map((_, i) => i).sort((x, y) => a.s[x] - a.s[y]);
    byChrom.set(chrom, { start: Uint32Array.from(order, i => a.s[i]), end: Uint32Array.from(order, i => a.e[i]), lr: Float32Array.from(order, i => a.v[i]) });
  }
  return { sample, build, n, byChrom };
}
