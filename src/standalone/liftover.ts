/**
 * The window's coordinates in the other human build (GRCh38 ⇄ GRCh37), for the second coordinate line of the ruler.
 * Ensembl's assembly map (GET /map/human/{from}/{region}/{to}, https://rest.ensembl.org/documentation/info/assembly_map)
 * gives the pieces of the region that map, each as an `original` stretch and its `mapped` one; the lift of a position
 * is then a shift within its piece, done here for any position of the window without asking again. The mapping
 * comes from the GRCh37 contigs kept in GRCh38 (exact sequence): a stretch with no counterpart (a gap, a fixed or
 * replaced contig) has no piece and lifts to nothing.
 *
 * Only the chromosome and the window's coordinates are sent, in 1 Mb pieces cached for the page's life; the
 * reads, the samples and the file names stay on the computer.
 */
import type { GenomeBuild } from './ensembl';

/** 1-based inclusive stretches: `from` in the build of the view, `to` in the other one */
export interface LiftSegment {
  from: { start: number; end: number };
  to: { chrom: string; start: number; end: number; strand: 1 | -1 };
}

const HOSTS = ['https://rest.ensembl.org', 'https://grch37.rest.ensembl.org'];
const CHUNK = 1_000_000;
const cache = new Map<string, Promise<LiftSegment[]>>();

/** Ensembl's name of a chromosome (no "chr", MT for the mitochondrion) and back in the style of the view. */
const toEnsembl = (chrom: string) => { const c = chrom.replace(/^chr/i, ''); return /^m$/i.test(c) ? 'MT' : c; };
const fromEnsembl = (name: string, like: string) => {
  const prefixed = /^chr/i.test(like);
  if (name === 'MT') return prefixed ? 'chrM' : 'MT';
  return prefixed ? `chr${name}` : name;
};

async function fetchChunk(from: GenomeBuild, to: GenomeBuild, chrom: string, start: number, end: number, signal?: AbortSignal): Promise<LiftSegment[]> {
  const path = `/map/human/${from}/${toEnsembl(chrom)}:${start}..${end}:1/${to}?content-type=application/json`;
  let last: unknown = null;
  // the main server maps both ways; the GRCh37 one is the fallback when it does not answer
  for (const host of HOSTS) {
    try {
      const r = await fetch(host + path, { signal, headers: { Accept: 'application/json' } });
      if (!r.ok) throw new Error(`Ensembl ${r.status}`);
      return parseMap(await r.json(), from, chrom);
    } catch (e) {
      if ((e as { name?: string })?.name === 'AbortError') throw e;
      last = e;
    }
  }
  throw last instanceof Error ? last : new Error('the Ensembl assembly map did not answer');
}

/**
 * The answer's pieces. Which side is which is read from each side's `assembly`, not from the field names: answers
 * with `original` and `mapped` swapped have been reported (ensembl-dev, 2019).
 */
export function parseMap(d: unknown, from: GenomeBuild, chrom: string): LiftSegment[] {
  const list = (d as { mappings?: unknown[] })?.mappings;
  if (!Array.isArray(list)) return [];
  const out: LiftSegment[] = [];
  for (const m of list as { original?: any; mapped?: any }[]) {
    if (!m?.original || !m?.mapped) continue;
    const [src, dst] = m.mapped.assembly === from && m.original.assembly !== from ? [m.mapped, m.original] : [m.original, m.mapped];
    if (![src.start, src.end, dst.start, dst.end].every(Number.isFinite)) continue;
    out.push({ from: { start: src.start, end: src.end }, to: { chrom: fromEnsembl(String(dst.seq_region_name), chrom), start: dst.start, end: dst.end, strand: dst.strand === -1 ? -1 : 1 } });
  }
  return out.sort((a, b) => a.from.start - b.from.start);
}

/** The pieces over [start, end] (1-based inclusive) of `chrom`, from the 1 Mb chunks that cover it. */
export async function liftSegments(from: GenomeBuild, to: GenomeBuild, chrom: string, start: number, end: number, signal?: AbortSignal): Promise<LiftSegment[]> {
  const parts: Promise<LiftSegment[]>[] = [];
  for (let c = Math.floor((Math.max(1, start) - 1) / CHUNK); c * CHUNK < end; c++) {
    const key = `${from}>${to}:${chrom}:${c}`;
    let p = cache.get(key);
    if (!p) {
      p = fetchChunk(from, to, chrom, c * CHUNK + 1, (c + 1) * CHUNK, signal);
      cache.set(key, p);
      p.catch(() => cache.delete(key));
    }
    parts.push(p);
  }
  const all = (await Promise.all(parts)).flat();
  return all.filter(s => s.from.end >= start && s.from.start <= end);
}

/** A position (1-based) in the other build, or null where nothing maps. */
export function liftPos(segs: LiftSegment[], pos: number): { chrom: string; pos: number } | null {
  for (const s of segs) {
    if (pos < s.from.start || pos > s.from.end) continue;
    const k = pos - s.from.start;
    return { chrom: s.to.chrom, pos: s.to.strand === 1 ? s.to.start + k : s.to.end - k };
  }
  return null;
}

/**
 * A stretch (1-based inclusive) in the other build: from its ends when both map, else from the pieces inside it.
 * `partial` when some of it has no counterpart, `chroms` > 1 when it lands on several sequences.
 */
export function liftRange(segs: LiftSegment[], start: number, end: number): { chrom: string; start: number; end: number; partial: boolean; split: boolean } | null {
  const inside = segs.filter(s => s.from.end >= start && s.from.start <= end);
  if (!inside.length) return null;
  const pts: { chrom: string; pos: number }[] = [];
  let covered = 0;
  for (const s of inside) {
    const a = Math.max(start, s.from.start), b = Math.min(end, s.from.end);
    covered += b - a + 1;
    pts.push(liftPos([s], a)!, liftPos([s], b)!);
  }
  const chrom = pts[0].chrom, same = pts.filter(p => p.chrom === chrom);
  return { chrom, start: Math.min(...same.map(p => p.pos)), end: Math.max(...same.map(p => p.pos)), partial: covered < end - start + 1, split: same.length < pts.length };
}
