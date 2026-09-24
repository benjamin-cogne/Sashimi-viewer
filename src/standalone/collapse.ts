/**
 * Browser-side port of backend/read_collapse.py (see that module's docstring for the
 * method). All coordinates are 0-based half-open.
 */
import type { AlignedRead, ReadGroup, VariantSite } from '../components/sashimi/types';
import { depthArray } from './alignments';

const BASE_CODE: Record<string, number> = { A: 0, C: 1, G: 2, T: 3, N: 4 };
const CODE_BASE = 'ACGTN';
/** position × INDEL_KEY + length, for insertions and deletions (positions up to 1e9 stay exact) */
const INDEL_KEY = 8_388_608;

/**
 * Variant sites of a set of reads over [start, end): an allele needs `minAlt` reads and `minVaf` of the depth. Every
 * alternate base counts, whatever its base quality, as in the full scan (alleles.ts): leaving the low-quality ones out
 * of the count while their reads stay in the depth under-read every allele fraction, by a third on ONT data, where it
 * took homozygous sites (~0.67) for heterozygous ones and dropped heterozygous ones (~0.33) under the phasing window.
 * `minBq` is kept for the callers' signature; the phasing (alleleAt) still trusts only alleles of that quality.
 * Keys are numbers, not strings: this runs several times per collapse over every mismatch of every read.
 */
export function callSites(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, _minBq = 20, minIndel = 1): VariantSite[] {
  const depth = depthArray(reads, start, end), w = Math.max(0, end - start);
  // SNV counts in one typed array (position × 5 bases): no Map, no key, nothing for the collector
  const snv = new Int32Array(w * 5), otherSnv = new Map<string, number>(), ins = new Map<number, number>(), del = new Map<number, number>();
  const bump = <K,>(m: Map<K, number>, k: K) => m.set(k, (m.get(k) ?? 0) + 1);
  for (const r of reads) {
    for (const m of r.m) {
      const pos = m[0];
      if (pos < start || pos >= end) continue;
      const c = BASE_CODE[m[1]];
      if (c === undefined) bump(otherSnv, `${pos}\t${m[1]}`); else snv[(pos - start) * 5 + c]++;
    }
    for (const [pos, len] of r.i) if (len >= minIndel && pos >= start && pos < end) bump(ins, pos * INDEL_KEY + Math.min(len, INDEL_KEY - 1));
    for (const [ds, de] of r.d) if (de - ds >= minIndel && ds >= start && ds < end) bump(del, ds * INDEL_KEY + Math.min(de - ds, INDEL_KEY - 1));
  }
  const sites: VariantSite[] = [];
  const addSnv = (pos: number, base: string, n: number) => {
    const d = depth[pos - start];
    if (n >= minAlt && d && n / d >= minVaf) {
      const rb = ref && pos - refStart >= 0 && pos - refStart < ref.length ? ref[pos - refStart] : '?';
      sites.push({ pos, kind: 'snv', ref: rb, alt: base, length: 0, alt_count: n, depth: d, vaf: n / d });
    }
  };
  for (let i = 0; i < snv.length; i++) if (snv[i] >= minAlt) addSnv(start + Math.floor(i / 5), CODE_BASE[i % 5], snv[i]);
  for (const [k, n] of otherSnv) { const [p, base] = k.split('\t'); addSnv(parseInt(p), base, n); }
  for (const [k, n] of ins) {
    const pos = Math.floor(k / INDEL_KEY), len = k % INDEL_KEY;
    const d = pos - start < depth.length ? depth[pos - start] : 0;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos, kind: 'ins', ref: '', alt: `+${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  for (const [k, n] of del) {
    const ds = Math.floor(k / INDEL_KEY), len = k % INDEL_KEY;
    const d = depth[ds - start] + n;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos: ds, kind: 'del', ref: '', alt: `-${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  return sites.sort((x, y) => x.pos - y.pos || x.kind.localeCompare(y.kind) || (x.alt < y.alt ? -1 : x.alt > y.alt ? 1 : 0));
}

/**
 * A read's mismatches by position, for lookups in O(log n): alleleAt and readAlleles used to scan them all for every
 * site, which on a long read (hundreds of sequencing errors) times every site of the window was most of a collapse.
 * Built once per read, kept with it.
 */
const mismatchIndex = new WeakMap<AlignedRead, { pos: Int32Array; at: Int32Array }>();
/** The mismatch of `r` at `pos` (base, quality), or null. */
export function mismatchAt(r: AlignedRead, pos: number): [string, number] | null {
  let ix = mismatchIndex.get(r);
  if (!ix) {
    const order = r.m.map((_, i) => i).sort((i, j) => r.m[i][0] - r.m[j][0]);
    ix = { pos: Int32Array.from(order, i => r.m[i][0]), at: Int32Array.from(order) };
    mismatchIndex.set(r, ix);
  }
  let lo = 0, hi = ix.pos.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (ix.pos[mid] < pos) lo = mid + 1; else hi = mid; }
  if (lo >= ix.pos.length || ix.pos[lo] !== pos) return null;
  const m = r.m[ix.at[lo]];
  return [m[1], m[2]];
}
/** Whether `pos` lies in one of the read's aligned blocks (sorted, disjoint): bs <= pos < be, or bs < pos < be when `inner`. */
export function inBlocks(r: AlignedRead, pos: number, inner = false): boolean {
  const b = r.b;
  let lo = 0, hi = b.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (b[mid][1] <= pos) lo = mid + 1; else hi = mid; }
  return lo < b.length && (inner ? b[lo][0] < pos : b[lo][0] <= pos) && pos < b[lo][1];
}

function readAlleles(r: AlignedRead, sites: VariantSite[], minBq: number): string[] {
  return sites.map(s => {
    if (s.kind === 'snv') {
      if (!inBlocks(r, s.pos)) return '.';
      const m = mismatchAt(r, s.pos);
      if (m) return m[1] < minBq ? '.' : (m[0] === s.alt ? 'alt' : m[0]);
      return 'ref';
    }
    if (s.kind === 'ins') {
      const alt = r.i.some(([p, l]) => p === s.pos && l === s.length);
      const covered = alt || inBlocks(r, s.pos, true) || r.i.some(([p]) => p === s.pos);
      return !covered ? '.' : (alt ? 'alt' : 'ref');
    }
    const has = r.d.some(([a, b]) => a === s.pos && b === s.pos + s.length);
    return has ? 'alt' : inBlocks(r, s.pos) ? 'ref' : '.';
  });
}

function spliceChain(r: AlignedRead): [number, number][] {
  const dels = new Set(r.d.map(([s, e]) => `${s}\t${e}`));
  const chain: [number, number][] = [];
  for (let k = 0; k + 1 < r.b.length; k++) {
    const be = r.b[k][1], ns = r.b[k + 1][0];
    if (!dels.has(`${be}\t${ns}`) && ns > be) chain.push([be, ns]);
  }
  return chain;
}

function mergeBlocks(blocks: [number, number][]): [number, number][] {
  const out: [number, number][] = [];
  for (const [s, e] of [...blocks].sort((a, b) => a[0] - b[0] || a[1] - b[1])) {
    if (out.length && s <= out[out.length - 1][1]) out[out.length - 1][1] = Math.max(out[out.length - 1][1], e);
    else out.push([s, e]);
  }
  return out;
}

interface Group { chain: [number, number][]; alleles: string[]; n: number; reads: AlignedRead[]; blocks: [number, number][]; span?: [number, number]; chainSet?: Set<string>; spec?: number }
const specificity = (g: Group) => (g.spec ??= g.chain.length + g.alleles.filter(a => a !== '.').length);
/** blocks are merged (sorted, disjoint) before any comparison: the span is their first start and last end */
const span = (g: Group): [number, number] => (g.span ??= g.blocks.length ? [g.blocks[0][0], g.blocks[g.blocks.length - 1][1]] : [0, 0]);
const chainKey = (c: [number, number][]) => c.map(j => `${j[0]}-${j[1]}`).join(',');

function compatible(g: Group, anchor: Group): boolean {
  const aj = (anchor.chainSet ??= new Set(anchor.chain.map(j => `${j[0]}-${j[1]}`)));
  for (const j of g.chain) if (!aj.has(`${j[0]}-${j[1]}`)) return false;
  for (const [js, je] of anchor.chain) for (const [bs, be] of g.blocks) if (bs < je && js < be) return false;
  let shared = g.chain.length;
  for (let i = 0; i < g.alleles.length; i++) {
    const a = g.alleles[i], b = anchor.alleles[i];
    if (a !== '.') { if (b === '.' || a !== b) return false; shared++; }
  }
  if (shared) return true;
  const [gs, ge] = span(g), [as, ae] = span(anchor);
  return specificity(g) === 0 && gs < ae && as < ge;
}

/**
 * Long reads, clustered with a tolerance: a read (or a set of identical ones) joins the group whose consensus it
 * agrees with at their shared sites, up to TOLERANT_MAX_DISAGREE of them disagreeing (a sequencing error, a
 * mistyped site), and any read may seed a group. Exact patterns cannot group long reads: each covers its own run of
 * sites and carries its own errors, so no two are identical and nearly all fell into "minor" (301 of 324 on a
 * synthetic 30x ONT window, two haplotypes). A read agreeing with two groups goes to the one it agrees with by
 * TOLERANT_MARGIN more sites, else stays ambiguous between them.
 */
const TOLERANT_MAX_DISAGREE = 0.2, TOLERANT_MARGIN = 2;

/** chain part of compatible(): the group's junctions all in the anchor's, and none of its blocks across an anchor junction */
function chainCompatible(g: Group, anchor: Group): boolean {
  const aj = (anchor.chainSet ??= new Set(anchor.chain.map(j => `${j[0]}-${j[1]}`)));
  for (const j of g.chain) if (!aj.has(`${j[0]}-${j[1]}`)) return false;
  for (const [js, je] of anchor.chain) for (const [bs, be] of g.blocks) if (bs < je && js < be) return false;
  return true;
}

/** Tolerant clustering (long reads): fills anchors / members / ambiguous / minor as the exact placement does. */
function clusterTolerant(ordered: Group[], nSites: number, minSupport: number, anchors: Group[], members: Group[][],
  ambiguous: Map<string, { hits: number[]; groups: Group[] }>, minor: Group[]): void {
  interface Cluster { seed: Group; members: Group[]; tally: Map<string, number>[]; cons: string[]; n: number; span: [number, number] }
  const clusters: Cluster[] = [];
  const pending: { g: Group; hits: number[] }[] = [];
  const join = (c: Cluster, g: Group) => {
    g.alleles.forEach((a, i) => {
      if (a === '.') return;
      const t = (c.tally[i] ??= new Map());
      t.set(a, (t.get(a) ?? 0) + g.n);
      let best = '.', bn = 0, tie = false;
      for (const [k, v] of t) { if (v > bn) { best = k; bn = v; tie = false; } else if (v === bn) tie = true; }
      c.cons[i] = tie ? '.' : best;
    });
    c.n += g.n;
    const [gs, ge] = span(g); c.span = [Math.min(c.span[0], gs), Math.max(c.span[1], ge)];
  };
  // left to right, so that each read meets the group of the reads just before it over the sites they share (by how
  // informative they are, groups were seeded all along the window and never met)
  const byPos = [...ordered].sort((a, b) => span(a)[0] - span(b)[0] || specificity(b) - specificity(a));
  for (const g of byPos) {
    const calls: number[] = [];
    g.alleles.forEach((a, i) => { if (a !== '.') calls.push(i); });
    const scored: { i: number; agree: number }[] = [];
    clusters.forEach((c, ci) => {
      if (!chainCompatible(g, c.seed)) return;
      let agree = 0, disagree = 0;
      for (const i of calls) { const b = c.cons[i]; if (b === '.' || b === undefined) continue; if (b === g.alleles[i]) agree++; else disagree++; }
      if (!calls.length) { const [gs, ge] = span(g); if (gs < c.span[1] && c.span[0] < ge) scored.push({ i: ci, agree: 0 }); return; }
      if (agree >= 1 && disagree <= TOLERANT_MAX_DISAGREE * (agree + disagree)) scored.push({ i: ci, agree: agree - disagree });
    });
    scored.sort((a, b) => b.agree - a.agree);
    if (scored.length === 1 || (scored.length > 1 && scored[0].agree - scored[1].agree >= TOLERANT_MARGIN)) join(clusters[scored[0].i], g), clusters[scored[0].i].members.push(g);
    else if (scored.length > 1) pending.push({ g, hits: scored.map(x => x.i) });
    else if (calls.length) clusters.push({ seed: g, members: [], tally: [], cons: new Array(nSites).fill('.'), n: 0, span: span(g) }), join(clusters[clusters.length - 1], g);
    else minor.push(g);
  }
  // groups that agree where they overlap are one (a stretch without shared sites split them): merged while any pair
  // shares at least TOLERANT_MARGIN called sites with at most TOLERANT_MAX_DISAGREE of them disagreeing
  const alive = clusters.map(() => true);
  const mergeInto = new Map<number, number>();
  for (let changed = true; changed;) {
    changed = false;
    for (let x = 0; x < clusters.length; x++) {
      if (!alive[x]) continue;
      let best = -1, bestAgree = 0;
      for (let y = 0; y < clusters.length; y++) {
        if (y === x || !alive[y] || !chainCompatible(clusters[x].seed, clusters[y].seed)) continue;
        let agree = 0, disagree = 0;
        const cx = clusters[x].cons, cy = clusters[y].cons;
        for (let i = 0; i < nSites; i++) { if (cx[i] === '.' || cy[i] === '.') continue; if (cx[i] === cy[i]) agree++; else disagree++; }
        if (agree >= TOLERANT_MARGIN && disagree <= TOLERANT_MAX_DISAGREE * (agree + disagree) && agree - disagree > bestAgree) { best = y; bestAgree = agree - disagree; }
      }
      if (best < 0) continue;
      const [keep, drop] = clusters[x].n >= clusters[best].n ? [x, best] : [best, x];
      const k = clusters[keep], d = clusters[drop];
      for (const g of [d.seed, ...d.members]) { join(k, g); k.members.push(g); }
      alive[drop] = false; mergeInto.set(drop, keep); changed = true;
    }
  }
  const root = (c: number): number => { while (mergeInto.has(c)) c = mergeInto.get(c)!; return c; };
  for (const p of pending) p.hits = [...new Set(p.hits.map(root))];
  // groups of fewer than minSupport reads are minor patterns; ambiguity between surviving groups only
  const index = new Map<number, number>();
  clusters.forEach((c, ci) => {
    if (!alive[ci]) return;
    if (c.n < minSupport) { minor.push(c.seed, ...c.members); return; }
    index.set(ci, anchors.length);
    anchors.push({ ...c.seed, alleles: c.cons, spec: undefined, span: undefined });
    members.push(c.members);
  });
  for (const { g, hits } of pending) {
    const live = hits.map(h => index.get(h)).filter((x): x is number => x != null);
    if (live.length === 1) members[live[0]].push(g);
    else if (!live.length) minor.push(g);
    else { const k = live.join(','); const e = ambiguous.get(k); if (e) e.groups.push(g); else ambiguous.set(k, { hits: live, groups: [g] }); }
  }
}

export function collapseReads(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20, minSupport = 3, minIndel = 1, tolerant = false): { sites: VariantSite[]; groups: ReadGroup[]; total: number } {
  const sites = callSites(reads, start, end, ref, refStart, minAlt, minVaf, minBq, minIndel);
  const exact = new Map<string, Group>();
  for (const r of reads) {
    const chain = spliceChain(r), alleles = readAlleles(r, sites, minBq);
    const key = `${chainKey(chain)}|${alleles.join('')}`;
    let g = exact.get(key);
    if (!g) { g = { chain, alleles, n: 0, reads: [], blocks: [] }; exact.set(key, g); }
    g.n++; g.reads.push(r); g.blocks.push(...r.b.map(b => [b[0], b[1]] as [number, number]));
  }
  for (const g of exact.values()) g.blocks = mergeBlocks(g.blocks);
  const ordered = [...exact.values()].sort((a, b) => specificity(b) - specificity(a) || b.n - a.n);
  const anchors: Group[] = [], members: Group[][] = [], ambiguous = new Map<string, { hits: number[]; groups: Group[] }>(), minor: Group[] = [];
  const place = (g: Group, allowNew: boolean) => {
    const hits = anchors.map((a, i) => (compatible(g, a) ? i : -1)).filter(i => i >= 0);
    if (hits.length === 1) members[hits[0]].push(g);
    else if (hits.length > 1) { const k = hits.join(','); const e = ambiguous.get(k); if (e) e.groups.push(g); else ambiguous.set(k, { hits, groups: [g] }); }
    else if (allowNew) { anchors.push(g); members.push([]); }
    else minor.push(g);
  };
  if (!tolerant) {
    for (const g of ordered) if (g.n >= minSupport) place(g, true);
    for (const g of ordered) if (g.n < minSupport) place(g, false);
  } else clusterTolerant(ordered, sites.length, minSupport, anchors, members, ambiguous, minor);
  const total = reads.length;
  const envelope = (gs: Group[]) => {
    const all = gs.flatMap(g => g.reads);
    const blocks = mergeBlocks(all.flatMap(r => r.b.map(b => [b[0], b[1]] as [number, number])));
    const n = all.length;
    const lo = blocks[0][0], hi = blocks[blocks.length - 1][1];
    const d = depthArray(all, lo, hi), thr = Math.max(1, n * 0.25);
    const dense: [number, number][] = [];
    let cur: [number, number] | null = null;
    for (let i = 0; i < d.length; i++) {
      if (d[i] >= thr) { if (!cur) cur = [lo + i, lo + i + 1]; else cur[1] = lo + i + 1; }
      else if (cur) { dense.push(cur); cur = null; }
    }
    if (cur) dense.push(cur);
    return { blocks, dense, n };
  };
  let groups: ReadGroup[] = anchors.map((a, i) => {
    const { blocks, dense, n } = envelope([a, ...members[i]]);
    return { id: `H${i + 1}`, kind: 'consensus' as const, n, frac: total ? n / total : 0, chain: a.chain, alleles: a.alleles, blocks, dense, absorbed: n - a.n };
  });
  const order = groups.map((_, i) => i).sort((i, j) => groups[j].n - groups[i].n);
  const rename = new Map(order.map((i, rank) => [`H${i + 1}`, `H${rank + 1}`]));
  groups = order.map(i => ({ ...groups[i], id: rename.get(groups[i].id)! }));
  for (const { hits, groups: gs } of [...ambiguous.values()].sort((a, b) => b.groups.reduce((s, g) => s + g.n, 0) - a.groups.reduce((s, g) => s + g.n, 0))) {
    const { blocks, dense, n } = envelope(gs);
    const ids = hits.map(h => rename.get(`H${h + 1}`)!);
    groups.push({ id: ids.join('|'), kind: 'ambiguous', n, frac: total ? n / total : 0, compatible: ids, chain: [], alleles: sites.map(() => '.'), blocks, dense, absorbed: 0 });
  }
  if (minor.length) {
    const n = minor.reduce((s, g) => s + g.n, 0);
    groups.push({ id: 'minor', kind: 'minor', n, frac: total ? n / total : 0, patterns: minor.length, chain: [], alleles: sites.map(() => '.'), blocks: mergeBlocks(minor.flatMap(g => g.blocks)), dense: [], absorbed: 0 });
  }
  return { sites, groups, total };
}
