/**
 * Browser-side port of backend/read_collapse.py (see that module's docstring for the
 * method). All coordinates are 0-based half-open.
 */
import type { AlignedRead, ReadGroup, VariantSite } from '../components/sashimi/types';
import { depthArray } from './alignments';

export function callSites(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20): VariantSite[] {
  const depth = depthArray(reads, start, end);
  const snv = new Map<string, number>(), ins = new Map<string, number>(), del = new Map<string, number>();
  const bump = (m: Map<string, number>, k: string) => m.set(k, (m.get(k) || 0) + 1);
  for (const r of reads) {
    for (const [pos, base, qual] of r.m) if (qual >= minBq && pos >= start && pos < end) bump(snv, `${pos}\t${base}`);
    for (const [pos, len] of r.i) if (pos >= start && pos < end) bump(ins, `${pos}\t${len}`);
    for (const [ds, de] of r.d) if (ds >= start && ds < end) bump(del, `${ds}\t${de}`);
  }
  const sites: VariantSite[] = [];
  for (const [k, n] of snv) {
    const [p, base] = k.split('\t'); const pos = parseInt(p);
    const d = depth[pos - start];
    if (n >= minAlt && d && n / d >= minVaf) {
      const rb = ref && pos - refStart >= 0 && pos - refStart < ref.length ? ref[pos - refStart] : '?';
      sites.push({ pos, kind: 'snv', ref: rb, alt: base, length: 0, alt_count: n, depth: d, vaf: n / d });
    }
  }
  for (const [k, n] of ins) {
    const [p, l] = k.split('\t'); const pos = parseInt(p), len = parseInt(l);
    const d = pos - start < depth.length ? depth[pos - start] : 0;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos, kind: 'ins', ref: '', alt: `+${len}`, length: len, alt_count: n, depth: d, vaf: n / d });
  }
  for (const [k, n] of del) {
    const [a, b] = k.split('\t'); const ds = parseInt(a), de = parseInt(b);
    const d = depth[ds - start] + n;
    if (n >= minAlt && d && n / d >= minVaf) sites.push({ pos: ds, kind: 'del', ref: '', alt: `-${de - ds}`, length: de - ds, alt_count: n, depth: d, vaf: n / d });
  }
  return sites.sort((x, y) => x.pos - y.pos || x.kind.localeCompare(y.kind));
}

function readAlleles(r: AlignedRead, sites: VariantSite[], minBq: number): string[] {
  const mism = new Map<number, [string, number]>();
  for (const [pos, base, qual] of r.m) mism.set(pos, [base, qual]);
  const insSet = new Set(r.i.map(([p, l]) => `${p}\t${l}`));
  const delSet = new Set(r.d.map(([s, e]) => `${s}\t${e}`));
  return sites.map(s => {
    if (s.kind === 'snv') {
      const covered = r.b.some(([bs, be]) => bs <= s.pos && s.pos < be);
      if (!covered) return '.';
      const m = mism.get(s.pos);
      if (m) return m[1] < minBq ? '.' : (m[0] === s.alt ? 'alt' : m[0]);
      return 'ref';
    }
    if (s.kind === 'ins') {
      const covered = r.b.some(([bs, be]) => bs < s.pos && s.pos < be) || r.i.some(([p]) => p === s.pos);
      return !covered ? '.' : (insSet.has(`${s.pos}\t${s.length}`) ? 'alt' : 'ref');
    }
    const has = delSet.has(`${s.pos}\t${s.pos + s.length}`);
    const covered = has || r.b.some(([bs, be]) => bs <= s.pos && s.pos < be);
    return !covered ? '.' : (has ? 'alt' : 'ref');
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

interface Group { chain: [number, number][]; alleles: string[]; n: number; reads: AlignedRead[]; blocks: [number, number][] }
const specificity = (g: Group) => g.chain.length + g.alleles.filter(a => a !== '.').length;
const span = (g: Group): [number, number] => g.blocks.length ? [Math.min(...g.blocks.map(b => b[0])), Math.max(...g.blocks.map(b => b[1]))] : [0, 0];
const chainKey = (c: [number, number][]) => c.map(j => `${j[0]}-${j[1]}`).join(',');

function compatible(g: Group, anchor: Group): boolean {
  const aj = new Set(anchor.chain.map(j => `${j[0]}-${j[1]}`));
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

export function collapseReads(reads: AlignedRead[], start: number, end: number, ref: string | null, refStart: number,
  minAlt = 3, minVaf = 0.05, minBq = 20, minSupport = 3): { sites: VariantSite[]; groups: ReadGroup[]; total: number } {
  const sites = callSites(reads, start, end, ref, refStart, minAlt, minVaf, minBq);
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
  for (const g of ordered) if (g.n >= minSupport) place(g, true);
  for (const g of ordered) if (g.n < minSupport) place(g, false);
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
