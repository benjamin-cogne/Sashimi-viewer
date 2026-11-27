// Every probe taken into account (invented data only): parsing (header, unsorted rows, NaN, duplicate positions, chrM /
// unknown chromosomes), DLRS, segmentation (every non-NaN probe in one segment, chrY rule, hybrid search), a call's
// probes (ends inclusive, a probe straddling an end counted once), the panel's window (edges, NaN runs), moving median ends.
// Run: npx esbuild scripts/cgh-probe-audit-check.ts --bundle --platform=node --outfile=/tmp/cgh-probe-audit-check.cjs && node /tmp/cgh-probe-audit-check.cjs
import { parseProbes, dlrsOf, meanIn, visible, runningMedian, type CghProbeChrom } from '../src/components/sashimi/cghProbes';
import { detect, segmentChrom, chrYPart } from '../src/components/sashimi/cghSegment';

let fail = 0;
const ok = (c: boolean, m: string) => { if (!c) { fail++; console.log('FAIL', m); } };
const head = 'TYPE\tx\nAMADID\t000001\thg19\nProbe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(S1)\n';
const rows = [
  ['p3', '1', '3000', '3059', '0.30'], ['p1', '1', '1000', '1059', '0.10'], ['p2', '1', '2000', '2059', 'NaN'],
  ['d1', '1', '4000', '4059', '0.5'], ['d2', '1', '4000', '4059', '0.7'], // two names at one position: two features, both kept
  ['p5', '1', '5000', '5059', '-0.2'], ['m1', 'M', '10', '69', '0.1'], ['u1', 'Un_gl000220', '10', '69', '0.1'],
  ['x1', 'X', '100', '159', '0.05'], ['y1', 'Y', '20000', '20059', '0.4'], ['y2', 'Y', '9000000', '9000059', '-0.3'],
];
const text = head + rows.map(r => [r[0], r[1], r[2], r[3], '1', 'd', 'g', r[4]].join('\t')).join('\n') + '\n\n';
const p = parseProbes(new TextEncoder().encode(text));
const c1 = p.byChrom.get('chr1')!;
ok(p.n === 9, `stored ${p.n} of 9 (chrM, unknown left out)`);
ok([...c1.start].join() === '1000,2000,3000,4000,4000,5000', `chr1 sorted, duplicates kept: ${[...c1.start]}`);
ok(c1.name.slice(3, 5).join() === 'd1,d2', 'duplicate position: file order kept');
ok(Number.isNaN(c1.lr[1]), 'NaN kept in place');
ok(!p.byChrom.has('chrM') && p.byChrom.size === 3, 'chromosomes');
// DLRS: on the probes in position order (rows given shuffled), NaN bridged, chromosomes not joined
ok(dlrsOf(p.byChrom) === null, 'DLRS null under 10 differences');
const vals = [0.1, -0.1, 0.25, NaN, 0.05, -0.3, 0.2, 0.0, 0.15, -0.05, 0.1, 0.3, -0.2, 0.05];
const order = vals.map((_, i) => i).sort(() => 0).reverse();
const t2 = head + order.map(i => [`q${i}`, i < 7 ? '2' : '3', String(1000 * (i + 1)), String(1000 * (i + 1) + 59), '1', 'd', 'g', String(vals[i])].join('\t')).join('\n');
const diffs: number[] = [];
for (const [a0, b0] of [[0, 7], [7, 14]]) { const v = vals.slice(a0, b0).filter(x => !Number.isNaN(x)); for (let i = 1; i < v.length; i++) diffs.push(v[i] - v[i - 1]); }
const mu = diffs.reduce((a, b) => a + b) / diffs.length, want = Math.sqrt(diffs.reduce((a, b) => a + (b - mu) ** 2, 0) / (diffs.length - 1)) / Math.SQRT2;
const got = dlrsOf(parseProbes(new TextEncoder().encode(t2)).byChrom)!;
ok(diffs.length === 11 && Math.abs(got - want) < 1e-6, `DLRS ${got} vs ${want} (${diffs.length} differences)`);

// segmentation: every non-NaN probe in exactly one segment, per chromosome; the hybrid path too
const rnd = (() => { let x = 7; return () => ((x = (x * 1103515245 + 12345) % 2147483648) / 2147483648); })();
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
const chromOf = (n: number, nanEvery: number, ev?: [number, number, number]): CghProbeChrom => {
  const lr = new Float32Array(n);
  for (let i = 0; i < n; i++) lr[i] = i % nanEvery === 0 || (i > n / 2 && i < n / 2 + 40) ? NaN : 0.15 * gauss() + (ev && i >= ev[0] && i < ev[1] ? ev[2] : 0);
  return { start: Uint32Array.from({ length: n }, (_, i) => 1000 + i * 1000), end: Uint32Array.from({ length: n }, (_, i) => 1059 + i * 1000), lr, name: [] };
};
for (const n of [3000, 15000]) {
  const q = chromOf(n, 17, [n / 3, n / 3 + 60, -0.8]), s = segmentChrom('chr1', q, 0.15, 8), nn = [...q.lr].filter(v => !Number.isNaN(v)).length;
  ok(s.reduce((a, x) => a + x.n, 0) === nn, `n=${n}: segments hold ${s.reduce((a, x) => a + x.n, 0)} of ${nn} non-NaN probes`);
  ok(s.every((x, k) => k === 0 || x.i0 > s[k - 1].i1), `n=${n}: segments disjoint and ordered`);
  ok(s.some(x => x.mean < -0.5 && x.n >= 50), `n=${n}: the 60-probe loss found`);
}
const male = detect({ ...p, byChrom: new Map([['chrY', p.byChrom.get('chrY')!]]) }, undefined, true);
ok((male.segs.get('chrY') ?? []).reduce((a, x) => a + x.n, 0) === 1, 'male chrY: the PAR1 probe left out, the other in');
ok((detect(p, undefined, false).segs.get('chrY') ?? []).length === 0, 'female chrY: none');
ok([...chrYPart(p.byChrom.get('chrY')!, true, 'GRCh37').lr].filter(Number.isNaN).length === 1, 'chrY PAR rule');

// a call's probes: ends inclusive, a straddling probe counted once (midpoint), adjacent calls share none
const q: CghProbeChrom = { start: Uint32Array.from([100, 200, 290, 400, 500]), end: Uint32Array.from([159, 259, 349, 459, 559]), lr: Float32Array.from([1, 2, 3, NaN, 5]), name: [] };
ok(meanIn(q, 100, 259).n === 2, 'probes exactly at the call\'s ends counted');
ok(meanIn(q, 100, 300).n === 2 && meanIn(q, 301, 600).n === 2, 'a probe straddling a boundary counted on one side only');
ok(meanIn(q, 120, 259).n === 2 && meanIn(q, 140, 259).n === 1, 'a probe starting before the call counted by its midpoint');
ok(meanIn(q, 100, 559).n === 4 && meanIn(q, 100, 559).mean === 11 / 4, 'NaN skipped, mean of the rest');
// the window: a probe starting before it but reaching in is drawn; k/2 non-NaN extras each side across NaN
const wv = visible(q, 150, 400, 3);
ok(wv.join() === '0,1,2,4', `window 150-400 with k 3: ${wv}`);
const big = chromOf(200, 5), v2 = visible(big, 50_000, 80_000, 7), inside = [...big.start].map((s, i) => [s, i]).filter(([s, i]) => s >= 50_000 && s <= 80_000 && !Number.isNaN(big.lr[i])).map(([, i]) => i);
ok(inside.every(i => v2.includes(i)), 'every non-NaN probe of the window listed');
ok(v2.filter(i => big.start[i] < 50_000).length === 3 && v2.filter(i => big.start[i] > 80_000).length === 4, 'k/2 extras each side (NaN skipped)');
ok(runningMedian([5, 1, 9, 2, 7], 3).join() === '3,5,2,7,4.5', `running median ends: ${runningMedian([5, 1, 9, 2, 7], 3)}`);

console.log(fail ? `${fail} failed` : 'all passed');
process.exit(fail ? 1 : 0);
