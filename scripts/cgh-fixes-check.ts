// Unit checks of the review fixes (invented data only): CBS hybrid search, chrX copy numbers without a sex, the
// header-less export's ISCN column, PARs per build, log-ratio column guessing, a report without a Genome line.
// Run: npx esbuild scripts/cgh-fixes-check.ts --bundle --platform=node --outfile=/tmp/cgh-fixes-check.cjs && node /tmp/cgh-fixes-check.cjs
import { detect, segmentChrom, agrees, chrYPart, EXACT_MAX } from '../src/components/sashimi/cghSegment';
import { parseCghFile, exportCghFile } from '../src/components/sashimi/cghIntervals';
import { estCopies, mosaicOf, parOf } from '../src/components/sashimi/cghMosaic';
import { guessLrMore } from '../src/components/sashimi/cghMapping';
import type { CghProbeChrom } from '../src/components/sashimi/cghProbes';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = JSON.stringify(got) === JSON.stringify(want); if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', JSON.stringify(got), ok ? '' : `| want ${JSON.stringify(want)}`); };
let seed = 3; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
const chromOf = (n: number, ev: [number, number, number][]): CghProbeChrom => {
  const lr = Float32Array.from({ length: n }, (_, i) => gauss() * 0.15 + ev.reduce((s, [a, b, v]) => s + (i >= a && i < b ? v : 0), 0));
  return { start: Uint32Array.from({ length: n }, (_, i) => i * 1000 + 1), end: Uint32Array.from({ length: n }, (_, i) => i * 1000 + 60), lr, name: [] };
};
// 2. a long stretch (hybrid search) finds a 300-probe loss in 40k probes and a 30-probe gain; a short one is searched exactly
{
  const p = chromOf(40000, [[12000, 12300, -0.9], [30000, 30030, 0.8]]);
  const t0 = performance.now(), r = detect({ sample: '', build: 'GRCh37', n: 40000, byChrom: new Map([['chr1', p]]), dlrs: null });
  eq('hybrid: events', r.events.map(e => [e.kind, Math.round(e.i0 / 10), Math.round(e.i1 / 10)]), [['loss', 1200, 1230], ['gain', 3000, 3003]]);
  console.log('     40k-probe stretch in', Math.round(performance.now() - t0), 'ms; EXACT_MAX', EXACT_MAX);
}
// 3. a typeless table read from copy numbers, no sex: chrX x2 (a male's gain) and x1 (a female's loss) both kept, flagged
{
  const t = ['Chromosome\tStart\tEnd\tCN', 'chrX\t10000000\t12000000\t2', 'chrX\t20000000\t21000000\t1', 'chr3\t1000000\t2000000\t2', 'chr3\t3000000\t4000000\t3', 'chrX\t100000\t200000\t2'].join('\n');
  const cols = { chrom: 0, start: 1, stop: 2, pos: -1, type: -1, copies: 3, iscn: -1, lr: -1, probes: -1, size: -1, band: -1, genes: -1, cls: -1, note: -1, sample: -1, sup: -1 };
  const f = parseCghFile(t, 'x.tsv', 'P1', { cols, opts: { delim: '\t', headerRow: 0, zeroBased: false, sizeUnit: 'kb', build: null, typeMap: {} } });
  eq('X no sex: kept', f.samples[0].cnvs.map(c => [c.chrom, c.copies, c.kind, !!c.cnSex]), [['chrX', 2, 'gain', true], ['chrX', 1, 'loss', true], ['chr3', 3, 'gain', false]]);
  eq('X no sex: warned', f.warnings.some(w => /chrX row.*no sex/.test(w)), true);
  const m = parseCghFile('Sex: male\n' + t, 'x.tsv', 'P1', { cols, opts: { delim: '\t', headerRow: 1, zeroBased: false, sizeUnit: 'kb', build: null, typeMap: {} } });
  eq('X male in file: x1 normal, x2 gain', m.samples[0].cnvs.map(c => [c.chrom, c.copies, c.kind, !!c.cnSex]), [['chrX', 2, 'gain', false], ['chr3', 3, 'gain', false]]);
}
// 6. header-less export with Classification / Comment columns inside the rows: the ISCN after the widest row, no data overwritten
{
  const t = ['chr1\t1000000\t2000000\tGain\tVUS\tnote a\t3\t0.5\textra1', 'chr2\t1000000\t2000000\tLoss\t\t\t5\t-0.9\textra2'].join('\n');
  const cols = { chrom: 0, start: 1, stop: 2, pos: -1, type: 3, copies: -1, iscn: -1, lr: 7, probes: 6, size: -1, band: -1, genes: -1, cls: 4, note: 5, sample: -1, sup: -1 };
  const f = parseCghFile(t, 'x.tsv', 'P1', { cols, opts: { delim: '\t', headerRow: -1, zeroBased: false, sizeUnit: 'kb', build: 'GRCh37', typeMap: {} } });
  const out = exportCghFile(f, c => c.classification, c => c.note, c => `ISCN-${c.chrom}`).split('\n');
  eq('header-less export', out, ['chr1\t1000000\t2000000\tGain\tVUS\tnote a\t3\t0.5\textra1\tISCN-chr1', 'chr2\t1000000\t2000000\tLoss\t\t\t5\t-0.9\textra2\tISCN-chr2']);
}
// 7. PARs per build: a GRCh38 chrX call at 2.75 Mb is in PAR1 (diploid in a male), in GRCh37 it is not
{
  const c = { id: 'x', sample: 'S', sampleKey: 'S', chrom: 'chrX', start: 2700000, end: 2750000, cytoband: '', genes: [], sizeKb: 50, kind: 'loss' as const, typeText: '', probes: 5, logRatio: -1, iscn: null, copies: null, copiesHi: null, amp: false, suppressed: false, classification: '', classRest: '', note: '', line: 0 };
  eq('GRCh38 PAR1: x1 (diploid)', estCopies(c, 'male', 'GRCh38')?.n, 1);
  eq('GRCh37 non-PAR: x0 (hemizygous)', estCopies(c, 'male', 'GRCh37')?.n, 0);
  eq('mosaic GRCh38 PAR: no sex needed', mosaicOf(c, undefined, 'GRCh38')?.ok, true);
  eq('GRCh38 PARs', parOf('GRCh38'), { chrX: [[10001, 2781479], [155701383, 156030895]], chrY: [[10001, 2781479], [56887903, 57217415]] });
  const y = { chrom: 'chrY', kind: 'gain', start: 56900000, end: 57000000 }, x = { chrom: 'chrX', kind: 'gain', start: 155713480, end: 155813480 };
  eq('GRCh38 PAR2: chrY call meets its chrX twin', agrees(y, x, 'GRCh38'), true);
  const py: CghProbeChrom = { start: Uint32Array.from([20000, 2700000, 2800000, 56900000]), end: Uint32Array.from([20060, 2700060, 2800060, 56900060]), lr: Float32Array.from([1, 1, 1, 1]), name: [] };
  eq('chrY PAR probes dropped, GRCh38', Array.from(chrYPart(py, true, 'GRCh38').lr).map(v => Number.isNaN(v) ? 'NaN' : v), ['NaN', 'NaN', 1, 'NaN']);
  eq('chrY PAR probes dropped, GRCh37', Array.from(chrYPart(py, true, 'GRCh37').lr).map(v => Number.isNaN(v) ? 'NaN' : v), ['NaN', 1, 1, 1]);
}
// 9. guessLrMore: a Gain column of counts is not a log ratio; LogRatio_gain is; Loss with values within ±5 is
{
  const cols = { chrom: 0, start: 1, stop: 2, pos: -1, type: -1, copies: -1, iscn: -1, lr: -1, probes: -1, size: -1, band: -1, genes: -1, cls: -1, note: -1, sample: -1, sup: -1 };
  const h = ['chr', 'start', 'end', 'Gain count', 'Loss', 'LogRatio gain'], rows = [['1', '1', '2', '120', '-0.8', '0'], ['2', '1', '2', '300', '0', '12']];
  const r = guessLrMore(h, rows, cols);
  eq('guessLrMore', [r.cols.lr, r.lrMore], [4, [5]]);
}
// 10. an interval-based report without a Genome line: a header line naming the build, else GRCh37 assumed and said
{
  const body = ['AberrationNo\tChr\tCytoband\tStart\tStop\t#Probes\tAmplification\tGain\tLoss\tDeletion\tpval\tGene Names', 'INV1', '1\tchr1\tp36\t1000000\t2000000\t10\t0\t0\t-0.8\t0\t1e-9\tA'];
  const a = parseCghFile(['Aberration Algorithm: ADM-2', 'Reference: GRCh38 design', '', ...body].join('\n'), 'r.xls');
  eq('report: build from another header line', [a.build, a.warnings.some(w => /no genome build/.test(w))], ['GRCh38', false]);
  const b = parseCghFile(['Aberration Algorithm: ADM-2', '', ...body].join('\n'), 'r.xls');
  eq('report: none, GRCh37 assumed and said', [b.build, b.warnings.some(w => /GRCh37 assumed/.test(w))], [null, true]);
}
// 2 bis: below EXACT_MAX the search is the exact one (same breakpoints as a brute force on a short stretch)
{
  const p = chromOf(600, [[200, 260, -0.7]]);
  eq('exact segments', segmentChrom('chr2', p, 0.15, 8).map(s => [s.i0, s.i1]), [[0, 199], [200, 259], [260, 599]]);
}
console.log(bad ? `${bad} FAILED` : 'all ok');
process.exit(bad ? 1 : 0);
