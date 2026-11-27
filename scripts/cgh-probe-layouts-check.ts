// Unit checks of the other probe-level layouts (cghProbes.ts `probeLayoutOf`, `parseProbes` with a spec, `guessProbeCols`)
// on the invented files of `scripts/cgh-probe-fixture.mjs --layouts DIR` (the same probes and events as the CytoGenomics fixture).
// Run: node scripts/cgh-probe-fixture.mjs /tmp/f.zip --split /tmp/po --layouts /tmp/lay && npx esbuild scripts/cgh-probe-layouts-check.ts --bundle --platform=node --outfile=/tmp/cgh-probe-layouts-check.cjs && node /tmp/cgh-probe-layouts-check.cjs /tmp/lay /tmp/po
import { readFileSync } from 'node:fs';
import { parseProbes, probeLayoutOf, lrColumnSample, guessProbeCols } from '../src/components/sashimi/cghProbes';
import { detect, SEG_DEFAULTS } from '../src/components/sashimi/cghSegment';
import { splitDelimited } from '../src/components/sashimi/cghMapping';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = JSON.stringify(got) === JSON.stringify(want); if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', JSON.stringify(got), ok ? '' : `| want ${JSON.stringify(want)}`); };
const [lay, po] = process.argv.slice(2);
const read = (p: string) => new Uint8Array(readFileSync(p));
// the reference: the CytoGenomics file of the same probes
const ref = parseProbes(read(`${po}/CGH_Log_Ratio_DEMO02.txt`));
const refStarts = [...ref.byChrom].map(([c, p]) => `${c}:${p.start.length}:${p.start[0]}`).join();
const want = ['chr7 loss', 'chr13 loss', 'chr15 gain', 'chrX gain'];
const cases: [string, string[], string, boolean][] = [
  ['DEMO03.cnr', ['DEMO03'], 'CNVkit (.cnr)', true],
  ['illumina_DEMO02.txt', ['DEMO02'], 'Illumina-like (Log R Ratio)', false],
  ['chas_DEMO05.txt', ['chas_DEMO05'], 'Affymetrix ChAS-like', false],
  ['two_samples.csv', ['PATA', 'PATB'], 'generic (chromosome, position, log ratio)', false],
];
for (const [file, samples, layout, zero] of cases) {
  const b = read(`${lay}/${file}`), specs = probeLayoutOf(b, `x/${file}`);
  eq(`${file}: samples`, specs?.map(s => s.sample), samples);
  eq(`${file}: layout, 0-based`, [specs?.[0].layout, specs?.[0].zeroBased], [layout, zero]);
  const p = parseProbes(b, specs![0]);
  eq(`${file}: build`, p.build, file.endsWith('.cnr') || file.startsWith('illumina') ? null : 'GRCh37');
  // same probes as the CytoGenomics file: chromosomes, counts, first starts (0-based CNVkit starts moved to 1-based; XY / MT / 0 left out)
  eq(`${file}: probes as the reference`, [...p.byChrom].map(([c, x]) => `${c}:${x.start.length}:${x.start[0]}`).join(), refStarts);
  eq(`${file}: sorted`, [...p.byChrom.values()].every(x => x.start.every((s, i) => !i || s >= x.start[i - 1])), true);
  const nan = [...p.byChrom.values()].reduce((k, x) => k + x.lr.filter(Number.isNaN).length, 0);
  eq(`${file}: missing values as NaN (~0.5 %)`, nan > 150 && nan < 500, true);
  const ev = detect(p, SEG_DEFAULTS, true).events.map(e => `${e.chrom} ${e.kind}`);
  eq(`${file}: probe calls`, ev, want);
  if (specs!.length > 1) { const q = parseProbes(b, specs![1]); eq(`${file}: 2nd sample, no event`, detect(q, SEG_DEFAULTS, true).events.length, 0); eq(`${file}: 2nd sample name`, q.sample, 'PATB'); }
}
// ChAS: Log2Ratio taken, not WeightedLog2Ratio
eq('ChAS: plain log2 ratio column', probeLayoutOf(read(`${lay}/chas_DEMO05.txt`), 'chas_DEMO05.txt')![0].lr, 3);
// odd layout: not recognised (the dialog), its columns guessed by content
const odd = read(`${lay}/DEMO06_probes.txt`);
eq('odd: not recognised', probeLayoutOf(odd, 'DEMO06_probes.txt'), null);
const ls = new TextDecoder().decode(odd).split('\n').slice(0, 300).map(l => splitDelimited(l, ';'));
eq('odd: columns guessed by content', guessProbeCols(ls[0], ls.slice(1)), { chrom: 1, start: 2, end: -1, name: -1, lr: [3] });
// tables of calls and short files are not probe files
eq('interval report: not probes', probeLayoutOf(read(`${po}/IntervalBasedReport_DEMO02.xls`), 'r.xls'), null);
const short = new TextEncoder().encode('chrom\tpos\tlog2\n' + Array.from({ length: 500 }, (_, i) => `1\t${i * 1000 + 1}\t0.1`).join('\n'));
eq('500 rows: not probes', probeLayoutOf(short, 's.txt'), null);
const callish = new TextEncoder().encode('chrom\tstart\tend\tnum.mark\tseg.mean\n' + Array.from({ length: 3000 }, (_, i) => `1\t${i * 1000 + 1}\t${i * 1000 + 900}\t5\t0.1`).join('\n'));
eq('segments (num.mark, seg.mean): not probes', probeLayoutOf(callish, 's.txt'), null);
const agil = new TextEncoder().encode('Genome build: hg38\nSample: AG01\nProbeName,ChrName,Start,Stop,LogRatio\n' + Array.from({ length: 2500 }, (_, i) => `A_${i},chr${(i % 22) + 1},${i * 1000 + 1},${i * 1000 + 60},${i % 5 ? '0.01' : 'NaN'}`).join('\n'));
const as = probeLayoutOf(agil, 'agilent.csv');
eq('Agilent-like: sample line, build, layout', [as?.[0].sample, as?.[0].build, as?.[0].layout, as?.[0].name], ['AG01', 'GRCh38', 'Agilent-like', 0]);
eq('Agilent-like: 2,500 probes, 500 NaN', (() => { const p = parseProbes(agil, as![0]); return [p.n, [...p.byChrom.values()].reduce((k, x) => k + x.lr.filter(Number.isNaN).length, 0)]; })(), [2500, 500]);
// header names
for (const [h, s] of [['Log R Ratio', ''], ['S1.Log R Ratio', 'S1'], ['log2', ''], ['Log2Ratio', ''], ['WeightedLog2Ratio', ''], ['Normalized log ratio(X 1)', 'X 1'], ['PAT-A_log2ratio', 'PAT-A'], ['logR', ''], ['LRR', ''], ['log_ratio', ''], ['Log Ratio (B)', 'B'], ['B Allele Freq', null], ['depth', null], ['gene', null]] as [string, string | null][])
  eq(`header "${h}"`, lrColumnSample(h), s);
console.log(bad ? `${bad} FAILED` : 'all ok');
process.exit(bad ? 1 : 0);
