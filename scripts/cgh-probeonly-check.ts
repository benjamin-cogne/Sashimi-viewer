// Unit checks of the probes-only mode (cghProbeOnly.ts) on the invented probe file of scripts/cgh-probe-fixture.mjs.
// Run: node scripts/cgh-probe-fixture.mjs /tmp/f.zip --split /tmp/po && npx esbuild scripts/cgh-probeonly-check.ts --bundle --platform=node --outfile=/tmp/cgh-probeonly-check.cjs && node /tmp/cgh-probeonly-check.cjs /tmp/po
import { readFileSync } from 'node:fs';
import { parseProbes, probeSampleOf } from '../src/components/sashimi/cghProbes';
import { detect, SEG_DEFAULTS } from '../src/components/sashimi/cghSegment';
import { probeOnlyFile, probeOnlyCnvs, probeBuildOf, exportProbeCalls } from '../src/components/sashimi/cghProbeOnly';
import { decodeText } from '../src/components/sashimi/cghFormats';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = got === want; if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', got, ok ? '' : `| want ${want}`); };
const bytes = new Uint8Array(readFileSync(`${process.argv[2]}/CGH_Log_Ratio_DEMO02.txt`));
eq('sample from the header', probeSampleOf(bytes), 'DEMO02');
eq('build from the AMADID line', probeBuildOf(decodeText(bytes.subarray(0, 8192))), 'GRCh37');
eq('no build', probeBuildOf('AMADID\t000000\nProbe Name\tChrName'), null);
eq('hg38', probeBuildOf('AMADID\t000000\thg38\nProbe Name\tChrName'), 'GRCh38');
const f = probeOnlyFile('x › demo2/CGH_Log_Ratio.txt', 'DEMO02', 'GRCh37');
eq('format', f.format, 'probes'); eq('export name', f.saveAs, 'DEMO02.probe-calls.txt'); eq('one sample', f.samples.length, 1);
eq('name from the file name', probeOnlyFile('CGH_Log_Ratio_X.txt', '', null).samples[0].name, 'CGH_Log_Ratio_X');
const male = detect(parseProbes(bytes), SEG_DEFAULTS, true);
const cs = probeOnlyCnvs(male.events, f.samples[0].key, 'DEMO02');
const want = [['chr7', 'loss', 5], ['chr13', 'loss', 30], ['chr15', 'gain', 200], ['chrX', 'gain', 30]];
eq('events (male)', cs.length, 4);
// probes counted without the NaN ones (1 % of the fixture's)
for (const [chr, kind, n] of want) { const c = cs.find(x => x.chrom === chr); eq(`${chr} ${kind} ~${n}p`, c && c.kind === kind && c.probes! >= +n * 0.95 && c.probes! <= +n, true); }
eq('ids by position', cs[0].id.startsWith(`pc|${f.samples[0].key}|chr7:`), true);
eq('same ids recomputed', probeOnlyCnvs(male.events, f.samples[0].key, 'DEMO02').map(c => c.id).join(), cs.map(c => c.id).join());
const fx = { ...f, samples: [{ ...f.samples[0], cnvs: cs }] };
const out = exportProbeCalls(fx, SEG_DEFAULTS, new Date('2026-11-25T10:00:00Z'), { band: () => 'q11', copies: c => (c.kind === 'loss' ? 'x1' : 'x3'), iscn: () => 'arr', cls: c => (c.chrom === 'chr7' ? 'Pathogenic' : ''), note: c => (c.chrom === 'chr7' ? 'see\tthis' : '') });
const ls = out.trimEnd().split('\n');
eq('header: source', ls[0], '# Source probe file: CGH_Log_Ratio.txt');
eq('header: build', ls[2], '# Build: GRCh37');
eq('header: settings', ls[3], '# Detection: CBS, sensitivity normal (t = 8), probes >= 3, |mean log ratio| >= 0.25');
eq('header: date', ls[4], '# Date: 2026-11-25 10:00:00');
eq('columns', ls[6], 'Chromosome\tStart\tStop\tCytoband\tSize(kb)\tType\t#Probes\tMean Log Ratio\tCopies\tISCN\tClassification\tComment');
eq('rows', ls.length, 7 + 4);
const r7 = ls[7].split('\t');
eq('row: 12 cells', r7.length, 12); eq('row: chr7 class', r7[10], 'Pathogenic'); eq('row: note flattened', r7[11], 'see this'); eq('row: type', r7[5], 'Loss'); eq('row: probes', r7[6], '5');
console.log(bad ? `${bad} FAILED` : 'all ok'); process.exit(bad ? 1 : 0);
