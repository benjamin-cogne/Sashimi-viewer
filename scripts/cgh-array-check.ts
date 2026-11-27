// Unit checks of the array format (cghArray.ts), invented data only: probe counts → nearest standard format, design ids
// from header lines, and invented probe files of several densities.
// Run: for n in 60000 180000 400000 1000000 110000; do node scripts/cgh-probe-fixture.mjs /tmp/fx.zip --density $n /tmp/dens/d$n.txt; done
//      npx esbuild scripts/cgh-array-check.ts --bundle --platform=node --outfile=/tmp/cgh-array-check.cjs && node /tmp/cgh-array-check.cjs /tmp/dens
import { readFileSync } from 'fs';
import { arrayInfo, designOfMeta, nearestFormat, medianSpacing } from '../src/components/sashimi/cghArray';
import { parseProbes } from '../src/components/sashimi/cghProbes';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = JSON.stringify(got) === JSON.stringify(want); if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', JSON.stringify(got), ok ? '' : `| want ${JSON.stringify(want)}`); };
eq('nearest', [61246, 45000, 59000, 170000, 240000, 411000, 963000, 1_300_000, 110000, 2_000_000, 0].map(nearestFormat), ['60K', '60K', '60K', '180K', '180K', '400K', '1M', '1M', 'other', 'other', 'other']);
eq('design from report', designOfMeta([['Genome', 'hg19'], ['Array Design Info', 'X01[000000_20000101]']]), '000000_20000101');
eq('design from AMADID meta', designOfMeta([['AMADID', '012345']]), '012345');
eq('no design', designOfMeta([['Genome', 'hg19']]), null);
eq('metadata only: no chip', (({ chip, text }) => ({ chip, text }))(arrayInfo({ designs: [{ id: '000000_20000101', from: 'the report' }] })!), { chip: null, text: 'design 000000_20000101 (no probe file: format not deduced)' });
eq('nothing known: null', arrayInfo({ designs: [{ id: null, from: 'x' }], probes: null }), null);
const sp = new Map([['chr1', { start: [1, 101, 301, 601] }], ['chrX', { start: [1, 2, 3] }]]);
eq('median spacing (autosomes only)', medianSpacing(sp), 200);
const dir = process.argv[2];
if (dir) for (const [n, chip] of [[60000, '60K'], [180000, '180K'], [400000, '400K'], [1000000, '1M'], [110000, 'other']] as const) {
  const p = parseProbes(new Uint8Array(readFileSync(`${dir}/d${n}.txt`)));
  const a = arrayInfo({ designs: [{ id: p.design, from: "the probe file's AMADID line" }], probes: p })!;
  eq(`density ${n}: chip, design`, [a.chip, p.design], [chip, '000000']);
  console.log('     ', a.text);
}
else { console.log('FAIL no density dir given'); bad++; }
console.log(bad ? `${bad} FAILED` : 'all ok');
process.exit(bad ? 1 : 0);
