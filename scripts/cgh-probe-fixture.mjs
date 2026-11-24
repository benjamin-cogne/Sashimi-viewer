// Writes an invented array-CGH export (zip: IntervalBasedReport.xls + CGH_Log_Ratio.txt) with known events, to test
// the probe calls (cghSegment.ts). Invented data only: a fake sample, about 60k probes, Gaussian noise (SD 0.13).
// Usage: node scripts/cgh-probe-fixture.mjs out.zip [--female]
import { writeFileSync } from 'node:fs';
import { zipSync, strToU8 } from 'fflate';

const out = process.argv[2] ?? 'cgh-probe-fixture.zip', female = process.argv.includes('--female');
let seed = 20261124;
const rnd = () => { seed = (seed + 0x6d2b79f5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
const gauss = () => Math.sqrt(-2 * Math.log(rnd() || 1e-12)) * Math.cos(2 * Math.PI * rnd());
// GRCh37 lengths, Mb
const LEN = { 1: 249, 2: 243, 3: 198, 4: 191, 5: 181, 6: 171, 7: 159, 8: 146, 9: 141, 10: 136, 11: 135, 12: 134, 13: 115, 14: 107, 15: 103, 16: 90, 17: 81, 18: 78, 19: 59, 20: 63, 21: 48, 22: 51, X: 155, Y: 59 };
const STEP = 50_000;
// the known events: chromosome, first probe index on it, probes, log ratio added, in the file's calls or not
export const EVENTS = [
  { chr: '7', at: 1400, n: 5, lr: -0.9, inFile: true },
  { chr: '15', at: 600, n: 200, lr: 0.55, inFile: true },
  { chr: '13', at: 900, n: 30, lr: -0.4, inFile: false },
  { chr: 'X', at: 2500, n: 30, lr: 0.9, inFile: false },
];
const S = 'DEMO02';
const lines = ['AMADID\t000000\thg19', `Probe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(${S})\t`];
const span = {};
let k = 0;
for (const [chr, mb] of Object.entries(LEN)) {
  const n = Math.floor((mb * 1e6 - 2e6) / STEP);
  for (let i = 0; i < n; i++) {
    const s = 1_000_000 + i * STEP + Math.floor(rnd() * 20_000);
    let v = gauss() * 0.13;
    if (chr === 'Y' && female) v = gauss() * 0.5;
    for (const e of EVENTS) if (e.chr === chr && i >= e.at && i < e.at + e.n) { v += e.lr; (span[EVENTS.indexOf(e)] ??= []).push(s); }
    k++;
    lines.push(`A_00_P${String(k).padStart(6, '0')}\t${chr}\t${s}\t${s + 59}\t${k}\tUnknown\t\t${rnd() < 0.01 ? 'NaN' : v.toFixed(6)}`);
  }
}
const H = 'AberrationNo\tChr\tCytoband\tStart\tStop\t#Probes\tAmplification\tGain\tLoss\tDeletion\tpval\tGene Names';
const rows = EVENTS.filter(e => e.inFile).map((e, q) => {
  const p = span[EVENTS.indexOf(e)], v = ['0.000000', '0.000000', '0.000000', '0.000000'];
  v[e.lr > 0 ? 1 : 2] = e.lr.toFixed(6);
  return `${q + 1}\tchr${e.chr}\tq11\t${p[0]}\t${p[p.length - 1] + 59}\t${p.length}\t${v.join('\t')}\t1.0e-12\t`;
});
// a one-probe call the probes do not support (the probe calls need 3)
rows.push(`${rows.length + 1}\tchr3\tp21\t60000000\t60000059\t1\t0.000000\t0.000000\t-1.200000\t0.000000\t1.0e-6\t`);
const rep = ['Genome: hg19', `Array Design Info: ${S}[000000_20000101]`, '', H, S, ...rows, '', 'Number of calls per sample:', `${S} =${rows.length}`, ''];
writeFileSync(out, zipSync({ 'demo2/IntervalBasedReport.xls': strToU8(rep.join('\r\n')), 'demo2/CGH_Log_Ratio.txt': strToU8(lines.join('\n') + '\n') }));
console.log(`${out}: ${k} probes; events`, EVENTS.map((e, q) => `${e.chr}:${span[q][0]}-${span[q][span[q].length - 1]} ${e.n}p ${e.lr}`).join(', '));
