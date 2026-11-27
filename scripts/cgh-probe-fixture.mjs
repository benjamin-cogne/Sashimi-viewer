// Writes an invented array-CGH export (zip: IntervalBasedReport.xls + CGH_Log_Ratio.txt) with known events, to test
// the probe calls (cghSegment.ts). Invented data only: a fake sample, about 60k probes, Gaussian noise (SD 0.13).
// Usage: node scripts/cgh-probe-fixture.mjs out.zip [--female] [--split DIR] [--layouts DIR] [--hg38 DIR] [--density N FILE]
// --split DIR also writes, for the probes-only mode, DIR/CGH_Log_Ratio_DEMO02.txt (the probe file alone),
// DIR/demo2-probes-only.zip (the probe file and two dummy files, no table) and DIR/IntervalBasedReport_DEMO02.xls (the table alone).
// --layouts DIR writes the same invented probes (same events) in other probe-level layouts (cghProbes.ts `probeLayoutOf`):
// DEMO03.cnr (CNVkit, 0-based starts), illumina_DEMO02.txt (Illumina-like final report, `DEMO02.Log R Ratio`: the table's
// sample), chas_DEMO05.txt (ChAS-like, `#%` lines, chromosomes 23 / 24, NA), two_samples.csv (PATA / PATB columns,
// unsorted, empty and null cells) and DEMO06_probes.txt (an odd layout: `;`, decimal commas, no known column name: the dialog).
import { writeFileSync, mkdirSync } from 'node:fs';
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
const probes = [];
for (const [chr, mb] of Object.entries(LEN)) {
  const n = Math.floor((mb * 1e6 - 2e6) / STEP);
  for (let i = 0; i < n; i++) {
    const s = 1_000_000 + i * STEP + Math.floor(rnd() * 20_000);
    let v = gauss() * 0.13;
    if (chr === 'Y' && female) v = gauss() * 0.5;
    for (const e of EVENTS) if (e.chr === chr && i >= e.at && i < e.at + e.n) { v += e.lr; (span[EVENTS.indexOf(e)] ??= []).push(s); }
    k++;
    probes.push({ chr, s, v, nan: false });
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
const split = process.argv.indexOf('--split');
if (split > 0) {
  const dir = process.argv[split + 1], probes = strToU8(lines.join('\n') + '\n');
  mkdirSync(dir, { recursive: true });
  writeFileSync(`${dir}/CGH_Log_Ratio_${S}.txt`, probes);
  writeFileSync(`${dir}/demo2-probes-only.zip`, zipSync({ 'demo2/CGH_Log_Ratio.txt': probes, 'demo2/QC_Metrics.xml': strToU8('<qc><dlrs>0.13</dlrs></qc>'), 'demo2/plot.png': new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10]) }));
  writeFileSync(`${dir}/IntervalBasedReport_${S}.xls`, strToU8(rep.join('\r\n')));
  console.log(`${dir}: probe file alone, zip without a table, table alone`);
}
const lay = process.argv.indexOf('--layouts');
if (lay > 0) {
  const dir = process.argv[lay + 1];
  mkdirSync(dir, { recursive: true });
  // the same probes, a little noise added for the second sample; 0.5 % of the log ratios left out (written per layout)
  const pr = probes.map((p, i) => ({ ...p, i, miss: rnd() < 0.005, v2: gauss() * 0.13 }));
  const f = (v, d = 4) => v.toFixed(d);
  const w = (name, ls) => writeFileSync(`${dir}/${name}`, ls.join('\n') + '\n');
  w('DEMO03.cnr', ['chromosome\tstart\tend\tgene\tdepth\tlog2\tweight', ...pr.map(p => `chr${p.chr}\t${p.s - 1}\t${p.s + 59}\t-\t${f(80 + rnd() * 40, 2)}\t${p.miss ? '' : f(p.v)}\t0.9`)]);
  w('illumina_DEMO02.txt', ['[Header]', 'GSGT Version\t2.0.5', 'Processing Date\t1/1/2000 0:00 AM', 'Num SNPs\t' + pr.length, '[Data]',
    'Name\tChr\tPosition\tDEMO02.GType\tDEMO02.B Allele Freq\tDEMO02.Log R Ratio', ...pr.map(p => `rs${100000 + p.i}\t${p.chr}\t${p.s}\tAB\t0.5\t${p.miss ? 'NaN' : f(p.v)}`),
    'rs1\tXY\t1000\tAB\t0.5\t0.1', 'rs2\tMT\t100\tAB\t0.5\t0.1', 'rs3\t0\t0\t--\t0\tNaN']);
  w('chas_DEMO05.txt', ['#%chip_type=CytoScanHD_Array', '#%genome-version=hg19', 'ProbeSetName\tChromosome\tPosition\tLog2Ratio\tWeightedLog2Ratio\tSmoothSignal',
    ...pr.map(p => `C-${p.i}\t${p.chr === 'X' ? 23 : p.chr === 'Y' ? 24 : p.chr}\t${p.s}\t${p.miss ? 'NA' : f(p.v)}\t${f(p.v * 0.8)}\t2.0`)]);
  const sh = [...pr]; for (let i = sh.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [sh[i], sh[j]] = [sh[j], sh[i]]; }
  w('two_samples.csv', ['# invented, two samples', 'genome: GRCh37', 'chrom,pos,PATA.log2ratio,PATB.log2ratio', ...sh.map(p => `chr${p.chr},${p.s},${p.miss ? 'null' : f(p.v)},${p.i % 97 ? f(p.v2) : ''}`)]);
  w('DEMO06_probes.txt', ['Feature;Seq;Locus_bp;Signal', ...pr.map(p => `F${p.i};chr${p.chr};${p.s};${p.miss ? 'NA' : f(p.v).replace('.', ',')}`)]);
  console.log(`${dir}: ${pr.length} probes in five other layouts`);
}
// --hg38 DIR: an invented GRCh38 male (AMADID line hg38): CGH_Log_Ratio_DEMO38.txt and IntervalBasedReport_DEMO38.xls (Genome: hg38);
// events: chr22 loss (40 probes, -0.9), a chrX PAR1 gain (25 probes, +0.55; chrY PAR1 probes carry it too: never a chrY call),
// a chrX gain outside the PARs (30 probes, +0.9)
const h38 = process.argv.indexOf('--hg38');
if (h38 > 0) {
  const dir = process.argv[h38 + 1];
  mkdirSync(dir, { recursive: true });
  const L38 = { ...LEN, X: 156, Y: 57 };
  const E38 = [{ chr: '22', at: 400, n: 40, lr: -0.9 }, { chr: 'X', at: 5, n: 25, lr: 0.55 }, { chr: 'Y', at: 5, n: 25, lr: 0.55, par: true }, { chr: 'X', at: 2500, n: 30, lr: 0.9 }];
  const out = ['AMADID\t000000\thg38', 'Probe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(DEMO38)\t'], sp = E38.map(() => []);
  let q = 0;
  for (const [chr, mb] of Object.entries(L38)) {
    const n = Math.floor((mb * 1e6 - 2e6) / STEP);
    for (let i = 0; i < n; i++) {
      const s = 1_000_000 + i * STEP + Math.floor(rnd() * 20_000);
      let v = gauss() * 0.13;
      E38.forEach((e, j) => { if (e.chr === chr && i >= e.at && i < e.at + e.n) { v += e.lr; sp[j].push(s); } });
      q++;
      out.push(`A_38_P${String(q).padStart(6, '0')}\t${chr}\t${s}\t${s + 59}\t${q}\tUnknown\t\t${v.toFixed(6)}`);
    }
  }
  writeFileSync(`${dir}/CGH_Log_Ratio_DEMO38.txt`, out.join('\n') + '\n');
  const rows = E38.map((e, j) => ({ e, p: sp[j] })).filter(x => !x.e.par).map(({ e, p }, j) => {
    const v = ['0.000000', '0.000000', '0.000000', '0.000000']; v[e.lr > 0 ? 1 : 2] = e.lr.toFixed(6);
    return `${j + 1}\tchr${e.chr}\tq11\t${p[0]}\t${p[p.length - 1] + 59}\t${p.length}\t${v.join('\t')}\t1.0e-12\t`;
  });
  writeFileSync(`${dir}/IntervalBasedReport_DEMO38.xls`, ['Genome: hg38', 'Array Design Info: DEMO38[000000_20000101]', '', H, 'DEMO38', ...rows, '', 'Number of calls per sample:', `DEMO38 =${rows.length}`, ''].join('\r\n'));
  console.log(`${dir}: GRCh38 male, ${q} probes; events`, E38.map((e, j) => `${e.chr}:${sp[j][0]}-${sp[j][sp[j].length - 1]}`).join(', '));
}
// --density N FILE: a probe file alone (no event, noise only) with about N probes spread over the genome (GRCh37 lengths),
// to test the array format deduced from the probe count and spacing (cghArray.ts): e.g. 60000, 180000, 1000000.
const dens = process.argv.indexOf('--density');
if (dens > 0) {
  const want = Number(process.argv[dens + 1]), file = process.argv[dens + 2];
  const total = Object.values(LEN).reduce((a, mb) => a + mb * 1e6 - 2e6, 0), step = total / want;
  const out = ['AMADID\t000000\thg19', 'Probe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(DEMO-DENSITY)\t'];
  let q = 0;
  for (const [chr, mb] of Object.entries(LEN)) {
    const n = Math.floor((mb * 1e6 - 2e6) / step);
    for (let i = 0; i < n; i++) {
      const s = 1_000_000 + Math.floor(i * step + rnd() * step * 0.4);
      q++;
      out.push(`A_00_D${String(q).padStart(7, '0')}\t${chr}\t${s}\t${s + 59}\t${q}\tUnknown\t\t${(gauss() * 0.13).toFixed(4)}`);
    }
  }
  writeFileSync(file, out.join('\n') + '\n');
  console.log(`${file}: ${q} probes, step ${Math.round(step)} bp`);
}
