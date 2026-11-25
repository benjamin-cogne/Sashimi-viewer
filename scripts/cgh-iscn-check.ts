// Unit checks of the array ISCN (cghIscn.ts) and the export's ISCN column, on invented calls.
// Run: npx esbuild scripts/cgh-iscn-check.ts --bundle --platform=node --outfile=/tmp/cgh-iscn-check.cjs && node /tmp/cgh-iscn-check.cjs
import { callIscn, patientIscn, iscnBands, bandSpan, iscnBandOf } from '../src/components/sashimi/cghIscn';
import { exportCghFile } from '../src/components/sashimi/cghIntervals';
const base = { id: 'x', sample: 'S', sampleKey: 'S', cytoband: '', genes: [], sizeKb: 1, typeText: '', probes: 10, iscn: null, copies: null, copiesHi: null, amp: false, suppressed: false, classification: '', classRest: '', note: '', line: 0 };
const mk = (o: any) => ({ ...base, ...o });
let bad = 0;
const eq = (name: string, got: any, want: any) => { const ok = got === want; if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', got, ok ? '' : `| want ${want}`); };
// p-arm two bands, file copies
const a = mk({ chrom: 'chr1', start: 25889422, end: 28493687, kind: 'gain', copies: 3, logRatio: 0.58 });
eq('two bands', callIscn(a, 'GRCh38', 'p36.11 - p35.3', undefined).text, 'arr[GRCh38] 1p36.11p35.3(25,889,422_28,493,687)x3');
// single band, copies from the log ratio (no file copies), homogeneous
const b = mk({ chrom: 'chr22', start: 18894835, end: 21505417, kind: 'loss', logRatio: -0.95 });
eq('single band inferred x1', callIscn(b, 'GRCh37', 'q11.21 - q11.21', undefined).text, 'arr[GRCh37] 22q11.21(18,894,835_21,505,417)x1');
// mosaic: log ratio -0.4 on an autosome: f = 2 - 2*2^-0.4 = 0.48
const m = mk({ chrom: 'chr11', start: 41444865, end: 49146935, kind: 'loss', logRatio: -0.4 });
eq('mosaic', callIscn(m, 'GRCh38', 'p12p11.12', undefined).text, 'arr[GRCh38] 11p12p11.12(41,444,865_49,146,935)x1[0.48]');
// X in a male: loss to x0, gain x2
const x = mk({ chrom: 'chrX', start: 10000000, end: 12000000, kind: 'gain', logRatio: 0.95 });
eq('male X gain', callIscn(x, 'GRCh37', 'p22.2', 'male').text, 'arr[GRCh37] Xp22.2(10,000,000_12,000,000)x2');
// homozygous loss, x4
eq('x0', callIscn(mk({ chrom: 'chr2', start: 1000, end: 50000, kind: 'loss', logRatio: -3 }), 'GRCh37', 'p25.3', undefined).text, 'arr[GRCh37] 2p25.3(1,000_50,000)x0');
eq('x4', callIscn(mk({ chrom: 'chr2', start: 1000, end: 50000, kind: 'gain', logRatio: 1.0 }), 'GRCh37', 'p25.3', undefined).text, 'arr[GRCh37] 2p25.3(1,000_50,000)x4');
// whole chromosome
eq('whole chr', callIscn(mk({ chrom: 'chr21', start: 9411193, end: 48119895, kind: 'gain', logRatio: 0.58 }), 'GRCh37', 'q11.1q22.3', undefined, true).text, 'arr[GRCh37] (21)x3');
// file ISCN, same and different
const f1 = mk({ chrom: 'chr22', start: 18894835, end: 21505417, kind: 'loss', logRatio: -0.95, copies: 1, iscn: 'arr[hg19] 22q11.21(18894835_21505417)x1' });
eq('file same', callIscn(f1, 'GRCh37', 'q11.21', undefined).differs, false);
eq('file differs', callIscn({ ...f1, iscn: 'arr[hg19] 22q11.21(18894835_21505418)x1' }, 'GRCh37', 'q11.21', undefined).differs, true);
// no band: null
eq('no band', callIscn(b, 'GRCh37', '', undefined).text, null);
// patient: ordered by chromosome then position, X after 22
const pa = [x, m, a, b].map((c, i) => ({ c, part: callIscn(c, 'GRCh37', ['p22.2', 'p12p11.12', 'p36.11p35.3', 'q11.21'][i], 'male').part }));
eq('patient', patientIscn('GRCh37', pa), 'arr[GRCh37] 1p36.11p35.3(25,889,422_28,493,687)x3,11p12p11.12(41,444,865_49,146,935)x1[0.48],22q11.21(18,894,835_21,505,417)x1,Xp22.2(10,000,000_12,000,000)x2');
eq('patient none', patientIscn('GRCh37', []), null);
eq('bands', iscnBands('q11.21 - q11.23'), 'q11.21q11.23');
// export: a table without ISCN gets one; one with ISCN keeps its own
const mkFile = (head: string, row: string, cnv: any) => ({ id: 'f', name: 'f.txt', format: 'interval', build: 'GRCh37', lines: [head, row], eol: '\n', classCol: null, classSuffix: '', headerLine: 0, commentCol: null, delim: '\t', samples: [{ key: 'S', name: 'S', cnvs: [{ ...cnv, line: 1 }] }] } as any);
const e1 = exportCghFile(mkFile('Chr\tStart\tStop', 'chr22\t18894835\t21505417', b), () => 'Pathogenic', () => '', () => 'ISCNX');
eq('export adds ISCN', e1, 'Chr\tStart\tStop\tClassification\tComment\tISCN (Sashimi viewer)\nchr22\t18894835\t21505417\tPathogenic\t\tISCNX');
const e2 = exportCghFile(mkFile('Chr\tStart\tStop\tISCN', 'chr22\t18894835\t21505417\tarr', f1), () => 'Pathogenic', () => '', () => 'ISCNX');
eq('export keeps ISCN', e2, 'Chr\tStart\tStop\tISCN\tClassification\tComment\nchr22\t18894835\t21505417\tarr\tPathogenic\t');
// the band at sub-band precision from the cytogenetic bands (invented bands, 0-based half-open), the file's as the fallback
const bands7 = [['p11.1', 57.0, 58.0], ['q11.1', 58.0, 61.7], ['q11.21', 61.7, 67.0], ['q11.22', 67.0, 72.2], ['q11.23', 72.2, 77.5], ['q21.11', 77.5, 86.4]].map(([name, a, z]) => ({ name: name as string, start: (a as number) * 1e6, end: (z as number) * 1e6 }));
const w = mk({ chrom: 'chr7', start: 72700001, end: 74100000, kind: 'loss', copies: 1, cytoband: 'q11' });
eq('band from cytobands', iscnBandOf(w, bandSpan(bands7, w.start - 1, w.end)), 'q11.23');
eq('ISCN sub-band, single once', callIscn(w, 'GRCh37', iscnBandOf(w, bandSpan(bands7, w.start - 1, w.end)), undefined).text, 'arr[GRCh37] 7q11.23(72,700,001_74,100,000)x1');
const w2 = mk({ chrom: 'chr7', start: 65000000, end: 74100000, kind: 'gain', copies: 3, cytoband: 'q11' });
eq('ISCN two sub-bands pter to qter', callIscn(w2, 'GRCh37', iscnBandOf(w2, bandSpan(bands7, w2.start - 1, w2.end)), undefined).text, 'arr[GRCh37] 7q11.21q11.23(65,000,000_74,100,000)x3');
const w3 = mk({ chrom: 'chr7', start: 57500000, end: 62000000, kind: 'gain', copies: 3, cytoband: 'p11 - q11' });
eq('ISCN across the centromere', callIscn(w3, 'GRCh37', iscnBandOf(w3, bandSpan(bands7, w3.start - 1, w3.end)), undefined).text, 'arr[GRCh37] 7p11.1q11.21(57,500,000_62,000,000)x3');
eq('fallback: no cytobands, the file band', callIscn(w, 'GRCh37', iscnBandOf(w, undefined), undefined).text, 'arr[GRCh37] 7q11(72,700,001_74,100,000)x1');
eq('outside the bands: the file band', iscnBandOf(w, bandSpan(bands7, 200e6, 201e6)), 'q11');
console.log(bad ? `${bad} FAILED` : 'all passed');
