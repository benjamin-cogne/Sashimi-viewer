/**
 * Array-CGH: built-in examples, three INVENTED patients (no real sample) written as CytoGenomics interval tables
 * (GRCh37), loaded through the normal import path. The public build ships them: never put patient data here.
 *  - EXAMPLE-DEL22Q11: 22q11.2 LCR22A-D deletion, its ends those of ClinGen ISCA-37446 (HIGH);
 *  - EXAMPLE-PIEV-DUP: 1q21.1 distal BP3-BP4 duplication, AChro-Puce PIEV, ClinGen TS 2 region (MODERATE);
 *  - EXAMPLE-DER: terminal 16q24.3 loss breaking in ANKRD11 intron 4 (NM_013275.6, HI 3, KBG syndrome) and terminal
 *    9p24.3 gain: the pattern of a der(16)t(9;16), as the rearrangement hint reads it;
 *  - EXAMPLE-PROBES-ONLY: a probe file alone (no table: a probes-only patient), a 7q11.23 loss and a 15q11.2 BP1-BP2 gain.
 * Each comes with invented probe log ratios (cghExampleProbes.ts, a CGH_Log_Ratio.txt written at load time) matching its
 * calls; EXAMPLE-DEL22Q11's probes also hold a small 3p14.2 loss the table does not list (a probe call *not in the file*).
 */
import { cgh } from './cghStore';
import { exampleProbeFile, type ProbeEvent } from './cghExampleProbes';

const HEAD = 'Chromosome\tStart\tStop\tCytoband\tGene Name\tSize(kb)\tType\t#Probes\tMean Log Ratio\tISCN\tSuppress\tClassification';
type Row = [chrom: string, start: number, stop: number, band: string, genes: string, kind: 'Deletion' | 'Duplication', probes: number, logR: number];
/** an interval table as CytoGenomics exports it: header lines, then one row per call (size and ISCN written from the row) */
function table(name: string, sex: 'Male' | 'Female', rows: Row[]): string {
  const line = ([c, s, e, band, genes, kind, n, lr]: Row) => [c, s, e, band, genes, ((e - s) / 1000).toFixed(3), kind, n, lr.toFixed(3),
    `arr[hg19] ${c.replace('chr', '')}${band}(${s}_${e})x${kind === 'Deletion' ? 1 : 3}`, 'FALSE', ''].join('\t');
  return ['Genome\thg19', `Record Name\t${name}`, `Gender\t${sex}`, '', HEAD, ...rows.map(line)].join('\n') + '\n';
}
/** the table's calls as events of the probes, plus `more` (events the table does not list) */
const probesOf = (name: string, sex: 'Male' | 'Female', rows: Row[], more: ProbeEvent[] = []) => () =>
  exampleProbeFile(name, [...rows.map(([c, s, e, , , , n, lr]): ProbeEvent => [c, s, e, n, lr]), ...more], sex === 'Male');

interface Example { id: string; name: string; label: string; hint: string; text: string | null; probes: () => string }
const ex = (id: string, name: string, label: string, hint: string, sex: 'Male' | 'Female', rows: Row[], more: ProbeEvent[] = []): Example =>
  ({ id, name, label, hint, text: table(name, sex, rows), probes: probesOf(name, sex, rows, more) });

export const CGH_EXAMPLES: Example[] = [
  ex('del22q11', 'EXAMPLE-DEL22Q11', '22q11.2 deletion (LCR22A-D)', 'a typical pathogenic deletion: HIGH', 'Female', [
      ['chr1', 152556014, 152586532, 'q21.3', 'LCE3B, LCE3C', 'Deletion', 6, -0.903],
      ['chr22', 18912231, 21465672, 'q11.21', 'DGCR6, PRODH, DGCR2, TBX1, COMT, CRKL', 'Deletion', 186, -0.948],
    ], [['chr3', 60512000, 60838000, 12, -0.62]]),
  ex('piev', 'EXAMPLE-PIEV-DUP', '1q21.1 distal duplication (PIEV)', 'a susceptibility CNV (AChro-Puce PIEV): MODERATE', 'Male', [
      ['chr1', 110226519, 110251661, 'p13.3', 'GSTM1', 'Deletion', 5, -0.882],
      ['chr1', 146580123, 147390876, 'q21.1', 'CHD1L, BCL9, ACP6, GJA5, GJA8', 'Duplication', 47, 0.561],
    ]),
  ex('der', 'EXAMPLE-DER', 'terminal 16q loss + 9p gain (der)', 'an unbalanced translocation, a breakpoint inside ANKRD11', 'Female', [
      ['chr9', 204193, 1146312, 'p24.3', 'DOCK8, KANK1, DMRT1, DMRT3, DMRT2', 'Duplication', 41, 0.552],
      ['chr16', 89362118, 90155062, 'q24.3', 'ANKRD11, SPG7, CDK10, FANCA, MC1R, TUBB3', 'Deletion', 63, -0.917],
    ]),
  // probes alone: no table, so its calls are the probe calls (a probes-only patient); the sex is asked
  { id: 'probes', name: 'EXAMPLE-PROBES-ONLY', label: 'probes only: 7q11.23 loss + 15q11.2 gain', hint: 'a CGH_Log_Ratio.txt without its table: calls detected from the probes', text: null,
    probes: () => exampleProbeFile('EXAMPLE-PROBES-ONLY', [['chr7', 72744000, 74142000, 70, -0.93], ['chr15', 22770000, 23090000, 22, 0.57]], false) },
];

/** the examples' sample names, to tell them from the user's patients */
export const EXAMPLE_NAMES = CGH_EXAMPLES.map(x => x.name);
/** add examples (one id, or all) as files picked: the import question then asks as for any file */
export function loadExamples(id: string | 'all') {
  const xs = CGH_EXAMPLES.filter(x => id === 'all' || x.id === id);
  // each table with its probe file, as an export gives them (the probe file goes to the sample of its header's name)
  return cgh.add(xs.flatMap(x => [...(x.text ? [new File([x.text], `${x.name}.txt`, { type: 'text/plain' })] : []), new File([x.probes()], `CGH_Log_Ratio_${x.name}.txt`, { type: 'text/plain' })]));
}
