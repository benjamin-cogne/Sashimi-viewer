/**
 * Array-CGH: built-in examples, three INVENTED patients (no real sample) written as CytoGenomics interval tables
 * (GRCh37), loaded through the normal import path. The public build ships them: never put patient data here.
 *  - EXAMPLE-DEL22Q11: 22q11.2 LCR22A-D deletion, its ends those of ClinGen ISCA-37446 (HIGH);
 *  - EXAMPLE-PIEV-DUP: 1q21.1 distal BP3-BP4 duplication, AChro-Puce PIEV, ClinGen TS 2 region (MODERATE);
 *  - EXAMPLE-DER: terminal 16q24.3 loss breaking in ANKRD11 intron 4 (NM_013275.6, HI 3, KBG syndrome) and terminal
 *    9p24.3 gain: the pattern of a der(16)t(9;16), as the rearrangement hint reads it.
 */
import { cgh } from './cghStore';

const HEAD = 'Chromosome\tStart\tStop\tCytoband\tGene Name\tSize(kb)\tType\t#Probes\tMean Log Ratio\tISCN\tSuppress\tClassification';
type Row = [chrom: string, start: number, stop: number, band: string, genes: string, kind: 'Deletion' | 'Duplication', probes: number, logR: number];
/** an interval table as CytoGenomics exports it: header lines, then one row per call (size and ISCN written from the row) */
function table(name: string, sex: 'Male' | 'Female', rows: Row[]): string {
  const line = ([c, s, e, band, genes, kind, n, lr]: Row) => [c, s, e, band, genes, ((e - s) / 1000).toFixed(3), kind, n, lr.toFixed(3),
    `arr[hg19] ${c.replace('chr', '')}${band}(${s}_${e})x${kind === 'Deletion' ? 1 : 3}`, 'FALSE', ''].join('\t');
  return ['Genome\thg19', `Record Name\t${name}`, `Gender\t${sex}`, '', HEAD, ...rows.map(line)].join('\n') + '\n';
}

export const CGH_EXAMPLES: { id: string; name: string; label: string; hint: string; text: string }[] = [
  {
    id: 'del22q11', name: 'EXAMPLE-DEL22Q11', label: '22q11.2 deletion (LCR22A-D)', hint: 'a typical pathogenic deletion: HIGH',
    text: table('EXAMPLE-DEL22Q11', 'Female', [
      ['chr1', 152556014, 152586532, 'q21.3', 'LCE3B, LCE3C', 'Deletion', 6, -0.903],
      ['chr22', 18912231, 21465672, 'q11.21', 'DGCR6, PRODH, DGCR2, TBX1, COMT, CRKL', 'Deletion', 186, -0.948],
    ]),
  },
  {
    id: 'piev', name: 'EXAMPLE-PIEV-DUP', label: '1q21.1 distal duplication (PIEV)', hint: 'a susceptibility CNV (AChro-Puce PIEV): MODERATE',
    text: table('EXAMPLE-PIEV-DUP', 'Male', [
      ['chr1', 110226519, 110251661, 'p13.3', 'GSTM1', 'Deletion', 5, -0.882],
      ['chr1', 146580123, 147390876, 'q21.1', 'CHD1L, BCL9, ACP6, GJA5, GJA8', 'Duplication', 47, 0.561],
    ]),
  },
  {
    id: 'der', name: 'EXAMPLE-DER', label: 'terminal 16q loss + 9p gain (der)', hint: 'an unbalanced translocation, a breakpoint inside ANKRD11',
    text: table('EXAMPLE-DER', 'Female', [
      ['chr9', 204193, 1146312, 'p24.3', 'DOCK8, KANK1, DMRT1, DMRT3, DMRT2', 'Duplication', 41, 0.552],
      ['chr16', 89362118, 90155062, 'q24.3', 'ANKRD11, SPG7, CDK10, FANCA, MC1R, TUBB3', 'Deletion', 63, -0.917],
    ]),
  },
];

/** the examples' sample names, to tell them from the user's patients */
export const EXAMPLE_NAMES = CGH_EXAMPLES.map(x => x.name);
/** add examples (one id, or all) as files picked: the import question then asks as for any file */
export function loadExamples(id: string | 'all') {
  const xs = CGH_EXAMPLES.filter(x => id === 'all' || x.id === id);
  return cgh.add(xs.map(x => new File([x.text], `${x.name}.txt`, { type: 'text/plain' })));
}
