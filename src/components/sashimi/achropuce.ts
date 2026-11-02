/**
 * The recurrent (and some rare) CNVs of the French AChro-Puce network's catalogue, "Fiches d'aide à la classification
 * et à l'interprétation de CNVs récurrents ou rares", version 6 (November 2025):
 * https://acpa-achropuce.com/wp-content/uploads/2026/01/catalogue_2025-VF.pdf (from the network's postnatal page,
 * https://acpa-achropuce.com/diagnostic-postnatal/).
 *
 * For each sheet: the class the network gives each CNV in its conclusion table (classes 1 to 5, or PIEV: a CNV of
 * susceptibility to neurodevelopmental disorders with incomplete penetrance and variable expressivity), the table's
 * wording, a short summary written for this viewer (not the sheet's text: read the sheet), and where the CNV lies:
 * the coordinates the sheet gives (in its build; mapped to the other build through Ensembl), else the span of the genes
 * the sheet names in the region, from the gene models in view. A position the sheet gives only to the Mb, or taken from
 * genes, is approximate and drawn so.
 */
export type AchroClass = '1' | '2' | '3' | '4' | '5' | 'PIEV';
export interface AchroRow {
  type: 'del' | 'dup';
  /** for X-linked CNVs, the class differs between sexes */
  sex?: 'male' | 'female';
  cls: AchroClass;
  /** the conclusion table's wording, in English */
  label: string;
  /** a second class the table gives for one aspect (17q12 deletion: class 5 for the renal and pancreatic disease) */
  also?: { cls: AchroClass; label: string };
}
export interface AchroRegion {
  id: string;
  /** as the sheet names it */
  name: string;
  /** the sheet's title and authors */
  sheet: string; authors: string;
  /** coordinates the sheet gives: 1-based inclusive as written; `approx` when given only to the Mb or as a range of positions */
  coords?: { build: 'GRCh37' | 'GRCh38'; chrom: string; start: number; end: number; approx?: boolean; note?: string }[];
  /** genes the sheet names in the region: the span of those found places it when no coordinates fit the build */
  genes: string[];
  rows: AchroRow[];
  summary: string;
}

export const ACHRO_SOURCE = {
  title: 'AChro-Puce: Fiches d’aide à la classification et à l’interprétation de CNVs récurrents ou rares',
  version: 'version 6, November 2025',
  url: 'https://acpa-achropuce.com/wp-content/uploads/2026/01/catalogue_2025-VF.pdf',
  page: 'https://acpa-achropuce.com/diagnostic-postnatal/',
};

const susceptibility = (pen: 'high' | 'low' | 'incomplete') => `susceptibility to neurodevelopmental disorders, ${pen === 'incomplete' ? 'incomplete' : pen} penetrance and variable expressivity`;

export const ACHRO_REGIONS: AchroRegion[] = [
  {
    id: '1q21.1-prox', name: '1q21.1 proximal (BP2-BP3, TAR region)', sheet: 'CNVs en 1q21.1q21.2', authors: 'M. Egloff & M. Reynaud, 2025',
    genes: ['HJV', 'TXNIP', 'POLR3GL', 'RBM8A', 'PEX11B', 'ITGA10', 'PIAS3', 'POLR3C', 'CD160', 'PDZK1'],
    rows: [{ type: 'del', cls: '3', label: 'uncertain significance' }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 330 kb between BP2 and BP3, including RBM8A. Both the deletion and the duplication are of uncertain significance. With a TAR syndrome (thrombocytopenia, absent radius), a deletion points to a hypomorphic RBM8A allele unmasked on the other chromosome.',
  },
  {
    id: '1q21.1-dist', name: '1q21.1q21.2 distal (BP3-BP4)', sheet: 'CNVs en 1q21.1q21.2', authors: 'M. Egloff & M. Reynaud, 2025',
    genes: ['PRKAB2', 'CHD1L', 'FMO5', 'BCL9', 'ACP6', 'GJA5', 'GJA8', 'GPR89B', 'NBPF11'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('high') }, { type: 'dup', cls: 'PIEV', label: susceptibility('low') }],
    summary: 'About 1.2 Mb between BP3 and BP4, 9 coding genes including GJA5 and GJA8. Both are susceptibility CNVs (penetrance about 37 % for the deletion, 29 % for the duplication, Rosenfeld 2013, any phenotypic trait counted). Parental segregation and genetic counselling recommended.',
  },
  {
    id: '2q11.2', name: '2q11.2', sheet: 'CNVs en 2q11.2', authors: 'G. Jedraszak, 2022',
    coords: [{ build: 'GRCh37', chrom: 'chr2', start: 96739012, end: 97671429 }],
    genes: ['ADRA2B', 'LMAN2L', 'SNRNP200', 'NCAPH'],
    rows: [{ type: 'del', cls: '3', label: 'uncertain significance' }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 930 kb, 16 OMIM genes including ADRA2B, LMAN2L, SNRNP200 and NCAPH. Too few reports to conclude: uncertain significance for both.',
  },
  {
    id: '2q13', name: '2q13', sheet: 'CNVs en 2q13', authors: 'S. Redon & K. Uguen, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr2', start: 111399243, end: 113098686 }],
    genes: ['BUB1', 'ACOX1', 'BCL2L11', 'ANAPC1', 'MERTK', 'TMEM87B', 'FBLN7', 'ZC3H8', 'ZC3H6'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 1.7 Mb, 9 genes from BUB1 to ZC3H6; other sizes exist (110–114 Mb). The deletion, de novo or inherited from an unaffected or mildly affected parent, is a susceptibility CNV: parental segregation and genetic counselling recommended.',
  },
  {
    id: '3q29', name: '3q29', sheet: 'CNVs en 3q29', authors: 'M. Massier, G. Delplancq & V. Malan, 2024',
    coords: [{ build: 'GRCh38', chrom: 'chr3', start: 195999964, end: 197617792 }, { build: 'GRCh37', chrom: 'chr3', start: 195726835, end: 197344663 }],
    genes: ['TFRC', 'PCYT1A', 'RNF168', 'CEP19', 'PAK2', 'DLG1', 'BDH1'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('high') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: '1.6 Mb, 23 coding genes. Deletion (MIM #609425, ClinGen ISCA-37443): feeding difficulties, low weight, developmental delay or intellectual disability, neuropsychiatric disorders; high penetrance, very variable expressivity; parental segregation recommended. Duplication (MIM #611936): heterogeneous, often inherited; keep looking for another cause.',
  },
  {
    id: '10q11.2', name: '10q11.21q11.23 (critical region LCR2C-LCR2D)', sheet: 'CNVs en 10q11.21q11.23', authors: 'E. Pisan & V. Malan, 2024',
    genes: ['MAPK8', 'ARHGAP22', 'WDFY4', 'LRRC18', 'FAM170B', 'DRGX', 'ERCC6', 'PGBD3', 'SLC18A3', 'CHAT', 'OGDHL'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'Sizes vary (deletions 1.9–10.9 Mb) but all include the 1.5 Mb critical region between LCR2C and LCR2D (11 genes, candidates SLC18A3 and CHAT), drawn here. Deletion: susceptibility CNV, parental segregation recommended.',
  },
  {
    id: '15q11.2', name: '15q11.2 BP1-BP2', sheet: 'CNVs en 15q11.2 (BP1-BP2)', authors: 'R. Mathevet & P. Kuentz, 2025',
    coords: [{ build: 'GRCh38', chrom: 'chr15', start: 22782171, end: 23040548, note: 'ClinGen ISCA-37448' }],
    genes: ['TUBGCP5', 'CYFIP1', 'NIPA2', 'NIPA1'],
    rows: [{ type: 'del', cls: '3', label: 'uncertain significance' }, { type: 'dup', cls: '1', label: 'benign' }],
    summary: 'About 443 kb, 4 genes (NIPA1, NIPA2, CYFIP1, TUBGCP5). The deletion, very often inherited from unaffected parents and frequent in the general population (0.14–0.39 %), with very low penetrance, is of uncertain significance; the duplication is benign.',
  },
  {
    id: '15q13.3', name: '15q13.3 BP4-BP5', sheet: 'CNVs en 15q13.3 (BP4/BP5)', authors: 'G. Delplancq & S. Romana, 2024',
    genes: ['ARHGAP11B', 'MTMR10', 'FAN1', 'TRPM1', 'KLF13', 'OTUD7A', 'CHRNA7'],
    rows: [{ type: 'del', cls: 'PIEV', label: `${susceptibility('high')} (mainly neurological)` }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 1.54 Mb between BP4 and BP5, 7 coding genes; candidates OTUD7A and CHRNA7. Deletion (MIM #612001, ClinGen ISCA-37411): developmental delay or intellectual disability, language disorders, hypotonia, epilepsy; parental segregation recommended. Duplication: uncertain significance, segregation not needed.',
  },
  {
    id: '15q13.3-chrna7', name: '15q13.3 CHRNA7 (CHRNA7-LCR/BP5)', sheet: 'CNVs en 15q13.3 (BP4/BP5)', authors: 'G. Delplancq & S. Romana, 2024',
    genes: ['CHRNA7'],
    rows: [{ type: 'del', cls: 'PIEV', label: `${susceptibility('high')} (neurobehavioural)` }, { type: 'dup', cls: '1', label: 'benign' }],
    summary: 'The smaller CNV between the CHRNA7 low-copy repeat and BP5, containing CHRNA7. Deletion: susceptibility CNV; duplication: benign.',
  },
  {
    id: '16p13.11', name: '16p13.11 (interval II)', sheet: 'CNVs en 16p13.11', authors: 'C. Pebrel-Richard, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr16', start: 15480000, end: 16320000, approx: true, note: 'interval II, 15.48–16.32 Mb' }],
    genes: ['NDE1', 'MYH11', 'ABCC1'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 1.2 Mb (interval II), 33 RefSeq genes including NDE1 (candidate for the neurodevelopmental disorders) and MYH11. Deletion: susceptibility CNV, parental segregation recommended. Duplication: still controversial (aortic anomalies reported); uncertain significance, the diagnostic work-up goes on.',
  },
  {
    id: '16p12.2', name: '16p12.2', sheet: 'CNVs en 16p12.2', authors: 'A.-C. Tabet & A. Vitrac, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr16', start: 21946524, end: 22467284 }],
    genes: [],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 520 kb. The deletion is most often inherited from unaffected parents; the phenotype may be modulated by another genomic variant, worth looking for. The duplication is rare and of uncertain significance.',
  },
  {
    id: '16p11.2-dist', name: '16p11.2 distal (BP2-BP3)', sheet: 'CNVs en 16p11.2 distal (BP2-BP3)', authors: 'C. Pebrel-Richard & J.-M. Dupont, 2025',
    genes: ['TUFM', 'SH2B1', 'ATP2A1', 'CD19', 'LAT'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('high'), also: { cls: '5', label: 'obesity' } }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'About 220 kb including SH2B1. The deletion is highly penetrant (62 %, Rosenfeld 2013; 33 %, Goh 2025) with variable expressivity, over-represented in intellectual disability, autism and developmental delay, often inherited, with a major risk of obesity (class 5 for obesity). The duplication has low penetrance (risk of scoliosis).',
  },
  {
    id: '16p11.2-prox', name: '16p11.2 proximal (BP4-BP5)', sheet: 'CNVs en 16p11.2 proximal', authors: 'A.-C. Tabet & A. Vitrac, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr16', start: 29500000, end: 30100000, approx: true, note: 'BP4-BP5, 29.5–30.1 Mb' }],
    genes: ['SEZ6L2', 'PRRT2', 'KIF22', 'TAOK2', 'KCTD13', 'TBX6', 'MAPK3'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('high') }, { type: 'dup', cls: 'PIEV', label: susceptibility('high') }],
    summary: 'About 600 kb between BP4 and BP5, including TBX6. The deletion, with a proven neurodevelopmental disorder, is to be considered as contributing to the phenotype; the duplication, less frequent and less penetrant, often inherited from unaffected parents, is a susceptibility factor.',
  },
  {
    id: '17q12', name: '17q12 (RCAD region, HNF1B)', sheet: 'CNVs en 17q12', authors: 'C. Rooryck-Thambo, 2024',
    genes: ['ACACA', 'SYNRG', 'GGNBP2', 'LHX1', 'HNF1B'],
    rows: [{ type: 'del', cls: '5', label: 'renal and pancreatic disease', also: { cls: 'PIEV', label: susceptibility('low') } }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: '1.3–1.9 Mb, 14 to 19 genes including HNF1B. Deletion (RCAD, MIM 614527): renal malformation and dysfunction, MODY5 diabetes, liver involvement, developmental delay; pathogenic for the renal and pancreatic disease, a susceptibility CNV for neurodevelopment; parental testing recommended. Duplication (MIM 614526): mostly inherited, penetrance about 20 % (Goh 2025); uncertain significance.',
  },
  {
    id: '22q11.2-prox', name: '22q11.2 proximal (LCR22A-LCR22D)', sheet: 'CNVs en 22q11.2 : duplications proximales', authors: 'F. Vialard & C. Dupont, 2024',
    genes: ['PRODH', 'COMT', 'TBX1', 'HIRA', 'DGCR8', 'ZDHHC8', 'PI4KA', 'CRKL'],
    rows: [{ type: 'dup', cls: 'PIEV', label: susceptibility('low') }],
    summary: 'About 2.9 Mb between LCR22A and LCR22D, the reciprocal of the 22q11.2 (DiGeorge) deletion, including TBX1. Very variable; 73 % inherited. A susceptibility CNV that cannot by itself explain a severe phenotype: a genome-wide analysis should then be considered. (The deletion itself is not among these sheets.)',
  },
  {
    id: '22q11.2-cent', name: '22q11.2 central (LCR22B/C-LCR22D)', sheet: 'CNVs en 22q11.2 : variations centrales', authors: 'C. Missirian, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr22', start: 20731987, end: 21465672, note: 'LCR22B-D; LCR22C-D: 21,092,339–21,465,672' }],
    genes: ['ZNF74', 'SCARF2', 'KLHL22', 'MED15', 'PI4KA'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'Between LCR22B and LCR22D (734 kb) or LCR22C and LCR22D (374 kb), 14 RefSeq genes. A deletion calls for a family study and genetic counselling: a susceptibility CNV, often inherited. The duplication, with few patients reported, is of uncertain significance.',
  },
  {
    id: '22q11.2-dist-I', name: '22q11.2 distal type I (LCR22D-E)', sheet: 'CNVs en 22q11.2 : variations distales', authors: 'C. Missirian, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr22', start: 21920000, end: 22960000, approx: true, note: 'LCR22D-E, 21.92–22.96 Mb; also C-E and D-F' }],
    genes: [],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'Type I (LCR22C-E, D-E or D-F): mostly de novo; prematurity, pre- and postnatal growth delay, heart defect, learning difficulties in 40 %: a susceptibility CNV (low penetrance for cognition); in prenatal diagnosis, a specific follow-up of the pregnancy. Distal duplications: uncertain significance.',
  },
  {
    id: '22q11.2-dist-II', name: '22q11.2 distal type II (LCR22E-F)', sheet: 'CNVs en 22q11.2 : variations distales', authors: 'C. Missirian, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr22', start: 23050000, end: 23650000, approx: true, note: 'LCR22E-F, 23.05–23.65 Mb' }],
    genes: [],
    rows: [{ type: 'del', cls: '3', label: 'uncertain significance' }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'Type II (LCR22E-F): few patients reported, uncertain significance. Distal duplications: uncertain significance.',
  },
  {
    id: '22q11.2-dist-III', name: '22q11.2 distal type III (to LCR22H, SMARCB1)', sheet: 'CNVs en 22q11.2 : variations distales', authors: 'C. Missirian, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chr22', start: 23820000, end: 24995000, approx: true, note: 'LCR22F-H, 23.82–24.995 Mb; also D-H and E-H' }],
    genes: ['SMARCB1'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('low') }, { type: 'dup', cls: '3', label: 'uncertain significance' }],
    summary: 'Type III (LCR22D-H, E-H or F-H) includes SMARCB1: a susceptibility CNV (low penetrance for cognition), and a high risk of malignant rhabdoid tumour calling for surveillance. Distal duplications: uncertain significance.',
  },
  {
    id: 'Xp22.31', name: 'Xp22.31 (STS)', sheet: 'CNVs en Xp22.31', authors: 'D. Trost, 2024',
    coords: [{ build: 'GRCh37', chrom: 'chrX', start: 6455812, end: 8124954 }],
    genes: ['VCX3A', 'HDHD1', 'STS', 'VCX', 'PNPLA4', 'VCX2'],
    rows: [
      { type: 'del', sex: 'male', cls: '3', label: 'uncertain significance for neurodevelopmental disorders (X-linked ichthyosis: pathogenic)' },
      { type: 'dup', sex: 'male', cls: '2', label: 'likely benign' },
      { type: 'del', sex: 'female', cls: '2', label: 'likely benign (carrier)' },
      { type: 'dup', sex: 'female', cls: '1', label: 'benign' },
    ],
    summary: 'The deletion causes X-linked ichthyosis in boys (STS, MIM 308100); its role in intellectual disability is not established (no ID gene in the usual deletion; larger deletions more likely). Carrier girls are asymptomatic.',
  },
  {
    id: '2p25.3', name: '2p25.3 (MYT1L)', sheet: 'CNVs en 2p25.3 (rare)', authors: 'M. Bouassida, B. Hervé & F. Vialard, 2022',
    coords: [{ build: 'GRCh37', chrom: 'chr2', start: 1716437, end: 1815909, note: 'smallest overlap of the duplications (PXDN intron 1 to MYT1L intron 22); deletions vary, all include MYT1L' }],
    genes: ['MYT1L'],
    rows: [{ type: 'del', cls: '5', label: 'pathogenic' }, { type: 'dup', cls: 'PIEV', label: susceptibility('low') }],
    summary: 'Not recurrent (no segmental duplications): sizes vary, deletions all include MYT1L. Deletion: pathogenic (mild to moderate intellectual disability, obesity or overweight, behavioural disorders). Duplication: susceptibility CNV (schizophrenia, developmental delay, autism), inherited in about half of cases. Parental segregation recommended for both.',
  },
  {
    id: '2p16.3', name: '2p16.3 (NRXN1)', sheet: 'CNVs en 2p16.3 (rare)', authors: 'N. Celton, 2022',
    genes: ['NRXN1'],
    rows: [{ type: 'del', cls: 'PIEV', label: susceptibility('incomplete') }],
    summary: 'Non-recurrent, heterogeneous deletions, mostly 5′ in NRXN1 (promoter and first exons, NRXN1-α); penetrance depends on the transcripts hit. A severe phenotype should prompt a search for a biallelic NRXN1 deficiency.',
  },
];

/** class colours: pathogenic red, PIEV purple, uncertain grey, (likely) benign green */
export const ACHRO_COLOR: Record<AchroClass, string> = { '5': '#b91c1c', '4': '#ea580c', PIEV: '#7c3aed', '3': '#64748b', '2': '#22a35a', '1': '#15803d' };
export const ACHRO_CLASS_LABEL: Record<AchroClass, string> = { '5': 'class 5, pathogenic', '4': 'class 4, likely pathogenic', PIEV: 'PIEV, susceptibility (incomplete penetrance, variable expressivity)', '3': 'class 3, uncertain significance', '2': 'class 2, likely benign', '1': 'class 1, benign' };
export const achroShort = (c: AchroClass) => (c === 'PIEV' ? 'PIEV' : `cl.${c}`);
