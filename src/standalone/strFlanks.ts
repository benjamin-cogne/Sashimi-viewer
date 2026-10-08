/**
 * Published flanking sequences of repeat loci, used to find a repeat's ends in the reads (repeatScan.ts anchors):
 * STRique's targets (configs/repeat_config.tsv of https://github.com/giesselmann/STRique, sha256 178ee119842d…, downloaded
 * 2026-10-08; Giesselmann et al. 2019, Nat Biotechnol 37:1478–1481), keyed by STRchive locus id. `prefix` ends where
 * STRique's tract starts and `suffix` begins where it ends, both in the reference orientation; a flank may begin or end
 * with repeat units STRique leaves outside its tract (FMR1's suffix starts with AGG (CGG)9): the anchors skip them.
 *
 * STRique is under the MIT License:
 * Copyright (c) 2018-2019, Pay Giesselmann, Max Planck Institute for Molecular Genetics. Permission is hereby granted,
 * free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall
 * be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY
 * OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

export interface PublishedFlanks { tool: string; name: string; motif: string; prefix: string; suffix: string }

export const STR_FLANKS: Record<string, PublishedFlanks> = {
  FTDALS1_C9orf72: { tool: 'STRique', name: 'c9orf72', motif: 'GGCCCC', prefix: 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC', suffix: 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC' },
  FXS_FMR1: { tool: 'STRique', name: 'fmr1', motif: 'CGG', prefix: 'GCGGGCCGGGGGTTCGGCCTCAGTCAGGCGCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGCG', suffix: 'AGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCCAATGGCGCTTTCTACAAG' },
};
