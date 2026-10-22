<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/logo/logo-lockup-dark.svg">
  <img src="docs/logo/logo-lockup.svg" alt="Sashimi viewer" height="64">
</picture>

**Sashimi plots and read-level evidence from your own BAM/CRAM files, in one HTML file, entirely in the browser.**
No server, no installation, no upload: the files stay on your computer, so it runs on a hospital PC with patient data.

[![build](https://github.com/benjamin-cogne/Sashimi-viewer/actions/workflows/build.yml/badge.svg)](https://github.com/benjamin-cogne/Sashimi-viewer/actions/workflows/build.yml)
[![release](https://img.shields.io/github/v/release/benjamin-cogne/Sashimi-viewer?label=release)](https://github.com/benjamin-cogne/Sashimi-viewer/releases/latest)
[![license: CC BY-NC 4.0](https://img.shields.io/badge/license-CC%20BY--NC%204.0-lightgrey.svg)](LICENSE)

> [!WARNING]
> Designed by a molecular geneticist and written with Anthropic's Claude models. Check the HGVS
> nomenclature, the predicted transcript and the NMD verdict before anything from it goes into a
> clinical report. Feedback and bug reports are welcome.
> — Benjamin Cogné, Nantes University Hospital (CHU de Nantes)

## Start

1. Open **https://benjamin-cogne.github.io/Sashimi-viewer/**, or download
   [`sashimi-viewer.html`](https://github.com/benjamin-cogne/Sashimi-viewer/releases/latest/download/sashimi-viewer.html)
   from the latest release and double-click it (Edge, Chrome or Firefox).
2. Drop BAM (+ `.bai`) or CRAM (+ `.crai`) files, or their run folder, on the page. The first is the primary
   sample, the others the comparison samples (controls, parents, other patients). Each chip shows the library
   type (RNA / DNA) read from the file header; click it to change it.
3. Type a gene (`NF1`, suggested from the third letter), an ENSG id or coordinates (`chr17:31,229,000-31,231,000`)
   and press **Open**.

GRCh38 or GRCh37, read from the files' headers. Gene models (RefSeq, MANE Select) come from the UCSC API, or from
Ensembl when it fails. Reference bases come from the same APIs, or from a FASTA (+ `.fai`) you add. No data is
kept when the tab closes. The browser only remembers recent gene names, favourite GTEx tissues and, in Chromium,
a bookmark to the run folder.

## What it answers

| Clinical question | What the viewer shows |
|---|---|
| Does this intronic or synonymous VUS create a cryptic splice site or exon skipping? | A **red dashed arc** for a junction absent from the comparison samples; click it for donor/acceptor **c. positions**, the predicted transcript in **r. notation** and the reading frame; **🎬 Cartoon** gives the **NMD verdict** (55-nt rule, last-exon and start-proximal escape). |
| Would skipping this exon keep the reading frame? | A small red **frameshift sign** above every coding exon whose coding length is not a multiple of three (skipping it alone shifts the frame). The first and last coding exons carry none: skipping them removes the start or stop codon. Hover or click an exon for its coding length and codon phases. |
| Is the aberrant junction in *cis* with the variant? | **Reads** with mismatches (*reads of* the primary sample by default, any sample, or *all*), then **Collapsed**: two haplotypes per phase block, from the file's haplotags (HP/PS) or the page's read-based phasing. Each haplotype has its own junctions, and those the two use differently are flagged (Fisher's exact test). Or click the variant, **Phase with** → *⌖ Pick*, and click the arc: the junction's usage on each copy, a Fisher test and a verdict. |
| Are two variants in *cis* or in *trans*? | Click a site of the variants track, then **Phase with** `chr:pos; …`. The reads holding both sites (one read, a spliced read, two mates, a split read) give a 2×2 table, a verdict and a LOD. When no read holds both, the haplotags decide, or the reads phased between the two sites. |
| How much of the transcript is affected? | Per-sample **ψ** (rMATS-style inclusion) of the junction against its canonical alternative; **exon usage** from read depth compared across the open files (DEXSeq-style relative usage, robust z-score). |
| Is the gene on the minus strand? | With an RNA sample shown, the axis is reversed so the transcript reads 5′→3′ left to right (positions decrease to the right, unlike IGV), and the gene model carries a red **antisense gene** warning. A page of DNA samples only keeps the genomic orientation. |
| Is there intron retention or a cryptic exon? | The plot is on **genomic scale** by default: look at the intronic coverage and the **IR** pills. **Equal introns** (MISO / ggsashimi convention) is for exon-focused review. |
| My report says c.2033dup — where is that in the reads? | Type **`c.2033`** in the view's search box (next to the gene name). c./n. positions, intronic offsets (`c.288+1`), ranges (`c.123_125`) and whole variant descriptions (`c.234A>G`) move the window there and mark it. **`12`** or `exon 12` goes to that exon with its splice sites; `exons 3-5` frames the three. Both are read on the transcript drawn. |
| Is this an exon-level deletion or duplication? | Click an exon: its median depth relative to the gene's other coding exons, per sample, against the other open files. |
| Is this junction normal in some tissues? | **GTEx tissue tracks** (v10, v8 fallback; GRCh38 only) as reference splicing profiles. |
| Can I open it straight from the variant page of my interpretation tool? | **Deep links** (`#variant=NC_000017.11:g.43094464G>A&pad=100`) open the viewer on the variant ±100 bp with the variant marked; the BAM/CRAM files added afterwards appear as tracks. |
| Is this "mismatch" a known polymorphism? | **SNPs** track (dbSNP 155 common; AF ≥ 1 % by default, 0.1 or 5 %) and the **variant sites** scanned from every read of the window (★ on RNA, allele bars on DNA). Click a site for its reads, allele fraction and quality. |
| How long is this repeat expansion? | **Repeat inspector**, for the STRchive loci (a pill on the ruler) or any repeat selected with Shift+drag. Each read is sized between locus-specific anchors. The result is a histogram and a waterfall of the reads unit by unit, with the alleles called and a consensus per allele. |
| What does the protein look like afterwards? | **Splicing cartoon** (experimental): animated pre-mRNA, spliced mRNA, translation with UniProt/Pfam domains, NMD verdict, exportable as SVG/PNG. |
| Is this read-through or fusion transcript in frame? | Click the fusion or long-range arc: the reading frame from the two partners' CDS and the fused peptide across the junction; **🎬 Cartoon** shows the fused transcript, its translation, both genes' domains and whether it is degraded by NMD or makes an in-frame fusion protein. |

## What it shows

**RNA-seq**
- Sashimi plot on the MANE Select model (or any transcript): coverage, junction arcs with read counts or
  **% usage**, intron retention pills. Read-through and fusion junctions are named by gene (`GENEA exon 2 → GENEB exon 1`).
- Novel junctions in **HGVS c. / r.** notation with the reading frame; the **NMD** verdict in the splicing cartoon and for retained introns.
- **Exon usage** across the open files, **ψ** of a junction against its canonical form, **GTEx** tissues next to
  your samples.
- Stranded libraries: the orientation is inferred from the spliced reads. **Strands** draws the gene's strand
  above the baseline and the opposite strand mirrored below it.
- **Reads** with mismatches at the called sites (*Consensus*, on by default), mates joined, clipped and inserted
  bases. **Collapsed** gives two haplotypes, each with its own junctions, or consensus groups, one per haplotype ×
  isoform the reads support (short and long reads).

**Genomic DNA** (short or long reads)
- **Layers**: one colour-coded control, **C · V · M · R**, switches them for every DNA sample: coverage (on by
  default), variants, methylation and reads, always drawn in that order.
- A **variants track**: each site is a bar of its allele fraction over four quality cells (base quality or
  homopolymer, mapping quality, strand bias, read-position bias; green / amber / red), from an exact scan of the
  window in the background. Click a site for the reads' distributions. Phased files get two lanes, H1 and H2, with
  their phase sets. The track label sums up the common SNPs: heterozygous count, VAF range, *allele imbalance?* or
  *no heterozygous SNP (LOH / UPD?)*.
- Coverage drawn as the exact mean depth of every pixel, with a light band up to its highest base and a tick down to
  any base below half the mean: a narrow dropout stays visible on a whole-gene view.
- Reads and pairs, with discordant pairs coloured by orientation. The **two haplotypes** are drawn as consensus
  rows, from the file's haplotags (HP/PS from WhatsHap, LongPhase, HiPhase or DRAGEN) or from read-based phasing,
  and reads can be grouped by haplotype. Split reads are joined on one row, clipped sequences shown on a click, and
  the consensus of a clip cluster is ready for BLAT.
- Structural hints: deletions inside reads, split reads, discordant pairs, soft-clip clusters, duplicated stretches.
- **CpG methylation** of long reads (ONT, PacBio; MM / ML tags) in a panel under the coverage. A 5mC ribbon splits
  into HP1 / HP2 where both haplotagged copies are read (≥ 3 calls per CpG each, the minor ≥ 25 %, within one phase
  set) and joins elsewhere. The panel also shows the haplotypes' difference, the allele-specific stretches and the CpG
  islands (or the islands only), each island tagged with its difference to the primary (e.g. −26 %). Each read's CpGs
  are coloured in the reads track. Counted in the background at the reference's CpGs only, the modkit way (strands
  combined, primary alignments, 10th-percentile confidence filter per 100 kb block).

**Organising the work**
- One **view** per gene or locus (Shift+**Open** for a second), samples pooled into **groups**, **known variants**
  on every view, panels reordered by dragging.
- **Filters**: minimum MAPQ (all, 1, 10, 20, 30) and secondary alignments, for coverage, junctions, reads and sites.
- The other build's coordinates over the ruler (GRCh38 ⇄ GRCh37).
- **Sessions** (JSON, files by name or run folder). The **HTML export** carries each view's coverage, junctions
  and filters, plus its reads, methylation and phasing when they are shown, for readers without the BAM files.
- **Deep links** from another tool: `#variant=NC_000017.11:g.43094464G>A&label=BRCA1`. Other parameters are
  `locus`, `pad`, `gene`, `build` and `reads=1`.

[Full manual](docs/manual.md) ·
[Methods and validation](docs/algorithms.md) ·
[Deep links](docs/manual.md#variants-of-interest-and-deep-links) ·
[Deploying in a laboratory](docs/manual.md#deploying-in-a-clinical-laboratory) ·
[Troubleshooting](docs/manual.md#troubleshooting) ·
[Exported-page format](docs/embedded-format.md)

## Privacy

Files are read in the browser and never uploaded. The network serves gene models, reference bases and protein
domains (UCSC, Ensembl as fallback), common SNPs, GTEx values, gene-name suggestions and the other build's
coordinates. Requests carry a gene name, an interval or the letters typed, never reads, sample names or file names.
Details in the [manual](docs/manual.md#privacy).

## Build

```
npm ci
npm run build      # → sashimi-viewer.html, one self-contained file
```

Requirements and the layout of the sources are in the [manual](docs/manual.md#build-from-source).

## Cite

Cogné B. *Sashimi viewer: in-browser Sashimi plots from BAM/CRAM files.* Version 1.2, 2026.
https://github.com/benjamin-cogne/Sashimi-viewer — see [`CITATION.cff`](CITATION.cff).

## Licence

[CC BY-NC 4.0](LICENSE) © 2026 Benjamin Cogné, CHU Nantes. Free for non-commercial use with
attribution. Not a medical device: every observation must be checked in the alignments before it
enters a report.
