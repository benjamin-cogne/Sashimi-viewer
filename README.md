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
2. Drop BAM (+ `.bai`) or CRAM (+ `.crai`) files on the page. The first is the primary sample, the
   others the comparison samples (controls, parents, other patients).
3. Type a gene (`NF1`), an ENSG id or coordinates (`chr17:31,229,000-31,231,000`) and press **Open**.

GRCh38 or GRCh37. Gene models (RefSeq, MANE Select) and reference bases come from the UCSC API, or
from a FASTA you add. Nothing is kept when the tab closes.

## What it answers

| Clinical question | What the viewer shows |
|---|---|
| Does this intronic or synonymous VUS create a cryptic splice site or exon skipping? | A **red dashed arc** for a junction absent from the comparison samples; click it for donor/acceptor **c. positions**, the predicted transcript in **r. notation**, the… |
| Would skipping this exon keep the reading frame? | A small red **frameshift sign** above every coding exon whose coding length is not a multiple of three (skipping it alone shifts the frame); the first and last coding… |
| Is the aberrant junction in *cis* with the variant? | **Reads track** with mismatches (primary sample by default, any sample, or *All samples* for one reads track under each coverage track), then **Collapse**: read-based… |
| How much of the transcript is affected? | Per-sample **ψ** (rMATS-style inclusion) of the junction against its canonical alternative; **exon usage** from read depth compared across the open files… |
| Is the gene on the minus strand? | The axis is reversed so the transcript reads 5′→3′ left to right (positions decrease to the right, unlike IGV); the transcript track then carries a **red antisense… |
| Is there intron retention or a cryptic exon? | Switch from **equal introns** (exon-focused review, MISO / ggsashimi convention) to **genomic scale** and look at the coverage. |
| My report says c.2033dup — where is that in the reads? | Type **`c.2033`** in the view's search box (next to the gene name): c./n. positions, intronic offsets (`c.288+1`), ranges (`c.123_125`) and whole variant descriptions (`c.234A>G`) move the window there and mark it. Typing **`12`** or `exon 12` goes to that exon with its splice sites; `exons 3-5` frames the three. Both are read on the transcript drawn, so they follow the model you chose. |
| Is this an exon-level deletion or duplication? | Exon usage panel: median depth of each coding MANE exon relative to the other exons of the gene, per sample. |
| Is this junction normal in some tissues? | **GTEx tissue tracks** (v10, v8 fallback, hg38) as reference splicing profiles. |
| Can I open it straight from the variant page of my interpretation tool? | **Deep links** (`#variant=NC_000017.11:g.43094464G>A&pad=100`) open the browser on the variant ±100 bp with the variant marked; the BAM/CRAM files added afterwards… |
| Is this "mismatch" a known polymorphism? | **Common SNPs** track (dbSNP 155, MAF ≥ 1 %) and **Variant sites (★)** called from the reads of the window. |
| What does the protein look like afterwards? | **Splicing cartoon** (experimental): animated pre-mRNA, spliced mRNA, translation with UniProt/Pfam domains, NMD verdict, exportable as SVG/PNG. |

## What it shows

**RNA-seq**
- Sashimi plot on the MANE Select model: coverage, junction arcs with read counts or **% usage** per
  intron, intron retention pills.
- Novel junctions in **HGVS c. / r.** notation, reading frame and **NMD** verdict; a frameshift sign
  above every exon whose skipping shifts the frame.
- **Exon usage** across the open files, **ψ** of a junction against its canonical form, **GTEx**
  tissues next to your samples.
- **Reads track** with mismatches, mates joined, clipped and inserted bases; **Collapse** into two
  phased haplotypes or consensus groups.

**Genomic DNA** (short or long reads)
- Coverage with **variant sites** as allele bars (option off by default), from the reads track or
  from an exact scan of the window run in the background; allele balance of the common SNPs.
- Coverage drawn as the exact mean depth of every pixel, with a tick down to any base below half of
  it: a narrow dropout stays visible on a whole-gene view.
- Reads, pairs, the **two haplotypes as consensus rows**, from the file's haplotags (HP/PS from
  WhatsHap, LongPhase, HiPhase or DRAGEN) or from read-based phasing; reads grouped by haplotype;
  split reads joined on one row, clipped sequences on a click,
  consensus of a clip cluster ready for BLAT.
- **CpG methylation** of long reads (ONT, PacBio; MM / ML tags) in a panel under the coverage: a 5mC
  density ribbon that splits into the two haplotypes where the reads are phased and joins where they
  are not, their difference and the allele-specific stretches, the CpG islands (or the islands only),
  the mean island difference of each sample with the primary (e.g. −26 %); each read's
  CpGs coloured in the reads track. Counted in the background at the reference's CpGs only, the
  modkit way (strands combined, 10th-percentile confidence filter).

**Organising the work**
- Several genes as **views**, samples pooled into **groups**, **known variants** drawn on every view.
- **Sessions** (JSON, files by name or run folder) and an **HTML export** that embeds the data of
  every view, reads included, for readers without the BAM files.
- **Deep links** from another tool: `#variant=NC_000017.11:g.43094464G>A&label=BRCA1`.

[Full manual](docs/manual.md) ·
[Deploying in a laboratory](docs/manual.md#deploying-in-a-clinical-laboratory) ·
[Troubleshooting](docs/manual.md#troubleshooting) ·
[Exported-page format](docs/embedded-format.md)

## Privacy

Files are read in the browser and never uploaded. The network serves gene models, reference bases,
common SNPs and GTEx values only, never your reads. Details in the [manual](docs/manual.md#privacy).

## Build

```
npm install
npm run build      # → sashimi-viewer.html, one self-contained file
```

Tests and the layout of the sources are in the [manual](docs/manual.md#build-from-source).

## Cite

Cogné B. *Sashimi viewer: in-browser Sashimi plots from BAM/CRAM files.* Version 1.1, 2026.
https://github.com/benjamin-cogne/Sashimi-viewer — see [`CITATION.cff`](CITATION.cff).

## Licence

[CC BY-NC 4.0](LICENSE) © 2026 Benjamin Cogné, CHU Nantes. Free for non-commercial use with
attribution. Not a medical device: every observation must be checked in the alignments before it
enters a report.
