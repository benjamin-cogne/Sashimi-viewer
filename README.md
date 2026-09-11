<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/logo/logo-lockup-dark.svg">
  <img src="docs/logo/logo-lockup.svg" alt="Sashimi viewer" height="64">
</picture>

**Sashimi plots from your own RNA-seq BAM/CRAM files, in one HTML file, entirely in the browser.**
No server, no Python, no installation, no upload: reads are decoded on your machine, so the viewer
can be used on a hospital PC with patient data.

[![build](https://github.com/benjamin-cogne/Sashimi-viewer/actions/workflows/build.yml/badge.svg)](https://github.com/benjamin-cogne/Sashimi-viewer/actions/workflows/build.yml)
[![release](https://img.shields.io/github/v/release/benjamin-cogne/Sashimi-viewer?label=release)](https://github.com/benjamin-cogne/Sashimi-viewer/releases/latest)
[![license: CC BY-NC 4.0](https://img.shields.io/badge/license-CC%20BY--NC%204.0-lightgrey.svg)](LICENSE)

Built for clinical geneticists and for bioinformaticians who need to look at a splice variant, who want a reviewable, reference-anchored view (MANE Select, HGVS c./r. notation, NMD prediction).

---

## Get started in one minute

1. **Get the viewer** (one file, about 1 MB):
   - [Download the latest release](https://github.com/benjamin-cogne/Sashimi-viewer/releases/latest/download/sashimi-viewer.html)
     (a versioned file you can validate and keep), or
   - open it online at **https://benjamin-cogne.github.io/Sashimi-viewer/** (GitHub Pages, same file), or
   - take [`sashimi-viewer.html`](sashimi-viewer.html) from this repository (*Download raw file*, keep the `.html` extension).
2. **Open it** in Edge, Chrome or Firefox (double-click; no server needed).
3. **Add files**: drag one or several BAM (+ `.bai`) or CRAM (+ `.crai`) files onto the page.
   The first file is the *primary sample*, the others are *comparison samples* (controls, parents, other patients).
   Double-click a sample chip (or its ✎) to rename it, for instance `proband`, `mother`, `control`; the track labels follow.
4. **Type a gene** (`NF1`), an ENSG id, or coordinates (`chr17:31,229,000-31,231,000`) and press **Open**.

Nothing is stored between sessions. Close the tab and the data is gone.

## What it answers

| Clinical question | What the viewer shows |
|---|---|
| Does this intronic or synonymous VUS create a cryptic splice site or exon skipping? | A **red dashed arc** for a junction absent from the comparison samples; click it for donor/acceptor **c. positions**, the predicted transcript in **r. notation**, the reading frame and the **NMD verdict** (55-nt rule, last-exon escape). |
| Is the aberrant junction in *cis* with the variant? | **Reads track** with mismatches (primary sample by default, any sample, or *All samples* for one reads track under each coverage track), then **Collapse** into consensus groups (local haplotype × splicing pattern): the alternate allele seen only in exon-skipping reads is explicit. |
| How much of the transcript is affected? | Per-sample **ψ** (rMATS-style inclusion) of the junction against its canonical alternative; **exon usage** from read depth compared across the open files (DEXSeq-style relative usage, robust z-score). |
| Is the gene on the minus strand? | The axis is reversed so the transcript reads 5′→3′ left to right (positions decrease to the right, unlike IGV); the transcript track then carries a **red antisense warning** so nobody misreads a coordinate. |
| Is there intron retention or a cryptic exon? | Switch from **equal introns** (exon-focused review, MISO / ggsashimi convention) to **genomic scale** and look at the coverage. |
| Is this an exon-level deletion or duplication? | Exon usage panel: median depth of each coding MANE exon relative to the other exons of the gene, per sample. |
| Is this junction normal in some tissues? | **GTEx tissue tracks** (v10, v8 fallback, hg38) as reference splicing profiles. |
| Can I open it straight from the variant page of my interpretation tool? | **Deep links** (`#variant=NC_000017.11:g.43094464G>A&pad=100&reads=1`) open the browser on the variant ±100 bp with the variant marked; the BAM/CRAM files added afterwards appear as tracks. See *Open on a variant from another tool*. |
| Is this "mismatch" a known polymorphism? | **Common SNPs** track (dbSNP 155, MAF ≥ 1 %) and **Variant sites (★)** called from the reads of the window. |
| What does the protein look like afterwards? | **Splicing cartoon** (experimental): animated pre-mRNA, spliced mRNA, translation with UniProt/Pfam domains, NMD verdict, exportable as SVG/PNG. |

Plots export as **SVG** (vector, publication-ready) for reports.

## Input requirements

- **Alignments from a spliced aligner** (STAR, HISAT2, `minimap2 -ax splice`, Dragen RNA…).
  Junctions are counted from CIGAR `N` gaps, so a DNA aligner run on RNA-seq (BWA, Bowtie2) shows
  coverage but no arcs.
- **Coordinate-sorted and indexed**: `sample.bam` + `sample.bam.bai` (or `sample.bai`),
  `sample.cram` + `sample.cram.crai`. Add the index in the same drop; files are paired by name.
- **Genome build** selected in the header (GRCh38/hg38 default, GRCh37/hg19) must match the alignment.
- **CRAM needs its reference.** Drop the indexed FASTA used at alignment time (`.fa` + `.fai`, or
  bgzipped `.fa.gz` + `.gzi`) with the CRAM files. Without it the viewer fetches the needed sequence
  from the UCSC API, which is slower and requires the network. A FASTA also makes the reads track
  (mismatches) work offline for BAM files.
- **Read filters** follow `samtools` defaults: unmapped, secondary, QC-fail and duplicate reads are
  ignored. **Unique reads** keeps `NH:1` reads (STAR/HISAT2) or MAPQ ≥ 30 when no `NH` tag is present.
- **Library strandness** (fr-firststrand / dUTP rule) is detected per file and used for the exon
  usage statistics.

## Deploying in a clinical laboratory

The viewer is a static file: deploy it the way you deploy a PDF.

| Option | How | When |
|---|---|---|
| **Shared drive / intranet** | Copy `sashimi-viewer.html` from a release to a network folder or any intranet web server. | Validated, versioned copy for the whole lab; no dependence on github.com. |
| **GitHub Pages** | *Settings → Pages*, then either *Source: Deploy from a branch* (the root page forwards to `sashimi-viewer.html`, no Actions needed) or *Source: GitHub Actions* (`.github/workflows/pages.yml` builds and publishes on every push to `main`). | Always-current URL to share with colleagues. |
| **Release asset** | `git tag v1.1.0 && git push --tags`: the build workflow attaches the viewer to a GitHub release. | Traceable versions for your quality system. |

**Network.** The reads never leave the machine. The browser only queries public annotation APIs
with a gene name or a genomic interval, so these hosts must be reachable through the hospital
proxy (proxies usually allow browsers while blocking servers):

| Data | Host | Notes |
|---|---|---|
| Transcript models (RefSeq `NM_`/`NR_`, MANE Select, RefSeq Select), gene search, neighbouring genes, protein domains, reference bases | `api.genome.ucsc.edu` | hg38 and hg19 |
| Fallback for the row above, aliases and ENSG ids | `rest.ensembl.org`, `grch37.rest.ensembl.org` | models are then Ensembl `ENST` |
| Common SNPs | `api.genome.ucsc.edu` (`dbSnp155Common`), Ensembl variation as fallback | |
| GTEx tissue profiles | `gtexportal.org` (API v2) | hg38 only |

**Browsers.** Current Microsoft Edge, Google Chrome or Firefox. Internet Explorer and very old
browsers show a plain notice instead of a blank page.

**Limits.** Coverage windows up to 5 Mb, reads track up to 250 kb, 10,000 drawn reads (40,000 when
collapsing). Deep genes over large windows take a few seconds because decoding runs in the page.
Files can be tens of gigabytes: only the indexed slices of the window are read.

**Scope.** This is a visualisation and review tool for trained users. Splice, ψ, exon usage and NMD
calls are computed from the reads and the MANE model in view; confirm findings with your
validated pipeline before reporting.

## Open on a variant from another tool (deep links)

Any tool that knows a genomic HGVS notation (a variant interpretation site such as
[MobiDetails](https://mobidetails.chu-montpellier.fr/), a report generator, a database export) can
open the viewer on that variant with a link. The parameters go in the URL fragment (`#…`), which
never reaches a server log; the query string (`?…`) is accepted too.

```
https://benjamin-cogne.github.io/Sashimi-viewer/#variant=NC_000017.11:g.43094464G>A&label=BRCA1%20c.5266dupC&pad=100&reads=1
```

The browser opens at once on the variant ±100 bp: gene model, reference bases, the variant drawn as
a labelled marker with a guide line through every track, the position pinned above the ruler, and
the reads track armed. The alignments still come from the user (a web page cannot fetch a BAM by
itself): the BAM or CRAM files added afterwards, with the button or by dropping them on the page,
appear as tracks in the same window. Pressing *Open* again reproduces the linked window.

The same works without a link: typing a gene or a locus and pressing *Open* before adding any
file shows the annotation alone (gene model, all transcripts, common SNPs, GTEx tissue tracks).

| Parameter | Meaning | Default |
|---|---|---|
| `variant` | HGVS g. with an `NC_` accession (`NC_000017.11:g.43094464G>A`; the accession version sets the build), `chr17:g.43094464G>A`, a pseudo-VCF `17-43094464-G-A` or `chr17:43094464:G:A`, or a bare `chr17:43094464`. Several separated by commas. Deletions and duplications of 50 bp or more draw as bands. | one of `variant` / `locus` |
| `locus` | Position or interval to show instead of a window derived from the variants (`chr17:43094464`, `chr17:43000000-43100000`). | |
| `pad` | Flank in bp shown on each side of the variant, locus or band. | 100 |
| `label` | Text drawn next to the marker, typically the c. or p. notation; one per variant, comma-separated. | the g. notation |
| `gene` | Symbol used only when no RefSeq gene covers the window (deep intergenic positions). | inferred from the position |
| `build` | `GRCh38` or `GRCh37` (`hg38` / `hg19` accepted). Only needed when no variant carries an accession; a notation that disagrees with it is flagged in the marker tooltip. | GRCh38 |
| `reads` | `1` opens with the reads track on, which is what a ±100 bp window is for. | off |

Example link for a MobiDetails variant page (Jinja-style template; MobiDetails holds the hg38 and
hg19 genomic HGVS with the `NC_` accession, the gene and the c./p. notations):

```html
<a target="_blank" rel="noopener"
   href="https://benjamin-cogne.github.io/Sashimi-viewer/#variant={{ hg38_g_hgvs | urlencode }}&label={{ c_hgvs | urlencode }}&pad=100&reads=1">
  Open in Sashimi viewer (RNA-seq)
</a>
```

Passing the hg38 notation is enough; a laboratory whose alignments are on GRCh37 can send the hg19
accession instead, or add `build=GRCh37`. For deep-intronic variants send a larger `pad`, or the exon
window as `locus`. The link works the same with the viewer opened from a file share: the parameters
are read by the page, not by a server.

## Privacy

BAM/CRAM/FASTA files are opened with the browser's File API and decoded with
[`@gmod/bam`](https://github.com/GMOD/bam-js), [`@gmod/cram`](https://github.com/GMOD/cram-js) and
[`@gmod/indexedfasta`](https://github.com/GMOD/indexedfasta-js). No byte of them leaves the machine.
The only outgoing requests are the annotation queries listed above; they carry a gene name or an
interval, never sample data, sample names or file names. Nothing is written to disk or to browser
storage.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| "No index found for X" | Drop the `.bai` / `.crai` in the same selection; names must match (`x.bam` + `x.bam.bai` or `x.bai`). |
| "Gene lookup failed" | Check the symbol and the build; check that `api.genome.ucsc.edu` or `rest.ensembl.org` is reachable from the browser. |
| CRAM is slow or fails | Add the reference FASTA (+ `.fai`, + `.gzi` if bgzipped) used for alignment. |
| Coverage but no arcs | The BAM comes from an unspliced aligner, or *Min reads* is above the junction's support. |
| "This file is the development entry" | You opened `index.html` from disk. Use `sashimi-viewer.html`. On a web server the page forwards to it by itself. |
| Exon usage panel is grey | Fewer than the required reads, or a single file open: usage is compared across the open files. |

## Build from source

Requires Node.js ≥ 20 (`.nvmrc` pins 22).

```bash
npm ci               # exact versions from package-lock.json
npm run build        # tsc + vite → dist/index.html, copied to dist/sashimi-viewer.html and ./sashimi-viewer.html
npm run dev          # development server at http://localhost:3000/
npm run typecheck    # tsc only
```

The build produces **one self-contained HTML file** (`vite-plugin-singlefile`: all JavaScript and
CSS inlined, no external assets) so that it can be opened from a file share without a web server.
`sashimi-viewer.html` at the repository root is tracked in git on purpose: users download it
directly, and the CI warns when it is stale relative to `src/`. Commit the rebuilt file together
with source changes.

Continuous integration:

- `.github/workflows/build.yml`: type-checks and builds on every push and pull request, uploads the
  viewer as a workflow artifact, and on a tag `v*` attaches it to a GitHub release.
- `.github/workflows/pages.yml`: builds and publishes the viewer on GitHub Pages from `main` when
  Pages is set to *Source: GitHub Actions*. With *Deploy from a branch* the repository root is served
  as is and `index.html` forwards to the tracked `sashimi-viewer.html`, so either setting works.
  Delete the file if you do not want a public URL.

<details>
<summary><strong>Repository layout</strong></summary>

```
index.html                       development entry (Vite input); shows a notice if opened directly
sashimi-viewer.html              built single-file viewer (tracked, what users download)
scripts/finish.mjs               copies the build output to dist/ and to the repository root
docs/logo/                       logo (SVG, PNG, favicon.ico, social preview) and the script that regenerates it
src/standalone/
  main.tsx                       page shell: file picker, gene box, build selector, mounts the viewer
  link.ts                        deep links: variant / locus / pad / label / build / reads from the URL
  localSource.ts                 SashimiDataSource over local files (@gmod/bam, @gmod/cram, @gmod/indexedfasta)
  alignments.ts                  coverage runs, junction counts, CIGAR/mismatch decoding, strandness
  collapse.ts                    variant-site calling and consensus-group collapsing of reads
  ucsc.ts                        UCSC Genome Browser API client (RefSeq / MANE models, sequence, domains)
  ensembl.ts                     Ensembl REST client (fallback, ENSG resolution, GRCh37)
  snps.ts                        dbSNP common variants (UCSC, Ensembl fallback)
  gtex.ts                        GTEx Portal API v2 client
src/components/
  SashimiViewer.tsx              the viewer (tracks, tooltips, panels, exon usage, reads)
  sashimi/geometry.ts            scales, arcs, HGVS cDNA positions, exon usage statistics
  sashimi/spliceModel.ts         splice event → mRNA / protein / NMD model
  sashimi/SpliceCartoon.tsx      animated splicing cartoon
  sashimi/knownVariants.ts       HGVS g. / ISCN parsing for the known-variants panel
  sashimi/datasource.ts, types.ts   interfaces shared with the parent application
```

</details>

## Citation

If the viewer contributes to a publication, please cite it
(see [`CITATION.cff`](CITATION.cff); GitHub's *Cite this repository* button formats it):

> Cogné B. *Sashimi viewer: in-browser Sashimi plots from BAM/CRAM files.* CHU Nantes, 2026.
> https://github.com/benjamin-cogne/Sashimi-viewer

## References

- Katz Y, Wang ET, Airoldi EM, Burge CB. *Analysis and design of RNA sequencing experiments for
  identifying isoform regulation.* Nat Methods 2010;7:1009–1015. doi:10.1038/nmeth.1528 (the
  Sashimi plot, MISO `sashimi_plot`).
- Garrido-Martín D, Palumbo E, Guigó R, Breschi A. *ggsashimi: Sashimi plot revised for browser- and
  annotation-independent splicing visualization.* PLoS Comput Biol 2018;14:e1006360.
- Morales J et al. *A joint NCBI and EMBL-EBI transcript set for clinical genomics and research
  (MANE).* Nature 2022;604:310–315.
- Anders S, Reyes A, Huber W. *Detecting differential usage of exons from RNA-seq data (DEXSeq).*
  Genome Res 2012;22:2008–2017.
- Shen S et al. *rMATS: robust and flexible detection of differential alternative splicing from
  replicate RNA-Seq data.* PNAS 2014;111:E5593–E5601.
- Kurosaki T, Popp MW, Maquat LE. *Quality and quantity control of gene expression by
  nonsense-mediated mRNA decay.* Nat Rev Mol Cell Biol 2019;20:406–420 (55-nt rule used for the
  NMD verdict).
- GTEx Consortium. *The GTEx Consortium atlas of genetic regulatory effects across human tissues.*
  Science 2020;369:1318–1330.
- GMOD JavaScript libraries: [bam-js](https://github.com/GMOD/bam-js),
  [cram-js](https://github.com/GMOD/cram-js), [indexedfasta-js](https://github.com/GMOD/indexedfasta-js).

## License

[CC BY-NC 4.0](LICENSE) © 2026 Benjamin Cogné, CHU Nantes. Free for non-commercial use with
attribution, which covers diagnostic and research use in public hospitals and academic
laboratories. Contact the author for commercial licensing.
