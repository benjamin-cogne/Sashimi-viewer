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

> [!WARNING]
> Artificial intelligence is changing the practice of genetics in many
> ways, including how visualization tools like this one get built. This viewer was designed by a
> molecular geneticist with some coding skills, but it was written with Anthropic's Claude models
> (Opus 5 and Fable 5.1). It's a first release so check the
> HGVS nomenclature, the predicted transcript and amino-acid changes and the NMD verdict before anything from it goes into a clinical report. Feedback and bug reports are
> welcome and appreciated.
>
> — Benjamin Cogné, Nantes University Hospital (CHU de Nantes), France
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
   A status pill next to the buttons shows what the page is doing with a drop or a folder before the
   chips exist (reading the drop, listing the folder); each chip then spins while the file's header is
   read (a few blocks from the start of the file, never the index), and its badge shows the library
   type. The index is read when the first region is requested.
4. **Type a gene** (`NF1`), an ENSG id, or coordinates (`chr17:31,229,000-31,231,000`) and press **Open**.

Nothing is stored between sessions. Close the tab and the data is gone.

## What it answers

| Clinical question | What the viewer shows |
|---|---|
| Does this intronic or synonymous VUS create a cryptic splice site or exon skipping? | A **red dashed arc** for a junction absent from the comparison samples; click it for donor/acceptor **c. positions**, the predicted transcript in **r. notation**, the reading frame and the **NMD verdict** (55-nt rule, last-exon escape). |
| Would skipping this exon keep the reading frame? | A small red **frameshift sign** above every coding exon whose coding length is not a multiple of three (skipping it alone shifts the frame); the first and last coding exons carry none, since skipping them removes the start or stop codon. Hover or click an exon for its coding length and codon phases. |
| Is the aberrant junction in *cis* with the variant? | **Reads track** with mismatches (primary sample by default, any sample, or *All samples* for one reads track under each coverage track), then **Collapse** into consensus groups (local haplotype × splicing pattern): the alternate allele seen only in exon-skipping reads is explicit. |
| How much of the transcript is affected? | Per-sample **ψ** (rMATS-style inclusion) of the junction against its canonical alternative; **exon usage** from read depth compared across the open files (DEXSeq-style relative usage, robust z-score). |
| Is the gene on the minus strand? | The axis is reversed so the transcript reads 5′→3′ left to right (positions decrease to the right, unlike IGV); the transcript track then carries a **red antisense warning** so nobody misreads a coordinate. |
| Is there intron retention or a cryptic exon? | Switch from **equal introns** (exon-focused review, MISO / ggsashimi convention) to **genomic scale** and look at the coverage. |
| Is this an exon-level deletion or duplication? | Exon usage panel: median depth of each coding MANE exon relative to the other exons of the gene, per sample. |
| Is this junction normal in some tissues? | **GTEx tissue tracks** (v10, v8 fallback, hg38) as reference splicing profiles. |
| Can I open it straight from the variant page of my interpretation tool? | **Deep links** (`#variant=NC_000017.11:g.43094464G>A&pad=100`) open the browser on the variant ±100 bp with the variant marked; the BAM/CRAM files added afterwards appear as tracks. See *Open on a variant from another tool*. |
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
- **Very deep libraries** (targeted RNA-seq, highly expressed genes) are read on a budget so the page
  never freezes: each track decodes at most 250,000 reads per request. The window is first sized from
  the index (the margins around the visible region shrink, the visible region is always read in
  full); if that region alone holds more reads, every 2nd, 4th, 8th… read is kept (a systematic
  sample in file order) and depths, junction and boundary counts are scaled back by that factor.
  Such a track shows **≈ 1 read in k** next to its name and its read counts start with ≈: they are
  estimates, exact to within a few reads for anything with hundreds of reads but coarse for rare
  junctions (a junction seen once in the sample shows k reads, one seen in none shows nothing).
  Zooming in reduces the window until the counts are exact again. Percentages (Groups view,
  *Arc labels: % usage*) are unaffected in expectation. The reads track and the exon-usage
  statistics use the same filter-then-sample order, so names and sequences are decoded only for
  the reads actually drawn.
- **Equal introns** draws every intron at the same width so exons and junctions dominate; the
  *Intron width* box that appears next to it sets that width in bp-equivalents (default: the
  model's median exon length, kept between 80 and 300; clear the box to go back to it).
- **Reference model**: the top track shows the MANE Select transcript (RefSeq Select or the
  longest CDS when there is none). Tick *All transcripts* and click any model in the list to make
  it the reference: exon numbering, junction classes, HGVS and usage percentages follow it, until
  the next gene search. The × at the end of a model's label removes it from the list (a model that
  clutters the panel, a predicted XM_ model…); the panel header then offers *undo*, which brings
  back the last one removed, and *show all*. Removed models are saved with the session.
- **Depth axis**: *relative* (the default: each sample drawn as a percentage of its own maximum in
  the current window, axis 0–100 %, the maximum shown next to the sample name, so profiles of a
  shallow and a deep library can be compared by shape), *shared* (one axis for every sample,
  heights comparable) or *per sample* (each sample scaled to its own maximum, in reads).
- **Groups (aggregate view)**: *Groups…* creates named sample groups (patients, controls, a
  tissue…); the *Samples | Groups* switch then draws one pooled track per group. Coverage and
  junction reads are summed over the group's samples, the coverage is drawn relative to its own
  maximum, and each arc is labelled with a **percentage** instead of a read count. At every intron
  of the reference model the events that use its donor or its acceptor compete, and each shows its
  share of their reads, so the labels of one intron add up to 100 %: the canonical junction C, an
  alternative 5′ or 3′ site n, a pseudo-exon (an alternative-3′ arc A into a cryptic exon of at
  most 500 bp and an alternative-5′ arc B out of it, drawn in purple) weighted (A + B) / 2 with the
  same value on both arcs, exon skipping S, and intron retention weighted (R5 + R3) / 2 from the
  reads that run unspliced through the two boundaries of the intron (one aligned block with at
  least 6 bases on the exon side and 10 on the intron side), shown as a teal *IR* pill on the intron
  baseline. An exon-skipping arc spans two introns and shows 2·S over the two totals, which is the
  rMATS value 2·S / (I₁ + I₂ + 2·S) when nothing else competes at those introns; its tooltip gives
  its share at each intron. Every tooltip also gives the event against the canonical junction
  alone, rMATS-style (n / (n + C), (A + B) / (A + B + 2·C), 2·S / (I₁ + I₂ + 2·S),
  (R5 + R3) / (R5 + R3 + 2·C)). The *Intron retention* box next to *Min %* removes retention from
  the shares. Junctions touching no annotated splice site show their pooled read count (`n=…`) and
  compete nowhere. Red arcs are events seen in the first group only. The same percentages are
  available per sample in the Samples view through *Arc labels: % usage* (computed on each
  sample's own reads); the Groups view always shows percentages. In that mode *Min reads* becomes
  *Min %*, hiding events below that usage; hidden and off-screen junctions still count in the
  denominators, so a canonical arc below 100 % lists in its tooltip the other events at its intron,
  shown or not.
  With two groups or more, each percentage is followed by its **difference with every other
  group**, in points and in that group's colour (`+10 %`, `−10 %`); the tooltip names the group.
  Group colours are chosen in the *Groups…* dialog (click the swatch; *default* returns to the
  palette) and saved with the session.
- **Hiding an arc**: hover an arc and click the × at the end of its pill (or *Hide arc* in the
  junction's detail panel) to remove it from every track, for instance the `n=…` junctions of
  another transcript that clutter a group view. Hidden arcs still count in the percentages, like
  arcs under the thresholds; the *hidden arcs · show* chip in the toolbar brings them back, and the
  list is saved with the session. The same panel has **A− / A+** buttons that shrink or enlarge that
  junction's label on every track (70–250 %), for the arcs a figure should emphasise; also saved
  with the session.

## RNA-seq and genomic DNA

Each sample is either **RNA-seq** or **genomic DNA** (exome, genome, long reads), and the viewer
decides which from the file itself: the aligner named in the header (`@PG` lines: STAR, HISAT2,
TopHat, `minimap2 -ax splice`, DRAGEN RNA mean RNA; BWA, bowtie2, minimap2 genome presets, Isaac,
pbmm2, DRAGEN mean DNA), then the reads of the first gene opened, which have the last word: at a
multi-exon gene with at least 200 reads, 2 % or more of them spliced (a CIGAR `N` gap) means
RNA-seq, fewer than 0.2 % means DNA. The sample chip carries the verdict as an **RNA**, **DNA** or
**?** badge whose tooltip gives the evidence; a click on the badge switches the type and the choice
is then kept (saved with sessions and carried by exports). An undetermined sample behaves as RNA.

A DNA track draws the coverage and, when asked, the reads, and nothing about splicing: no junction
arcs, pills, usage percentages or retention. Its label says *DNA*, and the track is shorter since
no arc space is needed. Options that only concern splicing (*Reads | Usage*, *Min reads*, *Min %*,
*Intron retention*, the splicing legend) stay while at least one RNA track is shown, and are put
away when every shown sample is DNA. RNA and DNA samples can share a page and a view (a proband's
RNA next to the parents' genomes), but not a group.

**Structural hints.** Where an RNA track shows junction arcs, a DNA track shows the structural
evidence of its reads, drawn as evidence and never as calls: **deletions inside reads** (a CIGAR `D`
run of 50 bp or more) as solid red arcs; **split reads**, read as chains: every part of a read (the
primary alignment and its supplementary alignments, from the SA tag, whichever of them fall in
the window) is placed along the read by its clips, the parts are ordered along the read and each
pair of adjacent parts is one breakpoint, each read counted once. The breakpoint is typed by where
the read continues: further along the same strand, a dashed purple *deletion-type* arc; backwards,
a dashed green *duplication-type* arc; on the other strand, a dashed blue *inversion* arc; on
another chromosome, a purple `→ chr` pill; and an unaligned stretch of the read between two parts
adjacent on the reference, a purple `ins` pill with its size. Breakpoints are rounded to 5 bp.
**Discordant pairs** (insert size above five times the window's median, at least 1 kb, or both
mates on the same strand) are dashed amber arcs between 500 bp bins, and **soft-clip clusters**
(at least 3 reads clipped by 20 bases or more at one position, split reads excluded) teal pills on
the baseline (⇤ clipped before, ⇥ clipped after). *Min supporting reads* sets the support needed to draw a
hint. An arc's pill shows the read count (≈ on sampled windows); clicking an arc opens a panel with
the count in every DNA sample shown, the median insert size, and a g. notation for a deletion.
Arcs can be hidden, resized and dragged like junction arcs, and DNA groups pool the evidence of
their members. Long reads carry deletions, split reads and clips; short-read pairs add the
discordant pairs (CRAM mate fields are read when the file records them).

**Consensus reads.** *Consensus* (on by default) draws mismatches and indels only where a variant
site is called (at least 3 reads and *Min VAF*), so sequencing errors do not paint every read. It
is offered for every genomic DNA reads track, short reads included, and for long reads whatever
the library; off, every mismatch and indel of every read is drawn.

**Long reads (ONT, PacBio).** A reads track whose median aligned length is above 1 kb gets two
more noise controls: *Min VAF (long)* (20 %) is the allele fraction a site needs on long reads,
above the short-read *Min VAF*; *Min indel* (10 bp) hides and leaves uncalled the shorter indels
typical of homopolymer errors. The controls appear next to *Collapse* when such reads are shown
and are saved with the session.
Deletions of 50 bp or more inside reads remain structural evidence whatever the setting. On DNA
tracks the reads carry no exon–intron boundary outline (that teal mark is an intron-retention
device for RNA).

**Variants and allele balance.** A DNA track opens as a plain coverage histogram: no read is
decoded until asked. Its variant sites are drawn as allele-fraction bars on the coverage (no star
strip on DNA tracks; the bar's label gives the fraction, its tooltip the site), from the reads track
when it is shown, or from the **variants** chip next to the sample name otherwise: the chip scans
every read of the current window, whatever its width (tile by tile, so memory stays bounded, with
the progress shown on the chip and a click to stop), calls every site above *Min VAF*, and then reads
*variants ✓ N*. From then on the variants follow the window: moving or widening it scans only the
part not scanned yet and merges it, so the bars are always those of the frame shown; another
chromosome, the unique-reads switch, a lower *Min VAF* or a changed long-read threshold start the
window over (a raised *Min VAF* filters at once). A click on the ✓ chip forgets the variants (plain
coverage until the next click). The **Variants** toggle of the options panel removes the bars (and the
chip) for a plain coverage. When every shown
sample is DNA the axis keeps the genomic orientation, coordinates increasing to the right, even
for a minus-strand gene. *Common SNPs* stay off by default on DNA tracks as on RNA ones; switch
them on to separate a known polymorphism from a novel change. The track label then
summarises the **allele balance** of the common SNPs covered: how many are heterozygous
(0.2 ≤ VAF ≤ 0.8) and the range of their fractions around 0.5. Heterozygous fractions far from 0.5
(median deviation above 0.15 over at least 5 SNPs) are flagged *allele imbalance?* (mosaic copy
change, loss of heterozygosity, contamination), and a window with at least 8 homozygous common SNPs
and none heterozygous is flagged *no heterozygous SNP (LOH / UPD?)*. Fractions come from the drawn
reads, up to 2,500 in the window, so they are approximate on very deep data.

## Views

Every gene or locus opened from the search box in the header becomes a **view**. Opening a second
one folds the first into a tab under the sample chips (a *Views* row appears from the second view
on); opening a third keeps the two. Click a tab to reopen that view exactly as you left it: the
gene, the window, the highlight, and every option of the viewer (depth axis, filters, groups,
reads track, chosen transcript, hidden arcs, label sizes), since each view keeps its own full
snapshot. × forgets a view. A new view starts with the options of the view you came from. Views
are saved with the session and are what the HTML export carries. **SVG · N views** (header) saves
every view on one SVG page, stacked vertically under their titles (gene, locus and window), each
drawn with its own options: each tab is shown in turn while its plot is captured, then the current
view comes back. The *SVG* button inside the plot still saves the current view alone.

## Sharing a view without the alignments

**Export HTML** (next to the session buttons) downloads a copy of the viewer itself with the data of
every registered view embedded: for each view, the gene models (reference transcript, every
transcript, neighbouring genes) and, for every loaded sample, the coverage, junctions and
intron-retention counts of a window around it chosen in the export dialog (see below), together
with the sample names, the groups and every option. The recipient opens the file in
any browser: no BAM, no server, nothing to install. They see an *Exported viewer* banner, the same
sample chips and view tabs, and can switch views, zoom, pan inside each exported window, toggle
groups or usage percentages, click arcs and exons for the HGVS and frame details. Coverage is
exact where you had it exact and marked ≈ where the source had sampled it.

**Reads.** For every view whose reads track is on, the export also embeds the reads of every loaded
sample, with the reference bases and the mismatches, so the recipient sees the same pile-up, can
switch between reads and the collapsed haplotype view and change *Min VAF*. Read names are
replaced by numbers. The dialog that opens on *Export HTML* chooses the window exported around each view, for the
coverage, junctions and retention counts as much as for the reads (the view as shown; the view with
a half-width margin on each side, the default; or the widest window the viewer itself loads, up to
2 Mb for coverage and 100 kb for reads) and the number of reads per sample (up to 2,500 as displayed, about 300 kB per sample and view; up to 20,000 for
zooming in; or every read of the window, which can reach tens of MB for a deep window). Views wider
than 100 kb have no reads track and are skipped.

What the export does not carry: the exon-depth statistics (they need the alignments; the recipient
can add the BAM/CRAM files and the page then reads them as usual), reads of views whose reads track
was off, and regions outside the exported windows (a track there says *not in this exported file*). Gene
lookups, common SNPs, GTEx and reference bases still come from the network when it is available.
The exported page is a normal viewer: it can save sessions and export again. Exporting needs the
built `sashimi-viewer.html` (or the published page), not the development entry. Sizes: about 1 MB
for the viewer plus the run-length coverage of each window, typically a few hundred kB per sample
and view for a gene, more for very deep or very wide windows.

## Sessions

*Save session* in the header downloads a JSON file (the name is editable, default
`sashimi-session-GENE-DATE.json`) recording the genome build, the run folder and every alignment
by its path inside that folder (with index name and size), the sample names and order (first =
primary), the FASTA, the gene and the window shown, and every option of the viewer: depth axis,
arc labels, thresholds, reads track, groups (by sample name), the chosen reference transcript…
The variants of interest of the header are saved too.

Give the files through **+ Run folder…** (or drop the folder on the page): the page lists the
BAM/CRAM files with their index and the FASTA found inside, without reading them, and remembers
their relative paths. **+ Files…** adds individual files as before. Some browsers word their folder
dialog as an *upload* ("this will upload all files from…"): that is the browser's generic
wording for letting a page read a folder. Nothing is sent anywhere; the viewer has no server and
reads the files on the computer, as for individually added files.

*Load session* reads the JSON, switches the build and gets the files back:

- **Chrome, Edge, Opera**: the browser keeps a bookmark of the run folder (and of individually
  added files) in its own storage, not a copy of the data. Loading a session reopens the files
  from it after a single permission click (remembered by the browser until revoked). Then the
  samples, gene, window and options come back at once.
- **Firefox, Safari** (no bookmarks): the page asks for the run folder; drop it on the page or
  pick it with *Choose the folder*, and the same happens.

A browser never reveals or stores file paths as text, which is why the folder is asked for once
rather than typed. *Open with the files present* starts with a subset. Options also carry over
from one gene search to the next within a page.

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

## Variants of interest (known variants)

The **Known variants** row of the header takes any number of variants to keep in sight: a locus
(`chr17:43,094,464`, or an interval `chr17:43,094,464-43,094,470`, an HGVS genomic notation such
as `chr17:g.43094464A>G` or `NC_000017.11:g.43094464A>G`, or a VCF-like line) and a free label
(`BRCA1 p.Glu23Asp`, a sample name, anything), then **+ Add**. Each becomes a chip (× removes it)
and is drawn on every view like the variants of a deep link: the **Known variants** option of the
viewer, on by default as soon as one exists, shows a panel under the transcript with the labels,
guide lines through every track, and a *go to…* menu to centre the view on one of them. A variant
on another chromosome can be chosen too: the gene at its position is opened (coding gene first, then
the largest overlap) and the view centred on the variant; when no RefSeq gene is there, the window
opens alone, in genomic (sense) orientation. They are saved with the session and carried by the
exported page.

## Open on a variant from another tool (deep links)

Any tool that knows a genomic HGVS notation (a variant interpretation site such as
[MobiDetails](https://mobidetails.chu-montpellier.fr/), a report generator, a database export) can
open the viewer on that variant with a link. The parameters go in the URL fragment (`#…`), which
never reaches a server log; the query string (`?…`) is accepted too.

```
https://benjamin-cogne.github.io/Sashimi-viewer/#variant=NC_000017.11:g.43094464G>A&label=BRCA1%20c.5266dupC&pad=100
```

The browser opens at once on the variant ±100 bp: gene model, reference bases, the variant drawn as
a labelled marker with a guide line through every track, the position pinned above the ruler, and
the reads track armed. The alignments still come from the user (a web page cannot fetch a BAM by
itself): the BAM or CRAM files added afterwards, with the button or by dropping them on the page,
appear as tracks in the same window. Pressing *Open* again reproduces the linked window.

The same works without a link: typing a gene or a locus and pressing *Open* before adding any
file shows the annotation alone (gene model, all transcripts, common SNPs, GTEx tissue tracks).
A locus searched by coordinates (in the header or the viewer's search box) is highlighted as a
dashed band or line with its position pinned above the ruler; right-click the band, the line or
the pinned label to remove the highlight (the next search draws a new one).

| Parameter | Meaning | Default |
|---|---|---|
| `variant` | HGVS g. with an `NC_` accession (`NC_000017.11:g.43094464G>A`; the accession version sets the build), `chr17:g.43094464G>A`, a pseudo-VCF `17-43094464-G-A` or `chr17:43094464:G:A`, or a bare `chr17:43094464`. Several separated by commas. Deletions and duplications of 50 bp or more draw as bands. | one of `variant` / `locus` |
| `locus` | Position or interval to show instead of a window derived from the variants (`chr17:43094464`, `chr17:43000000-43100000`). | |
| `pad` | Flank in bp shown on each side of the variant, locus or band. | 100 |
| `label` | Text drawn next to the marker, typically the c. or p. notation; one per variant, comma-separated. | the g. notation |
| `gene` | Symbol used only when no RefSeq gene covers the window (deep intergenic positions). | inferred from the position |
| `build` | `GRCh38` or `GRCh37` (`hg38` / `hg19` accepted). Only needed when no variant carries an accession; a notation that disagrees with it is flagged in the marker tooltip. | GRCh38 |
| `reads` | `1` opens with the reads track already on (add it when the link targets a ±100 bp window and the alignments are wanted at once). By default the reads track stays off until *Reads* is ticked in the toolbar or the *reads* chip of a sample is clicked. | off |

Example link for a MobiDetails variant page (Jinja-style template; MobiDetails holds the hg38 and
hg19 genomic HGVS with the `NC_` accession, the gene and the c./p. notations):

```html
<a target="_blank" rel="noopener"
   href="https://benjamin-cogne.github.io/Sashimi-viewer/#variant={{ hg38_g_hgvs | urlencode }}&label={{ c_hgvs | urlencode }}&pad=100">
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
| Counts start with ≈, "≈ 1 read in k" next to a sample | The window holds more than 250,000 reads for that sample: it was sampled 1 in k and scaled back. Zoom in for exact counts. |

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
- `.github/workflows/traffic.yml`: every Monday (or on demand), archives the repository traffic that
  GitHub keeps for only 14 days (views, clones, referrers, popular pages) and the adoption figures
  (stars, forks, release downloads) into `docs/stats/*.csv`, so usage can be reported over time.
- `.github/workflows/pages.yml`: builds and publishes the viewer on GitHub Pages from `main` when
  Pages is set to *Source: GitHub Actions*. It also publishes the `dev` branch, when it exists, as
  an unlisted preview under `/dev/` with a red **DEV MODE** banner, the branch and commit in the
  page, a `[DEV]` tab title and a red icon (`VITE_DEV_MODE=1` at build time). A push to either
  branch republishes both. Work on `dev`, check the preview, merge into `main`; commit the rebuilt
  `sashimi-viewer.html` on `main` only, so merges never conflict on it. With *Deploy from a branch* the repository root is served
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
