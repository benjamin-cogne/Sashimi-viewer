<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logo/logo-lockup-dark.svg">
  <img src="logo/logo-lockup.svg" alt="Sashimi viewer" height="64">
</picture>

# User manual · version 1.2

Sashimi plots and read-level evidence from your own BAM/CRAM files, in one HTML file that runs entirely in the
browser: no server, no installation, no upload. The [README](../README.md) is the short introduction.

> [!WARNING]
> A review tool for trained users, not a medical device. Check the HGVS nomenclature, the predicted transcript and
> the NMD verdict before anything from it goes into a clinical report.

This manual says what each feature shows and how to use it. **Every control and mark has a tooltip with the full
detail** (formulas, thresholds, what was counted): hover it when a number needs explaining. The methods, with their
validation and references, are in [Methods and validation](algorithms.md).

[Start](#start) · [The screen](#the-screen) · [RNA-seq](#rna-seq) · [Genomic DNA](#genomic-dna) ·
[Reads and phasing](#reads-haplotypes-and-phasing) · [Methylation](#methylation) ·
[Deep links](#variants-of-interest-and-deep-links) · [Sessions and export](#sessions-and-html-export) ·
[Builds](#reference-builds-and-lift-over) · [Privacy](#privacy) · [Deploying](#deploying-in-a-clinical-laboratory) ·
[Troubleshooting](#troubleshooting) · [Build](#build-from-source)

## Start

1. **Get the viewer**: open **https://benjamin-cogne.github.io/Sashimi-viewer/**, or download
   [`sashimi-viewer.html`](https://github.com/benjamin-cogne/Sashimi-viewer/releases/latest/download/sashimi-viewer.html)
   from the latest release (a versioned file to validate and keep) and double-click it. Edge, Chrome or Firefox.
2. **Add files**: drop BAM (+ `.bai`) or CRAM (+ `.crai`) files or a run folder on the page, or use **+ Add** ›
   *Run folder…* / *Files…*. The first file is the **primary sample** (★), the others comparison samples (controls,
   parents, other patients). Double-click a chip to rename it.
3. **Open a gene**: type a symbol (`NF1`), an ENSG id or coordinates (`chr17:31,229,000-31,231,000`) and press
   **Open**. Names are suggested from the third letter (HGNC symbols, aliases, previous symbols).

**Files.** Coordinate-sorted and indexed (the index is matched by name). RNA-seq must come from a spliced aligner
(STAR, HISAT2, `minimap2 -ax splice`, DRAGEN RNA): junctions are CIGAR `N` gaps. **CRAM needs its FASTA** (`.fa` +
`.fai`, or `.fa.gz` + `.gzi`) in the same drop; without one, bases come from the UCSC/Ensembl APIs. Unmapped,
secondary, QC-fail and duplicate reads are left out, as by `samtools`; see [Filters](#filters).

**Chips.** The badge says **RNA**, **DNA** or **?** (treated as RNA), read from the aligner in the header (`@PG`,
including DRAGEN and Oxford Nanopore basecalling models), else from the reads of the first multi-exon gene (≥ 2 %
spliced: RNA; < 0.2 %: DNA). RNA chips add the strandedness, **REV** (dUTP), **FWD** (QuantSeq FWD, ONT direct RNA),
**UNS**, or **±?** until a gene with ≥ 50 spliced reads is seen. Click a badge to set it; hover for the evidence.
Files of one name in different folders take the folders that tell them apart (`p1/Aligned.sortedByCoord.out`).

**Moving.** Drag to pan · Ctrl+drag zooms into a region · Ctrl+scroll zooms at the cursor · double-click resets ·
Shift+drag inspects a repeat · right-click a highlight removes it. The search box **inside the plot** (next to the
gene name) also takes positions on the transcript drawn: `c.234`, `c.-12`, `c.*30`, `c.288+1`, `c.234_267`, whole
descriptions (`c.2033dup`: the change is ignored), `n.412`, `12`, `exon 12`, `exons 3-5`. A position the transcript
does not have is refused with the reason (often a sign it was written against another transcript).

## The screen

**App bar**: logo · sample chips · **+ Add** · search box with the *GRCh 38 | 37* build switch · **Save session** /
**Load session** · **Export** (*Export HTML*, *SVG · N views*) · ⌃ folds it to one line. Second line: the **Views**
tabs, the **Known variants** menu, the reporting reminder.

**Views.** Each gene or locus opened is a tab with its own window, tracks and options. Searching one that has a view
switches to it; **Shift+Enter** (or Shift+click on *Open*) opens a second view on purpose. × forgets a view.

**Samples and groups.** **n / N shown** shows or hides each loaded sample (*Show all*, *Hide all*, *Primary only*).
A track's *reads* chip shows its reads; *★ make primary* promotes it. **⋯** › *Sample groups…* defines groups; the
*Samples | Groups* switch draws one pooled track per group. RNA and DNA can share a view, not a group.

| Toolbar section | Controls |
|---|---|
| **Arcs** (RNA) | *Reads \| Usage* labels · *min* reads (3) or *min* % (1 %) · *Retention* · *n hidden · show* |
| **Structure** (all DNA) | *min* reads for a structural hint: 3, or 20 on tracks deeper than 200× (*auto*) |
| **Show** | *Reads* (RNA) · *Strands* · *Variants* + *min VAF* (RNA) · *All transcripts* · *SNPs* + *AF ≥* · *Equal introns* · *Known variants* + *go to…* |
| **Layers** (DNA) | **C** Coverage · **V** Variants · **M** Methylation · **R** Reads, for all DNA samples; *CpG islands only*; *min VAF* |
| **Depth** | *shared* (heights comparable) · *own* · *% max* (default: shapes comparable) |
| **Filters** | *Minimum MAPQ*, *Secondary alignments* (a badge counts those on) |

Layers are drawn in the order C, V, M, R under each DNA track; a layer switched off is no longer read.

### Filters

- **Minimum MAPQ** (*all*, 1, 10, 20, 30; `samtools view -q`) applies to coverage, arcs, variant sites and scans,
  and reads. STAR scores 255 unique and 3/1/0 for 2, 3–4, ≥ 5 loci; HISAT2 60/1/0; bwa and minimap2 0–60. So **≥ 10
  keeps STAR's and HISAT2's unique reads**, ≥ 20–30 bwa's and minimap2's confident placements. Sessions saved with
  the former *Unique reads* switch open at ≥ 30.
- **Secondary alignments** also counts the other placements of multi-mapped reads (flag `0x100`), as IGV does, e.g.
  a variant of a paralogue then shows at both copies. Methylation stays on primary alignments. Aligners give
  secondaries MAPQ 0–3, so a threshold ≥ 1 removes most of them again.
- In the reads track, alignments with **MAPQ < 10 are drawn hollow**.

**Reads strip** (with the reads track): *reads of* (a sample or *all*) · *window ≤* 100 kb (default) to 2 Mb ·
*Raw | Collapsed* · *group*, or *haplotypes* and *phase* · *colour* (*strand*, *1st-of-pair*: IGV's pink and
lavender) · *Pairs* · *Clipped* · *Inserted* · *Consensus* · *min VAF* (10 %) · *long reads* (20 %, tracks whose
median aligned length exceeds 1 kb).

**Panels** under the ruler (gene model, sequence strip, known variants, SNPs, transcripts, tracks) move by their
grip ⋮⋮; double-click a grip to restore the order. The **sequence strip** under the gene model shows the bases and,
for every coding gene in view (the gene and its neighbours, one row each, named in the margin), its amino acids in
its own frame (start green, stop red), once for all tracks.

## RNA-seq

The model is **MANE Select** (else RefSeq Select or the longest CDS). With *All transcripts* on, click a model to
make it the reference for exon numbers, HGVS and usage. **A minus-strand gene is drawn 5′→3′**: positions decrease
to the right (unlike IGV), with a red *antisense gene* warning. A page of DNA samples keeps genomic orientation.

| Mark | Meaning |
|---|---|
| Arc with a count | Spliced reads of a junction |
| **Red dashed arc** | Junction of the primary sample absent from all comparison samples |
| Purple arcs | Pseudo-exon (cryptic exon ≤ 500 bp) |
| Teal **IR** pill | Intron retention (*Usage*) |
| Dashed fuchsia arc | Long-range junction, one end outside the gene (read-through, fusion after a deletion); shown from its reads whatever *min* % |
| Solid fuchsia arc · `→ chr` | Fusion junction from split and clipped reads (inversion, duplication-type, deletion-type between genes), with its reading frame · partner on another chromosome |
| Violet, under the baseline | Opposite strand (*Strands*) |
| Red frameshift sign on an exon | Coding length not a multiple of 3: skipping it alone shifts the frame |

**Click an arc**: donor and acceptor in **c.** notation, the transcript in **r.** notation (`r.X_Ydel` for a skip,
`r.X_X+1ins…` for a pseudo-exon), the reading frame, and each sample's share against the canonical alternative.
Buttons: **Supporting reads** (the reads track on that arc's reads and mates, up to 300), **Hide arc** (still counted
in percentages), **A− / A+** (label size) and **🎬 Cartoon** (experimental): pre-mRNA, mRNA and protein with
UniProt/Pfam domains and the **NMD verdict** (a stop > 55 nt before the last exon junction is degraded; in the last
exon, within 55 nt of the last junction, or within 150 nt of the start codon it escapes), as SVG/PNG.

**Click a fusion arc** (*GENEA exon 2 → GENEB exon 2*), or a dashed long-range arc whose two ends lie on exons of two
genes on one strand (a read-through such as *GLE1 exon 12 → SPTAN1 exon 2*): the supporting reads per sample, the fusion
transcript in HGVS notation (`NM_…:r.-87_1776::NM_…:r.-3_*328`) and the **reading frame**. From the two partners' coding
sequences it reads **in frame** (the 3′ partner's protein follows) or **frameshift**, or says the junction lies in a
UTR. The fused sequence is then translated across the junction from the reference: the peptide on each side
(`…QKLR[G]|VSCQ…`, the codon across the junction in brackets) and where the protein ends. That shows whether an
in-frame fusion carries the 3′ partner's protein to its own stop, or stops early (in its 5′ UTR, at the junction
codon). Both partners' models must be loaded: pan over the partner gene if the popover asks.

**🎬 Cartoon** on a fusion (experimental): the 5′ partner's transcript on top, the fused one below (its exons, then the
3′ partner's in teal), spliced, translated across the junction and judged like any transcript. The verdict is **in-frame
fusion protein** when translation reaches the 3′ partner's own stop codon; otherwise the stop is premature and the NMD
rules apply to the fused transcript. Both partners' UniProt/Pfam domains are drawn on the fusion protein: the 5′ one's
kept up to the junction, the 3′ one's carried when read in their frame. A junction in the 5′ partner's 5′ UTR uses the
3′ partner's start codon (promoter swap). The fusion protein's own stability is not predicted.

**% usage.** In *Usage*, an arc's label is its share of the reads competing at its intron, so one intron's labels sum
to 100 %: canonical, alternative 5′/3′ site, pseudo-exon (A + B)/2, skipping (2·S over the two introns), retention
(R5 + R3)/2 with *Retention* on. Hidden and off-screen arcs still count; junctions touching no annotated site show
`n=…`. The Groups view always uses percentages, and with ≥ 2 groups adds each difference in points.

**Click an IR pill**: the retained intron's HGVS (`NM_…:r.290_291ins290+1_291-1`), its effect on the protein (frame,
first stop, NMD likely unless last intron), and per sample the reads unspliced through each boundary, the canonical
junction, the rMATS-style ratio and an IRFinder-style IR ratio. Retention and cryptic exons are easiest to see at
genomic scale; *Equal introns* (off by default) compresses introns so exons and junctions dominate.

**Click an exon**: **ψ** (inclusion / (inclusion + skipping) junction reads) per sample, and the **depth-based
usage** (median exon depth over the gene's other coding exons) against the other samples, with a robust z-score from
5 controls. It points to exon-level deletions or duplications. It needs comparison samples.

**Strands** splits each stranded sample: the gene's strand above the baseline (arcs, usage, sites are its own), the
opposite strand in violet below. An opposite-strand arc's panel names the gene it joins, or says *likely the gene's
reads on the wrong strand* under 5 % of an intron with ≥ 20 reads (stranded kits leave about 0.5–3 % there).

**GTEx** tissue tracks (v10, v8 fallback; hg38 only): the tissue's median junction reads as arcs and the usage (0–1)
of each exon its junctions use, other transcripts' exons included.

**Variant sites** (*Variants*, on by default): a background scan of every read (views ≤ 3 Mb); ≥ 3 alternate reads
and *min VAF*. A star and an allele bar mark each site; click for reads by strand, the **allelic balance** (binomial
test against 50 %: *balanced*, *tilted*, near-monoallelic; in RNA, NMD, silencing or imprinting tilt it), median
qualities, the quality cells, whether it is a common SNP, **Copy** chips (`chr1:5291`, `1-5291-C-A`,
`chr1:g.5291C>A`, `NM_…:c.281C>A`; indels not 3′-shifted) and [Phase with](#phase-with).

**Long reads**: a junction within 6 bp of one with ≥ 20× its reads is counted with it; a real cryptic site holding
≥ 5 % as many reads keeps its own arc. Short-read junctions are never merged.

## Genomic DNA

**Coverage** is exact at any depth (nothing sampled); the view is counted first (*reading the view… N %*), then the
margins. Zoomed out, the line is each pixel's mean depth, a pale band its highest base, and a tick marks any base
below half the mean, so a narrow dropout stays visible on a whole gene.

**Variants track** (layer **V**): a background scan of every read (views ≤ 3 Mb). A site needs ≥ 3 alternate reads
and *min VAF* (10 %, *long reads* 20 %); every alternate base counts whatever its quality. Each site is a **bar** of
its allele fraction (50 % guide; base colour, purple insertion, black deletion), over four **quality cells**,
green pass · amber check · red likely artefact · grey too few reads:

| Cell | Measures | Amber | Red |
|---|---|---|---|
| **BQ** (SNV) | alternate bases with base quality < 20 | ≥ 25 % | ≥ 50 % |
| **HP** (indel) | reference homopolymer length | ≥ 6 (long reads ≥ 4) | ≥ 10 (≥ 7) |
| **MQ** | alternate reads with MAPQ < 20, vs other reads | ≥ 20 % and +15 points | ≥ 50 % and +30 points |
| **SB** | strand bias, alternate vs other reads (Fisher) | p < 0.001 | p < 1e-4 |
| **END** | alternate calls within 10 bases of a read end (Fisher) | p < 0.001 | p < 1e-4 |

They help judge a call (what GATK's hard filters and an IGV reviewer look at); they are not a caller's filter.
Click a site for the reads' distributions (strand, MAPQ, base quality, distance to end), Copy chips and *Phase with*.

**Haplotype lanes.** Haplotagged files (`HP`/`PS` from WhatsHap, LongPhase, HiPhase, DRAGEN) get lanes **H1** and
**H2**. A site goes on H1 when ≥ 70 % of HP 1 reads carry it and ≤ 20 % of HP 2 reads do (H2 likewise), on both when
both carry it at ≥ 70 %. **Dashed red circle**: *to check* (mosaic, artefact, mis-phasing, collapsed duplication);
**?**: on neither haplotype of its set; **grey dot**: not phased (< 3 reads per haplotype, or a red cell). Each phase
set is a band in its own hue; a red **⫽** means no read linked two sets. **HP 1 / HP 2 are the tool's labels, not a
parental origin**, arbitrary within each set, and need not match the `0|1` order of the pipeline's VCF.

**Common SNP balance.** With *SNPs* and **V** on, the track label gives the heterozygous common SNPs (0.2–0.8) and
flags *allele imbalance?* (median deviation from 0.5 > 0.15 over ≥ 5 SNPs) or *no heterozygous SNP (LOH / UPD?)*
(≥ 8 homozygous, none heterozygous). A lead to check, not a copy-number call.

**Structural evidence**, always shown on DNA tracks, is evidence and never a call:
deletions inside reads (CIGAR `D` ≥ 50 bp, solid red); split reads (dashed green duplication, dashed blue
inversion, `→ chr` and `ins` pills); discordant pairs (coloured as in the reads track; the breakpoint lies within an
insert size of the arc's ends); soft-clip clusters (≥ 3 reads clipped ≥ 20 bases: teal ⇤ ⇥ pills). Click an arc for
its count in every DNA sample and a g. notation; a duplication's tooltip says whether its whole exons keep the frame.
Hover a track's *DNA* tag to see why an arc is missing.

**Repeat inspector.** The 77 pathogenic repeat loci of STRchive get an orange ruler pill (*FMR1 CGG repeat ·
inspect ▸*); any other repeat opens with Shift+drag, then *⟲ Inspect repeat*. Long reads are sized between two
locus-specific anchors, drawn as a size **histogram** and a **waterfall** of units (pathogenic motif blue, benign
teal, interruption orange, e.g. FMR1's AGG), with up to two **alleles** (mode, P5–P95, category; *broad*: mosaicism
or stutter, not decided). Categories: FMR1 by ACMG (≤ 44, 45–54, 55–200, > 200), HTT by ACMG/ASHG (≤ 26, 27–35,
36–39, ≥ 40), others by STRchive. Under 20 spanning reads an amber banner warns and reads stopping inside count as
**lower bounds**: long alleles are spanned less often. Evidence, not a diagnostic call.

## Reads, haplotypes and phasing

**Reads track**: mismatches, indels and splice gaps, IGV-style, up to the *window ≤* width; up to 2,500 reads per
track, sampled evenly (*downsampled, zoom in for all*). **Consensus** (on) draws mismatches and indels only at called
sites; indels ≥ 50 bp are always drawn. Lower *min VAF* or turn *Consensus* off for changes in few reads (RNA editing,
subclonal variants).
- **Pairs** share a row. Discordant pairs, on DNA: green ← → (duplication type), red → ← far apart (deletion type),
  blue same strand (inversion type); amber: mate on another chromosome or pair not proper (the only class on RNA).
- **Clipped** draws soft-clipped bases, dimmed where they match the reference, and joins the parts of a split read on
  one row; **Inserted** writes inserted bases. **Click a read** for its clipped and inserted sequences, its other
  parts, and *Fetch from the primary record* for a supplementary one. **Click a soft-clip pill** for a consensus to
  paste into BLAT.

**Collapsed** folds the reads (up to 40,000 decoded, in the background):
- *haplotypes: 2 phased*: one consensus row per haplotype and phase set; coloured bases are variants carried by ≥ half
  its reads. **H1 of one phase set is unrelated to H1 of the next.** Sites the two do not split (80/20) are listed on a
  *not split* row (mis-phasing, mosaic, third haplotype, collapsed duplication).
- *phase: HP tags* uses the file's haplotags (untagged reads counted, left out); *phase: reads* phases the window in
  the page (heterozygous sites at 25–75 %, linked by ≥ 2 fragments, ≤ 20 % disagreeing; a read and its mate are one
  fragment).
- **RNA**: each haplotype row carries its junctions. A junction whose use differs between them (Fisher p < 0.001 and
  ≥ 20 points) is drawn in each haplotype's colour and counted in the header: a splice change **in cis** with that
  haplotype's variant.
- *haplotypes: any*: consensus groups, one per local haplotype × splice pattern; small groups fold into a minor
  bucket; ambiguous rows read `H1|H2 ?`.

**Group** (raw reads): *phased here* (the page's phasing, up to 10,000 reads) or *HP tags*. Each read takes its
**phase set's hue** (HP 1 darker, HP 2 lighter), the same as the variant lanes and methylation.

### Phase with

Cis or trans for two variants (a possible compound heterozygote), or which copy uses a junction (allele-specific
splicing). Click a site, then in **Phase with** type positions separated by `;` (`chr:pos`, `pos`, `chr:pos:A>G`,
`chr:pos:+2`, `1-55051215-G-GA`) and press **Phase**, or **⌖ Pick** and click sites (Esc stops). On RNA, pick a
junction arc or an *IR* label, or type the intron `chr:first-last` (1-based, as STAR's `SJ.out.tab`).
- **Evidence**, in order: the reads holding both sites, joined by name however far apart (long read, spliced read,
  mates, split read); else the haplotags, if both sites share a phase set; else, without haplotags, the reads of the
  stretch between the sites, phased. The verdict names its evidence; haplotags contradicting the reads give *conflicting*.
- **Verdict**: ≥ 3 informative fragments with ≤ 20 % against give **in cis** / **in trans**; more against,
  **conflicting** (mosaic, paralogue, third haplotype, chimeras); fewer, **undetermined**. The LOD is shown, not
  deciding. A site under 25 % VAF counts only the fragments carrying it; over 75 % it is homozygous.
- **Variant × junction**: the junction's usage on each copy and a Fisher test. **on the variant's copy** when ≤ 20 % of
  the fragments using it carry the reference allele, **on the other copy** when ≤ 20 % carry the variant, else **on
  both copies**. NMD under-counts the degraded copy; short reads link a variant and a junction only within a fragment.
- The panel draws the **map of the two copies**, the reads by haplotag and each fragment; several positions give a
  matrix of pairs. Switch the sample at the top: in a trio, the parents' reads answer for themselves.

## Methylation

Layer **M**, for long-read DNA with `MM`/`ML` tags (ONT, PacBio), views ≤ 2 Mb, reference needed. A panel under each
such track shows:
- a **header**: CpG islands (green) and the 5mC fraction per haplotype, the phased share, CpGs and calls;
- a **ribbon**, blue unmethylated to red methylated, that **splits into HP 1 / HP 2** where each haplotype has ≥ 3
  calls per CpG and the less covered holds ≥ 25 % of the tagged calls, and is one band elsewhere (no heterozygous
  variant, a missing copy, between phase sets). It works within **one phase set** at a time (a strip in its hue, a
  red **⫽** where the set changes); each pixel pools the calls of the CpGs under it, pale under 3 calls;
- the **HP 1 − HP 2 difference** (its sign means nothing across sets) and **allele-specific methylation** framed in
  purple: ≥ 5 consecutive CpGs differing by ≥ 50 points with ≥ 10 calls per haplotype (imprinting, X inactivation);
- **island differences** with ≥ 2 long-read samples: one tag per island, primary vs the others' mean (e.g. **−26 %**),
  blue less, red more, grey under 10 points, bold from 20.

*CpG islands only* restricts every figure to the islands. With the reads track on a view ≤ 30 kb, each read's CpGs
are red or blue, grey below the confidence threshold. Counting follows `modkit pileup --cpg --combine-strands`:
reference CpGs, strands combined, 5hmC not counted as 5mC, primary alignments, and calls under the **10th-percentile
confidence** filtered, the threshold estimated per **100 kb block** (a CpG's value never depends on the zoom).

## Variants of interest and deep links

**Known variants** (app bar): a locus, interval, HGVS g. (`chr17:g.43094464A>G`, `NC_000017.11:g.43094464A>G`) or
VCF-like line, with a label (`BRCA1 p.Glu23Asp`). Each is drawn on every view; *go to…* (Show) centres on one.
Saved with sessions and exports.

**Deep links** open the viewer on a variant from another tool (MobiDetails, a report). Parameters go in the fragment
(`#…`, never sent to a server) or the query string. The page opens on the variant with its marker; BAM/CRAM files
added afterwards appear as tracks there. It works the same from a file share.

```
https://benjamin-cogne.github.io/Sashimi-viewer/#variant=NC_000017.11:g.43094464G>A&label=BRCA1%20c.5266dupC&pad=100
```

| Parameter | Meaning | Default |
|---|---|---|
| `variant` | HGVS g. with `NC_` accession (its version sets the build), `chr17:g.43094464G>A`, `17-43094464-G-A`, `chr17:43094464:G:A` or `chr17:43094464`; comma-separated; ≥ 50 bp deletions/duplications draw as bands | `variant` or `locus` |
| `locus` | Position or interval to show instead | |
| `pad` | Flank in bp on each side | 100 |
| `label` | Text per marker, comma-separated. A `+` before a digit or `*` stays a plus (`c.157+5G>A`); other `+` read as spaces; `%2B` is always a plus | the g. notation |
| `gene` | Symbol used only when no RefSeq gene covers the window | from the position |
| `build` | `GRCh38` / `GRCh37` (`hg38` / `hg19`), when no accession says it | GRCh38 |
| `reads` | `1` opens with the reads track on | off |

Template: `href="…/#variant={{ hg38_g_hgvs | urlencode }}&label={{ c_hgvs | urlencode }}&pad=100"`. A GRCh37
laboratory sends the hg19 accession or `build=GRCh37`; deep-intronic variants need a larger `pad` or a `locus`.

## Sessions and HTML export

**Save session** downloads a JSON file (`sashimi-session-GENE-DATE.json`): the build, the run folder and each file's
path in it, sample names and order, the FASTA, every view with its options, groups, panel order and known variants.
**Load session**: Chrome and Edge reopen the files from the browser's bookmark of the folder after one click; Firefox
and Safari ask for the folder (*Choose the folder*, or drop it). *Open with the files present* starts with a subset;
moved files and a missing FASTA are named. Add files as a **run folder** so sessions can find them.

**Export › Export HTML** writes one page with the data of every view embedded, for readers without the alignments.
The dialog sets the **window** (the view · the view with half its width on each side, recommended · the widest the
viewer loads: 2 Mb of coverage, the reads window for reads) and the **reads per sample** (as displayed, up to 20,000 ·
dense, 100,000 · all). Each view carries its gene models; every sample's coverage, junctions and retention counts
**with the view's own filters** (MAPQ, secondary); the **reads** of views whose reads track is on (mates joined,
clips, haplotags; names replaced by numbers); the **methylation** and reads' CpG calls of views with **M** on; and the
**page's phasing** of views grouped *phased here*. The reader can switch views, zoom and pan within the windows, click
arcs and exons, collapse reads. Not carried: depth-based exon usage (add the BAM/CRAM), reads of views without the
reads track, anything outside the windows (*not in this exported file*). Keep the tab open during the export; it
needs the built page. Format: [`embedded-format.md`](embedded-format.md).

**SVG**: *Export › SVG · N views* stacks every view on one vector page; **⋯** › *Export this view as SVG* saves one.

## Reference builds and lift-over

The *GRCh 38 | 37* switch (GRCh38 by default) must match the alignments. The first file whose header gives its build
(`@SQ AS` or chromosome lengths) sets it while no view is open; a file on another build, or one without annotation
here (T2T-CHM13), raises an amber bar. With views open, switching asks:
- **Keep the views on GRCh38, add the GRCh37 coordinates** (default): a second coordinate line over the ruler, both
  builds at the pointer, and the other build's window in the header (click to copy; ≈ when partly unmapped);
- **Move everything to GRCh37**, for files aligned on it: every view is lifted and reopened.

Lift-over uses Ensembl's assembly map; only the window's coordinates are sent. GTEx is hg38 only.

## Privacy

BAM/CRAM/FASTA files are read with the browser's File API and decoded on the machine
([`@gmod/bam`](https://github.com/GMOD/bam-js), [`@gmod/cram`](https://github.com/GMOD/cram-js),
[`@gmod/indexedfasta`](https://github.com/GMOD/indexedfasta-js)); no byte of them leaves it. Requests go only to the
annotation hosts listed below and carry a gene name, a genomic interval or the letters typed in a search box, never
reads, sample names or file names. The browser keeps the last genes opened, the folded bar, favourite GTEx tissues
and, in Chrome and Edge, a bookmark (a permission, not a copy) of the folders and files added. **Session files and
exported pages are patient data**: a session names files and samples, an exported page holds reads and coverage.

## Deploying in a clinical laboratory

The viewer is a static file: deploy it as you would a PDF.

| Option | How | When |
|---|---|---|
| **Shared drive / intranet** | Copy `sashimi-viewer.html` from a release to a network folder or intranet server | A validated, versioned copy, independent of github.com |
| **GitHub Pages** | *Settings → Pages*: *Deploy from a branch*, or *GitHub Actions* (`pages.yml`) | An always-current URL |
| **Release asset** | `git tag v1.2.0 && git push --tags`: the build workflow attaches the viewer to a release | Traceable versions for a quality system |

**Network**: the browser must reach `api.genome.ucsc.edu` (models, MANE, sequence, protein domains, dbSNP 155 common
SNPs), `rest.ensembl.org` and `grch37.rest.ensembl.org` (fallback, aliases, ENSG ids, lift-over),
`clinicaltables.nlm.nih.gov` (gene-name suggestions) and `gtexportal.org` (GTEx).
**Browsers**: current Edge, Chrome or Firefox. **Limits**: coverage ≤ 5 Mb, variants ≤ 3 Mb, methylation ≤ 2 Mb,
SNPs ≤ 1 Mb, reads window ≤ 2 Mb (100 kb by default), 2,500 reads drawn per track.
**Scope**: splicing, ψ, exon usage, NMD, phasing and methylation figures come from the reads and the model in view;
confirm findings with your validated pipeline before reporting.

## Troubleshooting

| Symptom | Cause → fix |
|---|---|
| "No index found for X" | Drop the `.bai` / `.crai` with its file, named alike (`x.bam` + `x.bam.bai` or `x.bai`) |
| "Gene lookup failed" | Check the symbol and the build, and that the proxy lets the browser reach `api.genome.ucsc.edu` / `rest.ensembl.org` |
| CRAM slow or failing | No reference → add the alignment FASTA (+ `.fai`, + `.gzi` if bgzipped) |
| Coverage but no arcs | Unspliced aligner, sample taken as **DNA** (click its badge), or *min* reads too high |
| Variant sites missing, reference error on the track | No FASTA and a UCSC/Ensembl request failed → move the view to retry, or add the FASTA |
| "zoom in below … to load reads" | View wider than the reads window → zoom in, or raise *window ≤* |
| No methylation panel | Not long-read DNA with `MM`/`ML` tags, view > 2 Mb, or no reference |
| Exon usage without *vs controls* or z | One sample only, or fewer than 5 controls → add comparison samples |
| *Strands* inactive, badge **±?** | Orientation unknown → open a multi-exon gene with ≥ 50 spliced reads, or set the badge |
| *Phase with*: *undetermined* | < 3 informative fragments (sites beyond a fragment, no haplotags, two phase sets) |
| A session finds no files | Moved files, or Firefox/Safari → drop the run folder or *Choose the folder* |
| "not in this exported file" | Outside the exported window → re-export wider, or add the alignments |
| "≈ 1 read in k" on counts | A page exported by an older version that sampled deep windows → re-export from the files |
| No structural arcs | Fewer reads than the *Structure* minimum (3, or 20 on tracks deeper than 200×), or the sample taken as **RNA** → lower *min*, or click its badge |
| "This file is the development entry" | `index.html` opened from disk → open `sashimi-viewer.html` |

## Build from source

Node.js ≥ 20 (`.nvmrc` pins 22).

```bash
npm ci               # exact versions from package-lock.json
npm run build        # tsc + vite → one self-contained sashimi-viewer.html (in dist/ and at the root)
npm run dev          # development server at http://localhost:3000/
npm run typecheck    # tsc only
```

Everything is inlined in one file (`vite-plugin-singlefile`) so it opens from a file share. The root
`sashimi-viewer.html` is tracked: commit it rebuilt with source changes. `build.yml` warns when it is stale and
attaches it to a release on a `v*` tag; `pages.yml` publishes `main`, and `dev` under `/dev/` with a DEV MODE banner.
The page shell and data sources are in `src/standalone/`, the viewer in `src/components/`.

## Citation and licence

Cogné B. *Sashimi viewer: in-browser Sashimi plots from BAM/CRAM files.* CHU Nantes, 2026.
https://github.com/benjamin-cogne/Sashimi-viewer (see [`CITATION.cff`](../CITATION.cff)). The methods, their validation
and references (Sashimi plot, MANE, rMATS, IRFinder, NMD rules, GTEx, WhatsHap, modkit, STRchive, ACMG) are in
[Methods and validation](algorithms.md).

[CC BY-NC 4.0](../LICENSE) © 2026 Benjamin Cogné, CHU Nantes. Free for non-commercial use with attribution, which
covers diagnostic and research use in public hospitals and academic laboratories. Contact the author for commercial
licensing.
