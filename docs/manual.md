<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logo/logo-lockup-dark.svg">
  <img src="logo/logo-lockup.svg" alt="Sashimi viewer" height="64">
</picture>



# User manual

This is the full manual of Sashimi viewer. The short introduction is the [README](../README.md) at the
repository root.

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
   - take [`sashimi-viewer.html`](../sashimi-viewer.html) from this repository (*Download raw file*, keep the `.html` extension).
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
| Is the aberrant junction in *cis* with the variant? | **Reads track** with mismatches (primary sample by default, any sample, or *All samples* for one reads track under each coverage track), then **Collapse**: read-based **phasing** into two haplotypes per phase block (the alleles seen together in the same reads and mates), each drawn with its own junctions and those the two haplotypes use differently flagged (Fisher's exact test), or *Haplotypes: any* for consensus groups (local haplotype × splicing pattern) where the alternate allele seen only in exon-skipping reads is explicit. |
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
- **Very deep libraries** (capture panels at thousands of ×, targeted RNA-seq, highly expressed
  genes): coverage, junctions and intron-retention counts are **exact at any depth**. Every read of
  the window is counted straight from its alignment into per-position tallies, without building a
  read in memory, and nothing is sampled. The visible region is counted first and drawn as it fills
  (*reading the view… N %* next to the sample), then the margins around it (*reading the margins…*),
  as far as about 2 million reads of margin; past that they shrink, and panning out of them costs a
  new read. What has been counted stays with the sample while you stay in that part of the chromosome (a jump
  further than a window away starts over): panning
  back, *Unique reads*, another transcript model or its exon boundaries are answered from the
  tallies at once, without reading the file again (up to 8 Mb per sample and 256 MB for all samples;
  the least recently used are dropped first). A junction is the `N` operation of a read's CIGAR, as
  STAR's `SJ.out.tab`, regtools or pysam count it. The reads track and the exon-usage statistics still
  decode reads, on their own budget (every 2nd, 4th, 8th… read past it, drawn as *downsampled*).
- **Coverage drawing at any zoom**: when a pixel stands for more than one base (a whole gene on
  screen), the line is the exact **mean** depth of that pixel's bases (IGV's default windowing
  function for bigWig tracks; close to UCSC's *mean+whiskers*), over the usual shaded area; a pale
  band above it reaches the pixel's **highest** base.
  Where a base of the pixel falls **below half its mean**, a thin tick runs from the line down to
  that base's depth: a dropout narrower than a pixel stays visible on the whole gene, while the
  base-to-base noise of a shallow library (±20–40 % at 30×) is not drawn as a second line. All three
  come from a pyramid of per-bin minima, maxima and sums, so drawing costs the same whatever the
  zoom, the depth or the number of samples. Zoomed in to the base, the line is the depth itself.
  Reading the tallies back costs what the window holds, not its width: a zoomed-out view where
  only the gene is covered costs about as much as the gene. A wide window that is covered
  everywhere (a genome) is read back in short steps that hand the page back to the browser, and
  the file itself is decoded about 20 000 reads at a time, so zooming out past what has been read
  keeps the page moving while the new part is read.
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
  its share at each intron. When that skipping arc is the only competitor on both sides of the
  skipped exon, the two inclusion junctions show the inclusion level of the event pooled over the
  two introns, (I₁ + I₂) / (I₁ + I₂ + 2·S): the same value on both, the complement of the skipping
  label, rather than their own per-intron shares (which the tooltip still gives). As soon as
  something else competes at one of the two introns, the per-intron shares come back on the labels
  and the pooled inclusion level moves to the tooltip. Every tooltip also gives the event against the canonical junction
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

**Structural hints** (development builds only in 1.1: shown on the dev page and with `?sv=1` on the
URL, hidden on the release until they mature). Where an RNA track shows junction arcs, a DNA track shows the structural
evidence of its reads, drawn as evidence and never as calls: **deletions inside reads** (a CIGAR `D`
run of 50 bp or more) as solid red arcs; **split reads**, read as chains: every part of a read (the
primary alignment and its supplementary alignments, from the SA tag, whichever of them fall in
the window) is placed along the read by its clips, the parts are ordered along the read and each
pair of adjacent parts is one breakpoint, each read counted once. The breakpoint is typed by where
the read continues: further along the same strand, a *deletion*; backwards, a dashed green
*duplication*; on the other strand, a dashed blue *inversion*; on another chromosome, a purple
`→ chr` pill; and an unaligned stretch of the read between two parts adjacent on the reference, a
purple `ins` pill with its size.

**One arc per event.** The reads of one event rarely agree on its breakpoints to the base: long
reads scatter them by a few to a few tens of bases (more in repeats), an inversion shows two
junctions (+ → − and − → +) a few bases apart, and one deletion is a CIGAR `D` in some reads and a
split alignment (hard- or soft-clipped supplementary part) or a clipped read placed by realignment
in others. Arcs of one kind (deletion: CIGAR `D` together with deletion-type split reads and placed
clips; duplication; inversion) are therefore merged when every start and every end of the group
lies within 5 % of the event's length of each other, at least 20 and at most 100 bp. The ceiling
follows the long-read SV callers' ONT settings (cuteSV recommends a 100 bp cluster bias for ONT
deletions; Sniffles2 merges within 150 bp by default). The floor and the relative term keep small
distinct events apart: two 60 bp deletions 40 bp apart stay two arcs. The event's breakpoints are
the read-weighted medians of what it merged, and its pill counts distinct reads (a read showing
both junctions of an inversion counts once), so *Min supporting reads* applies to the event. The
tooltip and the panel give the evidence (CIGAR D, split reads, placed or rescued clipped reads,
the two inversion junctions) and the spread of the merged breakpoints. A deletion is a solid red
arc when only CIGARs carry it, dashed when split or clipped reads support it. The panel matches
the event in the other samples with the same tolerance, and the Groups view pools the members'
events the same way.
**Discordant pairs** (insert size above five times the window's median, at least 1 kb, or both
mates on the same strand) are dashed amber arcs between 500 bp bins, and **soft-clip clusters**
(at least 3 reads clipped by 20 bases or more at one position, split reads excluded) teal pills on
the baseline (⇤ clipped before, ⇥ clipped after). **Clipped reads without an SA tag are placed by
realignment**: the clipped consensus of each cluster is searched on both strands of the window's
reference (its 20 bases next to the breakpoint must occur once, and up to 60 bases must agree with
at most one mismatch per 20); a placed cluster becomes the arc a split read would give between the
clip and the placed sequence, deletion-type, duplication-type or inversion, merged with the
split-read arc at the same breakpoints. The arc shows one count; its tooltip and its panel give the
split reads and the placed clipped reads apart. Hard-clipped records without SA tag, which carry no
sequence, count in the cluster at their clip position and follow it into the arc. Clusters that
find no unique place stay pills. **Clipped reads are then rescued at the breakpoints already
seen**, the way a two-pass aligner uses the junctions of its first pass: a read soft-clipped at an
end of an arc (chain, deletion or placed cluster, within the 5 bp rounding) whose clipped bases,
8 or more (12 when several breakpoints share the position), match the reference at the other end
in the arc's orientation joins the arc; a hard clip without SA tag at an end is attached by
position. Rescued reads never make an arc on their own: the sample needs an aligned split read, a
deletion or a placed cluster there first. The breakpoints of the other DNA samples shown are used
as targets too; a sample's reads rescued at a breakpoint it does not carry itself are listed in the
arc's panel ("no arc") without being drawn, which is how a parent's short clips show up under a
child's breakpoint. The arc tooltip and panel give the aligned, placed and rescued reads apart.
The reference comes from the FASTA at any width, or from the web APIs for windows up to 500 kb. *Min supporting reads* sets the support needed to draw a
hint. Its default is 3 reads, and **20 on a deep track**: a DNA track whose median depth over the
covered bases of the model's coding exons is above 200× (measured once per gene, on the first
complete coverage; over the whole covered window when no model is open). At a few hundred × a handful of reads with a long
clip, a split alignment or an odd insert turns up at many places (PCR chimeras, adapter
read-through, capture edges), so 3 would draw arcs everywhere. The box then shows the value in
effect on the first DNA track, with *auto*; a number typed in it applies to every track, and the
*auto* link beside it returns to the default. An arc's pill shows the read count (≈ on sampled windows); clicking an arc opens a panel with
the count in every DNA sample shown, the median insert size, and a g. notation for a deletion.
Arcs can be hidden, resized and dragged like junction arcs, and DNA groups pool the evidence of
their members. Long reads carry deletions, split reads and clips; short-read pairs add the
discordant pairs (CRAM mate fields are read when the file records them).

**Pairs.** In the reads track the two mates of a pair share one row and are joined by a line
(*Pairs*, on by default, next to *Reads*); the tooltip of a read gives its mate's position and the
insert size. Reads of a **discordant pair** are amber: mate on another chromosome, pair not flagged
as proper by the aligner, or, on genomic DNA, an insert above five times the median of the window.
The pairs are kept in exported pages.

**Clipped, inserted and split-read sequences.** The reads keep the bases their alignment leaves
out: the soft-clipped bases at each end, the inserted bases, the lengths of the hard clips and the
SA tag of a split read, in the live viewer, in exported pages (sections `clips`, `inserts` and `sa`
of the reads stream, `docs/embedded-format.md`) and, once the converter writes them, in containers.
Two options next to *Pairs*, off by default and saved with the session: **Clipped** draws the
soft-clipped bases beyond the ends of the reads, as letters when the zoom allows and as
base-coloured bars otherwise, dimmed where they match the reference continuing past the
alignment, so a real breakpoint sequence stands out from a run of errors or an adapter; the parts
of a split read that fall in the window share one row, joined by a dashed purple line, as mates do,
each part's hard clip drawn as a dashed grey stub of the right length (a hard clip whose other
part is outside the window is only in the tooltip, so a lone supplementary record looks like any
read). **Inserted** writes the inserted bases
inside the insertion marks when there is room; they are always in the tooltip. Both work on RNA
and DNA tracks, on the raw reads only (not on the consensus or haplotype rows). **Click a read**
for its sequence panel: the clipped and inserted sequences with a copy button, the positions of the
other parts of a split read, and, for a supplementary record, a *Fetch from the primary record*
button that reads the read's primary record (the one that carries the whole sequence, soft-clipped)
from the BAM/CRAM and shows the hard-clipped bases, reverse-complemented when the two records are on
opposite strands. **Click a soft-clip cluster pill** (⇤ n, ⇥ n on a DNA coverage track) for the
consensus of the clipped sequences of its reads, anchored at the breakpoint and ready to paste into
BLAT to place the other side of the junction.

**Consensus reads.** *Consensus* (on by default) draws mismatches and indels only where a variant
site is called (at least 3 reads and *Min VAF*), so sequencing errors do not paint every read. It
is offered for every genomic DNA reads track, short reads included, and for long reads whatever
the library; off, every mismatch and indel of every read is drawn. Deletions and insertions of
50 bp or more are always drawn (a deletion as a black line across the gap), called or not: they
are structural evidence, not sequencing noise, and the reads of one large deletion often place
its breakpoint a few bases apart, so no single site would gather them.

**Long reads (ONT, PacBio).** A reads track whose median aligned length is above 1 kb gets one
more noise control. The aligned length counts aligned and deleted bases, not the skipped introns: a
2×100 RNA-seq read spliced over a 3 kb intron is a short read. Counting the intron wrongly took most
RNA-seq tracks of multi-exon genes for long reads, with the long-read thresholds and grouping.
The control is *Min VAF (long)* (20 %), the allele fraction a site needs on long reads, above the
short-read *Min VAF*. It appears next to *Collapse* when such reads are shown and is saved
with the session. Indels of every size are called. Random homopolymer indels stay out of sight
because *Consensus* draws indels only at called sites. A homopolymer indel common enough to be
called is flagged by the variants track's HP check. (An earlier *Min indel* option hid the indels
under 10 bp; it also hid every real short indel of long reads, and was removed.)

**Junction arcs of long reads.** Long reads place a junction a few bases off where the bases next to
it carry errors (ONT especially). One junction then drew as a fan of arcs: the true one, and
neighbours 1–6 bp away with a handful of reads each. On a track whose reads average more than 1 kb
aligned, a junction within 6 bp at both ends of one used by at least 20 times as many reads is
counted with it. Its tooltip says how many reads were counted that way. The ψ values and the
boundary counts use the merged junction. IsoQuant corrects to the annotation within 6 bp on ONT data
(4 on PacBio), FLAIR within 15; here the reads are the only reference.

The 1-in-20 rule is strict on purpose. A variant that creates a cryptic donor or acceptor a few bases
from the canonical one gives a junction at one exact place, used by a real share of the reads (often
a whole haplotype's). Alignment jitter spreads a few per cent over many offsets. A cryptic site used
by at least 5 % as many reads as its neighbour stays its own arc. On a simulated gene, an acceptor
4 bp into an exon, used by 40 % of one haplotype's molecules, stayed apart (137 reads next to 594). It
also formed its own consensus group and was flagged as used by that haplotype only. A looser 1-in-4
rule would have merged it. The known limit: a real site used by fewer than 1 in 20 of its
neighbour's reads is merged into it. Short reads (STAR, HISAT2 place junctions to the base) are
never merged.
Deletions of 50 bp or more inside reads remain structural evidence whatever the setting. On DNA
tracks the reads carry no exon–intron boundary outline (that teal mark is an intron-retention
device for RNA).

**Haplotypes (phasing).** *Collapse* with *Haplotypes: 2* shows the **consensus of each of the
two haplotypes**: one row per haplotype and phase set, drawn like a read. Grey marks the stretches
covered by at least 3 of the haplotype's reads (blank where fewer). Coloured bases, deletion lines
(large ones included) and insertion marks are the variants at least half of its reads carry, so a
heterozygous variant lands on one row and a homozygous one on both. A variant of a haplotype must
also be a site of the window over all reads. Otherwise, where a haplotype has only 3–4 reads (at
the window's edges), two sequencing errors of the same base would make a false variant. Phase sets side by side share a
pair of rows, separated by a dashed line: the haplotypes are linked within a set, not across it
(H1 of one set is unrelated to H1 of the next). Hover a row for its reads, the median assignment
confidence and its variants.

The haplotypes come from one of two sources, chosen with *Phase*:

- **The file's haplotags** (default, used whenever the window has tagged reads). A phasing tool
  wrote them on each read it could place: `HP` (haplotype 1 or 2), `PS` (phase set) and sometimes
  `PC` (Phred-scaled confidence). Tools that write them: WhatsHap or LongPhase `haplotag` (the ONT
  wf-human-variation outputs), PacBio HiPhase, Illumina DRAGEN (TruPath, where `HP` is a "copy
  label" and may go above 2). That phasing comes from the genome's variants, so it reaches across
  what a window's reads cannot link. Untagged reads (no phased variant under them, or fitting both
  haplotypes) are left out and counted.
- **The reads (in-page)**. The window is phased from its own reads, in the spirit of WhatsHap and
  HapCUT2 but in the browser, in milliseconds. Every heterozygous site (25–75 % alternate allele
  among the called sites, so at least 3 reads and *Min VAF*) is linked to the others through the
  **fragments** that cover both, a read and its mate counting as one fragment. Two sites are linked
  when at least 2 fragments cover both and at most 20 % of them disagree. The sites are walked in
  order, and each joins the current **phase block** when its trusted links agree, or starts a new
  one (no linking fragment, or contradicting links, which is what a third haplotype, mosaic alleles
  or errors look like). A single site that does not fit, while the next one does, is taken for an
  outlier and left unphased: an error-made site, or a homozygous one read as heterozygous. The block
  goes on past it instead of ending there. Each fragment then goes to the haplotype of its block that it matches at
  more sites, and the blocks act as phase sets. Short-read pairs phase sites within an insert of
  each other; long reads phase whole windows; RNA-seq phases too, with the fewer heterozygous sites
  its exons carry.

With either source, the heterozygous sites of the window are checked. The two haplotypes of a set
should split them (at least 80 % of one haplotype's reads carrying the alternate allele, at most
20 % of the other's). Sites that do not split are listed on a *not split* row: a mis-phased or
mosaic site, a third haplotype, or a collapsed duplication. Reads whose allele contradicts their
haplotype's are counted in the header (tagging errors, chimeras). With the in-page source, the
sites it could not phase are listed on that row too, with the reason.

**Junctions of each haplotype (RNA).** On RNA-seq, each haplotype row shows its splice junctions. A
line crosses each intron, and an arc spans a junction that skips over others (an exon skip). Each
carries the number of the haplotype's fragments that use it; fainter means used less where the
haplotype decides. For every junction carried by at least 3 fragments of a phase set, the viewer
counts each haplotype's fragments that go another way at one of its ends: another junction from its
donor or to its acceptor, or aligned bases across the exon–intron boundary. The two haplotypes' shares
are compared with Fisher's exact test (two-sided). A junction with p < 0.001 and at least 20 points
between the shares is drawn in the colour of each haplotype and counted in the header ("junctions used
differently by the haplotypes"). Its tooltip gives both shares and the p-value. That is the splice
change of a variant in *cis*: a donor or acceptor variant, or a created cryptic site, acts on its own
haplotype only. Counts are in fragments (a read and its mate are one molecule). On long reads, a
junction placed a few bases off is taken for the common one next to it, as in the consensus groups
below.

**Group by haplotype.** When the reads carry haplotags, *Group: haplotype (HP)* on the raw reads
track packs them per haplotag, like IGV's *Group alignments by tag HP*: HP 1, HP 2 (and higher
copy labels), then the untagged reads. Each group has a label row and a coloured band. Mates and the
parts of a split read stay together. A read's tooltip gives its HP, PS and PC.

*Haplotypes: any* switches the collapsed track to the **consensus groups** of the earlier collapse:
one row per local haplotype × splice pattern with its read count, groups below *Min reads* folded
into a minor bucket. A read and its mate count as one **fragment**, one molecule: together they see
sites and junctions further apart than either read, and *Min reads* counts fragments. Reads fitting
several groups make ambiguous rows, labelled with those groups ("H1|H2 ?"). Two kinds of grouping
are used:
- **Unspliced short reads (DNA)** are grouped by identical patterns of alleles, seeded by patterns
  seen in at least *Min reads* fragments. Every pattern is then placed on the one seed it fits,
  or made ambiguous when it fits several.
- **Long reads, and spliced short reads (RNA-seq)**, are clustered with a tolerance. No two long reads
  carry the same pattern: each covers its own run of sites, with its own sequencing errors. On
  RNA-seq, short fragments see different pieces of the same transcript (exons 2–3, 3–4, …), which as
  exact patterns split one isoform of one haplotype into dozens of rows.
  - The fragments are taken left to right. Each joins the group whose consensus it agrees with at
    their shared sites, up to 20 % of them disagreeing, and whose splicing it never contradicts:
    - a junction of one crossed by aligned bases of the other (a skipped exon, a retained intron);
    - two different junctions overlapping (another donor or acceptor).
  - A shared junction counts as agreement. A fragment with calls meets groups through its sites only,
    since junctions that both haplotypes splice say nothing of the phase.
  - A fragment agreeing with two groups goes to the one it agrees with by 2 more sites or junctions,
    and is ambiguous otherwise.
  - Groups are then merged when they agree where they overlap: at least 2 shared sites, or at least
    2 fragments linking them (a mate pair that agrees with each at a site of its own). A group that
    could join two groups that differ from each other (two isoforms, or two haplotypes where it sees
    no site of either) stays apart, because the reads cannot say which.
  - Every fragment is finally placed again against the final groups.
  - Homozygous sites (above 90 % alternate) are shown but not used, because every read agrees there.
  - On long reads, a junction within 6 bp of one used by at least 20 times as many reads is taken for
    it, as on the junction arcs (see *Long reads* above).
  - When most fragments agree with nobody (a wrong reference, a junk window), the exact patterns are
    used instead.

A diploid DNA window thus gives its two haplotypes as two groups, whatever the number of groups the
other alleles (mosaic, a paralogue, a third copy) add. An RNA window gives one group per haplotype ×
isoform where the reads show both. On a simulated gene (9 exons, one haplotype skipping exon 5 in
cis with its alleles, an alternative acceptor on both), every group was pure:
- **ONT-like cDNA** (4 % errors, 10 % of junctions shifted by 1–4 bp): the five haplotype ×
  isoform classes as five groups.
- **2×100 pairs**: the exon-skipping reads as one group of the alternate haplotype. The pieces the
  reads cannot assign to one isoform are rows of their own, or ambiguous rows.

**Allele fractions of the reads track.** Every alternate base counts in a site's fraction,
whatever its base quality, as in the variants track. Before, the bases under Q20 were left out of
the count but not of the depth. On ONT data, where a third of the bases can be under Q20, every
fraction came out a third too low. Homozygous sites (about 0.67) were then taken for heterozygous
ones, and heterozygous ones (about 0.33) fell under the phasing window. The phasing itself still
uses only alleles of Q20 or more.

**Speed.** A collapsed window (up to 40,000 reads decoded, then phased or grouped) is computed in
the background worker, and only its answer, not the reads, comes back to the page. A 1,000×
capture, which used to freeze the page for half a second while collapsing, no longer holds it up.
Alleles are kept per fragment only at the sites it covers, and each fragment is compared only with
the groups and phase blocks it overlaps. A window of thousands of apparent sites (reads against the
wrong reference) used to take 1–2 minutes and gigabytes to collapse, and over ten minutes on 1.5
million reads. It now takes about 6 s. These choices are saved with the session. Exported pages keep the haplotags of
their embedded reads and compute the haplotypes on the spot.

**Layers.** Four layers make a DNA track, switched in the toolbar by one colour-coded control,
**C · V · M · R**. A filled letter is a layer that is on; a click switches it for every DNA sample
at once. They are always drawn in this order, top to bottom:
- **C, coverage** (blue, on by default): the coverage histogram and the structural arcs over it.
  Switched off, a DNA track keeps its label band and the layers under it, for example to compare
  the methylation or the variants of several samples in little room. The coverage stays loaded,
  so switching it back on is immediate. RNA tracks always show their coverage, since the sashimi
  plot is drawn on it;
- **V, variants** (amber): the variants track described below;
- **M, methylation** (red): CpG methylation of long reads (below). *CpG islands only* appears next
  to the control while it is on;
- **R, reads** (slate): the alignments. The sample chip, or the selector next to the control,
  picks which sample's reads are shown, or all of them.

A page without a DNA sample shows a plain *Reads* option instead.

Why global toggles: the same question is usually asked of every sample at once (patient against
controls), and switching a layer off stops its reading and releases what it held, so the page only
ever carries what is shown.

**Variants track.** With *Variants* on, the reads of each DNA track's window are scanned in the
background (a Web Worker), for views up to 3 Mb. The scan counts every read, whatever the width and
depth. A site needs at least 3 supporting reads and an alternate-allele fraction of *Min VAF* or
more (*Min VAF (long)* for long reads). Every alternate base counts in that fraction, whatever its
base quality: the BQ check below says how much of it rests on low-quality bases, rather than
hiding them. The sites then
follow the window: moving or widening it scans only what is new, and zooming back in, going back,
or another *Min VAF* are answered at once. Each site is drawn in a track of its own under the
coverage, never on it:
- a **bar** as high as its alternate-allele fraction, on a 0–100 % axis with a 50 % guide, in the
  colour of the alternate base (purple for an insertion, black for a deletion). The grey above it
  is the reference share. The fraction is written over the bar where there is room;
- four **quality cells** under the bar, green (pass), amber (check), red (likely artefact) or grey
  (too few reads to judge). Where sites are too close for four cells, a single cell takes the
  colour of the worst check. From left to right:
  1. **BQ** (SNV): the share of the alternate bases seen with a base quality under 20. They count in
     the allele fraction; when they are most of the allele, it may be sequencing errors (a noisy
     cycle, a homopolymer tail) rather than an allele. Amber from 25 %, red from 50 %. For an indel the first cell is **HP**, the
     length of the reference homopolymer at the site, where polymerase slippage (short reads) and
     basecalling (ONT) make indel errors. Short reads: amber from 6, red from 10. Long reads: amber
     from 4, red from 7.
  2. **MQ**: the share of the alternate reads mapped with MAPQ < 20, against the other reads over
     the site. Red when half of them or more map poorly while the other reads do not (30 points
     more); amber from 20 % (15 points more), or when every read maps poorly there (a repeat).
  3. **SB** (strand bias): the alternate reads' strand split against the other reads' (Fisher's
     exact test, two-sided, as GATK's FisherStrand). Amber at p < 0.001, red at p < 1e-4. An allele
     carried on one strand only is the signature of oxidative damage (8-oxoG, G>T) or of a library
     or PCR artefact.
  4. **END** (read-position bias): the alternate calls within 10 bases of an alignment end against
     the other reads' (Fisher, one-sided). Amber at p < 0.001, red at p < 1e-4. It catches
     misaligned read ends near an indel, and adapter or clip artefacts.

  Both tests compare two sets of counts, so a reference side of a few reads, itself noisy, weighs
  as little as it should. The cut-offs are strict because a window holds hundreds of sites: at
  p < 0.01 a few would be coloured by chance alone. On clean synthetic data (199 sites, 40×) none
  is flagged. GATK's hard filter for SNVs (FS > 60) sits near p = 1e-6.

Testing each share against the same reads' reference side keeps a library whose reads all run one
way (amplicons), or a region where every read maps poorly, from being called an artefact. The
checks follow what GATK's hard filters (FisherStrand, MQ and MQRankSum, ReadPosRankSum) and a
reviewer in IGV look at. They are an aid to judge a call, not a variant caller's filter.

- under them, the **allele balance** strip: the major allele fraction, max(VAF, 1 − VAF), smoothed
  along the window. Each pixel pools the sites under it, or its 8 nearest within 100 kb, leaving out
  sites under 20 % and those flagged red. Where heterozygous sites (VAF 0.2–0.8) are among them, the
  strip takes their median:
  - 0.5, in green, is a balanced diploid region;
  - it turns amber, then red, as one allele takes over. A copy gain gives 0.67, and a mosaic change
    anything in between.

  Where hardly any site is heterozygous, and a wider stretch agrees (at least 15 of the 20 nearest
  sites within 500 kb, under 10 % heterozygous), the region is a **run of homozygosity** in red:
  loss of heterozygosity, uniparental disomy or identity by descent. Homozygous sites alone never
  make a region red. In a normal genome about 40 % of the variant sites are homozygous for the
  alternate allele, which is why the median of the heterozygous sites gives the balance, and not a
  mean over all. A thin line gives the value, 0.5 at the bottom to 1 at the top, and the hover card
  gives the sites pooled. Twenty sites is far fewer than PLINK's default run (100 SNPs over 1 Mb on
  arrays), so take a red stretch as a lead to check, for example against a copy-number view.

Hovering a site gives its values. A click opens the **distributions from the reads** over the
site, decoded there (up to 5,000):
- the alternate and reference reads by strand;
- their mapping qualities;
- the alternate bases' qualities;
- their distances to the nearer read end.

Distributions are compared as shares, so a skew shows whatever the depth. The track's header gives
the sites of the view and how many pass, need a check or are flagged.

The scan counts this evidence without holding a record per error. The first two calls of each
(position, base) are kept as flags in a byte, and only an allele seen a third time gets an entry.
On a noisy long-read library most mismatches are isolated errors, so the added memory stays
moderate.

**Allele balance.** *Common SNPs* stay off by default on DNA tracks as on RNA ones; switch them on
to separate a known polymorphism from a novel change. With Variants on, the track label then
summarises the allele balance of the common SNPs covered: how many are heterozygous
(0.2 ≤ VAF ≤ 0.8) and the range of their fractions around 0.5. The label flags two patterns:
- *allele imbalance?*: heterozygous fractions far from 0.5 (median deviation above 0.15 over at
  least 5 SNPs), as in a mosaic copy change, loss of heterozygosity or contamination;
- *no heterozygous SNP (LOH / UPD?)*: a window with at least 8 homozygous common SNPs and none
  heterozygous.

When every shown
sample is DNA the axis keeps the genomic orientation, coordinates increasing to the right, even
for a minus-strand gene.

**CpG methylation (long reads).** ONT (dorado, Guppy) and PacBio (jasmine, primrose) basecallers
write the 5-methylcytosine calls of each read in its base-modification tags, `MM` (which bases carry
a call) and `ML` (the probability of each call, 0–255), as the SAM specification describes (the
legacy `Mm` / `Ml` names are read too). Switch on **Methylation** in the options panel, shown when a
DNA track is open. Each long-read DNA track then gets a panel under its coverage:

- a **header** with the CpG islands of the reference (green; at least 200 bp with GC ≥ 50 % and an
  observed/expected CpG ratio ≥ 0.6, Gardiner-Garden & Frommer 1987) and, on the right, the 5mC
  fraction of the view per haplotype, the share of the view that is phased, the CpGs and the calls
  counted (its tooltip gives the filter);
- a **ribbon** coloured from blue (unmethylated) to red (methylated), the same convention as IGV.
  Where the reads are **phased** it splits in two, HP 1 on top and HP 2 below; where they are not,
  it joins into one band of **all reads**. A thin vertical line marks each split and join. Unphased
  stretches include regions with no heterozygous variant to tag the reads by, homozygous or
  hemizygous stretches, and the gaps between phase sets. A pixel is split when both haplotypes carry a
  fair share of its calls: each at least 20 % of them, the tagged reads together at least 60 %, and
  each haplotype 3 calls or more. Split stretches narrower than 4 pixels are joined.
  A file without haplotags gives one band all along;
- a **regional density at every zoom**. Each pixel shows the pooled fraction of methylated calls
  over the CpGs under it: the calls add up, and per-CpG fractions are not averaged, so better-covered
  sites weigh more. When fewer than 6 CpGs lie under the pixel, it pools the 6 nearest, as long as
  they are within 1 kb of it. Zoomed out, a pixel is exactly its stretch of genome; zoomed in, the
  window stays a few CpGs wide and narrows as you zoom, so the panel keeps reading as a density
  instead of isolated sites. A CpG desert wider than that stays a gap. A pixel with fewer than
  3 calls is drawn pale. The hover card gives the fractions of the pixel under the pointer, split or
  joined as drawn, with the CpGs and the span pooled;
- under the ribbon, the **difference HP 1 − HP 2** where it is split, as bars up (HP 1 more
  methylated, in its colour) or down (HP 2). **Allele-specific methylation** is framed in purple
  across the panel. A frame marks a run of at least 5 consecutive CpGs covered on both haplotypes
  whose fractions differ by 50 points or more, with at least 10 calls on each side, and its label
  gives the difference. Imprinted differentially methylated regions, the inactive X of a female
  sample, allele-specific promoters and cis-acting variants look like this.

The calls follow the conventions of `modkit pileup --cpg --combine-strands` (Oxford Nanopore's
reference tool; `pb-CpG-tools` for PacBio does the same per haplotype):
- only the **CpG sites of the reference** are counted, never every C of the reads. Both strands are
  combined: a − strand read's call, on the G of the CpG, is counted at the C of the + strand;
- a call is **methylated when 5mC is the likelier state**. 5hmC, when called, is not counted as 5mC:
  its probability is set aside and the other two renormalised (modkit's *traditional* preset);
- calls whose confidence is below the **10th percentile** of the window's calls are **filtered
  out**, as modkit does by default. The threshold and the calls left out are in the header's
  tooltip;
- the tags describe the read as sequenced: a hard-clipped record is used only when its `MN` tag
  confirms its sequence matches them, and a record whose sequence no longer matches is skipped.
  Secondary, duplicate and QC-failed records are left out as elsewhere.

*mosdepth*, often cited for long-read depth, has no methylation mode; its role here is the coverage
track. The counting runs **in the background** (the Web Worker of the variant scans): the page never
waits for it. The counts are kept per sample over a window of the view plus half of it on each
side, so panning and zooming inside it cost nothing new. Drawing is computed per pixel, not per
read or CpG, so it takes the same time at 1 kb or at 2 Mb. It is counted for views up to 2 Mb
(wider, the panel asks to zoom in) and needs the reference sequence.

With the **reads track** open on a view of 30 kb or less, each read's CpG calls are drawn on it:
red for methylated, blue for unmethylated. Calls below the confidence threshold are left uncoloured,
so the read's grey reads as "no confident call". *Group: haplotype (HP)* then shows the two
haplotypes' reads apart, and an allele-specific region stands out as a red block over a blue one.
Wider than about 15 kb (under 0.08 pixel per base), a read's consecutive calls of one kind (within
100 bp of each other) are drawn as one bar rather than one line each. The panel is not part of exported pages, which carry no
base-modification tags.

**CpG islands only.** Once *Methylation* is on, a second option, *CpG islands only*, restricts
the panel to the CpG islands of the reference. Only their CpGs are counted in the ribbon, the
haplotype difference, the allele-specific frames and the figures of the header. A pixel's smoothing
never reaches outside its island, and the panel stays empty between islands. Promoter islands are
where methylation carries most of its regulatory and diagnostic meaning: an imprinting centre, a
hypermethylated tumour-suppressor promoter, a fragile-site expansion such as *FMR1*. Leaving out
the gene-body CpGs (mostly methylated) makes the islands stand out.

**Island differences between samples.** With two or more long-read DNA samples open, each panel
gets a row of tags, one per CpG island, whichever mode the panel is in:
- on the **primary sample** (the first one), its mean island methylation minus the mean of the
  other samples, for example **−26 %** for an island 26 points less methylated than in the others;
- on **each other sample**, its own mean minus the primary's.

A tag is coloured blue (less methylated) or red (more), and grey below 10 points. From 20 points it
is bold, and on the primary it is filled. The mean is taken over the island's CpGs covered by at least
5 calls **in both samples**. These CpGs are paired, so a CpG covered in only one sample cannot tilt
the comparison, and an island needs at least 3 of them. Each CpG counts once, as the mean of the
CpGs' 5mC fractions (the usual mean β of array and bisulfite analyses), whatever its depth. When
tags would overlap, the largest differences are kept. The tooltip lists each sample's mean, its
difference and the CpGs compared. The differences are computed when the counts change, not while
the view moves.

**Staying fluid.** Methylation is designed to stay fluid:
- **Counting in the background.** It runs in the Web Worker, and the panel is computed per pixel.
- **Moving read calls.** Each read's calls are built once, as paths relative to the read, and a pan
  or zoom only moves them. With reads, variants and methylation all on, frames stay at the
  display's rate on 30× long reads, like frames with all three off.
- **Deep windows.** A deep window can hold more than 400 coloured reads, for example 120× of
  5 kb reads. There, the calls on the reads are left out while the view moves and come back
  200 ms after it stops.

Switching an option off releases what it holds:
- **Methylation:** the counts of every sample, here and in the worker, and the calls carried by the
  reads;
- **Variants:** the scanned sites and the allele counts of the worker;
- **Reads:** the reads of the reads track, and the records the alignment libraries decoded for it.

They are read again when switched back on.

**Memory.** The alignment libraries keep decoded records so that a genome browser can pan without
decompressing again. Their default is 1 GB per file. This page counts coverage, variants and
methylation once into compact states of its own, so it gives the libraries one budget of 128 MB
for all files together, in the page and in the background worker each. A record nothing has used
for 45 s is dropped, and the records are let go as soon as what they were decoded for is counted.
On a 1,000× capture, the page's memory with coverage alone went from 396 MB to 121 MB. The
decompressor keeps a working area as large as the largest block it has inflated (about 100 MB on
such data); this is a one-off ceiling, not a growth. The reads track draws mismatches as one shape
per base colour, not one per mismatch. A noisy long-read library with *Consensus* off, or a
reference of another genome build, used to put over 100,000 shapes on the page; it now stays
under 3,000.
## Moving inside a view

The search box next to the gene name moves the window without leaving the view. It takes, in
this order:

- **Coordinates** — `chr17:43,094,464` (a 1 kb window) or `chr17:43,000,000-43,100,000`. On
  another chromosome the gene at the locus is opened.
- **A c. or n. position** on the transcript drawn — `c.234`, `c.-12` (5′ UTR), `c.*30` (3′ UTR),
  `c.288+1` and `c.289-2` (intronic, counted from the splice site), `n.412` for a non-coding
  model. A range is written HGVS-style with an underscore: `c.234_267`. A whole variant
  description works too — `c.234A>G`, `c.123_125del`, `c.2033dup` — the window moves to the
  position and the change itself is ignored, so a c. copied from a report or from this viewer's
  own arc tooltips can be pasted straight in.
- **An exon number** — `12`, `exon 12`, or `exons 3-5` for a run. Exons are numbered in
  transcription order, as they are drawn, so exon 1 is the 5′ one on either strand. The window
  frames the exon with enough intron on each side to show both splice sites.
- **A gene symbol or an ENSG id**, which opens that gene in the view.

c. positions and exon numbers are read on the model currently displayed, so they follow the
transcript you chose in *All transcripts*, and a number that does not exist on it is refused with
the reason: a c. past the stop codon, an intronic offset on a base that is not a splice site
(usually the sign that the number was written against another transcript), an exon number the
model does not have. The position reached is highlighted like a searched locus — a dashed line
with its coordinates pinned above the ruler — and right-clicking the highlight removes it.

The header box above the panel opens genes and loci into **views** and does not take c. positions
or exon numbers; typing one there says which box to use.

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

**Hide panel** (next to *+ Files…*) folds the upper panel away: the notes, build, files, sample
chips, session buttons and known variants disappear and one bar stays, with the logo, a **Show
panel** button, the *Views* tabs and the search box, so a new gene or locus can still be opened
while the plot gets the rest of the window height. The choice is remembered by the browser.

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
sample, with the reference bases, the mismatches and the mate positions, so the recipient sees the
same pile-up with its pairs, can switch between reads, phased haplotypes and consensus groups and
change *Min VAF*. Read names are replaced by numbers.

**Size.** Reads and coverage are stored in a compact binary form (columns of small integers,
compressed with the browser's own deflate, written in base64): about 6 bytes per read with its
mate, 24 times less than the JSON of earlier exports, and 5 times less for the coverage. The
format is specified in [`docs/embedded-format.md`](embedded-format.md) so other tools can
write or read it. A page decodes a view the first time it is shown, off the main thread, the
active view first and the others in idle moments, so opening a page with many views costs
nothing up front. Browsers without a built-in decompressor (before Chrome 80, Firefox 113,
Safari 16.4) use a small JavaScript one bundled in the viewer. Pages exported by earlier versions
still open. While an export runs (HTML or SVG), a progress window names the step in progress (view, sample,
coverage or reads) with a bar over the total number of steps; the download starts when it
closes, so the tab should stay open. The dialog that opens on *Export HTML* chooses the window exported around each view, for the
coverage, junctions and retention counts as much as for the reads (the view as shown; the view with
a half-width margin on each side, the default; or the widest window the viewer itself loads, up to
2 Mb for coverage and 100 kb for reads) and the number of reads per sample (up to 20,000 as
displayed, about 150 kB per sample and view; up to 100,000 for deep windows, under 1 MB; or every
read of the window, a few MB for a deep window). Views wider than 100 kb have no reads track and
are skipped.

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
  junctionSnap.ts                long-read junctions a few bases off a much more common one, merged into it
  phasing.ts                     read-based phasing of the heterozygous sites into two-haplotype blocks
  haplotypes.ts                  haplotype consensus rows from the HP/PS haplotags or the in-page phasing
  svmerge.ts                     structural arcs with nearby breakpoints merged into events
  alleles.ts                     allele counts of the full variant scan, with each call's quality evidence
  methylation.ts                 CpG methylation from the MM / ML tags: counts per CpG × haplotype, filter threshold, CpG islands
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
  sashimi/siteQuality.ts         quality checks of a variant site (BQ / homopolymer, MQ, strand, read position)
  sashimi/datasource.ts, types.ts   interfaces shared with the parent application
```

</details>

## Citation

If the viewer contributes to a publication, please cite it
(see [`CITATION.cff`](../CITATION.cff); GitHub's *Cite this repository* button formats it):

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
- Gardiner-Garden M, Frommer M. *CpG islands in vertebrate genomes.* J Mol Biol 1987;196:261–282
  (CpG island criteria).
- The SAM/BAM Format Specification Working Group. *Sequence Alignment/Map Optional Fields
  Specification*, section 1.7 "Base modifications" (MM, ML, MN tags).
  https://samtools.github.io/hts-specs/SAMtags.pdf
- Oxford Nanopore Technologies. *modkit* documentation: `pileup`, `--cpg`, `--combine-strands`,
  filter threshold (10th percentile), 5hmC handling. https://nanoporetech.github.io/modkit/
- PacBio. *pb-CpG-tools*: CpG methylation probabilities from HiFi reads, per haplotype.
  https://github.com/PacificBiosciences/pb-CpG-tools
- GMOD JavaScript libraries: [bam-js](https://github.com/GMOD/bam-js),
  [cram-js](https://github.com/GMOD/cram-js), [indexedfasta-js](https://github.com/GMOD/indexedfasta-js).

## License

[CC BY-NC 4.0](../LICENSE) © 2026 Benjamin Cogné, CHU Nantes. Free for non-commercial use with
attribution, which covers diagnostic and research use in public hospitals and academic
laboratories. Contact the author for commercial licensing.
