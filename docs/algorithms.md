# Methods and validation

This document describes the methods behind the figures Sashimi viewer 1.2.0 draws, for reviewers, the validation and
quality files of a clinical laboratory, and bioinformaticians. It complements the user manual
([`docs/manual.md`](manual.md)), which says how to use each feature; here each section gives the rules, thresholds and
formulas, and names the source files that implement them. **The code is the reference**: where this text and the code
differ, the code is right, and every number below can be checked in the file named at the head of its section
(paths relative to the repository root, e.g. `src/standalone/methylation.ts`). Coordinates are 0-based half-open
unless written as `chr:pos` (1-based), as in the code.

## 1. Reads and filters

`src/standalone/alignments.ts`, `src/standalone/coverage.ts`, `src/standalone/localSource.ts`, `src/standalone/mates.ts`

**Flags.** BAM and CRAM files are decoded in the browser with the GMOD libraries [31]. As `samtools` [43] does by
default, unmapped (0x4), secondary (0x100), QC-failed (0x200) and duplicate (0x400) records are left out everywhere;
supplementary records (0x800) are kept. The *Secondary alignments* switch re-admits 0x100 records in coverage,
junctions, variant sites and the reads track. Methylation and the repeat inspector use primary records only.

**Minimum MAPQ.** The filter takes five values, **0 (all), 1, 10, 20, 30**, and keeps a record when its MAPQ is at or
above the threshold, as `samtools view -q`. A typed value is read as the largest of the five not above it. It applies to
coverage, junctions, boundary counts, variant sites, structural evidence, exon-depth usage and the reads track; it does
not apply to methylation. The scale is the aligner's: STAR [42] gives 255 to a unique read and 3, 1, 0 to a read placed
at 2, 3–4, 5 or more loci; HISAT2 60, 1, 0; bwa and minimap2 0–60, a Phred-scaled probability that the placement is
wrong. MAPQ 255 is therefore treated as ≥ 30. Changing the filter costs no new read of the file for coverage and
variant counts: each read under MAPQ 30 is also counted in its bucket (MAPQ 0, 1–9, 10–19, 20–29), and a
threshold subtracts the buckets under it from the total. In the reads track, a read under MAPQ 10 is drawn hollow
(display only).

**Mates and fragments** (`mates.ts`). A read and its mate form one *fragment* (one molecule) for phasing and consensus
groups. Among reads that point at each other's start with the other pair bit, the mate is (1) the read of the same name,
which makes a BAM's or CRAM's pairing exact, else (2) the one with the opposite template length, else (3) the first by
position. Exported pages number their reads (`read N`) and store mate links (`columnar.ts`), joined before any rule.

**Long reads.** A set of reads is taken for long reads when its **median** aligned length (aligned bases only, not
introns `N` nor deletions `D`) is above 1 kb (`readsWindow.ts`, `isLongRead`). A 2×100 RNA-seq read spliced over a 3 kb
intron, or a 2×150 read carrying an 8 kb deletion, stays a short read. The coverage counts use the **mean** aligned
bases per read over their tallies with the same 1 kb cut-off (`coverage.ts`, `LONG_READ_BASES`).

## 2. Library type and strandedness

`src/standalone/localSource.ts` (`classifyHeader`), `src/standalone/main.tsx`, `src/components/SashimiViewer.tsx`,
`src/standalone/alignments.ts`

**RNA or DNA**, decided in this order:

1. The `@PG` lines (a DRAGEN *Hash Table Build* line, which describes the reference index, is ignored). STAR, HISAT(2),
   TopHat, Subjunc, GSNAP, OLego, MapSplice, CRAC, `minimap2 -x splice`, DRAGEN `--enable-rna true` mean RNA; bwa(-mem2),
   bowtie(2), other minimap2 presets, Isaac, NovoAlign, NGMLR, Winnowmap, pbmm2 and other DRAGEN runs mean DNA.
2. An Oxford Nanopore basecalling model (read group `DS:basecall_model=…` or the command line): `rna00x` is direct RNA,
   `dna_r…` DNA. A cDNA run is basecalled with a DNA model and is told apart by its aligner or its reads.
3. A file whose `@PG` names MinKNOW, dorado or guppy and no spliced aligner is taken as DNA.
4. The reads of the first multi-exon gene opened with **≥ 200 reads**, which override the header: **≥ 2 %** of reads
   spliced (a CIGAR `N`) means RNA, **< 0.2 %** DNA, in between no change. A choice made by the user is never overridden.

An undetermined sample is handled as RNA.

**Strandedness (chip).** From the spliced reads of the gene's own annotated introns, once **≥ 50** of them are counted:
the share lying on the gene's strand under the reverse (dUTP, `fr-firststrand`) rule, where read 2 of a pair lies on
the transcript strand and read 1 or a single read on the other one. **≥ 90 %** is *reverse*, **≤ 10 %** *forward*,
anything between *unstranded*, as RSeQC's `infer_experiment.py` [36] reads it. The call made on the most reads stands
unless the user sets one.

**Strandedness (exon-depth usage).** The depth-based exon usage (§4.4) makes its own call from the reads of its exon
scans: the first 3,000 of them, ≥ 100 needed, all reads (spliced or not), with the same 90 % / 10 % cut-offs
(`detectStrandness`); a stranded library then counts only the gene-strand reads.

**Wrong-strand reads.** An opposite-strand junction holding **< 5 %** of the reads of an intron that has ≥ 20 reads on
the gene's strand is labelled as the gene's own reads on the wrong strand: stranded protocols leave about 0.5–3 % there
(0.47–0.63 % for the best protocols [12]; 0.6–3 % of the sense signal in dUTP libraries [13]).

## 3. Coverage, junctions and boundary counts

`src/standalone/coverage.ts`, `src/standalone/localSource.ts`, `src/standalone/junctionSnap.ts`,
`src/components/sashimi/geometry.ts`

### 3.1 Exact tallies

Coverage, junction counts and intron-retention counts are **exact at any depth**: no read is sampled. Each aligned block
(CIGAR `M`, `=`, `X`) of each record goes, from the packed CIGAR and without building a read object, into two sparse
histograms (4 kb chunks, allocated only where something lands): block starts *S* and block ends *E*. Then:

- depth(*x*) = *S*(≤ *x*) − *E*(≤ *x*). Deleted bases (`D`) and introns (`N`) are not depth, as in `samtools depth`
  without `-J`;
- reads unspliced through a window [*A*, *B*) = blocks with start ≤ *A* and end ≥ *B* = *S*(≤ *A*) − *E*(< *B*) + the
  blocks strictly inside (*A*, *B*); the last term involves only blocks of ≤ 16 bases, kept apart in a sorted list;
- **a junction is the CIGAR `N` operation itself**, as STAR's `SJ.out.tab`, regtools and pysam count it. A deletion
  next to an intron (`…M D N…` or `…N D M…`) is not part of the junction; two deletions in a row make no junction.
  The same rule holds for junctions read from decoded reads (`collapse.ts`, `readJunctions`).

A sample keeps its tallies for one stretch of a chromosome: reads *starting* in [*ps*, *pe*) are *owned*, reads
starting before *ps* and reaching into it form a *spill*, so the counts are exact over [*ps*, *pe*) whatever the read or
intron length. Panning back, another MAPQ threshold, another transcript model or its exon boundaries are answered from
the tallies without reading the file again. With *Strands*, the reads on + under the reverse rule are counted again in
a `plus` layer; the − strand is all − plus. The view is counted first and drawn as it fills, then the margins, up to
about 2 million reads by the index's estimate (halved while the window looks too deep, dropped under 2 kb).

### 3.2 Boundary-spanning reads

A read counts as unspliced through an exon–intron boundary when one aligned block covers at least **6 bases on the exon
side and 10 on the intron side** (`SPAN_EXON_ANCHOR`, `SPAN_INTRON_ANCHOR`): for an intron start *p*, a block over
[*p* − 6, *p* + 10); for an intron end *q*, over [*q* − 10, *q* + 6). These counts feed intron retention (§4.2).

### 3.3 Long-read junction merge

`junctionSnap.ts`. Long reads place a junction a few bases off where the bases next to it carry errors (ONT
especially), so one junction draws as a fan: the true one and neighbours 1–6 bp away with a handful of reads each. On
long reads only, junctions are taken in decreasing order of reads, and a junction **within 6 bp at both ends** of one
used by **at least 20 times as many reads** is counted with it (`JUNCTION_SNAP_BP = 6`, `JUNCTION_SNAP_RATIO = 20`).
The merged arc's tooltip gives how many reads were counted that way; ψ, usage and boundary counts use the merged
junction. For comparison, IsoQuant matches and corrects junctions to the annotation within 6 bp on ONT data (4 on
PacBio; `--delta`) [34] and FLAIR within 15 bp [35]; here there is no annotation, only the reads.

The ratio is strict on purpose: a cryptic site created by a variant gives a junction at one exact place used by a real
share of the reads (often a whole haplotype's), while jitter spreads a few per cent over many offsets. A site used by
**≥ 5 %** as many reads as its neighbour keeps its own arc; **known limit**: a real site under 1 in 20 is merged.
**Short-read junctions are never merged** (STAR and HISAT2 place junctions to the base), in the arcs, the consensus
groups, the haplotype junctions and *Phase with* alike.

### 3.4 Drawing: the min/max pyramid

`geometry.ts` (`DepthIndexBuilder`, `depthRange`, `buildCoveragePaths`). Depth is kept as runs of equal depth. A
pyramid holds, for bins of 16·2^*k* bases, the **minimum, maximum and sum** of depth over each bin, built in one pass as
the runs are swept. Any range then decomposes into O(log *n*) whole bins plus two ragged edges read from the runs, so
the exact minimum, maximum and mean of every pixel column cost the same at any zoom, depth or number of samples (the
principle of bigWig zoom levels and of the M4 aggregation [39]). When a pixel stands for more than one base:

- the line is the **mean** depth of the pixel's bases (IGV's default windowing for bigWig [38]);
- a pale band reaches the pixel's **highest** base;
- a tick runs down to the **lowest** base where it is below **half the mean**, so a dropout narrower than a pixel stays
  visible, while the base-to-base noise of a shallow library (±20–40 % at 30×) is not drawn as a second line.

## 4. Splicing measures

`src/components/sashimi/aggregate.ts`, `src/components/sashimi/geometry.ts`, `src/components/sashimi/retention.ts`,
`src/components/sashimi/spliceModel.ts`, `src/components/SashimiViewer.tsx`

The Sashimi plot follows MISO [1] and ggsashimi [2]. The reference model is MANE Select [3] (RefSeq Select, else the
longest CDS, when there is none); the user can make any listed transcript the reference.

### 4.1 % usage (Usage mode, Groups view)

At every intron of the reference model, the events using its donor or its acceptor compete, and each is labelled with
its weight over the sum of the weights at that intron, so the labels of one intron add up to 100 %:

| Event | Weight |
|---|---|
| Canonical junction | its reads *C* |
| Alternative 5′ or 3′ site | its reads *n* |
| Pseudo-exon: an alternative-3′ arc *A* into a cryptic exon of ≤ 500 bp and an alternative-5′ arc *B* out of it | (*A* + *B*) / 2, same value on both arcs |
| Exon skipping (spans two introns) | 2·*S* over the two introns' totals |
| Intron retention (with *Retention* on) | (*R5* + *R3*) / 2, reads unspliced through the two boundaries (§3.2) |

An exon-skipping label is the rMATS SE value [5] 2·*S* / (*I*₁ + *I*₂ + 2·*S*) when nothing else competes at those
introns; the two inclusion junctions then show the pooled inclusion level (*I*₁ + *I*₂) / (*I*₁ + *I*₂ + 2·*S*), the
complement of the skipping label. As soon as something else competes at one of the introns, the per-intron shares come
back on the labels. Tooltips also give each event against the canonical form alone, rMATS-style: *n* / (*n* + *C*),
(*A* + *B*) / (*A* + *B* + 2·*C*), 2·*S* / (*I*₁ + *I*₂ + 2·*S*), (*R5* + *R3*) / (*R5* + *R3* + 2·*C*). Hidden,
sub-threshold and off-screen junctions still count in the denominators; junctions touching no annotated site compete
nowhere. Groups sum their members' coverage and junction reads, and show each share's difference with the other groups.

**Share against the canonical alternative** (arc panel, `junctionAlternative`): a skipping junction against the mean of
the two inclusion junctions of the skipped block, *S* / (*S* + (*I*₁ + *I*₂)/2); an alternative donor or acceptor
against the canonical junction of the annotated exon on that side (rMATS A5SS / A3SS).

**Long-range junctions** (one end outside the gene, not annotated in another of its transcripts) are drawn from their
read count (*Min reads*) whatever *Min %*.

**RNA fusion junctions** come from the structural evidence of RNA tracks (§9: split reads and clipped reads placed by
realignment; pairs and CIGAR `D`/`I` are left out on RNA). Each end is drawn on the nearest exon boundary **within
5 bp**. A junction is kept when it spans **≥ 1 kb** and at least one end lies on an exon boundary; a deletion-type one
only when it joins two genes (inside one gene it is splicing). The 1 kb floor removes fold-back and template-switch
chimeras, which read as tiny inversions. Clipped bases are searched up to 1 Mb on each side of the window with a FASTA,
500 kb in all through the web APIs.

### 4.2 Intron retention

For a retained intron (IR pill), per sample:

- reads unspliced through the 5′ (donor) and 3′ (acceptor) boundary, *R5* and *R3*, and the canonical junction *C*;
- the share shown on the pill (§4.1) and the rMATS-style ratio (*R5* + *R3*) / (*R5* + *R3* + 2·*C*);
- an **IRFinder-style IR ratio** [16]: *d* / (*d* + *C*), with *d* the median depth inside the intron (the part in
  view, sampled at up to 201 evenly spaced positions). IRFinder takes the spliced reads at either end of the intron;
  here, the canonical junction's. Depth inside the intron, not only at its ends, so a pre-mRNA or a cryptic exon inside
  it raises the ratio too.

HGVS inserts the intron between its two exon bases (`r.<a>_<b>ins<a>+1_<b>-1`). Inside the CDS, the frame and the first
stop codon are read from the intron's bases, from the codon it splits; a stop with junctions downstream makes NMD likely.

### 4.3 NMD rules

`spliceModel.ts`, `nmdVerdict` (the splicing cartoon), after the EJC model [6]:

| Condition (first that applies) | Verdict |
|---|---|
| The annotated start codon is not in the transcript | start lost |
| No stop codon before the poly(A) | no stop (non-stop decay) |
| The canonical stop codon is reached | no PTC |
| The PTC lies in the last exon | escapes NMD |
| The PTC lies ≤ **55 nt** upstream of the last exon–exon junction | escapes NMD |
| The PTC lies < **150 nt** after the start codon | probable escape by re-initiation (uncertain) |
| Otherwise | degraded by NMD |

A 3′ UTR left longer than **1,000 nt** is noted as a possible trigger of EJC-independent NMD.

### 4.4 Exon inclusion and depth-based exon usage

**ψ of an exon** (exon popover, `exonPsi`): *I* = mean of the reads of the junctions ending at the exon's 5′ and 3′
boundaries (inclusion), *E* = reads of junctions spanning the whole exon (skipping); ψ = *I* / (*I* + *E*).

**Depth-based usage** (DEXSeq-style relative usage [4]):

- intervals: the coding part of each exon of a coding model (UTR parts left out: 5′ drop-off and 3′ poly(A) bias),
  the whole exon when it holds no coding base or the model is non-coding;
- usage of exon *i* = median depth of *i* / median of the median depths of the **other coding exons** (all other exons
  when fewer than two coding ones remain); its standard deviation by the delta method, usage·√(1/*n*_exon + 1/*n*_ref)
  from the read counts;
- controls = every other loaded sample; *vs controls* = usage / median usage of the controls;
- **robust z** = (usage − median) / max(1.4826·MAD, 0.02), only with **≥ 5 controls**; the 0.02 floor (2 points of
  usage) keeps a flat cohort from inflating z;
- exon depths are read from up to 100,000 reads per exon (scaled back when sampled), stranded as described in §2.

**GTEx tracks** [7] draw each exon's usage from the tissue's median junction reads alone, the only data GTEx exposes: at
each boundary, junctions ending there / (those + junctions spanning that boundary), averaged over the two boundaries.
Exons of other transcripts implied by the tissue's junctions (an acceptor followed by a donor within 3 kb) are added.

## 5. Variant sites and quality cells

`src/standalone/alleles.ts`, `src/standalone/collapse.ts` (`callSites`), `src/components/sashimi/siteQuality.ts`

### 5.1 Calling rule

The variants scan counts every read of the window (views ≤ 3 Mb) into per-position allele counts in a background
worker. A site needs **≥ 3 supporting reads** and an alternate fraction **≥ *min VAF*** (default 10 %; on long reads,
the larger of *min VAF* and *long reads*, default 20 %).

- **Every alternate base counts, whatever its base quality.** Bases under Q20 are flagged for the BQ cell instead of
  being dropped: leaving them out of the allele count while their reads stay in the depth under-reads every fraction (by
  a third on ONT data, where homozygous sites then read about 0.67 and heterozygous ones about 0.33).
- SNV denominator: reads with an aligned base at the position.
- **Indel denominator**: every read spanning the position, as a pileup counts it (`samtools depth -J`): the reads with a
  base there plus the reads deleting it, whatever their deletion. Insertions and deletions are keyed by position and
  length.
- The reads track calls its own sites from the reads it decoded with the same rule (`callSites`). The in-page phasing
  still uses only alleles of base quality ≥ 20 (§7).

### 5.2 Quality cells

Each site carries, for its alternate reads and for the other reads over the position, the share on the + strand, with
MAPQ < 20, with the call within 10 bases of an alignment end, and (SNVs) with base quality < 20. The counters are
unbounded (a packed counter past 8,191 continues in an overflow table), so deep amplicons are counted exactly.

| Cell | Measures | Amber (check) | Red (likely artefact) | Grey |
|---|---|---|---|---|
| **BQ** (SNV) | share of alternate bases with BQ < 20 | ≥ 25 % | ≥ 50 % | — |
| **HP** (indel, replaces BQ) | reference homopolymer length at the site | ≥ 6; long reads ≥ 4 | ≥ 10; long reads ≥ 7 | — |
| **MQ** | share of alternate reads with MAPQ < 20 (*a*) vs other reads (*o*) | *a* ≥ 20 % and *a* − *o* ≥ 15 points, or *o* ≥ 50 % | *a* ≥ 50 % and *a* − *o* ≥ 30 points | — |
| **SB** | strand split, alternate vs other reads, Fisher exact, two-sided (as GATK FisherStrand [37]) | p < 0.001 | p < 1e-4 | < 5 alternate or < 5 other reads |
| **END** | calls within 10 bases of an alignment end, alternate vs other reads, Fisher exact, one-sided | p < 0.001 | p < 1e-4 | < 5 alternate or < 5 other reads |

The cut-offs are strict because a window holds hundreds of sites (GATK's SNV hard filter FS > 60 sits near p = 1e-6).
Testing against the position's other reads keeps an amplicon library whose reads all run one way, or a region where every
read maps poorly, from being called an artefact. Sites too close for four cells get one cell, the worst colour. These
checks help judge a call; they are not a variant caller's filter.

### 5.3 Allelic balance and common-SNP balance

- **RNA site panel**: two-sided exact binomial test of alternate vs reference reads against 50 %. p ≥ 0.01: *balanced*;
  otherwise alternate share ≥ 0.9: near-monoallelic (alternate); ≤ 0.1: alternate barely expressed; else *tilted*. In
  RNA the fraction follows expression as well as genotype (NMD, allelic silencing, imprinting).
- **DNA track label** (common SNPs on, ≥ 3 common SNVs called): heterozygous = 0.2 ≤ VAF ≤ 0.8, homozygous = VAF > 0.8.
  *allele imbalance?* when the median |VAF − 0.5| of the heterozygous ones exceeds **0.15** over **≥ 5** of them; *no
  heterozygous SNP (LOH / UPD?)* when none is heterozygous and **≥ 8** are homozygous. A lead to check, not a
  copy-number call.

## 6. Haplotypes from haplotags

`src/components/sashimi/sitePhase.ts`, `src/standalone/alleles.ts` (`siteHaplotypesOf`), `SashimiViewer.tsx`

Phasing tools write `HP` (haplotype) and `PS` (phase set) on the reads they place: WhatsHap and LongPhase `haplotag`
[24, 25, 28], HiPhase [29], DRAGEN. For each site, the scan counts per (phase set, haplotype) the tagged reads whose
alignment spans the site and those carrying its allele (any base quality). **The site's phase set is the one with the
most HP 1 + HP 2 reads over it** (ties: the lower PS); reads of other sets are pooled apart. With *f*₁, *f*₂ the shares of
HP 1 and HP 2 reads carrying the allele:

| Call | Rule |
|---|---|
| not judged (grey dot) | fewer than 3 reads on either haplotype |
| both lanes (homozygous) | *f*₁ ≥ 0.7 and *f*₂ ≥ 0.7 |
| H1 | *f*₁ ≥ 0.7 and *f*₂ ≤ 0.2 |
| H2 | *f*₂ ≥ 0.7 and *f*₁ ≤ 0.2 |
| neither (**?**) | *f*₁ ≤ 0.2 and *f*₂ ≤ 0.2: the allele sits in untagged reads or another set's |
| to check (dashed red circle) | anything else: in part of one haplotype (mosaic, artefact, partial homopolymer indel) or on both without being homozygous (mis-phasing, collapsed duplication, third haplotype) |

A site whose quality cells include a red one is drawn as a grey dot, not phased. The panel adds a two-sided Fisher
test of HP 1 against HP 2. 70 % rather than the 80 % of the haplotype consensus (§7.3) keeps the indels long reads carry
only in part (homopolymers) on their haplotype. **HP numbers are the tool's labels**, arbitrary within each set: H1 of
one set is not known to be H1 of the next, and the numbering need not follow the `0|1` order of a VCF from another
phasing run.

## 7. Read-based phasing and consensus groups

`src/standalone/phasing.ts`, `src/standalone/haplotypes.ts`, `src/standalone/collapse.ts`, `src/standalone/readsWindow.ts`

### 7.1 In-page phasing

In the spirit of WhatsHap [24, 25] and HapCUT2 [26], in a greedy form that runs in milliseconds:

1. **Sites**: the window's called sites (§5.1); those at **25–75 %** alternate are phased, those above 75 % are
   homozygous (both haplotypes), those under 25 % are reported unphased (*low*).
2. **Fragments**: a read and its mate. A fragment's allele at a site is unknown when the base has quality < 20, is a
   third allele, or the two mates disagree.
3. **Links**: for each pair of heterozygous sites covered by the same fragments (each site linked to the next 32 it
   shares a fragment with), the fragments in phase (ref–ref, alt–alt) and in opposite phase. A link is **trusted with
   ≥ 2 fragments and ≤ 20 % of them disagreeing**.
4. **Blocks**: sites in order; a site joins the current block when its trusted links to the block's sites agree on its
   phase (a vote weighted by fragments, ≤ 20 % against); else a new block starts, with the reason kept (no link, or
   contradicting links). **A single site that does not fit while the next one does is an outlier**: it is left
   unphased and the block goes on. Single-site blocks are unphased. Blocks act as phase sets.
5. **Fragments to haplotypes**: per block, the haplotype the fragment matches at more sites. In collapsed mode a
   fragment goes to the block where it is most clearly on one side (largest margin); in *group: phased here*, to the
   block where it covers the most phased sites (blocks of ≥ 2 sites); ties stay unassigned.

Short-read pairs link sites within an insert; long reads phase whole windows; RNA-seq phases with the fewer
heterozygous sites its exons carry, spliced reads linking exons.

### 7.2 Haplotag source

With *phase: HP tags* (the default whenever the window holds tagged reads), each fragment takes the tag of the mate that
has one, and groups are (phase set, haplotype); untagged fragments are counted and left out.

### 7.3 Haplotype consensus and checks

For each (phase set, haplotype): **covered** where ≥ 3 of its reads align or delete bases; a variant is in the consensus
when **≥ 50 %** of the haplotype's reads carry it (at least 2 reads), inside covered stretches, and when it is also a site
of the window over all reads (two sequencing errors among the 3–4 reads of a haplotype at a window edge would otherwise
make a variant). **Split check**: at each heterozygous site (25–75 % over all reads) with ≥ 3 reads on each haplotype,
the two haplotypes should split it, **≥ 80 %** of one carrying the alternate and **≤ 20 %** of the other; sites that do
not are listed as *not split*. Reads whose allele contradicts their haplotype at a split site are counted.

### 7.4 Junctions per haplotype (RNA)

For every junction carried by **≥ 3 fragments** of a phase set, each haplotype's fragments are counted as *using* it or
*going another way* at one of its ends (another junction from its donor or to its acceptor, an exon skip included, or
aligned bases across the exon–intron boundary). With each haplotype holding ≥ 3 informative fragments, the two shares
are compared by a two-sided Fisher exact test; a junction with **p < 0.001 and ≥ 20 points** between the shares is
flagged as used differently by the haplotypes: the splice change of a variant in *cis*. Long-read junctions are merged
as in §3.3 first.

### 7.5 Consensus groups (*haplotypes: any*)

One row per local haplotype × splice pattern. A fragment's pattern is its junctions and its allele at each called site
it covers (unknown where its mates disagree).

- **Unspliced short reads (DNA)**: exact patterns. Seeds are the patterns of ≥ *Min reads* fragments, most specific first,
  that fit no seed already taken; every pattern is then placed on the one compatible seed (same junctions, no block
  across the seed's junctions, the same allele at every site it calls), ambiguous when several fit, *minor* when none;
  seeds left under *Min reads* are dropped and placement is redone.
- **Long reads, and spliced short reads (RNA-seq)**: tolerant clustering, since no two long reads share a pattern and
  RNA fragments see different pieces of one transcript. Fragments are taken left to right; each joins the group whose
  consensus it agrees with at their shared sites with **≤ 20 % disagreeing** and whose splicing it never contradicts (a
  junction of one crossed by aligned bases of the other; two different overlapping junctions). Shared junctions count
  as agreement; a fragment with calls meets groups through its sites only. A fragment agreeing with two groups goes to
  the one it agrees with by **≥ 2** more sites or junctions, else is ambiguous. Groups are then merged when they agree
  where they overlap with **≥ 2 shared sites** or **≥ 2 fragments linking them** (at a site of each), unless the group
  could join two groups that differ from each other (it then stays apart). Every fragment is finally placed again
  against the final groups. Sites above **90 %** alternate are shown but not used. When, after 300 patterns, more than
  30 % of them seed their own group (a wrong reference, a junk window), the exact patterns are used instead.

## 8. Pair phasing (*Phase with*)

`src/standalone/pairPhase.ts`, `src/components/sashimi/PairPhasePanel.tsx`

The pairwise form of read-backed phasing, for a possible compound heterozygote.

**Sites.** Typed as `chr:pos` (1-based) with an optional allele (`A>G`, VCF-style `TTTG>T` anchored on the base before,
`+2` / `-3`, gnomAD `1-55051215-G-GA`). The variant scan's site is used when there is one; without an allele, the most
supported allele at the position. Reads are read in a small window around each site (the base, or the deleted bases,
±1 base), never the stretch between them, **up to 20,000 reads per site**; past that, whole fragments are kept by name
hash at a power-of-two rate, so the fragments kept at the deeper site are kept at the other one too.

**Fragments** are joined by read name across the two windows: one read holding both sites (long reads), an RNA read
spliced across the introns between them (phASER's principle [27]), two mates, or the parts of a split read. Each site's
allele is taken from the copy of the read decoded at that site. An allele is unknown under base quality 20, for another
allele, or when two reads of the fragment disagree.

**Table and verdict.** Fragments by allele at A × allele at B. Two heterozygous sites: alt–alt and ref–ref are *cis*,
alt–ref and ref–alt *trans*. A site **under 25 % VAF** counts only the fragments carrying its allele (its reference reads
come from cells without it too). A site **over 75 %** is homozygous: no verdict needed.

| Informative fragments *n* | Against | Verdict |
|---|---|---|
| < 3 | — | undetermined |
| ≥ 3 | trans ≤ 20 % of *n* | in cis |
| ≥ 3 | cis ≤ 20 % of *n* | in trans |
| ≥ 3 | > 20 % against either | conflicting |

**LOD** = log10 L(cis) / L(trans) with an error rate ε = 1 % per allele read, i.e. (cis − trans) ×
log10[((1 − ε)² + ε²) / (2ε(1 − ε))] ≈ 1.69 per fragment (3:0 gives 5.1, 30:0 gives 50.8). It grows with depth
whatever the share against, so it is shown and does not decide (likelihood after HapCUT2 [26]).

**Haplotags.** Each site is placed in its phase set by the rule of §6; same HP is *cis*, different HP *trans*; two
different phase sets, or a site not placed on H1/H2, leave the tags undecided.

**Order of evidence.**

1. The reads holding both sites, when they give cis or trans. Haplotags that agree are named; haplotags that disagree
   make the pair **conflicting** (a switch error of the tool, or reads joined that are not one molecule).
2. The haplotags, when the reads cannot decide: fewer than 3 informative fragments, or *conflicting* with fewer than
   **10** informative fragments.
3. Without haplotags deciding, and no homozygous site: the reads of the whole stretch between the sites (±1 kb, up to
   3 Mb, up to 40,000 reads) are phased by §7.1 (sites ≥ 3 reads and ≥ 5 % VAF, 20 % on long reads; long reads on
   their SNVs alone, since their homopolymer indel errors broke the blocks every few kb). The pair is cis or trans when
   one block holds both sites; otherwise the panel says why (a site not phased, two blocks no read links).
4. Otherwise the reads' verdict stands (undetermined, conflicting).

**A variant and a splice junction** (allele-specific splicing). A junction is typed as `chr:first-last` (the intron's
first and last bases, 1-based) or picked. A fragment *uses* the junction when one of its reads is spliced exactly from
its first to its last intronic base (long reads' junctions merged as in §3.3; short reads never); it *passes it
otherwise* when a read spans ≥ 5 bases on each side of one of its splice sites without using it (unspliced, spliced
elsewhere, skipping). The table crosses allele × use; the **usage on each copy** is uses / (uses + other) per allele, and
a two-sided Fisher exact test measures the association. From the fragments using the junction (*u* ≥ 3 needed):

- **on the variant's copy** when ≤ 20 % of them carry the reference allele;
- **on the other copy** when ≤ 20 % of them carry the variant, provided ≥ 3 fragments carry the variant over the
  junction's splice sites (else undetermined);
- **on both copies** otherwise; when p < 0.01, the copy that uses it more is named.

A homozygous variant cannot tell the copies apart. Counts are of molecules that reached the library: a transcript
degraded by NMD is under-counted on its copy.

## 9. Structural evidence

`src/standalone/alignments.ts` (`structuralEvidence`, `placeClip`, `rescueClipEnds`), `src/standalone/svmerge.ts`,
`src/standalone/arcSupport.ts`, `src/components/SashimiViewer.tsx`

Structural hints are drawn on every DNA track as **evidence, never as calls**. The coverage scan keeps only the records
that can carry evidence (an SA tag, a clip of ≥ 8 bases, a `D` of ≥ 50 bp, an insertion of ≥ 50 bp, discordant mates)
and an exact histogram of the proper pairs' |TLEN|, whose median is the window's insert size.

### 9.1 Kinds of evidence

- **Deletions inside reads**: a CIGAR `D` of **≥ 50 bp**.
- **Split reads, read as chains.** Every part of a read (the primary record and the parts its SA tag lists, primary
  preferred as the chain's record) is placed on the read by its clips; parts are ordered along the read, and each pair
  of adjacent parts is one breakpoint, each read counted once. Where the read continues: on the other strand,
  *inversion*; ≥ 50 bp further on the reference, *deletion-type*; ≥ 50 bp back, *duplication-type*; on another
  chromosome, a `→ chr` pill; an unaligned stretch of ≥ 50 read bases between two parts adjacent on the reference, an
  `ins` pill with its length. Breakpoints are rounded to 5 bp.
- **Soft-clip clusters**: ≥ 3 reads without SA tag clipped by **≥ 20 bases** at one position (soft or hard clips).
- **Clipped reads placed by realignment.** The clipped consensus of a cluster (majority base where ≥ 60 % agree, N
  otherwise, while ≥ 2 reads cover the position) is searched on both strands of the window's reference: its 20 bases
  next to the breakpoint must occur **exactly once**, and up to 60 bases must then agree with **at most one mismatch per
  20**. A placed cluster becomes the breakpoint a split read would give, typed the same way. Hard-clipped records without
  SA tag count in the cluster at their clip position.
- **CIGAR insertions of ≥ 50 bp** (some realigners, such as ABRA2 [40], write a tandem duplication as an insertion of the
  copy's bases). Insertions are grouped (position within 5 bp, length within 10 %). When the first and last 24 inserted
  bases are the reference just before (or just after) the insertion point, **within 60 bp**, the reads count on a
  duplication arc over the copied stretch; otherwise they count as an insertion pill. The copy's far end is checked
  where the reference reaches it (the reference is read up to 50 kb further around the window); otherwise the
  insertion's length places it.
- **Discordant pairs** (§9.3).

### 9.2 One arc per event

Arcs of one family (deletion: CIGAR `D`, deletion-type split reads and placed clips; duplication; inversion) are
merged when **every start and every end** of the group lies within a tolerance of each other (complete linkage):
**tolerance = 5 % of the event's length, at least 20 and at most 100 bp**, the smaller event's tolerance applying. The
ceiling follows long-read SV callers' ONT settings (cuteSV recommends a 100 bp cluster bias for ONT deletions,
`--max_cluster_bias_DEL` [32]; Sniffles2 merges insertions and deletions within 150 bp by default outside tandem
repeats, `--cluster-merge-pos` [33]); the floor and the relative term keep two 60 bp deletions 40 bp apart as two
arcs. The event's breakpoints are the read-weighted medians of what it merged; its count is the number of **distinct
reads** (a read showing both junctions of an inversion counts once). The same tolerance matches an event across
samples and pools the members of a DNA group.

### 9.3 Discordant pairs

A pair is counted once, from its leftmost mate, or from the other mate when the leftmost one is not among the records;
from primary records only. Its extent comes from the **mates' positions, not TLEN** (an aligner may give the pairs of a
duplication a short TLEN and the proper-pair flag). Orientation is judged before the proper-pair flag. With *m* the
median insert:

| Class | Rule |
|---|---|
| inversion-type | both mates on one strand |
| deletion- or duplication-type, **pair across one junction** | opposite strands, the forward mate's leading clip and the reverse mate's trailing clip both ≥ 10 bases, and ≥ 50 bp between the reverse mate's end and the forward mate's start: deletion-type when the junction goes forward, duplication-type when it goes back |
| duplication-type (facing away, ← →) | the reverse mate **ends** before the forward mate starts (20 bp of slack), and the mates are > max(300 bp, 2·*m*) apart |
| deletion-type (facing each other, → ←) | span (minus the deleted bases the reads already carry) > max(5·*m*, 1 kb) |

A pair a mate of which carries a CIGAR `D` of ≥ 50 bp is not counted (the event is drawn from that alignment). A pair a
mate of which carries an insertion of ≥ 50 bp is kept only when it reads as duplication-type. Pairs are binned by 500 bp;
bins of one class within 500 bp of the strongest at both ends join it. The arc is drawn where the reads place the
breakpoints: for a duplication the **outermost** reads (start of the leftmost left mate, end of the rightmost right
mate), for a deletion the **innermost**, for an inversion the median; the true breakpoint lies within an insert size of
these ends. Discordant arcs are matched across samples with the same class and both ends within **1 kb**.

### 9.4 Rescue at known breakpoints

As a two-pass aligner uses the junctions of its first pass, a read **soft-clipped by ≥ 8 bases** (≥ 12 when several
breakpoints have an end at its position) at an end of an arc of the sample (chain, CIGAR deletion or placed cluster)
joins that arc when its clipped bases match the reference at the arc's other end in the arc's orientation (no mismatch
under 20 bases, one per 20 above) and fit exactly one breakpoint. The clip must sit within **2 bp** of the arc's
(rounded) end, or **at** the end for a CIGAR deletion. A hard clip of ≥ 20 bases without SA tag is attached by position
when exactly one breakpoint end can hold it. Rescued reads never make an arc on their own. The breakpoints of the other
DNA samples shown are targets too (their reads read ±400 bp around each end); a sample's reads rescued at a breakpoint
it does not carry are listed in the arc's panel (*no arc*), not drawn.

### 9.5 Support threshold and supporting reads

An arc is drawn from **3** reads by default, and **20 on a deep track**: a DNA track whose median depth over the covered
bases of the model's coding exons (≥ 100 such bases, else the covered bases of the window) is **above 200×**, measured
once per gene on its first complete coverage. At a few hundred ×, PCR chimeras, adapter read-through and capture edges
give a handful of odd reads at many places. A number typed by the user applies to every track.

*Supporting reads* (`arcSupport.ts`), tested before any downsampling: a junction's `N` exactly (within 6 bp on long
reads); a `D` ≥ 50 bp at a deletion, or a break at either end of the arc (an alignment end followed by a clip of ≥ 8
bases, an insertion ≥ 50 bp, an SA part's end) within the merge tolerance; for pairs, the arc's class, a mate within
1.5 kb of each end, the class's distance (1 kb; 300 bp facing away) and half the arc spanned. 300 kept, plus mates.

## 10. CpG methylation

`src/standalone/methylation.ts`, `src/standalone/localSource.ts`, `src/components/SashimiViewer.tsx`

### 10.1 Calls

Counted the way `modkit pileup --cpg --combine-strands` (traditional preset) builds a CpG table [10]; pb-CpG-tools does
the same per haplotype for PacBio [11].

- **Sites from the reference**: every CpG of the reference; a + strand read counts at the C, a − strand read at the G
  (the C of the other strand), both summed on the + strand C, and only when the read carries the canonical base there.
- **Tags**: `MM` / `ML` (legacy `Mm` / `Ml`) as the SAM specification describes [9]; ML byte *N* ↦ (*N* + 0.5)/256; a `?`
  entry says nothing of the bases it skips, a `.` entry (or none) says they are unmodified. A record is skipped when its
  sequence no longer matches the tags: hard-clipped without an `MN` tag, or `MN` ≠ SEQ length.
- **5hmC set aside**: its probability is removed and the 5mC and unmodified probabilities renormalised; the call is the
  likelier state, its probability the call's confidence.
- Primary alignments only (no secondary, supplementary, duplicate or QC-failed record); the MAPQ filter is not applied.

### 10.2 Confidence filter

Calls below the **10th percentile** of the confidences are filtered (modkit's default). modkit estimates one threshold
per file from a sample of reads; here it is estimated **per 100 kb block** of the chromosome (a fixed grid from its
start) from every call at the block's CpGs, so a CpG's value depends only on its reads and its block, never on the zoom,
what was viewed before, the layer being toggled or the file format. Confidences are kept in **16 bins** over [0.5, 1],
bin *k* starting at 1 − 2^−(1 + *k*/2) (0.5, 0.65, 0.75, … 0.992, 0.9945, 0.9961, 0.9972: near 1, about one ML step).
The threshold is the bin where the cumulative count reaches 10 %; the bins below it are dropped (≤ 10 % of the calls). A window is counted over whole blocks: the view plus half of it on each side (within 2 Mb in all), rounded out
to 100 kb blocks. Each read's ticks in the reads track are compared with its block's threshold by the same bins.

### 10.3 Ribbon, phase sets and smoothing

Counts are kept per CpG, per haplotag (untagged or HP > 2, HP 1, HP 2) and per confidence bin, and the tagged calls again
per phase set. **Each CpG takes the phase set holding most of its passing tagged calls** (ties: the lower PS); the other
sets' tagged calls join the untagged ones there. A pixel never pools two sets: its HP 1 and HP 2 come from the set
holding most of its tagged calls.

- **Smoothing (one column per pixel)**: the pooled fraction of methylated calls over the CpGs under the pixel (calls are
  summed, not per-CpG fractions averaged). When fewer than **6** CpGs lie under it, the nearest CpGs on either side are
  added up to 6, **within 1 kb** of the pixel. A CpG desert wider than that stays a gap. A column with **< 3 calls** is
  drawn pale.
- **Split into HP 1 / HP 2** where each haplotype holds **≥ 3 calls per CpG** of the column and the less covered one holds
  **≥ 25 %** of the tagged calls; elsewhere one band of all reads. Untagged reads do not count against the split (a tool
  tags reads by their alleles, so the tagged reads are a fair sample of each copy). Split runs narrower than 4 pixels are
  joined.
- **Difference** HP 1 − HP 2 per split column, in its phase set's hue (its sign has no meaning across sets).
- **Allele-specific methylation (ASM)**: a run of **5 consecutive CpGs** covered on both haplotypes, within one phase
  set, whose pooled fractions differ by **≥ 50 points** with **≥ 10 calls** on each side; overlapping runs of the same
  sign are merged and the difference recomputed over the stretch's own CpGs.
- **CpG islands** (Gardiner-Garden & Frommer [8]): 200 bp windows with GC ≥ 50 % and observed/expected CpG ≥ 0.6, merged
  where they overlap, found over the whole counted stretch so an island is never cut by the view. *CpG islands only*
  keeps the CpGs of the islands and confines the smoothing to their bounds.

### 10.4 Island differences between samples

With two or more long-read DNA samples, per island: over the CpGs covered by **≥ 5 calls in both samples** (paired), the
mean of the per-CpG 5mC fractions (mean β, all reads); an island needs **≥ 3** such CpGs. The primary sample's tag is its
mean minus the mean of the others'; each other sample's, its own mean minus the primary's. Grey under **10 points**, bold
from **20**.

## 11. Repeat inspector

`src/standalone/repeatScan.ts`, `src/standalone/repeatConsensus.ts`, `src/standalone/strCatalog.ts`,
`src/standalone/strFlanks.ts`, `src/components/sashimi/RepeatInspector.tsx`

**Loci.** The catalogue holds the **77** pathogenic tandem-repeat loci of STRchive [17] with motifs of ≤ 12 bp. The tract
is found in the reference around the catalogue's coordinates: from the seed giving the longest tract, motif units in
phase, extended while no two units in a row differ from the motif by more than one base, ending on exact units, at least
3 units. A stretch selected by hand takes the shortest motif (1–6 bp) whose rotations cover most of it (≥ 50 %; a longer
motif must cover 5 points more). Interruptions are the catalogue's, else the reference units one base off the motif.

**Reads.** Primary records over the locus (± 490 bp), all MAPQ, up to 20,000 (whole fragments, every *k*-th past that).
Each read is rebuilt over the locus from its soft clips, aligned bases with their mismatches, and insertions, then
searched as stored and reverse-complemented.

**Anchors.** Two 40 bp locus-specific stretches, one on each side, place the tract's ends. Published flanks are used where
a tool gives them (STRique's FMR1 and C9orf72 targets [21]), placed in the reference next to the tract; otherwise the
reference's, as STRipy [22] and ExpansionHunter [41] do. Each anchor is the 40 bp window **nearest the tract, within
150 bp, that differs from the repeat by ≥ 35 %** of its bases (semi-global edit distance), else the most different one.
In a read, exact 10 bp seeds (every 2 bases) propose places and the whole anchor must align there with **≤ 6 edits
(15 %)**.

**Classes of reads.**

| Class | Rule | Use |
|---|---|---|
| spanning | both anchors, once each, 5′ before 3′, tract **≥ 80 %** motif, benign-motif or known-interruption units | sized: tract length / motif length, rounded |
| impure | both anchors, tract < 80 % | counted apart |
| anchor twice / out of order | an anchor found twice, or 3′ before 5′ (concatemer, fold-back) | counted apart |
| stopping inside (lower bound) | one anchor; the motif-rich stretch next to it, followed while ≥ 60 % of the last 10 units are motif units; ≥ 3 units | lower bound |
| all repeat | no anchor, ≥ 10 motif units covering ≥ 80 % of the read | lower bound |

Units are read in the motif's phase; one or two bases breaking the phase (a sequencing indel) are skipped to the next
unit.

**Sizes and alleles.** With **≥ 20** spanning reads, sizes come from them alone; below, the reads stopping inside are
added as lower bounds. Alleles: up to two modes of a Gaussian kernel density on log(size) (Silverman bandwidth
0.9·min(sd, IQR/1.34)·*n*^−0.2, at least 0.04; grid step 0.005), so that a peak's width grows with its size as stutter
does; a peak is kept when it reaches 15 % of the highest and the valley towards any higher peak drops under 70 % of it;
reads are split at the valleys. Each allele gives its mode, median, P5–P95 and read count; it is **broad** when it holds
≥ 10 reads and (P90 − P10) / median > 0.25 (somatic mosaicism or PCR/sequencing stutter, not decided). A note warns that
long alleles are under-represented when ≥ 10 % of the reads stopping inside go beyond the spanning reads' P95.

**Categories.** FMR1 by the ACMG technical standard [18]: normal ≤ 44, intermediate 45–54, premutation 55–200, full
mutation > 200. HTT by the ACMG/ASHG guidelines [19]: normal ≤ 26, intermediate 27–35, reduced penetrance 36–39, full
penetrance ≥ 40. Other loci: STRchive's benign, intermediate and pathogenic ranges. Loss of the 5′ AGG interruptions of
an FMR1 allele, which makes expansion on maternal transmission more likely [20], reads as a missing orange stripe in the
waterfall (after TRGT's waterfall plots [23]).

**Consensus.** Up to 300 reads (evenly by size). Flanks base by base: each read's 120 bp flank aligned on the
reference's (fixed at the tract), the majority base, gap or insertion while half the reads reach the position. Tract unit
by unit: each read's units aligned on those of the median-size read (centre-star), the majority at each column, so
interruptions keep their place. Support < 60 % is lowercase. Reads stopping inside make the tract only when no listed
read spans it (a lower bound). Each allele's consensus uses its **P10–P90 core** (all its reads if the core has < 3).

The inspector reports evidence from the reads, not a diagnostic call.

## 12. Sampling, memory and speed

`src/standalone/localSource.ts`, `src/standalone/readsWindow.ts`, `src/components/SashimiViewer.tsx`,
`src/components/sashimi/datasource.ts`

| What | Budget |
|---|---|
| Coverage, junctions, boundary counts | exact, every read; view counted first, margins up to ~2 M reads; kept up to 8 Mb per sample, 256 MB for all samples (least recently used dropped) |
| Variant scan | every read; views ≤ 3 Mb; allele counts kept up to 4 Mb per sample, 384 MB for all |
| Methylation | views ≤ 2 Mb; kept up to 6 Mb per sample, 256 MB for all; read ticks on views ≤ 30 kb |
| Reads track | loaded for views ≤ the reads window (100 kb default; 250 kb, 500 kb, 1 Mb, 2 Mb); **2,500 reads drawn per track**, whole fragments kept (every *k*-th fragment past the cap) |
| Collapsed mode | up to 40,000 reads decoded, in the background worker |
| *group: phased here* | up to 10,000 reads phased, then 2,500 drawn |
| Supporting reads | every read tested; 300 kept, plus mates |
| *Phase with* | 20,000 reads per site; 40,000 for the stretch between sites |
| Repeat inspector | 20,000 reads over the locus |
| Exon-depth usage | 100,000 reads per exon |
| Record cache of the alignment libraries | one budget of 128 MB for all files, in the page and in the worker each; a record unused for 45 s is dropped |

**Sampling rules.** Where a cap applies, the decoded reads are thinned by name hash at a power-of-two rate (a read is
kept when the hash of its name is a multiple of 2^*k*, as `samtools view --subsample` does), so both mates of a pair
share their fate; the reads track then keeps every *k*-th fragment. Coverage, junctions, retention counts, variant counts
and methylation are never sampled.

**Speed.** Records are decoded in tiles of about 20,000 records sized from the index (16 kb floor, the BAI linear-index
interval; 250 kb ceiling), the thread handed back between tiles; a coverage slice yields every 8 ms. Variant scans,
methylation and collapsing run in a background worker. The scan gives an allele a table entry only at its third call
(the first two calls' flags fit in a 16-bit word), so the isolated errors of noisy long reads cost little memory. With
more than 400 reads carrying CpG ticks in view, the ticks are left out while the view moves and return 200 ms after.

## 13. Reference builds

`src/standalone/assembly.ts`, `src/standalone/liftover.ts`

A file's build is read from its header: an `@SQ AS` tag naming it, else the length of a primary chromosome, which
differs between GRCh38 and GRCh37 [15] (chr1: 248,956,422 vs 249,250,621 bp; T2T-CHM13 and NCBI36 are recognised and
reported as unannotated). Coordinates are lifted between GRCh38 and GRCh37 with Ensembl's assembly map
(`/map/human/…`) [14], which follows the GRCh37 contigs kept in GRCh38 to the base; only the chromosome and window
coordinates are sent, in 1 Mb pieces. A stretch with no counterpart lifts to nothing.

## 14. Validation

Measured on the 1.2.0 code during development, on public data and simulations; not re-run automatically (the
repository has no test suite). *Right*: agreeing with the truth set; *wrong*: a cis/trans call contradicting it.

**Real data: ONT GIAB HG002** (wf-human-variation, 60×, haplotagged BAM and phased VCF of the same pipeline [30]),
CFTR and HBB regions (500 kb each).

| Check | Result |
|---|---|
| Haplotype lanes (§6) vs the phased VCF, all phased heterozygous SNVs | 1,413 of 1,413 agree within every phase set (every pair in cis or trans as the VCF says); the HP numbering and the VCF's genotype order agree in some sets and are swapped in others |
| *Phase with*, pairs of phased heterozygous SNVs of one set, ≤ 40 kb apart (CFTR) | 294 of 301 right, none wrong; from the reads alone 232 right, none wrong (the rest had no long read spanning both) |
| *Phase with*, pairs ≤ 150 kb apart | CFTR 198 of 202 right, HBB 201 of 201, none wrong; most rest on the haplotags, which come from the same pipeline as the VCF, so this checks the panel's reading of the tags, not the tool's phasing |
| *Phase with* undetermined pairs | each holds a site at 28–32 % VAF that the tags do not place |
| *Phase with* through the stretch (§8, step 3), BAMs stripped of their haplotags, pairs ≤ 150 kb | CFTR 200 of 202 right (145 of them through the stretch), HBB 195 of 201 right, none wrong; the rest are breaks between blocks |
| Methylation filter (§10.2), CFTR 200 kb | 7.0 % of the calls filtered |
| Methylation phase sets (§10.3), CFTR 200 kb | 179 of 1,502 CpGs hold the reads of two sets; at each CpG a set's counts equal those of a file holding only that set's reads |
| Methylation split (§10.3), haplotagged ONT sample at *UBE3A*, 90 kb view | 52 % of the view drawn split; the rest is a gap between phase sets where no read is tagged |
| Variant scan time, 500 kb at 60× ONT | 5.0–5.3 s for a phased file, 4.1–4.6 s unphased |

**Simulations.**

| Check | Result |
|---|---|
| Haplotype lanes: two phase sets, a mosaic SNV, an allele only in untagged reads, short pairs and noisy long reads | every site on the expected lane; counts equal an independent pysam count |
| *Phase with*, simulated 2×150 pairs | every pair within an insert right, none wrong |
| *Phase with*, simulated RNA-seq (2×100, three exons over 25 kb) | every pair linked by spliced reads or mates right, including exons 15 kb apart, none wrong |
| Consensus groups (§7.5): 9-exon gene, one haplotype skipping exon 5 in cis with its alleles, an alternative acceptor on both; ONT-like cDNA (4 % errors, 10 % of junctions shifted 1–4 bp) | the five haplotype × isoform classes as five pure groups |
| Same gene, 2×100 pairs | the exon-skipping reads as one group of the alternate haplotype; pieces the reads cannot assign to one isoform as rows of their own or ambiguous rows; every group pure |
| Long-read junction merge (§3.3): acceptor 4 bp into an exon, used by 40 % of one haplotype's molecules | kept as its own arc (137 reads next to 594), its own consensus group, flagged as used by that haplotype only (a 1-in-4 rule would have merged it) |
| Quality cells (§5.2), clean synthetic data, 199 sites at 40× | none flagged |

**Structural evidence, cases behind the rules (§9).**

| Case | Observation |
|---|---|
| LDLR capture, a 4 kb tandem duplication | the arc falls at 11,229,087–11,233,051, at the feet of its supporting reads (bin edges would give 11,229,000–11,233,500); BWA gave the pairs a TLEN of about 110 and the proper-pair flag with mates 3.9 kb apart, hence the extent from positions |
| LDLR capture, reads written with a CIGAR insertion | 93 reads with a 3,964 bp insertion written at 40 places inside the repeat give one duplication arc next to the 55 pairs facing away; without the insertion-mate and crossing-pair rules, 33 pairs of an LDLR duplication read as mates facing each other far apart (a deletion) |
| VWF deletion | comparing mate starts instead of the reverse mate's end reads 103 pairs as a duplication over 288 deletion reads |
| Exons named in a duplication's tooltip | whole exons 13–14 of LDLR NM_000527.5 hold 295 coding bases, not a multiple of 3: their tandem copy shifts the frame |

## References

[1]–[31] are the user manual's former list, in its order; [32]–[43] are the tools the methods are compared with.

- [1] Katz Y, Wang ET, Airoldi EM, Burge CB. Analysis and design of RNA sequencing experiments for identifying isoform
  regulation. *Nat Methods* 2010;7:1009–1015. doi:10.1038/nmeth.1528 (the Sashimi plot, MISO `sashimi_plot`).
- [2] Garrido-Martín D, Palumbo E, Guigó R, Breschi A. ggsashimi: Sashimi plot revised for browser- and
  annotation-independent splicing visualization. *PLoS Comput Biol* 2018;14:e1006360.
- [3] Morales J et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature*
  2022;604:310–315.
- [4] Anders S, Reyes A, Huber W. Detecting differential usage of exons from RNA-seq data (DEXSeq). *Genome Res*
  2012;22:2008–2017.
- [5] Shen S et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq
  data. *PNAS* 2014;111:E5593–E5601.
- [6] Kurosaki T, Popp MW, Maquat LE. Quality and quantity control of gene expression by nonsense-mediated mRNA decay.
  *Nat Rev Mol Cell Biol* 2019;20:406–420 (the 50–55 nt rule of the NMD verdict).
- [7] GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*
  2020;369:1318–1330.
- [8] Gardiner-Garden M, Frommer M. CpG islands in vertebrate genomes. *J Mol Biol* 1987;196:261–282.
- [9] The SAM/BAM Format Specification Working Group. *Sequence Alignment/Map Optional Fields Specification*, §1.7
  Base modifications (MM, ML, MN). https://samtools.github.io/hts-specs/SAMtags.pdf
- [10] Oxford Nanopore Technologies. *modkit* documentation: `pileup`, `--cpg`, `--combine-strands`, filter threshold
  (10th percentile), 5hmC handling. https://nanoporetech.github.io/modkit/
- [11] PacBio. *pb-CpG-tools*: CpG methylation probabilities from HiFi reads, per haplotype.
  https://github.com/PacificBiosciences/pb-CpG-tools
- [12] Levin JZ et al. Comprehensive comparative analysis of strand-specific RNA sequencing methods. *Nat Methods*
  2010;7:709–715. doi:10.1038/nmeth.1491.
- [13] Mourão K, Schurch NJ, Lucoszek R et al. Detection and mitigation of spurious antisense expression with RoSA.
  *F1000Research* 2019;8:819.
- [14] Yates A et al. The Ensembl REST API: Ensembl data for any language. *Bioinformatics* 2015;31:143–145.
- [15] Schneider VA et al. Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring
  quality of the reference assembly. *Genome Res* 2017;27:849–864.
- [16] Middleton R et al. IRFinder: assessing the impact of intron retention on mammalian gene expression. *Genome
  Biol* 2017;18:51. doi:10.1186/s13059-017-1184-4.
- [17] Hiatt L, Weisburd B, Dolzhenko E et al. STRchive: a dynamic resource detailing population-level and
  locus-specific insights at tandem repeat disease loci. *Genome Med* 2025;17:29. doi:10.1186/s13073-025-01454-4
  (loci, motifs and ranges; data under the MIT License).
- [18] Spector E, Behlmann A, Kronquist K et al. Laboratory testing for fragile X, 2021 revision: a technical standard
  of the American College of Medical Genetics and Genomics (ACMG). *Genet Med* 2021;23:799–812.
  doi:10.1038/s41436-021-01115-y.
- [19] ACMG/ASHG Huntington Disease Genetic Testing Working Group. Laboratory guidelines for Huntington disease
  genetic testing. *Am J Hum Genet* 1998;62:1243–1247.
- [20] Nolin SL et al. Fragile X full mutation expansions are inhibited by one or more AGG interruptions in
  premutation carriers. *Genet Med* 2015;17:358–364.
- [21] Giesselmann P et al. Analysis of short tandem repeat expansions and their methylation state with nanopore
  sequencing. *Nat Biotechnol* 2019;37:1478–1481 (STRique; its published FMR1 and C9orf72 flanks; MIT License).
- [22] Halman A, Dolzhenko E, Oshlack A. STRipy: a graphical application for enhanced genotyping of pathogenic short
  tandem repeats in sequencing data. *Hum Mutat* 2022;43:859–868.
- [23] Dolzhenko E et al. Characterization and visualization of tandem repeats at genome scale. *Nat Biotechnol* 2024
  (TRGT and its read waterfall plots).
- [24] Martin M, Patterson M, Garg S et al. WhatsHap: fast and accurate read-based phasing. *bioRxiv* 2016:085050.
  doi:10.1101/085050 (its `haplotag` command writes HP / PS).
- [25] Patterson M, Marschall T, Pisanti N et al. WhatsHap: weighted haplotype assembly for future-generation
  sequencing reads. *J Comput Biol* 2015;22:498–509. doi:10.1089/cmb.2014.0157.
- [26] Edge P, Bafna V, Bansal V. HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies.
  *Genome Res* 2017;27:801–812. doi:10.1101/gr.213462.116 (the likelihood behind the pairwise LOD).
- [27] Castel SE, Mohammadi P, Chung WK, Shen Y, Lappalainen T. Rare variant phasing and haplotypic expression from
  RNA sequencing with phASER. *Nat Commun* 2016;7:12817. doi:10.1038/ncomms12817.
- [28] Lin JH, Chen LC, Yu SC, Huang YT. LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and
  large variants. *Bioinformatics* 2022;38(7):1816–1822.
- [29] Holt JM et al. HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing.
  *Bioinformatics* 2024;40(2):btae042. doi:10.1093/bioinformatics/btae042.
- [30] Oxford Nanopore Technologies. GIAB 2023.05 open data release: HG002 haplotagged BAM and phased VCF of
  wf-human-variation, `s3://ont-open-data/giab_2023.05`.
- [31] GMOD JavaScript libraries: bam-js (https://github.com/GMOD/bam-js), cram-js (https://github.com/GMOD/cram-js),
  indexedfasta-js (https://github.com/GMOD/indexedfasta-js).
- [32] Jiang T, Liu Y, Jiang Y et al. Long-read-based human genomic structural variation detection with cuteSV.
  *Genome Biol* 2020;21:189. doi:10.1186/s13059-020-02107-y.
- [33] Smolka M, Paulin LF, Grochowski CM et al. Detection of mosaic and population-level structural variants with
  Sniffles2. *Nat Biotechnol* 2024;42:1571–1580. doi:10.1038/s41587-023-02024-y.
- [34] Prjibelski AD, Mikheenko A, Joglekar A et al. Accurate isoform discovery with IsoQuant using long reads. *Nat
  Biotechnol* 2023;41:915–918. doi:10.1038/s41587-022-01565-y.
- [35] Tang AD, Soulette CM, van Baren MJ et al. Full-length transcript characterization of SF3B1 mutation in chronic
  lymphocytic leukemia reveals downregulation of retained introns. *Nat Commun* 2020;11:1438.
  doi:10.1038/s41467-020-15171-6 (FLAIR).
- [36] Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. *Bioinformatics* 2012;28:2184–2185.
  doi:10.1093/bioinformatics/bts356.
- [37] McKenna A et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA
  sequencing data. *Genome Res* 2010;20:1297–1303. doi:10.1101/gr.107524.110 (hard filters FisherStrand, MQ,
  MQRankSum, ReadPosRankSum).
- [38] Robinson JT et al. Integrative genomics viewer. *Nat Biotechnol* 2011;29:24–26. doi:10.1038/nbt.1754.
- [39] Jugel U, Jerzak Z, Hackenbroich G, Markl V. M4: a visualization-oriented time series data aggregation. *PVLDB*
  2014;7(10):797–808. doi:10.14778/2732951.2732953.
- [40] Mose LE, Perou CM, Parker JS. Improved indel detection in DNA and RNA via realignment with ABRA2.
  *Bioinformatics* 2019;35:2966–2973. doi:10.1093/bioinformatics/btz033.
- [41] Dolzhenko E, van Vugt JJFA, Shaw RJ et al. Detection of long repeat expansions from PCR-free whole-genome
  sequence data. *Genome Res* 2017;27:1895–1903. doi:10.1101/gr.225672.117 (ExpansionHunter).
- [42] Dobin A et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2013;29:15–21.
  doi:10.1093/bioinformatics/bts635.
- [43] Danecek P et al. Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10:giab008.
  doi:10.1093/gigascience/giab008.
