# Array-CGH module: what this branch adds, and how to merge it

Branch `dev-WD`. Base: `e0f80a4` (*SV browser: transcript HGVS placed by VariantValidator first, Ensembl as fallback*), which
`dev` contains. The module is one new page, **one level under the SV browser** (SV browser → **Array-CGH ▸**). It adds
files of its own and touches four shared files in a few lines. Nothing else of the viewer changes. With no array-CGH
file loaded, the only visible change is the **Array-CGH ▸** button in the SV browser's header.

## What the page does

- **Examples.** *Examples ▾* (empty page and frame bar, as the SV browser offers its examples) adds five invented patients (`cghExamples.ts`, CytoGenomics text, GRCh37) through the normal import question: 22q11.2 LCR22A-D deletion (HIGH), 1q21.1 distal BP3-BP4 duplication (AChro-Puce PIEV, ClinGen TS 2 region: MODERATE), terminal 16q24.3 loss breaking in ANKRD11 intron 4 + terminal 9p24.3 gain (der() hint, ANKRD11 5′ part cut, HIGH), and EXAMPLE-PROBES-ONLY (a probe file alone: 7q11.23 loss, 15q11.2 BP1-BP2 gain; a probes-only patient), EXAMPLE-MOSAIC-T12 (a mosaic trisomy 12 written as a 12p gain, a 12q gain and a gain nested in 12q, mean log ratio ~+0.25: one CNV group); *All (five)*. Coordinates checked on UCSC hg19 (ncbiRefSeqCurated, ClinGen tracks). Each table comes with invented probes (`cghExampleProbes.ts`): a `CGH_Log_Ratio.txt` written at load time (no embedded data; seeded mulberry32 per sample name, the same file each time; ~59.7k probes over the GRCh37 lengths, none in centromeres / acrocentric p arms / chrY heterochromatin, noise SD 0.15, 0.3 % NaN, a female's chrY noise SD 0.3), each call given exactly its #Probes evenly from its start to its stop with their mean set to its Mean Log Ratio (so `meanIn` gives the file's #Probes and mean; terminal calls have no probe beyond them), given as a second file to `cgh.add` (attached by its header's sample name). EXAMPLE-DEL22Q11's probes also hold a 12-probe 3p14.2 loss (−0.62) not in its table: one probe call *not in the file*. Checked in node (`detect`, defaults): every table call found with its probes and mean, the extra loss found, no other event; the probes-only example gives its two events. Public build: never real data there.
- **Quick guide.** *? Quick guide* in the page's bar: an interactive tour of 27 steps on the examples (`CghGuide.tsx`): input files (text only), examples, import question (DLRS, minimums, sex), tabs, filters, More filters, priority, why, browser below, the call's ISCN, the call's probe panel, der() hint, gene cut, PIEV, classify / note, the patient's ISCN card, a CNV group (EXAMPLE-MOSAIC-T12), whole genome, a chromosome (opened through a `cgh-chrom` window event `CghGenome` listens to), log-ratio range + show probes, probe calls (*n not in the file*), probes-only patient, lab list, retry, export, ⓘ methods (`data-cgh-disclaimer`: a complementary tool, tested on Agilent 60K), Esc. An outline on the target (the page dimmed around it, the target's box clipped by its scrolling ancestors and scrolled into view again while hidden; `maxH` outlines only a tall target's top) and a popover beside it (accent border and title band, a pointer toward the target, never over the target: the side kept while it fits, else the on-screen spot covering least of it; it slides (220 ms) rather than jumps, keeps its place up to 1.5 s while a target is missing, then goes bottom centre; draggable by its title); each step does its own action (load the examples, Analyse, open More filters, a call, a tab, the whole genome, a chromosome): the user only clicks *Next*; Esc / *Skip tour* ends it (the store's `help` keeps Esc from acting on the page).
- **Esc and leaving.** Esc on the page first clears the view (`cgh.escape()`: `shown` cleared, the store's `viewReset` counter folds the detail in `CghFrame`, the chromosome view in `CghGenome`, the call window and the whole genome in `CghPage`; filters, tab, sort, table scroll kept; the page tells the store what it holds open through `setViewBusy`); on that base view Esc, or the bar's *← SV browser*, opens *Leave Array-CGH…?* (`leaving`, `LeaveDialog` in `CghPage.tsx`: Cancel focused, Esc = Cancel, the number of patients in `dirtyKeys`; with no patient loaded it leaves at once). The store is module-level, so leaving keeps the patients for the tab's life. A menu open, the guide, the import / column questions keep their own Esc; Esc in a text field (search, note) only blurs it. Targets found by `data-cgh-*` attributes. The former screenshot guide (`src/assets/cgh-guide/`, `docs/cgh-demo/`) is removed.
- **Input.** Files are read in the page, never uploaded:
  - CytoGenomics interval tables (`.txt`), CSV / TSV of any layout (a column-mapping dialog when the columns are not
    recognised, remembered by column names, only when the header row looks like one);
  - `.xlsx` (first sheet), BED, SV VCF (`.vcf`, `.vcf.gz`);
  - the mapping dialog also asks, besides chromosome / start / end (or one region column), for a type column, OR one or
    several log-ratio columns (`MapSpec.lrMore`; guessed when named log / lr / ratio, or gain / loss with values within ±5: gain / loss columns merged as the non-zero one, typed by a column name
    saying gain / loss, else the sign), OR a copy-number column, and names what is missing; with no type: copy number
    (< 2 loss, > 2 gain; 1 normal on a male's X outside the PARs and on Y; chrX x1 / x2 outside the PARs always kept, x1 as a loss, x2 as a gain, `CghCnv.cnSex`, with or without a sex in the file (a warning in the question when none), the one normal for the current sex hidden by `passes` / `hiddenBySex`, so a corrected sex brings it back; a range (x2~3, a male's X x1~2) a gain when `copiesHi` > normal, else a loss, never normal), else the log ratio's sign; rows with neither skipped, counted
    in the question; specs saved before `lrMore` load as before;
  - CytoGenomics' automated interval-based report (`IntervalBasedReport.xls`, tab-separated text despite its name, or `.xlsx`), recognised by its header (`AberrationNo`, `Chr`, `Start`, `Stop`, `#Probes`, `Amplification` / `Gain` / `Loss` / `Deletion`) and read without the mapping dialog (`parseAberrations` in `cghIntervals.ts`): one log ratio from the four columns (the non-zero one, the largest |value| if several), the type from the column (else the sign), `pval` (`CghCnv.pval`), genes, the build from `Genome:` (else a header line naming one, else GRCh37 with a warning in the question), the sample from its name line (else `Array Design Info`); copies null (taken from the log ratio, see below; ISCN written by `cghIscn.ts`); exported as it came plus Classification / Comment columns;
  - zips holding any of these (one patient each): tables found by content, the other files of an automated export left aside and listed in one note (*ignored: n other files*); a table needing the mapping is asked only when the zip holds no table read by itself, and the zip's probe file then waits with it (`toMap[].probes`) and goes to its sample(s) in `mapped` (no probes-only patient; skipped: a probes-only patient); a table read again with other columns (*columns…*: new file id, new sample keys) has its probes moved to the samples of the same name (`moveProbes`);
  - probes' log ratios (`CGH_Log_Ratio.txt`, `cghProbes.ts`): kept as bytes with the zip's patient (alone: with the loaded sample of the `Normalized log ratio(<name>)` header, else a **probes-only** patient, below), read on first `cgh.probes(sampleKey)` into `Uint32Array` start / end and `Float32Array` log ratios per chromosome (NaN kept) and probe names; `cgh.hasProbes`; `dlrs` (SD of consecutive differences / √2) in the import question and the call card. **Log ratio (probes)** panel (`CghProbePanel.tsx`) in the call's card above the browser: a canvas of the window's probes (binary search on the sorted starts, plus those starting before the window but reaching into it; `visible` / `runningMedian` in `cghProbes.ts`); the call's probes' mean (`meanIn`) over the probes whose midpoint lies in the call, ends included (was: start in the call, which counted a probe straddling the call's end: 1 of 105 real calls; now the file's #Probes on all 105), colours black / red / green at ±0.2, moving median 3 / 7 / 15, call shaded with its mean, ±0.2 dashed, hover tip; its own zoom / pan (the browser's window is internal to `SvGenomeBrowser`: not synced, that file untouched). Both open centred on the call ± max(50 % of its size, 100 kb) (`callMargin`): `drawOf` in `CghPage.tsx` gives the browser that window (less the browser's own 12 % pad) and rebuilds `rec` over it, no change to `SvGenomeBrowser`; a gain's view also takes the tandem copy, so it is wider. Hovering the call's shading away from the dots: the call's tip (`callLines` in `CghGenome.tsx`).
  - **Whole genome** (`CghGenome.tsx`): one patient (the tab's; in All patients a list in the panel), chr1-22, X, Y by GRCh37 length; with no probe file (interval table, BED, VCF, examples) a short panel instead: the listed calls as bars over the chromosomes to scale (click: open the call), no log-ratio axis / range / layout / chromosome view (an overview of where the calls lie is worth its 40 px; empty axes were not); probes drawn one mark per pixel cell and colour (a bit mask: each colour present in a cell drawn, a red probe no longer hidden by a later black one; every non-NaN probe counted, `data-drawn` = 61,246 − NaN on the three real arrays; 61k probes: ~20-50 ms), the median of each pixel line; file calls listed by the filters (red / blue) and probe calls (purple) as bars; layout *horizontal* (side by side, bars above) or *vertical* (CytoGenomics-like columns: position down from pter, log ratio across, bars beside, one row or two under 1,000 px), `state.view.layout`; a chromosome clicked → its view (a call bar clicked, in any view: the call), always horizontal, which **is** `CghProbePanel` without a call (`at`, `calls`, `onCall` props: every probe as a dot, moving median 3 / 7 / 15, segments, ⊖ ⊕ ← →, slider, drag, Ctrl / ⌘ + wheel, call bars above, the chromosome's **ideogram** under the plot: `Ideogram` in `CghProbePanel.tsx`, the SV browser's drawing redone with its exported `STAIN` colours, bands from the loader's `getCytobands` (UCSC, cached), names where they fit, *cen*, the window boxed, click / drag moves it; the call's probe panel has it too; ← → for the neighbours), a call bar clicked there opens the call; **log-ratio range** (`LrScale`, `state.view.scale`, one setting for the genome, chromosome and probe panels, localStorage `sashimi.cghView`): small buttons stacked in the plot's left margin beside the axis labels (`L` 60; vertical layout: at the foot of the first column, by its log-ratio axis), ⊕ / ±2 / ±1 (default) / ±0.5 / ⊖ out to ±3 / 4 / 6 / 8; **show probes** (`DotsToggle`, `state.view.dots`, same storage): off, no dots (no probe hover), the median line (genome: per pixel only), segments, call shading and bars kept; a probe beyond it drawn as a triangle at the edge (`clipMark`), *clipped at ±x* in the captions, colours still at ±0.2; hover tips for probes and bars; calls only without a probe file. It takes the call window's place: open while no call is, a strip once one is; opened again, the call window folds to a strip (one at a time). Methods of DLRS, probe calls, estimated copies, mosaic and this view: the *Methods* part of the page's ⓘ panel (`ProbeMethods` in `CghProbeCalls.tsx`, saying when no probe file is loaded that the probe methods do not apply). Nothing empty is drawn: without probes no probe panel, *Probe calls* line, DLRS or range control; the *Probes* / *Log R* columns and the import question's minimum of probes only when a file gives them; `fmtBp` / `sizeText` shared from `cghIntervals.ts`.
  - **Probe calls** (`cghSegment.ts`, pure; `CghProbeCalls.tsx`): events detected from a patient's probes, a second opinion. Per chromosome, NaN skipped, circular binary segmentation (Olshen et al. 2004) with a fixed threshold `t` on the arc statistic, σ = MAD of consecutive differences / √2 taken once per sample (no permutations: O(n²) per stretch up to `EXACT_MAX` = 10,000 probes, about half a second for 61k probes; a longer stretch (180k / 1M arrays) searched on 2,500 blocks then exactly around the best arc, at the stretch's ends and over every arc of up to two blocks (O(n·b): a 5-probe loss inside a 32-probe block is no longer diluted; invented 80k-probe chromosome, 5 probes at −0.8, σ 0.15, 20 random places: 13 / 20 found before, 20 / 20 after, no extra; the three real 60K arrays, no chromosome over 10,000 probes, give identical JSON before / after): 1M invented probes in about a second; one patient at a time via `cgh.runProbeCalls` in a Web Worker, `cghProbeWorker.ts`, inlined by Vite, which also reads the probe file: `cgh.probes` returns null while it is read and the page redraws when done, `cgh.probeGen` versions a patient's probes for memos and results; new probes during a reading drop it (`dropProbes` clears `probeParsing`) and the new bytes are read; a worker error is tried again on the page, and a failure there is settled, said on the patient's *Probe calls* line (*probe calls failed: …*, `cgh.probeErrorOf`) and not retried until the settings, sex or probes change; without a worker, on the page); segments called with ≥ `minProbes` (3) probes and |mean| ≥ `minLr` (0.25), a < `minProbes` segment between two of one direction absorbed, adjacent same-direction merged; chrY PAR left out (the chrX PAR probes again; PARs per build, `PARS` in `cghMosaic.ts`, GRCh38 PAR1 chrX / chrY 10,001-2,781,479, PAR2 chrX 155,701,383-156,030,895 / chrY 56,887,903-57,217,415, checked by sequence identity on the UCSC API), chrY read for a male only (a female's gives many false events). Settings in `state.probeOpts` (localStorage `sashimi.cghProbeCalls`), sensitivity high / normal / low = t 6 / 8 / 10. Events not in the file (`agrees`: same type, reciprocal overlap ≥ 50 %, chrY PAR mapped to chrX; the file's calls themselves carry no mark, being made from the same probes) listed per patient (*▸ n not in the file*), drawn as a `CghCnv` with `fromProbes` (`probeCnv`), never exported; segments a purple line in the probe panel. Defaults calibrated locally on real arrays against CytoGenomics ADM-2 calls (data not in the repo): at t 8 every file call of ≥ 3 probes found, 1-2 extras per array, no 1-2-probe call reproduced; t 10 lost calls of 3-5 probes, t 6 tripled the extras. Invented recall (`scripts/cgh-probe-fixture.mjs`: del -0.9 × 5, dup +0.55 × 200, mosaic -0.4 × 30, male X dup +0.9 × 30): 4/4, no extra.
  - **Probes only** (`cghProbeOnly.ts`): a probe file with no matching patient (alone, or in a zip without a table) becomes a `CghFile` of `format: 'probes'` (no lines; name from the header, else the file name; build from the AMADID line, `probeBuildOf`); its sample's `cnvs` are its probe calls (`probeOnlyCnvs`, ids `pc|<key>|chr:start-end`, so classes / notes / annotations stay on a call whose position is unchanged), rewritten by `applyProbeOnly` in `cghStore.ts` whenever its probe calls are recomputed (settings, sex; also while its question is open, which shows *probes only*, the sex, `ProbeSettings` and the DLRS). Rows like an interval report's (copies null: estimated from the log ratio; band from UCSC); a *probes only* tab tag and a banner (`PROBE_ONLY_WARNING`); not in the grey *not in the file* list. Export: `exportProbeCalls`, a new `<name>.probe-calls.txt` (header lines: source, build, settings, date, and the rows left without a band when the UCSC cytobands could not be read; the bands are read for every call first, per chromosome, cached by the loader; Chromosome … Comment). Its calls are not drawn again as probe calls in the whole genome. A table of the same name added later (through `add`, not a mapped table) gets a copy of the probes and replaces it at `accept` (`pending.replaces`, then `cgh.remove`); classes / notes of the probe calls are not carried over, said in the question (simplest safe choice: Cancel lets them be exported first); a probes-only patient still waiting in the question is dropped at once. Checks: `scripts/cgh-probeonly-check.ts` on the fixture's `--split` files (probe file alone, zip with only probes and dummies, table alone). Real arrays alone (local, not in the repo): every ≥ 3-probe CytoGenomics call found (4/4, 10/10, 9/9), 1 extra each.
  - **Other probe files** (`probeLayoutOf`, `parseProbes(bytes, spec)` in `cghProbes.ts`; `ProbeSpec` per sample, `probeSpec` map beside `probeBytes` in the store, sent to the worker): Agilent FE-like, Illumina final report (`<sample>.Log R Ratio`), ChAS-like (`Log2Ratio` before `WeightedLog2Ratio`), CNVkit `.cnr` (0-based starts moved to 1-based; `.cnr` added to the pickers and zips), DNAcopy-like inputs, recognised by header (any of tab / `;` / `,`, `#` or `key: value` lines before; ≥ `PROBE_MIN_ROWS` = 2,000 rows; no type / size / probe-count / `num.mark` / `seg.mean` column); one patient per log-ratio column (`lrColumnSample`: the name in the header, else a `Sample:` line, else the file name); 23 / 24 → X / Y, MT / XY / 0 left out, empty / NA / NaN / null → NaN, chromosomes and rows sorted; the layout and the positions' base said in the question. Then exactly the CGH_Log_Ratio path (attached by sample name, else probes-only). Not recognised and ≥ 2,000 lines: the mapping dialog opens in **probe mode** (a first radio *Calls* / *Probe log ratios*; chromosome, position or start, end, log-ratio column(s) each named, probe name, 0-based, build; `guessProbeCols` by names then content; only the first 300 rows read; `cgh.mappedProbes`), remembered by header names in `sashimi.cghProbeMaps` (column indexes and options only, no sample name). Checks: `scripts/cgh-probe-layouts-check.ts` on `cgh-probe-fixture.mjs --layouts DIR` (same probes as the CytoGenomics fixture in five layouts: same probes, NaN, 4/4 events; a 2-sample CSV; the odd `;` / decimal-comma file not recognised, its columns guessed); real CGH_Log_Ratio.txt (local): probe calls byte-identical before / after, and the generic reader on them gives the same probes and calls. `--hg38 DIR`: an invented GRCh38 male (chr22 loss, chrX PAR1 gain with chrY PAR1 probes carrying it, chrX non-PAR gain) with its hg38 table. Fixed with it: the table's *Position* header names the files' build (was always GRCh37); a male's chrX / chrY call inside a PAR is no longer drawn with one copy (`one` in `CghPage.tsx`, PARs of the call's build). No-worker check (Chromium, `window.Worker` deleted or throwing): same probe calls as with the worker; the main thread blocks up to ~0.5 s per 60k-probe file (~1.3 s for three at once).
  - Encodings: UTF-8 with or without BOM, UTF-16, Windows-1252. A binary `.xls` is refused with a clear message.
- **Import question.** Minimum number of probes (3) and minimum size (log slider: any, 10 kb to 5 Mb, or a value typed in
  bp / kb / Mb, also in the filter bar); the 200 largest calls annotated by
  themselves, the others on demand (*annotate*); each patient's calls under the minimums revealed on demand (*+ N more*); sex per patient.
- **Table.** One tab per patient plus *All patients*, filters (series: the number of other patients carrying the call
  compared with = ≥ ≤ > < n; a sample name loaded twice is one patient, said at import), gene search (file genes and
  RefSeq genes read), sort by priority; *Probes* and *Log R* (with the mosaic tag) in two columns; column widths set by
  dragging a header's right border (`useColWidths` in `CghPanel.tsx`: pointer events, min 36 px, then `table-layout:
  fixed` with every column's measured width; double-click: automatic again), kept in localStorage `sashimi.cghColumns`
  (column ids and pixels only).
- **Per call, read on GRCh37** by the SV browser's own routes (UCSC, PanelApp, gnomAD):
  - RefSeq genes; PanelApp ratings and modes of inheritance (asked 25 genes a request);
  - ClinGen HI / TS of genes and of ISCA regions; gnomAD LOEUF;
  - AChro-Puce recurrent CNVs; DGV Gold (3 relations);
  - exons cut and reading frame; terminal / arm / chromosome labels; rearrangement hints;
  - mosaic fraction (AChro-Puce formula; no DLRS input); a copy number from the log ratio when the file gives none (`estCopies` in `cghMosaic.ts`: nearest log2(n/2), x0 at −2; male X / Y normal 1), written as any copy number (*del x1*, *copies x1*, no question mark), *mosaic ~n%* only when the AChro-Puce fraction of a one-copy change is < 80 %; **ISCN** (`cghIscn.ts`, array ISCN as in ISCN 2024: bands pter→qter, one band once, the band from the UCSC cytobands at sub-band precision (`bandSpan`, `iscnBandOf`: `CnvAnn.band` now always read; the file's band only without them), positions with commas, `x1[0.48]` mosaic, `(21)x3` whole chromosome; `CghIscnView.tsx`: on the call card and detail with *copy*, in the row tooltip, the file's own ISCN beside when it differs; an **ISCN — patient** card per patient tab (under the filters, above the table; `normalIscn`: `arr (1-22,X)x2` female, `arr (X,Y)x1,(1-22)x2` male, no line when sex not given) from the calls of the chosen classes (`state.view.iscnClasses`, default P / LP / VUS), 1-22, X, Y then position, the normal result by sex when none; checks in `scripts/cgh-iscn-check.ts`, invented calls); a disclaimer in the ⓘ panel (copy numbers inferred from the log ratio; everything to be checked in the analysis software); recurrence in the loaded series; the laboratory's benign
    list (BED, kept in `localStorage`).
- **CNV groups** (`cghGroups.ts`, pure; `groupsOf` / `activeGroups` in `cghStore.ts`). A patient's calls of one direction on one chromosome (analysed: over the minimums, not suppressed, not hidden by the sex; at least two) form a group when together they cover ≥ 90 % (`GROUP_COVER`) of an arm's probe-covered span, or of both arms (a whole chromosome; an acrocentric 13/14/15/21/22: its q arm alone). The span left out: the centromere and the UCSC `gvar` / `stalk` bands (embedded per build, read 2026-11-25, `BLIND`), the PARs on X / Y, and GRCh37's 9p11.2-p13.1 (from 38.7 Mb, as `cghExampleProbes.ts`). A call inside a larger one is a **nested** member (ADM-2 reports nested intervals): shown with `└╴`, not counted in the probes nor weighted in the mean log ratio (weighted by probes, by size when a count is missing). Copies: the members' own when they agree, else `estCopies` on the group as one call (mean log ratio), mosaic fraction from `mosaicOf` (< 80 %). Suggestion (`groupSuggestion`, words only, never a class): *trisomy 8*, *mosaic trisomy 8 (~29 %)*, *partial monosomy 4p*, *partial tetrasomy 12p*, x0 nullisomy … x6 hexasomy; X read with the sex (female loss: *monosomy X (45,X) — or a sex mismatch with the reference*; male loss: expected against a female reference; male gain: *XXY?*; female gain: trisomy / tetrasomy X; no sex: both readings), Y likewise (a Y gain ≥ +1.5 in a male: a female reference). ISCN (`groupIscn`): `(8)x3[0.29]` for a chromosome, an arm by the first member's first band and the last member's last band with the span (`callIscn` on the group as one call). Table: an outlined block (`GRP_*` classes): a summary row (GROUP tag, span, size, copies, probes, weighted log ratio, suggestion, arms covered, members, ISCN + copy, its own class `classes['grp|<key>|<chrom>|<kind>']` and note, *→ calls* gives its class to every member, *ungroup*), then its members as usual rows; it takes the place of its first member listed (by position or priority), all members kept together, shown when any member passes the filters. *ungroup*: `ungrouped[id]` (this tab), a line *group* to restore it. Patient ISCN: a group with a class ticked replaces its members' entries. Whole genome: a violet outline over the group's bars, the group in the bar's tip. Export: a **column** `CNV group (Sashimi viewer)` (`OWN_GROUP_HEADER`, never guessed when read back, rewritten at the next export) on each member's row: `label: suggestion [group class]`, added only when a call of the file is in a group; VCF: INFO `CGH_GROUP`; probe-calls export: a *CNV group* column. Chosen over extra group rows, which no format holds cleanly.
- **Priority.** HIGH / MODERATE / LOW from listed rules (`cghScore.ts`). It is not a classification. Per Riggs et al. 2020: a ClinGen HI 3 / TS 3 region is HIGH only wholly in the call (2A, share ≥ 0.999), MODERATE at 80-99 % (2B: check the critical region); an HI 3 gene losing only UTR exons (2C-2 / 2D-1) or only its last exon (2D-2 / 2D-3) is MODERATE (`sparing()`, from `CnvExons.cds` and the exon ranks; genes unknown to RefSeq stay HIGH).
- **Retry.** Failures kept per call and source (`CnvAnn.failed`; PanelApp per gene): *↻ retry (n)* in the column header (per patient; a menu in All patients) and a *retry* per source in the detail re-ask only what failed (`retry()` in `cghAnnotate.ts`); an unexpected error ends the call on what was read.
- **Classification, note, CNV-Hub link.** **Export classified table** writes each file back with Classification and
  Comment columns (a VCF: INFO `CGH_CLASS` / `CGH_NOTE`; an xlsx: a tab-separated `.classified.txt`), plus an ISCN column
  (VCF: INFO `CGH_ISCN`) only when the file has none: a CytoGenomics interval table keeps its own. The column added is named `ISCN (Sashimi viewer)` (`OWN_ISCN_HEADER`): read back, `guessCols` ignores it (not the file's ISCN) and the export rewrites it; a BED gets class, note and ISCN as three last columns.
  A class not changed in the page is written as the file had it (`s.classes[id] ?? c.classification`, not `classOf`, which
  translates *bénin* to *Benign* for display); a file mixing line endings (CR LF title, LF after) keeps each one (`CghFile.eols`).
  Several files: a list with one *export* per file and *export all* (one download at a time, each file marked as it starts).
- **The browser below.** The call clicked is drawn as DEL / DUP in the SV browser's own `SvGenomeBrowser`, with the
  same tracks, plus **ClinGen CNV regions** (only here, through the optional `extraTrack` prop): HI 3 / TS 3 dosage
  regions (red / blue, link to ClinGen), the ClinGen recurrent CNV list v2.0 (hatched orange) and its breakpoint / LCR
  segments (grey), folded by default (*▸ show recurrent CNVs / LCRs*), a *PIEV* tag where AChro-Puce classes the region PIEV (≥ 80 % reciprocal, same type). Embedded data
  (`cghClinGenData.ts`, read 2026-11-20 from the UCSC API and https://ftp.clinicalgenome.org/; no licence stated in the
  BED files), no request. The call's detail adds *ClinGen recurrent CNVs*: regions matched (share of region / of call)
  and whether each end lies in a breakpoint / LCR segment (± 50 kb). Priority rules unchanged.

Privacy: patient data stay in memory; `localStorage` only holds settings, column names of mapped tables and the lab list.

## Files

**New (the module, all in `src/components/sashimi/`):**

| File | Role |
|---|---|
| `CghPage.tsx` | the page (bar, frame, browser below); `CghOpenButton`, `CghPageSlot` used by the SV browser |
| `CghPanel.tsx` | the frame: tabs, filters, table, detail, import question |
| `CghMapDialog.tsx` | column-mapping dialog |
| `cghStore.ts` | state (files, classes, notes, sex, filters, annotations) |
| `cghIntervals.ts` | interval table / BED parsing, export, classes, CNV-Hub link |
| `cghMapping.ts` | delimiter and header detection, column guessing, remembered mappings |
| `cghFormats.ts` | decoding (BOM, UTF-16, cp1252), xlsx reading, `.xls` detection |
| `cghVcf.ts` | SV VCF reading and export |
| `cghAnnotate.ts` | per-call annotation queue (genes, dosage, DGV, AChro-Puce, bands, PanelApp batches) |
| `cghScore.ts` | priority rules |
| `cghMosaic.ts` | mosaic fraction (port of the AChro-Puce calculator) |
| `cghGroups.ts` | CNV groups: calls covering a whole arm / chromosome, their metrics, suggestion and ISCN |
| `cghSeries.ts` | recurrence across loaded patients |
| `cghTopo.ts` | terminal / arm / chromosome labels, rearrangement hints |
| `cghLabBenign.ts` | the laboratory's benign list |
| `cghExons.ts` | exons cut and reading frame |
| `cghExamples.ts` | built-in invented example patients |
| `cghExampleProbes.ts` | the examples' invented probe files, generated at load time (seeded) |
| `cghProbes.ts` | probes' log ratios (`CGH_Log_Ratio.txt` and other probe-level layouts): detection, lazy reading into typed arrays, DLRS, window search, remembered probe mappings |
| `CghProbePanel.tsx` | the *Log ratio (probes)* panel above the browser |
| `CghGenome.tsx` | the *Whole genome* panel (genome and chromosome views), the call tooltip |
| `cghSegment.ts` | probe calls: segmentation (CBS, fixed threshold), calling, agreement with the file's calls |
| `cghProbeOnly.ts` | probes-only patients: a probe file alone as a patient, its probe calls as its calls, the `.probe-calls.txt` export |
| `cghProbeWorker.ts` | the Web Worker reading probe files and computing probe calls off the page |
| `cghIscn.ts` | array ISCN (as in ISCN 2024) of a call and of a patient |
| `CghIscnView.tsx` | the ISCN shown: call card, detail, row tooltip, *ISCN — patient* card, copy buttons |
| `CghProbeCalls.tsx` | the *Probe calls* line: settings, per-patient lists, the *probes: seen / not seen* mark |
| `CghGuide.tsx` | the interactive tour (*? Quick guide*) |
| `CghClinGenTrack.tsx` | the *ClinGen CNV regions* track (`clinGenTrack(build)`) and the detail's recurrent-CNV lines |
| `cghRecurrent.ts` | ClinGen regions on a chromosome, AChro-Puce PIEV tags, a call's recurrent matches and LCR ends |
| `cghClinGenData.ts` | generated data: HI 3 / TS 3 ISCA regions, recurrent CNV list v2.0 (both builds) |

Also `scripts/cgh-probe-audit-check.ts` checks that every probe is taken into account (invented data: header, shuffled rows, NaN, duplicate positions kept, chrM / unknown chromosomes left out, DLRS in position order with NaN bridged, every non-NaN probe in one CBS segment, exact and hybrid, chrY rule, probes at and across a call's ends, the window's edges, moving-median ends); `scripts/cgh-fixes-check.ts` checks the review fixes on invented data (hybrid CBS, an 80k-probe chromosome with a 5-probe loss, chrX without sex and with a stated male, copy-number ranges, header-less ISCN column, the page's own ISCN column read back, BED's three columns, the no-band count of a probes-only export, PARs per build, log-ratio column guessing, report without `Genome:`); `scripts/cgh-store-check.ts` checks the store with a stand-in worker (a zip whose table waits for its columns keeps its probes, *columns…* moves them, new probes during a reading, a worker error); `scripts/cgh-iscn-check.ts` checks the ISCN and its export column on invented calls; `scripts/cgh-probe-fixture.mjs` writes an invented zip (report + probes, known events) to test the probe calls (`--split DIR`: also the probe file alone, a zip of it with dummies, the table alone); `scripts/cgh-probeonly-check.ts` checks the probes-only mode on them; `scripts/cgh-clingen-data.mjs` writes `cghClinGenData.ts` from the downloaded UCSC JSON and ClinGen BED files.

Also new: this file, `ARRAY-CGH.md`.

**Shared files changed (whole footprint: 4 files, about 75 lines):**

| File | Change |
|---|---|
| `src/components/sashimi/SvBrowserDialog.tsx` | imports `CghOpenButton`, `CghPageSlot`, `cgh`; the button in the header; the page slot; Esc handler: the array-CGH questions, guide (`help`) and leave question (`leaving`) answer it themselves, an open page menu closes first, a page text field is only blurred, else `cgh.escape()` (clear the view, then ask to leave); with the page closed, Esc closes the SV browser as before |
| `src/components/sashimi/SvGenomeBrowser.tsx` | `export const MOI` (was `const`); `export const STAIN` (was `const`; the array-CGH ideogram's colours); optional loader field `getDosageRegions` in `SvGeneLoaders`; optional prop `strongEvent` (off by default: shades the event, bolder breakpoints); optional prop `extraTrack` + exported `SvExtraTrack` interface (a caller's track under AChro-Puce: state `showExtra`, `extraLaid` / `yExtra` in the layout before `yRefBand`, its side label, its `<g>` before the reference band, a `Pill` after AChro-Puce's and its legend; absent in the SV browser, so nothing drawn there); **bug fix, also visible in the SV browser**: the SVG's ids (`svb-clip`, gradients, pattern) made unique per instance with `useId` (`sid(...)`): the array-CGH page's browser drawn over the SV browser's resolved `url(#svb-clip)` to the other, shorter one, so a grown ClinVar SV track (*▸ all SVs*) hid all below it; bug fix, also visible in the SV browser: in "▸ all SVs", ClinVar SVs that find no free lane among the 60 are no longer drawn over others but counted ("+n SV not drawn (zoom in)") |
| `src/standalone/ucsc.ts` | new `getDosageRegions` + `DosageRegion` (ClinGen ISCA regions; `getDosageScores` untouched) |
| `src/standalone/svLoaders.ts` | wires `getDosageRegions` |

**Docs and build:** `docs/manual.md` (an *Array-CGH* section in *Genomic DNA*, about 100 lines, and one sentence in *Privacy*), `README.md`
(one bullet), and the rebuilt `sashimi-viewer.html`. No dependency added: `fflate` was already in `package.json`, and
`package.json` / `package-lock.json` are unchanged.

## Checks done on this branch

- The SV browser, before (`e0f80a4`) and after, gives pixel-identical screenshots on three searches, apart from the
  button.
- The real CytoGenomics export (36 calls) parses identically through every change and exports byte-identically.
- Unit checks for formats, scores, mosaic, DGV, exons; Playwright runs on real data.
- A 120-call load test annotates in about 70 s with no request refused.

## Merging into `dev` (to do by hand)

**Update.** `dev` has since merged this branch (`a720409`). As of `origin/dev` = `3248eee`, merging the later commits
(size slider, series operator, DLRS removed, the SV browser's unique SVG ids) gives **one conflict**,
`sashimi-viewer.html` (take either side, rebuild); `SvGenomeBrowser.tsx`, `CghPanel.tsx` and `docs/manual.md` merge
automatically. The *ClinGen CNV regions* commit adds `extraTrack` lines in `SvGenomeBrowser.tsx` next to the
AChro-Puce state, layout, label, drawing, pill and legend lines: if `dev` changed those, keep both sides' lines. The
first merge, for the record:

As of `origin/dev` = `6515ae6` (checked in memory with `git merge-tree`, nothing written), a merge gives **two
conflicts**:

1. **`sashimi-viewer.html`.** It is a build output: take either side, then rebuild (step 5).
2. **`src/components/sashimi/SvGenomeBrowser.tsx`.** `dev` changed the same three places:
   - **Props interface.** `dev` added `onGeneSummary?`, this branch `strongEvent?` (later also `extraTrack?`, with the `SvExtraTrack` interface after `Props`). Keep **all** fields.
   - **The `MOI` constant.** `dev` moved it to the new `SvGeneSummary.tsx` (`export const MOI` there) and removed it from
     `SvGenomeBrowser.tsx`. Keep `dev`'s side: **drop** this branch's `export const MOI` block. Then, in
     `CghPanel.tsx`, change
     `import { MOI, type SvGeneLoaders } from './SvGenomeBrowser';` to
     `import { type SvGeneLoaders } from './SvGenomeBrowser';` plus `import { MOI } from './SvGeneSummary';`.
   - **Function signature.** Keep both destructured props: `…, variant, onGeneSummary, strongEvent, extraTrack }: Props)`.
   - The other two edits of this branch in that file (`getDosageRegions?` in `SvGeneLoaders`, and the `strongEvent`
     lines in the drawing, next to the comment *breakpoints through the reference*) merged without conflict. Check that
     they are still there.

`SvBrowserDialog.tsx`, `docs/manual.md`, `ucsc.ts` and `svLoaders.ts` merged automatically at that date. If `dev` has
moved further, re-run the dry run below.

### Steps

1. Look at what a merge would do, without changing anything:

```bash
git fetch origin dev dev-WD
```

```bash
git merge-tree --write-tree --name-only origin/dev origin/dev-WD
```

The lines after the first one, if any, are the conflicting files.

2. Merge on a local branch:

```bash
git checkout -b merge-cgh origin/dev
```

```bash
git merge --no-ff origin/dev-WD
```

3. Resolve `SvGenomeBrowser.tsx` and `CghPanel.tsx` as described above. For `sashimi-viewer.html`, take either side:

```bash
git checkout --theirs sashimi-viewer.html
```

4. Check that the module only touches what it should. The shared files should show only the lines listed above:

```bash
git diff origin/dev --stat
```

```bash
git diff origin/dev -- src/components/sashimi/SvBrowserDialog.tsx src/components/sashimi/SvGenomeBrowser.tsx src/standalone/ucsc.ts src/standalone/svLoaders.ts
```

5. Type-check and rebuild the single file:

```bash
npm ci
```

```bash
npm run build
```

6. Test by hand in `sashimi-viewer.html`:
   - Open the SV browser on `#sv=DEL%2020:10489766-10783340`: the browser is as before, with the **Array-CGH ▸** button,
     and *Genes on the rebuilt chromosome* still shows its gene summary (from `dev`'s `onGeneSummary`).
   - Click **Array-CGH ▸**, add an interval table, analyse it, click a call: it is drawn below; Esc clears the view, Esc
     again asks to leave (Cancel / Leave).
7. Commit the merge, then push the branch you choose.

### If a conflict appears elsewhere

- **The module's own files.** `dev` has no file named `Cgh*` / `cgh*`, so a conflict there means `dev` started its own
  version: compare before choosing.
- **`docs/manual.md`.** Keep both texts. This branch's text is the *Array-CGH* section at the end of *Genomic DNA*,
  plus one sentence in *Privacy*.
- **`SvBrowserDialog.tsx`.** This branch's lines are:
  - the imports `CghOpenButton, CghPageSlot` and `cgh`;
  - `<CghOpenButton />` in the title row;
  - `<CghPageSlot loadersFor={cghLoaders} />` before the closing `</div>`;
  - the `cghLoaders` memo;
  - the Esc handler (`cgh.get().pending`, `toMap`, `help`, `leaving`, `open`, the menu / text-field capture, `cgh.escape()`).
  Re-apply them onto `dev`'s version.

### To leave the module out of a build

Remove the button and the page slot in `SvBrowserDialog.tsx`. The new files are then unused, and the other shared
changes are inert: `MOI` was already exported, and the optional prop and loader are unused.
