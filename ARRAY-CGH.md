# Array-CGH module: what this branch adds, and how to merge it

Branch `dev-WD`. Base: `e0f80a4` (*SV browser: transcript HGVS placed by VariantValidator first, Ensembl as fallback*), which
`dev` contains. The module is one new page, **one level under the SV browser** (SV browser → **Array-CGH ▸**). It adds
files of its own and touches four shared files in a few lines. Nothing else of the viewer changes. With no array-CGH
file loaded, the only visible change is the **Array-CGH ▸** button in the SV browser's header.

## What the page does

- **Examples.** *Examples ▾* (empty page and frame bar, as the SV browser offers its examples) adds three invented patients (`cghExamples.ts`, CytoGenomics text, GRCh37) through the normal import question: 22q11.2 LCR22A-D deletion (HIGH), 1q21.1 distal BP3-BP4 duplication (AChro-Puce PIEV, ClinGen TS 2 region: MODERATE), terminal 16q24.3 loss breaking in ANKRD11 intron 4 + terminal 9p24.3 gain (der() hint, ANKRD11 5′ part cut, HIGH). Coordinates checked on UCSC hg19 (ncbiRefSeqCurated, ClinGen tracks). Public build: never real data there.
- **Quick guide.** *? Quick guide* in the page's bar: an interactive tour of 15 steps on the examples (`CghGuide.tsx`): an outline on the target and a popover beside it (kept on screen, following scroll, never over the target when it fits elsewhere, draggable by its title); each step does its own action (load the examples, Analyse, open More filters or a call): the user only clicks *Next*; Esc / *Skip tour* ends it (the store's `help` keeps Esc from closing the page). Targets found by `data-cgh-*` attributes. The former screenshot guide (`src/assets/cgh-guide/`, `docs/cgh-demo/`) is removed.
- **Input.** Files are read in the page, never uploaded:
  - CytoGenomics interval tables (`.txt`), CSV / TSV of any layout (a column-mapping dialog when the columns are not
    recognised, remembered by column names, only when the header row looks like one);
  - `.xlsx` (first sheet), BED, SV VCF (`.vcf`, `.vcf.gz`);
  - the mapping dialog also asks, besides chromosome / start / end (or one region column), for a type column, OR one or
    several log-ratio columns (`MapSpec.lrMore`: gain / loss columns merged as the non-zero one, typed by a column name
    saying gain / loss, else the sign), OR a copy-number column, and names what is missing; with no type: copy number
    (< 2 loss, > 2 gain; 1 normal on a male's X and on Y), else the log ratio's sign; rows with neither skipped, counted
    in the question; specs saved before `lrMore` load as before;
  - CytoGenomics' automated interval-based report (`IntervalBasedReport.xls`, tab-separated text despite its name, or `.xlsx`), recognised by its header (`AberrationNo`, `Chr`, `Start`, `Stop`, `#Probes`, `Amplification` / `Gain` / `Loss` / `Deletion`) and read without the mapping dialog (`parseAberrations` in `cghIntervals.ts`): one log ratio from the four columns (the non-zero one, the largest |value| if several), the type from the column (else the sign), `pval` (`CghCnv.pval`), genes, the build from `Genome:`, the sample from its name line (else `Array Design Info`); copies null (no ISCN written, mosaic from the log ratio); exported as it came plus Classification / Comment columns;
  - zips holding any of these (one patient each): tables found by content, the other files of an automated export left aside and listed in one note (*ignored: n other files*); a table needing the mapping is asked only when the zip holds no table read by itself;
  - probes' log ratios (`CGH_Log_Ratio.txt`, `cghProbes.ts`): kept as bytes with the zip's patient (alone: with the loaded sample of the `Normalized log ratio(<name>)` header, else a note *not matched*), read on first `cgh.probes(sampleKey)` into `Uint32Array` start / end and `Float32Array` log ratios per chromosome (NaN kept) and probe names; `cgh.hasProbes`; `dlrs` (SD of consecutive differences / √2) in the import question and the call card. **Log ratio (probes)** panel (`CghProbePanel.tsx`) in the call's card above the browser: a canvas of the window's probes (binary search on the sorted starts), colours black / red / green at ±0.2, moving median 3 / 7 / 15, call shaded with its mean, ±0.2 dashed, hover tip; its own zoom / pan (the browser's window is internal to `SvGenomeBrowser`: not synced, that file untouched).
  - **Probe calls** (`cghSegment.ts`, pure; `CghProbeCalls.tsx`): events detected from a patient's probes, a second opinion. Per chromosome, NaN skipped, circular binary segmentation (Olshen et al. 2004) with a fixed threshold `t` on the arc statistic, σ = MAD of consecutive differences / √2 taken once per sample (no permutations: O(n²) per chromosome, about half a second for 61k probes, one patient per tick via `cgh.runProbeCalls`, no worker); segments called with ≥ `minProbes` (3) probes and |mean| ≥ `minLr` (0.25), a < `minProbes` segment between two of one direction absorbed, adjacent same-direction merged; chrY PAR left out (the chrX PAR probes again), chrY read for a male only (a female's gives many false events). Settings in `state.probeOpts` (localStorage `sashimi.cghProbeCalls`), sensitivity high / normal / low = t 6 / 8 / 10. File calls marked *probes: seen / not seen* (`agrees`: same type, reciprocal overlap ≥ 50 %, chrY PAR mapped to chrX); events not in the file listed per patient (*▸ n not in the file*), drawn as a `CghCnv` with `fromProbes` (`probeCnv`), never exported; segments a purple line in the probe panel. Defaults calibrated locally on real arrays against CytoGenomics ADM-2 calls (data not in the repo): at t 8 every file call of ≥ 3 probes found, 1-2 extras per array, no 1-2-probe call reproduced; t 10 lost calls of 3-5 probes, t 6 tripled the extras. Invented recall (`scripts/cgh-probe-fixture.mjs`: del -0.9 × 5, dup +0.55 × 200, mosaic -0.4 × 30, male X dup +0.9 × 30): 4/4, no extra.
  - Encodings: UTF-8 with or without BOM, UTF-16, Windows-1252. A binary `.xls` is refused with a clear message.
- **Import question.** Minimum number of probes (3) and minimum size (log slider: any, 10 kb to 5 Mb, or a value typed in
  bp / kb / Mb, also in the filter bar); the 200 largest calls annotated by
  themselves, the others on demand (*annotate*); each patient's calls under the minimums revealed on demand (*+ N more*); sex per patient.
- **Table.** One tab per patient plus *All patients*, filters (series: the number of other patients carrying the call
  compared with = ≥ ≤ > < n; a sample name loaded twice is one patient, said at import), gene search (file genes and
  RefSeq genes read), sort by priority.
- **Per call, read on GRCh37** by the SV browser's own routes (UCSC, PanelApp, gnomAD):
  - RefSeq genes; PanelApp ratings and modes of inheritance (asked 25 genes a request);
  - ClinGen HI / TS of genes and of ISCA regions; gnomAD LOEUF;
  - AChro-Puce recurrent CNVs; DGV Gold (3 relations);
  - exons cut and reading frame; terminal / arm / chromosome labels; rearrangement hints;
  - mosaic fraction (AChro-Puce formula; no DLRS input); an estimated copy number from the log ratio when the file gives none (`estCopies` in `cghMosaic.ts`: nearest log2(n/2), x0 at −2; male X / Y normal 1), shown as *x1?* / *estimated x1 (or mosaic ~n%)*, never put in an ISCN; recurrence in the loaded series; the laboratory's benign
    list (BED, kept in `localStorage`).
- **Priority.** HIGH / MODERATE / LOW from listed rules (`cghScore.ts`). It is not a classification. Per Riggs et al. 2020: a ClinGen HI 3 / TS 3 region is HIGH only wholly in the call (2A, share ≥ 0.999), MODERATE at 80-99 % (2B: check the critical region); an HI 3 gene losing only UTR exons (2C-2 / 2D-1) or only its last exon (2D-2 / 2D-3) is MODERATE (`sparing()`, from `CnvExons.cds` and the exon ranks; genes unknown to RefSeq stay HIGH).
- **Retry.** Failures kept per call and source (`CnvAnn.failed`; PanelApp per gene): *↻ retry (n)* in the column header (per patient; a menu in All patients) and a *retry* per source in the detail re-ask only what failed (`retry()` in `cghAnnotate.ts`); an unexpected error ends the call on what was read.
- **Classification, note, CNV-Hub link.** **Export classified table** writes each file back with Classification and
  Comment columns (a VCF: INFO `CGH_CLASS` / `CGH_NOTE`; an xlsx: a tab-separated `.classified.txt`).
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
| `cghSeries.ts` | recurrence across loaded patients |
| `cghTopo.ts` | terminal / arm / chromosome labels, rearrangement hints |
| `cghLabBenign.ts` | the laboratory's benign list |
| `cghExons.ts` | exons cut and reading frame |
| `cghExamples.ts` | built-in invented example patients |
| `cghProbes.ts` | probes' log ratios (`CGH_Log_Ratio.txt`): detection, lazy reading into typed arrays, DLRS, window search |
| `CghProbePanel.tsx` | the *Log ratio (probes)* panel above the browser |
| `cghSegment.ts` | probe calls: segmentation (CBS, fixed threshold), calling, agreement with the file's calls |
| `CghProbeCalls.tsx` | the *Probe calls* line: settings, per-patient lists, the *probes: seen / not seen* mark |
| `CghGuide.tsx` | the interactive tour (*? Quick guide*) |
| `CghClinGenTrack.tsx` | the *ClinGen CNV regions* track (`clinGenTrack(build)`) and the detail's recurrent-CNV lines |
| `cghRecurrent.ts` | ClinGen regions on a chromosome, AChro-Puce PIEV tags, a call's recurrent matches and LCR ends |
| `cghClinGenData.ts` | generated data: HI 3 / TS 3 ISCA regions, recurrent CNV list v2.0 (both builds) |

Also `scripts/cgh-probe-fixture.mjs` writes an invented zip (report + probes, known events) to test the probe calls; `scripts/cgh-clingen-data.mjs` writes `cghClinGenData.ts` from the downloaded UCSC JSON and ClinGen BED files.

Also new: this file, `ARRAY-CGH.md`.

**Shared files changed (whole footprint: 4 files, about 65 lines):**

| File | Change |
|---|---|
| `src/components/sashimi/SvBrowserDialog.tsx` | imports `CghOpenButton`, `CghPageSlot`, `cgh`; the button in the header; the page slot; Esc goes up one level (question / mapping dialog → page → viewer) |
| `src/components/sashimi/SvGenomeBrowser.tsx` | `export const MOI` (was `const`); optional loader field `getDosageRegions` in `SvGeneLoaders`; optional prop `strongEvent` (off by default: shades the event, bolder breakpoints); optional prop `extraTrack` + exported `SvExtraTrack` interface (a caller's track under AChro-Puce: state `showExtra`, `extraLaid` / `yExtra` in the layout before `yRefBand`, its side label, its `<g>` before the reference band, a `Pill` after AChro-Puce's and its legend; absent in the SV browser, so nothing drawn there); **bug fix, also visible in the SV browser**: the SVG's ids (`svb-clip`, gradients, pattern) made unique per instance with `useId` (`sid(...)`): the array-CGH page's browser drawn over the SV browser's resolved `url(#svb-clip)` to the other, shorter one, so a grown ClinVar SV track (*▸ all SVs*) hid all below it; bug fix, also visible in the SV browser: in "▸ all SVs", ClinVar SVs that find no free lane among the 60 are no longer drawn over others but counted ("+n SV not drawn (zoom in)") |
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
   - Click **Array-CGH ▸**, add an interval table, analyse it, click a call: it is drawn below, and Esc goes back one
     level at a time.
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
  - the Esc handler (`cgh.get().pending`, `toMap`, `open`).
  Re-apply them onto `dev`'s version.

### To leave the module out of a build

Remove the button and the page slot in `SvBrowserDialog.tsx`. The new files are then unused, and the other shared
changes are inert: `MOI` was already exported, and the optional prop and loader are unused.
