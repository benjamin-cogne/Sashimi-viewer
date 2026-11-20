# Array-CGH module: what this branch adds, and how to merge it

Branch `dev-WD`. Base: `e0f80a4` (*SV browser: transcript HGVS placed by VariantValidator first, Ensembl as fallback*), which
`dev` contains. The module is one new page, **one level under the SV browser** (SV browser → **Array-CGH ▸**). It adds
files of its own and touches four shared files in a few lines. Nothing else of the viewer changes. With no array-CGH
file loaded, the only visible change is the **Array-CGH ▸** button in the SV browser's header.

## What the page does

- **Quick guide.** *? Quick guide* in the page's bar: six steps, one annotated screenshot each (`CghGuide.tsx`, images in `src/assets/cgh-guide/`, synthetic data, remade by `node docs/cgh-demo/make-guide-shots.cjs`).
- **Input.** Files are read in the page, never uploaded:
  - CytoGenomics interval tables (`.txt`), CSV / TSV of any layout (a column-mapping dialog when the columns are not
    recognised, remembered by column names);
  - `.xlsx` (first sheet), BED, SV VCF (`.vcf`, `.vcf.gz`);
  - zips holding any of these (one patient each).
  - Encodings: UTF-8 with or without BOM, UTF-16, Windows-1252. A binary `.xls` is refused with a clear message.
- **Import question.** Minimum number of probes (3) and minimum size (log slider: any, 10 kb to 5 Mb, or a value typed in
  bp / kb / Mb, also in the filter bar); the 200 largest calls annotated by
  themselves, the others on demand (*annotate*); each patient's calls under the minimums revealed on demand (*+ N more*); sex per patient.
- **Table.** One tab per patient plus *All patients*, filters (series: the number of other patients carrying the call
  compared with = ≥ ≤ > < n), sort by priority.
- **Per call, read on GRCh37** by the SV browser's own routes (UCSC, PanelApp, gnomAD):
  - RefSeq genes; PanelApp ratings and modes of inheritance (asked 25 genes a request);
  - ClinGen HI / TS of genes and of ISCA regions; gnomAD LOEUF;
  - AChro-Puce recurrent CNVs; DGV Gold (3 relations);
  - exons cut and reading frame; terminal / arm / chromosome labels; rearrangement hints;
  - mosaic fraction (AChro-Puce formula; no DLRS input); recurrence in the loaded series; the laboratory's benign
    list (BED, kept in `localStorage`).
- **Priority.** HIGH / MODERATE / LOW from listed rules (`cghScore.ts`). It is not a classification.
- **Classification, note, CNV-Hub link.** **Export classified table** writes each file back with Classification and
  Comment columns (a VCF: INFO `CGH_CLASS` / `CGH_NOTE`; an xlsx: a tab-separated `.classified.txt`).
- **The browser below.** The call clicked is drawn as DEL / DUP in the SV browser's own `SvGenomeBrowser`, with the
  same tracks.

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

Also new: this file, `ARRAY-CGH.md`.

**Shared files changed (whole footprint: 4 files, about 50 lines):**

| File | Change |
|---|---|
| `src/components/sashimi/SvBrowserDialog.tsx` | imports `CghOpenButton`, `CghPageSlot`, `cgh`; the button in the header; the page slot; Esc goes up one level (question / mapping dialog → page → viewer) |
| `src/components/sashimi/SvGenomeBrowser.tsx` | `export const MOI` (was `const`); optional loader field `getDosageRegions` in `SvGeneLoaders`; optional prop `strongEvent` (off by default: shades the event, bolder breakpoints); **bug fix, also visible in the SV browser**: the SVG's ids (`svb-clip`, gradients, pattern) made unique per instance with `useId` (`sid(...)`): the array-CGH page's browser drawn over the SV browser's resolved `url(#svb-clip)` to the other, shorter one, so a grown ClinVar SV track (*▸ all SVs*) hid all below it; bug fix, also visible in the SV browser: in "▸ all SVs", ClinVar SVs that find no free lane among the 60 are no longer drawn over others but counted ("+n SV not drawn (zoom in)") |
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
automatically. The first merge, for the record:

As of `origin/dev` = `6515ae6` (checked in memory with `git merge-tree`, nothing written), a merge gives **two
conflicts**:

1. **`sashimi-viewer.html`.** It is a build output: take either side, then rebuild (step 5).
2. **`src/components/sashimi/SvGenomeBrowser.tsx`.** `dev` changed the same three places:
   - **Props interface.** `dev` added `onGeneSummary?`, this branch `strongEvent?`. Keep **both** fields.
   - **The `MOI` constant.** `dev` moved it to the new `SvGeneSummary.tsx` (`export const MOI` there) and removed it from
     `SvGenomeBrowser.tsx`. Keep `dev`'s side: **drop** this branch's `export const MOI` block. Then, in
     `CghPanel.tsx`, change
     `import { MOI, type SvGeneLoaders } from './SvGenomeBrowser';` to
     `import { type SvGeneLoaders } from './SvGenomeBrowser';` plus `import { MOI } from './SvGeneSummary';`.
   - **Function signature.** Keep both destructured props: `…, variant, onGeneSummary, strongEvent }: Props)`.
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
