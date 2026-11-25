/**
 * The probe calls under the filters: events detected from each patient's probes (cghSegment.ts) that the file does not
 * call, a second opinion listed per patient on demand, drawn like a call when clicked, never exported. Their settings
 * (fewest probes, |mean log ratio|, sensitivity) in one popover.
 */
import { cgh, useCgh, probeCnv } from './cghStore';
import { agrees, SENSITIVITY, SEG_DEFAULTS, type ProbeCalls, type Sensitivity } from './cghSegment';
import { estCopies, estTip } from './cghMosaic';
import { Stepper, Popover } from './controls';
import { fmtBp as fmt, sizeText, type CghCnv, type CghFile, type CghSample } from './cghIntervals';
import type { GenomeBuild } from '../../standalone/ensembl';
import { isProbeOnly } from './cghProbeOnly';

const sensOf = (t: number): Sensitivity => (t <= SENSITIVITY.high ? 'high' : t >= SENSITIVITY.low ? 'low' : 'normal');
const HOW = 'Circular binary segmentation of the probes (Olshen et al. 2004) with a fixed threshold, the segments called by their probes and mean log ratio, adjacent ones of one direction merged; chrY only for a patient said male. A second opinion beside the file\'s calls.';

/** the detected events the file does not call (same type, reciprocal overlap ≥ 50 %) */
export const extrasOf = (r: ProbeCalls, x: CghSample, build: GenomeBuild | null = 'GRCh37') => r.events.filter(e => !x.cnvs.some(c => agrees(e, c, build)));


export function ProbeSettings() {
  const o = useCgh().probeOpts, sens = sensOf(o.t);
  const def = o.minProbes === SEG_DEFAULTS.minProbes && o.minLr === SEG_DEFAULTS.minLr && o.t === SEG_DEFAULTS.t;
  return (
    <Popover width={300} title={HOW} label="Probe calls settings"
      buttonClass="inline-flex items-center gap-1 h-[22px] px-1.5 rounded border border-slate-200 bg-white text-[11px] text-slate-700 hover:border-slate-300"
      button={<>settings{def ? '' : ' •'} ▾</>}>
      <div className="flex flex-col items-start gap-2 p-1 text-xs" data-cgh-probe-settings>
        <span className="[&>span]:!h-[24px]"><Stepper value={o.minProbes} onChange={v => cgh.setProbeOpts({ minProbes: v })} min={1} max={100} label="probes ≥" title="Fewest probes of an event (3 by default: events of 1-2 probes are mostly noise)" /></span>
        <span className="[&>span]:!h-[24px]"><Stepper value={o.minLr} onChange={v => cgh.setProbeOpts({ minLr: v })} min={0.1} max={2} step={0.05} width={44} label="|log ratio| ≥" title="Smallest |mean log ratio| of an event (0.25 by default: a one-copy change in about 35 % of cells)" /></span>
        <label className="inline-flex items-center gap-1 text-slate-600" title="How easily a run of probes is cut off from its neighbours (the threshold on the segmentation statistic: high 6, normal 8, low 10)">sensitivity
          <select value={sens} onChange={e => cgh.setProbeOpts({ t: SENSITIVITY[e.target.value as Sensitivity] })} className="h-[22px] rounded border border-slate-200 bg-white text-[11px] outline-none cursor-pointer">
            {(Object.keys(SENSITIVITY) as Sensitivity[]).map(k => <option key={k} value={k}>{k}</option>)}</select></label>
        {!def && <button className="text-accent-700 hover:underline" onClick={() => cgh.setProbeOpts(SEG_DEFAULTS)}>back to the defaults</button>}
        <span className="text-[10.5px] text-slate-500 leading-4">{HOW}</span>
      </div>
    </Popover>
  );
}

/** the line under the filters: per patient with probes, its probe calls not in the file, listed on demand */
export function ProbeCallsBar({ listed, shown, onShow }: { listed: { x: CghSample; f: CghFile }[]; shown: string | null; onShow: (c: CghCnv, f: CghFile) => void }) {
  const s = useCgh();
  const withProbes = listed.filter(({ x }) => cgh.hasProbes(x.key));
  if (!withProbes.length) return null;
  const many = withProbes.length > 1;
  return (
    <div className="px-5 pb-1.5 text-[11px] text-slate-600" data-cgh-probe-calls>
      <div className="flex flex-wrap items-center gap-x-2 gap-y-1">
        <span className="font-semibold text-slate-800" title={HOW}>Probe calls</span>
        {withProbes.map(({ x, f }) => {
          const err = cgh.probeErrorOf(x.key);
          if (err) return <span key={x.key} className="h-[22px] inline-flex items-center px-1.5 rounded border border-red-200 bg-red-50 text-red-700" data-cgh-probe-error={x.key}
            title="The probe calls could not be computed for these probes and settings; changing a setting or the sex tries again">{many ? `${x.name}: ` : ''}probe calls failed: {err.slice(0, 120)}</span>;
          // a probes-only patient: its probe calls are its table's rows
          if (isProbeOnly(f)) { const r = cgh.probeCallsOf(x.key); return <span key={x.key} className="h-[22px] inline-flex items-center px-1.5 rounded border border-amber-200 bg-amber-50 text-amber-800" data-cgh-probe-toggle={x.key}
            title="No call table: the probe calls are this patient's calls, listed in the table">{many ? `${x.name}: ` : ''}{r ? `${r.events.length} call${r.events.length === 1 ? '' : 's'} = the table's rows` : 'computing…'}</span>; }
          const r = cgh.probeCallsOf(x.key), n = r ? extrasOf(r, x, f.build).length : null, open = !!s.probeOpen[x.key];
          return <button key={x.key} disabled={n == null} onClick={() => cgh.setProbeOpen(x.key, !open)} data-cgh-probe-toggle={x.key}
            className={`h-[22px] px-1.5 rounded border text-[11px] ${open ? 'border-accent-300 bg-accent-50 text-accent-800' : 'border-slate-200 bg-white text-slate-700 hover:border-slate-300'} disabled:text-slate-400`}
            title={r ? `${r.events.length} events found in ${x.name}'s probes (${Math.round(r.ms)} ms; noise ${r.sigma.toFixed(2)}), ${n} not in the file` : 'Being computed from the probes…'}>
            {open ? '▾' : '▸'} {many ? `${x.name}: ` : ''}{n == null ? 'computing…' : `${n} not in the file`}</button>;
        })}
        <ProbeSettings />
        {withProbes.some(({ f }) => !isProbeOnly(f)) && <span className="text-slate-400">· probe calls not in a file are never exported</span>}
      </div>
      {withProbes.filter(({ x, f }) => s.probeOpen[x.key] && !isProbeOnly(f)).map(({ x, f }) => {
        const r = cgh.probeCallsOf(x.key); if (!r) return null;
        const ev = extrasOf(r, x, f.build);
        return <div key={x.key} className="mt-1 border border-slate-200 rounded-lg overflow-auto max-h-[180px]" data-cgh-probe-list={x.key}>
          <table className="w-full text-[11px] border-collapse italic text-slate-500">
            <thead className="sticky top-0 bg-slate-50 not-italic text-left whitespace-nowrap"><tr className="border-b border-slate-200">
              {many && <th className="px-1.5 py-0.5 font-medium">Sample</th>}
              <th className="px-1.5 py-0.5 font-medium">Position</th><th className="px-1.5 py-0.5 font-medium text-right">Size</th><th className="px-1.5 py-0.5 font-medium">Type</th>
              <th className="px-1.5 py-0.5 font-medium text-right">Probes</th><th className="px-1.5 py-0.5 font-medium text-right">Mean LogR</th><th className="px-1.5 py-0.5 font-medium" title="Estimated from the mean log ratio">Copies</th>
              <th className="px-1.5 py-0.5 font-medium w-full" /></tr></thead>
            <tbody>
              {ev.map(e => {
                const c = probeCnv(e, x.key, x.name), est = estCopies(c, s.sex[x.key], f.build), on = shown === c.id;
                return <tr key={c.id} data-cgh-probe-row onClick={() => { cgh.setShown(c.id); onShow(c, f); }} title="From the probes, not in the file. Click: draw it below. Not exported."
                  className={`border-b border-slate-100 cursor-pointer ${on ? 'bg-accent-50 text-slate-700' : 'hover:bg-slate-50'}`}>
                  {many && <td className="px-1.5 py-0.5 whitespace-nowrap">{x.name}</td>}
                  <td className="px-1.5 py-0.5 font-mono text-[10.5px] whitespace-nowrap">{c.chrom}:{fmt(c.start)}-{fmt(c.end)}</td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums whitespace-nowrap">{sizeText(c.sizeKb)}</td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap">{e.kind === 'loss' ? 'del' : 'dup'}</td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums">{e.n}</td>
                  <td className="px-1.5 py-0.5 text-right tabular-nums">{e.mean.toFixed(2)}</td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap" title={est ? estTip(est, e.mean) : undefined}>{est?.text ?? ''}</td>
                  <td className="px-1.5 py-0.5 whitespace-nowrap">from probes</td>
                </tr>;
              })}
              {!ev.length && <tr><td colSpan={8} className="px-2 py-1">Every event found in the probes is in the file.</td></tr>}
            </tbody>
          </table>
        </div>;
      })}
    </div>
  );
}

/** the methods of what is read from the probes, in the page's ⓘ panel (written from cghProbes / cghSegment / cghMosaic / CghGenome) */
export function ProbeMethods({ probes }: { probes: boolean }) {
  return <div className="mt-1.5" data-cgh-probe-methods>
    <div className="font-semibold text-slate-700">Methods: copy number, ISCN{probes ? ', and what is read from the probes (CGH_Log_Ratio.txt)' : ''}</div>
    {!probes && <div className="italic">No probe file loaded: the probe plots, DLRS and probe calls below apply only when a CGH_Log_Ratio.txt is added (with a patient's table, or alone: a probes-only patient).</div>}
    <div><b>DLRS</b> (derivative log ratio spread), per patient: the standard deviation of the differences between the log ratios of consecutive probes (sorted by position, within each chromosome, every chromosome of the file including X and Y, chrM left out), divided by √2; probes without a value skipped; no outlier removed (events and waves raise it a little). Under 10 differences: not given. Words from Agilent's guidance: below 0.2 excellent, 0.2 to 0.3 good, above 0.3 poor.</div>
    <div><b>Probe calls</b>, per patient and chromosome, probes without a value skipped: circular binary segmentation (Olshen AB et al., Biostatistics 2004;5:557-72). The arc whose mean differs most from the rest of the stretch, Z = (S<sub>j</sub> − S<sub>i</sub> − k·T/n)·√(n / (k(n − k))) / σ, cuts it when |Z| ≥ t, and each piece is cut again (a stretch of more than 10,000 probes, as on 180k or 1M arrays: the arcs between 2,500 blocks first, then exactly around the best one, at the stretch's ends and every arc of up to two blocks, so a short event inside a block is not missed); a fixed threshold t, without the permutation test (sensitivity high / normal / low: t = 6 / 8 / 10). σ, the noise per probe, taken once per patient: the median absolute difference between consecutive probes × 1.4826 / √2 (robust to the events). A segment is a gain (loss) when its mean log ratio is ≥ +threshold (≤ −threshold; 0.25 by default); a segment of fewer than the minimum probes (3 by default) lying between two of one direction takes that direction; neighbouring segments of one direction are merged (mean weighted by probes); kept when the merged mean still passes the threshold and it has at least the minimum probes. chrY: its pseudo-autosomal probes left out (the chrX PAR ones again), the rest read only for a patient said male (cross-hybridisation in a female). A probe call is <i>not in the file</i> when no file call of the same chromosome and type overlaps it reciprocally by ≥ 50 % (the overlap over the longer of the two; a chrY PAR call compared as its chrX twin); those are listed, never exported. The file's calls carry no mark (they are made from these same probes).</div>
    <div><b>Probes only</b>: a CGH_Log_Ratio.txt added without its patient's table (alone, or in a zip without one; a zip whose table waits for its columns keeps its probes for that table) makes a patient of its own, named from its <i>Normalized log ratio(name)</i> header (else the file name), its build from the AMADID line. Its calls are these probe calls (current settings and sex, recomputed when they change; a class or note stays on a call whose position is unchanged): Probes and Log R the segment's, copies, mosaic and ISCN as above, band from the UCSC bands, then the usual annotation and priority. They are not CytoGenomics calls: to confirm in the dedicated software. Its calls are drawn once in the whole genome (as its calls, not again as probe calls). Exported as a new tab-separated file, <i>name</i>.probe-calls.txt (source file, build, settings, date; the UCSC bands read for every call first, a header line counting rows left without one; Chromosome, Start, Stop, Cytoband, Size(kb), Type, #Probes, Mean Log Ratio, Copies, ISCN, Classification, Comment). A table of the same patient added later takes the probes and replaces the probe-only patient when the question is accepted (its classes and notes are not carried over: the question says so).</div>
    <div><b>Copy number</b>, when the file gives none: the n whose log2(n / 2) is nearest the call's mean log ratio (x0 taken as −2, so a loss is x0 below −1.5, else x1; a gain x3 up to +0.79, x4 up to +1.16…; a male's X outside the PARs and chrY: log2(n / 1)). <b>Mosaic</b>: for a one-copy change, the fraction of cells from the AChro-Puce calculator (loss 2 − 2·2<sup>LR</sup>, gain 2·2<sup>LR</sup> − 2; one-copy X / Y of a male 1 − 2<sup>LR</sup> and 2<sup>LR</sup> − 1); below 80 % the call is shown <i>mosaic ~n%</i> (an autosomal loss above −0.74, a gain below +0.49), else as a whole x1 / x3.</div>
    <div><b>ISCN</b> (ISCN 2024, chapter 8, microarray): <span className="font-mono">arr[build] chromosome, bands from pter to qter (one band once), (start_end) with thousands commas, x copies</span>, the mosaic fraction in brackets (<span className="font-mono">x1[0.48]</span>), a whole chromosome as <span className="font-mono">(21)x3</span>; the band is read from the UCSC cytogenetic bands at sub-band precision, the file's band only when they cannot be read, the copies the file's, else those above. The file's own ISCN is shown beside ours when they differ. The patient's ISCN joins the calls of the chosen classes by commas, chromosomes 1-22, X, Y, then position; with none, <i>no abnormal call selected</i> (no normal-result form written). Inheritance is not written.</div>
    <div className="text-amber-800"><b>Disclaimer.</b> Copy numbers are inferred from the log ratio with the thresholds above. All numbers, calls and nomenclature must be checked in the dedicated analysis software; the authors cannot be held responsible for any misinterpretation of these data.</div>
    <div><b>Whole genome view</b>: each probe drawn in its pixel cell (one mark per cell and colour, so 61k probes cost no more than the pixels), chromosomes side by side (horizontal) or as columns, position going down from pter and log ratio across (vertical); the dark line is the median of each pixel line. A chromosome opened is drawn as the call's probe panel (every probe, moving median of 3 / 7 / 15). Log-ratio range ±0.5 / ±1 (default) / ±2, further with ⊖ (buttons left of the log-ratio axis), shared by the three views; <i>show probes</i> off draws only the median line, segments and calls; every probe with a value is drawn (each colour present in a pixel cell marked); a call's probes' mean is over the probes whose midpoint lies within it, ends included; a probe beyond it is drawn as a triangle at the edge; the colours stay at ±0.2. A male's X and Y may lie near −1 (one copy against a reference of two): not normalised.</div>
  </div>;
}
