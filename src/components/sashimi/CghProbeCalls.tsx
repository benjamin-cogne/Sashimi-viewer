/**
 * The probe calls under the filters: events detected from each patient's probes (cghSegment.ts) that the file does not
 * call, a second opinion listed per patient on demand, drawn like a call when clicked, never exported. Their settings
 * (fewest probes, |mean log ratio|, sensitivity) in one popover.
 */
import { cgh, useCgh, probeCnv } from './cghStore';
import { agrees, SENSITIVITY, SEG_DEFAULTS, type ProbeCalls, type Sensitivity } from './cghSegment';
import { estCopies, estTip } from './cghMosaic';
import { Stepper, Popover } from './controls';
import type { CghCnv, CghFile, CghSample } from './cghIntervals';

const fmt = (n: number) => Math.round(n).toLocaleString('en-US');
const sizeText = (kb: number) => (kb >= 1000 ? `${(kb / 1000).toFixed(2)} Mb` : kb >= 1 ? `${kb.toFixed(kb >= 100 ? 0 : 1)} kb` : `${Math.round(kb * 1000)} bp`);
const sensOf = (t: number): Sensitivity => (t <= SENSITIVITY.high ? 'high' : t >= SENSITIVITY.low ? 'low' : 'normal');
const HOW = 'Circular binary segmentation of the probes (Olshen et al. 2004) with a fixed threshold, the segments called by their probes and mean log ratio, adjacent ones of one direction merged; chrY only for a patient said male. A second opinion beside the file\'s calls.';

/** the detected events the file does not call (same type, reciprocal overlap ≥ 50 %) */
export const extrasOf = (r: ProbeCalls, x: CghSample) => r.events.filter(e => !x.cnvs.some(c => agrees(e, c)));

/** a file call's mark: whether the probe calls see it (only once they are computed) */
export function ProbeSeenTag({ c }: { c: CghCnv }) {
  const r = cgh.probeCallsOf(c.sampleKey);
  if (!r) return null;
  const seen = r.events.some(e => agrees(e, c));
  return <span data-cgh-seen={seen ? 'seen' : 'not'} className={`ml-1 text-[9.5px] ${seen ? 'text-slate-400' : 'text-amber-700'}`}
    title={seen ? 'The probe calls find this call too (same type, reciprocal overlap ≥ 50 %)' : 'The probe calls do not find this call with their settings (often a call of 1-2 probes)'}>probes: {seen ? 'seen' : 'not seen'}</span>;
}

function Settings() {
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
        {withProbes.map(({ x }) => {
          const r = cgh.probeCallsOf(x.key), n = r ? extrasOf(r, x).length : null, open = !!s.probeOpen[x.key];
          return <button key={x.key} disabled={n == null} onClick={() => cgh.setProbeOpen(x.key, !open)} data-cgh-probe-toggle={x.key}
            className={`h-[22px] px-1.5 rounded border text-[11px] ${open ? 'border-accent-300 bg-accent-50 text-accent-800' : 'border-slate-200 bg-white text-slate-700 hover:border-slate-300'} disabled:text-slate-400`}
            title={r ? `${r.events.length} events found in ${x.name}'s probes (${Math.round(r.ms)} ms; noise ${r.sigma.toFixed(2)}), ${n} not in the file` : 'Being computed from the probes…'}>
            {open ? '▾' : '▸'} {many ? `${x.name}: ` : ''}{n == null ? 'computing…' : `${n} not in the file`}</button>;
        })}
        <Settings />
        <span className="text-slate-400">· file calls marked <i>probes: seen / not seen</i> · probe calls are never exported</span>
      </div>
      {withProbes.filter(({ x }) => s.probeOpen[x.key]).map(({ x, f }) => {
        const r = cgh.probeCallsOf(x.key); if (!r) return null;
        const ev = extrasOf(r, x);
        return <div key={x.key} className="mt-1 border border-slate-200 rounded-lg overflow-auto max-h-[180px]" data-cgh-probe-list={x.key}>
          <table className="w-full text-[11px] border-collapse italic text-slate-500">
            <thead className="sticky top-0 bg-slate-50 not-italic text-left whitespace-nowrap"><tr className="border-b border-slate-200">
              {many && <th className="px-1.5 py-0.5 font-medium">Sample</th>}
              <th className="px-1.5 py-0.5 font-medium">Position</th><th className="px-1.5 py-0.5 font-medium text-right">Size</th><th className="px-1.5 py-0.5 font-medium">Type</th>
              <th className="px-1.5 py-0.5 font-medium text-right">Probes</th><th className="px-1.5 py-0.5 font-medium text-right">Mean LogR</th><th className="px-1.5 py-0.5 font-medium" title="Estimated from the mean log ratio">Copies</th>
              <th className="px-1.5 py-0.5 font-medium w-full" /></tr></thead>
            <tbody>
              {ev.map(e => {
                const c = probeCnv(e, x.key, x.name), est = estCopies(c, s.sex[x.key]), on = shown === c.id;
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
