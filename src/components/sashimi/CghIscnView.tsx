/** The array ISCN shown: a call's (call card, detail, table tooltip) and a patient's line, each with a copy button. */
import { useState } from 'react';
import { cgh, useCgh, classOf, type CghState } from './cghStore';
import type { CghCnv, CghSample } from './cghIntervals';
import { callIscn, patientIscn, iscnBandOf, type CallIscn } from './cghIscn';

const ISCN_OFFERED = ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'PIEV', 'Likely benign', 'Benign'];
const SHORT: Record<string, string> = { Pathogenic: 'P', 'Likely pathogenic': 'LP', 'Uncertain significance': 'VUS', PIEV: 'PIEV', 'Likely benign': 'LB', Benign: 'B' };

/** a call's ISCN from the state (band: the cytogenetic bands read, at sub-band precision, else the file's; whole chromosome from the topology) */
export const iscnOfCall = (s: CghState, c: CghCnv, build: string): CallIscn =>
  callIscn(c, build, iscnBandOf(c, s.ann[c.id]?.band), s.sex[c.sampleKey], s.ann[c.id]?.topo?.span === 'chr');

/** the tooltip lines of a call's ISCN (ours, the file's when it differs) */
export const iscnTip = (i: CallIscn) => [i.text ?? i.file ?? '', i.differs ? `file: ${i.file} (differs)` : ''].filter(Boolean).join('\n');

export function CopyButton({ text, what }: { text: string; what: string }) {
  const [done, setDone] = useState(false);
  return <button className="ml-1 px-1 rounded border border-slate-200 bg-white text-[10px] font-sans font-medium text-slate-600 hover:border-accent-300 hover:bg-accent-50 align-middle" data-cgh-copy
    title={`Copy ${what}`} onClick={e => { e.stopPropagation(); void navigator.clipboard?.writeText(text).then(() => { setDone(true); setTimeout(() => setDone(false), 1200); }, () => {}); }}>{done ? 'copied' : 'copy'}</button>;
}

/** a call's ISCN with its copy button; the file's own ISCN beside when it differs (flagged) */
export function CallIscnText({ i }: { i: CallIscn }) {
  if (!i.text && !i.file) return null;
  return <span className="font-mono text-[12px] text-slate-800" data-cgh-iscn>
    {i.text ?? i.file}<CopyButton text={(i.text ?? i.file)!} what="the ISCN" />
    {i.differs && <span className="ml-1.5 font-sans text-[10.5px] text-amber-700" title="The file's own ISCN differs from the one written here (positions, bands, copies or mosaic fraction)" data-cgh-iscn-differs>
      file: <span className="font-mono">{i.file}</span> (differs)</span>}
  </span>;
}

/** a patient's ISCN from its calls of the chosen classes (a setting), ordered by chromosome then position */
export function PatientIscn({ x, build }: { x: CghSample; build: string }) {
  const s = useCgh(), want = s.view.iscnClasses;
  const calls = x.cnvs.filter(c => want.includes(classOf(s, c).trim())).map(c => ({ c, part: iscnOfCall(s, c, build).part }));
  const text = patientIscn(build, calls), missing = calls.filter(y => !y.part).length;
  const flip = (v: string) => cgh.setView({ iscnClasses: want.includes(v) ? want.filter(w => w !== v) : [...want, v] });
  return <div className="mx-5 mb-1.5 flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] text-slate-600" data-cgh-patient-iscn>
    <span className="font-semibold text-slate-900 text-xs">Patient ISCN</span>
    {text ? <span className="font-mono text-[12px] text-slate-800">{text}<CopyButton text={text} what="the patient's ISCN" /></span>
      : <span className="text-slate-500">no abnormal call selected</span>}
    {missing > 0 && <span className="text-amber-700" title="Calls of these classes without a band or a copy number">{missing} call{missing > 1 ? 's' : ''} left out (no band or copy number)</span>}
    <span className="inline-flex items-center gap-0.5" title="The classes whose calls are written in the patient's ISCN (kept as a setting)">from
      {ISCN_OFFERED.map(v => <button key={v} onClick={() => flip(v)} aria-pressed={want.includes(v)} title={v}
        className={`h-[20px] px-1 rounded border text-[10px] ${want.includes(v) ? 'border-accent-300 bg-accent-50 text-accent-700 font-semibold' : 'border-slate-200 bg-white text-slate-500'}`}>{SHORT[v]}</button>)}</span>
  </div>;
}
