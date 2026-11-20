/**
 * The array-CGH page's quick guide: an interactive tour on the built-in examples (cghExamples.ts). Each step outlines
 * its target and says what to do; a step advances by itself when the user does it, Next does it for them. The page
 * stays usable under it. Esc or Skip ends it; the store's `help` keeps Esc from closing the page too.
 */
import { useEffect, useLayoutEffect, useRef, useState, type ReactNode } from 'react';
import { cgh, useCgh } from './cghStore';
import { EXAMPLE_NAMES, loadExamples } from './cghExamples';

const B = ({ children }: { children: ReactNode }) => <b className="font-semibold text-slate-900">{children}</b>;
const q = <T extends Element = HTMLElement>(sel: string) => document.querySelector<T>(sel);
type Kind = 'loss' | 'gain';
/** an example's call (sample, chromosome, type), when loaded */
const callOf = (sample: string, chrom: string, kind: Kind) => cgh.get().files.flatMap(f => f.samples).find(x => x.name === sample)?.cnvs.find(c => c.chrom === chrom && c.kind === kind);
const loaded = () => EXAMPLE_NAMES.every(n => cgh.get().files.some(f => f.samples.some(x => x.name === n)));
const pendingExamples = () => !!cgh.get().pending?.files.some(f => f.samples.some(x => EXAMPLE_NAMES.includes(x.name)));
/** the call's row; All patients opened when another patient's tab hides it */
const rowOf = (sample: string, chrom: string, kind: Kind) => {
  const c = callOf(sample, chrom, kind), r = c ? q(`tr[data-cgh-call="${c.id}"]`) : null;
  if (c && !r && cgh.get().sample && cgh.get().sample !== c.sampleKey) cgh.setSample(null);
  return r;
};
const shownIs = (sample: string, chrom: string, kind: Kind) => { const c = callOf(sample, chrom, kind); return !!c && cgh.get().shown === c.id; };
/** click the call's row unless it is the one shown (a second click folds it) */
const open = (sample: string, chrom: string, kind: Kind) => { if (!shownIs(sample, chrom, kind)) rowOf(sample, chrom, kind)?.click(); };
const shownRow = () => q('tr[data-cgh-row][data-shown]');

interface Step {
  title: string; body: ReactNode;
  target: () => Element | null;
  /** the user did what the step asks: the tour moves on by itself */
  done?: () => boolean;
  /** Next does it for them */
  act?: () => void | Promise<void>;
}
const DEL = 'EXAMPLE-DEL22Q11', DER = 'EXAMPLE-DER', PIEV = 'EXAMPLE-PIEV-DUP';
const STEPS: Step[] = [
  { title: 'Add samples, or the examples', target: () => q('button[aria-label="Examples"]'),
    body: <>Your files go in with <B>Add samples</B> (or dropped on the page): interval tables, BED, VCF, zips; read here, never sent. For this tour, open <B>Examples ▾</B> and pick <B>All three</B> (invented patients).</>,
    done: () => pendingExamples() || loaded(), act: () => { if (!pendingExamples() && !loaded()) return loadExamples('all'); } },
  { title: 'The import question', target: () => q('[data-cgh-question]'),
    body: <>Calls with fewer <B>probes</B> than the minimum (3), or under the <B>minimum size</B>, are not analysed. Check each patient's <B>sex</B> (read from the file here). Then click <B>Analyse</B>.</>,
    done: () => !cgh.get().pending && loaded(), act: () => q('[data-cgh-analyse]')?.click() },
  { title: 'Patients', target: () => q('[role="tablist"][aria-label="Patients"]'),
    body: <>One tab per patient, and <B>All patients</B> in one table with a Sample column. <B>×</B> on a tab takes the patient off (the examples too).</> },
  { title: 'Filters', target: () => q('[data-cgh-filters]'),
    body: <><B>probes ≥</B>, <B>losses / gains</B>, <B>size</B>, <B>class</B>; search a gene or a band (16p11); order by <B>priority</B> to see HIGH first.</> },
  { title: 'More filters', target: () => q('button[aria-label="More filters"]'),
    body: <>Open <B>More filters</B>: <B>series</B> (how many other loaded patients carry the call), <B>lab list</B>, <B>hide suppressed</B>.</>,
    done: () => !!q('button[aria-label="More filters"][aria-expanded="true"]') },
  { title: 'Priority', target: () => rowOf(DEL, 'chr22', 'loss')?.querySelector('[data-cgh-priority]') ?? null,
    body: <>Each call gets <B>HIGH</B>, <B>MODERATE</B> or <B>LOW</B>: an order for the review, not a class. Click the <B>22q11.2</B> deletion of {DEL} to see why.</>,
    done: () => shownIs(DEL, 'chr22', 'loss'), act: () => open(DEL, 'chr22', 'loss') },
  // the detail itself while its genes are still read
  { title: 'Why this priority', target: () => q('[data-cgh-reasons]') ?? shownRow()?.nextElementSibling ?? null,
    body: <>The rules the call meets, after the ACMG / ClinGen CNV standards (here a whole ClinGen HI 3 region, 2A). Every rule is listed under the table, in <B>ⓘ How the calls are read…</B></> },
  { title: 'Interpret in the browser', target: () => q('[data-cgh-browser]'),
    body: <>The call is drawn below in the SV browser's own browser: genes, <B>DGV Gold</B>, <B>ClinVar</B>, <B>AChro-Puce</B>, <B>PanelApp</B>. Read the call there before classifying it.</> },
  { title: 'A possible der()', target: () => rowOf(DER, 'chr16', 'loss')?.querySelector('[data-cgh-topo]') ?? rowOf(DER, 'chr16', 'loss'),
    body: <>{DER} has a terminal 16q loss and a terminal 9p gain: <B>ter ⚠</B> hints a derivative of an unbalanced translocation, to check by parental karyotype / FISH. Click the <B>16q24.3</B> deletion.</>,
    done: () => shownIs(DER, 'chr16', 'loss'), act: () => open(DER, 'chr16', 'loss') },
  { title: 'A gene cut', target: () => q('[data-cgh-gene="ANKRD11"]') ?? shownRow()?.nextElementSibling ?? null,
    body: <>The deletion breaks inside <B>ANKRD11</B> (HI 3, KBG syndrome): the exons in the call and the effect on the frame are listed, read from the call's ends (theoretical: the true breakpoint lies between probes).</> },
  { title: 'A susceptibility CNV', target: () => rowOf(PIEV, 'chr1', 'gain'),
    body: <>{PIEV}: a 1q21.1 distal duplication. Click it.</>,
    done: () => shownIs(PIEV, 'chr1', 'gain'), act: () => open(PIEV, 'chr1', 'gain') },
  { title: 'PIEV', target: () => q('[data-cgh-achro]'),
    body: <>AChro-Puce classes this duplication <B>PIEV</B> (incomplete penetrance, variable expressivity): <B>MODERATE</B>, with the ClinGen TS 2 region. Parental segregation is advised.</> },
  { title: 'The lab benign list', target: () => q('[data-cgh-lab]'),
    body: <><B>Lab benign CNVs (BED)…</B> loads the laboratory's benign list, kept in this browser. Calls sharing 80 % with an entry are tagged, with a <B>→ Benign?</B> suggestion.</> },
  { title: 'Retry', target: () => q('[data-cgh-thead]'),
    body: <>When a source (UCSC, PanelApp, gnomAD, DGV) does not answer, <B>↻ retry (n)</B> appears in its column header, and <B>retry</B> in the call's detail: only what failed is asked again.</> },
  { title: 'Classify and note', target: () => shownRow()?.querySelector('[data-cgh-class]') ?? null,
    body: <>Choose the class in <B>Classification</B>; <B>✎</B> (or the box in the detail) adds a note. Both stay in this tab only.</> },
  { title: 'Export', target: () => q('[data-cgh-export-btn]'),
    body: <><B>Export classified table</B> writes each file back with Classification and Comment columns. The examples can then be closed with <B>×</B> on their tabs.</> },
];

/** "? Quick guide" in the page's bar: starts the tour */
export function CghGuideButton() {
  return <button onClick={() => cgh.setHelp(true)} className="ml-auto inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="An interactive tour of this page on invented example patients">? Quick guide</button>;
}

type Box = { x: number; y: number; w: number; h: number };
const PAD = 6, GAP = 12, W = 340;
/** where the popover goes: below the target, else above, right, left; kept on screen */
function place(t: Box | null, h: number): { left: number; top: number } {
  const vw = window.innerWidth, vh = window.innerHeight, clampX = (x: number) => Math.max(8, Math.min(vw - W - 8, x)), clampY = (y: number) => Math.max(56, Math.min(vh - h - 8, y));
  if (!t) return { left: (vw - W) / 2, top: Math.max(56, (vh - h) / 2) };
  if (t.y + t.h + GAP + h <= vh - 8) return { left: clampX(t.x), top: t.y + t.h + GAP };
  if (t.y - GAP - h >= 56) return { left: clampX(t.x), top: t.y - GAP - h };
  if (t.x + t.w + GAP + W <= vw - 8) return { left: t.x + t.w + GAP, top: clampY(t.y) };
  if (t.x - GAP - W >= 8) return { left: t.x - GAP - W, top: clampY(t.y) };
  return { left: vw - W - 8, top: vh - h - 8 };
}

export function CghGuide() {
  const on = useCgh().help;
  const [i, setI] = useState(0);
  const [box, setBox] = useState<Box | null>(null);
  const [h, setH] = useState(160);
  const card = useRef<HTMLDivElement>(null);
  const step = STEPS[i];
  const end = () => { cgh.setHelp(false); setI(0); };
  const go = (k: number) => setI(Math.max(0, Math.min(STEPS.length - 1, k)));
  // the step may move on by itself meanwhile (its action done): never skip one
  const next = async () => { if (i === STEPS.length - 1) { end(); return; } const at = i; if (!step.done?.()) await step.act?.(); setI(k => (k === at ? Math.min(STEPS.length - 1, k + 1) : k)); };
  useEffect(() => {
    if (!on) return;
    // capture: before the page's own Esc; arrows left to fields
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') { e.stopPropagation(); end(); return; }
      if (/^(INPUT|TEXTAREA|SELECT)$/.test((e.target as HTMLElement)?.tagName)) return;
      if (e.key === 'ArrowRight') setI(k => Math.min(STEPS.length - 1, k + 1));
      else if (e.key === 'ArrowLeft') setI(k => Math.max(0, k - 1));
    };
    window.addEventListener('keydown', onKey, true);
    return () => window.removeEventListener('keydown', onKey, true);
  }, [on]); // eslint-disable-line react-hooks/exhaustive-deps
  // each step: its target scrolled into view (under the page's sticky bar), then followed frame by frame
  useEffect(() => {
    if (!on) return;
    let raf = 0, scrolled = false, was = !!step.done?.();
    const tick = () => {
      const el = step.target();
      if (el && !scrolled) {
        scrolled = true;
        const r = el.getBoundingClientRect();
        // always: a row may lie in the page's view but outside its table's own scroll box
        (el as HTMLElement).style.scrollMarginTop = '64px';
        el.scrollIntoView({ block: r.height > window.innerHeight * 0.5 ? 'start' : 'center' });
      }
      const r = el?.getBoundingClientRect();
      // a tall target is outlined on its visible part only, so that the popover finds room beside it
      const nb = r && r.width > 0 ? (() => { const y = Math.max(r.top, 52), b = Math.min(r.bottom, y + window.innerHeight * 0.45); return { x: r.left, y, w: r.width, h: Math.max(8, b - y) }; })() : null;
      setBox(o => (o && nb && Math.abs(o.x - nb.x) + Math.abs(o.y - nb.y) + Math.abs(o.w - nb.w) + Math.abs(o.h - nb.h) < 1) || (!o && !nb) ? o : nb);
      const d = !!step.done?.();
      if (d && !was) { setI(k => (k === i ? Math.min(STEPS.length - 1, k + 1) : k)); return; }
      was = d;
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [on, i]); // eslint-disable-line react-hooks/exhaustive-deps
  useLayoutEffect(() => { if (card.current) setH(card.current.offsetHeight); });
  if (!on) return null;
  const pos = place(box && { x: box.x - PAD, y: box.y - PAD, w: box.w + 2 * PAD, h: box.h + 2 * PAD }, h);
  const own = !!cgh.get().files.length && cgh.get().files.some(f => f.samples.some(x => !EXAMPLE_NAMES.includes(x.name)));
  const btn = 'h-[26px] px-2.5 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed';
  return (
    <div className="fixed inset-0 z-[90] pointer-events-none" data-cgh-guide>
      {box
        ? <div className="fixed rounded-lg border-2 border-accent-500 shadow-[0_0_0_9999px_rgba(15,23,42,0.28)] transition-all duration-150" data-cgh-spot
            style={{ left: box.x - PAD, top: box.y - PAD, width: box.w + 2 * PAD, height: box.h + 2 * PAD }} />
        : <div className="fixed inset-0 bg-slate-900/20" />}
      <div ref={card} role="dialog" aria-label={`Quick guide: ${step.title}`} data-cgh-tour
        className="fixed pointer-events-auto bg-white rounded-xl shadow-2xl border border-slate-200 text-slate-900 text-[12.5px] leading-5"
        style={{ left: pos.left, top: pos.top, width: W }}>
        <div className="px-3.5 pt-2.5 flex items-baseline gap-2">
          <div className="font-bold text-[13.5px]">{step.title}</div>
          <div className="ml-auto text-[11px] text-slate-500 whitespace-nowrap">step {i + 1} of {STEPS.length}</div>
        </div>
        <div className="px-3.5 py-1.5 text-slate-700">{step.body}
          {i === 0 && own && <div className="mt-1 text-[11.5px] text-amber-800">The tour adds the examples beside your patients; close them afterwards with × on their tabs.</div>}
          {!box && <div className="mt-1 text-[11px] text-slate-500">(its place is not on the page yet: Next goes on)</div>}
        </div>
        <div className="px-3.5 pb-2.5 pt-1 flex items-center gap-1.5">
          <button className="text-[11.5px] text-slate-500 hover:text-slate-800 hover:underline" onClick={end} title="End the tour (Esc)">Skip tour</button>
          <span className="ml-auto" />
          <button className={btn} disabled={i === 0} onClick={() => go(i - 1)}>Back</button>
          <button className="h-[26px] px-3 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700" onClick={() => void next()} data-cgh-next>{i === STEPS.length - 1 ? 'Done' : 'Next'}</button>
        </div>
      </div>
    </div>
  );
}
