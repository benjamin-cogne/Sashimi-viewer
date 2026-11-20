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
  /** done by the tour as the step opens: the user only clicks Next */
  enter?: () => void | Promise<void>;
}
const DEL = 'EXAMPLE-DEL22Q11', DER = 'EXAMPLE-DER', PIEV = 'EXAMPLE-PIEV-DUP';
const moreBtn = () => q<HTMLButtonElement>('button[aria-label="More filters"]');
const moreOpen = () => moreBtn()?.getAttribute('aria-expanded') === 'true';
const closeMore = () => { if (moreOpen()) moreBtn()?.click(); };
const STEPS: Step[] = [
  { title: 'Add samples, or the examples', target: () => q('button[aria-label="Examples"]'),
    body: <>Your files go in with <B>Add samples</B> (or dropped on the page): interval tables, BED, VCF, zips; read here, never sent. <B>Examples ▾</B> adds invented patients: Next loads the three of them for this tour.</> },
  { title: 'The import question', target: () => q('[data-cgh-question]'),
    enter: () => { if (!pendingExamples() && !loaded()) return loadExamples('all'); },
    body: <>Calls with fewer <B>probes</B> than the minimum (3), or under the <B>minimum size</B>, are not analysed. Each patient's <B>sex</B> is set here (read from the file when it says so). Next clicks <B>Analyse</B>.</> },
  { title: 'Patients', target: () => q('[role="tablist"][aria-label="Patients"]'),
    enter: () => { if (pendingExamples()) q<HTMLButtonElement>('[data-cgh-analyse]')?.click(); },
    body: <>One tab per patient, and <B>All patients</B> in one table with a Sample column. <B>×</B> on a tab takes the patient off (the examples too).</> },
  { title: 'Filters', target: () => q('[data-cgh-filters]'), enter: closeMore,
    body: <><B>probes ≥</B>, <B>losses / gains</B>, <B>size</B>, <B>class</B>; search a gene or a band (16p11); order by <B>priority</B> to see HIGH first.</> },
  { title: 'More filters', target: () => moreBtn()?.parentElement?.querySelector('[role="dialog"]') ?? moreBtn(), enter: () => { if (!moreOpen()) moreBtn()?.click(); },
    body: <><B>More filters</B>: <B>series</B> (how many other loaded patients carry the call), <B>lab list</B>, <B>hide suppressed</B>.</> },
  { title: 'Priority', target: () => rowOf(DEL, 'chr22', 'loss')?.querySelector('[data-cgh-priority]') ?? null, enter: closeMore,
    body: <>Each call gets <B>HIGH</B>, <B>MODERATE</B> or <B>LOW</B>: an order for the review, not a class. Here the <B>22q11.2</B> deletion of {DEL}.</> },
  // the detail itself while its genes are still read
  { title: 'Why this priority', target: () => q('[data-cgh-reasons]') ?? shownRow()?.nextElementSibling ?? null,
    enter: () => open(DEL, 'chr22', 'loss'),
    body: <>A click on a call opens its detail: the rules it meets, after the ACMG / ClinGen CNV standards (here a whole ClinGen HI 3 region, 2A). Every rule is listed under the table, in <B>ⓘ How the calls are read…</B></> },
  { title: 'Interpret in the browser', target: () => q('[data-cgh-browser]'), enter: () => open(DEL, 'chr22', 'loss'),
    body: <>The call is drawn below in the SV browser's own browser: genes, <B>DGV Gold</B>, <B>ClinVar</B>, <B>AChro-Puce</B>, <B>PanelApp</B>. Read the call there before classifying it.</> },
  { title: 'A possible der()', target: () => rowOf(DER, 'chr16', 'loss')?.querySelector('[data-cgh-topo]') ?? rowOf(DER, 'chr16', 'loss'),
    body: <>{DER} has a terminal 16q loss and a terminal 9p gain: <B>ter ⚠</B> hints a derivative of an unbalanced translocation, to check by parental karyotype / FISH.</> },
  { title: 'A gene cut', target: () => q('[data-cgh-gene="ANKRD11"]') ?? shownRow()?.nextElementSibling ?? null,
    enter: () => open(DER, 'chr16', 'loss'),
    body: <>The 16q deletion breaks inside <B>ANKRD11</B> (HI 3, KBG syndrome): the exons in the call and the effect on the frame are listed, read from the call's ends (theoretical: the true breakpoint lies between probes).</> },
  { title: 'A susceptibility CNV', target: () => q('[data-cgh-achro]') ?? rowOf(PIEV, 'chr1', 'gain'),
    enter: () => open(PIEV, 'chr1', 'gain'),
    body: <>{PIEV}: a 1q21.1 distal duplication that AChro-Puce classes <B>PIEV</B> (incomplete penetrance, variable expressivity): <B>MODERATE</B>, with the ClinGen TS 2 region. Parental segregation is advised.</> },
  { title: 'The lab benign list', target: () => q('[data-cgh-lab]'),
    body: <><B>Lab benign CNVs (BED)…</B> loads the laboratory's benign list, kept in this browser. Calls sharing 80 % with an entry are tagged, with a <B>→ Benign?</B> suggestion.</> },
  { title: 'Retry', target: () => q('[data-cgh-thead]'),
    body: <>When a source (UCSC, PanelApp, gnomAD, DGV) does not answer, <B>↻ retry (n)</B> appears in its column header, and <B>retry</B> in the call's detail: only what failed is asked again.</> },
  { title: 'Classify and note', target: () => shownRow()?.querySelector('[data-cgh-class]') ?? null, enter: () => open(PIEV, 'chr1', 'gain'),
    body: <>Choose the class in <B>Classification</B>; <B>✎</B> (or the box in the detail) adds a note. Both stay in this tab only.</> },
  { title: 'Export', target: () => q('[data-cgh-export-btn]'),
    body: <><B>Export classified table</B> writes each file back with Classification and Comment columns. The examples can then be closed with <B>×</B> on their tabs.</> },
];

/** "? Quick guide" in the page's bar: starts the tour */
export function CghGuideButton() {
  return <button onClick={() => cgh.setHelp(true)} className="ml-auto inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="An interactive tour of this page on invented example patients">? Quick guide</button>;
}

type Box = { x: number; y: number; w: number; h: number };
const PAD = 6, GAP = 12, W = 340, TOP = 56;
const hit = (a: Box, b: Box) => a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h;
const boxOf = (el: Element | null | undefined): Box | null => {
  const r = el?.getBoundingClientRect();
  if (!r || r.width <= 0) return null;
  // the part on screen only
  const y = Math.max(r.top, TOP - 4), b = Math.min(r.bottom, window.innerHeight);
  return b > y ? { x: r.left, y, w: r.width, h: b - y } : null;
};
/** where the popover goes: beside the target (below, above, right, left), else a corner; never over what is to be
 * clicked, and over the target only when nothing else fits */
function place(t: Box | null, c: Box | null, h: number): { left: number; top: number } {
  const vw = window.innerWidth, vh = window.innerHeight;
  const cx = (x: number) => Math.max(8, Math.min(vw - W - 8, x)), cy = (y: number) => Math.max(TOP, Math.min(vh - h - 8, y));
  const corners = [{ left: vw - W - 16, top: vh - h - 16 }, { left: 16, top: vh - h - 16 }, { left: vw - W - 16, top: TOP + 8 }, { left: 16, top: TOP + 8 }];
  if (!t) return corners[0];
  const side = [
    { left: cx(t.x), top: t.y + t.h + GAP }, { left: cx(t.x), top: t.y - GAP - h },
    { left: t.x + t.w + GAP, top: cy(t.y) }, { left: t.x - GAP - W, top: cy(t.y) },
  ];
  const fits = (p: { left: number; top: number }) => p.left >= 8 && p.left + W <= vw - 8 && p.top >= TOP && p.top + h <= vh - 8;
  const pad = (b: Box) => ({ x: b.x - PAD, y: b.y - PAD, w: b.w + 2 * PAD, h: b.h + 2 * PAD });
  const free = (p: { left: number; top: number }, avoid: (Box | null)[]) => avoid.every(a => !a || !hit({ x: p.left, y: p.top, w: W, h }, pad(a)));
  return [...side, ...corners].find(p => fits(p) && free(p, [t, c]))
    ?? [...side, ...corners].find(p => fits(p) && free(p, [c]))
    ?? corners[0];
}

export function CghGuide() {
  const on = useCgh().help;
  const [i, setI] = useState(0);
  const [box, setBox] = useState<Box | null>(null);
  const [h, setH] = useState(160);
  // the popover dragged by its title, for this step
  const [drag, setDrag] = useState({ dx: 0, dy: 0 });
  useEffect(() => setDrag({ dx: 0, dy: 0 }), [i]);
  const grab = (e: React.PointerEvent) => {
    const x0 = e.clientX - drag.dx, y0 = e.clientY - drag.dy;
    const mv = (m: PointerEvent) => setDrag({ dx: m.clientX - x0, dy: m.clientY - y0 });
    const up = () => { window.removeEventListener('pointermove', mv); window.removeEventListener('pointerup', up); };
    window.addEventListener('pointermove', mv); window.addEventListener('pointerup', up);
  };
  const card = useRef<HTMLDivElement>(null);
  const step = STEPS[i];
  const end = () => { cgh.setHelp(false); setI(0); };
  const go = (k: number) => setI(Math.max(0, Math.min(STEPS.length - 1, k)));
  // the step may move on by itself meanwhile (its action done): never skip one
  const next = () => { if (i === STEPS.length - 1) end(); else go(i + 1); };
  // the step's own action, once, as it opens
  useEffect(() => { if (on) void STEPS[i].enter?.(); }, [on, i]);
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
    let raf = 0, scrolled = false;
    const tick = () => {
      const el = step.target();
      if (el && !scrolled) {
        scrolled = true;
        const r = el.getBoundingClientRect();
        // always: a row may lie in the page's view but outside its table's own scroll box
        (el as HTMLElement).style.scrollMarginTop = '64px';
        el.scrollIntoView({ block: r.height > window.innerHeight * 0.5 ? 'start' : 'center' });
      }
      const same = (o: Box | null, n: Box | null) => (!o && !n) || (!!o && !!n && Math.abs(o.x - n.x) + Math.abs(o.y - n.y) + Math.abs(o.w - n.w) + Math.abs(o.h - n.h) < 1);
      const nb = boxOf(el);
      setBox(o => (same(o, nb) ? o : nb));
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [on, i]); // eslint-disable-line react-hooks/exhaustive-deps
  useLayoutEffect(() => { if (card.current) setH(card.current.offsetHeight); });
  if (!on) return null;
  const at = place(box && { x: box.x - PAD, y: box.y - PAD, w: box.w + 2 * PAD, h: box.h + 2 * PAD }, null, h);
  const pos = { left: at.left + drag.dx, top: at.top + drag.dy };
  const own = !!cgh.get().files.length && cgh.get().files.some(f => f.samples.some(x => !EXAMPLE_NAMES.includes(x.name)));
  const btn = 'h-[26px] px-2.5 rounded-[6px] border border-slate-200 text-slate-700 hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed';
  return (
    <div className="fixed inset-0 z-[90] pointer-events-none" data-cgh-guide>
      {box
        ? <div className="fixed rounded-lg border-2 border-accent-500 shadow-[0_0_0_9999px_rgba(15,23,42,0.28)] transition-all duration-150" data-cgh-spot
            style={{ left: box.x - PAD, top: box.y - PAD, width: box.w + 2 * PAD, height: box.h + 2 * PAD }} />
        : null}
      <div ref={card} role="dialog" aria-label={`Quick guide: ${step.title}`} data-cgh-tour
        className="fixed pointer-events-auto bg-white rounded-xl shadow-2xl border border-slate-200 text-slate-900 text-[12.5px] leading-5"
        style={{ left: pos.left, top: pos.top, width: W }}>
        <div className="px-3.5 pt-2.5 flex items-baseline gap-2 cursor-move select-none" onPointerDown={grab} title="Drag to move this box">
          <div className="font-bold text-[13.5px]">{step.title}</div>
          <div className="ml-auto text-[11px] text-slate-500 whitespace-nowrap">step {i + 1} of {STEPS.length}</div>
        </div>
        <div className="px-3.5 py-1.5 text-slate-700">{step.body}
          {i === 0 && own && <div className="mt-1 text-[11.5px] text-amber-800">The tour adds the examples beside your patients; close them afterwards with × on their tabs.</div>}
        </div>
        <div className="px-3.5 pb-2.5 pt-1 flex items-center gap-1.5">
          <button className="text-[11.5px] text-slate-500 hover:text-slate-800 hover:underline" onClick={end} title="End the tour (Esc)">Skip tour</button>
          <span className="ml-auto" />
          <button className={btn} disabled={i === 0} onClick={() => go(i - 1)}>Back</button>
          <button className="h-[26px] px-3 rounded-[6px] bg-accent-600 text-white font-semibold hover:bg-accent-700" onClick={next} data-cgh-next>{i === STEPS.length - 1 ? 'Done' : 'Next'}</button>
        </div>
      </div>
    </div>
  );
}
