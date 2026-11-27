/**
 * The array-CGH page's quick guide: an interactive tour on the built-in examples (cghExamples.ts). Each step outlines
 * its target (the rest of the page dimmed) and a popover beside it, never over it, with a pointer toward it; each step
 * does its own action as it opens, so the user only clicks Next. The page stays usable under it. Esc or Skip ends it;
 * the store's `help` keeps Esc from closing the page too.
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
  /** outline only the target's top part (a tall target: the popover fits under it) */
  maxH?: number;
  /** done by the tour as the step opens: the user only clicks Next */
  enter?: () => void | Promise<void>;
}
const DEL = 'EXAMPLE-DEL22Q11', DER = 'EXAMPLE-DER', PIEV = 'EXAMPLE-PIEV-DUP', PO = 'EXAMPLE-PROBES-ONLY', MOS = 'EXAMPLE-MOSAIC-T12';
const keyOf = (name: string) => cgh.get().files.flatMap(f => f.samples).find(x => x.name === name)?.key ?? null;
/** a patient's own tab */
const tab = (name: string) => { const k = keyOf(name); if (k && cgh.get().sample !== k) cgh.setSample(k); };
const frame = () => new Promise(r => requestAnimationFrame(() => r(null)));
/** wait (up to `ms`) for an element to appear */
const waitFor = async (f: () => Element | null, ms = 3000) => { const t = performance.now(); while (!f() && performance.now() - t < ms) await frame(); return f(); };
/** the whole genome panel open (it folds to a strip while a call is open) */
const genomeOpen = async () => { if (!q('[data-cgh-genome]')) q<HTMLButtonElement>('[data-cgh-genome-strip]')?.click(); await waitFor(() => q('[data-cgh-genome]')); };
const moreBtn = () => q<HTMLButtonElement>('button[aria-label="More filters"]');
const moreOpen = () => moreBtn()?.getAttribute('aria-expanded') === 'true';
const closeMore = () => { if (moreOpen()) moreBtn()?.click(); };
const STEPS: Step[] = [
  { title: 'Input files', target: () => q('[data-cgh-input]')?.parentElement ?? null,
    body: <>Add <B>CytoGenomics interval tables</B> (.txt), its <B>automated workflow zip</B> (IntervalBasedReport with <B>CGH_Log_Ratio.txt</B>: the probes), other probe tables (Illumina, ChAS-like, CNVkit) or any table of calls, BED, VCF. An odd file opens a <B>column dialog</B>. Everything is read here, never sent.</> },
  { title: 'The examples', target: () => q('button[aria-label="Examples"]'),
    body: <><B>Examples ▾</B> adds invented patients with their probes: three tables and one probe file alone. Next loads the four of them for this tour.</> },
  { title: 'The import question', target: () => q('[data-cgh-question]'), maxH: 440,
    enter: () => { if (!pendingExamples() && !loaded()) return loadExamples('all'); },
    body: <>Each patient's probe file is listed with its <B>DLRS</B> (the array's noise: under 0.2 excellent). Calls under the minimum <B>probes</B> or <B>size</B> are not analysed; the <B>sex</B> is set per patient. Next clicks <B>Analyse</B>.</> },
  { title: 'Patients', target: () => q('[role="tablist"][aria-label="Patients"]'),
    enter: async () => { if (pendingExamples()) { q<HTMLButtonElement>('[data-cgh-analyse]')?.click(); await frame(); } cgh.setSample(null); },
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
  { title: 'Interpret in the browser', target: () => q('[data-cgh-svb]'), maxH: 340, enter: () => open(DEL, 'chr22', 'loss'),
    body: <>The call is drawn below in the SV browser's own browser: genes, <B>DGV Gold</B>, <B>ClinVar</B>, <B>AChro-Puce</B>, <B>PanelApp</B>, ClinGen regions. Read the call there before classifying it.</> },
  { title: "The call's ISCN", target: () => q('[data-cgh-browser] [data-cgh-iscn]'), enter: () => open(DEL, 'chr22', 'loss'),
    body: <>The call's card gives its <B>ISCN</B> line (bands, positions, copies), with <B>copy</B> to paste it in a report; the file's own ISCN is shown beside when it differs.</> },
  { title: "The call's probes", target: () => q('[data-cgh-probes]'), enter: () => open(DEL, 'chr22', 'loss'),
    body: <><B>Log ratio (probes)</B>: the call centred with a margin, each probe a dot, the moving median, the ideogram under it. Hover a probe for its tip; the header gives the <B>probes' mean</B> beside the <B>file's mean</B>.</> },
  { title: 'A possible der()', target: () => rowOf(DER, 'chr16', 'loss')?.querySelector('[data-cgh-topo]') ?? rowOf(DER, 'chr16', 'loss'),
    body: <>{DER} has a terminal 16q loss and a terminal 9p gain: <B>ter ⚠</B> hints a derivative of an unbalanced translocation, to check by parental karyotype / FISH.</> },
  { title: 'A gene cut', target: () => q('[data-cgh-gene="ANKRD11"]') ?? shownRow()?.nextElementSibling ?? null,
    enter: () => open(DER, 'chr16', 'loss'),
    body: <>The 16q deletion breaks inside <B>ANKRD11</B> (HI 3, KBG syndrome): the exons in the call and the effect on the frame are listed, read from the call's ends (theoretical: the true breakpoint lies between probes).</> },
  { title: 'A susceptibility CNV', target: () => q('[data-cgh-achro]') ?? rowOf(PIEV, 'chr1', 'gain'),
    enter: () => open(PIEV, 'chr1', 'gain'),
    body: <>{PIEV}: a 1q21.1 distal duplication that AChro-Puce classes <B>PIEV</B> (incomplete penetrance, variable expressivity): <B>MODERATE</B>, with the ClinGen TS 2 region. Parental segregation is advised.</> },
  { title: 'Classify and note', target: () => shownRow()?.querySelector('[data-cgh-class]') ?? null, enter: () => { tab(PIEV); open(PIEV, 'chr1', 'gain'); },
    body: <>Choose the class in <B>Classification</B>; <B>✎</B> (or the box in the detail) adds a note. Both stay in this tab only.</> },
  { title: "The patient's ISCN", target: () => q('[data-cgh-patient-iscn]'), enter: () => tab(PIEV),
    body: <>A patient's tab shows its <B>ISCN</B> line, built from the calls of the classes ticked (P / LP / VUS by default), with <B>copy</B>. With none, the normal result for the patient's sex: here a male.</> },
  { title: 'A CNV group', target: () => q('[data-cgh-group]'), enter: () => tab(MOS),
    body: <>Calls of one direction covering a whole arm or chromosome (here 12p + 12q, with a call nested in 12q) make a <B>CNV group</B>, outlined: its weighted log ratio, probes, copies, ISCN and a <B>suggestion</B> (mosaic trisomy 12), its own class (<B>→ calls</B> gives it to its calls). <B>ungroup</B> splits it back.</> },
  { title: 'The whole genome', target: () => q('[data-cgh-genome]'), enter: async () => { tab(DEL); await genomeOpen(); },
    body: <>Every probe of the patient over chromosomes 1-22, X, Y, the file's calls as red / blue bars above, the probe calls in <B>purple</B>. <B>horizontal</B> / <B>vertical</B> sets the layout.</> },
  { title: 'A chromosome', target: () => q('[data-cgh-chrom-probes]'),
    enter: async () => { tab(DEL); await genomeOpen(); window.dispatchEvent(new CustomEvent('cgh-chrom', { detail: 'chr22' })); },
    body: <>A click on a chromosome opens it: every probe, the moving median, the <B>ideogram</B> (click or drag it to move); <B>⊖ ⊕</B>, the slider or a drag zoom and move. A call bar clicked opens the call.</> },
  { title: 'Log-ratio range', target: () => q('[data-cgh-genome] [data-cgh-scale]'),
    enter: async () => { tab(DEL); await genomeOpen(); if (!q('[data-cgh-chrom-probes]')) window.dispatchEvent(new CustomEvent('cgh-chrom', { detail: 'chr22' })); },
    body: <>The small buttons left of the axis set the log-ratio range (<B>±0.5 / ±1 / ±2</B>, ⊕ ⊖), shared by the three views; <B>show probes</B> above hides the dots, keeping the median and the calls.</> },
  { title: 'Probe calls', target: () => q('[data-cgh-probe-calls]'),
    enter: () => { tab(DEL); const k = keyOf(DEL); if (k) cgh.setProbeOpen(k, true); },
    body: <>Events found in the probes by the page (segmentation), a second opinion: <B>n not in the file</B> lists those the file does not call, here a small 3p14.2 loss. Never exported; <B>settings</B> sets the sensitivity.</> },
  { title: 'A probes-only patient', target: () => q('[data-cgh-probeonly-banner]'), enter: () => tab(PO),
    body: <>{PO} came as a probe file alone: its calls are the probe calls, marked <B>probes only</B>, to confirm in the analysis software. Its export is a new <B>.probe-calls.txt</B>.</> },
  { title: 'The lab benign list', target: () => q('[data-cgh-lab]'),
    body: <><B>Lab benign CNVs (BED)…</B> loads the laboratory's benign list, kept in this browser. Calls sharing 80 % with an entry are tagged, with a <B>→ Benign?</B> suggestion.</> },
  { title: 'Retry', target: () => q('[data-cgh-thead]'),
    body: <>When a source (UCSC, PanelApp, gnomAD, DGV) does not answer, <B>↻ retry (n)</B> appears in its column header, and <B>retry</B> in the call's detail: only what failed is asked again.</> },
  { title: 'Export', target: () => q('[data-cgh-export-btn]'),
    body: <><B>Export classified table</B> writes each file back with Classification and Comment columns. The examples can then be closed with <B>×</B> on their tabs.</> },
  { title: 'Methods', target: () => q('[data-cgh-disclaimer]'), enter: () => { const d = q<HTMLDetailsElement>('[data-cgh-methods]'); if (d) d.open = true; },
    body: <><B>ⓘ</B> under the table gives the methods: <B>DLRS</B>, probe-call <B>detection</B>, copies, ISCN. A complementary tool, tested on Agilent 60K arrays: check everything in the reference software.</> },
  { title: 'Esc', target: () => q('[data-cgh-esc-hint]'),
    body: <>On the page, <B>Esc</B> first clears the view (call, detail, whole genome); a second Esc asks before going back to the SV browser. During this tour, Esc ends the tour.</> },
];

/** "? Quick guide" in the page's bar: starts the tour */
export function CghGuideButton() {
  return <button onClick={() => cgh.setHelp(true)} className="ml-auto inline-flex items-center gap-1.5 h-7 px-2.5 rounded-lg text-xs font-medium border bg-bar-ctrl border-bar-line text-bar-text hover:border-accent-300" title="An interactive tour of this page on invented example patients">? Quick guide</button>;
}

type Box = { x: number; y: number; w: number; h: number };
type Side = 'below' | 'above' | 'right' | 'left' | 'spot';
type At = { left: number; top: number; side: Side };
const PAD = 6, GAP = 14, W = 340, TOP = 56;
/** how long a target may vanish (a redraw, a step's action under way) before the popover gives up its place */
const HOLD_MS = 1500;
const overlap = (a: Box, b: Box) => Math.max(0, Math.min(a.x + a.w, b.x + b.w) - Math.max(a.x, b.x)) * Math.max(0, Math.min(a.y + a.h, b.y + b.h) - Math.max(a.y, b.y));
/** the part of an element on screen: clipped by its scrolling ancestors (a row hidden in the table's own scroll box is not
 * on screen), the page's sticky bar and the window; `maxH`: only its top part */
const boxOf = (el: Element | null | undefined, maxH = Infinity): Box | null => {
  const r = el?.getBoundingClientRect();
  if (!el || !r || r.width <= 0) return null;
  let x0 = r.left, x1 = r.right, y0 = r.top, y1 = Math.min(r.bottom, r.top + maxH);
  for (let a = el.parentElement; a && a !== document.body; a = a.parentElement) {
    const st = getComputedStyle(a);
    if (!/(auto|scroll|hidden|clip)/.test(st.overflowX + st.overflowY)) continue;
    const b = a.getBoundingClientRect();
    x0 = Math.max(x0, b.left); x1 = Math.min(x1, b.right); y0 = Math.max(y0, b.top); y1 = Math.min(y1, b.bottom);
  }
  // the page's sticky bar hides what scrolls under it, but its own buttons are targets too
  x0 = Math.max(x0, 0); x1 = Math.min(x1, window.innerWidth); y0 = r.bottom <= TOP ? Math.max(y0, 0) : Math.max(y0, TOP - 4); y1 = Math.min(y1, window.innerHeight);
  return x1 - x0 > 2 && y1 - y0 > 2 ? { x: x0, y: y0, w: x1 - x0, h: y1 - y0 } : null;
};
/** where the popover goes: beside the target (below, above, right, left; the side it had kept while it still fits), never
 * over it; when nothing beside fits, the on-screen spot covering least of it; no target: bottom centre, a fixed place */
function place(t: Box | null, h: number, prev: Side | null): At {
  const vw = window.innerWidth, vh = window.innerHeight;
  const bottom: At = { left: Math.round((vw - W) / 2), top: vh - h - 16, side: 'spot' };
  if (!t) return bottom;
  const cx = (x: number) => Math.max(8, Math.min(vw - W - 8, x)), cy = (y: number) => Math.max(TOP, Math.min(vh - h - 8, y));
  const side: At[] = [
    { left: cx(t.x + t.w / 2 - W / 2), top: t.y + t.h + GAP, side: 'below' }, { left: cx(t.x + t.w / 2 - W / 2), top: t.y - GAP - h, side: 'above' },
    { left: t.x + t.w + GAP, top: cy(t.y + t.h / 2 - h / 2), side: 'right' }, { left: t.x - GAP - W, top: cy(t.y + t.h / 2 - h / 2), side: 'left' },
  ];
  const fits = (p: At) => p.left >= 8 && p.left + W <= vw - 8 && p.top >= TOP && p.top + h <= vh - 8;
  const cover = (p: At) => overlap({ x: p.left, y: p.top, w: W, h }, t);
  const ok = (p: At) => fits(p) && cover(p) === 0;
  const kept = prev ? side.find(p => p.side === prev && ok(p)) : undefined;
  if (kept) return kept;
  const first = side.find(ok);
  if (first) return first;
  const spots: At[] = [...side.filter(fits), bottom, { left: vw - W - 16, top: vh - h - 16, side: 'spot' }, { left: 16, top: vh - h - 16, side: 'spot' },
    { left: vw - W - 16, top: TOP + 8, side: 'spot' }, { left: 16, top: TOP + 8, side: 'spot' }];
  return spots.reduce((a, b) => (cover(b) < cover(a) ? b : a));
}

export function CghGuide() {
  const on = useCgh().help;
  const [i, setI] = useState(0);
  // the target's box; kept a moment when the target vanishes (`found` false: no spotlight meanwhile)
  const [box, setBox] = useState<Box | null>(null);
  const [found, setFound] = useState(false);
  const [h, setH] = useState(160);
  // the side chosen for this step, kept while it fits (no flip as the page scrolls)
  const sideRef = useRef<Side | null>(null);
  // the popover dragged by its title, for this step
  const [drag, setDrag] = useState({ dx: 0, dy: 0, on: false });
  useEffect(() => { setDrag({ dx: 0, dy: 0, on: false }); sideRef.current = null; }, [i]);
  const grab = (e: React.PointerEvent) => {
    const x0 = e.clientX - drag.dx, y0 = e.clientY - drag.dy;
    const mv = (m: PointerEvent) => setDrag({ dx: m.clientX - x0, dy: m.clientY - y0, on: true });
    const up = () => { setDrag(d => ({ ...d, on: false })); window.removeEventListener('pointermove', mv); window.removeEventListener('pointerup', up); };
    window.addEventListener('pointermove', mv); window.addEventListener('pointerup', up);
  };
  const card = useRef<HTMLDivElement>(null);
  const step = STEPS[i];
  const end = () => { cgh.setHelp(false); setI(0); setBox(null); setFound(false); };
  const go = (k: number) => setI(Math.max(0, Math.min(STEPS.length - 1, k)));
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
    let raf = 0, scrolled = 0, first = true, seen = performance.now();
    const tick = () => {
      const el = step.target();
      let nb = boxOf(el, step.maxH);
      // brought into view as it appears (always: a row may lie in the page's view but outside its table's own scroll box), and again while it stays hidden (a table's own scroll box, a panel opening above it)
      if (el && (first || (!nb && performance.now() - scrolled > 400))) {
        first = false; scrolled = performance.now();
        const r = el.getBoundingClientRect();
        (el as HTMLElement).style.scrollMarginTop = '64px';
        el.scrollIntoView({ block: r.height > (step.maxH ?? Infinity) || r.height > window.innerHeight * 0.5 ? 'start' : 'center' });
        nb = boxOf(el, step.maxH);
      }
      const same = (o: Box | null, n: Box | null) => (!o && !n) || (!!o && !!n && Math.abs(o.x - n.x) + Math.abs(o.y - n.y) + Math.abs(o.w - n.w) + Math.abs(o.h - n.h) < 1);
      if (nb) { seen = performance.now(); setBox(o => (same(o, nb) ? o : nb)); setFound(true); }
      else { setFound(false); if (performance.now() - seen > HOLD_MS) setBox(null); }
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [on, i]); // eslint-disable-line react-hooks/exhaustive-deps
  useLayoutEffect(() => { if (card.current && Math.abs(card.current.offsetHeight - h) > 1) setH(card.current.offsetHeight); });
  if (!on) return null;
  const t = box && { x: box.x - PAD, y: box.y - PAD, w: box.w + 2 * PAD, h: box.h + 2 * PAD };
  const at = place(t, h, sideRef.current);
  if (t) sideRef.current = at.side;
  const pos = { left: at.left + drag.dx, top: at.top + drag.dy };
  // the pointer toward the target, on the popover's side facing it (none when moved by hand or in a fixed spot)
  const arrow = t && found && at.side !== 'spot' && !drag.dx && !drag.dy
    ? (at.side === 'below' || at.side === 'above'
      ? { side: at.side, off: Math.max(16, Math.min(W - 16, t.x + t.w / 2 - at.left)) }
      : { side: at.side, off: Math.max(16, Math.min(h - 16, t.y + t.h / 2 - at.top)) })
    : null;
  const own = !!cgh.get().files.length && cgh.get().files.some(f => f.samples.some(x => !EXAMPLE_NAMES.includes(x.name)));
  const btn = 'h-[26px] px-2.5 rounded-[6px] border border-slate-300 text-slate-700 bg-white hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed';
  const tri = 'absolute w-0 h-0 border-solid border-transparent';
  return (
    <div className="fixed inset-0 z-[90] pointer-events-none" data-cgh-guide>
      {box && found
        ? <div className="fixed rounded-lg border-2 border-accent-500 shadow-[0_0_0_9999px_rgba(15,23,42,0.42)] transition-all duration-200" data-cgh-spot
            style={{ left: box.x - PAD, top: box.y - PAD, width: box.w + 2 * PAD, height: box.h + 2 * PAD }} />
        : null}
      <div ref={card} role="dialog" aria-label={`Quick guide: ${step.title}`} data-cgh-tour data-side={at.side}
        className="fixed pointer-events-auto bg-white rounded-xl border-2 border-accent-600 text-slate-900 text-[12.5px] leading-5 shadow-[0_18px_48px_-6px_rgba(15,23,42,0.6),0_0_0_4px_rgba(255,255,255,0.85)]"
        style={{ left: pos.left, top: pos.top, width: W, transition: drag.on ? 'none' : 'left 220ms ease-out, top 220ms ease-out' }}>
        {arrow?.side === 'below' && <span className={`${tri} border-x-[9px] border-b-[10px] border-b-accent-600`} style={{ left: arrow.off - 9, top: -12 }} />}
        {arrow?.side === 'above' && <span className={`${tri} border-x-[9px] border-t-[10px] border-t-accent-600`} style={{ left: arrow.off - 9, bottom: -12 }} />}
        {arrow?.side === 'right' && <span className={`${tri} border-y-[9px] border-r-[10px] border-r-accent-600`} style={{ top: arrow.off - 9, left: -12 }} />}
        {arrow?.side === 'left' && <span className={`${tri} border-y-[9px] border-l-[10px] border-l-accent-600`} style={{ top: arrow.off - 9, right: -12 }} />}
        <div className="px-3.5 py-1.5 flex items-baseline gap-2 cursor-move select-none bg-accent-600 text-white rounded-t-[9px]" onPointerDown={grab} title="Drag to move this box">
          <div className="font-bold text-[13.5px]">{step.title}</div>
          <div className="ml-auto text-[11px] text-white/85 whitespace-nowrap">step {i + 1} of {STEPS.length}</div>
        </div>
        <div className="px-3.5 pt-2 pb-1.5 text-slate-700">{step.body}
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
