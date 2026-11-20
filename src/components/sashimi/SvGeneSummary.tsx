/**
 * What a CNV holds, at a glance, next to its nomenclature: its coding genes (whole or cut by a breakpoint), those
 * green in PanelApp with their mode of inheritance, those intolerant of loss of function (gnomAD LOEUF under the
 * recommended threshold). Each count opens its list of genes, which can be copied; a gene clicked is shown in the
 * browser.
 *
 * The SV browser fills it (SvGenomeBrowser's onGeneSummary): the genes the event changes (story), their constraint
 * (gnomAD) and their PanelApp entries (every panel of the source, asked for each changed gene).
 */
import { useState } from 'react';
import type { MoiClass } from '../../standalone/panelapp';
import { LOEUF_INTOLERANT, LOEUF_CONSTRAINED, type ImpactKind } from './svGenes';

/** Modes of inheritance (PanelApp), as colours and words. */
export const MOI: Record<MoiClass | 'none', { color: string; short: string; label: string }> = {
  AD: { color: '#dc2626', short: 'AD', label: 'monoallelic (AD)' },
  'AD+AR': { color: '#7c3aed', short: 'AD/AR', label: 'mono- or biallelic' },
  AR: { color: '#2563eb', short: 'AR', label: 'biallelic (AR)' },
  XLR: { color: '#db2777', short: 'XLR', label: 'X-linked, females biallelic' },
  XLD: { color: '#9d174d', short: 'XLD', label: 'X-linked, females may be affected' },
  MT: { color: '#059669', short: 'MT', label: 'mitochondrial' },
  '?': { color: '#94a3b8', short: '?', label: 'unknown / other' },
  none: { color: '#cbd5e1', short: '—', label: 'not in PanelApp' },
};
const MOI_ORDER: MoiClass[] = ['AD', 'AD+AR', 'AR', 'XLD', 'XLR', 'MT', '?'];

export interface GeneSummaryItem {
  symbol: string; start: number; end: number;
  kind: ImpactKind;
  /** the whole gene lost, gained (or moved), or only part of it (a breakpoint inside) */
  whole: boolean;
  coding: boolean;
  loeuf: number | null; pli: number | null; version: string | null;
  /** its mode of inheritance where it is green in PanelApp; null: green nowhere; 'pending': not known yet */
  green: MoiClass | null | 'pending';
  /** the green panels' names (a few) */
  greenPanels: string[];
}
export interface GeneSummary {
  /** a gain or loss of copies (CNV), or a rearrangement without one (inversion…) */
  cnv: boolean;
  items: GeneSummaryItem[];
  /** constraint still loading */
  constraintLoading: boolean;
  /** the PanelApp asked (England / Australia) */
  source: string;
}

/** gnomAD's recommended LOEUF threshold for the version a score comes from (v4.1.1: 0.45; v2.1.1: 0.35) */
const intolerantBelow = (v: string | null) => LOEUF_INTOLERANT[v ?? 'v4'] ?? LOEUF_INTOLERANT.v4;
const generalBelow = (v: string | null) => LOEUF_CONSTRAINED[v ?? 'v4'] ?? LOEUF_CONSTRAINED.v4;

type Cat = 'coding' | 'green' | 'constrained';

export default function SvGeneSummary({ s, onGene }: { s: GeneSummary; onGene?: (g: GeneSummaryItem) => void }) {
  const [open, setOpen] = useState<Cat | null>(null);
  const [copied, setCopied] = useState<string | null>(null);
  const copy = (key: string, text: string) => {
    navigator.clipboard?.writeText(text).then(() => { setCopied(key); setTimeout(() => setCopied(c => (c === key ? null : c)), 1500); }).catch(() => { /* clipboard blocked */ });
  };
  const bySym = (a: GeneSummaryItem, b: GeneSummaryItem) => a.symbol.localeCompare(b.symbol);
  const coding = s.items.filter(i => i.coding).sort(bySym);
  const nonCoding = s.items.length - coding.length;
  const whole = coding.filter(i => i.whole), part = coding.filter(i => !i.whole);
  const green = s.items.filter(i => i.green && i.green !== 'pending').sort(bySym);
  const pending = s.items.filter(i => i.green === 'pending').length;
  const greenBy = MOI_ORDER.map(c => ({ c, genes: green.filter(i => i.green === c) })).filter(g => g.genes.length);
  const known = coding.filter(i => i.loeuf != null);
  const intolerant = known.filter(i => i.loeuf! < intolerantBelow(i.version)).sort((a, b) => a.loeuf! - b.loeuf!);
  const between = known.filter(i => i.loeuf! >= intolerantBelow(i.version) && i.loeuf! < generalBelow(i.version)).sort((a, b) => a.loeuf! - b.loeuf!);
  const tolerant = known.length - intolerant.length - between.length, unknownC = coding.length - known.length;
  const where = s.cnv ? 'in the CNV' : 'in the rearrangement';

  /** a bar of parts, each a share of the whole */
  const bar = (parts: { n: number; color: string; title: string }[]) => {
    const total = parts.reduce((a, p) => a + p.n, 0);
    return (
      <div className="mt-1.5 flex h-[7px] rounded-full overflow-hidden bg-slate-100">
        {total > 0 && parts.filter(p => p.n > 0).map((p, i) => <div key={i} style={{ width: `${(p.n / total) * 100}%`, background: p.color }} title={p.title} />)}
      </div>
    );
  };
  const tile = (cat: Cat, n: number | string, label: string, sub: React.ReactNode, parts: { n: number; color: string; title: string }[], accent: string, title: string) => (
    <button type="button" onClick={() => setOpen(o => (o === cat ? null : cat))} title={title}
      className={`text-left rounded-lg border px-2.5 py-2 transition-colors ${open === cat ? 'border-accent-400 bg-accent-50/60 ring-1 ring-accent-200' : 'border-slate-200 bg-white hover:border-slate-300 hover:bg-slate-50'}`}>
      <div className="flex items-baseline gap-1.5">
        <span className="text-[22px] leading-none font-bold tabular-nums" style={{ color: accent }}>{n}</span>
        <span className="text-[11px] font-semibold text-slate-700 leading-tight">{label}</span>
      </div>
      {bar(parts)}
      <div className="mt-1 text-[10.5px] text-slate-500 leading-snug">{sub}</div>
    </button>
  );
  const chip = (g: GeneSummaryItem, extra?: React.ReactNode, color?: string) => (
    <button key={g.symbol} type="button" onClick={() => onGene?.(g)}
      title={`${g.symbol}: ${g.whole ? 'whole gene' : 'part of the gene'} (${g.kind})${g.loeuf != null ? ` · LOEUF ${g.loeuf.toFixed(2)} (gnomAD ${g.version})` : ''}${g.greenPanels.length ? `\ngreen in: ${g.greenPanels.join('; ')}` : ''}\nclick: show it in the browser`}
      className="inline-flex items-center gap-1 h-[22px] px-1.5 rounded-md border border-slate-200 bg-white hover:border-accent-300 hover:bg-accent-50 text-[11px] font-semibold text-slate-800">
      {color && <span className="inline-block w-2 h-2 rounded-full" style={{ background: color }} />}
      <span className={g.whole ? '' : 'italic'}>{g.symbol}</span>{extra}
    </button>
  );
  const copyBtn = (key: string, genes: GeneSummaryItem[], what = 'copy') => (
    <button type="button" onClick={() => copy(key, genes.map(g => g.symbol).join(', '))} disabled={!genes.length}
      className="text-[10.5px] px-1.5 py-0.5 rounded border border-slate-300 text-slate-600 hover:bg-slate-50 disabled:opacity-40" title="Copy the gene symbols, comma-separated">
      {copied === key ? '✓ copied' : what}
    </button>
  );
  const group = (key: string, head: React.ReactNode, genes: GeneSummaryItem[], render: (g: GeneSummaryItem) => React.ReactNode) => genes.length > 0 && (
    <div key={key} className="mt-1.5">
      <div className="flex items-center gap-2 text-[10.5px] font-semibold text-slate-600">{head}<span className="font-normal text-slate-400">{genes.length}</span>{copyBtn(key, genes)}</div>
      <div className="mt-1 flex flex-wrap gap-1">{genes.map(render)}</div>
    </div>
  );
  // every gene with its readings, as a table for a report or a spreadsheet
  const table = () => ['gene\tcoding\tin the event\timpact\tLOEUF\tgnomAD\tPanelApp green (mode)',
    ...[...s.items].sort(bySym).map(g => [g.symbol, g.coding ? 'yes' : 'no', g.whole ? 'whole' : 'part', g.kind, g.loeuf != null ? g.loeuf.toFixed(3) : '', g.version ?? '',
      g.green === 'pending' ? '(not known yet)' : g.green ? MOI[g.green].short : ''].join('\t'))].join('\n');

  return (
    <div className="min-w-0">
      <div className="flex items-center gap-2">
        <span className="text-[10.5px] font-semibold text-slate-700">Genes {where}</span>
        <button type="button" onClick={() => copy('table', table())} className="text-[10.5px] text-accent-700 hover:underline" title="Copy every gene with its readings (tab-separated: a spreadsheet or a report table)">{copied === 'table' ? 'copied' : 'copy table'}</button>
      </div>
      <div className="mt-1 grid grid-cols-3 gap-1.5">
        {tile('coding', coding.length, 'coding genes',
          <>{whole.length} whole · {part.length} cut{nonCoding ? <> · +{nonCoding} non-coding</> : null}</>,
          [{ n: whole.length, color: '#0f766e', title: `${whole.length} whole` }, { n: part.length, color: '#5eead4', title: `${part.length} cut by a breakpoint` }],
          '#0f766e', `Protein-coding genes ${where}: wholly ${s.cnv ? 'lost or gained' : 'inside it'}, or cut by a breakpoint. Click for the list.`)}
        {tile('green', green.length, 'PanelApp green',
          greenBy.length ? <span className="inline-flex flex-wrap gap-x-1.5">{greenBy.map(g => <span key={g.c} style={{ color: MOI[g.c].color }} className="font-semibold">{MOI[g.c].short} {g.genes.length}</span>)}</span>
            : pending ? <>checking {pending} gene{pending > 1 ? 's' : ''}…</> : <>none</>,
          greenBy.map(g => ({ n: g.genes.length, color: MOI[g.c].color, title: `${MOI[g.c].label}: ${g.genes.length}` })),
          '#15803d', `Genes green (diagnostic grade) in at least one ${s.source} panel, by mode of inheritance (that of their green entries)${pending ? `; ${pending} still being checked` : ''}. Click for the list.`)}
        {tile('constrained', intolerant.length, 'constrained',
          <>LOEUF &lt; {intolerantBelow(intolerant[0]?.version ?? known[0]?.version ?? null)}{between.length ? <> · +{between.length} under {generalBelow(null)}</> : null}{s.constraintLoading ? ' · loading…' : unknownC ? <> · {unknownC} unknown</> : null}</>,
          [{ n: intolerant.length, color: '#b91c1c', title: `${intolerant.length} LOEUF under the recommended threshold` }, { n: between.length, color: '#fdba74', title: `${between.length} between it and ${generalBelow(null)}` }, { n: tolerant, color: '#cbd5e1', title: `${tolerant} not constrained` }, { n: unknownC, color: '#eef2f6', title: `${unknownC} without a score` }],
          '#b91c1c', `Coding genes ${where} intolerant of loss of function: gnomAD LOEUF under the recommended threshold (${LOEUF_INTOLERANT.v4} for v4.1.1, ${LOEUF_INTOLERANT['v2.1.1']} for v2.1.1); then those under the v4.0 cut-off (${LOEUF_CONSTRAINED.v4}). Click for the list.`)}
      </div>
      {open && (
        <div className="mt-1.5 rounded-lg border border-slate-200 bg-slate-50/70 px-2.5 py-1.5 max-h-[220px] overflow-auto">
          {open === 'coding' && <>
            {group('whole', s.cnv ? 'whole genes' : 'whole genes inside', whole, g => chip(g))}
            {group('part', 'cut by a breakpoint (part of the gene)', part, g => chip(g))}
            {!coding.length && <div className="text-[11px] text-slate-500">No coding gene {where}.</div>}
            {coding.length > 0 && <div className="mt-1.5">{copyBtn('allcoding', coding, `copy all ${coding.length}`)}</div>}
          </>}
          {open === 'green' && <>
            {greenBy.map(({ c, genes }) => group(`g${c}`, <span style={{ color: MOI[c].color }}>{MOI[c].label}</span>, genes, g => chip(g, null, MOI[c].color)))}
            {!green.length && <div className="text-[11px] text-slate-500">{pending ? `Checking ${pending} gene${pending > 1 ? 's' : ''} in ${s.source}…` : `No gene ${where} is green in ${s.source}.`}</div>}
            {green.length > 0 && <div className="mt-1.5 flex items-center gap-2">{copyBtn('allgreen', green, `copy all ${green.length}`)}{pending > 0 && <span className="text-[10.5px] text-slate-400">{pending} still being checked</span>}</div>}
          </>}
          {open === 'constrained' && <>
            {group('intol', <>LOEUF under the recommended threshold</>, intolerant, g => chip(g, <span className="font-mono font-normal text-[10px] text-red-700">{g.loeuf!.toFixed(2)}</span>, '#b91c1c'))}
            {group('between', <>LOEUF {LOEUF_INTOLERANT.v4}–{LOEUF_CONSTRAINED.v4} (the v4.0 cut-off)</>, between, g => chip(g, <span className="font-mono font-normal text-[10px] text-orange-700">{g.loeuf!.toFixed(2)}</span>, '#fdba74'))}
            {!intolerant.length && !between.length && <div className="text-[11px] text-slate-500">{s.constraintLoading ? 'Loading gnomAD constraint…' : `No constrained coding gene ${where}.`}</div>}
          </>}
        </div>
      )}
    </div>
  );
}
