/**
 * Repeat inspector: the sizes of a tandem repeat in each sample's reads (repeatScan.ts), as a histogram over the size
 * categories of the locus and, on the same axis, a waterfall of reads drawn unit by unit (pathogenic motif, benign or
 * reference motif, known interruption, other unit), with the alleles called from the modes.
 *
 * Sizes come from the reads spanning the repeat (both flanks) when there are SPAN_MIN of them or more. Below that the
 * reads stopping inside the repeat are added at the length they reach, a lower bound, and the panel says so.
 * Dragging across the plot selects a size range: its reads, their spread and interruptions.
 */
import { useEffect, useMemo, useRef, useState } from 'react';
import type { SashimiDataSource } from './datasource';
import { Segmented } from './controls';
import {
  FLANK, SPAN_MIN, callAlleles, categoryOf, interruptionPattern, measureReads, quantile, revComp,
  type Allele, type CategoryTone, type ReadRepeat, type RepeatLocus, type RepeatReads,
} from '../../standalone/repeatScan';

/** reads asked of a sample over the locus (every k-th kept past it) */
const INSPECT_MAX_READS = 20000;
/** reads drawn in the waterfall, picked evenly by size */
const WATERFALL_ROWS = 160;
const UNIT_COLOR: Record<string, string> = { P: '#6f8fd6', B: '#0f9488', I: '#e07b1f', o: '#a21caf' };
const TONE: Record<CategoryTone, string> = { normal: '#0ca30c', intermediate: '#fab219', premutation: '#ec835a', reduced: '#ec835a', pathogenic: '#d03b3b' };

export interface InspectorSample { id: number; name: string; color: string }
interface SampleResult { reads: RepeatReads; sampled: boolean; total: number }

const fmt = (n: number) => n.toLocaleString('en-US');
const pct = (a: number, b: number) => (b ? `${Math.round((100 * a) / b)} %` : '—');
/** Sample names without what they all share (a run's prefix, a flow cell's suffix), cut at a separator: barcode66, barcode67. */
function shortNames(names: string[]): string[] {
  if (names.length < 2) return names;
  const sep = /[_\-.]/;
  let pre = names[0];
  for (const n of names) while (!n.startsWith(pre)) pre = pre.slice(0, -1);
  let suf = names[0];
  for (const n of names) while (!n.endsWith(suf)) suf = suf.slice(1);
  // back to a separator, so that names keep whole words
  while (pre && !sep.test(pre[pre.length - 1])) pre = pre.slice(0, -1);
  while (suf && !sep.test(suf[0])) suf = suf.slice(1);
  return names.map(n => { const s = n.slice(pre.length, n.length - suf.length); return s.length >= 2 ? s : n; });
}
const niceStep = (span: number) => [1, 2, 5, 10, 20, 25, 50, 100, 200, 500].find(s => span / s <= 90) ?? 1000;

export function RepeatInspector({ locus, samples, ds, initialSample, onClose }: {
  locus: RepeatLocus; samples: InspectorSample[]; ds: SashimiDataSource; initialSample?: number; onClose: () => void;
}) {
  const [sid, setSid] = useState<number>(initialSample ?? samples[0]?.id ?? 0);
  const [results, setResults] = useState<Record<number, { loading?: boolean; error?: string; r?: SampleResult }>>({});
  const [brush, setBrush] = useState<[number, number] | null>(null);
  const [hover, setHover] = useState<{ x: number; y: number; text: string } | null>(null);
  const drag = useRef<number | null>(null);
  const motifs = locus.catalog;
  const geneMotif = motifs?.strand === '-' ? revComp(locus.motif) : null;

  // the reads of the sample over the locus, measured once
  useEffect(() => {
    if (!sid || results[sid]) return;
    let cancelled = false;
    setResults(p => ({ ...p, [sid]: { loading: true } }));
    (async () => {
      const a = locus.start - FLANK - 300, b = locus.end + FLANK + 300;
      const resp = await ds.getReads(sid, locus.chrom, a, b, false, INSPECT_MAX_READS, 'reads', 0, 1, {});
      const reads = resp.reads;
      let lo = Math.min(a, ...reads.map(r => r.s)), hi = Math.max(b, ...reads.map(r => r.e));
      lo = Math.max(0, lo);
      let ref = resp.reference && resp.reference.start <= lo && resp.reference.start + resp.reference.seq.length >= hi ? resp.reference : null;
      if (!ref) {
        const s = await ds.getReference(locus.chrom, lo, hi).catch(() => null);
        ref = s ? { start: lo, seq: s } : resp.reference;
      }
      const measured = measureReads(reads, ref?.seq ?? '', ref?.start ?? lo, locus);
      if (!cancelled) setResults(p => ({ ...p, [sid]: { r: { reads: measured, sampled: resp.total > reads.length, total: resp.total } } }));
    })().catch(e => { if (!cancelled) setResults(p => ({ ...p, [sid]: { error: e?.message ?? String(e) } })); });
    return () => { cancelled = true; };
  }, [sid, locus]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => setBrush(null), [sid]);
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') { if (brush) setBrush(null); else onClose(); } };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [brush, onClose]);

  const cur = results[sid];
  const view = useMemo(() => {
    const r = cur?.r;
    if (!r) return null;
    const enough = r.reads.spanning.length >= SPAN_MIN;
    const used: ReadRepeat[] = enough ? r.reads.spanning : [...r.reads.spanning, ...r.reads.truncated];
    const units = used.map(x => x.units).sort((a, b) => a - b);
    const alleles = callAlleles(units);
    const spanSorted = r.reads.spanning.map(x => x.units).sort((a, b) => a - b);
    const p95 = quantile(spanSorted, 0.95);
    const beyond = r.reads.truncated.filter(x => spanSorted.length && x.units > p95).length;
    const xMaxRaw = Math.max(quantile(units, 0.99) * 1.12, locus.refUnits * 2.5, 30);
    const step = niceStep(xMaxRaw);
    const xMax = Math.ceil(xMaxRaw / step) * step;
    // waterfall rows: evenly by size
    const order = [...used].sort((a, b) => a.units - b.units);
    const rows = order.length <= WATERFALL_ROWS ? order : Array.from({ length: WATERFALL_ROWS }, (_, i) => order[Math.round((i * (order.length - 1)) / (WATERFALL_ROWS - 1))]);
    return { r, enough, used, units, alleles, beyond, p95, xMax, step, rows };
  }, [cur, locus.refUnits]);

  const W = 1000, L = 64, R = 24, top = 46, hH = 150;
  const x = (u: number) => L + (Math.min(u, view?.xMax ?? 1) / (view?.xMax ?? 1)) * (W - L - R);
  const ux = (px: number) => Math.max(0, ((px - L) / (W - L - R)) * (view?.xMax ?? 1));

  const svgPoint = (e: React.MouseEvent<SVGSVGElement>) => {
    const b = e.currentTarget.getBoundingClientRect();
    return { px: ((e.clientX - b.left) / b.width) * W, py: ((e.clientY - b.top) / b.height) * (e.currentTarget.viewBox.baseVal.height) };
  };

  const stats = (list: ReadRepeat[]) => {
    const u = list.map(r => r.units).sort((a, b) => a - b);
    const pats = new Map<string, number>();
    for (const r of list) { const p = interruptionPattern(r.tokens); pats.set(p, (pats.get(p) ?? 0) + 1); }
    const top3 = [...pats].sort((a, b) => b[1] - a[1]).slice(0, 3);
    return { n: list.length, median: quantile(u, 0.5), p5: quantile(u, 0.05), p95: quantile(u, 0.95), pats: top3 };
  };
  const patText = (p: string) => (p ? `${locus.interruptions[0] ?? 'interruption'} at ${p}` : `no ${locus.interruptions.length ? locus.interruptions[0] : 'interruption'}`);
  const alleleReads = (a: Allele) => (view ? view.used.filter(r => r.units > a.lo && r.units <= a.hi) : []);

  const body = (() => {
    if (!cur || cur.loading) return <div className="px-5 py-10 text-center text-sm text-slate-500">Reading the reads over the repeat…</div>;
    if (cur.error) return <div className="px-5 py-6 text-sm text-red-700">Could not read this sample: {cur.error}</div>;
    if (!view) return null;
    const { r, enough, alleles, xMax, step, rows, units, beyond, p95 } = view;
    if (!units.length) return <div className="px-5 py-6 text-sm text-slate-600">No read of this sample reaches the repeat with a flank on either side ({fmt(r.reads.total)} reads over the locus).</div>;
    const nb = Math.ceil(xMax / step);
    const span = new Array(nb).fill(0), trunc = new Array(nb).fill(0);
    for (const u of r.reads.spanning) span[Math.min(nb - 1, Math.floor(u.units / step))]++;
    if (!enough) for (const u of r.reads.truncated) trunc[Math.min(nb - 1, Math.floor(u.units / step))]++;
    const yMax = Math.max(1, ...span.map((v, i) => v + trunc[i]));
    const base = top + hH, wTop = base + 34;
    const rowH = Math.max(1.4, Math.min(4, 320 / Math.max(1, rows.length))), wH = rows.length * rowH;
    const H = wTop + wH + 10;
    const sel = brush ? view.used.filter(u => u.units >= Math.min(...brush) && u.units <= Math.max(...brush)) : null;
    const selStats = sel ? stats(sel) : null;
    return (
      <>
        {!enough && (
          <div className="mx-5 mt-3 rounded-lg border border-amber-300 bg-amber-50 px-3 py-2 text-[12px] text-amber-900" role="alert">
            <b>Few reads span the repeat:</b> {fmt(r.reads.spanning.length)} read{r.reads.spanning.length === 1 ? '' : 's'} with both flanks (at least {SPAN_MIN} are needed to size it from them alone).
            The {fmt(r.reads.truncated.length)} reads that stop inside the repeat are added at the length they reach, hatched: <b>lower bounds</b>, the alleles may be longer than shown.
          </div>
        )}
        {enough && beyond > 0 && beyond >= 0.1 * r.reads.truncated.length && (
          <div className="mx-5 mt-3 rounded-lg border border-amber-200 bg-amber-50/60 px-3 py-1.5 text-[11.5px] text-amber-900">
            {fmt(beyond)} reads stop inside the repeat beyond the longest spanning ones ({'>'} {p95} units): long alleles are under-represented among the reads that span it.
          </div>
        )}
        <div className="flex flex-wrap gap-2 px-5 pt-3">
          <Chip k="spanning reads (both flanks)" v={fmt(r.reads.spanning.length)} />
          <Chip k={enough ? 'stop inside the repeat · not counted' : 'stop inside the repeat · lower bounds'} v={fmt(r.reads.truncated.length)} />
          {alleles.map((a, i) => {
            const c = categoryOf(locus.categories, a.mode), st = stats(alleleReads(a));
            return (
              <Chip key={i} k={`allele ${i + 1} · ${fmt(a.n)} reads`} accent={c ? TONE[c.tone] : undefined}
                v={<><b>{a.mode}</b> <span className="font-normal text-slate-600">{locus.motif}</span>{c && <span className="ml-1.5 text-[11px] font-semibold px-1.5 py-0.5 rounded-full border" style={{ borderColor: TONE[c.tone], color: '#334155' }}>{c.name}</span>}</>}
                sub={`P5–P95 ${a.p5}–${a.p95}${a.broad ? ' · broad' : ''}${locus.interruptions.length && st.pats[0] ? ` · ${patText(st.pats[0][0])} (${pct(st.pats[0][1], st.n)})` : ''}`} />
            );
          })}
          {alleles.some(a => a.broad) && <Chip k="reading" v={<span className="text-red-700">broad spread</span>} sub="somatic mosaicism, or PCR / sequencing stutter" />}
        </div>
        <div className="px-3 pt-1 relative">
          <svg viewBox={`0 0 ${W} ${H}`} width="100%" style={{ maxWidth: W, cursor: 'crosshair', userSelect: 'none' }} fontFamily="Inter, system-ui, sans-serif"
            onMouseDown={e => { const { px } = svgPoint(e); if (px < L) return; drag.current = ux(px); setBrush([ux(px), ux(px)]); }}
            onMouseMove={e => {
              const { px, py } = svgPoint(e);
              if (drag.current != null) { setBrush([drag.current, ux(px)]); return; }
              if (py > base || px < L || px > W - R) { setHover(null); return; }
              const i = Math.floor(ux(px) / step);
              if (i < 0 || i >= nb) { setHover(null); return; }
              setHover({ x: px, y: py, text: `${i * step}–${(i + 1) * step - 1} ${locus.motif}: ${fmt(span[i])} spanning${!enough && trunc[i] ? ` + ${fmt(trunc[i])} at least` : ''}` });
            }}
            onMouseUp={() => { drag.current = null; if (brush && Math.abs(brush[1] - brush[0]) < step / 2) setBrush(null); }}
            onMouseLeave={() => { drag.current = null; setHover(null); }}>
            <defs>
              <pattern id="ri-hatch" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
                <rect width="5" height="5" fill="#dbe3f5" /><line x1="0" y1="0" x2="0" y2="5" stroke="#6f8fd6" strokeWidth="2" />
              </pattern>
            </defs>
            {/* size categories across both panels: a label inside its band when it fits, else above it, pushed right of the previous one */}
            {(() => {
              let aboveRight = -Infinity;
              return locus.categories.map(c => {
                if (c.min > xMax) return null;
                const a = x(Math.max(0, c.min - 0.5)), b = x(Math.min(xMax, c.max + 0.5));
                const label = `${c.name} ${c.max === Infinity ? `> ${c.min - 1}` : c.min === 0 ? `≤ ${c.max}` : `${c.min}–${c.max}`}`, tw = label.length * 6.3;
                const inside = tw + 8 <= b - a;
                let lx = inside ? a + 4 : Math.max((a + b) / 2 - tw / 2, aboveRight + 8);
                if (!inside) aboveRight = lx + tw;
                lx = Math.min(lx, W - R - tw);
                return (
                  <g key={c.name}>
                    <rect x={a} y={top - 18} width={b - a} height={H - top + 8} fill={TONE[c.tone]} opacity={0.07} />
                    <rect x={a} y={top - 18} width={b - a} height={3} fill={TONE[c.tone]} />
                    {!inside && <line x1={(a + b) / 2} y1={top - 18} x2={Math.min(Math.max((a + b) / 2, lx), lx + tw)} y2={top - 22} stroke={TONE[c.tone]} />}
                    <text x={lx} y={inside ? top - 5 : top - 25} fontSize={11} fontWeight={600} fill="#334155">{label}</text>
                  </g>
                );
              });
            })()}
            {/* histogram */}
            {span.map((v, i) => {
              const a = x(i * step) + 0.5, w = Math.max(0.8, x((i + 1) * step) - x(i * step) - 1.5), hs = (v / yMax) * (hH - 12), ht = (trunc[i] / yMax) * (hH - 12);
              return (
                <g key={i}>
                  {v > 0 && <rect x={a} y={base - hs} width={w} height={hs} rx={1.5} fill={UNIT_COLOR.P} />}
                  {trunc[i] > 0 && <rect x={a} y={base - hs - ht} width={w} height={ht} fill="url(#ri-hatch)" stroke="#6f8fd6" strokeWidth={0.5} />}
                </g>
              );
            })}
            <line x1={L} y1={base} x2={W - R} y2={base} stroke="#94a3b8" />
            <text x={L - 8} y={base - hH + 16} textAnchor="end" fontSize={10} fill="#64748b">{fmt(yMax)}</text>
            <text x={L - 8} y={base} textAnchor="end" fontSize={10} fill="#64748b">0</text>
            <text x={L - 8} y={base - hH / 2} textAnchor="end" fontSize={10} fill="#64748b">reads</text>
            {Array.from({ length: Math.floor(xMax / (step * 10)) + 1 }, (_, i) => i * step * 10).map(u => (
              <g key={u}>
                <line x1={x(u)} y1={base} x2={x(u)} y2={base + 4} stroke="#94a3b8" />
                <text x={x(u)} y={base + 15} textAnchor="middle" fontSize={10.5} fill="#475569">{u}</text>
              </g>
            ))}
            <text x={W - R} y={base + 28} textAnchor="end" fontSize={10.5} fill="#64748b">{locus.motif} units (tract length / {locus.k}) · one axis for both panels</text>
            {/* reference and alleles */}
            <line x1={x(locus.refUnits)} y1={top} x2={x(locus.refUnits)} y2={H - 6} stroke="#334155" strokeDasharray="4 3" />
            <text x={x(locus.refUnits) + 4} y={top + 10} fontSize={10.5} fontWeight={600} fill="#334155">reference {locus.refUnits}</text>
            {view.alleles.map((a, i) => (
              <g key={i}>
                <line x1={x(a.mode)} y1={top + 14} x2={x(a.mode)} y2={H - 6} stroke="#1e3a8a" strokeWidth={1.3} strokeDasharray="6 3" />
                <text x={x(a.mode) + 4} y={top + 26 + i * 13} fontSize={10.5} fontWeight={700} fill="#1e3a8a">allele {i + 1} · {a.mode}</text>
              </g>
            ))}
            {/* waterfall */}
            <text x={L - 8} y={wTop + 8} textAnchor="end" fontSize={10} fill="#64748b">reads</text>
            <text x={L - 8} y={wTop + 20} textAnchor="end" fontSize={10} fill="#64748b">by size</text>
            {rows.map((r, k) => {
              const y = wTop + k * rowH, runs: JSX.Element[] = [];
              const t = r.tokens;
              for (let i = 0; i < t.length;) {
                let j = i; while (j < t.length && t[j] === t[i]) j++;
                if (i < xMax) runs.push(<rect key={i} x={x(i)} y={y} width={Math.max(t[i] === 'P' ? 0.4 : 1.4, x(Math.min(j, xMax)) - x(i))} height={rowH - 0.35} fill={UNIT_COLOR[t[i]] ?? UNIT_COLOR.o} />);
                i = j;
              }
              // the read's size: its units may differ from its tokens by an indel's worth
              return (
                <g key={k} opacity={r.truncated ? 0.55 : 1}>
                  {runs}
                  {r.truncated && <path d={`M${x(r.units)},${y}l3,${rowH / 2}l-3,${rowH / 2}`} fill="none" stroke="#475569" strokeWidth={0.8} />}
                </g>
              );
            })}
            {brush && (
              <rect x={x(Math.min(...brush))} y={top - 18} width={Math.max(1, x(Math.max(...brush)) - x(Math.min(...brush)))} height={H - top + 8} fill="#4f46e5" opacity={0.1} stroke="#4f46e5" strokeDasharray="4 2" pointerEvents="none" />
            )}
          </svg>
          {hover && !drag.current && (
            <div className="pointer-events-none absolute z-10 rounded-md border border-slate-200 bg-white/95 px-2 py-1 text-[11px] shadow" style={{ left: `calc(${(hover.x / W) * 100}% + 18px)`, top: hover.y }}>{hover.text}</div>
          )}
        </div>
        <div className="px-5 pb-1 text-[11.5px] text-slate-600 min-h-[22px]">
          {selStats && brush ? (
            <span>
              <b>{Math.round(Math.min(...brush))}–{Math.round(Math.max(...brush))} units:</b> {fmt(selStats.n)} reads ({pct(selStats.n, view.used.length)}) · median {selStats.median} · P5–P95 {selStats.p5}–{selStats.p95}
              {locus.interruptions.length > 0 && selStats.pats.length > 0 && <> · {selStats.pats.map(([p, n]) => `${patText(p)} ${pct(n, selStats.n)}`).join(', ')}</>}
              {locus.categories.length > 0 && <> · {locus.categories.map(c => `${c.name.toLowerCase()} ${pct(sel!.filter(u => u.units >= c.min && u.units <= c.max).length, selStats.n)}`).join(', ')}</>}
              <button onClick={() => setBrush(null)} className="ml-2 text-indigo-700 underline">clear</button>
            </span>
          ) : <span className="text-slate-400">Drag across the plot to measure a size range: its reads, spread and interruptions.</span>}
        </div>
        <div className="flex flex-wrap items-center gap-x-4 gap-y-1 px-5 pt-1 text-[11.5px] text-slate-600">
          <Swatch c={UNIT_COLOR.P} t={`${locus.pathogenic.join(' / ')}${locus.benign.length ? ' (pathogenic motif)' : ''}`} />
          {locus.benign.length > 0 && <Swatch c={UNIT_COLOR.B} t={`${locus.benign.join(' / ')} (benign or reference motif)`} />}
          {locus.interruptions.length > 0 && <Swatch c={UNIT_COLOR.I} t={`${locus.interruptions.join(' / ')} (interruption)`} />}
          <Swatch c={UNIT_COLOR.o} t="other unit (variant or sequencing error)" />
          {!enough && <span className="inline-flex items-center gap-1.5"><svg width="14" height="10"><rect width="14" height="10" fill="url(#ri-hatch)" stroke="#6f8fd6" strokeWidth="0.5" /></svg>stops inside the repeat: at least that long</span>}
        </div>
        <div className="px-5 pt-2 pb-4 text-[11px] leading-relaxed text-slate-500">
          Each read is rebuilt over the locus from its soft clips, aligned bases and insertions; its repeat is what lies between the {FLANK}-bp flanks of the reference
          (found wherever they are in the read), its size the tract length divided by {locus.k}, so that sequencing indels inside it average out. Primary alignments only;
          {r.sampled ? ` ${fmt(INSPECT_MAX_READS)} of the ${fmt(r.total)} reads over the locus (every k-th kept);` : ` ${fmt(r.reads.total)} reads over the locus;`} {fmt(r.reads.skipped)} do not reach the repeat with a flank.
          Alleles: modes of a kernel density on the log of the size, at most two; "broad" when P10–P90 exceeds a quarter of the median.
          {locus.categorySource && <> Categories: {locus.categorySource}.</>} PCR amplification and nanopore sequencing both add stutter to GC-rich repeats: a broad allele is not by itself mosaicism.
          Evidence from the reads, not a diagnostic call.
        </div>
      </>
    );
  })();

  return (
    <div className="fixed inset-0 z-50 bg-black/35 flex items-start justify-center p-4 overflow-y-auto" onMouseDown={onClose}>
      <div className="bg-white rounded-xl shadow-2xl border border-slate-200 w-full max-w-[1040px] text-slate-900 mt-6" onMouseDown={e => e.stopPropagation()} role="dialog" aria-label="Repeat inspector">
        <div className="flex items-start justify-between gap-4 px-5 pt-4 pb-3 border-b border-slate-200">
          <div className="min-w-0">
            <div className="font-bold text-[15px]">Repeat inspector · {locus.label} · {locus.chrom}:{fmt(locus.start + 1)}-{fmt(locus.end)}</div>
            <div className="text-[12px] text-slate-500 mt-0.5">
              Reference: {locus.refStructure} = {locus.refUnits} units{geneMotif ? ` (${geneMotif} on ${motifs!.gene}'s strand)` : ''}
              {motifs && <> · {motifs.disease} ({motifs.inh}) · {motifs.where}</>}
            </div>
          </div>
          <div className="flex items-center gap-3 shrink-0">
            {samples.length > 1 && (
              <Segmented size="sm" label="Sample" title="The sample whose reads are measured" value={String(sid)} onChange={v => setSid(Number(v))}
                options={samples.map((s, i) => ({ value: String(s.id), label: shortNames(samples.map(x => x.name))[i], hint: s.name, dot: s.color }))} />
            )}
            <button onClick={onClose} className="text-slate-400 hover:text-slate-700 text-xl leading-none" title="Close (Esc)">×</button>
          </div>
        </div>
        {body}
      </div>
    </div>
  );
}

function Chip({ k, v, sub, accent }: { k: string; v: React.ReactNode; sub?: string; accent?: string }) {
  return (
    <div className="rounded-[10px] border border-slate-200 bg-slate-50 px-3 py-1.5 min-w-[150px]" style={accent ? { borderLeft: `3px solid ${accent}` } : undefined}>
      <div className="text-[10.5px] text-slate-500">{k}</div>
      <div className="text-[15px] font-semibold leading-tight">{v}</div>
      {sub && <div className="text-[10.5px] text-slate-500 mt-0.5">{sub}</div>}
    </div>
  );
}
function Swatch({ c, t }: { c: string; t: string }) {
  return <span className="inline-flex items-center gap-1.5"><span className="inline-block w-3 h-2.5 rounded-[2px]" style={{ background: c }} />{t}</span>;
}
