/**
 * The sequences of a repeat inspector's reads (RepeatInspector.tsx): for each read its 5′ flank aligned base to base on
 * the reference's, its tract as runs of units ("(CGG)10 AGG (CGG)9 …", every unit on demand) and its 3′ flank, with the
 * anchors underlined; reads that stop inside the repeat show what follows the stretch measured, reads set apart for
 * their anchors the whole read with every anchor found in it.
 */
import { useEffect, useMemo, useState } from 'react';
import { ANCHOR, FLANK_SHOWN, alignFlank, unitRuns, type ReadRepeat, type RepeatLocus } from '../../standalone/repeatScan';
import { CONSENSUS_LOW, CONSENSUS_MAX_READS, consensusFasta, consensusSequence, repeatConsensus, type RepeatConsensus } from '../../standalone/repeatConsensus';

const PAGE = 15;
const UNIT_TEXT: Record<string, string> = { P: '#3b5bb5', B: '#0f766e', I: '#c2410c', o: '#a21caf', x: '#94a3b8' };
const fmt = (n: number) => n.toLocaleString('en-US');

export function ReadSequences({ locus, title, reads, sample, onClose }: { locus: RepeatLocus; title: string; reads: ReadRepeat[]; sample: string; onClose: () => void }) {
  const [desc, setDesc] = useState(true);
  const [page, setPage] = useState(0);
  /** the consensus of the reads listed, once asked for */
  const [cons, setCons] = useState<{ busy: boolean; c?: RepeatConsensus | null } | null>(null);
  useEffect(() => { setCons(null); setPage(0); }, [reads]);
  const makeConsensus = () => {
    setCons({ busy: true });
    // after the "computing" state is drawn: a few hundred reads take a fraction of a second
    setTimeout(() => setCons({ busy: false, c: repeatConsensus(reads, locus) }), 30);
  };
  const sorted = useMemo(() => [...reads].sort((a, b) => (desc ? b.units - a.units : a.units - b.units)), [reads, desc]);
  const pages = Math.max(1, Math.ceil(sorted.length / PAGE));
  const shown = sorted.slice(page * PAGE, page * PAGE + PAGE);
  return (
    <div className="mx-5 mb-4 rounded-lg border border-slate-200">
      <div className="flex flex-wrap items-center gap-2 px-3 py-2 border-b border-slate-200 bg-slate-50 rounded-t-lg text-[12px]">
        <b>Read sequences · {title}</b>
        <button onClick={makeConsensus} disabled={cons?.busy} className="rounded-md bg-indigo-600 px-2.5 py-0.5 font-semibold text-white hover:bg-indigo-700 disabled:opacity-50"
          title={`Consensus of these reads: flanks base by base on the reference, the tract unit by unit (aligned on the read of median size), each position with its support${reads.length > CONSENSUS_MAX_READS ? `; ${CONSENSUS_MAX_READS} of the ${fmt(reads.length)} reads, evenly by size` : ''}`}>
          {cons?.busy ? 'Computing…' : 'Consensus'}
        </button>
        <span className="text-slate-500">{fmt(reads.length)} read{reads.length === 1 ? '' : 's'}</span>
        {reads.length > 1 && (
          <button onClick={() => { setDesc(d => !d); setPage(0); }} className="ml-2 rounded border border-slate-300 bg-white px-2 py-0.5 hover:bg-slate-100">{desc ? 'longest first ↓' : 'shortest first ↑'}</button>
        )}
        {pages > 1 && (
          <span className="flex items-center gap-1 ml-2">
            <button disabled={page === 0} onClick={() => setPage(p => p - 1)} className="rounded border border-slate-300 bg-white px-1.5 disabled:opacity-40">‹</button>
            <span className="tabular-nums text-slate-600">{fmt(page * PAGE + 1)}–{fmt(Math.min(sorted.length, page * PAGE + PAGE))} of {fmt(sorted.length)}</span>
            <button disabled={page >= pages - 1} onClick={() => setPage(p => p + 1)} className="rounded border border-slate-300 bg-white px-1.5 disabled:opacity-40">›</button>
          </span>
        )}
        <span className="ml-auto flex items-center gap-3 text-[11px] text-slate-500">
          <span><span className="px-1 rounded bg-red-100 text-red-800 font-mono">T</span> differs from the reference</span>
          <span><span className="px-1 rounded bg-amber-100 text-amber-900 font-mono">A</span> extra base in the read</span>
          <span><span className="font-mono text-red-700">-</span> base missing</span>
          <span><span className="font-mono underline decoration-2 decoration-indigo-500">anchor</span></span>
          <button onClick={onClose} className="text-slate-400 hover:text-slate-700 text-base leading-none" title="Close the sequences">×</button>
        </span>
      </div>
      {cons && !cons.busy && (cons.c ? <ConsensusCard c={cons.c} locus={locus} sample={sample} title={title} onClose={() => setCons(null)} />
        : <div className="px-3 py-2 text-[12px] text-slate-500 border-b border-slate-200">No consensus: none of these reads has a measured tract.</div>)}
      <div className="divide-y divide-slate-100">
        {shown.map(r => <ReadCard key={r.name} r={r} locus={locus} />)}
        {!shown.length && <div className="px-3 py-4 text-[12px] text-slate-500">No read here.</div>}
      </div>
    </div>
  );
}

/** One read: its flanks against the reference's, and its tract. */
function ReadCard({ r, locus }: { r: ReadRepeat; locus: RepeatLocus }) {
  const [bases, setBases] = useState(false);
  const [a, b] = r.tract;
  const chimeric = !!r.note;
  const leftAnchored = !r.truncated || r.from === 'left';
  const rightAnchored = !r.truncated || r.from === 'right';
  const tract = r.seq.slice(a, b);
  const runs = useMemo(() => (chimeric ? [] : unitRuns(tract, locus)), [tract, locus, chimeric]);
  const left = useMemo(() => {
    if (chimeric || !leftAnchored) return null;
    const from = Math.max(0, a - FLANK_SHOWN - 15), part = r.seq.slice(from, a);
    return { ...alignFlank(part, locus.refL, true), from, anchor: r.hitsL[0] != null ? [r.hitsL[0] - ANCHOR, r.hitsL[0]] as [number, number] : null };
  }, [r, a, locus, chimeric, leftAnchored]);
  const right = useMemo(() => {
    if (chimeric || !rightAnchored) return null;
    const part = r.seq.slice(b, b + FLANK_SHOWN + 15);
    return { ...alignFlank(part, locus.refR, false), from: b, anchor: r.hitsR[0] != null ? [r.hitsR[0], r.hitsR[0] + ANCHOR] as [number, number] : null };
  }, [r, b, locus, chimeric, rightAnchored]);
  const ed = (v?: number) => (v == null ? '—' : `${v} edit${v === 1 ? '' : 's'}`);
  const k = locus.k;
  return (
    <div className="px-3 py-2 text-[11.5px]">
      <div className="flex flex-wrap items-baseline gap-x-3 text-slate-600">
        <b className="text-slate-900">{r.name}</b>
        <span>{r.reverse ? 'reverse strand' : 'forward strand'}</span>
        {!chimeric && <span><b className="text-slate-900">{r.truncated ? '≥ ' : ''}{r.units}</b> {locus.motif} ({fmt(b - a)} bp)</span>}
        <span>5′ anchor {r.hitsL.length > 1 ? `×${r.hitsL.length}` : ed(r.edL)}</span>
        <span>3′ anchor {r.hitsR.length > 1 ? `×${r.hitsR.length}` : ed(r.edR)}</span>
        {!chimeric && <span>{Math.round(r.purity * 100)} % motif units</span>}
        {r.truncated && <span className="text-amber-700">stops inside the repeat ({r.from === 'left' ? 'no 3′ anchor' : r.from === 'right' ? 'no 5′ anchor' : 'no anchor'})</span>}
        {r.note && <span className="text-red-700">set apart: {r.note}</span>}
        <span className="text-slate-400">read {fmt(r.seq.length)} bp</span>
      </div>
      {chimeric ? (
        <WholeRead r={r} />
      ) : (
        <div className="mt-1 grid grid-cols-[72px_1fr] gap-x-2 gap-y-0.5 items-start">
          <span className="text-slate-500 pt-px">5′ flank</span>
          {left ? <Aligned al={left} /> : <Raw text={r.seq.slice(Math.max(0, a - 60), a)} note={a ? 'the read before the stretch measured (no 5′ anchor)' : 'the read starts here'} />}
          <span className="text-slate-500 pt-px">repeat</span>
          <div className="font-mono leading-5 break-all">
            {runs.map((u, i) => (
              <span key={i} style={{ color: UNIT_TEXT[u.cls] }} className={u.cls === 'x' ? 'text-[10px] lowercase' : 'font-semibold'}>
                {u.cls === 'x' ? u.unit.toLowerCase() : u.n > 1 ? `(${u.unit})${u.n}` : u.unit}{' '}
              </span>
            ))}
            <button onClick={() => setBases(v => !v)} className="ml-1 font-sans text-[11px] text-indigo-700 underline">{bases ? 'hide units' : 'every unit'}</button>
            {bases && <Units runs={runs} k={k} />}
          </div>
          <span className="text-slate-500 pt-px">3′ flank</span>
          {right ? <Aligned al={right} /> : <Raw text={r.seq.slice(b, b + 60)} note={b < r.seq.length ? 'the read after the stretch measured (no 3′ anchor)' : 'the read ends here'} />}
        </div>
      )}
    </div>
  );
}

/** A flank against the reference's: two rows, differences marked, the anchor underlined in the read. */
function Aligned({ al }: { al: { ref: string; read: string; from: number; anchor: [number, number] | null } }) {
  // read positions of the columns, to underline the anchor
  // every base of the read's part is in the row, in order: the first one is at al.from
  const cols: { r: string; q: string; pos: number | null }[] = [];
  let pos = al.from;
  for (let i = 0; i < al.read.length; i++) {
    const q = al.read[i], rr = al.ref[i];
    const here = q !== '-' && q !== ' ' ? pos++ : null;
    cols.push({ r: rr, q, pos: here });
  }
  const edits = cols.filter(c => c.q !== ' ' && c.r !== ' ' && c.q !== c.r).length;
  return (
    <div className="font-mono text-[11px] leading-[15px] overflow-x-auto whitespace-pre">
      <div className="text-slate-400">{cols.map((c, i) => <span key={i}>{c.r}</span>)} <span className="font-sans text-[10px]">reference</span></div>
      <div>
        {cols.map((c, i) => {
          const anchor = al.anchor && c.pos != null && c.pos >= al.anchor[0] && c.pos < al.anchor[1];
          const cls = c.q === ' ' ? '' : c.q === '-' ? 'text-red-700' : c.r === '-' ? 'bg-amber-100 text-amber-900' : c.r !== ' ' && c.r !== c.q ? 'bg-red-100 text-red-800' : 'text-slate-800';
          return <span key={i} className={`${cls} ${anchor ? 'underline decoration-2 decoration-indigo-500 underline-offset-2' : ''}`}>{c.q}</span>;
        })} <span className="font-sans text-[10px] text-slate-500">read · {edits} difference{edits === 1 ? '' : 's'}</span>
      </div>
    </div>
  );
}

function Raw({ text, note }: { text: string; note: string }) {
  return <div className="font-mono text-[11px] leading-[15px] break-all text-slate-500">{text || '—'} <span className="font-sans text-[10px] text-slate-400">{note}</span></div>;
}

/** Every unit of the tract, numbered every 10. */
function Units({ runs, k }: { runs: { unit: string; cls: string; n: number }[]; k: number }) {
  const units: { u: string; cls: string }[] = [];
  for (const r of runs) for (let i = 0; i < (r.cls === 'x' ? 1 : r.n); i++) units.push({ u: r.unit, cls: r.cls });
  let n = 0;
  return (
    <div className="mt-1 font-mono text-[11px] leading-5 break-all">
      {units.map((x, i) => {
        if (x.cls !== 'x') n++;
        return (
          <span key={i}>
            {x.cls !== 'x' && (n - 1) % 10 === 0 && <sup className="text-[8px] text-slate-400 mr-px">{n}</sup>}
            <span style={{ color: UNIT_TEXT[x.cls] }} className={x.cls === 'x' ? 'text-[10px]' : ''}>{x.cls === 'x' ? x.u.toLowerCase() : x.u}</span>
            {x.cls !== 'x' && <span className="text-slate-300">{k > 1 ? '·' : ''}</span>}
          </span>
        );
      })}
    </div>
  );
}

/** A read set apart for its anchors: the whole read, every anchor found underlined (5′ indigo, 3′ teal). */
function WholeRead({ r }: { r: ReadRepeat }) {
  const max = 3000, seq = r.seq.slice(0, max);
  const inL = (p: number) => r.hitsL.some(h => p >= h - ANCHOR && p < h), inR = (p: number) => r.hitsR.some(h => p >= h && p < h + ANCHOR);
  const spans: JSX.Element[] = [];
  let i = 0;
  while (i < seq.length) {
    const kind = inL(i) ? 'L' : inR(i) ? 'R' : '';
    let j = i; while (j < seq.length && (inL(j) ? 'L' : inR(j) ? 'R' : '') === kind) j++;
    spans.push(<span key={i} className={kind === 'L' ? 'bg-indigo-100 text-indigo-900' : kind === 'R' ? 'bg-teal-100 text-teal-900' : 'text-slate-600'}>{seq.slice(i, j)}</span>);
    i = j;
  }
  return (
    <div className="mt-1">
      <div className="text-[10.5px] text-slate-500">5′ anchor at {r.hitsL.map(h => fmt(h - ANCHOR)).join(', ') || 'none'} · 3′ anchor at {r.hitsR.map(h => fmt(h)).join(', ') || 'none'} (read positions; <span className="bg-indigo-100 px-0.5">5′</span> <span className="bg-teal-100 px-0.5">3′</span>)</div>
      <div className="font-mono text-[10.5px] leading-[14px] break-all max-h-48 overflow-y-auto">{spans}{r.seq.length > max ? ` … (+${fmt(r.seq.length - max)} bp)` : ''}</div>
    </div>
  );
}

/** The consensus of the reads listed: a summary line, its flanks against the reference, its tract by units, its support. */
function ConsensusCard({ c, locus, sample, title, onClose }: { c: RepeatConsensus; locus: RepeatLocus; sample: string; title: string; onClose: () => void }) {
  const [copied, setCopied] = useState('');
  const left = useMemo(() => ({ ...alignFlank(c.flankL.seq, locus.refL, true), from: 0, anchor: null }), [c, locus]);
  const right = useMemo(() => ({ ...alignFlank(c.flankR.seq, locus.refR, false), from: 0, anchor: null }), [c, locus]);
  const intr = locus.interruptions[0];
  const summary = `${locus.label} ${locus.chrom}:${fmt(locus.start + 1)}-${fmt(locus.end)} · ${sample} · consensus of ${c.used} read${c.used === 1 ? '' : 's'} (${title}): ${c.structure} · ${c.lowerBound ? '≥ ' : ''}${c.size} ${locus.motif}` +
    `${intr ? ` · ${c.interruptions.length ? `${intr} at ${c.interruptions.join(', ')}` : `no ${intr}`}` : ''} · mean support ${Math.round(c.meanSupport * 100)} %`;
  const copy = (text: string, what: string) => { navigator.clipboard?.writeText(text).then(() => { setCopied(what); setTimeout(() => setCopied(''), 1200); }).catch(() => {}); };
  const download = () => {
    const blob = new Blob([consensusFasta(c, locus, sample, title)], { type: 'text/plain' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `${locus.label.replace(/\s+/g, '_')}_${sample.replace(/[^\w.-]+/g, '_')}_consensus.fa`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  // runs of the same unit, a low-support unit apart
  const runs: { unit: string; cls: string; n: number; low: boolean }[] = [];
  for (const u of c.units) {
    const low = u.support < CONSENSUS_LOW, last = runs[runs.length - 1];
    if (last && last.unit === u.unit && last.low === low) last.n++; else runs.push({ unit: u.unit, cls: u.cls, n: 1, low });
  }
  return (
    <div className="px-3 py-2 text-[11.5px] border-b-2 border-indigo-200 bg-indigo-50/40">
      <div className="flex flex-wrap items-baseline gap-x-3 text-slate-600">
        <b className="text-indigo-900">Consensus</b>
        <span>{fmt(c.used)} of the {fmt(c.listed)} reads{c.listed > CONSENSUS_MAX_READS ? ` (${fmt(CONSENSUS_MAX_READS)} at most, evenly by size)` : c.listed > c.used ? ' (reads without a measured tract left out)' : ''}</span>
        <span><b className="text-slate-900">{c.lowerBound ? '≥ ' : ''}{c.size}</b> {locus.motif} · median of the reads {c.medianSize}</span>
        <span>mean support {Math.round(c.meanSupport * 100)} %</span>
        {c.lowerBound && <span className="text-amber-700">from reads stopping inside the repeat: a lower bound</span>}
        <span className="ml-auto flex gap-2">
          <button onClick={() => copy(summary, 'summary')} className="rounded border border-slate-300 bg-white px-2 py-0.5 hover:bg-slate-50">{copied === 'summary' ? 'copied' : 'copy summary'}</button>
          <button onClick={() => copy(consensusSequence(c), 'sequence')} className="rounded border border-slate-300 bg-white px-2 py-0.5 hover:bg-slate-50">{copied === 'sequence' ? 'copied' : 'copy sequence'}</button>
          <button onClick={download} className="rounded border border-slate-300 bg-white px-2 py-0.5 hover:bg-slate-50">FASTA</button>
          <button onClick={onClose} className="text-slate-400 hover:text-slate-700" title="Close the consensus">×</button>
        </span>
      </div>
      <div className="mt-1 font-mono text-[11.5px] text-slate-800 break-all">{summary}</div>
      <div className="mt-1 grid grid-cols-[72px_1fr] gap-x-2 gap-y-0.5 items-start">
        <span className="text-slate-500 pt-px">5′ flank</span>
        {c.flankL.seq ? <Aligned al={left} /> : <Raw text="" note="no read of the list has its 5′ flank" />}
        <span className="text-slate-500 pt-px">repeat</span>
        <div className="font-mono leading-5 break-all">
          {runs.map((u, i) => (
            <span key={i} style={{ color: UNIT_TEXT[u.cls] }} className={u.low ? 'opacity-70' : 'font-semibold'} title={u.low ? `support under ${Math.round(CONSENSUS_LOW * 100)} %` : undefined}>
              {u.n > 1 ? `(${u.low ? u.unit.toLowerCase() : u.unit})${u.n}` : u.low ? u.unit.toLowerCase() : u.unit}{' '}
            </span>
          ))}
        </div>
        <span className="text-slate-500 pt-px">3′ flank</span>
        {c.flankR.seq ? <Aligned al={right} /> : <Raw text="" note="no read of the list has its 3′ flank" />}
        <span className="text-slate-500 pt-px">support</span>
        <SupportTrack c={c} />
      </div>
    </div>
  );
}

/** Support along the consensus (5′ flank, repeat, 3′ flank), with the lowercase threshold. */
function SupportTrack({ c }: { c: RepeatConsensus }) {
  const parts = [
    { label: "5′ flank", sup: c.flankL.support },
    { label: 'repeat', sup: c.units.map(u => u.support) },
    { label: "3′ flank", sup: c.flankR.support },
  ];
  const total = parts.reduce((n, p) => n + p.sup.length, 0) || 1, W = 900, H = 34;
  let x0 = 0;
  return (
    <svg viewBox={`0 0 ${W} ${H + 12}`} width="100%" style={{ maxWidth: W }} className="block">
      {parts.map(p => {
        const w = (p.sup.length / total) * W, start = x0;
        x0 += w;
        const d = p.sup.map((v, i) => `${i ? 'L' : 'M'}${(start + ((i + 0.5) / p.sup.length) * w).toFixed(1)},${(H - v * H).toFixed(1)}`).join('');
        return (
          <g key={p.label}>
            <rect x={start} y={0} width={w} height={H} fill={p.label === 'repeat' ? '#eef2ff' : '#f8fafc'} />
            {p.sup.length > 0 && <path d={d} fill="none" stroke="#4f46e5" strokeWidth={1} />}
            <text x={start + 3} y={H + 10} fontSize={9.5} fill="#64748b">{p.label} · {p.label === 'repeat' ? `${p.sup.length} units` : `${p.sup.length} bp`}</text>
          </g>
        );
      })}
      <line x1={0} x2={W} y1={H - CONSENSUS_LOW * H} y2={H - CONSENSUS_LOW * H} stroke="#d97706" strokeDasharray="3 3" />
      <text x={W - 2} y={H - CONSENSUS_LOW * H - 2} textAnchor="end" fontSize={9} fill="#b45309">{Math.round(CONSENSUS_LOW * 100)} %: lowercase under it</text>
    </svg>
  );
}
