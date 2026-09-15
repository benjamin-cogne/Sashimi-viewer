/**
 * Standalone Sashimi viewer: a single HTML page. The user adds BAM/CRAM files (with
 * their .bai/.crai), optionally a reference FASTA, and types a gene. Everything is
 * decoded in the browser; only gene lookups (RefSeq models from the UCSC API, Ensembl REST as
 * fallback) and reference sequence (when no FASTA is given) are fetched from the network.
 */
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import ReactDOM from 'react-dom/client';
import SashimiViewer, { type ViewerSettings, type ViewerState } from '../components/SashimiViewer';
import { buildSession, defaultSessionName, matchSession, parseSession, viewerSettingsOf, type SessionFile } from './session';
import { LocalDataSource, type LocalSample } from './localSource';
import type { GenomeBuild } from './ensembl';
import { parseLocus } from '../components/sashimi/geometry';
import { describeLink, parseLink } from './link';
import '../index.css';

/** Unreleased build (branch dev published under /dev/): banner, tab title and a red favicon, so it is never mistaken for the stable page. */
const DEV = import.meta.env.VITE_DEV_MODE === '1'
  ? { branch: import.meta.env.VITE_DEV_BRANCH || 'dev', sha: (import.meta.env.VITE_DEV_SHA || '').slice(0, 7), date: import.meta.env.VITE_DEV_DATE || '', stable: import.meta.env.VITE_STABLE_URL || '../' }
  : null;
if (DEV) {
  document.title = `[DEV] ${document.title}`;
  const icon = document.querySelector<HTMLLinkElement>('link[rel="icon"]');
  if (icon) icon.href = 'data:image/svg+xml,' + encodeURIComponent('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64"><rect width="64" height="64" rx="14" fill="#dc2626"/><path d="M11 44 Q32 -22 53 44" fill="none" stroke="#fff" stroke-width="6" stroke-linecap="round"/><rect x="2" y="44" width="16" height="14" rx="2" fill="#fff"/><rect x="24" y="44" width="16" height="14" rx="2" fill="#fff"/><rect x="46" y="44" width="16" height="14" rx="2" fill="#fff"/></svg>');
}

/** Request carried by the URL (deep link from another tool), read once at start-up. */
const LINK = parseLink(window.location.hash, window.location.search);
const LINK_TEXT = LINK ? `${LINK.mark.chrom}:${LINK.mark.start.toLocaleString('en-US')}${LINK.mark.end > LINK.mark.start ? `-${LINK.mark.end.toLocaleString('en-US')}` : ''}` : '';

interface Pending { file: File }

const ALIGN_EXT = /\.(bam|cram)$/i;
const INDEX_EXT = /\.(bai|crai)$/i;
const FASTA_EXT = /\.(fa|fasta|fna)(\.gz)?$/i;

/** Pair alignment files with their indexes by name (case-insensitive, Windows tools often shout): x.bam + x.bam.bai or x.bai. */
function pairFiles(files: File[]): { samples: { name: string; kind: 'bam' | 'cram'; file: File; index: File }[]; unmatched: string[]; fasta?: { fa: File; fai: File; gzi?: File }; fastaMissing?: string } {
  const byName = new Map(files.map(f => [f.name.toLowerCase(), f]));
  const lookup = (n: string) => byName.get(n.toLowerCase());
  const samples: { name: string; kind: 'bam' | 'cram'; file: File; index: File }[] = [];
  const unmatched: string[] = [];
  for (const f of files) {
    if (!ALIGN_EXT.test(f.name)) continue;
    const kind = f.name.toLowerCase().endsWith('.cram') ? 'cram' : 'bam';
    const stem = f.name.replace(ALIGN_EXT, '');
    const idx = lookup(`${f.name}.${kind === 'bam' ? 'bai' : 'crai'}`) || lookup(`${stem}.${kind === 'bam' ? 'bai' : 'crai'}`);
    if (idx) samples.push({ name: stem, kind, file: f, index: idx });
    else unmatched.push(f.name);
  }
  let fasta: { fa: File; fai: File; gzi?: File } | undefined, fastaMissing: string | undefined;
  const fa = files.find(f => FASTA_EXT.test(f.name));
  if (fa) {
    const fai = lookup(`${fa.name}.fai`);
    const gzi = lookup(`${fa.name}.gzi`);
    if (fai && (!fa.name.toLowerCase().endsWith('.gz') || gzi)) fasta = { fa, fai, gzi };
    else fastaMissing = fa.name;
  }
  return { samples, unmatched, fasta, fastaMissing };
}

/** The project mark (docs/logo): three exons, canonical arcs, and a dashed red arc skipping the middle exon. */
function Logo({ size = 28 }: { size?: number }) {
  return (
    <svg width={size} height={size} viewBox="0 0 64 64" aria-hidden="true" className="shrink-0">
      <path d="M11 46 Q21.5 14 32 46 M32 46 Q42.5 14 53 46" fill="none" stroke="#3730a3" strokeWidth="3.5" strokeLinecap="round" />
      <path d="M11 46 Q32 -26 53 46" fill="none" stroke="#dc2626" strokeWidth="3.5" strokeLinecap="round" strokeDasharray="5 4" />
      <path d="M18 51 H25 M39 51 H46" stroke="#3730a3" strokeWidth="2.5" strokeLinecap="round" />
      <rect x="4" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
      <rect x="25" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
      <rect x="46" y="46" width="14" height="10" rx="1.5" fill="#3730a3" />
    </svg>
  );
}

function App() {
  const [build, setBuild] = useState<GenomeBuild>(LINK?.build ?? 'GRCh38');
  const [samples, setSamples] = useState<LocalSample[]>([]);
  const [renaming, setRenaming] = useState<{ id: number; value: string } | null>(null);
  const [fasta, setFasta] = useState<{ fa: File; fai: File; gzi?: File } | undefined>();
  const [notes, setNotes] = useState<string[]>([]);
  const [gene, setGene] = useState(LINK_TEXT);
  const [opened, setOpened] = useState<{ geneName: string; chrom: string; start: number; end: number; view?: { start: number; end: number }; mark?: { start: number; end: number }; reads?: boolean } | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // ---- Sessions: the viewer's latest state (options + navigation), a loaded session waiting for its files ----
  const viewerStateRef = useRef<ViewerState | null>(null);
  const [sessionName, setSessionName] = useState(() => defaultSessionName());
  const [sessionNameEdited, setSessionNameEdited] = useState(false);
  const [pendingSession, setPendingSession] = useState<{ file: SessionFile; name: string } | null>(null);
  const [viewerSettings, setViewerSettings] = useState<Partial<ViewerSettings> | undefined>();
  const [sessionSeq, setSessionSeq] = useState(0);   // bumps the viewer key so a loaded session always remounts the viewer
  const nextId = useRef(1);
  const dsRef = useRef<LocalDataSource>();
  if (!dsRef.current) { dsRef.current = new LocalDataSource({ build }); if (LINK) dsRef.current.knownVariants = LINK.variants; }
  const ds = dsRef.current;
  (window as any).__sashimiDs = ds; // for debugging from the console

  const addFiles = useCallback((list: FileList | File[]) => {
    const files = Array.from(list);
    const { samples: found, unmatched, fasta: fa, fastaMissing } = pairFiles(files);
    const msgs: string[] = [];
    const added: LocalSample[] = found.map(s => ({ id: nextId.current++, ...s }));
    for (const s of added) ds.addSample(s);
    if (added.length) setSamples(prev => [...prev, ...added]);
    if (unmatched.length) msgs.push(`No index found for ${unmatched.join(', ')} (add the .bai / .crai file together with it)`);
    if (fa) { setFasta(fa); ds.setReference({ build, fasta: fa }); msgs.push(`Reference FASTA: ${fa.fa.name}`); }
    if (fastaMissing) msgs.push(`${fastaMissing} needs its .fai index${fastaMissing.toLowerCase().endsWith('.gz') ? ' and .gzi' : ''}`);
    setNotes(msgs);
  }, [ds, build]);

  const renameSample = useCallback((id: number, raw: string) => {
    const name = raw.trim();
    setRenaming(null);
    if (!name) return;
    ds.renameSample(id, name);
    setSamples(prev => prev.map(s => (s.id === id ? { ...s, name } : s)));
  }, [ds]);
  const sampleNames = useMemo(() => Object.fromEntries(samples.map(s => [s.id, s.name])) as Record<number, string>, [samples]);

  const removeSample = useCallback((id: number) => { ds.removeSample(id); setSamples(prev => prev.filter(s => s.id !== id)); }, [ds]);

  const changeBuild = (b: GenomeBuild) => { setBuild(b); ds.setReference({ build: b, fasta }); };

  /** Open the window a deep link asked for: the gene at the locus (coding first), else the `gene` parameter, with the variant pinned. */
  const openLink = useCallback(async () => {
    if (!LINK) return;
    setBusy(true); setError(null);
    try {
      const { view, mark } = LINK;
      const genes = await ds.getRegionGenes(view.chrom, view.start, view.end);
      const ov = (g: { start: number; end: number }) => Math.min(g.end, view.end) - Math.max(g.start, view.start);
      let best = [...genes].sort((a, b) => Number(b.biotype === 'protein_coding') - Number(a.biotype === 'protein_coding') || ov(b) - ov(a))[0];
      if (!best && LINK.gene) {
        const t = await ds.getTranscript(LINK.gene, LINK.gene.toUpperCase().startsWith('ENSG') ? LINK.gene : undefined);
        best = { gene_name: t.gene_name, start: t.start, end: t.end } as any;
      }
      if (!best) throw new Error(`no RefSeq gene at ${view.chrom}:${view.start.toLocaleString()}-${view.end.toLocaleString()} (add gene=SYMBOL to the link)`);
      setOpened({ geneName: best.gene_name, chrom: view.chrom, start: best.start, end: best.end, view: { start: view.start, end: view.end }, mark: { start: mark.start, end: mark.end }, reads: LINK.reads });
    } catch (e: any) {
      setError(`Could not open the linked locus: ${e.message}. Check the genome build (${LINK.build}) and that api.genome.ucsc.edu (or rest.ensembl.org) is reachable.`);
    }
    setBusy(false);
  }, [ds]);

  // A link opens the browser right away; files added afterwards become tracks
  const autoOpened = useRef(false);
  useEffect(() => {
    if (LINK && !autoOpened.current) { autoOpened.current = true; openLink(); }
  }, [openLink]);

  const open = useCallback(async () => {
    const q = gene.trim();
    if (!q) return;
    if (LINK && q === LINK_TEXT) return openLink();
    setBusy(true); setError(null);
    try {
      const locus = parseLocus(q);
      if (locus) {
        // coordinates: open the gene at the locus (coding first, largest overlap) and show the locus itself
        const genes = await ds.getRegionGenes(locus.chrom, locus.start, locus.end);
        const ov = (g: { start: number; end: number }) => Math.min(g.end, locus.end) - Math.max(g.start, locus.start);
        const best = [...genes].sort((a, b) => Number(b.biotype === 'protein_coding') - Number(a.biotype === 'protein_coding') || ov(b) - ov(a))[0];
        if (!best) throw new Error(`no RefSeq gene at ${q}`);
        const point = locus.start === locus.end;
        setOpened({ geneName: best.gene_name, chrom: locus.chrom, start: best.start, end: best.end, view: { start: point ? Math.max(1, locus.start - 500) : locus.start, end: point ? locus.start + 500 : locus.end } });
      } else {
        const t = await ds.getTranscript(q, q.toUpperCase().startsWith('ENSG') ? q : undefined);
        setOpened({ geneName: t.gene_name, chrom: t.chrom, start: t.start, end: t.end });
      }
    } catch (e: any) {
      setError(`Gene lookup failed: ${e.message}. Check the symbol, the genome build and that api.genome.ucsc.edu (or rest.ensembl.org) is reachable.`);
    }
    setBusy(false);
  }, [gene, ds, openLink]);

  // Options the viewer starts with: a loaded session, else the previous viewer's options (kept across gene searches;
  // a chosen reference transcript only when it is the same gene)
  const viewerInit: Partial<ViewerSettings> | undefined = (() => {
    if (viewerSettings) return viewerSettings;
    const prev = viewerStateRef.current;
    if (!prev) return undefined;
    return { ...prev, transcriptId: prev.gene.name === opened?.geneName ? prev.transcriptId : undefined };
  })();
  // order-independent: promoting another sample to primary keeps the viewer (and its view) mounted
  const viewerKey = useMemo(() => `${opened?.geneName}|${opened?.view ? `${opened.view.start}-${opened.view.end}` : ''}|${opened?.mark ? `${opened.mark.start}-${opened.mark.end}` : ''}|${[...samples.map(s => s.id)].sort((a, b) => a - b).join(',')}|${build}|${fasta?.fa.name || ''}|${sessionSeq}`, [opened, samples, build, fasta, sessionSeq]);
  const makePrimary = useCallback((id: number) => setSamples(prev => [...prev.filter(s => s.id === id), ...prev.filter(s => s.id !== id)]), []);

  // ---- Save / load a session (JSON) ----
  const saveSession = useCallback(() => {
    const session = buildSession({ build, samples, fasta, state: opened ? viewerStateRef.current : null });
    const name = (sessionName.trim() || defaultSessionName(session.gene?.name)).replace(/\.json$/i, '') + '.json';
    const url = URL.createObjectURL(new Blob([JSON.stringify(session, null, 2)], { type: 'application/json' }));
    const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    setNotes([`Session saved as ${name} (${session.samples.length} sample${session.samples.length === 1 ? '' : 's'}${session.gene ? `, ${session.gene.name}` : ''}). Files are recorded by name: add them again when loading.`]);
  }, [build, samples, fasta, opened, sessionName]);

  /** Applies a loaded session with the files present: names, order, options, then the gene and window. */
  const applySession = useCallback((session: SessionFile) => {
    const { matched } = matchSession(session, samples);
    for (const { entry, sample } of matched) if (sample.name !== entry.name) ds.renameSample(sample.id, entry.name);
    const renamed = samples.map(s => { const m = matched.find(x => x.sample.id === s.id); return m ? { ...s, name: m.entry.name } : s; });
    const ordered = [...matched.map(m => renamed.find(s => s.id === m.sample.id)!), ...renamed.filter(s => !matched.some(m => m.sample.id === s.id))];
    setSamples(ordered);
    setViewerSettings(viewerSettingsOf(session, ordered));
    setSessionSeq(n => n + 1);
    setPendingSession(null);
    if (session.gene) {
      setGene(session.gene.name);
      setOpened({ geneName: session.gene.name, chrom: session.gene.chrom, start: session.gene.start, end: session.gene.end, view: session.gene.view, mark: session.gene.mark ?? undefined });
    }
    setNotes([`Session loaded: ${matched.length}/${session.samples.length} sample${session.samples.length === 1 ? '' : 's'}${session.gene ? `, ${session.gene.name}` : ''}.`]);
  }, [samples, ds]);

  const loadSession = useCallback(async (file: File) => {
    setError(null);
    try {
      const session = parseSession(await file.text());
      setSessionName(file.name); setSessionNameEdited(true);
      if (session.build !== build) { setBuild(session.build); ds.setReference({ build: session.build, fasta }); }
      setPendingSession({ file: session, name: file.name });
    } catch (e: any) {
      setError(`Could not load the session ${file.name}: ${e.message}`);
    }
  }, [build, ds, fasta]);

  // A pending session applies itself as soon as every one of its files is present (or right away when it needs none)
  useEffect(() => {
    if (!pendingSession) return;
    const { missing } = matchSession(pendingSession.file, samples);
    if (!missing.length) applySession(pendingSession.file);
  }, [pendingSession, samples, applySession]);

  // Default file name follows the gene until the user types one
  useEffect(() => { if (!sessionNameEdited) setSessionName(defaultSessionName(opened?.geneName)); }, [opened?.geneName, sessionNameEdited]);

  const dropRef = useRef<HTMLDivElement>(null);
  const onDrop = (e: React.DragEvent) => { e.preventDefault(); addFiles(e.dataTransfer.files); };

  return (
    <div className="min-h-screen bg-gray-100 text-gray-900">
      {DEV && (
        <div className="px-5 py-1.5 text-xs font-semibold text-white flex flex-wrap items-center gap-x-3 gap-y-1"
          style={{ background: 'repeating-linear-gradient(135deg, #b91c1c 0 14px, #dc2626 14px 28px)' }}>
          <span className="px-1.5 py-0.5 rounded bg-white text-red-700 tracking-wider">DEV MODE</span>
          <span>Unreleased build from branch <code>{DEV.branch}</code>{DEV.sha ? <> · commit <code>{DEV.sha}</code></> : null}{DEV.date ? ` · built ${DEV.date}` : ''} · not for clinical use.</span>
          <a href={DEV.stable} className="underline ml-auto">Go to the stable version →</a>
        </div>
      )}
      <header className="bg-white border-b border-gray-200 px-5 py-3 space-y-2">
        {/* Row 1: title, build, gene search. Row 2: files and sample chips, which may wrap over several lines without moving the search. */}
        <div className="flex flex-wrap items-center gap-x-6 gap-y-2">
        <div className="flex items-center gap-2.5 flex-1 min-w-[320px]">
          <Logo size={34} />
          <div className="min-w-0">
            <h1 className="text-lg font-bold leading-tight">Sashimi <span className="font-normal">viewer</span></h1>
            <p className="text-xs text-gray-500">Files are read in your browser and never uploaded. Gene models (RefSeq, UCSC API) and reference bases come from the network unless you add a FASTA.</p>
            <p className="text-xs font-medium text-amber-700" role="note">⚠ Check the HGVS nomenclature, the predicted transcript, amino-acid changes and the NMD verdict before anything from it goes into a clinical report.</p>
          </div>
        </div>
        <div className="flex items-center gap-4 shrink-0 ml-auto">
        <label className="flex items-center gap-1 text-xs text-gray-600">Build
          <select value={build} onChange={e => changeBuild(e.target.value as GenomeBuild)} className="border border-gray-300 rounded px-1 py-0.5 text-xs bg-white">
            <option value="GRCh38">GRCh38 / hg38</option>
            <option value="GRCh37">GRCh37 / hg19</option>
          </select>
        </label>
        <form onSubmit={e => { e.preventDefault(); open(); }} className="flex items-center gap-1">
          <input value={gene} onChange={e => setGene(e.target.value)} placeholder="Gene, ENSG or chr:pos…" title="A gene symbol, an ENSG id, or coordinates (chr17:43,094,464 or chr17:43,000,000-43,100,000: the gene at the locus is opened)" className="border border-gray-300 rounded px-2 py-1 text-sm w-48 bg-white" />
          <button type="submit" disabled={busy || !gene.trim()} className="px-3 py-1 text-sm rounded bg-indigo-600 text-white disabled:opacity-40 hover:bg-indigo-700">{busy ? '…' : 'Open'}</button>
        </form>
        </div>
        </div>
        <div className="flex flex-wrap items-center gap-x-3 gap-y-1.5">
        <label className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 cursor-pointer font-medium">
          + Add BAM / CRAM (+ index) · FASTA
          <input type="file" multiple className="hidden" accept=".bam,.bai,.cram,.crai,.fa,.fasta,.fna,.gz,.fai,.gzi" onChange={e => { if (e.target.files) addFiles(e.target.files); e.target.value = ''; }} />
        </label>
        <div className="flex flex-wrap items-center gap-1.5">
          {samples.map((s, i) => (
            <span key={s.id} className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border ${i === 0 ? 'bg-indigo-50 border-indigo-300 text-indigo-800' : 'bg-gray-50 border-gray-300 text-gray-700 hover:border-indigo-300 cursor-pointer'}`}
              title={`${s.file.name} · ${(s.file.size / 1e9).toFixed(2)} GB · ${s.kind.toUpperCase()}${i === 0 ? ' · primary sample' : ' · click to make it the primary sample'} · double-click to rename`}
              onClick={() => { if (i !== 0 && renaming?.id !== s.id) makePrimary(s.id); }}
              onDoubleClick={e => { e.stopPropagation(); setRenaming({ id: s.id, value: s.name }); }}>
              {i === 0 && <span title="primary sample">★</span>}
              {renaming?.id === s.id ? (
                <input autoFocus value={renaming.value} onChange={e => setRenaming({ id: s.id, value: e.target.value })}
                  onBlur={() => renameSample(s.id, renaming.value)}
                  onKeyDown={e => { if (e.key === 'Enter') renameSample(s.id, renaming.value); else if (e.key === 'Escape') setRenaming(null); }}
                  onClick={e => e.stopPropagation()} className="w-32 px-1 py-0 text-xs rounded border border-indigo-300 bg-white text-gray-900" title="Enter to confirm, Esc to cancel" />
              ) : s.name}
              <button onClick={e => { e.stopPropagation(); setRenaming({ id: s.id, value: s.name }); }} className="text-gray-400 hover:text-indigo-600" title="Rename">✎</button>
              <button onClick={e => { e.stopPropagation(); removeSample(s.id); }} className="text-gray-400 hover:text-red-500" title="Remove">×</button>
            </span>
          ))}
          {fasta && <span className="px-2 py-0.5 rounded-full text-xs border bg-emerald-50 border-emerald-300 text-emerald-800" title={fasta.fa.name}>FASTA · {fasta.fa.name}</span>}
        </div>
        <div className="flex items-center gap-1.5 ml-auto" title="A session file (JSON) records the alignment files by name, the sample names and order, the FASTA, the gene and window, and every option of the viewer. Load it later and add the same files again.">
          <span className="text-xs text-gray-600">Session</span>
          <input value={sessionName} onChange={e => { setSessionName(e.target.value); setSessionNameEdited(true); }} spellCheck={false}
            className="border border-gray-300 rounded px-2 py-0.5 text-xs w-64 bg-white font-mono" title="File name of the session to save (.json)" />
          <button onClick={saveSession} disabled={!samples.length && !opened} className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 font-medium disabled:opacity-40" title="Download the session as a JSON file">Save session</button>
          <label className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 cursor-pointer font-medium" title="Load a session JSON file, then add the alignment files it names">
            Load session
            <input type="file" accept=".json,application/json" className="hidden" onChange={e => { const f = e.target.files?.[0]; if (f) loadSession(f); e.target.value = ''; }} />
          </label>
        </div>
        </div>
      </header>
      {pendingSession && (() => {
        const { matched, missing } = matchSession(pendingSession.file, samples);
        return (
          <div className="mx-5 mt-2 px-3 py-2 text-xs rounded-lg border border-indigo-200 bg-indigo-50 text-indigo-900 flex flex-wrap items-center gap-x-3 gap-y-1">
            <span className="font-semibold">Session {pendingSession.name}</span>
            <span>{matched.length}/{pendingSession.file.samples.length} files present · add (or drop on the page): {missing.map(m => m.file).join(', ')}</span>
            <button onClick={() => applySession(pendingSession.file)} className="px-2 py-0.5 rounded border border-indigo-300 bg-white hover:bg-indigo-100 font-medium">Open with the files present</button>
            <button onClick={() => setPendingSession(null)} className="text-indigo-500 hover:text-indigo-800" title="Forget this session">×</button>
          </div>
        );
      })()}
      {(notes.length > 0 || error || (opened && !samples.length)) && (
        <div className="px-5 py-2 text-xs space-y-0.5">
          {opened && !samples.length && (
            <div className="text-indigo-800">
              {LINK ? `${describeLink(LINK)} · ` : ''}No alignment yet: add the BAM or CRAM files (with their index) with the button above, or drop them on the page. They stay on this computer.
            </div>
          )}
          {notes.map((n, i) => <div key={i} className="text-gray-600">{n}</div>)}
          {error && <div className="text-red-600">{error}</div>}
        </div>
      )}
      {!opened ? (
        <div ref={dropRef} onDragOver={e => e.preventDefault()} onDrop={onDrop}
          className="m-6 p-10 border-2 border-dashed border-indigo-300 rounded-2xl bg-white text-center">
          <div className="text-xl font-semibold text-indigo-700">Drop BAM or CRAM files here</div>
          <div className="text-sm text-gray-600 mt-2 max-w-2xl mx-auto">
            Add each alignment with its index (<code>.bam</code> + <code>.bai</code>, or <code>.cram</code> + <code>.crai</code>). The first file is the primary sample, the others are comparison samples; click a sample chip (or "make primary" on its track) to switch.
            CRAM needs the reference: add an indexed FASTA (<code>.fa</code> + <code>.fai</code>, bgzipped with <code>.gzi</code>) or let the page fetch it from the UCSC API.
            Then type a gene and press Open.
          </div>
          <div className="text-xs text-gray-500 mt-4">
            {samples.length ? `${samples.length} sample${samples.length > 1 ? 's' : ''} ready` : 'No files yet'} · reads are decoded locally with the GMOD BAM/CRAM libraries
          </div>
        </div>
      ) : (
        <div className="p-3" onDragOver={e => e.preventDefault()} onDrop={onDrop}>
          <SashimiViewer key={viewerKey} geneName={opened.geneName} chrom={opened.chrom} geneStart={opened.start} geneEnd={opened.end}
            sampleId={samples[0]?.id ?? 0} sampleName={samples[0]?.name ?? ''} runId={0} darkMode={false} onClose={() => setOpened(null)} embedded dataSource={ds} allowPrimarySwitch onPrimaryChange={makePrimary} initialView={opened.view} initialMark={opened.mark} initialReads={opened.reads} sampleNames={sampleNames}
            initialSettings={viewerInit} onStateChange={s => { viewerStateRef.current = s; if (viewerSettings) setViewerSettings(undefined); }} />
        </div>
      )}
      <footer className="px-5 py-3 text-[11px] text-gray-500 flex flex-wrap gap-x-3 gap-y-1">
        <span>Benjamin Cogné (CHU Nantes, 2026)</span>
        <span>·</span>
        <span>made with Claude Opus 5 and Fable 5.1</span>
        <span>·</span>
        <a href="https://github.com/benjamin-cogne/Sashimi-viewer" target="_blank" rel="noopener noreferrer" className="hover:text-indigo-700 underline decoration-dotted">github.com/benjamin-cogne/Sashimi-viewer</a>
        <span>·</span>
        <span>CC BY-NC 4.0</span>
        {DEV && <span className="text-red-700 font-semibold">DEV build {DEV.branch}{DEV.sha ? ` @ ${DEV.sha}` : ''}</span>}
      </footer>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById('root')!).render(<App />);
