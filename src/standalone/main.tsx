/**
 * Standalone Sashimi viewer: a single HTML page. The user adds BAM/CRAM files, optionally a
 * reference FASTA, and types a gene. Indexes are matched by name (sample.bam.bai or sample.bai,
 * sample.cram.crai or sample.crai) in whatever order the files arrive; giving the page the folder
 * (picker or drop) lets it find them itself. Everything is decoded in the browser; only gene
 * lookups (RefSeq models from the UCSC API, Ensembl REST as fallback) and reference sequence (when
 * no FASTA is given) are fetched from the network.
 */
import React, { useCallback, useMemo, useRef, useState } from 'react';
import ReactDOM from 'react-dom/client';
import SashimiViewer from '../components/SashimiViewer';
import { LocalDataSource, type LocalSample } from './localSource';
import type { GenomeBuild } from './ensembl';
import { parseLocus } from '../components/sashimi/geometry';
import { filesFromDirectoryHandle, filesFromDrop, pairPool, supportsDirectoryPicker, WANTED_EXT, type PairedFasta, type PendingFile } from './files';
import '../index.css';

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

/** A file input styled as a small button; the value is cleared so the same file can be chosen again. */
function Picker({ label, accept, directory, className, inputRef, onFiles }: { label: React.ReactNode; accept?: string; directory?: boolean; className?: string; inputRef?: React.RefObject<HTMLInputElement>; onFiles: (files: FileList) => void }) {
  return (
    <label className={className}>
      {label}
      <input ref={inputRef} type="file" multiple className="hidden" accept={accept} {...(directory ? ({ webkitdirectory: '', directory: '' } as any) : {})}
        onChange={e => { if (e.target.files) onFiles(e.target.files); e.target.value = ''; }} />
    </label>
  );
}

const ACCEPT_ANY = '.bam,.bai,.csi,.cram,.crai,.fa,.fasta,.fna,.gz,.fai,.gzi';
const ACCEPT_INDEX = '.bai,.csi,.crai,.fai,.gzi';

function App() {
  const [build, setBuild] = useState<GenomeBuild>('GRCh38');
  const [samples, setSamples] = useState<LocalSample[]>([]);
  const [pending, setPending] = useState<PendingFile[]>([]);
  const [fasta, setFasta] = useState<PairedFasta | undefined>();
  const [fastaPending, setFastaPending] = useState<PendingFile | undefined>();
  const [notes, setNotes] = useState<string[]>([]);
  const [gene, setGene] = useState('');
  const [opened, setOpened] = useState<{ geneName: string; chrom: string; start: number; end: number; view?: { start: number; end: number } } | null>(null);
  const [busy, setBusy] = useState(false);
  const [dragging, setDragging] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const nextId = useRef(1);
  /** Every file added so far, by lower-cased name (Windows tools often shout): later files replace earlier ones. */
  const pool = useRef(new Map<string, File>());
  /** Lower-cased names of alignment files already registered as samples. */
  const registered = useRef(new Map<string, number>());
  const dsRef = useRef<LocalDataSource>();
  if (!dsRef.current) dsRef.current = new LocalDataSource({ build });
  const ds = dsRef.current;
  (window as any).__sashimiDs = ds; // for debugging from the console

  /** Put files in the pool and re-pair everything: new samples get an id, alignments without index wait, the FASTA is picked up. */
  const addFiles = useCallback((list: FileList | File[]) => {
    const files = Array.from(list).filter(f => WANTED_EXT.test(f.name));
    const skipped = Array.from(list).length - files.length;
    for (const f of files) pool.current.set(f.name.toLowerCase(), f);
    const r = pairPool(pool.current);
    const added: LocalSample[] = [];
    for (const s of r.samples) {
      const key = s.file.name.toLowerCase();
      if (registered.current.has(key)) continue;
      const sample: LocalSample = { id: nextId.current++, ...s };
      registered.current.set(key, sample.id);
      ds.addSample(sample);
      added.push(sample);
    }
    if (added.length) setSamples(prev => [...prev, ...added]);
    setPending(r.pending);
    setFastaPending(r.fastaPending);
    const msgs: string[] = [];
    if (r.fasta && r.fasta.fa !== fasta?.fa) { setFasta(r.fasta); ds.setReference({ build, fasta: r.fasta }); msgs.push(`Reference FASTA: ${r.fasta.fa.name}`); }
    if (skipped) msgs.push(`${skipped} file${skipped > 1 ? 's' : ''} ignored (not BAM/CRAM, index or FASTA)`);
    if (!files.length && !skipped) msgs.push('No files received');
    setNotes(msgs);
  }, [ds, build, fasta]);

  const removeSample = useCallback((id: number) => {
    ds.removeSample(id);
    setSamples(prev => prev.filter(s => s.id !== id));
    for (const [key, sid] of registered.current) if (sid === id) { registered.current.delete(key); pool.current.delete(key); }
  }, [ds]);

  const forgetPending = useCallback((name: string) => {
    pool.current.delete(name.toLowerCase());
    setPending(prev => prev.filter(p => p.file.name !== name));
  }, []);

  const changeBuild = (b: GenomeBuild) => { setBuild(b); ds.setReference({ build: b, fasta }); };

  /** Folder: File System Access API where available (Edge, Chrome), otherwise the directory input. Only wanted extensions are read. */
  const folderInput = useRef<HTMLInputElement>(null);
  const openFolder = useCallback(async () => {
    if (!supportsDirectoryPicker()) { folderInput.current?.click(); return; }
    try {
      const dir = await (window as any).showDirectoryPicker({ mode: 'read' });
      addFiles(await filesFromDirectoryHandle(dir));
    } catch (e: any) {
      if (e?.name !== 'AbortError') setError(`Could not read the folder: ${e?.message || e}`);
    }
  }, [addFiles]);

  const open = useCallback(async () => {
    const q = gene.trim();
    if (!q || !samples.length) return;
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
  }, [gene, samples.length, ds]);

  // order-independent: promoting another sample to primary keeps the viewer (and its view) mounted
  const viewerKey = useMemo(() => `${opened?.geneName}|${opened?.view ? `${opened.view.start}-${opened.view.end}` : ''}|${[...samples.map(s => s.id)].sort((a, b) => a - b).join(',')}|${build}|${fasta?.fa.name || ''}`, [opened, samples, build, fasta]);
  const makePrimary = useCallback((id: number) => setSamples(prev => [...prev.filter(s => s.id === id), ...prev.filter(s => s.id !== id)]), []);

  const onDrop = async (e: React.DragEvent) => {
    e.preventDefault(); setDragging(false);
    try { addFiles(await filesFromDrop(e.dataTransfer)); }
    catch (err: any) { setError(`Could not read the dropped items: ${err?.message || err}`); }
  };
  const dropProps = {
    onDragOver: (e: React.DragEvent) => { e.preventDefault(); if (!dragging) setDragging(true); },
    onDragLeave: () => setDragging(false),
    onDrop,
  };

  const btn = 'px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 cursor-pointer font-medium';

  return (
    <div className="min-h-screen bg-gray-100 text-gray-900" {...dropProps}>
      <header className="bg-white border-b border-gray-200 px-5 py-3 flex flex-wrap items-center gap-x-6 gap-y-2">
        <div className="flex items-center gap-2.5">
          <Logo size={34} />
          <div>
            <h1 className="text-lg font-bold leading-tight">Sashimi <span className="font-normal">viewer</span></h1>
            <p className="text-xs text-gray-500">Files are read in your browser and never uploaded. Gene models (RefSeq, UCSC API) and reference bases come from the network unless you add a FASTA.</p>
          </div>
        </div>
        <label className="flex items-center gap-1 text-xs text-gray-600">Build
          <select value={build} onChange={e => changeBuild(e.target.value as GenomeBuild)} className="border border-gray-300 rounded px-1 py-0.5 text-xs bg-white">
            <option value="GRCh38">GRCh38 / hg38</option>
            <option value="GRCh37">GRCh37 / hg19</option>
          </select>
        </label>
        <div className="flex items-center gap-1.5">
          <Picker label="+ Add BAM / CRAM" accept={ACCEPT_ANY} className={btn} onFiles={addFiles} />
          <button type="button" onClick={openFolder} className={btn} title="Read a folder: every BAM/CRAM in it is paired with its index by name; a FASTA (+ .fai) is used as reference. Nothing is uploaded.">Open folder…</button>
          <Picker label="" directory inputRef={folderInput} className="hidden" onFiles={addFiles} />
        </div>
        <div className="flex flex-wrap items-center gap-1.5">
          {samples.map((s, i) => (
            <span key={s.id} className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border ${i === 0 ? 'bg-indigo-50 border-indigo-300 text-indigo-800' : 'bg-gray-50 border-gray-300 text-gray-700 hover:border-indigo-300 cursor-pointer'}`}
              title={`${s.file.name} + ${s.index.name} · ${(s.file.size / 1e9).toFixed(2)} GB · ${s.kind.toUpperCase()}${i === 0 ? ' · primary sample' : ' · click to make it the primary sample'}`}
              onClick={() => { if (i !== 0) makePrimary(s.id); }}>
              {i === 0 && <span title="primary sample">★</span>}{s.name}<button onClick={e => { e.stopPropagation(); removeSample(s.id); }} className="text-gray-400 hover:text-red-500" title="Remove">×</button>
            </span>
          ))}
          {pending.map(p => (
            <span key={p.file.name} className="flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border bg-amber-50 border-amber-300 text-amber-900" title={`${p.file.name} has no index yet. Expected next to it: ${p.wanted.join(' or ')}`}>
              {p.file.name}
              <Picker label={<span className="underline decoration-dotted cursor-pointer">add {p.wanted[0]}</span>} accept={ACCEPT_INDEX} className="contents" onFiles={addFiles} />
              <button onClick={() => forgetPending(p.file.name)} className="text-amber-500 hover:text-red-500" title="Forget this file">×</button>
            </span>
          ))}
          {fasta && <span className="px-2 py-0.5 rounded-full text-xs border bg-emerald-50 border-emerald-300 text-emerald-800" title={`${fasta.fa.name} + ${fasta.fai.name}${fasta.gzi ? ' + ' + fasta.gzi.name : ''}`}>FASTA · {fasta.fa.name}</span>}
          {fastaPending && (
            <span className="flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border bg-amber-50 border-amber-300 text-amber-900" title={`${fastaPending.file.name} needs ${fastaPending.wanted.join(' and ')}`}>
              {fastaPending.file.name}
              <Picker label={<span className="underline decoration-dotted cursor-pointer">add {fastaPending.wanted.join(' + ')}</span>} accept={ACCEPT_INDEX} className="contents" onFiles={addFiles} />
            </span>
          )}
        </div>
        <form onSubmit={e => { e.preventDefault(); open(); }} className="flex items-center gap-1 ml-auto">
          <input value={gene} onChange={e => setGene(e.target.value)} placeholder="Gene, ENSG or chr:pos…" title="A gene symbol, an ENSG id, or coordinates (chr17:43,094,464 or chr17:43,000,000-43,100,000: the gene at the locus is opened)" className="border border-gray-300 rounded px-2 py-1 text-sm w-48 bg-white" />
          <button type="submit" disabled={busy || !samples.length || !gene.trim()} className="px-3 py-1 text-sm rounded bg-indigo-600 text-white disabled:opacity-40 hover:bg-indigo-700">{busy ? '…' : 'Open'}</button>
        </form>
      </header>
      {(notes.length > 0 || error) && (
        <div className="px-5 py-2 text-xs space-y-0.5">
          {notes.map((n, i) => <div key={i} className="text-gray-600">{n}</div>)}
          {error && <div className="text-red-600">{error}</div>}
        </div>
      )}
      {!opened ? (
        <div className={`m-6 p-10 border-2 border-dashed rounded-2xl bg-white text-center transition-colors ${dragging ? 'border-indigo-500 bg-indigo-50' : 'border-indigo-300'}`}>
          <div className="text-xl font-semibold text-indigo-700">Drop BAM or CRAM files, or their folder, here</div>
          <div className="text-sm text-gray-600 mt-2 max-w-2xl mx-auto space-y-1">
            <p>Each alignment is paired with its index by name (<code>sample.bam</code> + <code>sample.bam.bai</code> or <code>sample.bai</code>; <code>sample.cram</code> + <code>sample.cram.crai</code> or <code>sample.crai</code>), in any order. Drop the folder, or use <em>Open folder…</em>, and the page finds the indexes itself. A browser cannot read a file that was not given to it, so a BAM added alone waits until its index arrives.</p>
            <p>The first file is the primary sample, the others are comparison samples; click a sample chip (or "make primary" on its track) to switch. CRAM needs the reference: add an indexed FASTA (<code>.fa</code> + <code>.fai</code>, bgzipped with <code>.gzi</code>) or let the page fetch it from the UCSC API. Then type a gene and press Open.</p>
          </div>
          <div className="text-xs text-gray-500 mt-4">
            {samples.length ? `${samples.length} sample${samples.length > 1 ? 's' : ''} ready` : 'No files yet'}{pending.length ? ` · ${pending.length} waiting for an index` : ''} · reads are decoded locally with the GMOD BAM/CRAM libraries
          </div>
        </div>
      ) : (
        <div className="p-3">
          <SashimiViewer key={viewerKey} geneName={opened.geneName} chrom={opened.chrom} geneStart={opened.start} geneEnd={opened.end}
            sampleId={samples[0].id} sampleName={samples[0].name} runId={0} darkMode={false} onClose={() => setOpened(null)} embedded dataSource={ds} allowPrimarySwitch onPrimaryChange={makePrimary} initialView={opened.view} />
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
      </footer>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById('root')!).render(<App />);
