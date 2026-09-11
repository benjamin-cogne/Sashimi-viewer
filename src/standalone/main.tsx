/**
 * Standalone Sashimi viewer: a single HTML page. The user adds BAM/CRAM files (with
 * their .bai/.crai), optionally a reference FASTA, and types a gene. Everything is
 * decoded in the browser; only gene lookups (RefSeq models from the UCSC API, Ensembl REST as
 * fallback) and reference sequence (when no FASTA is given) are fetched from the network.
 */
import React, { useCallback, useMemo, useRef, useState } from 'react';
import ReactDOM from 'react-dom/client';
import SashimiViewer from '../components/SashimiViewer';
import { LocalDataSource, type LocalSample } from './localSource';
import type { GenomeBuild } from './ensembl';
import { parseLocus } from '../components/sashimi/geometry';
import '../index.css';

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
  const [build, setBuild] = useState<GenomeBuild>('GRCh38');
  const [samples, setSamples] = useState<LocalSample[]>([]);
  const [fasta, setFasta] = useState<{ fa: File; fai: File; gzi?: File } | undefined>();
  const [notes, setNotes] = useState<string[]>([]);
  const [gene, setGene] = useState('');
  const [opened, setOpened] = useState<{ geneName: string; chrom: string; start: number; end: number; view?: { start: number; end: number } } | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const nextId = useRef(1);
  const dsRef = useRef<LocalDataSource>();
  if (!dsRef.current) dsRef.current = new LocalDataSource({ build });
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

  const removeSample = useCallback((id: number) => { ds.removeSample(id); setSamples(prev => prev.filter(s => s.id !== id)); }, [ds]);

  const changeBuild = (b: GenomeBuild) => { setBuild(b); ds.setReference({ build: b, fasta }); };

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

  const dropRef = useRef<HTMLDivElement>(null);
  const onDrop = (e: React.DragEvent) => { e.preventDefault(); addFiles(e.dataTransfer.files); };

  return (
    <div className="min-h-screen bg-gray-100 text-gray-900">
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
        <label className="px-3 py-1 text-xs rounded border border-gray-300 bg-white hover:bg-indigo-50 cursor-pointer font-medium">
          + Add BAM / CRAM (+ index) · FASTA
          <input type="file" multiple className="hidden" accept=".bam,.bai,.cram,.crai,.fa,.fasta,.fna,.gz,.fai,.gzi" onChange={e => { if (e.target.files) addFiles(e.target.files); e.target.value = ''; }} />
        </label>
        <div className="flex flex-wrap items-center gap-1.5">
          {samples.map((s, i) => (
            <span key={s.id} className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-xs border ${i === 0 ? 'bg-indigo-50 border-indigo-300 text-indigo-800' : 'bg-gray-50 border-gray-300 text-gray-700 hover:border-indigo-300 cursor-pointer'}`}
              title={`${s.file.name} · ${(s.file.size / 1e9).toFixed(2)} GB · ${s.kind.toUpperCase()}${i === 0 ? ' · primary sample' : ' · click to make it the primary sample'}`}
              onClick={() => { if (i !== 0) makePrimary(s.id); }}>
              {i === 0 && <span title="primary sample">★</span>}{s.name}<button onClick={e => { e.stopPropagation(); removeSample(s.id); }} className="text-gray-400 hover:text-red-500" title="Remove">×</button>
            </span>
          ))}
          {fasta && <span className="px-2 py-0.5 rounded-full text-xs border bg-emerald-50 border-emerald-300 text-emerald-800" title={fasta.fa.name}>FASTA · {fasta.fa.name}</span>}
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
