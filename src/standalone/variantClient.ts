/**
 * The page's side of the variant worker (variantWorker.ts): a VariantScanner that a LocalDataSource hands its
 * full variant scans to. The worker is built into the page (Vite's inline worker, so the single HTML file
 * still works offline); samples and the reference are sent to it once, when a scan first needs them.
 * If the worker cannot be created, or dies, the data source is left to scan on the page as before.
 */
import VariantWorker from './variantWorker?worker&inline';
import type { LocalDataSource, VariantScanner, ReferenceChoice } from './localSource';
import type { VariantScan, VariantScanOptions } from '../components/sashimi/datasource';
import type { MethylWindow } from './methylation';

export function attachVariantWorker(ds: LocalDataSource): boolean {
  let worker: Worker;
  try { worker = new VariantWorker(); } catch (e) { console.warn('[sashimi] variant scans stay on the page (no worker):', e); return false; }
  let req = 0;
  const pending = new Map<number, { resolve: (r: any) => void; reject: (e: Error) => void; onProgress?: (f: number) => void }>();
  /** the sample files and the reference the worker has, by identity */
  const sent = new Map<number, File>();
  let refSent: ReferenceChoice | null = null;
  const abortError = () => new DOMException('Variant scan cancelled', 'AbortError');

  worker.onmessage = (e: MessageEvent<any>) => {
    const m = e.data, p = pending.get(m.req);
    if (!p) return;
    if (m.type === 'progress') p.onProgress?.(m.fraction);
    else if (m.type === 'done') { pending.delete(m.req); p.resolve(m.result); }
    else if (m.type === 'error') { pending.delete(m.req); p.reject(m.name === 'AbortError' ? abortError() : new Error(m.message)); }
  };
  worker.onerror = (e: ErrorEvent) => {
    // a worker that fails is dropped: what is running fails, the next scans run on the page
    console.warn('[sashimi] variant worker failed, scans go back to the page:', e.message);
    for (const p of pending.values()) p.reject(new Error(`variant worker failed: ${e.message}`));
    pending.clear();
    if (ds.variantScanner === scanner) ds.variantScanner = undefined;
    worker.terminate();
  };

  const scanner: VariantScanner = {
    scan(sampleId: number, chrom: string, start: number, end: number, uniqueOnly: boolean, minVaf: number, opts?: VariantScanOptions): Promise<VariantScan> {
      const s = ds.sampleFiles(sampleId);
      if (!s) return Promise.reject(new Error('Sample not found'));
      if (opts?.signal?.aborted) return Promise.reject(abortError());
      if (refSent !== ds.reference) { worker.postMessage({ type: 'reference', reference: ds.reference }); refSent = ds.reference; }
      if (sent.get(sampleId) !== s.file) {
        worker.postMessage({ type: 'sample', sample: { id: s.id, name: s.name, kind: s.kind, file: s.file, index: s.index } });
        sent.set(sampleId, s.file);
      }
      const id = ++req;
      return new Promise<VariantScan>((resolve, reject) => {
        pending.set(id, { resolve, reject, onProgress: opts?.onProgress });
        opts?.signal?.addEventListener('abort', () => worker.postMessage({ type: 'cancel', req: id }), { once: true });
        worker.postMessage({ type: 'scan', req: id, sampleId, chrom, start, end, uniqueOnly, minVaf,
          opts: { longReadMinIndel: opts?.longReadMinIndel, longReadMinVaf: opts?.longReadMinVaf, haplotypes: opts?.haplotypes, phaseSource: opts?.phaseSource } });
      });
    },
    methyl(sampleId: number, chrom: string, start: number, end: number, opts?: { signal?: AbortSignal; onProgress?: (fraction: number) => void }): Promise<MethylWindow> {
      const s = ds.sampleFiles(sampleId);
      if (!s) return Promise.reject(new Error('Sample not found'));
      if (opts?.signal?.aborted) return Promise.reject(abortError());
      if (refSent !== ds.reference) { worker.postMessage({ type: 'reference', reference: ds.reference }); refSent = ds.reference; }
      if (sent.get(sampleId) !== s.file) {
        worker.postMessage({ type: 'sample', sample: { id: s.id, name: s.name, kind: s.kind, file: s.file, index: s.index } });
        sent.set(sampleId, s.file);
      }
      const id = ++req;
      return new Promise<MethylWindow>((resolve, reject) => {
        pending.set(id, { resolve, reject, onProgress: opts?.onProgress });
        opts?.signal?.addEventListener('abort', () => worker.postMessage({ type: 'cancel', req: id }), { once: true });
        worker.postMessage({ type: 'methyl', req: id, sampleId, chrom, start, end });
      });
    },
    release(what: 'methylation' | 'variants') { worker.postMessage({ type: 'release', what }); },
    forget(sampleId: number) { if (sent.delete(sampleId)) worker.postMessage({ type: 'forget', id: sampleId }); },
    reference(reference: ReferenceChoice) { worker.postMessage({ type: 'reference', reference }); refSent = reference; },
  };
  ds.variantScanner = scanner;
  return true;
}
