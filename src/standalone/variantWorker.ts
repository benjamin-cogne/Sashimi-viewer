/**
 * Web Worker running full variant scans off the page's thread.
 *
 * A scan decodes every read of a window with its sequence (at 1 000× over a gene, over a million reads,
 * seconds of work however it is counted). On the page that work took the frames with it: a zoom-out that
 * widened the window froze the interface. Here it runs in its own LocalDataSource, on the same files
 * (File objects cross to a worker) and the same reference (a local FASTA, or the web APIs, fetched from
 * here), with the same code as the page would run — countVariants, allele counts kept per sample.
 *
 * Messages in: reference, sample, forget, scan {req, …}, methyl {req, …}, cancel {req}. Out: progress, done, error.
 * The CpG methylation of long reads (methylation.ts) is counted here too, for the same reason.
 */
import { LocalDataSource, type LocalSample, type ReferenceChoice } from './localSource';
import type { VariantScanOptions } from '../components/sashimi/datasource';

type In =
  | { type: 'reference'; reference: ReferenceChoice }
  | { type: 'sample'; sample: LocalSample }
  | { type: 'forget'; id: number }
  | { type: 'cancel'; req: number }
  | { type: 'methyl'; req: number; sampleId: number; chrom: string; start: number; end: number }
  | { type: 'scan'; req: number; sampleId: number; chrom: string; start: number; end: number; uniqueOnly: boolean; minVaf: number; opts: Omit<VariantScanOptions, 'signal' | 'onProgress'> };

const ds = new LocalDataSource({ build: 'GRCh38' });
const running = new Map<number, AbortController>();
const post = (m: unknown) => (self as unknown as Worker).postMessage(m);

self.onmessage = (e: MessageEvent<In>) => {
  const m = e.data;
  if (m.type === 'reference') ds.setReference(m.reference);
  else if (m.type === 'sample') ds.addSample(m.sample);
  else if (m.type === 'forget') ds.removeSample(m.id);
  else if (m.type === 'cancel') running.get(m.req)?.abort();
  else if (m.type === 'methyl') {
    // CpG methylation (methylation.ts): the window comes back as typed arrays, handed over without a copy
    const ctl = new AbortController();
    running.set(m.req, ctl);
    ds.countMethylation(m.sampleId, m.chrom, m.start, m.end, { signal: ctl.signal, onProgress: fraction => post({ type: 'progress', req: m.req, fraction }) })
      .then(result => (self as unknown as Worker).postMessage({ type: 'done', req: m.req, result }, [result.pos.buffer, ...result.mod.map(a => a.buffer), ...result.total.map(a => a.buffer)]))
      .catch((err: any) => post({ type: 'error', req: m.req, name: err?.name ?? 'Error', message: err?.message ?? String(err) }))
      .finally(() => running.delete(m.req));
  }
  else if (m.type === 'scan') {
    const ctl = new AbortController();
    running.set(m.req, ctl);
    ds.countVariants(m.sampleId, m.chrom, m.start, m.end, m.uniqueOnly, m.minVaf,
      { ...m.opts, signal: ctl.signal, onProgress: fraction => post({ type: 'progress', req: m.req, fraction }) })
      .then(result => post({ type: 'done', req: m.req, result }))
      .catch((err: any) => post({ type: 'error', req: m.req, name: err?.name ?? 'Error', message: err?.message ?? String(err) }))
      .finally(() => running.delete(m.req));
  }
};
