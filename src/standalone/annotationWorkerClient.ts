/**
 * The page's side of the annotation worker (annotationWorker.ts): a name scanner that runs there, or here when the
 * browser cannot start it.
 */
import AnnotationWorker from './annotationWorker?worker&inline';
import { scanNames, type AnnotationFormat, type NameIndex } from './annotationScan';

export function scanNamesInWorker(file: Blob, gz: boolean, format: AnnotationFormat, onProgress?: (fraction: number) => void): Promise<NameIndex> {
  let worker: Worker;
  try { worker = new AnnotationWorker(); } catch (e) {
    console.warn('[sashimi] annotation names read on the page (no worker):', e);
    return scanNames(file, gz, format, onProgress);
  }
  return new Promise<NameIndex>((resolve, reject) => {
    worker.onmessage = (e: MessageEvent<{ progress?: number; index?: NameIndex; error?: string }>) => {
      const m = e.data;
      if (m.progress != null) { onProgress?.(m.progress); return; }
      worker.terminate();
      if (m.index) resolve(m.index); else reject(new Error(m.error || 'the annotation worker failed'));
    };
    worker.onerror = e => { worker.terminate(); reject(new Error(e.message || 'the annotation worker failed')); };
    worker.postMessage({ file, gz, format });
  });
}
