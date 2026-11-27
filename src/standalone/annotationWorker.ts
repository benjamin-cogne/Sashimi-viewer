/**
 * The name index of a tabix-indexed annotation file (annotationScan.ts scanNames), built away from the page's thread:
 * a whole human annotation is millions of lines. Asked with { file, gz, format }; answers { progress } as it reads,
 * then { index } or { error }.
 */
import { scanNames, type AnnotationFormat } from './annotationScan';

self.onmessage = async (e: MessageEvent<{ file: Blob; gz: boolean; format: AnnotationFormat }>) => {
  const { file, gz, format } = e.data;
  try {
    const index = await scanNames(file, gz, format, progress => self.postMessage({ progress }));
    self.postMessage({ index });
  } catch (err) {
    self.postMessage({ error: (err as Error)?.message || String(err) });
  }
};
