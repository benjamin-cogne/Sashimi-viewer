/**
 * Web Worker reading a probe log-ratio file (cghProbes.ts `parseProbes`) and computing the probe calls (cghSegment.ts
 * `detect`) off the page's thread: a 1M-probe array takes seconds each, which on the page froze it. Same code as on the
 * page. In: { id, parse: bytes, spec? } or { id, byChrom, build, opts, male } (starts, ends, log ratios per chromosome, names
 * left out). Out: { id, result } (the probes with each chromosome's names joined by newlines, or the calls) or { id, error }.
 */
import { detect, type SegOpts } from './cghSegment';
import { parseProbes, type CghProbeChrom, type ProbeSpec } from './cghProbes';
import type { GenomeBuild } from '../../standalone/ensembl';

type In = { id: number; parse: Uint8Array; spec?: ProbeSpec } | { id: number; parse?: undefined; byChrom: Map<string, CghProbeChrom>; build: GenomeBuild | null; opts: SegOpts; male: boolean };
const post = (m: unknown) => (self as unknown as Worker).postMessage(m);
self.onmessage = (e: MessageEvent<In>) => {
  const d = e.data;
  try {
    if (d.parse) {
      const p = parseProbes(d.parse, d.spec);
      post({ id: d.id, result: { ...p, byChrom: new Map([...p.byChrom].map(([c, x]) => [c, { ...x, name: x.name.join('\n') }])) } });
    } else post({ id: d.id, result: detect({ sample: '', build: d.build, n: 0, byChrom: d.byChrom, dlrs: null }, d.opts, d.male) });
  } catch (x) { post({ id: d.id, error: (x as Error).message }); }
};
