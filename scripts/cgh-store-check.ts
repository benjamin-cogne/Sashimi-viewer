// Store checks of the review fixes (invented data only): a zip whose table waits for its columns keeps its probes for
// it, a table read again with other columns keeps its probes, new probes during a reading, a worker error.
// The worker is a stand-in (globalThis.TestWorker). Run:
// node -e "require('esbuild').build({entryPoints:['scripts/cgh-store-check.ts'],bundle:true,platform:'node',outfile:'/tmp/cgh-store-check.cjs',logLevel:'warning',plugins:[{name:'w',setup(b){b.onResolve({filter:/[?]worker/},()=>({path:'w',namespace:'w'}));b.onLoad({filter:/.*/,namespace:'w'},()=>({contents:'export default class{constructor(){return new globalThis.TestWorker()}}',loader:'js'}))}}]})" && node /tmp/cgh-store-check.cjs
import { zipSync, strToU8 } from 'fflate';
import { cgh } from '../src/components/sashimi/cghStore';
import { parseCghFile } from '../src/components/sashimi/cghIntervals';
import { parseProbes } from '../src/components/sashimi/cghProbes';
import { detect } from '../src/components/sashimi/cghSegment';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = JSON.stringify(got) === JSON.stringify(want); if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', JSON.stringify(got), ok ? '' : `| want ${JSON.stringify(want)}`); };
const wait = (ms: number) => new Promise(r => setTimeout(r, ms));
/** the stand-in worker: reads after `delay` ms; `mode` 'err' answers the probe calls with an error, 'bad' also gives probes the page cannot segment */
const W = { delay: 30, mode: 'ok' as 'ok' | 'err' | 'bad' };
(globalThis as any).TestWorker = class {
  onmessage: ((e: { data: unknown }) => void) | null = null; onerror: (() => void) | null = null;
  postMessage(m: any) {
    setTimeout(() => {
      let data: unknown;
      if (m.parse) {
        const p = parseProbes(m.parse);
        data = { id: m.id, result: { ...p, byChrom: new Map([...p.byChrom].map(([c, x]) => [c, W.mode === 'bad' ? { ...x, lr: undefined, name: '' } : { ...x, name: x.name.join('\n') }])) } };
      } else data = W.mode === 'ok' ? { id: m.id, result: detect({ sample: '', build: m.build, n: 0, byChrom: m.byChrom, dlrs: null }, m.opts, m.male) } : { id: m.id, error: 'boom' };
      this.onmessage?.({ data });
    }, W.delay);
  }
  terminate() { /* */ }
};
let seed = 7; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
/** an invented probe file: n probes on chr1, a loss of 20 probes at -0.9 */
const probeText = (sample: string, n: number) => ['AMADID\t000000\thg19', `Probe Name\tChrName\tStart\tStop\tFeatureNum\tDescription\tName of Gene\tNormalized log ratio(${sample})\t`,
  ...Array.from({ length: n }, (_, i) => `P${i}\t1\t${1e6 + i * 5000}\t${1e6 + i * 5000 + 59}\t${i}\t\t\t${(gauss() * 0.13 + (i >= 400 && i < 420 ? -0.9 : 0)).toFixed(4)}`)].join('\n') + '\n';
const file = (name: string, b: Uint8Array | string) => new File([typeof b === 'string' ? strToU8(b) : b], name);
/** a table whose columns are not recognised (the mapping dialog) */
const odd = 'Alpha\tBeta\tGamma\tDelta\nchr1\t3000000\t3100000\tfoo\nchr2\t5000000\t5200000\tbar\n';
const spec = { cols: { chrom: 0, start: 1, stop: 2, pos: -1, type: 3, copies: -1, iscn: -1, lr: -1, probes: -1, size: -1, band: -1, genes: -1, cls: -1, note: -1, sample: -1, sup: -1 },
  opts: { delim: '\t', headerRow: 0, zeroBased: false, sizeUnit: 'kb' as const, build: null, typeMap: { foo: 'loss' as const, bar: 'gain' as const } } };

(async () => {
  // review 4: a zip whose table waits for its columns, with a probe file: no probes-only patient; the probes go to the mapped table
  await cgh.add([file('DEMOZ.zip', zipSync({ 'z/calls.txt': strToU8(odd), 'z/CGH_Log_Ratio.txt': strToU8(probeText('DEMOZ', 1000)) }))]);
  const s1 = cgh.get();
  eq('zip: table asked, probes held with it', [s1.toMap.length, s1.toMap[0]?.probes?.length, s1.pending?.files.length ?? 0], [1, 1, 0]);
  const m1 = parseCghFile(s1.toMap[0].text, s1.toMap[0].name, s1.toMap[0].fallback, spec);
  cgh.mapped(m1);
  const s2 = cgh.get(), k1 = m1.samples[0].key;
  eq('mapped: one patient, with the probes', [s2.pending?.files.length, s2.pending?.files[0].format, cgh.hasProbes(k1)], [1, 'interval', true]);
  // review 1: read again with other columns: new keys, the probes moved to them
  cgh.remap(s2.pending!.files[0]);
  const h = cgh.get().toMap[0], m2 = parseCghFile(h.text, h.name, h.fallback, spec);
  cgh.mapped(m2);
  const k2 = m2.samples[0].key;
  eq('remapped: probes under the new key only', [k1 !== k2, cgh.hasProbes(k2), cgh.hasProbes(k1)], [true, true, false]);
  cgh.cancel();
  // review 2: new probes while the former ones are being read: the new ones read, the probe calls made
  W.delay = 60;
  await cgh.add([file('DEMOR.zip', zipSync({ 'r/calls.txt': strToU8('Chromosome\tStart\tEnd\tType\nchr1\t3000000\t3100000\tDeletion\n'), 'r/CGH_Log_Ratio.txt': strToU8(probeText('DEMOR', 1000)) }))]);
  const k3 = cgh.get().pending!.files[0].samples[0].key;
  cgh.probes(k3);
  await wait(10);
  await cgh.add([file('CGH_Log_Ratio_DEMOR.txt', probeText('DEMOR', 800))]);
  cgh.runProbeCalls([k3]);
  await wait(400);
  eq('new probes read, calls made', [cgh.probes(k3)?.n, !!cgh.probeCallsOf(k3), cgh.probeCallsOf(k3)?.events.length], [800, true, 1]);
  cgh.cancel();
  // review 5: the worker answers with an error and the page's own try fails too: settled, said, the queue goes on
  W.delay = 5; W.mode = 'bad';
  await cgh.add([file('CGH_Log_Ratio_DEMOE.txt', probeText('DEMOE', 600))]);
  const k4 = cgh.get().pending!.files[0].samples[0].key;
  await wait(200);
  eq('error said', [cgh.probeErrorOf(k4)?.length ? 'set' : null, cgh.probeCallsOf(k4)], ['set', null]);
  W.mode = 'err';
  await cgh.add([file('CGH_Log_Ratio_DEMOF.txt', probeText('DEMOF', 600))]);
  const k5 = cgh.get().pending!.files.find(f => f.samples[0].name === 'DEMOF')!.samples[0].key;
  await wait(200);
  eq('worker error, page fallback: calls made', [cgh.probeErrorOf(k5), cgh.probeCallsOf(k5)?.events.length], [null, 1]);
  console.log(bad ? `${bad} FAILED` : 'all ok'); process.exit(bad ? 1 : 0);
})();
