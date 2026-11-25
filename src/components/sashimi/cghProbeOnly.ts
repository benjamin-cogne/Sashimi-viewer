/**
 * "Probes only": a CGH_Log_Ratio.txt added with no matching patient becomes a patient of its own, whose calls are the
 * probe calls (cghSegment.ts, the current settings) made into ordinary calls (annotated, classified, exported to a new
 * tab-separated file `<name>.probe-calls.txt`). Call ids are position-based (`pc|<key>|chr:start-end`): a class or a
 * note stays on a call whose position is unchanged when the calls are computed again (settings, sex). Pure.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { cghBuild, nextFileId, leafOf, type CghCnv, type CghFile } from './cghIntervals';
import type { ProbeEvent, SegOpts } from './cghSegment';

/** the banner of a probe-only patient (its tab and the export's header) */
export const PROBE_ONLY_WARNING = 'Calls detected from the probes by this page (CBS), not by CytoGenomics: to confirm in the dedicated software.';
export const isProbeOnly = (f: CghFile) => f.format === 'probes';

/** the build an AMADID / header line of a probe file names (hg19 / hg38), else null */
export function probeBuildOf(head: string): GenomeBuild | null {
  for (const l of head.split(/\r\n|\n|\r/).slice(0, 30)) { if (/^probe\s*name\t/i.test(l)) break; const b = cghBuild(l); if (b) return b; }
  return null;
}

/** a patient made from a probe file alone: no line, no call until the probe calls are computed */
export function probeOnlyFile(fileName: string, sample: string, build: GenomeBuild | null): CghFile {
  const id = nextFileId('probes'), name = sample || leafOf(fileName).replace(/\.(txt|tsv)(\.gz)?$/i, '');
  return {
    id, name: leafOf(fileName), saveAs: `${name.replace(/[\\/:*?"<>|]+/g, '_')}.probe-calls.txt`, format: 'probes', build, buildText: build ? (build === 'GRCh37' ? 'hg19' : 'hg38') : null,
    meta: [['Source probe file', leafOf(fileName)]], lines: [], eol: '\n', classCol: null, classSuffix: '', headerLine: null, commentCol: null,
    samples: [{ key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs: [], sex: null }],
    warnings: build ? [] : ['the probe file names no build (AMADID line): read as GRCh37'], delim: '\t', spec: null, skipped: 0, firstSkip: null,
  };
}

/** the probe calls as the patient's calls (copies, band and ISCN come from the log ratio and the UCSC bands, as for any call without them) */
export const probeOnlyCnvs = (events: ProbeEvent[], key: string, sample: string): CghCnv[] => events.map((e, i) => ({
  id: `pc|${key}|${e.chrom}:${e.start}-${e.end}`, sample, sampleKey: key, chrom: e.chrom, start: e.start, end: e.end, cytoband: '', genes: [],
  sizeKb: (e.end - e.start + 1) / 1000, kind: e.kind, typeText: e.kind === 'loss' ? 'Loss (probes)' : 'Gain (probes)', probes: e.n, logRatio: e.mean,
  iscn: null, copies: null, copiesHi: null, amp: false, suppressed: false, classification: '', classRest: '', note: '', line: i }));

const SENS: Record<number, string> = { 6: 'high', 8: 'normal', 10: 'low' };
export const optsText = (o: SegOpts) => `CBS, sensitivity ${SENS[o.t] ?? `t ${o.t}`} (t = ${o.t}), probes >= ${o.minProbes}, |mean log ratio| >= ${o.minLr}`;

/** the export: header lines (#), then one row per call */
export function exportProbeCalls(file: CghFile, o: SegOpts, at: Date, cell: { band: (c: CghCnv) => string; copies: (c: CghCnv) => string; iscn: (c: CghCnv) => string; cls: (c: CghCnv) => string; note: (c: CghCnv) => string }): string {
  const flat = (t: string) => t.replace(/[\t\r\n]+/g, ' ').trim();
  const out = [
    `# Source probe file: ${leafOf(file.name)}`,
    `# Sample: ${file.samples.map(x => x.name).join(', ')}`,
    `# Build: ${file.build ?? 'GRCh37 (none stated)'}`,
    `# Detection: ${optsText(o)}`,
    `# Date: ${at.toISOString().slice(0, 19).replace('T', ' ')}`,
    `# ${PROBE_ONLY_WARNING}`,
    ['Chromosome', 'Start', 'Stop', 'Cytoband', 'Size(kb)', 'Type', '#Probes', 'Mean Log Ratio', 'Copies', 'ISCN', 'Classification', 'Comment'].join('\t'),
  ];
  const cs = file.samples.flatMap(x => x.cnvs), noBand = cs.filter(c => !cell.band(c)).length;
  // bands read from UCSC when the page exported: a row left without one is counted
  if (noBand) out.splice(out.length - 1, 0, `# Rows without a cytoband (UCSC cytobands not read): ${noBand}`);
  for (const c of cs) out.push([c.chrom, c.start, c.end, cell.band(c), c.sizeKb.toFixed(2), c.kind === 'loss' ? 'Loss' : 'Gain', c.probes ?? '', c.logRatio?.toFixed(3) ?? '',
    cell.copies(c), cell.iscn(c), flat(cell.cls(c)), flat(cell.note(c))].join('\t'));
  return out.join('\n') + '\n';
}
