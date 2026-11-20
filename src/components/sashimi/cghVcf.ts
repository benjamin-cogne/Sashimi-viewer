/**
 * An SV / CNV VCF read as array-CGH calls: the records whose SVTYPE (or symbolic ALT) is DEL, DUP or CNV, placed by POS
 * and END (or SVLEN); the copy number from the first sample's CN, else <CNn>. BND, INV, INS and the other types are
 * skipped and counted. One sample read (the first); the build from ##reference or the ##contig lines. No probe count,
 * as for a BED. Classes and notes come back from INFO CGH_CLASS / CGH_NOTE (written by the export).
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { cghChrom, cghBuild, nextFileId, leafOf, type CghCnv, type CghFile, type CghKind, type CghSample } from './cghIntervals';
import { eolOf } from './cghMapping';

export const isVcf = (text: string) => /^﻿?##fileformat=VCF/i.test(text);
/** chr1's length in each build: a ##contig line with no assembly still tells it */
const CHR1: Record<number, GenomeBuild> = { 249250621: 'GRCh37', 248956422: 'GRCh38' };
const pctDecode = (t: string) => { try { return decodeURIComponent(t); } catch { return t; } };
const count = (n: number, one: string, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

export function parseVcf(text: string, fileName: string, fallbackSample?: string): CghFile {
  const eol = eolOf(text), lines = text.replace(/^﻿/, '').split(/\r\n|\n|\r/);
  const id = nextFileId('vcf'), leaf = leafOf(fileName);
  let build: GenomeBuild | null = null, buildText: string | null = null, samples: string[] = [];
  const meta: [string, string][] = [];
  const warnings: string[] = [];
  const other = new Map<string, number>();
  let bad = 0, neutral = 0, absent = 0, mito = 0, noCn = 0, firstSkip: string | null = null;
  const skip = (why: string) => { firstSkip ??= why; };
  let s: CghSample | null = null;
  const cnvs: CghCnv[] = [];
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    if (!l.trim()) continue;
    if (l.startsWith('##')) {
      const ref = l.match(/^##reference=(.*)$/i)?.[1], asm = l.match(/^##contig=<.*\bassembly=([^,>]+)/i)?.[1];
      const len = l.match(/^##contig=<ID=(?:chr)?1,.*\blength=(\d+)/i)?.[1];
      const b = cghBuild(ref) ?? cghBuild(asm) ?? (len ? CHR1[+len] ?? null : null);
      if (b && !build) { build = b; buildText = (ref ?? asm ?? `chr1 length ${len}`).trim(); }
      if (ref) meta.push(['reference', ref.trim()]);
      continue;
    }
    if (l.startsWith('#CHROM')) { samples = l.split('\t').slice(9); continue; }
    const f = l.split('\t');
    if (f.length < 8) { bad++; skip(`line ${i + 1}: fewer than 8 columns`); continue; }
    const info = new Map(f[7].split(';').map(x => { const k = x.indexOf('='); return k < 0 ? [x, ''] : [x.slice(0, k), x.slice(k + 1)]; }));
    const alt = f[4] ?? '', type = (info.get('SVTYPE') ?? alt.match(/^<([A-Z0-9]+)/i)?.[1] ?? '').toUpperCase();
    const fmt = f[8]?.split(':') ?? [], smp = f[9]?.split(':') ?? [];
    const fv = (k: string) => { const j = fmt.indexOf(k); return j >= 0 ? smp[j] : undefined; };
    const chrom = cghChrom(f[0]);
    if (chrom === 'chrM') { mito++; skip(`line ${i + 1}: chrM`); continue; }
    // the first sample does not carry it (0/0)
    const gt = fv('GT');
    if (gt && /^0([/|]0)*$/.test(gt)) { absent++; skip(`line ${i + 1}: not carried (GT ${gt})`); continue; }
    const cnText = fv('CN') ?? info.get('CN') ?? alt.match(/^<CN(\d+)>$/i)?.[1];
    const cn = cnText != null && /^\d+(\.\d+)?$/.test(cnText) ? Math.round(+cnText) : null;
    let kind: CghKind | null = null;
    if (/^(DEL|LOSS)/.test(type)) kind = 'loss';
    else if (/^(DUP|GAIN|AMP)/.test(type)) kind = 'gain';
    else if (type === 'CNV' || /^CN\d+$/.test(type)) {
      if (cn == null) { noCn++; skip(`line ${i + 1}: CNV without a copy number`); continue; }
      if (cn === 2) { neutral++; skip(`line ${i + 1}: copy number 2`); continue; }
      kind = cn < 2 ? 'loss' : 'gain';
    } else { const k = type || 'no SVTYPE'; other.set(k, (other.get(k) ?? 0) + 1); skip(`line ${i + 1}: ${k}`); continue; }
    const pos = Number(f[1]), endInfo = Number(info.get('END')), svlen = Math.abs(Number(info.get('SVLEN')));
    // a symbolic allele's POS is the base before the event
    const end = Number.isFinite(endInfo) && endInfo > 0 ? endInfo : Number.isFinite(svlen) && svlen > 0 ? pos + svlen : NaN;
    if (!chrom || !Number.isInteger(pos) || !Number.isFinite(end) || end < pos) { bad++; skip(`line ${i + 1}: no position, END or SVLEN`); continue; }
    const start = end > pos ? pos + 1 : pos;
    const name = samples[0] || fallbackSample || leaf.replace(/\.vcf(\.gz)?$/i, '');
    s ??= { key: `${id}:${name}`, name, fileId: id, declaredCalls: null, cnvs, sex: null };
    cnvs.push({
      id: `${id}:${i}`, sample: name, sampleKey: s.key, chrom, start, end, cytoband: '', genes: [], sizeKb: (end - start + 1) / 1000, kind,
      typeText: type === 'CNV' ? `CNV CN=${cn}` : type, probes: null, logRatio: null, iscn: null, copies: cn, copiesHi: null, amp: false,
      suppressed: false, classification: pctDecode(info.get('CGH_CLASS') ?? ''), classRest: '', note: pctDecode(info.get('CGH_NOTE') ?? ''), line: i,
    });
  }
  warnings.push('VCF: no probe count, the minimum number of probes does not apply');
  if (samples.length > 1) warnings.push(`${samples.length} samples in the VCF: only the first (${samples[0]}) is read`);
  if (other.size) warnings.push(`${[...other].map(([k, n]) => `${n} ${k}`).join(', ')} skipped (not a deletion, duplication or CNV)`);
  if (neutral) warnings.push(`${count(neutral, 'CNV')} of copy number 2 skipped`);
  if (noCn) warnings.push(`${count(noCn, 'CNV')} without a copy number skipped`);
  if (absent) warnings.push(`${count(absent, 'record')} not carried by ${samples[0] ?? 'the sample'} (GT 0/0) skipped`);
  if (mito) warnings.push(`${count(mito, 'record')} on chrM skipped`);
  if (bad) warnings.push(`${count(bad, 'record')} without a readable position skipped`);
  if (!s) throw new Error(`${fileName}: no deletion, duplication or CNV record (SVTYPE DEL / DUP / CNV with END or SVLEN) could be read.`);
  const skipped = bad + neutral + absent + mito + noCn + [...other.values()].reduce((a, b) => a + b, 0);
  return {
    id, name: fileName, saveAs: `${leaf.replace(/\.vcf(\.gz)?$/i, '').replace(/\.gz$/i, '')}.classified.vcf`, format: 'vcf', build, buildText, meta, lines, eol,
    classCol: null, classSuffix: '', headerLine: null, commentCol: null, samples: [s], warnings, delim: '\t', spec: null, fallback: fallbackSample, skipped, firstSkip,
  };
}
