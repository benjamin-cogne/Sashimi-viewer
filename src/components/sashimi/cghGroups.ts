/**
 * CNV groups: several calls of one patient, of one direction (loss / gain), that together cover a whole chromosome arm or a
 * whole chromosome. A scanner often cuts an aneuploidy at the centromere (8p + 8q), and ADM-2 also reports smaller calls
 * nested in a larger one. The group gives the event's own metrics (log ratio weighted by probes, nested calls left out of
 * the weights and of the counts), copies / mosaic fraction (cghMosaic.ts), an ISCN and a suggestion in words
 * (trisomy 8, partial monosomy 4p, monosomy X or a sex mismatch…). A suggestion only, never a class.
 */
import type { GenomeBuild } from '../../standalone/ensembl';
import { chromLength } from '../../standalone/assembly';
import type { CghCnv, CghKind, CghSex } from './cghIntervals';
import { estCopies, mosaicOf, parOf, MOSAIC_TAG_BELOW } from './cghMosaic';
import { callIscn, iscnBands } from './cghIscn';

/** the share of an arm's probe-covered span (centromere, heterochromatin, acrocentric p arms, PARs left out) the calls must cover */
export const GROUP_COVER = 0.9;
const ACROCENTRIC = new Set(['13', '14', '15', '21', '22']);
/**
 * Per build, per chromosome, in units of 100 kb: centromere start, centromere end, then the spans without probes on arrays
 * (UCSC cytoBand `acen`, `gvar`, `stalk` merged; `E` = the chromosome's end). Read from the UCSC API on 2026-11-25. Used when
 * the patient has no probe file; with one, an arm's span is that of its probes (`GroupOpts.probes`).
 */
const BLIND: Record<GenomeBuild, string> = {
  GRCh37: '1:1215,1289,1215-1426 2:905,968,905-968 3:879,939,879-983 4:482,527,482-527 5:461,507,461-507 6:587,633,587-633 7:580,617,580-617 8:431,481,431-481 9:473,507,473-659 10:380,423,380-423 11:516,557,516-557 12:333,382,333-382 13:163,195,0-195 14:161,191,0-191 15:158,207,0-207 16:346,386,346-470 17:222,258,222-258 18:154,190,154-190 19:244,286,200-324 20:256,294,256-294 21:109,143,0-143 22:122,179,0-179 X:581,630,581-630 Y:116,134,116-134,288-E',
  GRCh38: '1:1217,1251,1217-1432 2:918,960,918-960 3:878,940,878-986 4:482,518,482-518 5:461,514,461-514 6:585,626,585-626 7:581,621,581-621 8:432,472,432-472 9:422,455,422-615 10:380,416,380-416 11:510,558,510-558 12:332,378,332-378 13:165,189,0-189 14:161,182,0-182 15:175,205,0-205 16:353,384,353-470 17:227,274,227-274 18:154,215,154-215 19:242,281,199-319 20:257,304,257-304 21:109,130,0-130 22:137,174,0-174 X:581,638,581-638 Y:103,106,103-106,266-E',
};
interface ChromInfo { len: number; cen: [number, number]; blind: [number, number][] }
const infoCache = new Map<string, ChromInfo | null>();
/** a chromosome's centromere and blind spans (0-based half-open, the PARs added on X and Y); null when unknown */
export function chromInfo(build: GenomeBuild, chrom: string): ChromInfo | null {
  const k = `${build}|${chrom}`;
  if (infoCache.has(k)) return infoCache.get(k)!;
  const len = chromLength(build, chrom), n = chrom.replace(/^chr/, '');
  const row = BLIND[build]?.split(' ').find(r => r.startsWith(`${n}:`));
  let out: ChromInfo | null = null;
  if (len && row) {
    const v = (t: string) => (t === 'E' ? len : Number(t) * 1e5);
    const [cs, ce, ...spans] = row.slice(n.length + 1).split(',');
    const blind = spans.map(x => x.split('-').map(v) as [number, number]);
    for (const [a, b] of parOf(build)[`chr${n}`] ?? []) blind.push([a - 1, b]);
    out = { len, cen: [v(cs), v(ce)], blind };
  }
  infoCache.set(k, out);
  return out;
}

/** the length of [a, z) outside the blind spans */
const seen = (a: number, z: number, blind: [number, number][]) => {
  let n = Math.max(0, z - a);
  // blind spans may overlap (a PAR in a gvar): merged first
  const m = blind.map(([x, y]) => [Math.max(a, x), Math.min(z, y)]).filter(([x, y]) => y > x).sort((p, q) => p[0] - q[0]);
  let end = -Infinity;
  for (const [x, y] of m) { const s = Math.max(x, end); if (y > s) n -= y - s; end = Math.max(end, y); }
  return n;
};
/** fewer probes than this on an arm: its span is taken from the cytogenetic bands instead */
const MIN_ARM_PROBES = 10;
/** the probes of a chromosome (sorted by start, NaN kept): the positions only are read */
export interface GroupProbes { start: ArrayLike<number>; end: ArrayLike<number>; lr: ArrayLike<number> }
/**
 * The share of [a, z)'s probe-covered span covered by the intervals (0-based half-open): with the patient's probes, the
 * share of the arm's probes (non-NaN, outside the blind spans, by midpoint) lying in an interval; else the share of the
 * arm's length outside the blind spans.
 */
function coverOf(a: number, z: number, blind: [number, number][], iv: [number, number][], pr?: GroupProbes): number {
  if (pr) {
    let tot = 0, got = 0;
    for (let i = 0; i < pr.start.length; i++) {
      const m = (pr.start[i] + pr.end[i]) / 2 - 1;
      if (m < a || m >= z || Number.isNaN(pr.lr[i]) || blind.some(([x, y]) => m >= x && m < y)) continue;
      tot++;
      if (iv.some(([x, y]) => m >= x && m < y)) got++;
    }
    if (tot >= MIN_ARM_PROBES) return got / tot;
  }
  const tot = seen(a, z, blind);
  if (tot <= 0) return 0;
  const m = iv.map(([x, y]) => [Math.max(a, x), Math.min(z, y)] as [number, number]).filter(([x, y]) => y > x).sort((p, q) => p[0] - q[0]);
  let got = 0, cur: [number, number] | null = null;
  for (const r of m) {
    if (cur && r[0] <= cur[1]) cur[1] = Math.max(cur[1], r[1]);
    else { if (cur) got += seen(cur[0], cur[1], blind); cur = [r[0], r[1]]; }
  }
  if (cur) got += seen(cur[0], cur[1], blind);
  return got / tot;
}

export interface CnvGroup {
  /** `grp|<sample key>|<chrom>|<kind>`: its class and note are kept under it */
  id: string;
  sampleKey: string; sample: string; chrom: string; kind: CghKind;
  /** the whole chromosome, or one arm */
  scope: 'chr' | 'p' | 'q';
  /** every member, by position; `nested`: the ids of those inside a larger member (not weighted, not counted) */
  members: CghCnv[]; nested: Set<string>;
  start: number; end: number; spanKb: number;
  /** the probes of the members not nested (null when a member has no count) */
  probes: number | null;
  /** the members' mean log ratios weighted by their probes (by their size when a count is missing), nested left out */
  logRatio: number | null;
  /** the arms covered (8p, 8q; an acrocentric's q arm alone) and the share of each */
  arms: { arm: string; cover: number }[];
  /** the copy number (the members' own when they agree, else from the log ratio) and the mosaic fraction (< 80 %); a
   * reference mismatch: the stated sex's normal copies (x1), no mosaic */
  copies: number | null; mosaic: number | null;
  /** a male's whole X near −1, or his Y far above 0: what a normal male gives against a female reference (`MISMATCH_TIP`) */
  mismatch: boolean;
  /** the copies as the table writes them: x3, x3 mosaic ~29% */
  copiesText: string;
  /** `8 (whole chromosome)`, `4p (whole arm)` */
  label: string;
  /** in words, a suggestion only */
  suggestion: string;
  /** the group as one call: chromosome, span, mean log ratio, probes (for the ISCN, the mosaic estimate, the tooltips) */
  as: CghCnv;
}

const WORD: Record<number, string> = { 0: 'nullisomy', 1: 'monosomy', 3: 'trisomy', 4: 'tetrasomy', 5: 'pentasomy', 6: 'hexasomy' };
const pct = (f: number) => `~${Math.round(f * 100)} %`;
/** what a group (and its calls) read as a reference mismatch show in place of copy and mosaic estimates */
export const MISMATCH_TIP = 'A male\'s whole X near log ratio −1 (or his Y far above 0) is what a normal male gives against a female reference. Read against the stated sex, the copy number is x1 (normal): no loss or gain as such, and no mosaic estimate (it assumes a sex-matched reference). Check the sex given and the reference used.';
const MISMATCH = 'or a sex mismatch with the reference: check the sex given';

/** the words of a group (a suggestion: the lab decides) */
export function groupSuggestion(g: Pick<CnvGroup, 'chrom' | 'kind' | 'scope' | 'copies' | 'mosaic' | 'logRatio'>, sex: CghSex | undefined): string {
  const n = g.chrom.replace(/^chr/, ''), where = g.scope === 'chr' ? n : `${n}${g.scope}`, k = g.copies;
  const mos = g.mosaic != null ? `mosaic ` : '', tail = g.mosaic != null ? ` (${pct(g.mosaic)})` : '';
  const auto = () => {
    const w = k != null ? WORD[k] : null;
    if (!w) return `${g.scope === 'chr' ? 'whole chromosome' : 'whole arm'} ${where} ${g.kind}${k != null ? ` (x${k})` : ''}`;
    return `${mos}${g.scope === 'chr' ? '' : 'partial '}${w} ${where}${tail}`;
  };
  if (n !== 'X' && n !== 'Y') return auto();
  const loss = g.kind === 'loss';
  // chrY, whole or its q arm: read against the sex (a female reference gives a male's Y far above 0)
  if (n === 'Y') {
    const far = (g.logRatio ?? 0) >= 1.5;
    if (sex === 'female') return `Y ${g.kind} in a female: a sex mismatch with the reference or a wrong sex given — check the sex`;
    if (sex === 'male') return loss ? `${mos}loss of Y${tail} (45,X or mosaic loss of Y) — ${MISMATCH}`
      : far ? `Y far above the reference: expected for a male against a female reference (sex-mismatched hybridisation), or a wrong sex given — check the sex and the reference; not a gain of Y as such`
      : `${g.mosaic != null ? `mosaic XYY${tail}` : 'XYY? (47,XYY)'} — ${MISMATCH}`;
    return loss ? `sex not given: Y loss — a female against a male reference, or a loss of Y in a male — set the patient's sex`
      : `sex not given: Y gain — a male against a female reference${far ? '' : ', or XYY in a male'}, or Y material in a female — set the patient's sex`;
  }
  if (g.scope !== 'chr') return `${auto()} (sex chromosome: read with the patient's sex)`;
  // chrX, whole
  {
    if (sex === 'female') {
      if (loss) return g.mosaic != null ? `mosaic monosomy X (${pct(g.mosaic)}; 45,X/46,XX)` : k === 0 ? `whole X loss (x0) — ${MISMATCH}` : `monosomy X (45,X) — ${MISMATCH}`;
      return g.mosaic != null ? `mosaic trisomy X (${pct(g.mosaic)})` : `${k === 4 ? 'tetrasomy X (48,XXXX)' : k != null && k >= 5 ? `X x${k}` : 'trisomy X (47,XXX)'} — ${MISMATCH}`;
    }
    if (sex === 'male') {
      if (loss) return `whole X lower in a male: expected for a male against a female reference (sex-mismatched hybridisation), or a wrong sex given — check the sex and the reference; not a loss of the X as such`;
      return g.mosaic != null ? `mosaic XXY (${pct(g.mosaic)}; 47,XXY/46,XY)?` : `${k === 3 ? 'XXXY? (48,XXXY)' : 'XXY? (47,XXY)'} — ${MISMATCH} (a female sample?)`;
    }
    return loss
      ? `sex not given: monosomy X (45,X) in a female, or a male against a female reference (sex mismatch) — set the patient's sex`
      : `sex not given: XXY (47,XXY) in a male, or trisomy X (47,XXX) in a female, or a sex mismatch with the reference — set the patient's sex`;
  }
}

export interface GroupOpts {
  build: GenomeBuild; sex: CghSex | undefined; eligible: (c: CghCnv) => boolean;
  /** the patient's probes by chromosome, when read: an arm's span is then that of its probes */
  probes?: Map<string, GroupProbes> | null;
}
/** the groups of one patient's calls */
export function cnvGroups(cnvs: CghCnv[], o: GroupOpts): CnvGroup[] {
  const out: CnvGroup[] = [];
  const by = new Map<string, CghCnv[]>();
  for (const c of cnvs) if (!c.fromProbes && o.eligible(c)) { const k = `${c.chrom}|${c.kind}`; const l = by.get(k) ?? []; l.push(c); by.set(k, l); }
  for (const [k, list] of by) {
    if (list.length < 2) continue;
    const [chrom, kind] = k.split('|') as [string, CghKind];
    const ci = chromInfo(o.build, chrom); if (!ci) continue;
    const n = chrom.replace(/^chr/, ''), acro = ACROCENTRIC.has(n);
    // nested: inside a member that starts before it (sorted by start, the longest first)
    const sorted = list.slice().sort((a, b) => a.start - b.start || b.end - a.end), nested = new Set<string>();
    let reach = -Infinity;
    for (const c of sorted) { if (c.end <= reach) nested.add(c.id); else reach = c.end; }
    const top = sorted.filter(c => !nested.has(c.id)), iv = top.map(c => [c.start - 1, c.end] as [number, number]);
    const pr = o.probes?.get(chrom);
    const pc = acro ? 0 : coverOf(0, ci.cen[0], ci.blind, iv, pr), qc = coverOf(ci.cen[1], ci.len, ci.blind, iv, pr);
    const p = !acro && pc >= GROUP_COVER, q = qc >= GROUP_COVER;
    const scope: CnvGroup['scope'] | null = (p || acro) && q ? 'chr' : p ? 'p' : q ? 'q' : null;
    if (!scope) continue;
    // the members: on the arm(s) covered, with the calls nested in them
    const onArm = (c: CghCnv) => scope === 'chr' || (scope === 'p' ? c.start - 1 < ci.cen[0] : c.end > ci.cen[1]);
    const tops = top.filter(onArm), topIds = new Set(tops.map(c => c.id));
    const members = sorted.filter(c => topIds.has(c.id) || (nested.has(c.id) && tops.some(t => c.start >= t.start && c.end <= t.end)));
    if (members.length < 2) continue;
    const start = Math.min(...tops.map(c => c.start)), end = Math.max(...tops.map(c => c.end));
    const probes = tops.every(c => c.probes != null) ? tops.reduce((a, c) => a + c.probes!, 0) : null;
    const withLr = tops.filter(c => c.logRatio != null), w = (c: CghCnv) => (probes != null ? c.probes! : c.sizeKb);
    const wsum = withLr.reduce((a, c) => a + w(c), 0);
    const logRatio = withLr.length && wsum > 0 ? withLr.reduce((a, c) => a + c.logRatio! * w(c), 0) / wsum : null;
    const cps = new Set(tops.map(c => (c.copies != null && c.copiesHi == null ? c.copies : null)));
    const common = cps.size === 1 ? [...cps][0] : null;
    const as: CghCnv = { ...tops[0], id: `grp|${tops[0].sampleKey}|${chrom}|${kind}`, start, end, sizeKb: (end - start + 1) / 1000, probes, logRatio,
      copies: common, copiesHi: null, amp: tops.some(c => c.amp), iscn: null, cytoband: '', genes: [], note: '', classification: '', classRest: '' };
    let copies: number | null = common, mosaic: number | null = null;
    // the cases groupSuggestion reads as expected for a male against a female reference: the stated sex's copies, x1
    const mismatch = o.sex === 'male' && ((n === 'X' && scope === 'chr' && kind === 'loss') || (n === 'Y' && kind === 'gain' && (logRatio ?? 0) >= 1.5));
    if (mismatch) copies = 1;
    else if (common != null) { const m = mosaicOf(as, o.sex, o.build); mosaic = m?.ok && m.f < MOSAIC_TAG_BELOW ? m.f : null; }
    else { const e = estCopies(as, o.sex, o.build); copies = e?.n ?? null; mosaic = e?.mosaic ?? null; }
    const arms = [...(!acro && (scope === 'chr' || scope === 'p') ? [{ arm: `${n}p`, cover: pc }] : []), ...(scope !== 'p' ? [{ arm: `${n}q`, cover: qc }] : [])];
    const g = { chrom, kind, scope, copies, mosaic, logRatio };
    out.push({ id: as.id, sampleKey: as.sampleKey, sample: as.sample, chrom, kind, scope, members, nested: new Set([...nested].filter(id => members.some(c => c.id === id))),
      start, end, spanKb: as.sizeKb, probes, logRatio, arms, copies, mosaic, mismatch,
      copiesText: mismatch ? 'x1 = normal male (female reference?)' : copies == null ? (as.amp ? 'amp' : '') : `x${copies}${mosaic != null ? ` mosaic ${pct(mosaic).replace(' ', '')}` : ''}`,
      label: `${scope === 'chr' ? n : `${n}${scope}`} (${scope === 'chr' ? 'whole chromosome' : 'whole arm'})`, suggestion: groupSuggestion(g, o.sex), as });
  }
  const ORDER = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY'];
  return out.sort((a, b) => ORDER.indexOf(a.chrom) - ORDER.indexOf(b.chrom) || a.start - b.start);
}

/** the group's ISCN: `(8)x3[0.29]` for a whole chromosome; an arm by its bands (first member's first band, last member's last) */
export function groupIscn(g: CnvGroup, build: string, sex: CghSex | undefined, bandOf: (c: CghCnv) => string) {
  // a normal male's X / Y against a female reference: nothing abnormal to write
  if (g.mismatch) return { text: null, part: null, file: null, differs: false };
  const first = g.members.find(c => c.start === g.start) ?? g.members[0], last = g.members.reduce((a, c) => (c.end > a.end ? c : a), g.members[0]);
  const band = iscnBands(`${bandOf(first)} ${bandOf(last)}`);
  const iscn = callIscn(g.as, build, band, sex, g.scope === 'chr');
  return iscn;
}
