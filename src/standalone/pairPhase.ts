/**
 * Phase of two variant sites from the reads: in cis (on the same molecule), in trans (on the two haplotypes), or
 * undetermined, with the fragments that say so. The pairwise form of read-backed phasing (WhatsHap, HapCUT2), the
 * question asked of a possible compound heterozygote.
 *
 *  1. Reads: a small window around each site, never the stretch between them, so two exons 200 kb apart cost what two
 *     neighbours do. The windows' reads are joined by read name into fragments: one read spanning both sites (long
 *     reads), an RNA-seq read spliced across the introns between them (phASER's principle), the two mates of a pair,
 *     the parts of a split read. A source that samples a deep window keeps whole fragments by name hash at a power-of-two
 *     rate, so the fragments kept at the deeper site are kept at the other one too.
 *  2. Alleles: each fragment's at each site (alleleAt, the page's phasing rule): unknown when not covered, under base
 *     quality MIN_BQ, another allele, or when two of its reads disagree.
 *  3. The 2×2 table, fragments by allele at A × allele at B. Which cells say something depends on the sites:
 *     - two heterozygous sites: alt–alt and ref–ref say cis, alt–ref and ref–alt trans;
 *     - a site under HET_MIN (mosaic, subclonal): its reference reads come from cells without it too, so only the
 *       fragments carrying its allele count (the other site's allele on them: alt cis, ref trans);
 *     - a homozygous site (over HET_MAX) is on both haplotypes: the question does not arise.
 *  4. The reads' verdict: at least MIN_FRAGMENTS informative fragments, at most MAX_AGAINST of them against. The LOD
 *     weighs the same counts, log10 L(cis) / L(trans) with an error of EPS per allele read, so that 30 to 0 reads
 *     stronger than 3 to 0; it is shown, not used to decide (it grows with depth whatever the share against).
 *  5. Haplotags (HP / PS from WhatsHap, LongPhase, HiPhase…), when the reads carry them: each site's haplotype within its
 *     phase set (sitePhase), then same HP cis, different HP trans; two phase sets are not linked. The tool phased through
 *     the heterozygous sites between the two, so it links sites no read spans (two exons 100 kb apart in long reads).
 *  6. Without haplotags deciding, and no read holding both sites: the reads of the whole stretch between them phased
 *     by the page (phasing.ts, WhatsHap's principle in a greedy form), the two sites in cis or in trans when one phase
 *     block holds both, through the heterozygous sites between them (chainPair).
 *  7. The verdict: the reads' when they decide (the haplotags agreeing said so, disagreeing make it conflicting: a
 *     switch error of the tool, or reads joined that are not one molecule); the haplotags' when the reads cannot
 *     decide (too few fragments, or under MIN_CONFLICT disagreeing among themselves); otherwise the reads' (homozygous,
 *     conflicting fragments).
 */
import type { AlignedRead, PhaseBlock, PhaseResult, VariantSite } from '../components/sashimi/types';
import { alleleAt, HET_MAX, HET_MIN, type Allele } from './phasing';
import { callSites, inBlocks, readJunctions } from './collapse';
import { siteHaplotypesOf } from './alleles';
import { sitePhase, type SitePhase } from '../components/sashimi/sitePhase';

export const MIN_FRAGMENTS = 3, MAX_AGAINST = 0.2, MIN_BQ = 20, EPS = 0.01;
/** Informative fragments under which the reads disagreeing among themselves leave the verdict to the haplotags. */
export const MIN_CONFLICT = 10;
/** Reads asked for around each site; past it the source keeps whole fragments by name hash. */
export const PAIR_READS_CAP = 20_000;

// ---------------------------------------------------------------- the positions typed

/** An allele as typed: an SNV's base, or an insertion or deletion of a length. */
export type AlleleSpec = { kind: 'snv'; alt: string } | { kind: 'ins' | 'del'; length: number };
export interface PairTarget {
  chrom: string;
  /** 0-based: the base typed (1-based) minus one */
  pos: number;
  /** a VCF-style REF>ALT: an indel's position is its anchor base, the base before the event */
  vcf?: boolean;
  allele?: AlleleSpec;
  /** as typed, for the messages */
  text: string;
}

const sameChrom = (a: string, b: string) => a.replace(/^chr/i, '').toUpperCase() === b.replace(/^chr/i, '').toUpperCase();

/**
 * Positions typed as `chr:pos`, separated by `;` (or new lines): `pos` alone is on `chrom`; an allele may follow, as
 * `chr:pos:A>G`, `chr:pos A>G`, `chr:pos:G`, `chr:pos:+2` / `-3` (the viewer's labels), or gnomAD's `1-55051215-G-GA`.
 * A chromosome written with or without `chr` is the view's when it names the same one.
 */
export function parseTargets(text: string, chrom: string): { targets: PairTarget[]; errors: string[] } {
  const targets: PairTarget[] = [], errors: string[] = [];
  for (const raw of text.split(/[;\n]+/)) {
    const t = raw.trim();
    if (!t) continue;
    const m = /^(?:([^\s:]+?)[:\s-])?\s*([\d,]+)(?:[\s:-]+([ACGTN]+)\s*[>/-]\s*([ACGTN]+)|[\s:-]+([ACGTN])|[\s:]+([+-])(\d+))?\s*$/i.exec(t);
    if (!m) { errors.push(`${t}: not a position (chr:pos, chr:pos:A>G, chr:pos:+2)`); continue; }
    const pos1 = parseInt(m[2].replace(/,/g, ''), 10);
    if (!Number.isFinite(pos1) || pos1 < 1) { errors.push(`${t}: not a position`); continue; }
    const c = m[1] ? (sameChrom(m[1], chrom) ? chrom : m[1]) : chrom;
    const target: PairTarget = { chrom: c, pos: pos1 - 1, text: t };
    if (m[3]) {
      const ref = m[3].toUpperCase(), alt = m[4].toUpperCase();
      if (ref.length === 1 && alt.length === 1) target.allele = { kind: 'snv', alt };
      else if (ref[0] !== alt[0]) { errors.push(`${t}: an indel's REF and ALT start with the same (anchor) base, as in a VCF`); continue; }
      else if (ref.length === alt.length) { errors.push(`${t}: a multi-base substitution; give its bases one by one`); continue; }
      else { target.vcf = true; target.allele = { kind: alt.length > ref.length ? 'ins' : 'del', length: Math.abs(alt.length - ref.length) }; }
    } else if (m[5]) target.allele = { kind: 'snv', alt: m[5].toUpperCase() };
    else if (m[6]) target.allele = { kind: m[6] === '+' ? 'ins' : 'del', length: parseInt(m[7], 10) };
    targets.push(target);
  }
  return { targets, errors };
}

/** The window of reads a site needs: the base, or the deleted bases, with one base on each side. */
export function siteWindow(s: { pos: number; kind: VariantSite['kind']; length: number }): [number, number] {
  return [Math.max(0, s.pos - 1), s.pos + (s.kind === 'del' ? s.length : 1) + 1];
}
/** The window to read for a position typed, before its allele is known: room for the indels anchored at it. */
export function targetWindow(t: PairTarget): [number, number] {
  return [Math.max(0, t.pos - 1), t.pos + 2 + (t.allele?.kind === 'del' ? t.allele.length : 0)];
}

const matches = (s: VariantSite, a: AlleleSpec) => s.kind === a.kind && (a.kind === 'snv' ? s.alt === a.alt : s.length === a.length);

/**
 * The site at a position typed, from the reads over it (`known`: the variant scan's sites, preferred when one is
 * there). An SNV or an indel starting at the base typed; an indel anchored at it (VCF style) starts at the next base.
 * Without an allele, the most supported one, the others named.
 */
export function resolveTarget(t: PairTarget, reads: AlignedRead[], ref: string | null, refStart: number, known: VariantSite[] = []):
  { site: VariantSite | null; note?: string } {
  const [ws, we] = targetWindow(t);
  const called = callSites(reads, ws, we, ref, refStart, 1, 0, 0, 1);
  const at = (s: VariantSite) => t.allele ? (t.vcf && s.kind !== 'snv' ? s.pos === t.pos + 1 : s.pos === t.pos) && matches(s, t.allele)
    : s.pos === t.pos || (s.kind !== 'snv' && s.pos === t.pos + 1);
  const scan = known.filter(at), here = called.filter(at);
  const pick = (l: VariantSite[]) => [...l].sort((a, b) => b.alt_count - a.alt_count || (a.pos === t.pos ? -1 : 1))[0];
  const best = pick(scan.length ? scan : here);
  if (!best) {
    const depth = reads.filter(r => inBlocks(r, t.pos)).length;
    if (!depth) return { site: null, note: `${t.text}: no read over the position${reads.length ? '' : ` (is ${t.chrom} in the file, the position on it?)`}` };
    return { site: null, note: t.allele ? `${t.text}: no read carries that allele (${depth} reads over the position)` : `${t.text}: no alternate allele in the ${depth} reads over the position` };
  }
  // the counts of the reads read here, the scan's quality and haplotype evidence kept
  const mine = here.find(s => s.pos === best.pos && s.kind === best.kind && s.alt === best.alt);
  const site: VariantSite = { ...best, ...(mine ? { alt_count: mine.alt_count, depth: mine.depth, vaf: mine.vaf } : {}) };
  // the other alleles worth naming: a single noisy read's indel (long reads) is not one
  const others = (scan.length ? scan : here).filter(s => s !== best && s.alt_count >= Math.max(2, 0.05 * s.depth));
  return { site, note: !t.allele && others.length ? `${t.text}: ${others.length} other allele${others.length > 1 ? 's' : ''} there (${others.map(s => `${s.alt} ×${s.alt_count}`).join(', ')}); the most supported taken, add it to the position to choose (chr:pos:A>G)` : undefined };
}

// ---------------------------------------------------------------- fragments

/**
 * The reads of a window by name: a fragment's records there. A source decodes a read's mismatches over the window asked
 * for (and its margin) only, so a long read read at A says nothing of its bases at B: each site's allele is taken from
 * the copy read at that site, and the two copies joined by name.
 */
function byName(reads: AlignedRead[]): Map<string, AlignedRead[]> {
  const by = new Map<string, AlignedRead[]>();
  for (const r of reads) { const l = by.get(r.n); if (l) l.push(r); else by.set(r.n, [r]); }
  return by;
}

/** How a fragment links two sites: one read, one read spliced between them, two mates, the parts of a split read. */
export type LinkKind = 'read' | 'spliced' | 'pair' | 'split';

export interface PairFragment {
  name: string;
  a: Allele; b: Allele;
  link: LinkKind;
  /** the fragment's haplotag, when its reads carry one */
  hp?: number; ps?: number;
  /** its reads disagree at a site (overlapping mates, a split read's parts): left unknown there */
  conflict?: boolean;
}

function alleleOf(fr: AlignedRead[], s: VariantSite): { a: Allele; by: AlignedRead | null; conflict: boolean } {
  let a: Allele = -1, by: AlignedRead | null = null, conflict = false;
  for (const r of fr) {
    const x = alleleAt(r, s, MIN_BQ);
    if (x < 0) continue;
    if (a < 0) { a = x; by = r; } else if (a !== x) { conflict = true; a = -1; by = null; break; }
  }
  return { a, by, conflict };
}

function linkOf(ra: AlignedRead, rb: AlignedRead, pa: number, pb: number): LinkKind {
  if (ra === rb || (ra.s === rb.s && ra.e === rb.e && ra.f === rb.f)) {
    let spliced = false;
    const lo = Math.min(pa, pb), hi = Math.max(pa, pb);
    readJunctions(ra, (_, be, ns) => { if (be > lo && ns <= hi) spliced = true; });
    return spliced ? 'spliced' : 'read';
  }
  return (ra.f & 192) && (rb.f & 192) && (ra.f & 192) !== (rb.f & 192) ? 'pair' : 'split';
}

// ---------------------------------------------------------------- the verdict

export type SiteClass = 'het' | 'low' | 'hom';
export const siteClass = (vaf: number): SiteClass => vaf > HET_MAX ? 'hom' : vaf < HET_MIN ? 'low' : 'het';

export type PairCall = 'cis' | 'trans' | 'conflict' | 'few' | 'hom';
export interface PairHaplotags {
  call: 'cis' | 'trans' | 'unlinked' | 'hom' | 'undecided';
  a: SitePhase; b: SitePhase;
  note: string;
}
/** What the verdict rests on: the fragments holding both sites, the file's haplotags, or both agreeing. */
export type PairSource = 'reads' | 'haplotags' | 'both' | 'chain';

/** The two sites in a phase block of the reads between them: the block, and its sites between them on A's copy (0) or the other (1). */
export interface PairChain {
  block: PhaseBlock;
  /** heterozygous sites of the block strictly between the two, with the copy carrying their allele */
  between: { pos: number; copy: 0 | 1 }[];
  /** A's allele on haplotype 1 or 2 of the block */
  hapA: 1 | 2;
}
export interface PairResult {
  /** fragments by allele at A × at B: [alt–alt, alt–ref, ref–alt, ref–ref] */
  table: [number, number, number, number];
  /** fragments over both sites whose allele is unknown at one of them (quality, another allele, reads that disagree) */
  unknown: number;
  /** the cells counted, as indices into `table`, and the informative fragments for cis and for trans */
  counted: number[];
  cis: number; trans: number;
  /** the verdict, from `source` */
  call: PairCall;
  source: PairSource;
  /** the reads' own verdict and its words, whatever the haplotags say */
  readCall: PairCall;
  readNote: string;
  lod: number;
  note: string;
  classes: [SiteClass, SiteClass];
  fragments: PairFragment[];
  /** fragments by how they link the sites, among those with both alleles known */
  links: Record<LinkKind, number>;
  haplotags: PairHaplotags | null;
  /** each site's reads by haplotag lane, when the reads carry haplotags */
  tagReads: [SiteTagReads | null, SiteTagReads | null];
  /** the phase block of the reads between the sites that decides (source 'chain'), or why none does */
  chain?: PairChain;
  chainNote?: string;
}

const CELL_CIS = [true, false, false, true];
const LOD_STEP = Math.log10(((1 - EPS) ** 2 + EPS ** 2) / (2 * EPS * (1 - EPS)));

/**
 * The reads over a site by haplotag lane, with their allele: H1 and H2 of the site's phase set (the set holding the most
 * tagged reads over it), the reads tagged in another set, the untagged ones. What the panel draws when the haplotags
 * phase the pair.
 */
export type TagLane = 'h1' | 'h2' | 'other' | 'none';
export interface SiteTagReads {
  /** the site's phase set, null when its tagged reads carry no PS */
  ps: number | null;
  lanes: Record<TagLane, { alt: number; ref: number }>;
}

/** A site's haplotype within its phase set and its reads by lane, from the reads over it (every base quality, as the scan counts). */
export function siteTags(reads: AlignedRead[], s: VariantSite): { phase: SitePhase | null; tags: SiteTagReads | null } {
  const per = new Map<string, { ps: number | null; hp: number; depth: number; alt: number }>();
  const seen: { ps: number | null; hp: number | null; a: Allele }[] = [];
  for (const r of reads) {
    const x = alleleAt(r, s, 0);
    if (x < 0) continue;
    seen.push({ ps: r.ps ?? null, hp: r.hp ?? null, a: x });
    if (r.hp == null) continue;
    const k = `${r.ps ?? ''}\t${r.hp}`;
    let e = per.get(k);
    if (!e) per.set(k, e = { ps: r.ps ?? null, hp: r.hp, depth: 0, alt: 0 });
    e.depth++; if (x === 1) e.alt++;
  }
  const hap = siteHaplotypesOf([...per.values()]);
  if (!hap) return { phase: null, tags: null };
  const lanes: SiteTagReads['lanes'] = { h1: { alt: 0, ref: 0 }, h2: { alt: 0, ref: 0 }, other: { alt: 0, ref: 0 }, none: { alt: 0, ref: 0 } };
  for (const x of seen) {
    const lane: TagLane = x.hp == null ? 'none' : x.ps === hap.ps && (x.hp === 1 || x.hp === 2) ? (x.hp === 1 ? 'h1' : 'h2') : 'other';
    if (x.a === 1) lanes[lane].alt++; else lanes[lane].ref++;
  }
  return { phase: sitePhase({ ...s, hap }), tags: { ps: hap.ps, lanes } };
}

export function haplotagPair(a: SitePhase | null, b: SitePhase | null, [nA, nB] = ['A', 'B']): PairHaplotags | null {
  if (!a || !b) return null;
  const word = (p: SitePhase) => p.note.split(' · ')[0];
  if (a.call === 'both' || b.call === 'both') return { call: 'hom', a, b, note: `${a.call === 'both' ? nA : nB} is on both haplotypes (homozygous)` };
  if (a.ps !== b.ps) return { call: 'unlinked', a, b, note: `different phase sets (PS ${a.ps ?? 'none'} and ${b.ps ?? 'none'}): the phasing tool did not link the two sites` };
  const on = (p: SitePhase) => p.call === 'h1' ? 1 : p.call === 'h2' ? 2 : 0;
  if (!on(a) || !on(b)) return { call: 'undecided', a, b, note: `not decided by the haplotags: ${!on(a) ? `${nA} ${word(a)}` : ''}${!on(a) && !on(b) ? '; ' : ''}${!on(b) ? `${nB} ${word(b)}` : ''}` };
  return on(a) === on(b)
    ? { call: 'cis', a, b, note: `both on H${on(a)} of phase set ${a.ps ?? '(no PS)'}` }
    : { call: 'trans', a, b, note: `${nA} on H${on(a)}, ${nB} on H${on(b)} of phase set ${a.ps ?? '(no PS)'}` };
}

/**
 * The phase of sites A and B (same chromosome, either order) from the reads read around each (siteWindow), joined by
 * name; `names` call them in the notes.
 */
export function phasePair(A: VariantSite, B: VariantSite, readsA: AlignedRead[], readsB: AlignedRead[], names: [string, string] = ['A', 'B']): PairResult {
  const table: [number, number, number, number] = [0, 0, 0, 0];
  const links: Record<LinkKind, number> = { read: 0, spliced: 0, pair: 0, split: 0 };
  const fragments: PairFragment[] = [];
  let unknown = 0;
  const atB = byName(readsB);
  for (const [name, fa] of byName(readsA)) {
    const fb = atB.get(name);
    if (!fb) continue;
    const xa = alleleOf(fa, A), xb = alleleOf(fb, B);
    const overA = xa.by || xa.conflict || fa.some(r => inBlocks(r, A.pos));
    const overB = xb.by || xb.conflict || fb.some(r => inBlocks(r, B.pos));
    if (!overA || !overB) continue;
    if (xa.a < 0 || xb.a < 0) { unknown++; continue; }
    // one record holding both alleles links them most directly (a foldback read's primary and supplementary can both
    // cover the two sites, in an order that differs between sources)
    const one = fa.find(r => alleleAt(r, A, MIN_BQ) >= 0 && fb.some(q => q.s === r.s && q.e === r.e && q.f === r.f && alleleAt(q, B, MIN_BQ) >= 0));
    const link = one ? linkOf(one, one, A.pos, B.pos) : linkOf(xa.by!, xb.by!, A.pos, B.pos);
    links[link]++;
    table[(xa.a ? 0 : 2) + (xb.a ? 0 : 1)]++;
    const tagged = fa.find(r => r.hp != null) ?? fb.find(r => r.hp != null);
    fragments.push({ name, a: xa.a, b: xb.a, link, ...(tagged ? { hp: tagged.hp, ps: tagged.ps } : {}), ...(xa.conflict || xb.conflict ? { conflict: true } : {}) });
  }
  const classes: [SiteClass, SiteClass] = [siteClass(A.vaf), siteClass(B.vaf)];
  // the cells that say something (see the head of the file): [alt–alt, alt–ref, ref–alt, ref–ref]
  const counted = classes[0] === 'low' && classes[1] === 'low' ? [0, 1, 2]
    : classes[0] === 'low' ? [0, 1] : classes[1] === 'low' ? [0, 2] : [0, 1, 2, 3];
  const cis = counted.filter(i => CELL_CIS[i]).reduce((n, i) => n + table[i], 0);
  const trans = counted.filter(i => !CELL_CIS[i]).reduce((n, i) => n + table[i], 0);
  const n = cis + trans, lod = (cis - trans) * LOD_STEP;
  const pct = (x: number) => `${Math.round((100 * x) / n)} %`;
  let call: PairCall, note: string;
  if (classes.includes('hom')) {
    call = 'hom';
    note = `${classes[0] === 'hom' ? names[0] : names[1]} is homozygous (VAF over ${HET_MAX * 100} %): on both haplotypes, so with the other variant whatever its phase`;
  } else if (n < MIN_FRAGMENTS) {
    call = 'few';
    note = n ? `only ${n} informative fragment${n > 1 ? 's' : ''} (${MIN_FRAGMENTS} needed)` : 'no fragment holds both sites with their alleles known';
  } else if (trans <= MAX_AGAINST * n) {
    call = 'cis'; note = `${cis} of ${n} informative fragments carry them together (${pct(cis)})`;
  } else if (cis <= MAX_AGAINST * n) {
    call = 'trans'; note = `${trans} of ${n} informative fragments carry one without the other (${pct(trans)})`;
  } else {
    call = 'conflict';
    note = `${cis} fragments say cis, ${trans} trans: more than ${MAX_AGAINST * 100} % against either. A mosaic variant, a mis-mapped paralog or collapsed duplication, a third haplotype, or reads joined that are not one molecule?`;
  }
  if (classes.includes('low') && call !== 'hom') note += ` · ${classes[0] === 'low' && classes[1] === 'low' ? 'both sites' : classes[0] === 'low' ? names[0] : names[1]} under ${HET_MIN * 100} % VAF: only the fragments carrying ${classes[0] === 'low' && classes[1] === 'low' ? 'an' : 'its'} allele counted`;
  const ta = siteTags(readsA, A), tb = siteTags(readsB, B);
  const haplotags = haplotagPair(ta.phase, tb.phase, names);
  const readCall = call, readNote = note;
  let source: PairSource = 'reads';
  const tagged = haplotags && (haplotags.call === 'cis' || haplotags.call === 'trans') ? haplotags.call : null;
  if (tagged && (call === 'cis' || call === 'trans')) {
    if (tagged === call) { source = 'both'; note += ` · the haplotags agree (${haplotags!.note})`; }
    else {
      call = 'conflict';
      note = `the reads say ${readCall} (${readNote}), the haplotags ${tagged} (${haplotags!.note}): a switch error of the phasing tool, or reads joined that are not one molecule?`;
    }
  } else if (tagged && (call === 'few' || (call === 'conflict' && n < MIN_CONFLICT))) {
    note = call === 'conflict'
      ? `${haplotags!.note} · the ${n} reads holding both disagree (${cis} cis, ${trans} trans), too few to weigh against the phasing`
      : `${haplotags!.note}: the phasing tool linked them through the heterozygous sites between them (${readNote.replace(/^no fragment/, 'no read')})`;
    call = tagged; source = 'haplotags';
  }
  return { table, unknown, counted, cis, trans, call, source, readCall, readNote, lod, note, classes, fragments, links, haplotags, tagReads: [ta.tags, tb.tags] };
}

/** Reads a stretch is phased over (whole fragments past it, as every reads window keeps them). */
export const CHAIN_READS_CAP = 40_000;
/**
 * Long reads phase a stretch on its SNVs alone: their indel errors in homopolymers, called at a third of the reads, broke
 * the blocks every few kb on ONT HG002 (no pair over 20 kb phased); on SNVs alone, one block runs through CFTR.
 */
export const CHAIN_LONG_MIN_INDEL = 1_000_000;

/** Whether the reads between two sites should be phased for this pair: the reads holding both cannot decide, nor the haplotags. */
export function needsChain(r: PairResult): boolean {
  return (r.call === 'few' || (r.call === 'conflict' && r.source === 'reads' && r.cis + r.trans < MIN_CONFLICT)) && r.haplotags?.call !== 'cis' && r.haplotags?.call !== 'trans' && !r.classes.includes('hom');
}

/**
 * Sites A and B in the page's phasing of the reads of the stretch holding them (`phase`, its blocks): in cis when the
 * block holding both has their alleles on the same haplotype, in trans on the two; else why not (a site not phased, two
 * blocks the reads do not link).
 */
export function chainPair(phase: PhaseResult, A: VariantSite, B: VariantSite, [nA, nB] = ['A', 'B']): { call: 'cis' | 'trans' | null; note: string; chain?: PairChain } {
  const at = (x: VariantSite) => phase.sites.findIndex(s => s.pos === x.pos && s.kind === x.kind && s.alt === x.alt);
  const ia = at(A), ib = at(B);
  const where = (i: number, n: string) => {
    if (i < 0) return `${n} is not among the sites phased (fewer than 3 reads with it in the stretch read)`;
    const u = phase.unphased.find(x => x.site === i);
    if (u) return `${n} is not phased: ${u.reason === 'low' ? 'its allele fraction is outside the heterozygous range' : u.reason === 'conflict' ? 'the reads linking it disagree' : 'no read links it to another heterozygous site'}`;
    return null;
  };
  const why = where(ia, nA) ?? where(ib, nB);
  if (why) return { call: null, note: why };
  const blockOf = (i: number) => phase.blocks.find(b => b.sites.includes(i));
  const ba = blockOf(ia), bb = blockOf(ib);
  if (!ba || !bb) return { call: null, note: 'a site outside the phase blocks' };
  const span = (b: PhaseBlock) => `${(b.start + 1).toLocaleString('en-US')}–${b.end.toLocaleString('en-US')}`;
  if (ba !== bb) {
    const later = ba.start < bb.start ? bb : ba;
    return { call: null, note: `the reads phase ${nA} in the block ${span(ba)} and ${nB} in ${span(bb)}: no read links the two blocks${later.breakBefore === 'conflict' ? ' (the reads between them disagree)' : ''}` };
  }
  const hap = (i: number): 1 | 2 => { const k = ba.sites.indexOf(i); return ba.h1[k] === 'alt' ? 1 : 2; };
  const hapA = hap(ia), hapB = hap(ib);
  const lo = Math.min(A.pos, B.pos), hi = Math.max(A.pos, B.pos);
  const between = ba.sites.filter(i => phase.sites[i].pos > lo && phase.sites[i].pos < hi).map(i => ({ pos: phase.sites[i].pos, copy: (hap(i) === hapA ? 0 : 1) as 0 | 1 }));
  const call = hapA === hapB ? 'cis' : 'trans';
  return {
    call, chain: { block: ba, between, hapA },
    note: `phased through the reads between them: one phase block (${span(ba)}, ${ba.sites.length} heterozygous sites, ${between.length} between the two) holds both, ${call === 'cis' ? 'on the same haplotype' : 'on the two haplotypes'}; ${ba.support[0] + ba.support[1]} fragments assigned, ${ba.conflicting} disagreeing`,
  };
}
