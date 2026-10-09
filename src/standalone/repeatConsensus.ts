/**
 * Consensus of a set of reads at a repeat locus (the repeat inspector's "Consensus"), repeat-aware:
 *
 * - **Flanks**, base by base: each read's flank aligned on the reference's (alignFlank, the read view's alignment), the
 *   majority base or gap at each reference position among the reads that reach it, and the bases a majority inserts
 *   between two positions; outward, as long as half the reads reach the position.
 * - **Tract**, unit by unit: a base-level alignment of reads of different lengths inside a pure repeat is arbitrary
 *   (every CGG matches every other), so each read's tract is read as units in the motif's phase (unitRuns) and the
 *   units are aligned instead, every read on the read of median size (center-star): the majority unit, or a gap, at
 *   each of its columns, plus the units a majority inserts. Interruptions (AGG) keep their place; the consensus has
 *   about the median size. Reads that stop inside the repeat are aligned on their anchored end only.
 *
 * Every position carries its support: the share of the reads covering it that agree with the consensus.
 */
import { FLANK_SHOWN, alignFlank, unitRuns, type ReadRepeat, type RepeatLocus } from './repeatScan';

/** reads used at most (evenly by size past that): the unit alignment is quadratic in the tract's units */
export const CONSENSUS_MAX_READS = 300;
/** a position whose support is under this is shown in lowercase */
export const CONSENSUS_LOW = 0.6;

export interface ConsensusPart { seq: string; support: number[] }
export interface RepeatConsensus {
  /** reads used, and listed (some left out past CONSENSUS_MAX_READS, or for having no tract) */
  used: number; listed: number;
  flankL: ConsensusPart; flankR: ConsensusPart;
  /** the tract, as units (their support each) */
  units: { unit: string; cls: 'P' | 'B' | 'I' | 'o'; support: number }[];
  /** consensus size in units, and the reads' median size */
  size: number; medianSize: number;
  /** the tract built from reads stopping inside the repeat: a lower bound */
  lowerBound: boolean;
  /** "(CGG)10 AGG (CGG)9 AGG (CGG)131" */
  structure: string;
  /** interruptions' positions (1-based units) */
  interruptions: number[];
  /** mean support over the whole sequence */
  meanSupport: number;
}

/** Majority base per reference position of a flank, from the reads' alignments on it (`atEnd`: the flank ends at the tract). */
function flankConsensus(parts: string[], refFlank: string, atEnd: boolean): ConsensusPart {
  const n = refFlank.length;
  // per reference position (in tract-outward order): votes; per gap before it: inserted strings
  const votes = Array.from({ length: n }, () => new Map<string, number>());
  const ins = Array.from({ length: n + 1 }, () => new Map<string, number>());
  const cover = new Array<number>(n).fill(0);
  for (const p of parts) {
    const al = alignFlank(p, refFlank, atEnd);
    // walk from the tract outward, so that position 0 is next to the tract
    const ref = atEnd ? [...al.ref].reverse() : [...al.ref], read = atEnd ? [...al.read].reverse() : [...al.read];
    let j = 0, pending = '';
    for (let c = 0; c < ref.length && j < n; c++) {
      const r = ref[c], q = read[c];
      if (r === '-') { if (q !== ' ') pending += q; continue; }
      if (r === ' ') continue;
      if (q === ' ') break; // the read stops here
      if (pending) { const s = atEnd ? pending.split('').reverse().join('') : pending; ins[j].set(s, (ins[j].get(s) ?? 0) + 1); pending = ''; }
      votes[j].set(q, (votes[j].get(q) ?? 0) + 1);
      cover[j]++;
      j++;
    }
  }
  const outSeq: string[] = [], outSup: number[] = [];
  // a position counts while half the reads at least reach it: further out, the few reads still going are past the
  // amplicon's end (primer, adapter, barcode)
  const minCover = Math.max(Math.min(3, parts.length), parts.length / 2);
  for (let j = 0; j < n; j++) {
    if (cover[j] < minCover) break;
    // insertions a majority of the reads reaching here carry
    let insN = 0, best = '', bestN = 0;
    for (const [s, c] of ins[j]) { insN += c; if (c > bestN) { best = s; bestN = c; } }
    if (insN > cover[j] / 2) for (const ch of atEnd ? best.split('').reverse() : best.split('')) { outSeq.push(ch); outSup.push(bestN / cover[j]); }
    let base = '', top = 0;
    for (const [b, c] of votes[j]) if (c > top) { base = b; top = c; }
    if (base !== '-') { outSeq.push(base); outSup.push(top / cover[j]); }
  }
  return atEnd ? { seq: outSeq.reverse().join(''), support: outSup.reverse() } : { seq: outSeq.join(''), support: outSup };
}

type Unit = { u: string; cls: 'P' | 'B' | 'I' | 'o' };
const unitsOf = (tract: string, locus: RepeatLocus): Unit[] => {
  const out: Unit[] = [];
  for (const r of unitRuns(tract, locus)) if (r.cls !== 'x') for (let i = 0; i < r.n; i++) out.push({ u: r.unit, cls: r.cls });
  return out;
};

/**
 * A read's units aligned on the center's: for each center column the read's unit or null (a gap), and the units it
 * inserts before each column (index n: after the last). `free`: the read may stop early (its end costs nothing).
 */
function alignUnits(center: string[], read: string[], free: boolean): { col: (string | null | undefined)[]; ins: string[][] } {
  const n = center.length, m = read.length;
  const D: Uint32Array[] = Array.from({ length: n + 1 }, () => new Uint32Array(m + 1));
  for (let i = 0; i <= n; i++) D[i][0] = i;
  for (let j = 0; j <= m; j++) D[0][j] = j;
  for (let i = 1; i <= n; i++) {
    const c = center[i - 1], Di = D[i], Dp = D[i - 1];
    for (let j = 1; j <= m; j++) Di[j] = Math.min(Dp[j - 1] + (c === read[j - 1] ? 0 : 1), Dp[j] + 1, Di[j - 1] + 1);
  }
  // a read stopping early: its best end anywhere on the center (the rest of the center uncovered)
  let bi = n;
  if (free) { let best = Infinity; for (let i = 0; i <= n; i++) if (D[i][m] < best) { best = D[i][m]; bi = i; } }
  const col: (string | null | undefined)[] = new Array(n).fill(undefined);
  const ins: string[][] = Array.from({ length: n + 1 }, () => []);
  let i = bi, j = m;
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && D[i][j] === D[i - 1][j - 1] + (center[i - 1] === read[j - 1] ? 0 : 1)) { col[i - 1] = read[j - 1]; i--; j--; }
    else if (i > 0 && D[i][j] === D[i - 1][j] + 1) { col[i - 1] = null; i--; }
    else { ins[i].unshift(read[j - 1]); j--; }
  }
  // past a free end: not covered (undefined), unlike a gap (null)
  return { col, ins };
}

const reverseUnits = (u: Unit[]) => [...u].reverse();

/** The consensus of reads (those with a tract; at most CONSENSUS_MAX_READS, evenly by size). */
export function repeatConsensus(reads: ReadRepeat[], locus: RepeatLocus): RepeatConsensus | null {
  const withTract = reads.filter(r => !r.note && r.tract[1] > r.tract[0]);
  if (!withTract.length) return null;
  const sorted = [...withTract].sort((a, b) => a.units - b.units);
  const pick = sorted.length <= CONSENSUS_MAX_READS ? sorted : Array.from({ length: CONSENSUS_MAX_READS }, (_, i) => sorted[Math.round((i * (sorted.length - 1)) / (CONSENSUS_MAX_READS - 1))]);
  const full = pick.filter(r => !r.truncated);
  // the tract from the reads spanning it; from reads stopping inside only when there is none (a lower bound then)
  const tractReads = full.length ? full : pick;
  const lowerBound = !full.length;

  // flanks: the reads anchored on that side
  const leftParts = pick.filter(r => !r.truncated || r.from === 'left').map(r => r.seq.slice(Math.max(0, r.tract[0] - FLANK_SHOWN - 15), r.tract[0]));
  const rightParts = pick.filter(r => !r.truncated || r.from === 'right').map(r => r.seq.slice(r.tract[1], r.tract[1] + FLANK_SHOWN + 15));
  const flankL = leftParts.length ? flankConsensus(leftParts, locus.refL, true) : { seq: '', support: [] };
  const flankR = rightParts.length ? flankConsensus(rightParts, locus.refR, false) : { seq: '', support: [] };

  // tract: every read on the read of median size; reads anchored at the 3′ end only are aligned backwards
  const tracts = tractReads.map(r => ({ r, u: unitsOf(r.seq.slice(r.tract[0], r.tract[1]), locus) }));
  const sizes = tracts.map(t => t.u.length).sort((a, b) => a - b);
  const medianSize = sizes[sizes.length >> 1];
  const centerT = tracts.reduce((b, t) => (Math.abs(t.u.length - medianSize) < Math.abs(b.u.length - medianSize) ? t : b));
  const center = centerT.u;
  const cs = center.map(x => x.u);
  const n = cs.length;
  const votes = Array.from({ length: n }, () => new Map<string, number>());
  const cover = new Array<number>(n).fill(0);
  const insVotes = Array.from({ length: n + 1 }, () => ({ reads: 0, byLen: new Map<number, number>(), units: [] as Map<string, number>[] }));
  const clsOf = new Map<string, Unit['cls']>(center.map(x => [x.u, x.cls]));
  for (const t of tracts) {
    for (const x of t.u) if (!clsOf.has(x.u)) clsOf.set(x.u, x.cls);
    const backwards = t.r.truncated && t.r.from === 'right';
    const seq = backwards ? reverseUnits(t.u).map(x => x.u) : t.u.map(x => x.u);
    const al = alignUnits(backwards ? [...cs].reverse() : cs, seq, t.r.truncated);
    const col = backwards ? [...al.col].reverse() : al.col;
    const ins = backwards ? [...al.ins].reverse().map(a => [...a].reverse()) : al.ins;
    col.forEach((v, i) => { if (v === undefined) return; cover[i]++; const k = v ?? '-'; votes[i].set(k, (votes[i].get(k) ?? 0) + 1); });
    ins.forEach((a, i) => {
      if (!a.length) return;
      const slot = insVotes[i];
      slot.reads++; slot.byLen.set(a.length, (slot.byLen.get(a.length) ?? 0) + 1);
      a.forEach((u, k) => { (slot.units[k] ??= new Map()).set(u, (slot.units[k].get(u) ?? 0) + 1); });
    });
  }
  const units: RepeatConsensus['units'] = [];
  const coverAt = (i: number) => (i < n ? cover[i] : cover[n - 1] ?? 0);
  for (let i = 0; i <= n; i++) {
    const slot = insVotes[i], c = Math.max(1, coverAt(i));
    if (slot.reads > c / 2) {
      // the inserted length most reads carry, its units by majority
      let len = 0, lenN = 0; for (const [l, k] of slot.byLen) if (k > lenN) { len = l; lenN = k; }
      for (let k = 0; k < len; k++) {
        let best = '', bn = 0; for (const [u, v] of slot.units[k] ?? []) if (v > bn) { best = u; bn = v; }
        if (best) units.push({ unit: best, cls: clsOf.get(best) ?? 'o', support: bn / c });
      }
    }
    if (i === n || !cover[i]) continue;
    let best = '', top = 0;
    for (const [u, v] of votes[i]) if (v > top) { best = u; top = v; }
    if (best !== '-') units.push({ unit: best, cls: clsOf.get(best) ?? 'o', support: top / cover[i] });
  }

  // the structure: runs of the same unit
  const runs: [string, number][] = [];
  for (const u of units) { const last = runs[runs.length - 1]; if (last && last[0] === u.unit) last[1]++; else runs.push([u.unit, 1]); }
  const structure = runs.map(([u, k]) => (k > 1 ? `(${u})${k}` : u)).join(' ');
  const interruptions = units.map((u, i) => (u.cls === 'I' ? i + 1 : 0)).filter(Boolean);
  const all = [...flankL.support, ...units.flatMap(u => new Array(u.unit.length).fill(u.support)), ...flankR.support];
  return {
    used: pick.length, listed: reads.length, flankL, flankR, units, size: units.length, medianSize, lowerBound, structure, interruptions,
    meanSupport: all.length ? all.reduce((a, b) => a + b, 0) / all.length : 0,
  };
}

/** The consensus as one sequence, positions under CONSENSUS_LOW in lowercase. */
export function consensusSequence(c: RepeatConsensus): string {
  const low = (s: string, sup: number[]) => s.split('').map((b, i) => (sup[i] < CONSENSUS_LOW ? b.toLowerCase() : b)).join('');
  return low(c.flankL.seq, c.flankL.support) + c.units.map(u => (u.support < CONSENSUS_LOW ? u.unit.toLowerCase() : u.unit)).join('') + low(c.flankR.seq, c.flankR.support);
}

/** FASTA of a consensus, 60 bases a line, its header describing it. */
export function consensusFasta(c: RepeatConsensus, locus: RepeatLocus, sample: string, title: string): string {
  const id = `${locus.label.replace(/\s+/g, '_')}_consensus`;
  const head = `>${id} ${locus.chrom}:${locus.start + 1}-${locus.end} sample=${sample.replace(/\s+/g, '_')} reads=${c.used}/${c.listed} set="${title}" size=${c.size}${locus.motif}${c.lowerBound ? '(lower_bound)' : ''} median=${c.medianSize} flanks=${c.flankL.seq.length}+${c.flankR.seq.length}bp support=${Math.round(c.meanSupport * 100)}% lowercase=support<${Math.round(CONSENSUS_LOW * 100)}% structure="${c.structure}"`;
  const seq = consensusSequence(c);
  const lines: string[] = [];
  for (let i = 0; i < seq.length; i += 60) lines.push(seq.slice(i, i + 60));
  return `${head}\n${lines.join('\n')}\n`;
}
