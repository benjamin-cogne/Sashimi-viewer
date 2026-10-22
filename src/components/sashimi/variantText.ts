/**
 * A variant site's coordinates as text to copy, in the forms other tools take:
 *  - the position, chr:pos (a deletion's first and last bases, an insertion's two neighbours);
 *  - gnomAD / VCF style, chrom-pos-ref-alt with the anchor base before an indel (gnomAD, VarSome, Franklin take it);
 *  - HGVS genomic, chr:g.… on the chromosome's own name;
 *  - HGVS coding on the shown transcript, NM_…:c.…, bases on the transcript's strand.
 * Indels are placed where the reads put them, not shifted 3′ as HGVS asks: a repeat may move them.
 */
import type { VariantSite } from './types';
import { cdnaPosition, type TxModel } from './geometry';

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const revComp = (s: string) => { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] ?? 'N'; return o; };

export interface CoordinateText { label: string; text: string; title: string }

export function variantCoordinates(s: VariantSite, chrom: string, ref: { start: number; seq: string } | null, insSeq: string | undefined, tx: TxModel | null): CoordinateText[] {
  const p = s.pos + 1;    // 1-based
  const base = (pos0: number) => { const i = pos0 - (ref?.start ?? 0); return ref && i >= 0 && i < ref.seq.length ? ref.seq[i].toUpperCase() : null; };
  const bases = (from0: number, to0: number) => { let o = ''; for (let i = from0; i < to0; i++) { const b = base(i); if (!b) return null; o += b; } return o; };
  const out: CoordinateText[] = [];
  const shift = ' (placed where the reads put it, not shifted 3′)';
  // the position
  out.push({ label: 'position', title: 'Position', text: s.kind === 'snv' ? `${chrom}:${p}` : s.kind === 'del' ? `${chrom}:${p}${s.length > 1 ? `-${p + s.length - 1}` : ''}` : `${chrom}:${p - 1}-${p}` });
  // gnomAD / VCF: the chromosome without "chr" (MT for chrM)
  const c = chrom.replace(/^chr/i, '').replace(/^M$/i, 'MT');
  const refBase = s.kind === 'snv' ? (s.ref && s.ref !== '?' ? s.ref : base(s.pos)) : null;
  if (s.kind === 'snv' && refBase) out.push({ label: 'vcf', title: 'gnomAD / VCF style (chrom-pos-ref-alt)', text: `${c}-${p}-${refBase}-${s.alt}` });
  else if (s.kind === 'del') {
    const r = bases(s.pos - 1, s.pos + s.length);
    if (r) out.push({ label: 'vcf', title: `gnomAD / VCF style, the base before the deletion as anchor${shift}`, text: `${c}-${p - 1}-${r}-${r[0]}` });
  } else if (s.kind === 'ins') {
    const a = base(s.pos - 1);
    if (a && insSeq) out.push({ label: 'vcf', title: `gnomAD / VCF style, the base before the insertion as anchor, the inserted bases most reads carry${shift}`, text: `${c}-${p - 1}-${a}-${a}${insSeq}` });
  }
  // HGVS genomic
  const g = s.kind === 'snv' ? (refBase ? `${chrom}:g.${p}${refBase}>${s.alt}` : null)
    : s.kind === 'del' ? `${chrom}:g.${p}${s.length > 1 ? `_${p + s.length - 1}` : ''}del`
    : insSeq ? `${chrom}:g.${p - 1}_${p}ins${insSeq}` : `${chrom}:g.${p - 1}_${p}ins(${s.length})`;
  if (g) out.push({ label: 'g', title: `HGVS genomic, on the chromosome's name${s.kind === 'snv' ? '' : shift}`, text: g });
  // HGVS coding, on the transcript shown, when the site lies in it
  if (tx && s.pos >= tx.start - 5000 && s.pos < tx.end + 5000) {
    const plus = tx.strand > 0, cp = (pos0: number) => cdnaPosition(pos0, tx).bare;
    const pre = `${tx.transcriptId}:${tx.cdsStart == null ? 'n' : 'c'}.`;
    let t: string | null = null;
    if (s.kind === 'snv' && refBase) t = `${pre}${cp(s.pos)}${plus ? refBase : COMP[refBase] ?? 'N'}>${plus ? s.alt : COMP[s.alt] ?? 'N'}`;
    else if (s.kind === 'del') {
      const first = plus ? s.pos : s.pos + s.length - 1, last = plus ? s.pos + s.length - 1 : s.pos;
      t = `${pre}${cp(first)}${s.length > 1 ? `_${cp(last)}` : ''}del`;
    } else if (s.kind === 'ins') {
      const left = plus ? s.pos - 1 : s.pos, right = plus ? s.pos : s.pos - 1;
      t = `${pre}${cp(left)}_${cp(right)}ins${insSeq ? (plus ? insSeq : revComp(insSeq)) : `(${s.length})`}`;
    }
    if (t) out.push({ label: 'c', title: `HGVS on ${tx.transcriptId} (${tx.geneName}), bases on its strand${s.kind === 'snv' ? '' : shift}`, text: t });
  }
  return out;
}
