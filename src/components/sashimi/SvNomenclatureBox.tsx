/**
 * The draft HGVS and ISCN of a rebuilt structural variant, with a copy button and the warning that goes with them.
 */
import { useEffect, useState } from 'react';
import type { SvReconstruction } from './svModel';
import { nomenclatureOf, NOMENCLATURE_WARNING, type Cytoband } from './svNomenclature';

/** The bands of a chromosome, fetched once per chromosome (undefined while missing or unavailable). */
export function useCytobands(chrom: string, get?: (chrom: string) => Promise<Cytoband[]>): Cytoband[] | undefined {
  const [bands, setBands] = useState<{ chrom: string; list: Cytoband[] } | null>(null);
  useEffect(() => {
    let alive = true;
    if (!get) return;
    get(chrom).then(list => { if (alive) setBands({ chrom, list }); }).catch(() => { /* bands unavailable: ISCN without bands */ });
    return () => { alive = false; };
  }, [chrom, get]);
  return bands?.chrom === chrom && bands.list.length ? bands.list : undefined;
}

export default function SvNomenclatureBox({ rec, build, bands }: { rec: SvReconstruction; build: string; bands?: Cytoband[] }) {
  const n = nomenclatureOf(rec, build, c => (c === rec.chrs[0] ? bands : undefined));
  const [copied, setCopied] = useState(false);
  if (!n.hgvs.length && !n.iscnSeq.length) return null;
  const unc = n.uncertainty ? ` (breakpoints ${n.uncertainty})` : '';
  const lines = [...n.hgvs.map(h => `${h}${unc}`), ...(n.iscn ? [n.iscn] : []), ...n.iscnSeq];
  const copy = () => {
    const text = `${lines.join('\n')}\n${NOMENCLATURE_WARNING}`;
    navigator.clipboard?.writeText(text).then(() => { setCopied(true); setTimeout(() => setCopied(false), 1500); }).catch(() => { /* clipboard blocked */ });
  };
  const row = (label: string, text: string, title: string) => (
    <div className="mt-0.5" title={title}>
      <span className="text-[10px] font-semibold text-gray-500 uppercase tracking-wide mr-1">{label}</span>
      <span className="font-mono text-[11px] bg-accent-50 text-accent-900 rounded px-1.5 py-0.5 select-all break-all">{text}</span>
    </div>
  );
  return (
    <div className="mt-1">
      <div className="flex items-center gap-2">
        <span className="text-[10.5px] font-semibold text-gray-700">Nomenclature (draft)</span>
        <button onClick={copy} className="text-[10.5px] text-accent-700 hover:underline">{copied ? 'copied' : 'copy'}</button>
      </div>
      {n.hgvs.map(h => row('HGVS', h, 'HGVS allele form: the changes on the reference from pter to qter (HGVS recommendations, complex HGVS/ISCN). Best positions from the reads, not shifted by the 3′ rule.'))}
      {n.uncertainty && <div className="text-[10.5px] text-gray-500">breakpoints {n.uncertainty} around the positions given</div>}
      {n.iscn && row('ISCN', n.iscn, 'ISCN: the rebuilt chromosome, with the bands of its breakpoints (UCSC cytoBand)')}
      {n.iscnSeq.map(s => row('ISCN seq', s, 'ISCN sequence form: the rebuilt chromosome from pter to qter, :: at each junction'))}
      {!n.iscn && <div className="text-[10.5px] text-gray-500">bands unavailable: ISCN with bands not given</div>}
      <div className="text-[11px] font-semibold text-amber-800 mt-1">⚠ {NOMENCLATURE_WARNING}</div>
    </div>
  );
}
