// Unit checks of the CNV groups (cghGroups.ts), invented data only: arms and whole chromosomes, chrY's arms (PARs, Yp11.2
// TSPY array to the centromere, Yq12 left out), calls added by hand.
// Run: npx esbuild scripts/cgh-groups-check.ts --bundle --platform=node --outfile=/tmp/cgh-groups-check.cjs && node /tmp/cgh-groups-check.cjs
import { cnvGroups, type GroupOpts } from '../src/components/sashimi/cghGroups';
import type { CghCnv } from '../src/components/sashimi/cghIntervals';
let bad = 0;
const eq = (name: string, got: unknown, want: unknown) => { const ok = JSON.stringify(got) === JSON.stringify(want); if (!ok) bad++; console.log(ok ? 'ok  ' : 'FAIL', name, '|', JSON.stringify(got), ok ? '' : `| want ${JSON.stringify(want)}`); };
let n = 0;
const C = (chrom: string, start: number, end: number, kind: 'loss' | 'gain', probes: number, lr: number): CghCnv => ({ id: `c${n++}`, sample: 'P1', sampleKey: 'P1', chrom, start, end, cytoband: '', genes: [], sizeKb: (end - start + 1) / 1000, kind, typeText: '', probes, logRatio: lr, iscn: null, copies: null, copiesHi: null, amp: false, suppressed: false, classification: '', classRest: '', note: '', line: n });
const opts = (o: Partial<GroupOpts> = {}): GroupOpts => ({ build: 'GRCh37', sex: 'female', eligible: c => (c.probes ?? 0) >= 3, ...o });
const brief = (cs: CghCnv[], o?: Partial<GroupOpts>) => cnvGroups(cs, opts(o)).map(g => [g.label, g.members.length, g.arms.map(a => `${a.arm} ${Math.round(a.cover * 100)}`).join(',')]);

// 1. other chromosomes unchanged
eq('trisomy 8 (p + q)', brief([C('chr8', 200000, 43000000, 'gain', 800, 0.56), C('chr8', 47000000, 146200000, 'gain', 1700, 0.58)]), [['8 (whole chromosome)', 2, '8p 99,8q 100']]);
eq('4p arm', brief([C('chr4', 70000, 25000000, 'loss', 500, -0.95), C('chr4', 26500000, 48100000, 'loss', 450, -1.0)]), [['4p (whole arm)', 2, '4p 97']]);
eq('21 (acrocentric, q only)', brief([C('chr21', 14400000, 33000000, 'gain', 300, 0.57), C('chr21', 33200000, 48000000, 'gain', 280, 0.6)]), [['21 (whole chromosome)', 2, '21q 99']]);
eq('partial cover: no group', brief([C('chr8', 200000, 30000000, 'gain', 500, 0.58), C('chr8', 50000000, 100000000, 'gain', 900, 0.58)]), []);

// 2. chrY (GRCh37): before the Yp11.2 blind span, Yp read ~64 % here (no Yp group, the Yq gains alone a Yq group); two Yp gains over its probe-covered span (PAR1 end 2.65 Mb to ~8.9 Mb) and two Yq gains up to Yq12
const yp = [C('chrY', 2660000, 5600000, 'gain', 70, 1.4), C('chrY', 6150000, 8870000, 'gain', 100, 1.5)];
const yq = [C('chrY', 14550000, 26270000, 'gain', 320, 1.45), C('chrY', 26350000, 28600000, 'gain', 8, 1.3)];
eq('Y: Yp + Yq gains → whole Y', brief([...yp, ...yq], { sex: 'male' }), [['Y (whole chromosome)', 4, 'Yp 90,Yq 91']]);
eq('Y: Yp gains alone → Yp', brief(yp, { sex: 'male' }), [['Yp (whole arm)', 2, 'Yp 90']]);
eq('Y: Yp half → no Yp arm', brief([C('chrY', 2660000, 5600000, 'gain', 70, 1.4), C('chrY', 5700000, 5900000, 'gain', 5, 1.4)], { sex: 'male' }), []);
// GRCh38 (Yp ~2.79-9.0 Mb, cen 10.3-10.6, Yq12 from 26.6 Mb)
eq('Y GRCh38: whole Y', brief([C('chrY', 2790000, 5750000, 'gain', 70, 1.4), C('chrY', 6280000, 9000000, 'gain', 100, 1.5), C('chrY', 10700000, 24000000, 'gain', 320, 1.45), C('chrY', 24100000, 26500000, 'gain', 8, 1.3)], { sex: 'male', build: 'GRCh38' }),
  [['Y (whole chromosome)', 4, 'Yp 91,Yq 98']]);

// 3. a call added by hand: counted and weighted, the arms' cover computed again; ignored when of the other direction or no such group
const t8 = [C('chr8', 200000, 40000000, 'gain', 700, 0.2), C('chr8', 47000000, 146200000, 'gain', 1700, 0.2)];
const small = C('chr8', 40500000, 43000000, 'gain', 2, 0.6), loss = C('chr8', 44000000, 44500000, 'loss', 20, -0.5);
const gid = 'grp|P1|chr8|gain';
const [g0] = cnvGroups([...t8, small, loss], opts());
const [g1] = cnvGroups([...t8, small, loss], opts({ manual: new Map([[small.id, gid], [loss.id, gid]]) }));
eq('manual: before', [g0.members.length, g0.probes, g0.logRatio?.toFixed(4), g0.manual.size], [2, 2400, '0.2000', 0]);
eq('manual: after (2-probe gain under the minimums added, the loss refused)', [g1.members.map(c => c.id), g1.probes, g1.logRatio?.toFixed(4), [...g1.manual]], [[t8[0].id, small.id, t8[1].id], 2402, ((0.2 * 2400 + 0.6 * 2) / 2402).toFixed(4), [small.id]]);
eq('manual: 8p cover goes up', g1.arms[0].cover > g0.arms[0].cover, true);
eq('manual: unknown group ignored', cnvGroups([...t8, small], opts({ manual: new Map([[small.id, 'grp|P1|chr9|gain']]) }))[0].members.length, 2);
// a call added that holds an auto member: the auto member becomes nested, not counted twice
const big = C('chr8', 100000, 42000000, 'gain', 2, 0.25);
const [g2] = cnvGroups([...t8, big], opts({ manual: new Map([[big.id, gid]]) }));
eq('manual: holds a member → it is nested', [g2.nested.has(t8[0].id), g2.probes], [true, 2 + 1700]);

console.log(bad ? `${bad} FAILED` : 'all ok');
process.exit(bad ? 1 : 0);
