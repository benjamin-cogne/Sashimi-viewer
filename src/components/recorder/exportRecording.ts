/**
 * A recording as one animated SVG (or an HTML page holding it, with play / pause and a time bar).
 *
 * The drawings copied during the recording are merged across frames: each element is written once, inline, shown
 * during the frames that hold it by a CSS animation (visibility, step-end); no <use> references, which the browser
 * would expand into a full copy of the drawing for every frame. A dialog's HTML blocks are drawn once each
 * (foreignObject), shown while they are in the dialog, and moved with the dialog's card. A camera above everything
 * zooms on areas (around the clicks, or drawn by hand) and back, with eased transforms: the drawing is vector, so a zoom
 * stays sharp. The pointer follows its path, smoothed; each click is a ring that widens and fades. Every animation has
 * the length of the whole recording, so they stay in step when the player pauses or seeks them, and when looping.
 */
import type { Recording, Layer } from './recorder';
import { renderCaptions, CAPTION_FADE_MS, type Caption, type Chapter } from './captions';

/** A zoom: the area shown (page coordinates) from `start` to `end` (recorded time, ms). */
export interface Zoom { start: number; end: number; rect: [number, number, number, number]; source: 'auto' | 'drawn' }
/** How recorded time plays: idle stretches cut, stretches cut by hand, stretches sped up or slowed down (recorded ms). */
export interface TimeEdits {
  /** idle stretches (nothing drawn changes, the pointer still) longer than this are cut down to it; 0 keeps them */
  maxIdleMs: number;
  cuts: { start: number; end: number }[];
  speeds: { start: number; end: number; factor: number }[];
}
export interface ExportOptions {
  time: TimeEdits;
  loop: boolean;
  /** text drawn in the recording replaced at export: sample names → "Sample 1"… */
  replace: { from: string; to: string }[];
  title: string;
  zooms: Zoom[];
  pointer: { show: boolean; smooth: boolean; scale: number; clicks: boolean };
  captions: Caption[];
  chapters: Chapter[];
}
/** The idle stretch kept by default, in ms. */
export const MAX_IDLE_DEFAULT_MS = 1000;
/** Auto-zoom: clicks closer than this in time (ms) share one zoom; the zoom starts ZOOM_IN_MS before the first and holds ZOOM_HOLD_MS after the last. */
export const ZOOM_GROUP_MS = 3000, ZOOM_IN_MS = 700, ZOOM_HOLD_MS = 1800, ZOOM_OUT_MS = 800;
/** Auto-zoom magnification by default. */
export const ZOOM_LEVEL_DEFAULT = 2;
/** Pause at the end before a loop starts again, in ms. */
const END_HOLD_MS = 1500;
const CLICK_MS = 600;
const EASE = 'cubic-bezier(.65,0,.35,1)';

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const pct = (t: number, d: number) => `${Math.max(0, Math.min(100, (t / d) * 100)).toFixed(3).replace(/\.?0+$/, '')}%`;
const n1 = (v: number) => String(Math.round(v * 10) / 10);

/** The idle stretches the pause shortening removes (recorded ms): what lies past the first `maxIdle` of each gap between events. */
export function idleCuts(rec: Recording, maxIdle: number): { start: number; end: number }[] {
  if (!(maxIdle > 0)) return [];
  const ev = [0, ...rec.frames.map(f => f.t), ...rec.pointer.map(p => p.t), ...rec.clicks.map(c => c.t), rec.duration].sort((a, b) => a - b);
  const out: { start: number; end: number }[] = [];
  for (let i = 1; i < ev.length; i++) if (ev[i] - ev[i - 1] > maxIdle) out.push({ start: ev[i - 1] + maxIdle, end: ev[i] });
  return out;
}

/**
 * Recorded time to played time: piecewise linear, with the rate 0 in a cut (idle or by hand) and 1 / factor in a
 * stretch played faster or slower (a cut wins over a speed).
 */
export function timeMap(rec: Recording, edits: TimeEdits): (t: number) => number {
  const cuts = [...idleCuts(rec, edits.maxIdleMs), ...edits.cuts.filter(c => c.end > c.start)];
  const speeds = edits.speeds.filter(s => s.end > s.start && s.factor > 0);
  if (!cuts.length && !speeds.length) return t => t;
  const cut = (t: number) => cuts.some(c => t >= c.start && t < c.end);
  // in overlapping speed stretches, the one set last applies
  const rate = (t: number) => { if (cut(t)) return 0; const sp = speeds.filter(s => t >= s.start && t < s.end); return sp.length ? 1 / sp[sp.length - 1].factor : 1; };
  const xs = [...new Set([0, rec.duration, ...cuts.flatMap(c => [c.start, c.end]), ...speeds.flatMap(s => [s.start, s.end])])].filter(x => x >= 0 && x <= rec.duration).sort((a, b) => a - b);
  const acc: number[] = [0], rates: number[] = [];
  for (let i = 0; i + 1 < xs.length; i++) { const r = rate((xs[i] + xs[i + 1]) / 2); rates.push(r); acc.push(acc[i] + r * (xs[i + 1] - xs[i])); }
  return t => {
    if (t <= 0) return t;
    if (t >= rec.duration) return acc[acc.length - 1] + (t - rec.duration);
    let lo = 0, hi = xs.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= t) lo = m; else hi = m; }
    return acc[lo] + rates[lo] * (t - xs[lo]);
  };
}
/** Played time back to recorded time (the map is monotone). */
export function invertTime(map: (t: number) => number, played: number, duration: number): number {
  let lo = 0, hi = duration;
  for (let k = 0; k < 40; k++) { const m = (lo + hi) / 2; if (map(m) < played) lo = m; else hi = m; }
  return hi;
}

/** The text between tags with the replacements made (attributes and tag names untouched). */
export function replaceText(s: string, rep: ExportOptions['replace']): string {
  if (!rep.length) return s;
  const list = rep.filter(r => r.from).map(r => ({ from: esc(r.from), to: esc(r.to) })).sort((a, b) => b.from.length - a.from.length);
  return s.replace(/>([^<]+)</g, (_m, txt: string) => `>${list.reduce((acc, r) => acc.split(r.from).join(r.to), txt)}<`);
}

/** The area the recording shows: every drawing and dialog card, clipped to what showed, with a margin. */
export function boundsOf(rec: Recording): [number, number, number, number] {
  // the whole page was recorded: the whole window
  if (rec.frames.some(f => f.layers[0]?.blocks.length || f.layers[0]?.fixed?.length)) return [0, 0, rec.width, rec.height];
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const add = (r: [number, number, number, number]) => { x0 = Math.min(x0, r[0]); y0 = Math.min(y0, r[1]); x1 = Math.max(x1, r[0] + r[2]); y1 = Math.max(y1, r[1] + r[3]); };
  const clamp = (r: [number, number, number, number]): [number, number, number, number] => { const a = Math.max(0, r[0]), b = Math.max(0, r[1]); return [a, b, Math.min(rec.width, r[0] + r[2]) - a, Math.min(rec.height, r[1] + r[3]) - b]; };
  for (const f of rec.frames) f.layers.forEach((l, li) => {
    for (const s of l.svgs) add(s.clip ?? clamp([s.x, s.y, s.w, s.h]));
    if (li > 0) for (const p of l.panels.slice(1, 2)) add(p.clip ?? clamp([p.x, p.y, p.w, p.h]));
  });
  if (!Number.isFinite(x0)) return [0, 0, rec.width, rec.height];
  const M = 8;
  x0 = Math.max(0, x0 - M); y0 = Math.max(0, y0 - M); x1 = Math.min(rec.width, x1 + M); y1 = Math.min(rec.height, y1 + M);
  return [Math.round(x0), Math.round(y0), Math.round(x1 - x0), Math.round(y1 - y0)];
}

/** A rectangle of the screen's proportions holding `r`, kept inside the bounds. */
function fitRect(r: [number, number, number, number], b: [number, number, number, number]): [number, number, number, number] {
  const ar = b[2] / b[3];
  let [x, y, w, h] = r;
  if (w / h > ar) { const nh = w / ar; y -= (nh - h) / 2; h = nh; } else { const nw = h * ar; x -= (nw - w) / 2; w = nw; }
  if (w >= b[2] || h >= b[3]) return [...b];
  x = Math.min(Math.max(x, b[0]), b[0] + b[2] - w); y = Math.min(Math.max(y, b[1]), b[1] + b[3] - h);
  return [x, y, w, h];
}

/** Zooms on the clicks inside the recorded area: clicks close in time share one, the area holding them all. */
export function autoZooms(rec: Recording, level = ZOOM_LEVEL_DEFAULT): Zoom[] {
  const b = boundsOf(rec);
  const inside = rec.clicks.filter(c => c.x >= b[0] && c.x <= b[0] + b[2] && c.y >= b[1] && c.y <= b[1] + b[3]);
  const groups: (typeof inside)[] = [];
  for (const c of inside) {
    const g = groups[groups.length - 1];
    if (g && c.t - g[g.length - 1].t < ZOOM_GROUP_MS) g.push(c); else groups.push([c]);
  }
  const w = b[2] / level, h = b[3] / level, m = 60;
  return groups.map(g => {
    const xs = g.map(c => c.x), ys = g.map(c => c.y);
    const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2;
    const rw = Math.max(w, Math.max(...xs) - Math.min(...xs) + 2 * m), rh = Math.max(h, Math.max(...ys) - Math.min(...ys) + 2 * m);
    return { start: g[0].t, end: g[g.length - 1].t + ZOOM_HOLD_MS, rect: fitRect([cx - rw / 2, cy - rh / 2, rw, rh], b), source: 'auto' as const };
  });
}

/** The pointer's path, in played time: still stretches kept still, moves smoothed (simplified, then through a Catmull-Rom curve). */
function pointerPath(rec: Recording, map: (t: number) => number, smooth: boolean): { t: number; x: number; y: number }[] {
  const pts = rec.pointer.map(p => ({ t: map(p.t), x: p.x, y: p.y }));
  if (!smooth || pts.length < 3) return pts;
  // runs of movement: a gap in the samples means the pointer stood still
  const runs: (typeof pts)[] = [];
  for (const p of pts) { const r = runs[runs.length - 1]; if (r && p.t - r[r.length - 1].t < 150) r.push(p); else runs.push([p]); }
  const out: typeof pts = [];
  for (const run of runs) {
    // Ramer–Douglas–Peucker on the positions (2.5 px)
    const keep = new Array(run.length).fill(false); keep[0] = keep[run.length - 1] = true;
    const rdp = (a: number, b: number) => {
      let best = -1, dmax = 2.5;
      const dx = run[b].x - run[a].x, dy = run[b].y - run[a].y, len = Math.hypot(dx, dy) || 1;
      for (let i = a + 1; i < b; i++) { const d = Math.abs(dy * run[i].x - dx * run[i].y + run[b].x * run[a].y - run[b].y * run[a].x) / len; if (d > dmax) { dmax = d; best = i; } }
      if (best >= 0) { keep[best] = true; rdp(a, best); rdp(best, b); }
    };
    if (run.length > 2) rdp(0, run.length - 1);
    const k = run.filter((_, i) => keep[i]);
    if (out.length) out.push({ ...out[out.length - 1], t: k[0].t - 1 }); // still until this run starts
    for (let i = 0; i < k.length - 1; i++) {
      const p0 = k[Math.max(0, i - 1)], p1 = k[i], p2 = k[i + 1], p3 = k[Math.min(k.length - 1, i + 2)];
      const steps = Math.max(1, Math.round((p2.t - p1.t) / 50));
      for (let s = 0; s < steps; s++) {
        const u = s / steps, u2 = u * u, u3 = u2 * u;
        const cr = (a: number, b: number, c: number, d: number) => 0.5 * (2 * b + (-a + c) * u + (2 * a - 5 * b + 4 * c - d) * u2 + (-a + 3 * b - 3 * c + d) * u3);
        out.push({ t: p1.t + (p2.t - p1.t) * u, x: cr(p0.x, p1.x, p2.x, p3.x), y: cr(p0.y, p1.y, p2.y, p3.y) });
      }
    }
    out.push(k[k.length - 1]);
  }
  return out;
}

export interface Exported { svg: string; css: string; hasHtml: boolean; duration: number; frames: number; bytes: number; bounds: [number, number, number, number]; chapters: { t: number; name: string }[] }

export function recordingToSvg(rec: Recording, opts: ExportOptions): Exported {
  return build(rec, opts, false) as Exported;
}

/**
 * The recording for the HTML player: the same scene, but nothing in it is a CSS animation. The player (playerHtml)
 * reads `timeline` and shows, hides and moves the elements itself, only when they change; the pointer, the clicks and
 * the camera are layers of the page above the scene, so moving them does not redraw it.
 */
export function recordingToPlayer(rec: Recording, opts: ExportOptions): PlayerExport {
  return build(rec, opts, true) as PlayerExport;
}

/** What the player runs (times in played ms, positions in the view's coordinates). */
export interface Timeline {
  W: number; H: number; D: number; loop: boolean;
  /** window classes: the stretches during which their elements show */
  wins: Record<string, [number, number][]>;
  /** move classes: the position of their elements from each time on */
  moves: Record<string, [number, number, number][]>;
  /** the camera: the area shown from each time on ([x, y, w, h] in view coordinates, null = all), eased between */
  cam: [number, [number, number, number, number] | null][];
  ptr: [number, number, number][]; ptrScale: number;
  clicks: [number, number, number][];
  fades: [string, number, number][];
}
export interface PlayerExport extends Omit<Exported, 'svg'> { scene: string; screen: string; timeline: Timeline }

function build(rec: Recording, opts: ExportOptions, player: boolean): Exported | PlayerExport {
  const map = timeMap(rec, opts.time);
  const end = map(rec.duration);
  const D = Math.max(1, end + (opts.loop ? END_HOLD_MS : 0));
  const B = boundsOf(rec);
  const [x0, y0, W, H] = B;
  const iter = opts.loop ? 'infinite' : '1';
  const tl: Timeline = { W: 0, H: 0, D, loop: opts.loop, wins: {}, moves: {}, cam: [], ptr: [], ptrScale: opts.pointer.scale, clicks: [], fades: [] };
  const css: string[] = [
    player ? '.fr{visibility:hidden}' : `.fr{visibility:hidden;animation:${D}ms step-end ${iter} forwards}`,
    player ? '' : `.a{animation-duration:${D}ms;animation-iteration-count:${iter};animation-fill-mode:forwards}`,
    // what a clip, mask, gradient or pattern holds is used, not drawn: it must stay visible whatever hides the group
    // around it (an invisible clip path clips everything away)
    'clipPath,clipPath *,mask *,pattern *,marker *,linearGradient *,radialGradient *,filter *{visibility:visible!important}',
  ];

  // frames shown for some time once the pauses are cut, with their window
  const shown = rec.frames.map((f, i) => ({ f, i, a: map(f.t), b: i + 1 < rec.frames.length ? map(rec.frames[i + 1].t) : Infinity })).filter(x => x.b > x.a);
  const order = new Map(shown.map((x, k) => [x.i, k]));
  // a class showing an element during a set of frames (runs of consecutive shown frames), one per distinct set
  const winClass = new Map<string, string>();
  const win = (frames: number[]): string => {
    const key = frames.join(',');
    let c = winClass.get(key);
    if (!c) {
      c = `w${winClass.size}`;
      let kf = '0%{visibility:hidden}', k = 0;
      const runs: [number, number][] = [];
      while (k < frames.length) {
        let e = k;
        while (e + 1 < frames.length && order.get(frames[e + 1]) === order.get(frames[e])! + 1) e++;
        const a = shown[order.get(frames[k])!].a, z = shown[order.get(frames[e])!].b;
        kf += `${pct(a, D)}{visibility:visible}${z === Infinity ? '' : `${pct(z, D)}{visibility:hidden}`}`;
        runs.push([Math.round(a), z === Infinity ? Math.round(D) + 1 : Math.round(z)]);
        k = e + 1;
      }
      if (player) tl.wins[c] = runs; else css.push(`@keyframes ${c}{${kf}}.${c}{animation-name:${c}}`);
      winClass.set(key, c);
    }
    return c;
  };
  /** A step-end animation through per-frame positions; null when the position never changes. */
  let nMoves = 0;
  const moveClass = (pos: { i: number; x: number; y: number }[]): string | null => {
    if (pos.every(q => q.x === pos[0].x && q.y === pos[0].y)) return null;
    const id = `m${nMoves++}`;
    if (player) { tl.moves[id] = pos.map(q => [Math.round(shown[order.get(q.i)!].a), q.x, q.y]); return id; }
    css.push(`@keyframes ${id}{${pos.map(q => `${pct(shown[order.get(q.i)!].a, D)}{transform:translate(${n1(q.x)}px,${n1(q.y)}px)}`).join('')}}.${id}{animation-name:${id};animation-timing-function:step-end}`);
    return id;
  };

  const clips = new Map<string, string>();
  const clipId = (c: [number, number, number, number]) => { const k = c.join(','); if (!clips.has(k)) clips.set(k, `k${clips.size}`); return clips.get(k)!; };
  const clipped = (c: [number, number, number, number] | null | undefined, inner: string) => (c ? `<g clip-path="url(#${clipId(c)})">${inner}</g>` : inner);
  const panel = (p: Layer['panels'][number]) => clipped(p.clip, `<rect x="${p.x}" y="${p.y}" width="${p.w}" height="${p.h}" rx="${p.radius}" fill="${p.fill}"${p.stroke ? ` stroke="${p.stroke}" stroke-width="${p.strokeWidth}"` : ''}/>`);

  /**
   * The drawings, written once: every SVG seen in one place (size, view box, style, clip) is one slot, moved to where
   * it is in each frame; inside, the elements of every frame are merged into one tree, an element with the same
   * content (or a group with the same tag) at the same rank kept once, shown during the frames that hold it. A
   * drawing that does not change is thus written once for the whole recording, and a change adds only what changed.
   */
  interface TreeNode { chunk: number; open?: string; frames: number[]; kids: TreeNode[][] }
  const insert = (alts: TreeNode[][], pos: number, id: number, frame: number) => {
    const c = rec.chunks[id];
    const list = (alts[pos] ??= []);
    let node = list.find(n => ('leaf' in c ? n.chunk === id : n.open === c.open));
    if (!node) { node = { chunk: id, open: 'leaf' in c ? undefined : c.open, frames: [], kids: [] }; list.push(node); }
    if (node.frames[node.frames.length - 1] !== frame) node.frames.push(frame);
    if (!('leaf' in c)) c.kids.forEach((k, j) => insert(node!.kids, j, k, frame));
  };
  const emit = (n: TreeNode, parentFrames: number): string => {
    const c = rec.chunks[n.chunk];
    const own = n.frames.length !== parentFrames;
    const inner = 'leaf' in c ? replaceText(c.leaf, opts.replace) : `${replaceText(n.open!, opts.replace)}${n.kids.map(list => list.map(k => emit(k, n.frames.length)).join('')).join('')}${c.close}`;
    return own ? `<g class="fr ${win(n.frames)}">${inner}</g>` : inner;
  };

  let nSlots = 0;
  // layer by layer: the page, then each dialog (its backdrop and card, its HTML, its drawings)
  const nLayers = Math.max(1, ...shown.map(x => x.f.layers.length));
  const body: string[] = [];
  let nBlocks = 0;
  for (let L = 0; L < nLayers; L++) {
    const here = shown.filter(x => x.f.layers[L]);
    const part: string[] = [];
    // the stacking levels of this layer (see zKey), lowest first: at each level its backgrounds, its HTML, its drawings
    const zs = new Set<string>(['']);
    for (const { f } of here) { const l = f.layers[L]; for (const it of [...l.panels, ...l.blocks, ...l.svgs]) zs.add(it.z ?? ''); }
    const levels = [...zs].sort();
    const zIs = (z: string) => (it: { z?: string }) => (it.z ?? '') === z;
    // the HTML blocks: each placement once, shown while a frame holds it; the page's (or the dialog card's) blocks
    // move with its origin, the fixed ones stay on the screen
    const blockSet = (pick: (l: Layer) => Layer['blocks'] | undefined): string[] => {
      const windows = new Map<string, { b: Layer['blocks'][number]; frames: number[] }>();
      for (const { f, i } of here) for (const b of pick(f.layers[L]) ?? []) {
        const k = `${b.id}|${b.x}|${b.y}|${b.w}|${b.h}|${b.clip?.join(',') ?? ''}`;
        if (!windows.has(k)) windows.set(k, { b, frames: [] });
        windows.get(k)!.frames.push(i);
      }
      const blocks: string[] = [];
      for (const { b, frames } of windows.values()) {
        nBlocks++;
        const h = rec.html[b.id];
        blocks.push(clipped(b.clip, `<foreignObject class="fr ${win(frames)}" x="${b.x}" y="${b.y}" width="${b.w}" height="${b.h}"><div xmlns="http://www.w3.org/1999/xhtml" style="${esc(h.style)};width:${b.w}px;height:${b.h}px;overflow:hidden">${replaceText(h.html, opts.replace)}</div></foreignObject>`));
      }
      return blocks;
    };
    let originMove: string | null | undefined;
    for (const z of levels) {
      // backgrounds: each distinct set once
      const panelSets = new Map<string, { ps: Layer['panels']; frames: number[] }>();
      for (const { f, i } of here) { const ps = f.layers[L].panels.filter(zIs(z)); if (!ps.length) continue; const k = JSON.stringify(ps); if (!panelSets.has(k)) panelSets.set(k, { ps, frames: [] }); panelSets.get(k)!.frames.push(i); }
      for (const { ps, frames } of panelSets.values()) part.push(`<g class="fr ${win(frames)}">${ps.map(panel).join('')}</g>`);
      const moving = blockSet(l => l.blocks.filter(zIs(z)));
      if (moving.length) {
        if (originMove === undefined) originMove = moveClass(here.map(({ f, i }) => ({ i, x: f.layers[L].origin[0], y: f.layers[L].origin[1] })));
        const o = here[0].f.layers[L].origin;
        part.push(originMove ? `<g class="a ${originMove}">${moving.join('\n')}</g>` : `<g transform="translate(${o[0]} ${o[1]})">${moving.join('\n')}</g>`);
      }
      // the drawings, slot by slot
      const slots = new Map<string, { s: Layer['svgs'][number]; frames: number[]; pos: { i: number; x: number; y: number }[]; kids: TreeNode[][] }>();
      for (const { f, i } of here) {
        const seen = new Map<string, number>();
        for (const sv of f.layers[L].svgs.filter(zIs(z))) {
          const base = `${sv.w}|${sv.h}|${sv.viewBox}|${sv.par}|${sv.style}|${sv.clip ? sv.clip.map((v, k) => v - (k === 0 ? sv.x : k === 1 ? sv.y : 0)).join(',') : ''}`;
          // two copies of one drawing in a frame are two slots
          const nth = seen.get(base) ?? 0; seen.set(base, nth + 1);
          const key = `${base}#${nth}`;
          let slot = slots.get(key);
          if (!slot) { slot = { s: sv, frames: [], pos: [], kids: [] }; slots.set(key, slot); }
          slot.frames.push(i); slot.pos.push({ i, x: sv.x, y: sv.y });
          sv.kids.forEach((k, j) => insert(slot!.kids, j, k, i));
        }
      }
      for (const { s: sv, frames, pos, kids } of slots.values()) {
        // the ids a drawing defines (clip paths, gradients, patterns) are its own: the page reuses one id while its
        // geometry changes (a clip that moves when a track is added), and the first in the file would serve them all
        const sfx = `_s${nSlots++}`;
        const inner = kids.map(list => list.map(k => emit(k, frames.length)).join('')).join('')
          .replace(/\sid="([^"]+)"/g, ` id="$1${sfx}"`)
          .replace(/url\(\s*(['"]?)#([^'")\s]+)\1\s*\)/g, `url(#$2${sfx})`)
          .replace(/(\s(?:xlink:)?href)="#([^"]+)"/g, `$1="#$2${sfx}"`);
        const el = `<svg width="${sv.w}" height="${sv.h}"${sv.viewBox ? ` viewBox="${sv.viewBox}"` : ''}${sv.par ? ` preserveAspectRatio="${sv.par}"` : ''} style="${esc(sv.style)}">${inner}</svg>`;
        // the clip moves with the drawing: it is given relative to the drawing's place
        const rel: [number, number, number, number] | null = sv.clip ? [sv.clip[0] - sv.x, sv.clip[1] - sv.y, sv.clip[2], sv.clip[3]] : null;
        const mv = moveClass(pos);
        const placed = mv ? `<g class="a ${mv}">${clipped(rel, el)}</g>` : `<g transform="translate(${sv.x} ${sv.y})">${clipped(rel, el)}</g>`;
        part.push(`<g class="fr ${win(frames)}">${placed}</g>`);
      }
    }
    // fixed and sticky elements of the page, lowest level first
    const fixedLevels = [...new Set(here.flatMap(({ f }) => (f.layers[L].fixed ?? []).map(b => b.z ?? '')))].sort();
    for (const z of fixedLevels) part.push(...blockSet(l => l.fixed?.filter(zIs(z))));
    body.push(`<g clip-path="url(#${clipId([0, 0, rec.width, rec.height])})">${part.join('\n')}</g>`);
  }

  // the pointer and the clicks
  let pointer = '';
  if (opts.pointer.show && rec.pointer.length && player) tl.ptr = pointerPath(rec, map, opts.pointer.smooth).map(p => [Math.round(p.t), Math.round(p.x - x0), Math.round(p.y - y0)]);
  else if (opts.pointer.show && rec.pointer.length) {
    const path = pointerPath(rec, map, opts.pointer.smooth);
    const kf = path.map(p => `${pct(p.t, D)}{transform:translate(${Math.round(p.x)}px,${Math.round(p.y)}px)}`);
    css.push(`@keyframes ptr{0%{transform:translate(${Math.round(path[0].x)}px,${Math.round(path[0].y)}px)}${kf.join('')}}#ptr{animation-name:ptr;animation-timing-function:linear}`);
    pointer = `<g id="ptr" class="a"><path transform="scale(${opts.pointer.scale})" d="M0 0 L0 16 L4.5 12.2 L7.6 19 L10.2 17.8 L7.2 11.2 L13 11.2 Z" fill="#111827" stroke="#ffffff" stroke-width="1.3" stroke-linejoin="round"/></g>`;
  }
  if (player && opts.pointer.clicks) tl.clicks = rec.clicks.map(c => [Math.round(map(c.t)), Math.round(c.x - x0), Math.round(c.y - y0)]);
  const clickEls = !opts.pointer.clicks || player ? '' : rec.clicks.map((c, i) => {
    const a = map(c.t), b = a + CLICK_MS;
    css.push(`@keyframes ck${i}{0%,${pct(a, D)}{opacity:0;transform:scale(.3)}${pct(a + 40, D)}{opacity:.9;transform:scale(.5)}${pct(b, D)},100%{opacity:0;transform:scale(1.6)}}#ck${i}{animation-name:ck${i};animation-timing-function:ease-out}`);
    return `<circle id="ck${i}" class="a" cx="${Math.round(c.x)}" cy="${Math.round(c.y)}" r="${14 * opts.pointer.scale}" fill="none" stroke="#0f766e" stroke-width="3" style="opacity:0;transform-box:fill-box;transform-origin:center"/>`;
  }).join('');

  // the camera: eased zooms on the chosen areas and back
  const zooms = [...opts.zooms].sort((a, b) => a.start - b.start);
  if (zooms.length) {
    const tf = (r: [number, number, number, number] | null) => {
      if (!r) return 'transform:translate(0px,0px) scale(1)';
      const s = W / r[2];
      return `transform:translate(${n1(-s * (r[0] - x0))}px,${n1(-s * (r[1] - y0))}px) scale(${s.toFixed(4)})`;
    };
    const kfs: { t: number; r: [number, number, number, number] | null }[] = [{ t: 0, r: null }];
    let last: [number, number, number, number] | null = null;
    zooms.forEach((z, i) => {
      const rect = fitRect(z.rect, B);
      const prevT = kfs[kfs.length - 1].t;
      const s = map(z.start), a = Math.max(prevT, s - ZOOM_IN_MS), b = Math.max(s, a + 200);
      kfs.push({ t: a, r: last }, { t: b, r: rect });
      const e = Math.max(map(z.end), b), next = zooms[i + 1];
      if (next && map(next.start) - ZOOM_IN_MS <= e + ZOOM_OUT_MS) { kfs.push({ t: e, r: rect }); last = rect; }
      else { kfs.push({ t: e, r: rect }, { t: e + ZOOM_OUT_MS, r: null }); last = null; }
    });
    if (player) tl.cam = kfs.map(k => [Math.round(k.t), k.r && [k.r[0] - x0, k.r[1] - y0, k.r[2], k.r[3]]]);
    else css.push(`@keyframes cam{${kfs.map(k => `${pct(k.t, D)}{${tf(k.r)};animation-timing-function:${EASE}}`).join('')}}#cam{animation-name:cam;transform-box:view-box;transform-origin:0 0}`);
  }

  const cap = renderCaptions(opts.captions.map(c => ({ ...c, text: replaceText(`>${c.text}<`, opts.replace).slice(1, -1) })), opts.chapters, map, D, { x0, y0, W, H }, rec.fontFamily, player ? tl.fades : undefined);
  css.push(...cap.css);
  const chapters = opts.chapters.map(c => ({ t: map(c.t), name: c.name })).sort((a, b) => a.t - b.t);

  const clipDefs = [...clips].map(([k, id]) => { const [x, y, w, h] = k.split(','); return `<clipPath id="${id}"><rect x="${x}" y="${y}" width="${w}" height="${h}"/></clipPath>`; }).join('');
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${esc(rec.fontFamily)}">
<title>${esc(opts.title)}</title>
<style>${rec.fonts?.length ? `/* Inter 5.3: Copyright 2016 The Inter Project Authors (https://github.com/rsms/inter). SIL Open Font License 1.1 (https://openfontlicense.org). */\n${rec.fonts.join('\n')}` : ''}
${css.join('\n')}</style>
<defs>${clipDefs}</defs>
<rect width="${W}" height="${H}" fill="${rec.background}"/>
<g id="cam" class="a"><g transform="translate(${-x0} ${-y0})">
${body.join('\n')}
${cap.content}
${clickEls}${pointer}
</g></g>
${player ? '' : cap.screen}
</svg>`;

  // the page's CSS rules the HTML blocks use
  let pageCss = '';
  if (nBlocks) {
    const classes = new Set<string>();
    for (const h of rec.html) for (const m of h.html.matchAll(/\sclass="([^"]*)"/g)) for (const c of m[1].split(/\s+/)) if (c) classes.add(c.replace(/&amp;/g, '&'));
    pageCss = rec.css.filter(r => r.classes.every(c => classes.has(c))).map(r => r.text).join('\n');
  }
  const meta = { css: pageCss, hasHtml: nBlocks > 0, duration: D, frames: shown.length, bounds: B, chapters };
  if (player) {
    tl.W = W; tl.H = H;
    const screen = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" font-family="${esc(rec.fontFamily)}">${cap.screen}</svg>`;
    return { ...meta, scene: svg, screen, timeline: tl, bytes: utf8Length(svg) + utf8Length(screen) + utf8Length(pageCss) + JSON.stringify(tl).length };
  }
  return { ...meta, svg, bytes: utf8Length(svg) + utf8Length(pageCss) };
}

/**
 * The HTML page of a recording: the scene, the camera, the pointer and the captions, run by a small player that only
 * touches what changes, when it changes (an element shown or hidden, a block moved), and moves the pointer, the click
 * rings and the camera as layers of their own, so the scene is not redrawn for them. Play / pause, a time bar with the
 * chapters marked, the chapter list; with "pause at chapters" the playback stops at each chapter, for a talk: → /
 * PageDown (a presentation clicker) goes on, or jumps to the next chapter while playing; ← / PageUp goes back to the
 * start of the chapter (the one before when at its start). Space toggles play / pause. `bare`: the drawing alone,
 * filling the window, still at the start until recPlayer.play() (the video export).
 */
export function playerHtml(out: PlayerExport, title: string, extra: { bare?: boolean; project?: string } = {}): string {
  const tl = out.timeline, D = Math.round(out.duration);
  const ch = extra.bare ? [] : out.chapters.filter(c => c.t >= 0 && c.t < D).map(c => ({ t: Math.round(c.t), name: c.name }));
  const json = (v: unknown) => JSON.stringify(v).replace(/</g, '\\u003c');
  const bar = extra.bare ? 0 : ch.length ? 118 : 80;
  return `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>${esc(title)}</title>
${out.css ? `<style>${out.css.replace(/<\/style/gi, '<\\/style')}</style>` : ''}
<style>
body{margin:0;background:${extra.bare ? '#fff' : '#f4f6f8'};font:13px 'Inter Viewer',Inter,system-ui,-apple-system,'Segoe UI',sans-serif;color:#1f2937${extra.bare ? ';overflow:hidden' : ''}}
.rec-wrap{max-width:100%;padding:${extra.bare ? '0' : '16px'};box-sizing:border-box}
.rec-stage{background:#fff;${extra.bare ? '' : 'border-radius:12px;box-shadow:0 1px 3px rgba(0,0,0,.12);'}overflow:hidden}
.rec-view{position:relative;overflow:hidden;margin:0 auto;aspect-ratio:${tl.W}/${tl.H};width:min(100%,calc((100vh - ${bar}px) * ${tl.W} / ${tl.H}))}
#camv{position:absolute;inset:0;transform-origin:0 0}
#camv>svg,.rec-screen{position:absolute;left:0;top:0;width:100%;height:100%;display:block}
.rec-screen{pointer-events:none}
#ptr{position:absolute;left:0;top:0;will-change:transform;pointer-events:none}
#ptr svg{display:block;overflow:visible}
.rec-click{position:absolute;left:0;top:0;box-sizing:border-box;border:3px solid #0f766e;border-radius:50%;opacity:0;pointer-events:none}
.rec-bar{display:flex;align-items:center;gap:10px;margin-top:10px}
.rec-bar button{width:34px;height:30px;border-radius:8px;border:1px solid #d3dbe5;background:#fff;cursor:pointer;font-size:14px}
.rec-track{position:relative;flex:1;display:flex;align-items:center}
.rec-track input{width:100%;accent-color:#0f766e;margin:0}
.rec-tick{position:absolute;top:-7px;width:2px;height:8px;background:#0f766e;border-radius:1px;pointer-events:none}
.rec-time{font-variant-numeric:tabular-nums;color:#5b6b7e;min-width:90px;text-align:right}
.rec-ch{display:flex;flex-wrap:wrap;align-items:center;gap:6px;margin-top:8px}
.rec-ch button{border:1px solid #d3dbe5;background:#fff;border-radius:999px;padding:3px 10px;font:inherit;font-size:12px;cursor:pointer;color:#334155}
.rec-ch button.on{background:#0f766e;border-color:#0f766e;color:#fff}
.rec-ch label{margin-left:auto;color:#5b6b7e;font-size:12px;display:flex;align-items:center;gap:4px}
</style></head><body><div class="rec-wrap">
<div class="rec-stage" id="stage"><div class="rec-view" id="view">
<div id="camv">${out.scene}<div id="clicks"></div><div id="ptr"></div></div>
<div class="rec-screen">${out.screen}</div>
</div></div>
${extra.bare ? '' : `<div class="rec-bar"><button id="play" title="Play / pause (space)">❚❚</button><div class="rec-track" id="track"><input id="seek" type="range" min="0" max="${D}" step="10" value="0"></div><span class="rec-time" id="time"></span></div>
${ch.length ? `<div class="rec-ch" id="chapters"><label title="Stop at each chapter; → or PageDown goes on"><input type="checkbox" id="pauseAt" checked> pause at chapters</label></div>` : ''}`}
</div>
<script type="application/json" id="rec-timeline">${json(tl)}</script>
<script>
(function(){
  var T=JSON.parse(document.getElementById('rec-timeline').textContent), D=T.D, CH=${json(ch)}, BARE=${extra.bare ? 'true' : 'false'};
  var view=document.getElementById('view'), camv=document.getElementById('camv'), ptr=document.getElementById('ptr'), clickBox=document.getElementById('clicks');
  var list=function(c){return Array.prototype.slice.call(document.getElementsByClassName(c));};
  // what the timeline drives, looked up once
  var wins=[], moves=[], fades=[], clicks=[];
  for(var c in T.wins) wins.push({els:list(c),runs:T.wins[c],on:null});
  for(var m in T.moves) moves.push({els:list(m),pts:T.moves[m],k:-2});
  T.fades.forEach(function(f){var el=document.getElementById(f[0]);if(el)fades.push({el:el,a:f[1],b:f[2],o:-1});});
  var NS='http://www.w3.org/2000/svg', sc=T.ptrScale||1;
  if(T.ptr.length){ptr.innerHTML='<svg width="'+(14*sc)+'" height="'+(20*sc)+'" viewBox="0 0 14 20"><path d="M0 0 L0 16 L4.5 12.2 L7.6 19 L10.2 17.8 L7.2 11.2 L13 11.2 Z" fill="#111827" stroke="#ffffff" stroke-width="1.3" stroke-linejoin="round"/></svg>';}
  T.clicks.forEach(function(k){var d=document.createElement('div');d.className='rec-click';clickBox.appendChild(d);clicks.push({el:d,t:k[0],x:k[1],y:k[2],on:false});});
  var ease=function(p){return p<.5?4*p*p*p:1-Math.pow(-2*p+2,3)/2;};
  var at=function(runs,t){for(var i=0;i<runs.length;i++)if(t>=runs[i][0]&&t<runs[i][1])return true;return false;};
  var camAt=function(t){
    var K=T.cam, full=[0,0,T.W,T.H];
    if(!K.length)return full;
    var r=function(x){return x||full;};
    if(t<=K[0][0])return r(K[0][1]);
    for(var i=0;i+1<K.length;i++){var a=K[i],b=K[i+1];if(t<b[0]){var p=b[0]>a[0]?ease((t-a[0])/(b[0]-a[0])):1,ra=r(a[1]),rb=r(b[1]);return [ra[0]+(rb[0]-ra[0])*p,ra[1]+(rb[1]-ra[1])*p,ra[2]+(rb[2]-ra[2])*p,ra[3]+(rb[3]-ra[3])*p];}}
    return r(K[K.length-1][1]);
  };
  var ptrAt=function(t){var P=T.ptr;if(!P.length)return null;if(t<=P[0][0])return P[0];var lo=0,hi=P.length-1;if(t>=P[hi][0])return P[hi];while(hi-lo>1){var md=(lo+hi)>>1;if(P[md][0]<=t)lo=md;else hi=md;}var a=P[lo],b=P[hi],p=(t-a[0])/((b[0]-a[0])||1);return [t,a[1]+(b[1]-a[1])*p,a[2]+(b[2]-a[2])*p];};
  var camNow=[0,0,T.W,T.H], lastCam='', still=0;
  var render=function(t){
    var i,j;
    for(i=0;i<wins.length;i++){var w=wins[i],on=at(w.runs,t);if(on!==w.on){w.on=on;for(j=0;j<w.els.length;j++)w.els[j].style.visibility=on?'visible':'hidden';}}
    for(i=0;i<moves.length;i++){var mv=moves[i],P=mv.pts,k=-1;for(j=0;j<P.length&&P[j][0]<=t;j++)k=j;if(k<0)k=0;if(k!==mv.k){mv.k=k;for(j=0;j<mv.els.length;j++)mv.els[j].style.transform='translate('+P[k][1]+'px,'+P[k][2]+'px)';}}
    for(i=0;i<fades.length;i++){var f=fades[i],F=${CAPTION_FADE_MS},o=t<f.a||t>f.b+F?0:t<f.a+F?(t-f.a)/F:t<=f.b?1:1-(t-f.b)/F;o=Math.round(o*100)/100;if(o!==f.o){f.o=o;f.el.style.opacity=o;}}
    var fx=view.clientWidth/T.W;
    // the camera: a layer of its own while it moves, drawn sharp again once it stops
    camNow=camAt(t);var s=T.W/camNow[2],key=camNow.map(function(v){return v.toFixed(2);}).join(',')+fx;
    if(key!==lastCam){lastCam=key;still=0;camv.style.willChange='transform';camv.style.transform='translate('+(-camNow[0]*s*fx)+'px,'+(-camNow[1]*s*fx)+'px) scale('+s+')';}
    else if(++still===20)camv.style.willChange='auto';
    var q=ptrAt(t);if(q)ptr.style.transform='translate('+(q[1]*fx)+'px,'+(q[2]*fx)+'px) scale('+fx+')';
    for(i=0;i<clicks.length;i++){var ck=clicks[i],p=(t-ck.t)/${CLICK_MS};if(p>=0&&p<=1){var sz=28*sc*fx,sc2=.5+1.1*p;ck.el.style.width=ck.el.style.height=sz+'px';ck.el.style.transform='translate('+(ck.x*fx-sz/2)+'px,'+(ck.y*fx-sz/2)+'px) scale('+sc2+')';ck.el.style.opacity=(p<.07?p/.07:1-p)*.9;ck.on=true;}else if(ck.on){ck.on=false;ck.el.style.opacity=0;}}
  };
  ptr.style.transformOrigin='0 0';
  // the clock
  var cur=0, paused=BARE, t0=0, last=0, dragging=false;
  var now=function(){return performance.now();};
  var time=function(){return paused?cur:Math.min(D,now()-t0);};
  var set=function(ms){cur=Math.max(0,Math.min(D,ms));t0=now()-cur;last=cur;render(cur);};
  var play=BARE?null:document.getElementById('play'), seek=document.getElementById('seek'), timeEl=document.getElementById('time');
  var setPaused=function(p){if(p){cur=time();paused=true;}else{if(cur>=D)cur=0;t0=now()-cur;paused=false;}if(play)play.textContent=paused?'▶':'❚❚';};
  var toggle=function(){setPaused(!paused);};
  var pauseAt=function(){var b=document.getElementById('pauseAt');return !!(b&&b.checked);};
  var chList=document.getElementById('chapters'), chBtns=[];
  CH.forEach(function(c,i){
    var b=document.createElement('button');b.textContent=(i+1)+'. '+c.name;b.onclick=function(){setPaused(pauseAt());set(c.t);};
    chList.insertBefore(b,chList.lastElementChild);chBtns.push(b);
    var k=document.createElement('div');k.className='rec-tick';k.style.left='calc('+(c.t/D*100)+'% - 1px)';k.title=c.name;document.getElementById('track').appendChild(k);
  });
  var chapterAt=function(t){var k=-1;CH.forEach(function(c,i){if(c.t<=t+1)k=i;});return k;};
  var next=function(){var t=time();if(paused){setPaused(false);return;}for(var i=0;i<CH.length;i++)if(CH[i].t>t+50){if(pauseAt())setPaused(true);set(CH[i].t);return;}};
  var prev=function(){var t=time(),k=chapterAt(t);if(k<0){setPaused(true);set(0);return;}if(t-CH[k].t<1000&&k>0)k--;setPaused(true);set(CH[k].t);};
  /** a point of the window (client coordinates) in the view's coordinates, under the current camera */
  var toView=function(x,y){var r=view.getBoundingClientRect(),u=(x-r.left)/r.width,v=(y-r.top)/r.height;return [camNow[0]+u*camNow[2],camNow[1]+v*camNow[3]];};
  window.recPlayer={pause:function(){setPaused(true);},play:function(){if(BARE){set(0);}setPaused(false);},seek:function(ms){set(ms);},time:time,paused:function(){return paused;},duration:D,chapters:CH,toView:toView};
  if(play)play.onclick=toggle;
  if(!BARE)document.addEventListener('keydown',function(e){
    if(e.code==='Space'){e.preventDefault();toggle();}
    else if(e.key==='ArrowRight'||e.key==='PageDown'){e.preventDefault();next();}
    else if(e.key==='ArrowLeft'||e.key==='PageUp'){e.preventDefault();prev();}
  });
  if(seek){seek.oninput=function(){dragging=true;set(Number(seek.value));};seek.onchange=function(){dragging=false;};}
  set(0);
  (function tick(){
    if(!paused){
      var t=now()-t0;
      if(t>=D){if(T.loop){t=t%D;t0=now()-t;}else{t=D;cur=D;setPaused(true);}}
      // stop at a chapter crossed while playing
      if(!dragging&&pauseAt()&&t>=last){for(var i=0;i<CH.length;i++){if(CH[i].t>0&&CH[i].t>last&&CH[i].t<=t){setPaused(true);t=CH[i].t;cur=t;break;}}}
      last=t;render(t);
    } else render(cur);
    if(seek&&!dragging)seek.value=String(Math.round(time()));
    if(timeEl){var f=function(ms){var s=Math.max(0,ms)/1000;return Math.floor(s/60)+':'+('0'+Math.floor(s%60)).slice(-2);};timeEl.textContent=f(time())+' / '+f(D);}
    var k=chapterAt(time());chBtns.forEach(function(b,i){b.className=i===k?'on':'';});
    requestAnimationFrame(tick);
  })();
})();
</script>${extra.project ? `\n<script type="application/json" id="sashimi-recording-project">${extra.project.replace(/</g, '\\u003c')}</script>` : ''}</body></html>`;
}

/** The size of a string in UTF-8, without copying it. */
export function utf8Length(s: string): number {
  let n = s.length;
  for (let i = 0; i < s.length; i++) { const c = s.charCodeAt(i); if (c > 0x7f) { n += c > 0x7ff ? (c >= 0xd800 && c <= 0xdbff ? (i++, 2) : 2) : 1; } }
  return n;
}

/** The id of the script that carries a recording and its edits inside an exported page, to edit it again. */
export const PROJECT_SCRIPT_ID = 'sashimi-recording-project';
