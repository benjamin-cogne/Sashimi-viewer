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
import { renderCaptions, type Caption, type Chapter } from './captions';

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
  const map = timeMap(rec, opts.time);
  const end = map(rec.duration);
  const D = Math.max(1, end + (opts.loop ? END_HOLD_MS : 0));
  const B = boundsOf(rec);
  const [x0, y0, W, H] = B;
  const iter = opts.loop ? 'infinite' : '1';
  const css: string[] = [
    `.fr{visibility:hidden;animation:${D}ms step-end ${iter} forwards}`,
    `.a{animation-duration:${D}ms;animation-iteration-count:${iter};animation-fill-mode:forwards}`,
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
      while (k < frames.length) {
        let e = k;
        while (e + 1 < frames.length && order.get(frames[e + 1]) === order.get(frames[e])! + 1) e++;
        const a = shown[order.get(frames[k])!].a, z = shown[order.get(frames[e])!].b;
        kf += `${pct(a, D)}{visibility:visible}${z === Infinity ? '' : `${pct(z, D)}{visibility:hidden}`}`;
        k = e + 1;
      }
      css.push(`@keyframes ${c}{${kf}}.${c}{animation-name:${c}}`);
      winClass.set(key, c);
    }
    return c;
  };
  /** A step-end animation through per-frame positions; null when the position never changes. */
  let nMoves = 0;
  const moveClass = (pos: { i: number; x: number; y: number }[]): string | null => {
    if (pos.every(q => q.x === pos[0].x && q.y === pos[0].y)) return null;
    const id = `m${nMoves++}`;
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

  // layer by layer: the page, then each dialog (its backdrop and card, its HTML, its drawings)
  const nLayers = Math.max(1, ...shown.map(x => x.f.layers.length));
  const body: string[] = [];
  let nBlocks = 0;
  for (let L = 0; L < nLayers; L++) {
    const here = shown.filter(x => x.f.layers[L]);
    const part: string[] = [];
    // backgrounds: each distinct set once
    const panelSets = new Map<string, { ps: Layer['panels']; frames: number[] }>();
    for (const { f, i } of here) { const ps = f.layers[L].panels; if (!ps.length) continue; const k = JSON.stringify(ps); if (!panelSets.has(k)) panelSets.set(k, { ps, frames: [] }); panelSets.get(k)!.frames.push(i); }
    for (const { ps, frames } of panelSets.values()) part.push(`<g class="fr ${win(frames)}">${ps.map(panel).join('')}</g>`);
    // the HTML blocks of this layer: each placement once, shown while a frame holds it; the page's (or the dialog
    // card's) blocks move with its origin, the fixed ones stay on the screen
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
    const moving = blockSet(l => l.blocks);
    if (moving.length) {
      const mv = moveClass(here.map(({ f, i }) => ({ i, x: f.layers[L].origin[0], y: f.layers[L].origin[1] })));
      const o = here[0].f.layers[L].origin;
      part.push(mv ? `<g class="a ${mv}">${moving.join('\n')}</g>` : `<g transform="translate(${o[0]} ${o[1]})">${moving.join('\n')}</g>`);
    }
    // the drawings, slot by slot
    const slots = new Map<string, { s: Layer['svgs'][number]; frames: number[]; pos: { i: number; x: number; y: number }[]; kids: TreeNode[][] }>();
    for (const { f, i } of here) {
      const seen = new Map<string, number>();
      for (const sv of f.layers[L].svgs) {
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
      const inner = kids.map(list => list.map(k => emit(k, frames.length)).join('')).join('');
      const el = `<svg width="${sv.w}" height="${sv.h}"${sv.viewBox ? ` viewBox="${sv.viewBox}"` : ''}${sv.par ? ` preserveAspectRatio="${sv.par}"` : ''} style="${esc(sv.style)}">${inner}</svg>`;
      // the clip moves with the drawing: it is given relative to the drawing's place
      const rel: [number, number, number, number] | null = sv.clip ? [sv.clip[0] - sv.x, sv.clip[1] - sv.y, sv.clip[2], sv.clip[3]] : null;
      const mv = moveClass(pos);
      const placed = mv ? `<g class="a ${mv}">${clipped(rel, el)}</g>` : `<g transform="translate(${sv.x} ${sv.y})">${clipped(rel, el)}</g>`;
      part.push(`<g class="fr ${win(frames)}">${placed}</g>`);
    }
    part.push(...blockSet(l => l.fixed));
    body.push(`<g clip-path="url(#${clipId([0, 0, rec.width, rec.height])})">${part.join('\n')}</g>`);
  }

  // the pointer and the clicks
  let pointer = '';
  if (opts.pointer.show && rec.pointer.length) {
    const path = pointerPath(rec, map, opts.pointer.smooth);
    const kf = path.map(p => `${pct(p.t, D)}{transform:translate(${Math.round(p.x)}px,${Math.round(p.y)}px)}`);
    css.push(`@keyframes ptr{0%{transform:translate(${Math.round(path[0].x)}px,${Math.round(path[0].y)}px)}${kf.join('')}}#ptr{animation-name:ptr;animation-timing-function:linear}`);
    pointer = `<g id="ptr" class="a"><path transform="scale(${opts.pointer.scale})" d="M0 0 L0 16 L4.5 12.2 L7.6 19 L10.2 17.8 L7.2 11.2 L13 11.2 Z" fill="#111827" stroke="#ffffff" stroke-width="1.3" stroke-linejoin="round"/></g>`;
  }
  const clickEls = !opts.pointer.clicks ? '' : rec.clicks.map((c, i) => {
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
    css.push(`@keyframes cam{${kfs.map(k => `${pct(k.t, D)}{${tf(k.r)};animation-timing-function:${EASE}}`).join('')}}#cam{animation-name:cam;transform-box:view-box;transform-origin:0 0}`);
  }

  const cap = renderCaptions(opts.captions.map(c => ({ ...c, text: replaceText(`>${c.text}<`, opts.replace).slice(1, -1) })), opts.chapters, map, D, { x0, y0, W, H }, rec.fontFamily);
  css.push(...cap.css);
  const chapters = opts.chapters.map(c => ({ t: map(c.t), name: c.name })).sort((a, b) => a.t - b.t);

  const clipDefs = [...clips].map(([k, id]) => { const [x, y, w, h] = k.split(','); return `<clipPath id="${id}"><rect x="${x}" y="${y}" width="${w}" height="${h}"/></clipPath>`; }).join('');
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${esc(rec.fontFamily)}">
<title>${esc(opts.title)}</title>
<style>${css.join('\n')}</style>
<defs>${clipDefs}</defs>
<rect width="${W}" height="${H}" fill="${rec.background}"/>
<g id="cam" class="a"><g transform="translate(${-x0} ${-y0})">
${body.join('\n')}
${cap.content}
${clickEls}${pointer}
</g></g>
${cap.screen}
</svg>`;

  // the page's CSS rules the HTML blocks use
  let pageCss = '';
  if (nBlocks) {
    const classes = new Set<string>();
    for (const h of rec.html) for (const m of h.html.matchAll(/\sclass="([^"]*)"/g)) for (const c of m[1].split(/\s+/)) if (c) classes.add(c.replace(/&amp;/g, '&'));
    pageCss = rec.css.filter(r => r.classes.every(c => classes.has(c))).map(r => r.text).join('\n');
  }
  return { svg, css: pageCss, hasHtml: nBlocks > 0, duration: D, frames: shown.length, bytes: utf8Length(svg) + utf8Length(pageCss), bounds: B, chapters };
}

/**
 * An HTML page around the SVG: play / pause, a time bar with the chapters marked, the chapter list. With "pause at
 * chapters" the playback stops at each chapter, for a talk: → / PageDown (a presentation clicker) goes on, or jumps to
 * the next chapter while playing; ← / PageUp goes back to the start of the chapter (the one before when at its start).
 * Space toggles play / pause.
 */
export function svgToHtml(out: Exported, title: string, extra: { bare?: boolean; project?: string } = {}): string {
  const D = Math.round(out.duration);
  if (extra.bare) {
    // for the video export: the drawing alone, filling the window, paused at the start until told to play
    return `<!doctype html><html><head><meta charset="utf-8">${out.css ? `<style>${out.css.replace(/<\/style/gi, '<\\/style')}</style>` : ''}
<style>html,body{margin:0;height:100%;background:#fff;overflow:hidden}body>svg{display:block;width:100vw;height:100vh}</style></head><body>${out.svg}
<script>(function(){var svg=document.querySelector('body>svg');var a=function(){return svg.getAnimations({subtree:true});};a().forEach(function(x){x.pause();x.currentTime=0;});
window.recPlayer={play:function(){a().forEach(function(x){x.currentTime=0;x.play();});},duration:${D}};})();</script></body></html>`;
  }
  const ch = out.chapters.filter(c => c.t >= 0 && c.t < D);
  const chJson = JSON.stringify(ch.map(c => ({ t: Math.round(c.t), name: c.name }))).replace(/</g, '\\u003c');
  return `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>${esc(title)}</title>
${out.css ? `<style>${out.css.replace(/<\/style/gi, '<\\/style')}</style>` : ''}
<style>
body{margin:0;background:#f4f6f8;font:13px Inter,system-ui,-apple-system,'Segoe UI',sans-serif;color:#1f2937}
.rec-wrap{max-width:100%;padding:16px;box-sizing:border-box}
.rec-stage{background:#fff;border-radius:12px;box-shadow:0 1px 3px rgba(0,0,0,.12);overflow:hidden}
.rec-stage>svg{display:block;width:100%;height:auto;max-height:calc(100vh - ${ch.length ? 118 : 80}px)}
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
<div class="rec-stage" id="stage">${out.svg}</div>
<div class="rec-bar"><button id="play" title="Play / pause (space)">❚❚</button><div class="rec-track" id="track"><input id="seek" type="range" min="0" max="${D}" step="10" value="0"></div><span class="rec-time" id="time"></span></div>
${ch.length ? `<div class="rec-ch" id="chapters"><label title="Stop at each chapter; → or PageDown goes on"><input type="checkbox" id="pauseAt" checked> pause at chapters</label></div>` : ''}
</div>
<script>
(function(){
  var D=${D}, CH=${chJson}, svg=document.querySelector('#stage>svg'), play=document.getElementById('play'), seek=document.getElementById('seek'), time=document.getElementById('time');
  var anims=function(){return svg.getAnimations({subtree:true});}, paused=false, dragging=false, last=0;
  var fmt=function(ms){var s=Math.max(0,ms)/1000;return Math.floor(s/60)+':'+('0'+Math.floor(s%60)).slice(-2);};
  var cur=function(){var a=anims()[0];return a?(Number(a.currentTime)||0)%D:0;};
  var set=function(ms){anims().forEach(function(a){a.currentTime=ms;});last=ms;};
  var setPaused=function(p){paused=p;anims().forEach(function(a){p?a.pause():a.play();});play.textContent=p?'▶':'❚❚';};
  var toggle=function(){setPaused(!paused);};
  var pauseAt=function(){var b=document.getElementById('pauseAt');return !!(b&&b.checked);};
  var chList=document.getElementById('chapters'), chBtns=[];
  CH.forEach(function(c,i){
    var b=document.createElement('button');b.textContent=(i+1)+'. '+c.name;b.onclick=function(){setPaused(pauseAt());set(c.t);};
    chList.insertBefore(b,chList.lastElementChild);chBtns.push(b);
    var k=document.createElement('div');k.className='rec-tick';k.style.left='calc('+(c.t/D*100)+'% - 1px)';k.title=c.name;document.getElementById('track').appendChild(k);
  });
  var chapterAt=function(t){var k=-1;CH.forEach(function(c,i){if(c.t<=t+1)k=i;});return k;};
  var next=function(){var t=cur();if(paused){setPaused(false);return;}for(var i=0;i<CH.length;i++)if(CH[i].t>t+50){if(pauseAt())setPaused(true);set(CH[i].t);return;}};
  var prev=function(){var t=cur(),k=chapterAt(t);if(k<0){setPaused(true);set(0);return;}if(t-CH[k].t<1000&&k>0)k--;setPaused(true);set(CH[k].t);};
  window.recPlayer={pause:function(){setPaused(true);},play:function(){setPaused(false);},seek:function(ms){set(ms);},time:cur,paused:function(){return paused;},duration:D,chapters:CH};
  play.onclick=toggle;
  document.addEventListener('keydown',function(e){
    if(e.code==='Space'){e.preventDefault();toggle();}
    else if(e.key==='ArrowRight'||e.key==='PageDown'){e.preventDefault();next();}
    else if(e.key==='ArrowLeft'||e.key==='PageUp'){e.preventDefault();prev();}
  });
  seek.oninput=function(){dragging=true;set(Number(seek.value));};
  seek.onchange=function(){dragging=false;};
  (function tick(){
    var t=cur();
    // stop at a chapter crossed while playing (not on a loop back to the start)
    if(!paused&&!dragging&&pauseAt()&&t>=last){for(var i=0;i<CH.length;i++){if(CH[i].t>0&&CH[i].t>last&&CH[i].t<=t){setPaused(true);set(CH[i].t);t=CH[i].t;break;}}}
    last=t;
    if(!dragging)seek.value=String(Math.round(t));
    time.textContent=fmt(t)+' / '+fmt(D);
    var k=chapterAt(t);chBtns.forEach(function(b,i){b.className=i===k?'on':'';});
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
