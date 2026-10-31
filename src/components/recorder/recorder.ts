/**
 * Live recording of what the viewer draws, as vectors.
 *
 * Every SVG on screen (the tracks, the cartoons, the SV browser; not the icons) is copied whenever the page changes,
 * at most RECORD_FPS times a second, with its place on the screen. Each drawing element is serialised once: a
 * MutationObserver marks what changed since the last copy, so an unchanged element is a reference to the copy already
 * kept, and a frame is a list of references. Full-screen dialogs are kept as their backdrop and their card (HTML is not
 * copied), the pointer as its positions and clicks. The export (exportRecording.ts) turns this into one animated SVG.
 */

/** Copies per second, at most: enough for the cartoons' movements, light on long recordings. */
export const RECORD_FPS = 12;
/** Pointer positions kept per second, at most (and only when it moved). */
export const POINTER_HZ = 20;
/** An element with more descendants than this is kept as a group of its children, so a change deep inside it only renews the part that changed. */
const LEAF_MAX_DESCENDANTS = 40;
/** SVGs this small are icons and logos: left out. */
const ICON_MAX_PX = 48;

export type Chunk = { leaf: string } | { open: string; close: string; kids: number[] };
export interface Placement { x: number; y: number; w: number; h: number; viewBox: string | null; par: string | null; style: string; clip: [number, number, number, number] | null; kids: number[] }
export interface Panel { x: number; y: number; w: number; h: number; fill: string; radius: number }
/** `order`: per layer (the page, then each dialog above it), how many panels then SVGs it holds, in drawing order. */
export interface Frame { t: number; panels: Panel[]; svgs: Placement[]; order: [number, number][] }
export interface Recording {
  width: number; height: number; background: string; fontFamily: string;
  frames: Frame[]; chunks: Chunk[];
  pointer: { t: number; x: number; y: number }[]; clicks: { t: number; x: number; y: number }[];
  duration: number;
}
export interface Recorder { stop(): Recording; elapsed(): number }

const SVG_NS_ATTR = / xmlns="http:\/\/www\.w3\.org\/2000\/svg"/g;
const skipRoot = (el: Element) => !!el.closest('[data-recorder="skip"]');

export function startRecording(root: Document = document): Recorder {
  const t0 = performance.now();
  const now = () => performance.now() - t0;
  const chunks: Chunk[] = [];
  const index = new Map<string, number>();
  const cache = new WeakMap<Element, number>();
  const boxCache = new WeakMap<Element, { clips: boolean; fixed: boolean }>();
  const ser = new XMLSerializer();
  const frames: Frame[] = [];
  const pointer: Recording['pointer'] = [];
  const clicks: Recording['clicks'] = [];
  let lastKey = '';

  const intern = (key: string, make: () => Chunk): number => {
    let id = index.get(key);
    if (id == null) { id = chunks.length; chunks.push(make()); index.set(key, id); }
    return id;
  };
  const skipped = (el: Element) => el.tagName.toLowerCase() === 'title' || el.getAttribute('data-export') === 'skip';
  /** The chunk of an element: the whole element when small, else its own tag around its children's chunks. */
  const chunkOf = (el: Element): number => {
    const hit = cache.get(el);
    if (hit != null) return hit;
    let id: number;
    const many = el.getElementsByTagName('*').length > LEAF_MAX_DESCENDANTS;
    const hasSkip = !!el.querySelector('title, [data-export="skip"]');
    if (!el.childElementCount || (!many && !hasSkip)) {
      const s = ser.serializeToString(el).replace(SVG_NS_ATTR, '');
      id = intern(`L${s}`, () => ({ leaf: s }));
    } else {
      const tag = el.tagName;
      const open = ser.serializeToString(el.cloneNode(false)).replace(SVG_NS_ATTR, '').replace(/\s*\/>$/, '>').replace(new RegExp(`></${tag}>$`), '>');
      const kids: number[] = [];
      for (const c of Array.from(el.children)) if (!skipped(c)) kids.push(chunkOf(c));
      id = intern(`C${open}|${kids.join(',')}`, () => ({ open, close: `</${tag}>`, kids }));
    }
    cache.set(el, id);
    return id;
  };
  /** The part of an element's box its scrolling ancestors leave visible. */
  const clipOf = (el: Element, r: DOMRect): [number, number, number, number] | null => {
    let x0 = r.left, y0 = r.top, x1 = r.right, y1 = r.bottom, clipped = false;
    for (let a = el.parentElement; a && a !== root.body; a = a.parentElement) {
      let box = boxCache.get(a);
      if (!box) { const cs = getComputedStyle(a); box = { clips: cs.overflowX !== 'visible' || cs.overflowY !== 'visible', fixed: cs.position === 'fixed' }; boxCache.set(a, box); }
      if (box.clips) {
        const b = a.getBoundingClientRect();
        if (b.left > x0 || b.top > y0 || b.right < x1 || b.bottom < y1) clipped = true;
        x0 = Math.max(x0, b.left); y0 = Math.max(y0, b.top); x1 = Math.min(x1, b.right); y1 = Math.min(y1, b.bottom);
      }
      // a fixed box (a dialog) is not clipped by what holds it
      if (box.fixed) break;
    }
    x0 = Math.max(x0, 0); y0 = Math.max(y0, 0); x1 = Math.min(x1, innerWidth); y1 = Math.min(y1, innerHeight);
    if (x1 <= x0 || y1 <= y0) return [0, 0, 0, 0];
    return clipped || r.left < 0 || r.top < 0 || r.right > innerWidth || r.bottom > innerHeight ? [x0, y0, x1 - x0, y1 - y0] : null;
  };
  const r1 = (v: number) => Math.round(v * 10) / 10;

  /** A copy of one SVG, or null when it is not drawn. */
  const place = (svg: SVGSVGElement): Placement | null => {
    const r = svg.getBoundingClientRect();
    if (r.width <= ICON_MAX_PX && r.height <= ICON_MAX_PX) return null;
    if (r.width < 1 || r.height < 1 || r.bottom < 0 || r.top > innerHeight || r.right < 0 || r.left > innerWidth) return null;
    const cs = getComputedStyle(svg);
    if (cs.visibility === 'hidden' || cs.display === 'none' || Number(cs.opacity) === 0) return null;
    const clip = clipOf(svg, r);
    if (clip && clip[2] === 0) return null;
    const kids: number[] = [];
    for (const c of Array.from(svg.children)) if (!skipped(c)) kids.push(chunkOf(c));
    return {
      x: r1(r.left), y: r1(r.top), w: r1(r.width), h: r1(r.height), viewBox: svg.getAttribute('viewBox'), par: svg.getAttribute('preserveAspectRatio'),
      style: `color:${cs.color};font-family:${cs.fontFamily.replace(/"/g, "'")};font-size:${cs.fontSize}`, clip: clip && clip.map(r1) as Placement['clip'], kids,
    };
  };

  const capture = () => {
    // full-screen dialogs, in page order: each above what comes before it
    const overlays = Array.from(root.querySelectorAll<HTMLElement>('.fixed.inset-0')).filter(ov => {
      if (skipRoot(ov)) return false;
      const cs = getComputedStyle(ov);
      return cs.display !== 'none' && cs.visibility !== 'hidden';
    });
    const layers: { panels: Panel[]; svgs: Placement[] }[] = [{ panels: [], svgs: [] }, ...overlays.map(ov => {
      const panels: Panel[] = [{ x: 0, y: 0, w: innerWidth, h: innerHeight, fill: getComputedStyle(ov).backgroundColor, radius: 0 }];
      const card = ov.firstElementChild as HTMLElement | null;
      if (card) {
        const r = card.getBoundingClientRect(), cc = getComputedStyle(card);
        if (r.width > 0 && cc.backgroundColor !== 'rgba(0, 0, 0, 0)') panels.push({ x: r1(r.left), y: r1(r.top), w: r1(r.width), h: r1(r.height), fill: cc.backgroundColor, radius: parseFloat(cc.borderTopLeftRadius) || 0 });
      }
      return { panels, svgs: [] as Placement[] };
    })];
    for (const svg of Array.from(root.querySelectorAll('svg'))) {
      if (svg.parentElement?.closest('svg') || skipRoot(svg)) continue;
      const p = place(svg);
      if (!p) continue;
      // the innermost dialog holding it
      let layer = 0;
      overlays.forEach((ov, i) => { if (ov.contains(svg)) layer = i + 1; });
      layers[layer].svgs.push(p);
    }
    const panels = layers.flatMap(l => l.panels);
    const svgs = layers.flatMap(l => l.svgs);
    const order = layers.map(l => [l.panels.length, l.svgs.length] as [number, number]);
    const key = JSON.stringify(layers);
    if (key === lastKey) return;
    lastKey = key;
    frames.push({ t: now(), panels, svgs, order });
  };

  // ---- scheduling: a copy after each burst of changes, at most RECORD_FPS a second ----
  let timer: number | null = null, lastCapture = -Infinity;
  const schedule = () => {
    if (timer != null) return;
    const wait = Math.max(0, lastCapture + 1000 / RECORD_FPS - now());
    timer = window.setTimeout(() => { timer = null; lastCapture = now(); requestAnimationFrame(capture); }, wait);
  };
  const dirty = (n: Node | null) => {
    for (let e: Element | null = n instanceof Element ? n : n?.parentElement ?? null; e; e = e.parentElement) {
      if (!cache.has(e)) { if (e.tagName.toLowerCase() === 'svg') break; continue; }
      cache.delete(e);
    }
  };
  const mo = new MutationObserver(records => {
    let any = false;
    for (const m of records) {
      const target = m.target instanceof Element ? m.target : m.target.parentElement;
      if (target && skipRoot(target)) continue;
      dirty(m.target);
      if (m.type === 'childList') m.removedNodes.forEach(n => { if (n instanceof Element) cache.delete(n); });
      any = true;
    }
    if (any) schedule();
  });
  mo.observe(root.body, { subtree: true, childList: true, attributes: true, characterData: true });
  const onScroll = () => schedule();
  let lastPtr = -Infinity;
  const onMove = (e: PointerEvent) => {
    const t = now();
    if (t - lastPtr < 1000 / POINTER_HZ) return;
    const p = pointer[pointer.length - 1];
    if (p && Math.abs(p.x - e.clientX) + Math.abs(p.y - e.clientY) < 2) return;
    lastPtr = t;
    pointer.push({ t, x: e.clientX, y: e.clientY });
  };
  const onDown = (e: PointerEvent) => {
    if (e.target instanceof Element && skipRoot(e.target)) return;
    const t = now();
    pointer.push({ t, x: e.clientX, y: e.clientY });
    clicks.push({ t, x: e.clientX, y: e.clientY });
  };
  addEventListener('scroll', onScroll, true);
  addEventListener('resize', onScroll);
  addEventListener('pointermove', onMove, { capture: true, passive: true });
  addEventListener('pointerdown', onDown, { capture: true, passive: true });
  capture();
  lastCapture = 0;

  return {
    elapsed: now,
    stop() {
      mo.disconnect();
      removeEventListener('scroll', onScroll, true);
      removeEventListener('resize', onScroll);
      removeEventListener('pointermove', onMove, true);
      removeEventListener('pointerdown', onDown, true);
      if (timer != null) clearTimeout(timer);
      capture();
      const bs = getComputedStyle(root.body);
      return {
        width: innerWidth, height: innerHeight, background: bs.backgroundColor === 'rgba(0, 0, 0, 0)' ? '#ffffff' : bs.backgroundColor, fontFamily: bs.fontFamily.replace(/"/g, "'"),
        frames, chunks, pointer, clicks, duration: now(),
      };
    },
  };
}
