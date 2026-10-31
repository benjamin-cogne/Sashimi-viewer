/**
 * Live recording of what the viewer draws, as vectors.
 *
 * Every SVG on screen (the tracks, the cartoons, the SV browser; not the icons) is copied whenever the page changes,
 * at most RECORD_FPS times a second, with its place on the screen. Each drawing element is serialised once: a
 * MutationObserver marks what changed since the last copy, so an unchanged element is a reference to the copy already
 * kept, and a frame is a list of references. The page around them is kept whole too: its HTML (toolbars, buttons,
 * text) cut into blocks, each block copied once while it does not change, so a ticking clock renews a few words and not
 * the whole page; the blocks move with the page as it scrolls (fixed and sticky ones stay on the screen). Full-screen
 * dialogs are layers above it: their backdrop, their card, and their HTML cut the same way. The page's CSS rules those
 * blocks use go with the recording. The pointer is kept as its positions and clicks. The export (exportRecording.ts) turns this into one animated SVG.
 */

/** Copies per second, at most: enough for the cartoons' movements, light on long recordings. */
export const RECORD_FPS = 12;
/** Pointer positions kept per second, at most (and only when it moved). */
export const POINTER_HZ = 20;
/** An element with more descendants than this is kept as a group of its children, so a change deep inside it only renews the part that changed. */
const LEAF_MAX_DESCENDANTS = 40;
/** An HTML element of a dialog with more elements than this (or holding a drawing) is cut into its children. */
const HTML_BLOCK_MAX = 60;
/** SVGs this small are icons and logos: left out (inside a dialog, they go with the HTML around them). */
const ICON_MAX_PX = 48;

export type Chunk = { leaf: string } | { open: string; close: string; kids: number[] };
type Clip = [number, number, number, number] | null;
export interface Placement { x: number; y: number; w: number; h: number; viewBox: string | null; par: string | null; style: string; clip: Clip; kids: number[] }
export interface Panel { x: number; y: number; w: number; h: number; fill: string; radius: number; stroke?: string; strokeWidth?: number; clip?: Clip }
/** A piece of a dialog's HTML: `id` into Recording.html, drawn in a box of its size; position and clip relative to the dialog's card. */
export interface Block { x: number; y: number; w: number; h: number; id: number; clip: Clip }
/**
 * What is drawn at one level: the page (layer 0), then each dialog above it. A dialog's blocks are placed relative to
 * its card (`origin`), so scrolling the dialog moves them together instead of copying them again.
 */
export interface Layer {
  panels: Panel[]; blocks: Block[]; svgs: Placement[]; origin: [number, number];
  /** blocks that stay put on the screen (fixed or sticky elements of the page), in screen coordinates */
  fixed?: Block[];
}
export interface Frame { t: number; layers: Layer[] }
export interface CssRule { classes: string[]; text: string }
export interface Recording {
  width: number; height: number; background: string; fontFamily: string;
  frames: Frame[]; chunks: Chunk[];
  /** the HTML blocks: the wrapper's inherited style and the element */
  html: { style: string; html: string }[];
  /** the page's CSS rules, with the classes each needs (kept at export when the blocks use them all) */
  css: CssRule[];
  pointer: { t: number; x: number; y: number }[]; clicks: { t: number; x: number; y: number }[];
  duration: number;
}
export interface Recorder { stop(): Recording; elapsed(): number }

const SVG_NS_ATTR = / xmlns="http:\/\/www\.w3\.org\/2000\/svg"/g;
const skipRoot = (el: Element) => !!el.closest('[data-recorder="skip"]');
const TRANSPARENT = 'rgba(0, 0, 0, 0)';
const INHERITED = ['color', 'font-family', 'font-size', 'font-weight', 'font-style', 'line-height', 'letter-spacing', 'text-align', 'white-space', 'text-transform'];

/** The page's CSS rules, each with the class names its selector needs. */
function pageCss(doc: Document): CssRule[] {
  const out: CssRule[] = [];
  const classesOf = (sel: string) => [...new Set([...sel.matchAll(/\.((?:\\.|[\w-])+)/g)].map(m => m[1].replace(/\\(.)/g, '$1')))];
  const walk = (rules: CSSRuleList, wrap: (t: string) => string) => {
    for (const r of Array.from(rules)) {
      if (r instanceof CSSStyleRule) {
        // a rule for several selectors is needed when any of them is
        const alts = r.selectorText.split(',').map(classesOf);
        out.push({ classes: alts.reduce((best, a) => (a.length < best.length ? a : best), alts[0] ?? []), text: wrap(r.cssText) });
      } else if (r instanceof CSSMediaRule) walk(r.cssRules, t => wrap(`@media ${r.conditionText}{${t}}`));
      else if (typeof CSSSupportsRule !== 'undefined' && r instanceof CSSSupportsRule) walk(r.cssRules, t => wrap(`@supports ${r.conditionText}{${t}}`));
      else out.push({ classes: [], text: wrap(r.cssText) });
    }
  };
  for (const sheet of Array.from(doc.styleSheets)) { try { walk(sheet.cssRules, t => t); } catch { /* a sheet from elsewhere: not readable */ } }
  return out;
}

export function startRecording(root: Document = document): Recorder {
  const t0 = performance.now();
  const now = () => performance.now() - t0;
  const chunks: Chunk[] = [];
  const index = new Map<string, number>();
  const cache = new WeakMap<Element, number>();
  const html: Recording['html'] = [];
  const htmlIndex = new Map<string, number>();
  const htmlCache = new WeakMap<Element, { id: number; w: number; h: number }>();
  const boxCache = new WeakMap<Element, { clips: boolean; fixed: boolean }>();
  const ser = new XMLSerializer();
  const frames: Frame[] = [];
  const pointer: Recording['pointer'] = [];
  const clicks: Recording['clicks'] = [];
  let lastKey = '';
  const r1 = (v: number) => Math.round(v * 10) / 10;
  const bigSvg = (s: Element) => { const r = s.getBoundingClientRect(); return !(r.width <= ICON_MAX_PX && r.height <= ICON_MAX_PX); };

  // ---- SVG drawings ----
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
  /** The part of an element's box its scrolling ancestors (up to a fixed one) and the window leave visible. */
  const clipOf = (el: Element, r: DOMRect): Clip => {
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
    return clipped || r.left < 0 || r.top < 0 || r.right > innerWidth || r.bottom > innerHeight ? [r1(x0), r1(y0), r1(x1 - x0), r1(y1 - y0)] : null;
  };
  /** The part of an element left visible by the scrolling boxes between it and its dialog's card (not the window). */
  const innerClip = (el: Element, r: DOMRect): Clip => {
    let x0 = r.left, y0 = r.top, x1 = r.right, y1 = r.bottom, clipped = false;
    for (let a = el.parentElement; a && !(a.parentElement?.matches('.fixed.inset-0')); a = a.parentElement) {
      let box = boxCache.get(a);
      if (!box) { const cs = getComputedStyle(a); box = { clips: cs.overflowX !== 'visible' || cs.overflowY !== 'visible', fixed: cs.position === 'fixed' }; boxCache.set(a, box); }
      if (box.fixed) break;
      if (!box.clips) continue;
      const b = a.getBoundingClientRect();
      if (b.left > x0 || b.top > y0 || b.right < x1 || b.bottom < y1) clipped = true;
      x0 = Math.max(x0, b.left); y0 = Math.max(y0, b.top); x1 = Math.min(x1, b.right); y1 = Math.min(y1, b.bottom);
    }
    if (x1 <= x0 || y1 <= y0) return [0, 0, 0, 0];
    return clipped ? [r1(x0), r1(y0), r1(x1 - x0), r1(y1 - y0)] : null;
  };
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
      style: `color:${cs.color};font-family:${cs.fontFamily.replace(/"/g, "'")};font-size:${cs.fontSize}`, clip, kids,
    };
  };

  // ---- the HTML of dialogs ----
  const htmlCount = (el: Element) => {
    let n = el.getElementsByTagName('*').length;
    for (const s of Array.from(el.getElementsByTagName('svg'))) if (!s.parentElement?.closest('svg')) n -= s.getElementsByTagName('*').length;
    return n;
  };
  const directText = (el: Element) => Array.from(el.childNodes).some(n => n.nodeType === 3 && n.textContent!.trim());
  const splitCache = new WeakMap<Element, boolean>();
  const split = (el: Element) => {
    let v = splitCache.get(el);
    if (v == null) {
      v = el.childElementCount > 0 && !directText(el) && (htmlCount(el) > HTML_BLOCK_MAX || Array.from(el.getElementsByTagName('svg')).some(s => !s.parentElement?.closest('svg') && bigSvg(s)));
      splitCache.set(el, v);
    }
    return v;
  };
  const SKIP_TAGS = new Set(['SCRIPT', 'STYLE', 'NOSCRIPT', 'TEMPLATE', 'LINK', 'META']);
  /** One block: the element as it shows now (form states written into it, drawings left to the SVG layer). */
  const blockId = (el: HTMLElement, r: DOMRect): number => {
    const hit = htmlCache.get(el);
    // a block keeps its copy while it neither changes nor resizes (a layout shift around it)
    if (hit && hit.w === r1(r.width) && hit.h === r1(r.height)) return hit.id;
    const clone = el.cloneNode(true) as HTMLElement;
    const live = [el, ...Array.from(el.querySelectorAll('input, select, textarea'))], copy = [clone, ...Array.from(clone.querySelectorAll('input, select, textarea'))];
    live.forEach((n, i) => {
      const c = copy[i];
      if (n instanceof HTMLInputElement && c instanceof HTMLInputElement) {
        if (n.type === 'checkbox' || n.type === 'radio') { if (n.checked) c.setAttribute('checked', ''); else c.removeAttribute('checked'); } else c.setAttribute('value', n.value);
      } else if (n instanceof HTMLSelectElement && c instanceof HTMLSelectElement) {
        Array.from(c.options).forEach((o, k) => { if (k === n.selectedIndex) o.setAttribute('selected', ''); else o.removeAttribute('selected'); });
      } else if (n instanceof HTMLTextAreaElement && c instanceof HTMLTextAreaElement) c.textContent = n.value;
    });
    // the drawings inside are drawn by the SVG layer: kept here as empty boxes of their size
    for (const s of Array.from(clone.getElementsByTagName('svg'))) if (!s.parentElement?.closest('svg')) {
      const orig = el.querySelectorAll('svg')[Array.from(clone.querySelectorAll('svg')).indexOf(s)];
      if (orig && bigSvg(orig)) s.replaceChildren();
    }
    clone.style.margin = '0'; clone.style.position = 'static'; clone.style.boxSizing = 'border-box';
    clone.style.width = `${r1(r.width)}px`; clone.style.height = `${r1(r.height)}px`;
    const pcs = el.parentElement ? getComputedStyle(el.parentElement) : null;
    const style = pcs ? INHERITED.map(p => `${p}:${pcs.getPropertyValue(p).replace(/"/g, "'")}`).join(';') : '';
    const text = ser.serializeToString(clone);
    const key = `${style}|${text}`;
    let id = htmlIndex.get(key);
    if (id == null) { id = html.length; html.push({ style, html: text }); htmlIndex.set(key, id); }
    htmlCache.set(el, { id, w: r1(r.width), h: r1(r.height) });
    return id;
  };
  /**
   * The blocks and backgrounds of the page's or a dialog's HTML, cutting large elements and those holding drawings into
   * their children. Dialogs met in the page are left to their own layer; fixed and sticky elements of the page (and
   * what they hold) go to the layer's `fixed` blocks, in screen coordinates.
   */
  const blocksOf = (el: Element, out: Layer, page = false, stuck = false) => {
    for (const c of Array.from(el.children)) {
      const small = c instanceof SVGSVGElement && !bigSvg(c);
      if (!(c instanceof HTMLElement || small) || skipRoot(c) || SKIP_TAGS.has(c.tagName)) continue;
      if (page && c.matches('.fixed.inset-0')) continue;
      if (c instanceof SVGSVGElement && !small) continue;
      const r = c.getBoundingClientRect();
      if (r.width < 1 || r.height < 1 || r.bottom < 0 || r.top > innerHeight || r.right < 0 || r.left > innerWidth) continue;
      const cs = getComputedStyle(c);
      if (cs.visibility === 'hidden' || cs.display === 'none' || Number(cs.opacity) === 0) continue;
      const clip = clipOf(c, r);
      if (clip && clip[2] === 0) continue;
      const pinned = stuck || (page && (cs.position === 'fixed' || cs.position === 'sticky'));
      if (split(c)) {
        const bw = parseFloat(cs.borderTopWidth) || 0;
        if (cs.backgroundColor !== TRANSPARENT || (bw > 0 && cs.borderTopStyle !== 'none'))
          out.panels.push({ x: r1(r.left), y: r1(r.top), w: r1(r.width), h: r1(r.height), fill: cs.backgroundColor === TRANSPARENT ? 'none' : cs.backgroundColor, radius: parseFloat(cs.borderTopLeftRadius) || 0, ...(bw > 0 ? { stroke: cs.borderTopColor, strokeWidth: bw } : {}), clip });
        blocksOf(c, out, page, pinned);
      } else if (pinned) {
        (out.fixed ??= []).push({ x: r1(r.left), y: r1(r.top), w: r1(r.width), h: r1(r.height), id: blockId(c as HTMLElement, r), clip });
      } else {
        const [ox, oy] = out.origin;
        // inside the card: only the dialog's own scrolling boxes clip (the window is clipped once for the whole dialog)
        const inner = innerClip(c, r);
        out.blocks.push({ x: r1(r.left - ox), y: r1(r.top - oy), w: r1(r.width), h: r1(r.height), id: blockId(c as HTMLElement, r), clip: inner && [r1(inner[0] - ox), r1(inner[1] - oy), inner[2], inner[3]] });
      }
    }
  };

  const capture = () => {
    // full-screen dialogs, in page order: each above what comes before it
    const overlays = Array.from(root.querySelectorAll<HTMLElement>('.fixed.inset-0')).filter(ov => {
      if (skipRoot(ov)) return false;
      const cs = getComputedStyle(ov);
      return cs.display !== 'none' && cs.visibility !== 'hidden';
    });
    // the page: its blocks placed in page coordinates, moved by the scroll
    const page: Layer = { panels: [], blocks: [], svgs: [], origin: [r1(-scrollX), r1(-scrollY)] };
    blocksOf(root.body, page, true);
    const layers: Layer[] = [page, ...overlays.map(ov => {
      const layer: Layer = { panels: [{ x: 0, y: 0, w: innerWidth, h: innerHeight, fill: getComputedStyle(ov).backgroundColor, radius: 0 }], blocks: [], svgs: [], origin: [0, 0] };
      const card = ov.firstElementChild as HTMLElement | null;
      if (card) {
        const r = card.getBoundingClientRect(), cc = getComputedStyle(card);
        layer.origin = [r1(r.left), r1(r.top)];
        if (r.width > 0 && cc.backgroundColor !== TRANSPARENT) layer.panels.push({ x: r1(r.left), y: r1(r.top), w: r1(r.width), h: r1(r.height), fill: cc.backgroundColor, radius: parseFloat(cc.borderTopLeftRadius) || 0, clip: clipOf(card, r) });
        blocksOf(card, layer);
      }
      return layer;
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
    const key = JSON.stringify(layers);
    if (key === lastKey) return;
    lastKey = key;
    frames.push({ t: now(), layers });
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
      cache.delete(e); htmlCache.delete(e); splitCache.delete(e);
      if (e.tagName.toLowerCase() === 'svg' && bigSvg(e)) break;
    }
  };
  const mo = new MutationObserver(records => {
    let any = false;
    for (const m of records) {
      const target = m.target instanceof Element ? m.target : m.target.parentElement;
      if (target && skipRoot(target)) continue;
      dirty(m.target);
      if (m.type === 'childList') m.removedNodes.forEach(n => { if (n instanceof Element) { cache.delete(n); htmlCache.delete(n); } });
      any = true;
    }
    if (any) schedule();
  });
  mo.observe(root.body, { subtree: true, childList: true, attributes: true, characterData: true });
  const onScroll = () => schedule();
  // form controls change without touching the DOM
  const onInput = (e: Event) => { if (e.target instanceof Element && !skipRoot(e.target)) { dirty(e.target); schedule(); } };
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
  addEventListener('input', onInput, true);
  addEventListener('change', onInput, true);
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
      removeEventListener('input', onInput, true);
      removeEventListener('change', onInput, true);
      removeEventListener('pointermove', onMove, true);
      removeEventListener('pointerdown', onDown, true);
      if (timer != null) clearTimeout(timer);
      capture();
      const bs = getComputedStyle(root.body);
      return {
        width: innerWidth, height: innerHeight, background: bs.backgroundColor === TRANSPARENT ? '#ffffff' : bs.backgroundColor, fontFamily: bs.fontFamily.replace(/"/g, "'"),
        frames, chunks, html, css: html.length ? pageCss(root) : [], pointer, clicks, duration: now(),
      };
    },
  };
}
