/**
 * Captions and chapters of a recording.
 *
 * - title card: a centred card with a title and a subtitle, over the whole view (chapters can show their name so);
 * - lower third: a dark bar at the bottom left, one to three lines of narration;
 * - callout: a label with a leader line pointing at a place of the drawing, which follows the camera when it zooms;
 * - highlight: a ring around an area of the drawing, optionally dimming the rest (spotlight).
 * Title cards and lower thirds are fixed on the screen; callouts and highlights belong to the drawing. All are text and
 * shapes, faded in and out. Times are recorded times (ms), so cutting the pauses keeps them on their moment.
 */

export type CaptionKind = 'title' | 'lower' | 'callout' | 'highlight';
export interface Caption {
  kind: CaptionKind; start: number; end: number;
  text: string; sub?: string;
  /** callout: the place pointed at (page coordinates) */
  at?: [number, number];
  /** highlight: the area (page coordinates), and whether the rest is dimmed */
  rect?: [number, number, number, number]; dim?: boolean;
}
export interface Chapter { t: number; name: string; card: boolean }

export const CAPTION_DEFAULT_MS: Record<CaptionKind, number> = { title: 3000, lower: 4000, callout: 3000, highlight: 3000 };
export const CHAPTER_CARD_MS = 2200;
const FADE_MS = 250;
const NAVY = 'rgba(15,27,43,0.88)', ACCENT = '#14857a';

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const pct = (t: number, d: number) => `${Math.max(0, Math.min(100, (t / d) * 100)).toFixed(3).replace(/\.?0+$/, '')}%`;

let ctx: CanvasRenderingContext2D | null = null;
/** Lines of `text` no wider than `max` px in the given font (measured by the browser; a rough estimate elsewhere). */
function wrap(text: string, px: number, weight: number, family: string, max: number, maxLines = 3): string[] {
  try { ctx ??= document.createElement('canvas').getContext('2d'); } catch { ctx = null; }
  const width = (s: string) => { if (!ctx) return s.length * px * 0.55; ctx.font = `${weight} ${px}px ${family}`; return ctx.measureText(s).width; };
  const lines: string[] = [];
  for (const para of text.split(/\n/)) {
    let line = '';
    for (const w of para.split(/\s+/).filter(Boolean)) {
      const next = line ? `${line} ${w}` : w;
      if (line && width(next) > max) { lines.push(line); line = w; } else line = next;
    }
    lines.push(line);
  }
  if (lines.length > maxLines) { const kept = lines.slice(0, maxLines); kept[maxLines - 1] += ' …'; return kept; }
  return lines;
}

export interface CaptionLayers { css: string[]; screen: string; content: string }

/** The fade of a caption: its element id, shown from `a` to `b` (played ms), fading in and out over FADE_MS. */
export type CaptionFade = [string, number, number];
export const CAPTION_FADE_MS = FADE_MS;

/**
 * The captions as SVG: `screen` goes above the camera (view coordinates 0..W, 0..H), `content` inside it (page
 * coordinates, so it zooms with the drawing). `map` turns recorded into played time, `D` is the animation length.
 * With `fades`, the fades are listed there for a player to run instead of written as CSS animations.
 */
export function renderCaptions(captions: Caption[], chapters: Chapter[], map: (t: number) => number, D: number,
  view: { x0: number; y0: number; W: number; H: number }, family: string, fades?: CaptionFade[]): CaptionLayers {
  const { x0, y0, W, H } = view;
  const css: string[] = [];
  let screen = '', content = '', n = 0;
  const fade = (a: number, b: number) => {
    const id = `cp${n++}`;
    if (fades) { fades.push([id, Math.round(a), Math.round(b)]); return id; }
    css.push(`@keyframes ${id}{0%,${pct(a, D)}{opacity:0}${pct(a + FADE_MS, D)},${pct(b, D)}{opacity:1}${pct(b + FADE_MS, D)},100%{opacity:0}}#${id}{animation-name:${id};animation-timing-function:linear}`);
    return id;
  };
  const card = (title: string, sub: string, a: number, b: number) => {
    const tp = Math.max(18, Math.round(W / 30)), sp = Math.max(13, Math.round(W / 62));
    const maxW = W * 0.62;
    const tl = wrap(title, tp, 700, family, maxW, 2), sl = sub ? wrap(sub, sp, 400, family, maxW, 2) : [];
    const pad = tp * 0.9, h = pad * 2 + tl.length * tp * 1.2 + (sl.length ? sp * 0.6 + sl.length * sp * 1.35 : 0), w = W * 0.7;
    const x = (W - w) / 2, y = (H - h) / 2;
    let ty = y + pad + tp;
    const lines = tl.map(l => { const s = `<text x="${W / 2}" y="${ty.toFixed(1)}" font-size="${tp}" font-weight="700" fill="#fff" text-anchor="middle">${esc(l)}</text>`; ty += tp * 1.2; return s; });
    ty += sp * 0.6 - tp * 0.2;
    const subs = sl.map(l => { const s = `<text x="${W / 2}" y="${ty.toFixed(1)}" font-size="${sp}" fill="#cbe9e4" text-anchor="middle">${esc(l)}</text>`; ty += sp * 1.35; return s; });
    screen += `<g id="${fade(a, b)}" class="a" style="opacity:0"><rect width="${W}" height="${H}" fill="rgba(15,27,43,0.35)"/><rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${w.toFixed(1)}" height="${h.toFixed(1)}" rx="${(tp * 0.6).toFixed(1)}" fill="${NAVY}"/><rect x="${(W / 2 - tp).toFixed(1)}" y="${(y + pad * 0.45).toFixed(1)}" width="${2 * tp}" height="3" rx="1.5" fill="${ACCENT}"/>${lines.join('')}${subs.join('')}</g>`;
  };
  // lower thirds shown at the same time stack upwards: each takes the lowest free slot
  const slots = new Map<Caption, number>();
  const busy: number[] = [];
  for (const c of captions.filter(x => x.kind === 'lower').sort((x, y) => x.start - y.start)) {
    let k = busy.findIndex(end => end <= c.start);
    if (k < 0) { k = busy.length; busy.push(0); }
    busy[k] = c.end + FADE_MS;
    slots.set(c, k);
  }
  const lowerHeight = (text: string) => { const fp = Math.max(13, Math.round(W / 58)), pad = fp * 0.7; return wrap(text, fp, 500, family, W * 0.6 - 2 * pad).length * fp * 1.35 + 2 * pad - fp * 0.35; };
  for (const c of captions) {
    const a = map(c.start), b = Math.max(a + 300, map(c.end));
    if (c.kind === 'title') card(c.text, c.sub ?? '', a, b);
    else if (c.kind === 'lower') {
      const fp = Math.max(13, Math.round(W / 58)), pad = fp * 0.7, maxW = W * 0.6;
      const lines = wrap(c.text, fp, 500, family, maxW - 2 * pad);
      try { ctx!.font = `500 ${fp}px ${family}`; } catch { /* no canvas */ }
      const tw = Math.min(maxW - 2 * pad, Math.max(...lines.map(l => (ctx ? ctx.measureText(l).width : l.length * fp * 0.55))));
      const w = tw + 2 * pad + 6, h = lines.length * fp * 1.35 + 2 * pad - fp * 0.35;
      // above the lower thirds shown at the same time in the slots below
      const k = slots.get(c) ?? 0;
      let below = 0;
      for (let j = 0; j < k; j++) {
        const under = captions.filter(o => o.kind === 'lower' && slots.get(o) === j && o.start < c.end && o.end > c.start);
        below += Math.max(0, ...under.map(o => lowerHeight(o.text))) + fp * 0.5;
      }
      const x = W * 0.03, y = H - H * 0.05 - h - below;
      const texts = lines.map((l, i) => `<text x="${(x + 6 + pad).toFixed(1)}" y="${(y + pad + fp * 0.95 + i * fp * 1.35).toFixed(1)}" font-size="${fp}" font-weight="500" fill="#fff">${esc(l)}</text>`).join('');
      screen += `<g id="${fade(a, b)}" class="a" style="opacity:0"><rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${w.toFixed(1)}" height="${h.toFixed(1)}" rx="${(fp * 0.4).toFixed(1)}" fill="${NAVY}"/><rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="6" height="${h.toFixed(1)}" rx="3" fill="${ACCENT}"/>${texts}</g>`;
    } else if (c.kind === 'callout' && c.at) {
      const fp = Math.max(12, Math.round(W / 75)), pad = fp * 0.55;
      const lines = wrap(c.text, fp, 600, family, W * 0.3);
      try { ctx!.font = `600 ${fp}px ${family}`; } catch { /* no canvas */ }
      const tw = Math.max(...lines.map(l => (ctx ? ctx.measureText(l).width : l.length * fp * 0.55)));
      const w = tw + 2 * pad, h = lines.length * fp * 1.3 + 2 * pad - fp * 0.3;
      const [px, py] = c.at;
      // the label above and to the right of the point, or to the left near the right edge
      const left = px + 40 + w > x0 + W;
      const bx = left ? px - 40 - w : px + 40, by = Math.max(y0 + 4, py - 40 - h);
      const lx = left ? bx + w : bx, ly = by + h;
      const texts = lines.map((l, i) => `<text x="${(bx + pad).toFixed(1)}" y="${(by + pad + fp * 0.9 + i * fp * 1.3).toFixed(1)}" font-size="${fp}" font-weight="600" fill="#fff">${esc(l)}</text>`).join('');
      content += `<g id="${fade(a, b)}" class="a" style="opacity:0"><line x1="${px}" y1="${py}" x2="${lx.toFixed(1)}" y2="${ly.toFixed(1)}" stroke="${ACCENT}" stroke-width="2"/><circle cx="${px}" cy="${py}" r="5" fill="${ACCENT}" stroke="#fff" stroke-width="2"/><rect x="${bx.toFixed(1)}" y="${by.toFixed(1)}" width="${w.toFixed(1)}" height="${h.toFixed(1)}" rx="${(fp * 0.4).toFixed(1)}" fill="${ACCENT}"/>${texts}</g>`;
    } else if (c.kind === 'highlight' && c.rect) {
      const [rx, ry, rw, rh] = c.rect;
      const dim = c.dim ? `<path d="M${x0} ${y0}h${W}v${H}h${-W}Z M${rx} ${ry}v${rh}h${rw}v${-rh}Z" fill-rule="evenodd" fill="rgba(15,27,43,0.45)"/>` : '';
      const label = c.text ? `<text x="${rx}" y="${(ry - 8).toFixed(1)}" font-size="${Math.max(12, Math.round(W / 75))}" font-weight="700" fill="${ACCENT}" stroke="#fff" stroke-width="3" paint-order="stroke">${esc(c.text)}</text>` : '';
      content += `<g id="${fade(a, b)}" class="a" style="opacity:0">${dim}<rect x="${rx}" y="${ry}" width="${rw}" height="${rh}" rx="8" fill="none" stroke="${ACCENT}" stroke-width="3"/>${label}</g>`;
    }
  }
  // chapters shown as a title card, faded in before their mark so a pause there shows it whole
  for (const ch of chapters) if (ch.card && ch.name) {
    const t = map(ch.t);
    card(ch.name, '', Math.max(0, t - FADE_MS), t + CHAPTER_CARD_MS);
  }
  return { css, screen, content };
}

/** Captions the recording suggests: each dialog opened, by the first words of its header. */
export function suggestCaptions(rec: import('./recorder').Recording): Caption[] {
  const out: Caption[] = [];
  let prev = 1;
  for (const f of rec.frames) {
    const n = f.layers.length;
    if (n > prev) {
      const layer = f.layers[n - 1];
      const first = layer.blocks.slice().sort((a, b) => a.y - b.y || a.x - b.x)[0];
      const text = first ? rec.html[first.id].html.replace(/<[^>]+>/g, ' ').replace(/&amp;/g, '&').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/\s+/g, ' ').trim().split(' · ')[0].slice(0, 70) : '';
      if (text) out.push({ kind: 'lower', start: f.t + 300, end: f.t + 300 + CAPTION_DEFAULT_MS.lower, text });
    }
    prev = n;
  }
  return out;
}
