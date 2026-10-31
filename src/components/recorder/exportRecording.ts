/**
 * A recording as one animated SVG (or an HTML page holding it, with play / pause and a time bar).
 *
 * The elements copied during the recording go once into <defs>; each frame is a group of <use> references, shown
 * from its time to the next frame's by a CSS animation (visibility, step-end). The pointer moves by a transform
 * animation through its positions; each click is a ring that widens and fades. Every animation has the length of the
 * whole recording, so they stay in step when the player pauses or seeks them (Web Animations API), and when looping.
 */
import type { Recording, Chunk } from './recorder';

export interface ExportOptions {
  /** idle stretches (nothing drawn changes, the pointer still) longer than this are cut down to it; 0 keeps them */
  maxIdleMs: number;
  loop: boolean;
  /** text drawn in the recording replaced at export: sample names → "Patient 1"… */
  replace: { from: string; to: string }[];
  title: string;
}
/** The idle stretch kept by default, in ms. */
export const MAX_IDLE_DEFAULT_MS = 1000;
/** Pause at the end before a loop starts again, in ms. */
const END_HOLD_MS = 1500;
const CLICK_MS = 600;

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const pct = (t: number, d: number) => `${Math.min(100, (t / d) * 100).toFixed(3).replace(/\.?0+$/, '')}%`;

/** A function from recorded time to played time, cutting idle stretches to `maxIdle`. */
function timeMap(rec: Recording, maxIdle: number): (t: number) => number {
  if (!(maxIdle > 0)) return t => t;
  const ev = [0, ...rec.frames.map(f => f.t), ...rec.pointer.map(p => p.t), ...rec.clicks.map(c => c.t), rec.duration].sort((a, b) => a - b);
  // each gap longer than maxIdle loses what lies past its first maxIdle
  const cuts: { s: number; e: number }[] = [];
  for (let i = 1; i < ev.length; i++) if (ev[i] - ev[i - 1] > maxIdle) cuts.push({ s: ev[i - 1] + maxIdle, e: ev[i] });
  return t => t - cuts.reduce((acc, c) => acc + Math.min(Math.max(t - c.s, 0), c.e - c.s), 0);
}

/** The text between tags with the replacements made (attributes and tag names untouched). */
function replaceText(s: string, rep: ExportOptions['replace']): string {
  if (!rep.length) return s;
  const list = rep.filter(r => r.from).map(r => ({ from: esc(r.from), to: esc(r.to) })).sort((a, b) => b.from.length - a.from.length);
  return s.replace(/>([^<]+)</g, (_m, txt: string) => `>${list.reduce((acc, r) => acc.split(r.from).join(r.to), txt)}<`);
}

export interface Exported { svg: string; duration: number; frames: number; bytes: number }

export function recordingToSvg(rec: Recording, opts: ExportOptions): Exported {
  const map = timeMap(rec, opts.maxIdleMs);
  const end = map(rec.duration);
  const D = Math.max(1, end + (opts.loop ? END_HOLD_MS : 0));
  // the area drawn over the recording: every copied SVG (clipped to what showed), with a margin
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const f of rec.frames) for (const s of f.svgs) {
    const [x, y, w, h] = s.clip ?? [s.x, s.y, s.w, s.h];
    x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x + w); y1 = Math.max(y1, y + h);
  }
  if (!Number.isFinite(x0)) { x0 = 0; y0 = 0; x1 = rec.width; y1 = rec.height; }
  const M = 8;
  x0 = Math.max(0, x0 - M); y0 = Math.max(0, y0 - M); x1 = Math.min(rec.width, x1 + M); y1 = Math.min(rec.height, y1 + M);
  const W = Math.round(x1 - x0), H = Math.round(y1 - y0);

  // definitions: every chunk a frame reaches, once
  const used = new Set<number>();
  const visit = (id: number) => { if (used.has(id)) return; used.add(id); const c = rec.chunks[id]; if ('kids' in c) c.kids.forEach(visit); };
  rec.frames.forEach(f => f.svgs.forEach(s => s.kids.forEach(visit)));
  const withId = (tagText: string, id: string): string | null => (/^<[^>]*\sid="/.test(tagText) ? null : tagText.replace(/^<([\w:-]+)/, `<$1 id="${id}"`));
  const defOf = (id: number, c: Chunk): string => {
    if ('leaf' in c) { const s = replaceText(c.leaf, opts.replace); return withId(s, `c${id}`) ?? `<g id="c${id}">${s}</g>`; }
    const inner = c.kids.map(k => `<use href="#c${k}"/>`).join('');
    const open = replaceText(`${c.open}`, opts.replace);
    const tagged = withId(open, `c${id}`);
    return tagged ? `${tagged}${inner}${c.close}` : `<g id="c${id}">${open}${inner}${c.close}</g>`;
  };
  const defs = [...used].sort((a, b) => a - b).map(id => defOf(id, rec.chunks[id])).join('\n');

  // clip rectangles, shared
  const clips = new Map<string, string>();
  const clipId = (c: [number, number, number, number]) => {
    const k = c.join(',');
    if (!clips.has(k)) clips.set(k, `k${clips.size}`);
    return clips.get(k)!;
  };
  const css: string[] = [
    `.f{visibility:hidden;animation:${D}ms step-end ${opts.loop ? 'infinite' : '1'} forwards}`,
    `.a{animation-duration:${D}ms;animation-iteration-count:${opts.loop ? 'infinite' : '1'};animation-fill-mode:forwards}`,
  ];
  const frames: string[] = [];
  rec.frames.forEach((f, i) => {
    const a = map(f.t), b = i + 1 < rec.frames.length ? map(rec.frames[i + 1].t) : Infinity;
    if (b <= a && i + 1 < rec.frames.length) return; // shown for no time once the pauses are cut
    css.push(`@keyframes f${i}{0%{visibility:hidden}${pct(a, D)}{visibility:visible}${b === Infinity ? '' : `${pct(b, D)}{visibility:hidden}`}}#f${i}{animation-name:f${i}}`);
    const panel = (p: Recording['frames'][number]['panels'][number]) => `<rect x="${p.x}" y="${p.y}" width="${p.w}" height="${p.h}" rx="${p.radius}" fill="${p.fill}"/>`;
    const svgEl = (s: Recording['frames'][number]['svgs'][number]) => {
      const inner = `<svg x="${s.x}" y="${s.y}" width="${s.w}" height="${s.h}"${s.viewBox ? ` viewBox="${s.viewBox}"` : ''}${s.par ? ` preserveAspectRatio="${s.par}"` : ''} style="${esc(s.style)}">${s.kids.map(k => `<use href="#c${k}"/>`).join('')}</svg>`;
      return s.clip ? `<g clip-path="url(#${clipId(s.clip)})">${inner}</g>` : inner;
    };
    // layer by layer: the page, then each dialog (its backdrop, its card, its drawings)
    let pi = 0, si = 0, body = '';
    for (const [np, ns] of f.order) {
      body += f.panels.slice(pi, pi + np).map(panel).join('') + f.svgs.slice(si, si + ns).map(svgEl).join('');
      pi += np; si += ns;
    }
    frames.push(`<g id="f${i}" class="f">${body}</g>`);
  });
  const clipDefs = [...clips].map(([k, id]) => { const [x, y, w, h] = k.split(','); return `<clipPath id="${id}"><rect x="${x}" y="${y}" width="${w}" height="${h}"/></clipPath>`; }).join('');

  // pointer: a transform animation through its positions
  let pointer = '';
  if (rec.pointer.length) {
    const kf = rec.pointer.map(p => `${pct(map(p.t), D)}{transform:translate(${Math.round(p.x)}px,${Math.round(p.y)}px)}`);
    const first = rec.pointer[0];
    css.push(`@keyframes ptr{0%{transform:translate(${Math.round(first.x)}px,${Math.round(first.y)}px)}${kf.join('')}}#ptr{animation-name:ptr;animation-timing-function:linear}`);
    pointer = `<g id="ptr" class="a"><path d="M0 0 L0 16 L4.5 12.2 L7.6 19 L10.2 17.8 L7.2 11.2 L13 11.2 Z" fill="#111827" stroke="#ffffff" stroke-width="1.3" stroke-linejoin="round"/></g>`;
  }
  const clickEls = rec.clicks.map((c, i) => {
    const a = map(c.t), b = a + CLICK_MS;
    css.push(`@keyframes ck${i}{0%,${pct(a, D)}{opacity:0;transform:scale(.3)}${pct(a + 40, D)}{opacity:.9;transform:scale(.5)}${pct(b, D)},100%{opacity:0;transform:scale(1.6)}}#ck${i}{animation-name:ck${i};animation-timing-function:ease-out}`);
    return `<circle id="ck${i}" class="a" cx="${Math.round(c.x)}" cy="${Math.round(c.y)}" r="14" fill="none" stroke="#0f766e" stroke-width="3" style="opacity:0;transform-box:fill-box;transform-origin:center"/>`;
  }).join('');

  const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${Math.round(x0)} ${Math.round(y0)} ${W} ${H}" width="${W}" height="${H}" font-family="${esc(rec.fontFamily)}">
<title>${esc(opts.title)}</title>
<style>${css.join('\n')}</style>
<defs>${clipDefs}
${defs}</defs>
<rect x="${Math.round(x0)}" y="${Math.round(y0)}" width="${W}" height="${H}" fill="${rec.background}"/>
${frames.join('\n')}
${clickEls}${pointer}
</svg>`;
  return { svg, duration: D, frames: frames.length, bytes: new Blob([svg]).size };
}

/** An HTML page around the SVG: play / pause, a time bar, the time; space toggles. */
export function svgToHtml(svg: string, title: string, durationMs: number): string {
  return `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>${esc(title)}</title>
<style>
body{margin:0;background:#f4f6f8;font:13px Inter,system-ui,-apple-system,'Segoe UI',sans-serif;color:#1f2937}
.wrap{max-width:100%;padding:16px;box-sizing:border-box}
.stage{background:#fff;border-radius:12px;box-shadow:0 1px 3px rgba(0,0,0,.12);overflow:hidden}
.stage svg{display:block;width:100%;height:auto}
.bar{display:flex;align-items:center;gap:10px;margin-top:10px}
.bar button{width:34px;height:30px;border-radius:8px;border:1px solid #d3dbe5;background:#fff;cursor:pointer;font-size:14px}
.bar input{flex:1;accent-color:#0f766e}
.time{font-variant-numeric:tabular-nums;color:#5b6b7e;min-width:90px;text-align:right}
</style></head><body><div class="wrap">
<div class="stage" id="stage">${svg}</div>
<div class="bar"><button id="play" title="Play / pause (space)">❚❚</button><input id="seek" type="range" min="0" max="${Math.round(durationMs)}" step="10" value="0"><span class="time" id="time"></span></div>
</div>
<script>
(function(){
  var D=${Math.round(durationMs)}, svg=document.querySelector('#stage svg'), play=document.getElementById('play'), seek=document.getElementById('seek'), time=document.getElementById('time');
  var anims=function(){return svg.getAnimations({subtree:true});}, paused=false, dragging=false;
  var fmt=function(ms){var s=Math.max(0,ms)/1000;return Math.floor(s/60)+':'+('0'+Math.floor(s%60)).slice(-2);};
  var cur=function(){var a=anims()[0];return a?(Number(a.currentTime)||0)%D:0;};
  var set=function(ms){anims().forEach(function(a){a.currentTime=ms;});};
  var toggle=function(){paused=!paused;anims().forEach(function(a){paused?a.pause():a.play();});play.textContent=paused?'▶':'❚❚';};
  play.onclick=toggle;
  document.addEventListener('keydown',function(e){if(e.code==='Space'){e.preventDefault();toggle();}});
  seek.oninput=function(){dragging=true;set(Number(seek.value));};
  seek.onchange=function(){dragging=false;};
  (function tick(){var t=cur();if(!dragging)seek.value=String(Math.round(t));time.textContent=fmt(t)+' / '+fmt(D);requestAnimationFrame(tick);})();
})();
</script></body></html>`;
}
