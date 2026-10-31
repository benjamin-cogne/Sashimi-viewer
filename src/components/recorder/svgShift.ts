/**
 * An SVG fragment up to a translation: the same drawing moved (a chevron following the view, a label sliding along an
 * arc) is written once in an export, with its positions in the timeline, instead of once per position.
 *
 * `shiftCanon(s)` rewrites the fragment as if moved so that its first x and first y are 0, and returns that offset; two
 * fragments with the same `key` are one drawing at two places. Only what moves with a translation is shifted: the
 * absolute coordinates (x, y, x1, cx, points, the absolute path commands); relative ones (dx, dy, the lowercase path
 * commands) stay. A leading translate() in a transform takes the whole shift for its element and everything inside it.
 * Anything a translation would not move rigidly, or that this cannot read for sure, gives null and the fragment is
 * kept as it is: a transform not starting with translate (rotate, scale, matrix), a coordinate in % or em, and the
 * gradients, patterns, masks, markers, filters and symbols (their coordinates are relative to something else).
 */

const NO_SHIFT = /<(?:linearGradient|radialGradient|pattern|mask|marker|filter|symbol)\b/;
const TAG = /<([a-zA-Z][\w:.-]*)((?:\s+[^\s=/>]+="[^"]*")*)\s*(\/?)>|<\/([a-zA-Z][\w:.-]*)\s*>/g;
const ATTR = /(\s)([^\s=/>]+)="([^"]*)"/g;
const NUM = /^-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const X_ATTRS = new Set(['x', 'x1', 'x2', 'cx']);
const Y_ATTRS = new Set(['y', 'y1', 'y2', 'cy']);
/** numbers to 0.01 px, as text (no "-0") */
export const r2 = (v: number): string => String(Math.round(v * 100) / 100 || 0);

/** Path commands: how many numbers each takes, and which of them are x (1) or y (2) when absolute. */
const PATH_ARGS: Record<string, number[]> = {
  M: [1, 2], L: [1, 2], T: [1, 2], H: [1], V: [2], C: [1, 2, 1, 2, 1, 2], S: [1, 2, 1, 2], Q: [1, 2, 1, 2], A: [0, 0, 0, 0, 0, 1, 2], Z: [],
};

interface Ref { x: number | null; y: number | null }

function shiftPath(d: string, dx: number, dy: number, ref: Ref): string | null {
  const toks = d.match(/[MmLlHhVvCcSsQqTtAaZz]|-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g);
  if (!toks) return d.trim() ? null : d;
  // anything else than commands, numbers, spaces and commas: not a path this reads
  if (d.replace(/[MmLlHhVvCcSsQqTtAaZz]|-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|[\s,]/g, '')) return null;
  const out: string[] = [];
  let cmd = '', k = 0, first = true, firstPair = false;
  for (const t of toks) {
    if (/[a-zA-Z]/.test(t)) {
      cmd = t; k = 0; out.push(t);
      // the first moveto of a path is absolute even when written lowercase
      firstPair = first && t === 'm'; first = false;
      continue;
    }
    if (!cmd) return null;
    const up = cmd.toUpperCase(), args = PATH_ARGS[up];
    if (!args || !args.length) return null;
    let v = Number(t);
    const role = args[k % args.length];
    const abs = cmd === up || (firstPair && k < 2);
    if (abs && role === 1) { if (ref.x == null) ref.x = v; v += dx; }
    else if (abs && role === 2) { if (ref.y == null) ref.y = v; v += dy; }
    out.push((out.length && /[\d.]$/.test(out[out.length - 1]) ? ' ' : '') + r2(v));
    k++;
    if (firstPair && k >= 2) firstPair = false;
  }
  return out.join('');
}

function shiftList(v: string, d: number, first: (n: number) => void): string | null {
  const parts = v.trim().split(/[\s,]+/);
  if (!parts.every(p => NUM.test(p))) return null;
  first(Number(parts[0]));
  return parts.map(p => r2(Number(p) + d)).join(' ');
}

/** The fragment moved by (dx, dy); the first absolute x and y met go to `ref`. null when it cannot be moved for sure. */
function shift(s: string, dx: number, dy: number, ref: Ref, absorbOnly: boolean): string | null {
  if (NO_SHIFT.test(s) || /style="[^"]*transform|\stransform-origin="/.test(s)) return null;
  // a clip, mask or paint defined elsewhere stays where it is while the element moves: only the fragment's own
  const own = new Set(Array.from(s.matchAll(/\sid="([^"]+)"/g), m => m[1]));
  for (const m of s.matchAll(/url\(\s*['"]?#([^'")\s]+)|href="#([^"]+)"/g)) if (!own.has(m[1] ?? m[2])) return null;
  // depth of the elements whose own translate took the shift (their content stays put)
  let absorbedAt = -1, depth = 0, rootAbsorbed = false, failed = false;
  const out = s.replace(TAG, (m, name: string | undefined, attrs: string, selfClose: string, closing: string | undefined) => {
    if (failed) return m;
    if (closing) { depth--; if (absorbedAt === depth) absorbedAt = -1; return m; }
    const here = depth;
    if (!selfClose) depth++;
    if (absorbedAt >= 0) return m;
    let absorbs = false;
    const tf = /\stransform="([^"]*)"/.exec(attrs)?.[1];
    if (tf != null && tf.trim()) {
      const t = /^\s*translate\(\s*(-?[\d.eE+-]+)(?:[\s,]+(-?[\d.eE+-]+))?\s*\)/.exec(tf);
      if (!t || !NUM.test(t[1]) || (t[2] != null && !NUM.test(t[2]))) { failed = true; return m; }
      absorbs = true;
    }
    if (here === 0 && absorbOnly && !absorbs) { failed = true; return m; }
    const newAttrs = attrs.replace(ATTR, (a, sp: string, k: string, v: string) => {
      if (failed) return a;
      if (absorbs) {
        if (k !== 'transform') return a;
        const nv = v.replace(/^(\s*translate\(\s*)(-?[\d.eE+-]+)(?:([\s,]+)(-?[\d.eE+-]+))?(\s*\))/, (_m, pre: string, tx: string, sep: string | undefined, ty: string | undefined, post: string) => {
          const x = Number(tx), y = ty != null ? Number(ty) : 0;
          if (ref.x == null) ref.x = x;
          if (ref.y == null) ref.y = y;
          return `${pre}${r2(x + dx)} ${r2(y + dy)}${post}`;
        });
        return `${sp}${k}="${nv}"`;
      }
      if (X_ATTRS.has(k) || Y_ATTRS.has(k)) {
        const isX = X_ATTRS.has(k);
        const nv = shiftList(v, isX ? dx : dy, n => { if (isX) { if (ref.x == null) ref.x = n; } else if (ref.y == null) ref.y = n; });
        if (nv == null) { failed = true; return a; }
        return `${sp}${k}="${nv}"`;
      }
      if (k === 'points') {
        const nums = v.trim().split(/[\s,]+/);
        if (nums.length % 2 || !nums.every(p => NUM.test(p))) { failed = true; return a; }
        if (ref.x == null) ref.x = Number(nums[0]);
        if (ref.y == null) ref.y = Number(nums[1]);
        return `${sp}${k}="${nums.map((p, i) => r2(Number(p) + (i % 2 ? dy : dx))).join(' ')}"`;
      }
      if (k === 'd' && name === 'path') {
        const nv = shiftPath(v, dx, dy, ref);
        if (nv == null) { failed = true; return a; }
        return `${sp}${k}="${nv}"`;
      }
      return a;
    });
    if (absorbs) { if (here === 0) rootAbsorbed = true; if (!selfClose) absorbedAt = here; }
    return `<${name}${newAttrs}${selfClose ? '/' : ''}>`;
  });
  if (failed || (absorbOnly && !rootAbsorbed)) return null;
  return out;
}

export interface Canon { key: string; x: number; y: number }

/**
 * The fragment moved to its own origin, and where that origin was. `absorbOnly`: for an element's opening tag alone
 * (its children kept apart), accepted only when its own translate takes the whole shift.
 */
export function shiftCanon(s: string, absorbOnly = false): Canon | null {
  const ref: Ref = { x: null, y: null };
  // a first pass finds the origin (the first x and y), the second moves the fragment there
  if (shift(s, 0, 0, ref, absorbOnly) == null) return null;
  const x = ref.x ?? 0, y = ref.y ?? 0;
  const key = shift(s, -x, -y, { x: null, y: null }, absorbOnly);
  return key == null ? null : { key, x: Math.round(x * 100) / 100, y: Math.round(y * 100) / 100 };
}
