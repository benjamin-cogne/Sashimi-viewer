import opentype from 'opentype.js';
import { writeFileSync, mkdirSync } from 'node:fs';

const OUT = new URL('.', import.meta.url).pathname; // run from a folder containing plex-400.ttf / plex-600.ttf (IBM Plex Sans, Google Fonts) with opentype.js installed
mkdirSync(OUT, { recursive: true });
const semi = opentype.loadSync('plex-600.ttf');
const reg = opentype.loadSync('plex-400.ttf');

// ---- palette (light ground / dark ground) ----
const P = {
  light: { ink: '#3730a3', tint: '#c7d2fe', red: '#dc2626', text: '#1c1b33', muted: '#5f5e7a', white: '#ffffff' },
  dark:  { ink: '#a5b4fc', tint: '#3730a3', red: '#f87171', text: '#eceaf6', muted: '#9b99b5', white: '#ffffff' },
};

// ---- marks: 64 x 64 viewBox, inner markup only ----
const marks = {
  'a-junction': (c) => `
  <path d="M14 42 Q32 -12 50 42" fill="none" stroke="${c.ink}" stroke-width="4.5" stroke-linecap="round"/>
  <path d="M24 48 H40" stroke="${c.ink}" stroke-width="3" stroke-linecap="round"/>
  <rect x="4" y="42" width="20" height="12" rx="2" fill="${c.ink}"/>
  <rect x="40" y="42" width="20" height="12" rx="2" fill="${c.ink}"/>`,
  'b-tile': (c, onDark) => `
  <rect width="64" height="64" rx="14" fill="${onDark ? '#4338ca' : c.ink}"/>
  <path d="M17 28 Q32 2 47 28" fill="none" stroke="${onDark ? '#c7d2fe' : c.tint}" stroke-width="3" stroke-linecap="round"/>
  <rect x="8" y="28" width="18" height="22" rx="2.5" fill="${c.white}"/>
  <rect x="38" y="34" width="18" height="16" rx="2.5" fill="${c.white}"/>
  <rect x="8" y="52" width="48" height="3" rx="1.5" fill="${c.white}" opacity="0.7"/>`,
  'c-skip': (c) => `
  <path d="M11 46 Q21.5 14 32 46" fill="none" stroke="${c.ink}" stroke-width="3.5" stroke-linecap="round"/>
  <path d="M32 46 Q42.5 14 53 46" fill="none" stroke="${c.ink}" stroke-width="3.5" stroke-linecap="round"/>
  <path d="M11 46 Q32 -26 53 46" fill="none" stroke="${c.red}" stroke-width="3.5" stroke-linecap="round" stroke-dasharray="5 4"/>
  <path d="M18 51 H25 M39 51 H46" stroke="${c.ink}" stroke-width="2.5" stroke-linecap="round"/>
  <rect x="4" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>
  <rect x="25" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>
  <rect x="46" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>`,
};

const svg = (w, h, inner) => `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}">${inner}\n</svg>\n`;

// ---- wordmark: outlined IBM Plex Sans, "Sashimi" 600 + "viewer" 400 ----
function word(font, text, x, y, size, fill) {
  const p = font.getPath(text, x, y, size);
  const adv = font.getAdvanceWidth(text, size);
  return { d: p.toPathData(2), adv, fill };
}

const exported = {};
for (const [name, fn] of Object.entries(marks)) {
  for (const theme of ['light', 'dark']) {
    const c = P[theme];
    const inner = fn(c, theme === 'dark');
    const markSvg = svg(64, 64, inner);
    writeFileSync(`${OUT}/${name}${theme === 'dark' ? '-dark' : ''}.svg`, markSvg);
    // lockup: mark 48px tall at x=0, text baseline aligned to mark's optical centre
    const size = 30;
    const w1 = word(semi, 'Sashimi', 0, 0, size, c.text);
    const w2 = word(reg, 'viewer', 0, 0, size, c.muted);
    const gap = 8, tx = 64 + 16;
    const tx2 = tx + w1.adv + gap;
    const totalW = Math.ceil(tx2 + w2.adv + 4);
    const baseline = 43;
    const lock = svg(totalW, 64, `
  <g>${inner}\n  </g>
  <path transform="translate(${tx} ${baseline})" d="${w1.d}" fill="${w1.fill}"/>
  <path transform="translate(${tx2.toFixed(1)} ${baseline})" d="${w2.d}" fill="${w2.fill}"/>`);
    writeFileSync(`${OUT}/${name}-lockup${theme === 'dark' ? '-dark' : ''}.svg`, lock);
    exported[`${name}-${theme}`] = { lockupWidth: totalW };
  }
}
console.log(JSON.stringify(exported));
