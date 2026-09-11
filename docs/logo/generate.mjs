// Rebuilds the Sashimi viewer logo assets from the mark geometry and the wordmark font.
//
//   cd docs/logo && npm i opentype.js@1.3.4 sharp@0.33 && node generate.mjs
//
// Needs IBM Plex Sans 400 and 600 as plex-400.ttf / plex-600.ttf next to this file (Google Fonts,
// SIL Open Font License). The wordmark is outlined to paths so the SVGs render without the font.
import opentype from 'opentype.js';
import sharp from 'sharp';
import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = HERE;
mkdirSync(join(OUT, 'png'), { recursive: true });
const semi = opentype.loadSync(join(HERE, 'plex-600.ttf'));
const reg = opentype.loadSync(join(HERE, 'plex-400.ttf'));

// Palette: the viewer's indigo one step darker; the red is the viewer's "junction absent from the
// comparison samples" colour, so the mark uses the plot's own semantics.
const P = {
  light: { ink: '#3730a3', red: '#dc2626', text: '#1c1b33', muted: '#5f5e7a' },
  dark:  { ink: '#a5b4fc', red: '#f87171', text: '#eceaf6', muted: '#9b99b5' },
};

// The mark, 64 x 64: three exons, two canonical arcs, one dashed aberrant arc skipping the middle exon.
const mark = (c) => `
  <path d="M11 46 Q21.5 14 32 46" fill="none" stroke="${c.ink}" stroke-width="3.5" stroke-linecap="round"/>
  <path d="M32 46 Q42.5 14 53 46" fill="none" stroke="${c.ink}" stroke-width="3.5" stroke-linecap="round"/>
  <path d="M11 46 Q32 -26 53 46" fill="none" stroke="${c.red}" stroke-width="3.5" stroke-linecap="round" stroke-dasharray="5 4"/>
  <path d="M18 51 H25 M39 51 H46" stroke="${c.ink}" stroke-width="2.5" stroke-linecap="round"/>
  <rect x="4" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>
  <rect x="25" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>
  <rect x="46" y="46" width="14" height="10" rx="1.5" fill="${c.ink}"/>`;

// Favicon cut for 16-32 px: the canonical arcs and the dashes go, the exons and the aberrant arc stay,
// heavier so they survive downsampling.
const favicon = (c) => `
  <path d="M11 44 Q32 -22 53 44" fill="none" stroke="${c.red}" stroke-width="6" stroke-linecap="round"/>
  <rect x="2" y="44" width="16" height="14" rx="2" fill="${c.ink}"/>
  <rect x="24" y="44" width="16" height="14" rx="2" fill="${c.ink}"/>
  <rect x="46" y="44" width="16" height="14" rx="2" fill="${c.ink}"/>`;

const svg = (w, h, inner) => `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}">${inner}\n</svg>\n`;
const word = (font, text, size) => ({ d: font.getPath(text, 0, 0, size).toPathData(2), adv: font.getAdvanceWidth(text, size) });

const files = {};
for (const theme of ['light', 'dark']) {
  const c = P[theme], sfx = theme === 'dark' ? '-dark' : '';
  files[`logo${sfx}.svg`] = svg(64, 64, mark(c));
  files[`favicon${sfx}.svg`] = svg(64, 64, favicon(c));
  const size = 30, w1 = word(semi, 'Sashimi', size), w2 = word(reg, 'viewer', size);
  const tx = 80, tx2 = tx + w1.adv + 8, totalW = Math.ceil(tx2 + w2.adv + 4);
  files[`logo-lockup${sfx}.svg`] = svg(totalW, 64, `
  <g>${mark(c)}\n  </g>
  <path transform="translate(${tx} 43)" d="${w1.d}" fill="${c.text}"/>
  <path transform="translate(${tx2.toFixed(1)} 43)" d="${w2.d}" fill="${c.muted}"/>`);
}
for (const [name, body] of Object.entries(files)) writeFileSync(join(OUT, name), body);

// Rasters: PNG at the sizes browsers and stores ask for, plus a favicon.ico wrapping the 16/32/48 PNGs.
const png = async (svgText, size) => sharp(Buffer.from(svgText), { density: 72 * size / 64 }).resize(size, size).png().toBuffer();
const outPng = {};
for (const size of [16, 32, 48]) outPng[size] = await png(files['favicon.svg'], size);
for (const size of [64, 128, 180, 256, 512]) outPng[size] = await png(files['logo.svg'], size);
for (const [size, buf] of Object.entries(outPng)) writeFileSync(join(OUT, 'png', `logo-${size}.png`), buf);

// ICO container with PNG entries (supported since Windows Vista, all current browsers).
const entries = [16, 32, 48].map(s => outPng[s]);
const header = Buffer.alloc(6); header.writeUInt16LE(0, 0); header.writeUInt16LE(1, 2); header.writeUInt16LE(entries.length, 4);
let offset = 6 + 16 * entries.length; const dirs = [], datas = [];
entries.forEach((buf, i) => {
  const s = [16, 32, 48][i], d = Buffer.alloc(16);
  d.writeUInt8(s, 0); d.writeUInt8(s, 1); d.writeUInt8(0, 2); d.writeUInt8(0, 3);
  d.writeUInt16LE(1, 4); d.writeUInt16LE(32, 6); d.writeUInt32LE(buf.length, 8); d.writeUInt32LE(offset, 12);
  offset += buf.length; dirs.push(d); datas.push(buf);
});
writeFileSync(join(OUT, 'favicon.ico'), Buffer.concat([header, ...dirs, ...datas]));
console.log('written:', Object.keys(files).join(', '), '+ png/ + favicon.ico');
