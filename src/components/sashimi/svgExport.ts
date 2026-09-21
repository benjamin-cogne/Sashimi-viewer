/**
 * The plot's SVG element as a standalone document: interactive helpers marked data-export="skip" removed,
 * namespace declared, an explicit white background, and colours written the way every renderer accepts:
 * rgba() (undefined in SVG 1.1, refused by Inkscape, Office and librsvg) becomes rgb() + opacity.
 */
export function serializePlotSvg(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.querySelectorAll('[data-export="skip"]').forEach(n => n.remove());
  clone.removeAttribute('data-sashimi-plot'); clone.removeAttribute('data-loading'); clone.removeAttribute('style');
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  const w = clone.getAttribute('width') ?? '1200', h = clone.getAttribute('height') ?? '600';
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  bg.setAttribute('x', '0'); bg.setAttribute('y', '0'); bg.setAttribute('width', w); bg.setAttribute('height', h); bg.setAttribute('fill', '#ffffff');
  clone.insertBefore(bg, clone.firstChild);
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + portableColors(new XMLSerializer().serializeToString(clone));
}

/** fill="rgba(r, g, b, a)" → fill="rgb(r, g, b)" fill-opacity="a" (same for stroke); rgb(a) inside style attributes left as they are. */
export function portableColors(svg: string): string {
  return svg.replace(/\b(fill|stroke)="rgba\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*([\d.]+)\s*\)"/g, (_, attr, r, g, b, a) => `${attr}="rgb(${r}, ${g}, ${b})" ${attr}-opacity="${a}"`);
}

/** A file name safe on every system, from a view label. */
export const safeFileName = (s: string) => s.replace(/[^\w.\-–·]+/g, '_').replace(/[·–]/g, '-').replace(/_+/g, '_').replace(/^_|_$/g, '').slice(0, 80) || 'view';

/** The ids of one SVG document made unique (id="x", url(#x), href="#x") so several plots can share a page. */
export function suffixIds(svg: string, suffix: string): string {
  return svg
    .replace(/\bid="([^"]+)"/g, (_, id) => `id="${id}-${suffix}"`)
    .replace(/url\(#([^)]+)\)/g, (_, id) => `url(#${id}-${suffix})`)
    .replace(/\b(xlink:)?href="#([^"]+)"/g, (_, x, id) => `${x ?? ''}href="#${id}-${suffix}"`);
}

/**
 * Several plots on one page, stacked vertically under their titles. Each plot's content is inlined in a
 * translated, clipped group (no nested <svg>, which several renderers mishandle); its ids are suffixed.
 */
export function stackSvgs(plots: { title: string; svg: string }[], gap = 28, titleH = 22): string {
  const items = plots.map((p, i) => {
    const doc = suffixIds(p.svg.replace(/^<\?xml[^>]*>\s*/, ''), `v${i + 1}`);
    const open = /<svg\b[^>]*>/.exec(doc);
    const openTag = open?.[0] ?? '<svg>';
    const attr = (name: string) => new RegExp(`\\s${name}="([^"]*)"`).exec(openTag)?.[1];
    const w = Number(attr('width') ?? 1200), h = Number(attr('height') ?? 600);
    const font = attr('font-family');
    const inner = doc.slice((open?.index ?? 0) + openTag.length, doc.lastIndexOf('</svg>'));
    return { inner, w, h, font, title: p.title };
  });
  const width = Math.max(...items.map(x => x.w));
  let y = 0;
  const parts: string[] = [];
  items.forEach((it, i) => {
    const clipId = `page-clip-v${i + 1}`;
    parts.push(`<text x="8" y="${y + 15}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="13" font-weight="700" fill="#111827">${escapeXml(it.title)}</text>`);
    parts.push(`<defs><clipPath id="${clipId}"><rect x="0" y="0" width="${it.w}" height="${it.h}"/></clipPath></defs>`);
    parts.push(`<g transform="translate(0, ${y + titleH})" clip-path="url(#${clipId})"${it.font ? ` font-family="${it.font}"` : ''}>${it.inner}</g>`);
    y += titleH + it.h + gap;
  });
  const height = Math.max(0, y - gap);
  return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n<rect x="0" y="0" width="${width}" height="${height}" fill="#ffffff"/>\n${parts.join('\n')}\n</svg>\n`;
}
const escapeXml = (s: string) => s.replace(/[<>&"']/g, c => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;', "'": '&apos;' }[c]!));
