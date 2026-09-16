/** The plot's SVG element as a standalone document: interactive helpers marked data-export="skip" removed, namespace declared. */
export function serializePlotSvg(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.querySelectorAll('[data-export="skip"]').forEach(n => n.remove());
  clone.removeAttribute('data-sashimi-plot'); clone.removeAttribute('data-loading');
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
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

/** Several plots on one page: stacked vertically, each under its title, as nested <svg> elements. */
export function stackSvgs(plots: { title: string; svg: string }[], gap = 28, titleH = 22): string {
  const items = plots.map((p, i) => {
    const body = suffixIds(p.svg.replace(/^<\?xml[^>]*>\s*/, ''), `v${i + 1}`);
    const w = Number(/<svg[^>]*\swidth="([\d.]+)"/.exec(body)?.[1] ?? 1200), h = Number(/<svg[^>]*\sheight="([\d.]+)"/.exec(body)?.[1] ?? 600);
    return { body, w, h, title: p.title };
  });
  const width = Math.max(...items.map(x => x.w));
  let y = 0;
  const parts: string[] = [];
  for (const it of items) {
    parts.push(`<text x="8" y="${y + 15}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="13" font-weight="700" fill="#111827">${escapeXml(it.title)}</text>`);
    parts.push(it.body.replace(/^<svg\b/, `<svg y="${y + titleH}"`));
    y += titleH + it.h + gap;
  }
  const height = Math.max(0, y - gap);
  return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n<rect width="${width}" height="${height}" fill="#ffffff"/>\n${parts.join('\n')}\n</svg>\n`;
}
const escapeXml = (s: string) => s.replace(/[<>&"']/g, c => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;', "'": '&apos;' }[c]!));
