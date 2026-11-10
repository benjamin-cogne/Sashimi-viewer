// Runs after `vite build`: dist/index.html is the single-file page. Keep it there for GitHub Pages,
// copy it as dist/sashimi-viewer.html (release asset) and as ./sashimi-viewer.html (tracked in git so
// the file can be downloaded straight from the repository).
//
// The bytes are written, not copied with copyFileSync: a copy also sets the source's permission bits on the
// destination, which only the file's owner may do, so on a server where the checkout belongs to another user
// (cloned with sudo, group-writable) it failed with EPERM, and the failed copy deleted the tracked page. The copy
// at the root is a convenience for the repository: a build that cannot write it still succeeds, dist/ is the page.
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const built = join('dist', 'index.html');
if (!existsSync(built)) { console.error('finish: dist/index.html not found, run `vite build` first'); process.exit(1); }
const page = readFileSync(built);
writeFileSync(join('dist', 'sashimi-viewer.html'), page);
const kb = (page.length / 1024).toFixed(0);
try {
  writeFileSync('sashimi-viewer.html', page);
  console.log(`sashimi-viewer.html updated (${kb} kB)`);
} catch (e) {
  console.warn(`finish: ./sashimi-viewer.html not updated (${e.code ?? e.message}: this user may not write it); the page is in dist/index.html and dist/sashimi-viewer.html (${kb} kB)`);
}
