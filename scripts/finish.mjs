// Runs after `vite build`: dist/index.html is the single-file page. Keep it there for GitHub Pages,
// copy it as dist/sashimi-viewer.html (release asset) and as ./sashimi-viewer.html (tracked in git so
// the file can be downloaded straight from the repository).
import { copyFileSync, existsSync, statSync } from 'node:fs';
import { join } from 'node:path';

const built = join('dist', 'index.html');
if (!existsSync(built)) { console.error('finish: dist/index.html not found, run `vite build` first'); process.exit(1); }
copyFileSync(built, join('dist', 'sashimi-viewer.html'));
copyFileSync(built, 'sashimi-viewer.html');
console.log(`sashimi-viewer.html updated (${(statSync(built).size / 1024).toFixed(0)} kB)`);
