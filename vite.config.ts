import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { viteSingleFile } from 'vite-plugin-singlefile'
import { readFileSync } from 'node:fs'

const pkg = JSON.parse(readFileSync(new URL('./package.json', import.meta.url), 'utf8')) as { version: string }

// Builds the viewer as ONE self-contained HTML file: dist/index.html, then copied by
// scripts/finish.mjs to dist/sashimi-viewer.html and to ./sashimi-viewer.html (tracked in git).
export default defineConfig({
  plugins: [react(), viteSingleFile()],
  define: { __APP_VERSION__: JSON.stringify(pkg.version) },
  server: { port: 3000 },
  build: {
    outDir: 'dist',
    rollupOptions: { input: 'index.html' },
    target: 'esnext',
    assetsInlineLimit: 100_000_000,
    chunkSizeWarningLimit: 100_000,
    sourcemap: false,
  },
  // the variant worker (src/standalone/variantWorker.ts) is inlined in the page as a classic script, the widest
  // browser support; its libraries' dynamic imports are folded into that one script, which cannot load chunks
  worker: { format: 'iife', rollupOptions: { output: { inlineDynamicImports: true } } },
  esbuild: { drop: ['debugger'] },
})
