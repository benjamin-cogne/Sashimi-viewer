import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { viteSingleFile } from 'vite-plugin-singlefile'

// Builds the viewer as ONE self-contained HTML file: dist/index.html, then copied by
// scripts/finish.mjs to dist/sashimi-viewer.html and to ./sashimi-viewer.html (tracked in git).
export default defineConfig({
  plugins: [react(), viteSingleFile()],
  server: { port: 3000 },
  build: {
    outDir: 'dist',
    rollupOptions: { input: 'index.html' },
    target: 'esnext',
    assetsInlineLimit: 100_000_000,
    chunkSizeWarningLimit: 100_000,
    sourcemap: false,
  },
  esbuild: { drop: ['debugger'] },
})
