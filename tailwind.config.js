/** @type {import('tailwindcss').Config} */
// The viewer's theme ("Ardoise"): slate-blue neutrals, a teal accent, a dark navy app bar. Gray and slate share one
// neutral scale; `accent` is the colour of buttons, active toggles and links; `bar` the app bar and its controls.
const neutral = {
  50: '#f5f7fa', 100: '#eceff4', 200: '#d3dbe5', 300: '#bcc6d3', 400: '#8796a8', 500: '#5b6b7e',
  600: '#475569', 700: '#334155', 800: '#1e293b', 900: '#0f1b2b', 950: '#0a1220',
}
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        gray: neutral,
        slate: neutral,
        accent: {
          50: '#ddf1ee', 100: '#cbe9e4', 200: '#b0ddd6', 300: '#93d2c9', 400: '#4fb3a6', 500: '#14857a',
          600: '#0f766e', 700: '#0c625b', 800: '#0a4f4a', 900: '#083d3a', 950: '#052b29',
        },
        bar: { DEFAULT: '#13253b', ctrl: '#1f3249', line: '#34475f', text: '#f1f5f9', muted: '#a3b5c9' },
        page: '#ebe8e3',
      },
      borderRadius: { lg: '6px', xl: '9px' },
    },
  },
  plugins: [],
}
