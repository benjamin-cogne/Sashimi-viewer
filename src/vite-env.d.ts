/// <reference types="vite/client" />

interface ImportMetaEnv {
  /** "1" in the unreleased build published under /dev/ (set by .github/workflows/pages.yml). */
  readonly VITE_DEV_MODE?: string;
  readonly VITE_DEV_BRANCH?: string;
  readonly VITE_DEV_SHA?: string;
  readonly VITE_DEV_DATE?: string;
  /** "1" shows the structural-variant hints on a build without the DEV MODE banner (they otherwise need ?sv=1). */
  readonly VITE_SV_HINTS?: string;
  /** URL of the stable build, linked from the DEV MODE banner. */
  readonly VITE_STABLE_URL?: string;
}

/** version from package.json, injected by the build */
declare const __APP_VERSION__: string;
