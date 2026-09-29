/// <reference types="vite/client" />

interface ImportMetaEnv {
  /** "1" in the unreleased build published under /dev/ (set by .github/workflows/pages.yml). */
  readonly VITE_DEV_MODE?: string;
  readonly VITE_DEV_BRANCH?: string;
  readonly VITE_DEV_SHA?: string;
  readonly VITE_DEV_DATE?: string;
  /** URL of the stable build, linked from the DEV MODE banner. */
  readonly VITE_STABLE_URL?: string;
}

/** version from package.json, injected by the build */
declare const __APP_VERSION__: string;
