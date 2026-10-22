/**
 * Loads every plugin of this folder: each src/plugins/<name>/index.ts, whose top level registers what it adds (a file
 * kind, fileKinds.ts). The page and the variant worker both import this module, so a plugin is present in both. With
 * no plugin folder, as in this repository, the glob matches nothing and the build is unchanged.
 */
import.meta.glob('./*/index.ts', { eager: true });

export {};
