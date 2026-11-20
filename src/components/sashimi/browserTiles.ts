/**
 * Data of a genome browser by tiles: the chromosome is cut into tiles of a size chosen from the span in view (1, 4, 16
 * or 64 Mb), each fetched once and kept while it is among the most recently used (a bounded cache, so a long walk
 * along a chromosome does not keep everything it saw). While a tile loads, the tiles of another size already loaded
 * over the same stretch stand in for it, so zooming never shows an empty track.
 */

export const TILE_SIZES = [1_000_000, 4_000_000, 16_000_000, 64_000_000] as const;
/** The tile size for a span in view: the smallest tile holding it (the view then touches two or three tiles). */
export const tileSizeFor = (span: number): number => TILE_SIZES.find(s => s >= span) ?? TILE_SIZES[TILE_SIZES.length - 1];

interface Entry<T> { size: number; i: number; items: T[] | null; error: string | null; used: number }

export class TileCache<T> {
  private entries = new Map<string, Entry<T>>();
  private tick = 0;
  /** bumped whenever a tile arrives or fails: what a React memo can depend on */
  version = 0;
  constructor(
    private fetcher: (start: number, end: number, size: number) => Promise<T[]>,
    private length: number,
    private cap: number,
    private onChange: () => void,
  ) {}

  private key = (size: number, i: number) => `${size}:${i}`;

  /** The tiles of `size` over [start, end) and one more on each side (fetched ahead for panning). */
  want(start: number, end: number, size: number): void {
    const i0 = Math.max(0, Math.floor(start / size) - 1), i1 = Math.min(Math.floor((this.length - 1) / size), Math.floor(end / size) + 1);
    for (let i = i0; i <= i1; i++) {
      const k = this.key(size, i);
      const e = this.entries.get(k);
      if (e) { e.used = ++this.tick; continue; }
      const entry: Entry<T> = { size, i, items: null, error: null, used: ++this.tick };
      this.entries.set(k, entry);
      this.fetcher(i * size, Math.min(this.length, (i + 1) * size), size)
        .then(items => { entry.items = items; }, err => { entry.error = String(err?.message || err); })
        .finally(() => { if (this.entries.get(k) === entry) { this.version++; this.evict(); this.onChange(); } });
    }
  }

  /** Drop the least recently used tiles past the cap (never one still loading). */
  private evict(): void {
    if (this.entries.size <= this.cap) return;
    const done = [...this.entries.entries()].filter(([, e]) => e.items || e.error).sort((a, b) => a[1].used - b[1].used);
    for (const [k] of done.slice(0, this.entries.size - this.cap)) this.entries.delete(k);
  }

  /**
   * The items over [start, end) at `size`, a tile not loaded yet replaced by the loaded tiles of other sizes over it.
   * `loading`: some of the stretch is still on its way; `error`: the last failure met.
   */
  collect(start: number, end: number, size: number): { items: T[]; loading: boolean; error: string | null } {
    const out: T[] = [];
    let loading = false, error: string | null = null;
    const last = Math.floor((this.length - 1) / size);
    for (let i = Math.max(0, Math.floor(start / size)); i * size < end && i <= last; i++) {
      const e = this.entries.get(this.key(size, i));
      if (e?.items) { out.push(...e.items); continue; }
      if (e?.error) error = e.error; else loading = true;
      const a = i * size, b = (i + 1) * size;
      for (const o of this.entries.values()) if (o.items && o.size !== size && o.i * o.size < b && (o.i + 1) * o.size > a) out.push(...o.items);
    }
    return { items: out, loading, error };
  }
}
