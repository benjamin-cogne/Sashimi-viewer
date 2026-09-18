# Columnar storage of reads and coverage (exported pages, version 2)

An exported page (`Export HTML`) carries, for every registered view and sample, the coverage,
the junctions and, when the reads track was on, the reads of the window. Since version 2 of the
export these are stored in a compact binary form instead of JSON: about 24 times smaller for
reads, 5 times for coverage, and decoded only when a view is shown. This document is the
specification, written so that another program (the command-line converter of a whole BAM, a
reader in another language) can produce or read the same bytes. The reference implementation is
`src/standalone/columnar.ts`; a round-trip test lives next to the simulation scripts.

The page payload itself stays a JSON object in a `<script id="sashimi-embedded"
type="application/json">` tag (`EmbeddedExport` in `src/standalone/embedded.ts`). In version 2,
`views[i].coverage[sampleId]` and `views[i].reads[sampleId]` hold `{ "bin": "<base64>", ... }`
objects; version 1 pages held the JSON forms and are still read.

## 1. Primitives

- **Byte order**: little-endian for the one fixed-width integer (the directory length).
- **Unsigned varint** (`u`): 7 bits per byte, least significant group first, high bit set on
  every byte but the last. Values are non-negative integers below 2^53.
- **Signed varint** (`s`): the value is folded onto the unsigned integers first
  (`zigzag`: 0, -1, 1, -2, 2 … become 0, 1, 2, 3, 4 …), then written as `u`.
- **String**: `u` byte length, then UTF-8 bytes.
- **Byte**: one raw byte.
- **Compression**: raw deflate (RFC 1951, no zlib or gzip header), what the browser's
  `CompressionStream('deflate-raw')` writes and `DecompressionStream('deflate-raw')` reads.
- **Base64**: standard alphabet with padding, when a stream is written into a text file.

## 2. Stream layout

```
u32 little-endian   length L of the directory
L bytes             directory, JSON (UTF-8)
section 1 bytes     as many bytes as the directory says
section 2 bytes
…
```

Directory:

```json
{ "v": 1, "kind": "reads" | "coverage",
  "sections": [ { "name": "core", "block": 0, "start": 5001, "end": 9222, "n": 1000, "bytes": 4381 }, … ],
  "meta": { … } }
```

- `v` is the format version of this document (1). A reader refuses a higher version.
- `sections` lists the compressed sections in stream order with their compressed byte length.
  A reader that does not know a section name skips it by its length. This is how fields are
  added later without breaking older readers.
- `meta` holds facts of the whole stream, listed per kind below.

## 3. Reads stream (`kind: "reads"`)

`meta`: `window {start, end}` (0-based half-open window the reads were taken from), `total`
(reads passing the filters in the window before any cap), `reads` (reads in the stream),
`reference_source` (`"fasta"`, `"ensembl"`, `"browser"` or null).

Reads are sorted by start, then end, and cut into **blocks** of at most 1,000 reads. Each block
has a `core` section followed by its companion sections, each present only when the block needs
it: `pairs` (any read has a mate), `clips` (soft-clipped bases, hard-clipped lengths), `inserts`
(inserted bases) and `sa` (SA tags of split reads). A reader takes every section that follows a
`core` section with the same block index, in any order, and skips the names it does not know. Block sections carry `block` (index), `start` and `end` (genomic span covered by the block's
reads, 0-based half-open) and `n` (reads in the block), so a reader can decode only the blocks
overlapping a window. Read names are not stored: a reader numbers the reads `read 1`,
`read 2`, … in stream order.

### 3.1 `core` section (inflated content)

All positions are 0-based. "Relative" means relative to the read's own start.

```
u   n                         reads in the block
n × u   start delta           start minus the previous read's start (the first: minus 0)
n × u   aligned length        end - start on the reference
n × u   flags                 SAM flags
n × u   MAPQ
n × s   NH                    NH tag, -1 when absent
n × ( u k, k × (u gap, u length) )
                              aligned blocks after the first: gap from the previous block's end
                              (a splice gap or a deletion), then the block's length
n × u   first block length
n × ( u k, k × (u offset, u length) )   deletions (CIGAR D): relative start, length
n × ( u k, k × (u offset, u length) )   insertions (CIGAR I): relative position, inserted bases
n × ( u k, k × (u delta, byte base, byte quality) )
                              mismatches: position as delta from the previous mismatch (the
                              first: from the read start), base code (A 0, C 1, G 2, T 3, N 4),
                              base quality
n × ( u left, u right )       soft-clipped bases at each end
```

The read's strand is bit 0x10 of the flags. A read with one aligned block has `k = 0` in the
"blocks after the first" column.

### 3.2 `pairs` section

```
u   c                         number of distinct mate chromosomes other than the read's own
c × string                    their names
n × s   code                  -1 no mate stored, 0 mate on the same chromosome, k+1 = names[k]
n × s   mate start delta      mate start minus the read start (0 when no mate)
n × s   template length       TLEN as the aligner set it (0 when no mate)
```

### 3.3 `clips` section (per block, present when any read of the block has a soft or hard clip)

```
n × ( string left, string right, u hard left, u hard right )
```

The soft-clipped bases at each end of the alignment, as the record stores them (reference
strand), empty when the sequence was not available; then the hard-clipped lengths (bases the
record does not carry: they sit in the read's primary record).

### 3.4 `inserts` section (per block, present when any read of the block has an inserted sequence)

```
for each read, for each insertion of the core section, in order:   string bases
```

### 3.5 `sa` section (per block, present when any read of the block is split)

```
n × string     the SA tag as the aligner wrote it ("rname,pos,strand,CIGAR,mapQ,NM;" per part), empty otherwise
```

### 3.6 `reference` section

```
u   start                     0-based start of the reference window
string                        bases, upper case
```

## 4. Coverage stream (`kind: "coverage"`)

Two sections, `runs` and `junctions`; `meta` is empty. The rest of a coverage answer (window,
boundary-spanning counts, sampling facts, structural evidence, error) stays in the JSON next to
the stream.

`runs` (run-length depth, 0-based half-open runs that abut):

```
u   n           runs
u   start       start of the first run
n × u  length
n × s  depth delta   depth minus the previous run's depth (the first: minus 0)
```

`junctions` (introns of spliced reads, sorted by start then end):

```
u   n
n × u  start delta   start minus the previous junction's start
n × u  length        end - start
n × u  count         reads
```

## 5. Sizes and timings (5,000 simulated paired 150-bp reads, 20 kb window)

| | bytes |
|---|---|
| reads as JSON (version 1) | 742,763 |
| reads stream with pairs and reference, before base64 | 40,140 |
| the same in the page (base64) | 53,520 |
| coverage as JSON (7,060 runs) | 35,508 |
| coverage stream in the page | 6,692 |

Encoding took 52 ms and decoding 36 ms in Node with the streaming codec; the page decodes each
view the first time it is shown, off the main thread, the active view first and the others in
idle moments.

## 6. Extending the format

Add a section with a new name; never change the meaning of an existing one. A field that
concerns every read (names, base qualities, full sequence, tags) becomes a per-block section
next to `core` and `pairs`, with one entry per read in block order, so the reader can decode it
only when a feature needs it. Bump `v` only when a reader of the previous version would
misread an existing section.
