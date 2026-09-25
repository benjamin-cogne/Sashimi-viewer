# Columnar storage of reads and coverage (exported pages)

An exported page (`Export HTML`) carries, for every registered view and sample, the coverage,
the junctions and, when the reads track was on, the reads of the window. Since version 2 of the
export these are stored in a compact binary form instead of JSON: about 24 times smaller for
reads, 5 times for coverage, and decoded only when a view is shown. This document is the
specification, written so that another program (a writer of whole files, a reader in another
language) can produce or read the same bytes. The reference implementation is
`src/standalone/columnar.ts`; a round-trip test lives next to the simulation scripts.

The page payload itself stays a JSON object in a `<script id="sashimi-embedded"
type="application/json">` tag (`EmbeddedExport` in `src/standalone/embedded.ts`). In version 2,
`views[i].coverage[sampleId]` and `views[i].reads[sampleId]` hold `{ "bin": "<base64>", ... }`
objects; version 1 pages held the JSON forms and are still read. Version 3 of the export writes
the streams of columnar version 4 (below), version 4 of the export those of columnar version 5, and
both read every earlier one.

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
L bytes             directory: varints from version 4, JSON (UTF-8) up to version 3
section 1 bytes     as many bytes as the directory says
section 2 bytes
…
```

The directory says what the stream holds. Its content, whatever the encoding:

```json
{ "v": 5, "kind": "reads" | "coverage",
  "sections": [ { "name": "core", "block": 0, "start": 5001, "end": 9222, "n": 1000, "bytes": 4381 }, … ],
  "meta": { … } }
```

From version 4 it is written as varints (§1). A reader tells the two apart by the first byte: a
JSON directory starts with `{` (0x7b), a binary one with its version, which is 4 or more.

```
u   v
u   kind                       0 reads, 1 coverage
u   number of sections
per section:
  u   name code                0 core, 1 pairs, 2 clips, 3 inserts, 4 sa, 5 mods, 6 reference,
                               7 runs, 8 junctions, 9 methylation, 11 hap, 12 seq;
                               10 = escape: a string with the name follows
  [string name]                only for code 10
  u   fields                   bit 1 block, 2 start, 4 end, 8 n: which of them follow
  u   bytes
  [u block]
  [s start]                    minus the start of the previous section that has one (0 at first)
  [s end]                      minus this section's start
  [u n]
u   meta fields                bit 1 window, 2 total, 4 reads, 8 reference_source
  [u window start, s window end - window start]
  [u total] [u reads]
  [string reference_source]    empty for null
string  the other meta keys as JSON, empty when there are none
```

The codes are append-only: a code never changes meaning, and the escape stays 10 however long
the list grows, so a section a reader does not know still travels under code 10 with its name.
A reader that meets a code it does not know skips that section by its length, like an unknown
name. `mods` (5), `methylation` (9) and `seq` (12, whole read sequences) are sections written by other
producers of these streams; this viewer neither writes nor reads them. JSON cost about 226 bytes a stream, the varints about
25: most of a small window's stream, and it adds up over the many windows of a page.

- `v` is the format version of this document (5). A reader refuses a higher version, and must
  decode each section the way the version it reads asks for:
  - **5** stores a read's alignment shape as a stream of CIGAR operations (§3.1) and a mate link
    once for both reads when they are in the same block (§3.2);
  - **4** writes the directory as varints (above); the sections are unchanged;
  - **3** links a read's mate to another read of the stream and stores the template length as a
    residual (§3.2);
  - **2** changed only the `mods` section of other producers; the sections of this document
    read as in version 1;
  - **1** is the original.
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
(inserted bases), `sa` (SA tags of split reads) and `hap` (haplotags of a phased file). A reader takes every section that follows a
`core` section with the same block index, in any order, and skips the names it does not know. Block sections carry `block` (index), `start` and `end` (genomic span covered by the block's
reads, 0-based half-open) and `n` (reads in the block), so a reader can decode only the blocks
overlapping a window. Read names are not stored: a reader numbers the reads `read 1`,
`read 2`, … in stream order.

### 3.1 `core` section (inflated content)

All positions are 0-based. "Relative" means relative to the read's own start. Version 5:

```
u   n                         reads in the block
n × u   start delta           start minus the previous read's start (the first: minus 0)
n × u   flags                 SAM flags
n × u   MAPQ
n × s   NH                    NH tag, -1 when absent
n × u   op count              CIGAR operations of the read, clips left out
K × byte   op code            0 M (=, X included), 1 I, 2 D, 3 N; K the sum of the counts, read by read
M × u   M length              the length of each M op, in op order
K-M × u   other length        the length of each I, D and N op, in op order
n × ( u k, k × (u delta, byte base, byte quality) )
                              mismatches: position as delta from the previous mismatch (the
                              first: from the read start), base code (A 0, C 1, G 2, T 3, N 4),
                              base quality
n × ( u left, u right )       soft-clipped bases at each end
```

The ops rebuild the read as a BAM record would: from the read start, an M op is an aligned block,
an I op an insertion before the current position, a D op a deletion and an N op skipped bases (an
intron), each of the last two moving the position on; the read ends where the ops end. A read
with no op (an unmapped read placed next to its mate) has no block and ends where it starts. The
strand is bit 0x10 of the flags.

The shape used to be stored in parts (versions 1–4, below): the end, the first block, each
further block as gap and length, then each deletion and each insertion again with its offset
from the read start. For a long read, whose hundreds of indels each split a block, that was three
numbers per indel where one op does, and offsets that grow along the read: 1,726 bytes per ONT
read against 562 in the four op columns (measured on 12,075 reads of a 30× simulated ONT locus,
each block of 1,000 reads deflated). A paired RNA-seq file goes from 1.3 to 0.81 bytes per read
for the shape, uniform short reads from 0.70 to 0.66. Keeping codes, M lengths and the other
lengths in columns of their own is what makes deflate find them: interleaved as `length × 4 + op`
the same ops cost 662 bytes per ONT read.

Versions 1–4:

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
n × ( u k, k × (u delta, byte base, byte quality) )   mismatches, as in version 5
n × ( u left, u right )       soft-clipped bases at each end
```

### 3.2 `pairs` section

Version 5:

```
u   c                         number of distinct mate chromosomes other than the read's own
c × string                    their names
n' × u  link                  0 none, else zigzag(signed distance to the mate's read in this stream) × 2
                              + 1 when the link is mutual inside the block (below); the reads a mutual
                              link already points back to have no entry
p × s   code                  for the reads with no link: -1 no mate stored, 0 mate on the same
                              chromosome, k+1 = names[k]
m × s   mate start delta      mate start minus the read start, for those with code ≥ 0
k × s   template length       for the paired reads (linked, or code ≥ 0): TLEN minus the derived
                              span when linked, else TLEN
```

A read with a link has its mate in the stream, so on its chromosome: it needs no code. When a
read links to a later read of the same block that links back to it, which is the case of nearly
every pair, the low bit says so and the later read's link is taken as the way back; it writes
neither link nor code. That halves the link column, the bulk of the section (1.20 bytes per read
of paired RNA-seq). A link into another block is written by both reads, so that each block
decodes without the others.

Versions 3–4 write every read's code, then every read's link as a signed varint (0 none);
versions 1–2 have no link column and write the mate start and template length of every paired
read.

Version 3 links a read to its mate instead of repeating where it is. The link is the signed
distance, in reads of this stream, from the read to one whose start equals the mate start; the
stream's reads are numbered from 0 across all its blocks, so `block × reads_per_block + i`, and
a link may point into another block. Where there is a link, the mate position is read back from
that read and the template length is stored as the residual against the span the aligner would
have written,

```
derived TLEN = (max(read end, mate end) − min(read start, mate start)) × (read start ≤ mate start ? 1 : −1)
```

which is zero for all but a handful of reads (3 in 245,613 measured). A read with no link keeps
the version 1–2 columns: the mate start as a delta and the template length verbatim. That is the
escape for a mate outside the stream's window, on another chromosome, or simply not written —
about 2% of the reads of a paired file cut into 100 kb windows.

For the position and the template length, the link only has to point at a read whose start
equals the mate start: that reproduces the mate position exactly and the residual absorbs the
rest. But a link is also how a reader tells which read *is* the mate, so the writer links the
true one when it is in the stream — the read with the same name, the other pair bit (64/128) and
a primary record — and only otherwise the read whose residual vanishes, then the nearest.

A reader takes two reads **linked to each other**, one of each pair bit, as a pair: the viewer
gives both the same mate key (`mk`, the stream's window start and the lower of the two indices,
so it is the same every time a stream is decoded and never shared across streams) and joins them
on it instead of by position. A one-way link is not a pairing and is left to the positional rule.
This matters where fragments share their starts: position alone then joins the first read it
finds, and a stream sorts reads by start *and end*, so that is often another fragment's mate
— 37 % of the pairs on a fixture made of such look-alikes, against none with the links. A pair
split across two streams has no link and is still joined by position (`src/standalone/mates.ts`:
same name first, then the opposite template length).

Measured on 300,000 paired 71 bp reads, the section falls from 3.54 to 1.35 bytes per read (−62%).
A reader that filters blocks out has to inflate the `core` of a block it is skipping when a link
points into it; `core` is the cheapest section to inflate, and a reader decoding a whole stream —
the usual case — never pays for it at all.

### 3.3 `clips` section (per block, present when any read of the block has a soft or hard clip)

```
n × ( string left, string right, u hard left, u hard right )
```

The soft-clipped bases at each end of the alignment, as the record stores them (reference
strand), empty when the sequence was not available; then the hard-clipped lengths (bases the
record does not carry: they sit in the read's primary record).

A string may be shorter than the clip's length in `core`: a writer may keep only the bases next
to the alignment, which are the *end* of a left clip and the *start* of a right one. A reader
places them from the alignment outwards, and treats the rest of the clip as bases not stored (the
viewer draws that part as a plain bar). The viewer's page export stores clips whole.

### 3.4 `inserts` section (per block, present when any read of the block has an inserted sequence)

```
for each read, for each insertion of the core section, in order:   string bases
```

### 3.5 `sa` section (per block, present when any read of the block is split)

```
n × string     the SA tag as the aligner wrote it ("rname,pos,strand,CIGAR,mapQ,NM;" per part), empty otherwise
```

### 3.6 `hap` section (per block, present when any read of the block carries a haplotag)

```
per read:  varint  HP     haplotype (1, 2, …) written by the phasing tool; 0 = untagged, nothing follows
           varint  PS+1   phase set + 1, 0 when the record has no PS tag        (tagged reads only)
           varint  PC+1   assignment confidence (Phred) + 1, 0 when absent       (tagged reads only)
```

### 3.7 `reference` section

```
u   start                     0-based start of the reference window
string                        bases, upper case
```

## 4. Coverage stream (`kind: "coverage"`)

Two sections, `runs` and `junctions`; `meta` is empty.
The rest of a coverage answer (window, boundary-spanning counts, sampling facts, structural
evidence, error) stays in the JSON next to the stream.

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

Version 3's mate links take about 8 % off a window of this shape, the `pairs` sections going from
1.25 to 0.96 bytes per read. That is the small end of the range: in a simulation every pair has
the same shape, so the mate deltas and template lengths are drawn from a narrow distribution and
deflate had already squeezed them to 1.25 bytes. Real mates vary, and there the links matter
much more — 3.54 to 1.35 bytes per read on a real paired RNA-seq file.

Encoding took 52 ms and decoding 36 ms in Node with the streaming codec; the page decodes each
view the first time it is shown, off the main thread, the active view first and the others in
idle moments.

## 6. Extending the format

Add a section with a new name, and give it the next code in the list (13); never change the
meaning of an existing one. A field that
concerns every read (names, base qualities, full sequence, tags) becomes a per-block section
next to `core` and `pairs`, with one entry per read in block order, so the reader can decode it
only when a feature needs it. Bump `v` only when a reader of the previous version would
misread an existing section.
