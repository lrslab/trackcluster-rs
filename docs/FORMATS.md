# File formats

See [`INTERCHANGE.md`](INTERCHANGE.md) for BAM/GFF3/GTF conversion and transcript exports.

## 0.2.0 output migration

Version 0.2.0 deliberately changes several
previously ambiguous output contracts. Count CSV is now RFC 4180 with the
header `gene,isoform_id,count`; novel transcript IDs use the `tc_novel_v1:`
structural namespace; and full `name2` payloads use the percent-encoded
`tc_name2_v1:` codec. Description files move from the legacy headerless layout
to `trackcluster-description-v2`. Consumers that assumed two-column count CSV,
representative-read novel IDs, raw comma-split `name2`, or headerless
description files must migrate using the contracts below. Readers retain
support for legacy unescaped `name2` payloads.

Description files use schema `trackcluster-description-v2`. Each starts with a
`#schema` line followed by an explicit column header, replacing the legacy
headerless files. The Figure 2 category values, UTR threshold, and class overlap
priority are defined in [`behavior/desc.md`](behavior/desc.md).

- `_desc.txt`: `isoform_id`, `reference_id`, `gene_id`, `missing_features`,
  `extra_features`
- `_class4.txt`: `isoform_id`, `class`
- `_fusion.txt`: `isoform_id`, `gene_ids`
- `_class12.txt`: `isoform_id`, `class`

The Rust rewrite primarily works with BED12 and bigGenePred-like files.

## BED coordinate conventions
- Coordinates are **0-based, half-open**: `[start, end)`.
- Two intervals overlap iff `a.start < b.end && b.start < a.end`.

## Required fields
At minimum, inputs should include BED12 columns:
- `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`
- `thickStart`, `thickEnd`, `itemRgb`
- `blockCount`, `blockSizes`, `blockStarts`

## bigGenePred extras
Many TrackCluster datasets use bigGenePred extra columns (e.g., gene name / annotation metadata).
The Rust parser preserves extra trailing fields and writes them back out when present.

### Converter-produced BED12+8

`bam2bigg` and `gff2bigg` write tab-separated, bigGenePred-compatible
**BED12+8 text**. They do not create the indexed binary `.bb` representation.
Both use BED's zero-based, half-open coordinates and emit all eight TrackCluster
extension columns in this order:

| Extra index | bigGenePred name | `bam2bigg` | `gff2bigg` |
| ---: | --- | --- | --- |
| 0 | `name2` | `none` | `none` |
| 1 | `cdsStartStat` | `none` | `none` |
| 2 | `cdsEndStat` | `none` | `none` |
| 3 | `exonFrames` | one `-1` per block | one `-1` per block |
| 4 | `type` | `nanopore_read` | `isoform_anno` |
| 5 | `geneName` / TrackCluster gene ID | `none` | resolved annotation gene(s), otherwise `none` |
| 6 | `geneName2` / TrackCluster sample group | `--group` or BAM filename stem | `none` |
| 7 | `geneType` / TrackCluster reserved | `none` | `none` |

The standard BED fields differ as follows:

- `bam2bigg` uses the BAM query name, MAPQ as BED score, flag-derived strand,
  strand color (`250,128,114` for `+`, `64,224,208` for `-`), and CIGAR-derived
  blocks. CIGAR `N` alone splits blocks. `thickStart` and `thickEnd` are both
  zero.
- `gff2bigg` uses the annotation transcript identity, score `100`, resolved
  exon strand, `itemRgb=0`, and exon-derived blocks. GFF/GTF coordinates
  `[start,end]` become BED `[start-1,end)`. CDS/UTR/phase are not transferred,
  so `thickStart` and `thickEnd` are both zero.

`gff2bigg` sorts output deterministically by transcript structure. `bam2bigg`
preserves retained BAM record order and preserves separate alignment instances,
including repeated query names.

### Gene name field (TrackCluster convention)
This Rust rewrite follows the legacy TrackCluster convention of storing a gene name in an extra field:
- Gene name lives at **extra field index 5** (0-based within `extra_fields`).
- Unassigned is `none`.
- Multi-gene values are joined with `||` (example: `GENE1||GENE2`).

### Biological gene IDs and filesystem keys

Per-gene paths never interpolate an unchecked biological gene ID. Gene IDs must be non-empty,
single path components and cannot be absolute paths, `.`/`..`, contain `/` or `\\`, or contain
NUL/control characters. The maximum accepted length is **4096 UTF-8 bytes**; longer IDs fail
validation before any per-gene path is constructed.

TrackCluster derives a deterministic filesystem key for each valid ID. ASCII letters, digits, `_`,
`-`, and `.` keep their historical spelling (for example, `GENE-1`); Unicode and other punctuation
are percent-encoded when the result fits, and oversized path keys use a fixed-size stable hash.
Per-gene directory and artifact filenames use this key. `<prefix>_gene_paths.tsv` and the batch-level
`clusterj_batch_gene_paths.tsv`/`cluster_batch_gene_paths.tsv` are versioned mappings from the
biological ID to its key. Each gene directory also contains `.trackcluster_gene_id`, so metadata
validation can round-trip hashed IDs without treating arbitrary directories as selected genes.
Before any preparation or batch publication, TrackCluster rejects a key that
equals one of its prefix-scoped merged/preparation filenames or fixed batch
report filenames. This keeps gene directories structurally disjoint from every
top-level pipeline artifact; changing the output prefix resolves a
prefix-derived collision, while a collision with a fixed
`clusterj_batch_*`/`cluster_batch_*` name requires a different biological gene
ID. These run-scoped namespace restrictions are checked in addition to the
context-free `GeneId` syntax and path-key rules above.

Gene lists and selectors stay in the biological namespace:
`<prefix>_gene.txt`, `clusterj_batch --gene-list`, and `--downsample-gene`
contain or accept biological gene IDs, not encoded keys or directory names.
The prefix-scoped gene list is also the preparation commit marker: preparation
publishes it last. An empty list means a replacement failed after publication
began, so the prepared generation is incomplete and must not be clustered or
counted until preparation succeeds again.

### Subread list / `name2` field (TrackCluster convention)
Isoforms produced by `clusterj`/`cluster` use extra field index `0` for the TrackCluster-style `name2` payload.

By default (`--name2-mode coverage`), this stores only a coverage value:
- `extra_fields[0] = "|<coverage>"` (no read IDs)
- `--name2-mode none`: `extra_fields[0] = "none"` (no payload)
- `--name2-mode full`: `extra_fields[0] = "tc_name2_v1:<encoded-read1>,<encoded-read2>,...,|<coverage>"`
  - Example: `tc_name2_v1:readA,readB,readC,|2.5`
  - Read IDs are UTF-8 percent encoded. Commas, pipes, percent signs, whitespace,
    and other reserved bytes therefore round-trip without being mistaken for
    payload delimiters.

Readers still accept the pre-v1 unescaped comma-separated payload for backward
compatibility. That legacy form cannot represent a comma inside a read ID;
re-emit the catalog or use the mapping TSV to migrate such data. All newly
written full payloads use `tc_name2_v1`.

When read IDs are omitted from `name2`, use the `*_read_to_isoform.tsv` mapping written by `clusterj`/`cluster`/`flow`. For ordinary recounts, prefer `trackcluster count --output-root <out> --prefix <prefix>`; it reads each gene folder directly, so unique assignment and retained-intron checks stay gene-local. The legacy standalone BED mode can still take `--read-to-isoform` (or auto-discover it next to the isoform BED), but its unique assignment scope is the supplied merged input. Counting defaults to unique best assignment: it expands compatible candidates against the isoform catalog, then selects the closest isoform per read before counting. Pass `--assignment-mode fractional` for compatibility with split multi-mapped counts from the mapping file.

Mapping files contain exactly two raw TSV fields, `read_id` and `isoform_id`,
without a header. Leading and trailing spaces are significant identity bytes and
round-trip unchanged; empty fields, tabs within a field, and embedded line
breaks are invalid.

In `flow` unique assignment mode, `<prefix>_read_to_isoform.tsv` remains the raw merged mapping from per-gene clustering. The selected mapping actually used for final counts is written separately as `<prefix>_read_to_isoform.unique.tsv`; use that file when auditing or reproducing unique-mode counts. `<prefix>_unique_assignment.provenance.tsv` records the effective junction tolerance, ordered one-to-one matcher, and explicit no-collapse policy for microfeatures.

## Stable identity contract

- Reference isoform IDs are preserved. Empty IDs, duplicate reference IDs, and
  references that claim the reserved `tc_novel_v1:` namespace are rejected.
- Novel isoforms use `tc_novel_v1:<gene-hex>:<chromosome-hex>:<strand>:<exons>`.
  This is a lossless structural namespace, not a representative read name or a
  truncated hash, so distinct `(gene, chromosome, strand, exon-chain)` tuples
  cannot silently collide.
- Catalog IDs are validated globally before a merged catalog is atomically
  published and again at counting boundaries.
- A read ID is an abundance **molecule ID**. Multiple BED rows with that ID are
  alignment instances during clustering, but repeated identical mapping rows
  are idempotent during counting. Fractional counting splits one molecule over
  its distinct candidate isoforms. Unique assignment accepts identical duplicate
  structures but rejects conflicting alignments for one molecule ID because
  there is no unambiguous structure to score.
- A read ID equal to a reference transcript ID remains a read; source provenance
  is structural and is never inferred from matching strings.

### Rejected-read diagnostics

`flow`, `preparedir`, `clusterj_batch`, `clusterj`, and `cluster` use `--invalid-read-policy skip` by default. This recovery
is deliberately record-local: a malformed read BED row or a parsed read row whose ID is empty is
excluded, while the remaining valid read tracks continue. Preparation publishes
`<prefix>_rejected_reads.tsv`; per-gene clustering publishes `rejected_reads.tsv` in each gene
directory for read tracks rejected while that gene is loaded. Header-only files mean that no read
was rejected at that boundary.

The direct single-gene commands derive the same diagnostic from `--out`: for example,
`--out isoform.bed` writes `isoform.rejected_reads.tsv`.

Both files have this stable header:

```tsv
source_path	line	read_id	kind	reason
```

- `source_path` and the 1-based physical `line` locate the input record.
- `read_id` is populated when an identifier is available; it is empty when parsing failed before a
  usable ID could be recovered.
- `kind` is `parse` for malformed BED and `identity` for an empty read ID.
- `reason` contains the parser or identity diagnostic. Backslashes, tabs, carriage returns, and
  newlines inside text fields are escaped as `\\`, `\t`, `\r`, and `\n`.

`--invalid-read-policy fail` restores strict read-track parsing: the first malformed or empty-ID
read fails the enclosing preparation/gene stage. Read-file I/O errors, reference parse/identity
errors, invalid configuration, and algorithm/integrity failures are fatal at their existing
boundary under both policies and never appear in these TSVs.

The rejected-read TSV is an additive audit artifact. Existing input BED, isoform BED, mapping,
count, and description/classification schemas are unchanged, as are the clustering and isoform
classification rules. However, rejected reads provide no junction, terminal, coverage, mapping, or
counting evidence, so filtering them can change which isoforms are called and consequently the
contents (but not the schema or category definitions) of classification outputs.

### Per-gene `run.json`

`flow` and `clusterj_batch` publish schema-versioned JSON completion manifests in each gene
directory. Source-aware tool identity uses manifest schema version 3; older per-gene manifests are
rebuilt rather than reused. The manifest has `status: "complete"`, a request fingerprint, input SHA-256 hashes,
effective options (including invalid-read policy, assignment mode, and unique-assignment junction
tolerance), package version, Git commit, deterministic source fingerprint, the per-gene seed, and
output SHA-256 hashes, byte sizes, and record counts. Clean Git checkouts use
`clean`. Cargo package builds and dirty source checkouts use a SHA-256
fingerprint of the actual build-source snapshot; rebuilding after editing an
unpacked package therefore cannot retain the official package's cache
identity. The
manifest is a completion marker rather than a user-editable configuration file: changing it or any
recorded output makes that gene stale on the next normal run. `flow --count-only` and
`count --output-root` require this marker for every selected gene and revalidate its self-fingerprint,
gene/mode/tool identity, current prepared-input contents, and every output hash, size, and record
count before publishing merged results. Legacy folders without a valid manifest must first be
rebuilt by a normal flow run.

Artifact publication is atomic per file, using a synced temporary sibling and rename. The
per-gene manifest is invalidated before replacement and published last, so it is the only reusable
completion marker. A set of flat merged/count/description files is not a cross-file filesystem
transaction; after an interruption during their publish phase, rerun the final stage. Those derived
files are outputs only and are never accepted as evidence that per-gene work completed.

### Batch summary status and recoverable gene errors

Batch summaries use `status=complete` when every gene succeeds, `status=partial` when the default
continue policy excludes one or more failed genes but retains at least one verified result, and
`status=failed` for strict mode, infrastructure failure, or an all-gene failure. The fields
`gene_error_policy`, `mergeable_genes`, `excluded_failed_genes`, and `infrastructure_errors` make
that distinction machine-readable. Detailed gene diagnostics are written to
`*_errors.txt` only when errors are recorded; a clean run removes any stale
error report.

Only processed genes, hash-verified cache hits, and semantic-empty genes with a valid empty
completion manifest are mergeable. Failed genes, their old artifacts, and their old downsampling
scale factors are excluded. This recovery policy applies to normal per-gene execution; count-only
manifest/hash verification remains all-or-nothing.

## Multi-sample manifest TSV
`count-multi` and `flow --manifest` expect a tab-separated manifest with:

- Required columns:
  - `sample`: unique sample name
  - `reads`: BED path (absolute or relative to manifest file)
- Optional columns:
  - `group`: condition/group label for pseudo-bulk summaries

Example:
```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

Constraints:
- `sample` must be unique.
- `sample` cannot contain `::` (reserved delimiter for pooled read IDs).
- Missing files fail fast with an error.

## Pooled read IDs
When pooling manifest reads, trackcluster rewrites read IDs as:
`<sample>::<orig_read_id>`

This guarantees sample identity for downstream per-sample counting.

## `count-multi` outputs

### Aggregate count table (`*.isoform_count.csv`)
Columns:
- `gene`
- `isoform_id`
- `count`

Semantics:
- The file is RFC 4180-compatible CSV written with field escaping; commas and
  quotes in identifiers do not create extra columns.
- `count` is exactly the sum of the sample columns in the matching
  `*.isoform_counts.matrix.tsv`.
- In `flow --manifest`, the main `<prefix>_isoform_count.csv` is synchronized
  from the same aggregate table, so total counts and sample counts share one
  assignment result.

### Long table (`*.isoform_usage.long.tsv`)
Columns:
- `gene`
- `isoform_id`
- `sample`
- `group` (present when any sample has a group)
- `count`
- `proportion`
- `gene_total`

Semantics:
- One row per `(gene, isoform_id, sample)` with non-zero count.
- `proportion` is within-gene usage:
  `count / sum(counts for all isoforms in the same gene and sample)`.

### Matrix table (`*.isoform_counts.matrix.tsv`)
Columns:
- `gene`
- `isoform_id`
- one column per sample (manifest order)

Semantics:
- Missing isoforms in a sample are represented as `0`.

### Group table (`*.isoform_usage.group.tsv`)
Only emitted when at least one sample has a non-empty `group` value.

Columns:
- `gene`
- `isoform_id`
- `group`
- `count`
- `proportion`
- `gene_total`

### Unique-assignment provenance (`*.unique_assignment.provenance.tsv`)

Emitted only in unique assignment mode. The table records the effective
`unique_assignment_junction_offset` and the one-molecule/one-isoform matching
policy used to generate the count tables. Fractional mode removes a stale
unique-mode provenance file when the same output prefix is rerun successfully.
