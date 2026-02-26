# File formats

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

### Gene name field (TrackCluster convention)
This Rust rewrite follows the legacy TrackCluster convention of storing a gene name in an extra field:
- Gene name lives at **extra field index 5** (0-based within `extra_fields`).
- Unassigned is `none`.
- Multi-gene values are joined with `||` (example: `GENE1||GENE2`).

### Subread list / `name2` field (TrackCluster convention)
Isoforms produced by `clusterj`/`cluster` use extra field index `0` for the TrackCluster-style `name2` payload.

By default (`--name2-mode coverage`), this stores only a coverage value:
- `extra_fields[0] = "|<coverage>"` (no read IDs)
- `--name2-mode none`: `extra_fields[0] = "none"` (no payload)
- `--name2-mode full`: `extra_fields[0] = "<read1>,<read2>,...,|<coverage>"`
  - Example: `readA,readB,readC,|2.5`

When read IDs are omitted from `name2`, use the `*_read_to_isoform.tsv` mapping written by `clusterj`/`cluster`/`flow` and pass it to `count` / `count-multi` via `--read-to-isoform` (or keep it next to the isoform BED for auto-discovery).

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
Only emitted when manifest has `group`.

Columns:
- `gene`
- `isoform_id`
- `group`
- `count`
- `proportion`
- `gene_total`
