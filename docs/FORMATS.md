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
Isoforms produced by `clusterj`/`cluster` store supporting read IDs in the first extra field:
- `extra_fields[0] = "<read1>,<read2>,...,|<coverage>"`
- Example: `readA,readB,readC,|2.5`

The `count` command parses the `<read1>,<read2>,...` portion and uses it to compute isoform counts.
