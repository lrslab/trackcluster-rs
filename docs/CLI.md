# CLI reference

This is a reference for the Rust binaries shipped in this repo:
- `trackcluster` (main CLI)
- `clusterj_batch` (gene-batched runner)

## Common conventions
- Inputs are BED12 / bigGenePred-like.
- Coordinates are 0-based, half-open.
- Most commands are strand-aware (the legacy TrackCluster behavior).
- **Gene name field**: stored in extra field index `5` (value like `GENE1` or `GENE1||GENE2`; `none` means unassigned).
- **Subread list (`name2`)**: stored in extra field index `0` for isoforms produced by `clusterj`/`cluster` (value like `read1,read2,|3.5`).

## `trackcluster`

### `trackcluster validate-bed`
Validate that a BED12/bigGenePred parses cleanly and print basic counts.

Example:
```bash
trackcluster validate-bed --input reads.bed
```

### `trackcluster flow`
Run the full pipeline as a single command:
1) `preparedir` (dedup + gene assignment into per-gene folders)
2) `clusterj_batch` (cluster per gene in parallel)
3) merge per-gene outputs into `<prefix>_isoform.bed` and `<prefix>_unused.bed`
4) `count` and `desc` on the merged isoforms

Key flags:
- `--reads/-s`: reads BED
- `--reference/-r`: reference BED
- `--output-root/-o`: output directory (created if missing)
- `--prefix`: output prefix (used for merged outputs like `<prefix>_isoform.bed`)
- `--threads/-t`: number of worker threads (parallel across genes)
- `--sw-score`: set to `-1` to disable 5' truncation collapsing
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--prepare-fraction-read`, `--prepare-fraction-ref`: overlap thresholds for gene assignment
- `--force`: overwrite existing per-gene outputs (otherwise genes with existing outputs are skipped)

Outputs (under `--output-root`):
- `<prefix>_isoform.bed`
- `<prefix>_unused.bed`
- `<prefix>_isoform_count.csv`
- `<prefix>_desc.txt`, `<prefix>_class4.txt`, `<prefix>_class12.txt`, `<prefix>_fusion.txt`
- `clusterj_batch_summary.txt`, `clusterj_batch_errors.txt`

Example:
```bash
trackcluster flow \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --threads 32 \
  --sw-score -1
```

### `trackcluster preparedir`
Split one reads BED + one reference BED into per-gene folders.

When to use:
- For manual gene-batched mode before running `clusterj_batch`.

Key flags:
- `--reads/-s`: reads BED
- `--reference/-r`: reference BED
- `--output-root/-o`: directory to create (gene folders written here)
- `--prefix`: prefix for summary outputs
- `--fraction-read`, `--fraction-ref`: overlap thresholds for gene assignment

Example:
```bash
trackcluster preparedir \
  --reads reads.bed \
  --reference ref.bed \
  --output-root tracktest \
  --prefix sample
```

### `trackcluster addgene`
Assign gene names to reads by overlap with reference and write `reads_gene.bed`.

Notes:
- Uses default overlap thresholds (same as `preparedir` defaults).
- If you need custom thresholds, use `preparedir` (it exposes them via flags).

Example:
```bash
trackcluster addgene \
  --reads reads.bed \
  --reference ref.bed \
  --out reads_gene.bed
```

### `trackcluster clusterj`
Junction-chain clustering (fast mode).

Outputs:
- `<out>`: isoform BED
- `<out>.read_to_isoform.tsv`: read -> isoform mapping
- `<out>.unused.bed`: rare/filtered reads

Key flags:
- `--reads/-s`, `--reference/-r`, `--out/-o`
- `--threads/-t`: number of worker threads
- `--sw-score`: set to `-1` to disable 5' truncation collapsing
- `--batch-size`, `--batch-rounds`: bounds for very large genes

Example:
```bash
trackcluster clusterj \
  --reads reads.bed \
  --reference ref.bed \
  --out isoform.bed \
  --threads 16 \
  --sw-score -1
```

### `trackcluster cluster`
Overlap-based clustering (slower, more permissive).

Outputs are the same as `clusterj`:
- isoform BED, mapping TSV, unused BED (derived from the `--out` prefix).

Example:
```bash
trackcluster cluster \
  --reads reads.bed \
  --reference ref.bed \
  --out isoform.bed \
  --threads 8
```

### `trackcluster count`
Compute isoform counts from the isoform BED produced by `clusterj`/`cluster`.

Output:
- CSV with header `isoform_id,count`

Example:
```bash
trackcluster count \
  --reads reads.bed \
  --reference ref.bed \
  --isoform isoform.bed \
  --out isoform_count.csv
```

### `trackcluster desc`
Describe/classify isoforms vs reference, and detect fusions.

Outputs (`--out <prefix>`):
- `<prefix>_desc.txt`
- `<prefix>_class4.txt`
- `<prefix>_class12.txt`
- `<prefix>_fusion.txt`

Important:
- `desc` groups isoforms by gene name field. If your isoforms do not have gene names, run `addgene` (or use `flow` / the `preparedir` -> `clusterj_batch` workflow).

Example:
```bash
trackcluster desc \
  --isoform isoform.bed \
  --reference ref.bed \
  --out sample
```

## `clusterj_batch`
Run `clusterj` per gene folder in parallel (optionally run `preparedir` first).

Typical usage (after `preparedir`):
```bash
clusterj_batch \
  --input-root tracktest \
  --output-root trackout \
  --threads 32 \
  --sw-score -1 \
  --force
```

One-command convenience (prepare + batch):
```bash
clusterj_batch \
  --prepare-reads reads.bed \
  --prepare-reference ref.bed \
  --prepare-prefix sample \
  --input-root tracktest \
  --output-root trackout \
  --threads 32 \
  --sw-score -1 \
  --force
```
