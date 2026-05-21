# Pipeline tutorial (Rust)

This document walks through the **end-to-end** Rust workflow, step by step, with examples.

You can run the whole pipeline in three ways:

1. **One-command mode (recommended)**: run `trackcluster flow` (does `preparedir` + per-gene clustering + merge + `count` + `desc`).
2. **Single-file mode**: run `trackcluster clusterj|cluster` directly on `reads.bed` + `ref.bed` (good for small inputs).
3. **Manual gene-batched mode**: run `preparedir` -> `clusterj_batch` -> combine outputs -> `count` + `desc` (junction-mode batch runner).

For multi-sample/condition studies, `flow` and `count-multi` support pooled isoform discovery with per-sample usage outputs.

## Prerequisites

Build the binaries:

```bash
cargo build --release
```

You should have:

- `./target/release/trackcluster`
- `./target/release/clusterj_batch`

## Inputs

- **Reads**: BED12 / bigGenePred-like (0-based, half-open).
- **Reference**: BED12 / bigGenePred-like, ideally with gene names populated in the extra fields.

Notes:

- Most commands say "sorted recommended". The current implementation will work on unsorted inputs, but performance and determinism are best on sorted inputs.
- Gene names in this rewrite are stored in extra field index `5` (the 6th "extra" column beyond BED12). Multiple genes can be joined with `||`. Unassigned is `none`.

## Option A (recommended) - One-command mode (`trackcluster flow`)

### What this does

`trackcluster flow` runs the legacy-style end-to-end pipeline in one command:

1. prepare per-gene folders (dedup + gene assignment)
2. run per-gene clustering internally (`clusterj` by default, or overlap-mode `cluster` with `--cluster-mode cluster`)
3. merge per-gene isoforms/unused reads into `<prefix>_isoform.bed` and `<prefix>_unused.bed`
4. run `count` and `desc` on the merged isoforms

### Command

```bash
trackcluster flow \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --threads 8
```

Defaults:
- `--cluster-mode clusterj` (junction mode). Use `--cluster-mode cluster` for the legacy overlap/intersection path.
- `--sw-score 11` (collapse 5' truncations); use `--sw-score -1` to disable.
- `--name2-mode coverage` (memory-friendly). Use `--name2-mode full` to embed read IDs in the isoform BED (larger outputs); otherwise rely on `*_read_to_isoform.tsv` for counting.
- SL-supported junction-mode 5' merge controls default to `--sl-partial-5prime-offset 15`, `--sl-same-junction-5prime-offset 25`, `--sl-5prime-cluster-offset 15`, and `--sl-5prime-min-support 2`.
- `--max-reads-per-gene 50000` (memory-friendly cap; set `0` to disable). Counts/usage tables are scaled when downsampling occurs.

Overlap-mode example:

```bash
trackcluster flow \
  --cluster-mode cluster \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --sw-score 11 \
  --batch-size 500 \
  --batch-rounds 100
```

Notes for overlap mode:
- `flow` writes per-gene overlap isoforms as `*_simple_coverage.bed` under the gene folders, then merges them into `<prefix>_isoform.bed`.
- Run summaries switch to `cluster_batch_summary.txt`, `cluster_batch_errors.txt`, and `cluster_batch_downsample.tsv`.
- In the second overlap pass, a short read is collapsed only when `score < --sw-score`; `score == --sw-score` stays as its own track.

### Multi-sample pooled command

Manifest TSV (`samples.tsv`):

```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

Run pooled clustering once, then compute per-sample usage:

```bash
trackcluster flow \
  --manifest samples.tsv \
  --reference ref.bed \
  --output-root out \
  --prefix pooled
```

### Outputs

Under `--output-root`:

- `sample_gene.txt`, `sample_dedup.bed`, `sample_novel.bed`
- per-gene folders with per-gene clustering outputs
- `sample_isoform.bed`, `sample_unused.bed`
- `sample_read_to_isoform.tsv`
- `sample_isoform_count.csv`
- `sample_desc.txt`, `sample_class4.txt`, `sample_class12.txt`, `sample_fusion.txt`

In manifest mode:

- `pooled_pooled_reads.bed` (only if `--emit-pooled-reads`; sample-tagged IDs like `S1::read123`)
- `pooled.isoform_count.csv` (aggregate counts derived from the per-sample matrix)
- `pooled.isoform_usage.long.tsv`
- `pooled.isoform_counts.matrix.tsv`
- `pooled.isoform_usage.group.tsv` (if manifest includes `group`)

In manifest mode, `pooled_isoform_count.csv` is synchronized from `pooled.isoform_count.csv`. This keeps the main total-count output consistent with the per-sample matrix and usage tables.

The batch runner also writes:

- `clusterj_batch_summary.txt` / `cluster_batch_summary.txt`
- `clusterj_batch_errors.txt` / `cluster_batch_errors.txt`
- `clusterj_batch_downsample.tsv` / `cluster_batch_downsample.tsv` (only when downsampling occurs)

## Step 0 - Validate inputs (optional but recommended)

Before running the pipeline, sanity-check your files parse cleanly:

```bash
trackcluster validate-bed --input reads.bed
trackcluster validate-bed --input ref.bed
```

If parsing fails, fix the BED12/bigGenePred formatting first (see `docs/FORMATS.md`).

## Option B - Manual gene-batched mode (`clusterj_batch`, junction mode)

This manual batched workflow is for `clusterj`. There is currently no separate `cluster_batch` binary; use `trackcluster flow --cluster-mode cluster` when you want the overlap-mode batched path.

### Step 1 - Prepare per-gene folders (`preparedir`)

#### What this step does

`preparedir` creates a "run folder" that mirrors the legacy Python pipeline's layout:

- Deduplicate reads by read name (keeps the read with the longest total exon length).
- Assign each read to one or more genes by span overlap against reference transcripts (strand-aware).
- Split reads and reference records into one folder per gene.

This bounds the clustering problem to gene-sized chunks (huge speed/memory win).

#### Command

```bash
trackcluster preparedir \
  --reads reads.bed \
  --reference ref.bed \
  --output-root tracktest \
  --prefix sample \
  --fraction-read 0.01 \
  --fraction-ref 0.05
```

#### Outputs

After running, `tracktest/` contains:

- `tracktest/sample_dedup.bed`
- `tracktest/sample_gene.txt`
- `tracktest/sample_novel.bed`
- `tracktest/<GENE>/<GENE>_nano.bed`
- `tracktest/<GENE>/<GENE>_gff.bed`

### Step 2 - Cluster per gene in parallel (`clusterj_batch`)

#### What this step does

`clusterj_batch` iterates over gene folders in `--input-root` and runs junction-mode clustering for each gene.

Internally, 5' truncation collapsing uses a junction-suffix index to avoid quadratic scans on large loci.
Low-support splice-junction sites are corrected before merging with `--junction-correction-offset` and `--junction-correction-min-support`.
Supported alternative SL 5' clusters can be protected from merging with the separate `--sl-*` terminal controls exposed by `trackcluster clusterj`, `trackcluster flow`, and `clusterj_batch`.
Use `--platform-preset rna002` for RNA002 reads (junction offset `15`; SL 5' offsets `20/25/20`) and `--platform-preset rna004` for RNA004 reads (conservative defaults: junction offset `10`; SL 5' offsets `15/25/15`). Explicit junction correction and SL options override the preset.
SL evidence is optional; reads without SL information still use junction correction and normal 5' truncation collapsing, but they are not protected by the SL-cluster rules.
Counting defaults to `--assignment-mode unique`, so each read contributes to only the closest compatible catalog isoform. Pass `--assignment-mode fractional` for compatibility with split multi-mapped counts from the mapping file.

For each gene, it writes:

- `*_simple_coveragej.bed` (isoforms)
- `*_unused.bed` (rare/filtered reads)
- `*_read_to_isoform.tsv` (mapping)

It also writes run summaries in the output root (e.g. `clusterj_batch_summary.txt`).

#### Command (gene-batched)

```bash
clusterj_batch \
  --input-root tracktest \
  --output-root trackout \
  --threads 8 \
  --force
```

### Step 3 - Combine per-gene outputs into a single isoform file

`clusterj_batch` produces one isoform BED per gene. Many downstream tools expect a single isoform file, so you usually concatenate them:

```bash
find trackout -mindepth 2 -maxdepth 2 -name '*_simple_coveragej.bed' -print0 \
  | sort -z \
  | xargs -0 cat > sample_isoform.bed
```

Optional: combine unused reads too:

```bash
find trackout -mindepth 2 -maxdepth 2 -name '*_unused.bed' -print0 \
  | sort -z \
  | xargs -0 cat > sample_unused.bed
```

Recommended: combine read-to-isoform mappings (required for default `--name2-mode coverage|none` counting):

```bash
find trackout -mindepth 2 -maxdepth 2 -name '*_read_to_isoform.tsv' -print0 \
  | sort -z \
  | xargs -0 cat > sample_read_to_isoform.tsv
```

## Step 4 - Count isoforms (`count`)

```bash
trackcluster count \
  --reads reads.bed \
  --reference ref.bed \
  --isoform sample_isoform.bed \
  --read-to-isoform sample_read_to_isoform.tsv \
  --out sample_isoform_count.csv
```

## Step 4b - Multi-sample usage (`count-multi`)

Use when isoforms were called from pooled reads and you need per-sample usage.

```bash
trackcluster count-multi \
  --manifest samples.tsv \
  --reference ref.bed \
  --isoform sample_isoform.bed \
  --read-to-isoform sample_read_to_isoform.tsv \
  --out sample
```

This writes:

- `sample.isoform_count.csv` (aggregate counts derived from the per-sample matrix)
- `sample.isoform_usage.long.tsv`
- `sample.isoform_counts.matrix.tsv`
- `sample.isoform_usage.group.tsv` (if groups are defined)

The aggregate count for each isoform is the sum of the corresponding row in `sample.isoform_counts.matrix.tsv`.

## Step 5 - Describe/classify isoforms (`desc`)

Critical requirement: isoforms must have gene names for per-gene comparisons.

If you used `preparedir` + `clusterj_batch`, gene names are naturally assigned within the per-gene workflow.

```bash
trackcluster desc \
  --isoform sample_isoform.bed \
  --reference ref.bed \
  --out sample
```

## Option C - Single-file mode (no `preparedir`)

If your dataset is small, you can cluster without the gene-batched folder structure:

```bash
trackcluster clusterj --reads reads.bed --reference ref.bed --out isoform.bed --threads 8
trackcluster cluster --reads reads.bed --reference ref.bed --out isoform.bed --threads 8
trackcluster count --reads reads.bed --reference ref.bed --isoform isoform.bed --read-to-isoform isoform.read_to_isoform.tsv --out isoform_count.csv
```

If you also want `desc`, make sure isoforms have gene names first:

```bash
trackcluster addgene --reads reads.bed --reference ref.bed --out reads_gene.bed
trackcluster clusterj --reads reads_gene.bed --reference ref.bed --out isoform.bed
trackcluster cluster --reads reads_gene.bed --reference ref.bed --out isoform.bed
trackcluster desc --isoform isoform.bed --reference ref.bed --out sample
```

## Reference script

See `scripts/run_full_flow_rust.sh` for an example end-to-end driver.
