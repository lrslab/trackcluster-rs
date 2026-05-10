# CLI reference

This is a reference for the Rust binaries shipped in this repo:
- `trackcluster` (main CLI)
- `clusterj_batch` (gene-batched runner)

## Common conventions
- Inputs are BED12 / bigGenePred-like.
- Coordinates are 0-based, half-open.
- Most commands are strand-aware (the legacy TrackCluster behavior).
- **Gene name field**: stored in extra field index `5` (value like `GENE1` or `GENE1||GENE2`; `none` means unassigned).
- **Subread payload (`name2`)**: stored in extra field index `0` for isoforms produced by `clusterj`/`cluster`.
  - Default (`--name2-mode coverage`): `|3.5` (no read IDs)
  - `--name2-mode full`: `read1,read2,|3.5`
  - `--name2-mode none`: `none` (no payload)
  - Use the `*_read_to_isoform.tsv` mapping for counting when read IDs are not embedded (and it is auto-discovered by `count` / `count-multi` when present next to the isoform BED).

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
2) internal per-gene clustering (`clusterj` by default, or overlap-mode `cluster` with `--cluster-mode cluster`)
3) merge per-gene outputs into `<prefix>_isoform.bed` and `<prefix>_unused.bed`
4) `count` and `desc` on the merged isoforms
5) when `--manifest` is used: run per-sample `count-multi` outputs from the pooled isoforms

Key flags:
- `--cluster-mode`: `clusterj` (default) or `cluster` (overlap/intersection mode)
- `--reads/-s`: reads BED (single-sample mode; mutually exclusive with `--manifest`)
- `--manifest`: sample manifest TSV for pooled clustering + per-sample usage
- `--reference/-r`: reference BED
- `--output-root/-o`: output directory (created if missing)
- `--prefix`: output prefix (used for merged outputs like `<prefix>_isoform.bed`)
- `--threads/-t`: number of worker threads (parallel across genes)
- `--sw-score`: Smith-Waterman cutoff for 5' truncation collapsing (default: `11`; set to `-1` to disable). In overlap mode, the second pass only collapses a short read when `score < --sw-score`; a read at the exact cutoff remains its own track.
- `--batch-size`, `--batch-rounds`: bounds for very large genes; in overlap mode these control iterative pre-merging rounds before the final two-pass overlap clustering
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size; mapping TSVs are used for counting)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction and SL 5' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `2`): internal junction-site correction controls used by junction-mode clustering.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): junction-mode SL 5' merge controls
- `--overlap-cutoff1`, `--overlap-cutoff2`, `--overlap-intron-weight`: overlap-mode controls used when `--cluster-mode cluster`
- `--prepare-fraction-read`, `--prepare-fraction-ref`: overlap thresholds for gene assignment
- `--emit-pooled-reads`: when using `--manifest`, also write `<prefix>_pooled_reads.bed`
- `--max-reads-per-gene` (default: `50000`; set `0` to disable), `--downsample-gene`, `--downsample-seed`: per-gene downsampling (writes `clusterj_batch_downsample.tsv` or `cluster_batch_downsample.tsv` and scales counts)
- `--force`: overwrite existing per-gene outputs (otherwise genes with existing outputs are skipped)
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)

Outputs (under `--output-root`):
- `<prefix>_isoform.bed`
- `<prefix>_unused.bed`
- `<prefix>_read_to_isoform.tsv`
- `<prefix>_isoform_count.csv`
- `<prefix>_desc.txt`, `<prefix>_class4.txt`, `<prefix>_class12.txt`, `<prefix>_fusion.txt`
- `clusterj_batch_summary.txt` / `cluster_batch_summary.txt` and matching `*_errors.txt`
- `clusterj_batch_downsample.tsv` / `cluster_batch_downsample.tsv` (only when downsampling occurs)
- `<prefix>_pooled_reads.bed` (manifest mode + `--emit-pooled-reads`)
- `<prefix>.isoform_usage.long.tsv` (manifest mode only)
- `<prefix>.isoform_counts.matrix.tsv` (manifest mode only)
- `<prefix>.isoform_usage.group.tsv` (manifest mode only; only when manifest has `group`)

#### Platform presets and SL/no-SL data

`--platform-preset` changes only the junction correction and SL 5' merge/protection defaults. The `--sw-score` default remains `11` for all presets.

| Preset | Recommended use | Junction correction offset | Junction min support | SL partial 5' offset | SL same-junction 5' offset | SL 5' cluster offset | SL min support |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `generic` | Default; conservative setting for mixed or unknown platforms | `10` | `2` | `15` | `25` | `15` | `2` |
| `rna002` | RNA002 direct RNA reads; more tolerant of junction and 5' end wobble | `15` | `2` | `20` | `25` | `20` | `2` |
| `rna004` | RNA004 direct RNA reads; use the conservative/default cutoffs | `10` | `2` | `15` | `25` | `15` | `2` |

SL information is optional. Reads with no SL evidence, or reads whose BED score is not greater than `--sw-score`, are treated as non-SL-supported reads: they still go through junction correction and normal 5' truncation collapsing, but they do not receive SL-cluster protection as alternative 5' isoforms. For datasets where many reads have no SL information, keep the default `--sw-score 11` and use `--platform-preset rna004` for RNA004 or `--platform-preset rna002` for RNA002. Use `--sw-score -1` only when you want to disable 5' truncation collapsing entirely.

Example:
```bash
trackcluster flow \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --threads 8
```

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

Manifest example:
```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

Manifest-mode example:
```bash
trackcluster flow \
  --manifest samples.tsv \
  --reference ref.bed \
  --output-root out \
  --prefix pooled
```

### `trackcluster preparedir`
Split one reads BED + one reference BED into per-gene folders.

When to use:
- For manual gene-batched mode before running `clusterj_batch`.
- Overlap-mode batching is not exposed as a separate `cluster_batch` binary; use `trackcluster flow --cluster-mode cluster` for that path.

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
If `--out isoform.bed`, this command writes:
- `isoform.bed`: isoform BED
- `isoform.read_to_isoform.tsv`: read -> isoform mapping
- `isoform.unused.bed`: rare/filtered reads

Key flags:
- `--reads/-s`, `--reference/-r`, `--out/-o`
- `--threads/-t`: number of worker threads
- `--sw-score`: Smith-Waterman cutoff for 5' truncation collapsing (default: `11`; set to `-1` to disable)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction and SL 5' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `2`): internal junction-site correction controls. This offset is separate from the SL/5' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls

Performance note:
- 5' truncation collapsing is implemented via a junction-suffix index to avoid quadratic scans on large loci; `--batch-size`/`--batch-rounds` remain useful as hard caps for extreme genes.
- SL-supported reads with enough nearby 5' support can be protected as alternative isoforms; singleton likely-degradation reads can still merge into compatible longer/reference tracks.

Example:
```bash
trackcluster clusterj \
  --reads reads.bed \
  --reference ref.bed \
  --out isoform.bed \
  --threads 16
```

### `trackcluster cluster`
Overlap-based clustering (slower, more permissive; intended to mimic the legacy two-round TrackCluster overlap/intersection path).

Outputs are the same as `clusterj`:
- isoform BED, mapping TSV, unused BED (derived from the `--out` prefix).

Key flags:
- `--reads/-s`, `--reference/-r`, `--out/-o`
- `--threads/-t`: number of worker threads
- `--batch-size`, `--batch-rounds`: optional overlap batching for large loci (`--batch-size 0` disables intermediate batching)
- `--sw-score`: Smith-Waterman cutoff for 5' truncation collapsing in pass 2 (default: `11`; set to `-1` to disable). In pass 2, a short read is only collapsed when `score < --sw-score`; `score == --sw-score` remains a separate track.
- `--cutoff1`, `--cutoff2`: overlap pass 1 / pass 2 cutoffs (default: `0.05`, `0.01`)
- `--intron-weight`: intron contribution to the combined overlap distance (default: `0.5`)
- `--name2-mode`: `coverage` (default), `full`, or `none`

Behavior summary:
- Pass 1 uses the `ratio` distance (`1 - overlap / union_len`) with `--cutoff1`.
- Pass 2 uses the `ratio_short` distance (`1 - overlap / min_len`) with `--cutoff2`.
- Exon and intron distances are combined with `--intron-weight`.

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

Tip:
- A mapping TSV is written by `clusterj`/`cluster`/`flow`. `count` will auto-discover it when it lives next to the isoform BED (recommended; required when `name2` does not embed read IDs, e.g. default `--name2-mode coverage`).

Output:
- CSV with header `isoform_id,count`

Example:
```bash
trackcluster count \
  --reads reads.bed \
  --reference ref.bed \
  --isoform isoform.bed \
  --read-to-isoform isoform.read_to_isoform.tsv \
  --out isoform_count.csv
```

### `trackcluster count-multi`
Compute per-sample isoform counts/proportions from pooled isoforms using a sample manifest.

Input:
- `--manifest`: TSV with required columns `sample`, `reads`; optional `group`
- `--reference/-r`: reference BED
- `--isoform/-i`: pooled isoform BED (typically from `flow --manifest` or pooled `clusterj`)
- `--read-to-isoform`: optional mapping TSV (recommended; required when isoform `name2` does not embed read IDs; auto-discovered when next to the isoform BED)
- `--out/-o`: output prefix

Outputs (`--out <prefix>`):
- `<prefix>.isoform_usage.long.tsv`
- `<prefix>.isoform_counts.matrix.tsv`
- `<prefix>.isoform_usage.group.tsv` (only when manifest includes `group`)

Long-table semantics:
- One row per `(gene, isoform, sample)` with non-zero count.
- `proportion` is within-gene usage for that `(gene, sample)`:
  `proportion = count / sum(count of all isoforms for the gene+sample)`.

Matrix-table semantics:
- Rows are `(gene, isoform)`.
- Columns are samples in manifest order.
- Missing isoforms in a sample appear as `0`.

Example:
```bash
trackcluster count-multi \
  --manifest samples.tsv \
  --reference ref.bed \
  --isoform pooled_isoform.bed \
  --read-to-isoform pooled_read_to_isoform.tsv \
  --out out/pooled
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

Scope:
- This binary is junction-mode only.
- For overlap-mode batched clustering, use `trackcluster flow --cluster-mode cluster` or `trackcluster cluster` directly on a single reads/reference pair.

Useful flags:
- `--sw-score`: Smith-Waterman cutoff for 5' truncation collapsing (default: `11`; set to `-1` to disable)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction and SL 5' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `2`): internal junction-site correction controls. This offset is separate from the SL/5' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)
- `--max-reads-per-gene` (default: `50000`; set `0` to disable), `--downsample-gene`, `--downsample-seed`: per-gene downsampling (writes `clusterj_batch_downsample.tsv`)

Typical usage (after `preparedir`):
```bash
clusterj_batch \
  --input-root tracktest \
  --output-root trackout \
  --threads 8 \
  --force
```

`clusterj_batch` writes the active platform preset, junction correction settings, and SL settings into `clusterj_batch_summary.txt` for reproducibility.

One-command convenience (prepare + batch):
```bash
clusterj_batch \
  --prepare-reads reads.bed \
  --prepare-reference ref.bed \
  --prepare-prefix sample \
  --input-root tracktest \
  --output-root trackout \
  --threads 8 \
  --force
```
