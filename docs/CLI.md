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
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (default: `-1`, no SW/SL 5' signal). Pass a non-negative cutoff such as `11` only when BED score is valid SL/SW 5' evidence. In overlap mode, the second pass protects a short read only when its score is at or above the cutoff; with `-1`, ordinary short-read merging still runs.
- `--batch-size`, `--batch-rounds`: bounds for very large genes; in overlap mode these control iterative pre-merging rounds before the final two-pass overlap clustering
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size; mapping TSVs are used for counting)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls used by junction-mode clustering.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): junction-mode SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): junction-mode same-junction 3' terminal retention controls.
- Junction-mode clustering retains supported same-junction 3' terminal clusters as isoforms; this is strand-aware, so minus-strand 3' early stops are recognized at the lower genomic coordinate.
- `--overlap-cutoff1`, `--overlap-cutoff2`, `--overlap-intron-weight`: overlap-mode controls used when `--cluster-mode cluster`
- `--prepare-fraction-read`, `--prepare-fraction-ref`: overlap thresholds for gene assignment
- `--assignment-mode`: final counting mode, `unique` (default) or `fractional`; unique mode expands candidates against the isoform catalog before choosing the closest compatible isoform, including retained 3' early-stop isoforms.
- `--emit-pooled-reads`: when using `--manifest`, also write `<prefix>_pooled_reads.bed`
- `--max-reads-per-gene` (default: `50000`; set `0` to disable), `--downsample-gene`, `--downsample-seed`: per-gene downsampling (writes `clusterj_batch_downsample.tsv` or `cluster_batch_downsample.tsv` and scales counts)
- `--force`: overwrite existing per-gene outputs (otherwise genes with existing outputs are skipped)
- `--count-only`: reuse existing per-gene outputs and run only merge, unique/fractional count, multi-sample count, and desc outputs
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)

Outputs (under `--output-root`):
- `<prefix>_isoform.bed`
- `<prefix>_unused.bed`
- `<prefix>_read_to_isoform.tsv`
- `<prefix>_read_to_isoform.unique.tsv` (unique assignment mode only; exact selected mapping used for final counts)
- `<prefix>_isoform_count.csv`
- `<prefix>_desc.txt`, `<prefix>_class4.txt`, `<prefix>_class12.txt`, `<prefix>_fusion.txt`
- `clusterj_batch_summary.txt` / `cluster_batch_summary.txt` and matching `*_errors.txt`
- `clusterj_batch_downsample.tsv` / `cluster_batch_downsample.tsv` (only when downsampling occurs)
- `<prefix>_pooled_reads.bed` (manifest mode + `--emit-pooled-reads`)
- `<prefix>.isoform_count.csv` (manifest mode only; aggregate counts derived from the per-sample matrix)
- `<prefix>.isoform_usage.long.tsv` (manifest mode only)
- `<prefix>.isoform_counts.matrix.tsv` (manifest mode only)
- `<prefix>.isoform_usage.group.tsv` (manifest mode only; only when manifest has `group`)

#### Platform presets and SL/no-SL data

`--platform-preset` changes only the junction correction, SL 5', and same-junction 3' merge/protection defaults. The junction-mode `--sw-score` default remains `-1` for all presets; pass `--sw-score 11` when the BED score should be used as valid SL/SW 5' evidence.
Junction min support is weighted site support: read sites contribute `1`, reference sites contribute `5`.

| Preset | Recommended use | Junction correction offset | Junction min support | SL partial 5' offset | SL same-junction 5' offset | SL 5' cluster offset | SL min support | 3' same-junction offset | 3' cluster offset | 3' min support |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `generic` | Default; conservative setting for mixed or unknown platforms | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |
| `rna002` | RNA002 direct RNA reads; more tolerant of junction and 5' end wobble | `15` | `5` | `20` | `25` | `20` | `2` | `50` | `15` | `5` |
| `rna004` | RNA004 direct RNA reads; use the conservative/default cutoffs | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |

SL information is optional. With the default `--sw-score -1`, all reads are treated as non-SL-supported and ordinary junction correction/truncation merging still runs. When `--sw-score` is non-negative, only reads whose BED score is greater than the cutoff receive SL-cluster protection as alternative 5' isoforms. For no-SL datasets, including human data where BED score is MAPQ or another non-SL value, keep `--sw-score -1`. For SL/SW-scored datasets, pass an explicit cutoff such as `--sw-score 11` together with the appropriate platform preset.

Same-junction 3' terminal clusters with nearby read support are retained as isoforms independently of SL evidence. By default, protection requires at least `5` reads within the active junction correction offset and a 3' end more than `50` bp from the merge target. The rule is strand-aware: on minus-strand transcripts, a 3' early stop appears as a higher `tx_start`/lower genomic terminal boundary relative to the full-length isoform.

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

Count-only rerun example:
```bash
trackcluster flow \
  --count-only \
  --reference ref.bed \
  --output-root out \
  --prefix sample
```

`--count-only` expects completed per-gene output folders under `--output-root`. It reads `<prefix>_gene.txt` when present, otherwise it discovers gene folders in the output root. In unique assignment mode, it selects reads directly from each gene folder using `{gene}_nano.bed`, the per-gene isoform BED, and `{gene}_read_to_isoform.tsv`; missing per-gene count inputs are skipped. Add `--manifest samples.tsv` to a count-only rerun when you need manifest-mode `*.isoform_usage.*` outputs regenerated.

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
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (default: `-1`, no SW/SL 5' signal; pass a non-negative cutoff such as `11` only when BED score is valid SL/SW 5' evidence)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls. This offset is separate from the SL/5' and 3' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): same-junction 3' terminal retention controls

Performance note:
- 5' truncation collapsing is implemented via a junction-suffix index to avoid quadratic scans on large loci; `--batch-size`/`--batch-rounds` remain useful as hard caps for extreme genes.
- SL-supported reads with enough nearby 5' support can be protected as alternative isoforms; singleton likely-degradation reads can still merge into compatible longer/reference tracks.
- Supported same-junction 3' terminal clusters are retained as isoforms and remain compatible with unique counting, which assigns reads to the closest terminal structure.

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
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection in pass 2 (default: `11`; set to `-1` to treat reads as having no SW 5' signal). In pass 2, a short read is protected only when its score is at or above the cutoff; with `-1`, ordinary short-read merging still runs.
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
Compute isoform counts from an existing `flow`/`clusterj-batch` output directory.
This is the recommended mode because unique assignment is rerun inside each
gene folder using `{gene}_nano.bed`, the per-gene isoform BED, and
`{gene}_read_to_isoform.tsv` before merged counts are written.
When using manual `preparedir` -> `clusterj_batch`, use the same directory for
`clusterj_batch --input-root` and `--output-root` if you want to recount with
`trackcluster count --output-root`; a separate cluster output directory does not
contain `{gene}_nano.bed` unless you copy the prepared inputs there.

Input:
- `--output-root/-o`: existing output directory containing per-gene folders
- `--prefix`: prefix for merged outputs
- `--reference/-r`: reference BED
- `--cluster-mode`: `clusterj` (default) or `cluster`
- `--assignment-mode`: `unique` (default) or `fractional`

Output:
- `<prefix>_isoform.bed`
- `<prefix>_read_to_isoform.tsv`
- `<prefix>_read_to_isoform.unique.tsv` in unique mode
- `<prefix>_isoform_count.csv`
- description files under `<prefix>_*`

Example:
```bash
trackcluster count \
  --reference ref.bed \
  --output-root out \
  --prefix sample
```

Legacy low-level mode can still count a standalone isoform BED. In this mode,
the selector sees the supplied BED/mapping as its full scope, so prefer
`--output-root` when continuing a per-gene cluster run.

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
- `--assignment-mode`: `unique` (default; expand candidates against the isoform catalog and assign each read to the closest compatible isoform using read/isoform structure) or `fractional` (split multi-mapped reads across mapped candidates)
- `--out/-o`: output prefix

Outputs (`--out <prefix>`):
- `<prefix>.isoform_count.csv` (aggregate counts derived from the per-sample matrix)
- `<prefix>.isoform_usage.long.tsv`
- `<prefix>.isoform_counts.matrix.tsv`
- `<prefix>.isoform_usage.group.tsv` (only when manifest includes `group`)

Aggregate count semantics:
- `count` is exactly the sum of the sample columns in `<prefix>.isoform_counts.matrix.tsv` for the same isoform.

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
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (default: `-1`, no SW/SL 5' signal; pass a non-negative cutoff such as `11` only when BED score is valid SL/SW 5' evidence)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls. This offset is separate from the SL/5' and 3' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): same-junction 3' terminal retention controls
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)
- `--max-reads-per-gene` (default: `50000`; set `0` to disable), `--downsample-gene`, `--downsample-seed`: per-gene downsampling (writes `clusterj_batch_downsample.tsv`)

Typical usage (after `preparedir`):
```bash
clusterj_batch \
  --input-root tracktest \
  --output-root tracktest \
  --threads 8 \
  --force
```

`clusterj_batch` writes the active platform preset, junction correction settings, SL settings, and same-junction 3' settings into `clusterj_batch_summary.txt` for reproducibility.

One-command convenience (prepare + batch):
```bash
clusterj_batch \
  --prepare-reads reads.bed \
  --prepare-reference ref.bed \
  --prepare-prefix sample \
  --input-root tracktest \
  --output-root tracktest \
  --threads 8 \
  --force
```

Use a distinct `--output-root` only when you intend to consume the per-gene
cluster outputs directly or concatenate them manually. For
`trackcluster count --output-root` unique assignment, the count root must also
contain the prepared `{gene}_nano.bed` files.
