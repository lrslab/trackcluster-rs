# Pipeline tutorial (Rust)

This document walks through the **end-to-end** Rust workflow, step by step, with examples.

Use a dedicated output directory whose complete tree is owned by TrackCluster.
The reads, reference, sample manifest, and every reads path named by a manifest
must remain outside `--output-root`. `flow`, `preparedir`, and output-root
recounting verify this boundary before mutation and reject pre-existing
symlink or hard-link aliases in the managed tree.

You can run the whole pipeline in three ways:

1. **One-command mode (recommended)**: run `trackcluster flow` (does `preparedir` + per-gene clustering + merge + `count` + `desc`).
2. **Single-file mode**: run `trackcluster clusterj|cluster` directly on `reads.bed` + `ref.bed` (good for small inputs).
3. **Manual gene-batched mode**: run `preparedir` -> `clusterj_batch` -> combine outputs -> `count` + `desc` (junction-mode batch runner).

If per-gene clustering already completed, `trackcluster flow --count-only` can rerun just the final merge/count/description stage.

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
- Gene-selection inputs such as `clusterj_batch --gene-list` and
  `--downsample-gene` use biological gene IDs. Per-gene directories and files
  use encoded `<gene-path-key>` values recorded in `<prefix>_gene_paths.tsv`;
  never substitute a folder name for a biological ID.

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
- `--sw-score -1` (the `flow` default in both cluster modes; treat reads as
  having no SW/SL 5' signal while keeping ordinary merging). Pass
  `--sw-score 11` only when BED score is valid SL/SW 5' evidence. The direct
  `trackcluster cluster` command separately defaults to `11` for legacy
  overlap-mode compatibility.
- `--name2-mode coverage` (memory-friendly). Use `--name2-mode full` to embed read IDs in the isoform BED (larger outputs); otherwise rely on `*_read_to_isoform.tsv` for counting.
- SL-supported junction-mode 5' merge controls default to `--sl-partial-5prime-offset 15`, `--sl-same-junction-5prime-offset 25`, `--sl-5prime-cluster-offset 15`, and `--sl-5prime-min-support 2`.
- Supported same-junction 3' terminal clusters are retained as isoforms with `--same-junction-3prime-offset 50`, `--3prime-min-support 5`, and a default `--3prime-cluster-offset` equal to the active junction correction offset.
- `--max-reads-per-gene 5000` (runtime- and memory-bounded cap; set `0` to disable). This lower default protects high-expression, single-exon loci such as mitochondrial `cox1`, whose diverse read endpoints can otherwise drive near-quadratic merge work. With the default 500-read merge batch it admits at most ten input batches per gene. Counts/usage tables are scaled when downsampling occurs. If an actually downsampled gene shares a molecule ID with another selected gene, final flow/count-only processing fails before publishing merged abundance outputs because independent per-gene reservoirs cannot preserve cross-gene candidate probabilities. Disable the cap or exclude every affected gene from downsampling.
- `--invalid-read-policy skip`: exclude only an individual malformed read track or a track with an empty read ID and continue with the valid tracks. Use `--invalid-read-policy fail` for strict read parsing.

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
- The run summary switches to `cluster_batch_summary.txt`.
  `cluster_batch_errors.txt` is present only when errors are recorded, and
  `cluster_batch_downsample.tsv` only when downsampling occurs.
- In the second overlap pass, a short read is protected only when its score is at or above `--sw-score`; with `--sw-score -1`, ordinary short-read merging still runs.

Count-only rerun:

```bash
trackcluster flow \
  --count-only \
  --reference ref.bed \
  --output-root out \
  --prefix sample
```

Use this when the per-gene folders already contain completed clustering outputs and you only need to regenerate merged isoforms, counts, multi-sample usage, or desc tables. `--count-only` uses `<prefix>_gene.txt` and `<prefix>_gene_paths.tsv` when they exist. A standalone batch tree without prefix-scoped metadata falls back to its versioned `clusterj_batch_gene_paths.tsv` or `cluster_batch_gene_paths.tsv`; arbitrary directories are never treated as the active gene set. It verifies every selected gene's completion manifest, current prepared-input hashes, cluster mode/tool identity, and every recorded output hash/count before publishing a merged file. A missing or stale artifact fails the entire rerun with a request to rebuild the producing flow or batch run; genes are never silently skipped. In unique assignment mode, the selected mapping used for counts is written to `<prefix>_read_to_isoform.unique.tsv`.

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
- `sample_gene_paths.tsv` (biological gene ID to filesystem-key mapping)
- `sample_rejected_reads.tsv` (read tracks rejected during preparation)
- per-gene folders with per-gene clustering outputs
- `<gene-path-key>/rejected_reads.tsv` (read tracks rejected while loading that gene)
- `sample_isoform.bed`, `sample_unused.bed`
- `sample_read_to_isoform.tsv`
- `sample_read_to_isoform.unique.tsv` (unique assignment mode; selected mapping used for final counts)
- `sample_unique_assignment.provenance.tsv` (unique mode; effective junction tolerance and one-to-one/no-collapse matching policy)
- `sample_isoform_count.csv`
- `sample_desc.txt`, `sample_class4.txt`, `sample_class12.txt`, `sample_fusion.txt`

In manifest mode:

- `pooled_pooled_reads.bed` (only if `--emit-pooled-reads`; sample-tagged IDs like `S1::read123`)
- `pooled.isoform_count.csv` (aggregate counts derived from the per-sample matrix)
- `pooled.isoform_usage.long.tsv`
- `pooled.isoform_counts.matrix.tsv`
- `pooled.isoform_usage.group.tsv` (when at least one sample has a non-empty `group`)

In manifest mode, `pooled_isoform_count.csv` is synchronized from `pooled.isoform_count.csv`. This keeps the main total-count output consistent with the per-sample matrix and usage tables.

The batch runner also writes:

- `clusterj_batch_summary.txt` / `cluster_batch_summary.txt`
- `clusterj_batch_errors.txt` / `cluster_batch_errors.txt` (only when errors
  are recorded)
- `clusterj_batch_downsample.tsv` / `cluster_batch_downsample.tsv` (only when downsampling occurs)

By default, an error isolated to one gene is recorded in the matching `*_errors.txt`, that gene is
excluded, and the flow continues through merge/count/description using only complete verified gene
artifacts. The summary status is `partial`. Use `--strict-gene-errors` when any gene failure must
stop the flow before downstream publication. Invalid global inputs/options, output failures,
infrastructure failures, an all-gene failure, and count-only integrity failures remain fatal.

Read-track recovery happens one level below this gene policy. With the default
`--invalid-read-policy skip`, only the bad row is excluded when a read BED record is malformed or
has an empty ID; the gene continues with its usable reads, and the rejection is auditable in the
TSVs above. If every read in a gene is rejected, the gene publishes valid empty artifacts rather
than poisoning unrelated genes. Read-file I/O failures, malformed references, invalid options, and
algorithm/integrity failures still fail their enclosing stage and are not converted into rejected
reads. `--invalid-read-policy fail` restores strict read parsing; combine it with
`--strict-gene-errors` when any resulting gene failure must prevent downstream publication.

This feature adds diagnostic TSVs but does not change the BED, isoform, mapping, count, or
description/classification schemas, nor the clustering/classification algorithms. Because excluded
reads no longer contribute splice-junction, terminal, support, mapping, or abundance evidence, the
resulting isoform calls and classification rows can still differ from a run whose input reads are
all valid.

The direct single-gene `clusterj` and `cluster` commands use the same policy. With
`--out isoform.bed`, their diagnostic is `isoform.rejected_reads.tsv`.

### Optional isoform-level modification post-processing

Modification import is intentionally separate from isoform discovery. First
normalize each caller output. For a primary genome-aligned Dorado/modBAM, for
example:

```bash
trackcluster mod-import-dorado \
  --sample S1 \
  --assay-id dorado_rna004_m6a \
  --bam S1.aligned.bam \
  --mod-code A+a \
  --model-id rna004_m6a_model \
  --candidate-rule all-target-canonical-bases \
  --source-emission-threshold 0.05 \
  --out S1.mod
```

This writes `S1.mod.observations.tsv`, `S1.mod.assay.json`, and
`S1.mod.import_qc.tsv`. Build one modification manifest row per sample and
compatible assay; `coverage_bam` can be `NA` when exact base-level coverage is
not required:

```tsv
sample	assay_id	observations	assay_metadata	coverage_bam
S1	dorado_rna004_m6a	S1.mod.observations.tsv	S1.mod.assay.json	S1.aligned.bam
S2	dorado_rna004_m6a	S2.mod.observations.tsv	S2.mod.assay.json	S2.aligned.bam
```

Then attach modification aggregation to the manifest-mode flow. It runs only
after the final unique read-to-isoform assignment has been published:

```bash
trackcluster flow \
  --manifest samples.tsv \
  --reference ref.bed \
  --output-root out \
  --prefix pooled \
  --assignment-mode unique \
  --mod-manifest mod_samples.tsv \
  --mod-reference-fasta genome.fa \
  --mod-analysis-threshold dorado_rna004_m6a=0.5 \
  --mod-eligibility-profile strict
```

The optional stage completes a hash-verified generation under
`<prefix>.mod.generations/<run_id>/`, synchronizes the flat compatibility files
`<prefix>.mod_join_qc.tsv`, `<prefix>.mod_site_join_qc.tsv`,
`<prefix>.isoform_mod_sites.tsv`, and `<prefix>.isoform_mod_design.tsv`, then
publishes `<prefix>.mod.current.json` last. Pass `--mod-contrasts contrasts.tsv`
to also write `<prefix>.isoform_mod_contrasts.tsv`. A rerun invalidates the
current pointer before core outputs can change; if the modification stage then
fails, surviving flat files are stale and downstream site-summary/contrast
commands reject them. Historical generation directories remain for audit.
The same aggregation can be run standalone with `trackcluster mod-aggregate`
when clustering/counting is already complete.

V1 requires exact unique assignment and keeps caller/model/chemistry strata
separate. Missing or unknown observations are not silently classified as
unmodified, and effect-only contrasts report `p_value` and `q_value` as `NA`.
The default `exploratory` eligibility profile is intended for QC and method
development. The `strict` profile additionally requires exact BAM coverage,
FASTA validation, Dorado version/model/threshold provenance verified from one
coherent source `@PG` record, and configurable covering, callable,
candidate-rate, and callable-rate guardrails; its default minima are
20, 10, 0.8, and 0.8 respectively and require caller/model-specific
calibration.
See the [CLI reference](CLI.md#trackcluster-mod-import-dorado), the exact
[modification formats and denominator rules](FORMATS.md#isoform-level-modification-formats),
and the [validation strategy and scientific limitations](MODIFICATION_VALIDATION.md).

### Safe resume and per-gene manifests

Every successfully clustered gene directory contains `run.json`. This versioned completion
manifest records SHA-256 hashes of the prepared reads and reference, all effective clustering and
assignment options, the deterministic per-gene seed, package version, Git commit, source
fingerprint, and hashes, byte sizes, and record counts for each per-gene output. Clean Git builds
use the `clean` marker; Cargo-packaged and dirty-checkout builds hash their actual source tree, so
each edited snapshot has its own SHA-256 identity. A normal rerun reuses a gene only when the
request fingerprint matches exactly and every recorded output still verifies. Existing filenames
without a valid manifest, changed inputs/options/tool versions, and missing or modified outputs are
treated as stale and rebuilt. `--force` always rebuilds and cannot be combined
with `--count-only`, which never executes the per-gene rebuild stage.

Per-gene output files are written and synced under temporary sibling names, then atomically
renamed. `run.json` is invalidated before rebuilding and published last, so an interrupted run
cannot advertise partially replaced files as complete. The three merged flow outputs use the same
temporary-write and atomic-rename protocol. Count/description and multi-sample files are all staged
successfully before their publish phase. Atomicity is per destination file: ordinary filesystems do
not provide a single transaction spanning the full flat output set. If interruption occurs during a
multi-file publish, rerun the command; only a valid per-gene `run.json` is a reusable completion
marker, and derived merged/count/description files are never treated as resumable state. Batch
summaries include one
`resume_decision\t<biological-gene-id>\t<action>\t<reason>` row per gene for
auditability.
The summary also records `gene_error_policy`, `mergeable_genes`,
`excluded_failed_genes`, and `infrastructure_errors`. Failed genes are never read during merging,
even if an older artifact remains in their directory; aggregate downsampling state is rebuilt only
from the mergeable subset.

## Step 0 - Validate inputs (optional but recommended)

Before running the pipeline, sanity-check your files parse cleanly:

```bash
trackcluster validate-bed --input reads.bed
trackcluster validate-bed --input ref.bed
```

If parsing fails, fix the BED12/bigGenePred formatting first (see
[`FORMATS.md`](FORMATS.md)). Use `--lenient` only for the explicitly documented
legacy-repair subset; `trackcluster validate-bed --help` and the
[CLI reference](CLI.md#trackcluster-validate-bed) describe the repair and
report contracts.

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
- `tracktest/sample_gene_paths.tsv`
- `tracktest/sample_novel.bed`
- `tracktest/sample_rejected_reads.tsv`
- `tracktest/<gene-path-key>/<gene-path-key>_nano.bed`
- `tracktest/<gene-path-key>/<gene-path-key>_gff.bed`

`sample_gene.txt` contains biological gene IDs. Use
`sample_gene_paths.tsv` to map them to the encoded keys used in these paths.
Preparation stages bucketed and pooled work before publishing, replaces each
artifact atomically, and publishes `sample_gene.txt` last as the generation
commit marker. If a late filesystem error leaves that file empty, rerun
`preparedir`; batch and count-only paths reject the incomplete generation.

### Step 2 - Cluster per gene in parallel (`clusterj_batch`)

#### What this step does

`clusterj_batch` runs junction-mode clustering for an authoritative gene selection. When preparation is inline, that selection is the just-prepared generation; otherwise `--gene-list` is required and contains one biological gene ID per line. Repeatable `--downsample-gene` values are biological IDs too; neither option accepts encoded directory keys. Arbitrary directory discovery is intentionally disabled so stale gene folders from an earlier preparation are never selected silently.
If you plan to run `trackcluster count --output-root` afterwards, write the batch outputs back into the prepared root so each `<gene-path-key>/` folder still contains `<gene-path-key>_nano.bed`, the key-named per-gene isoform BED, and `<gene-path-key>_read_to_isoform.tsv`. A separate output root is fine for manual concatenation, but it is not a complete per-gene count root unless you also copy the prepared gene inputs.

Internally, 5' truncation collapsing uses a junction-suffix index to avoid quadratic scans on large loci.
Low-support splice-junction sites are corrected before merging with `--junction-correction-offset` and `--junction-correction-min-support`.
Supported alternative SL 5' clusters can be protected from merging with the separate `--sl-*` terminal controls exposed by `trackcluster clusterj`, `trackcluster flow`, and `clusterj_batch`; pass a non-negative `--sw-score` such as `11` when BED score should be used as valid SL/SW 5' evidence.
Supported same-junction 3' terminal clusters can be protected with `--same-junction-3prime-offset`, `--3prime-cluster-offset`, and `--3prime-min-support`.
Use `--platform-preset rna002` for RNA002 reads (junction offset `15`; SL 5' offsets `20/25/20`; 3' cluster offset `15`) and `--platform-preset rna004` for RNA004 reads (conservative defaults: junction offset `10`; SL 5' offsets `15/25/15`; 3' cluster offset `10`). Explicit junction correction, SL, and 3' options override the preset.
SL evidence is optional; reads without SL information still use junction correction and normal 5' truncation collapsing, but they are not protected by the SL-cluster rules. For no-SL human data, keep the junction-mode default `--sw-score -1`.
Same-junction 3' retention is independent of SL evidence and is strand-aware. On the minus strand the biological 3' end lies on the lower-coordinate side, but an early stop truncates that side and therefore has a higher `tx_start` than the full-length isoform. Counting defaults to `--assignment-mode unique`, so each read contributes to only the closest compatible catalog isoform; retained 3' early-stop reads therefore count to the closer terminal isoform instead of a longer same-junction reference. Pass `--assignment-mode fractional` for compatibility with split multi-mapped counts from the mapping file.

For each gene, it writes:

- `*_simple_coveragej.bed` (isoforms)
- `*_unused.bed` (rare/filtered reads)
- `*_read_to_isoform.tsv` (mapping)
- `rejected_reads.tsv` (malformed/empty-ID input read tracks excluded for this gene)

It also writes run summaries in the output root (e.g. `clusterj_batch_summary.txt`).

#### Command (gene-batched)

```bash
clusterj_batch \
  --input-root tracktest \
  --gene-list tracktest/sample_gene.txt \
  --output-root tracktest \
  --threads 8 \
  --force
```

### Step 3 - Combine per-gene outputs into a single isoform file

`clusterj_batch` produces one isoform BED per gene. Many downstream tools expect a single isoform file, so you usually concatenate them:

```bash
find tracktest -mindepth 2 -maxdepth 2 -name '*_simple_coveragej.bed' -print0 \
  | sort -z \
  | xargs -0 cat > sample_isoform.bed
```

Optional: combine unused reads too:

```bash
find tracktest -mindepth 2 -maxdepth 2 -name '*_unused.bed' -print0 \
  | sort -z \
  | xargs -0 cat > sample_unused.bed
```

The recommended count path reads the completed per-gene output folders directly, so manual
concatenation is not needed for counting:

```bash
trackcluster count \
  --reference ref.bed \
  --output-root tracktest \
  --prefix sample
```

## Step 4 - Count isoforms (`count`)

`trackcluster count --output-root` reuses the same per-gene count boundary as `flow --count-only`.
In unique assignment mode, each `<gene-path-key>/` folder is counted using its
own `<gene-path-key>_nano.bed`, key-named per-gene isoform BED, and
`<gene-path-key>_read_to_isoform.tsv`; retained-intron searches do not cross
between gene folders.

```bash
trackcluster count \
  --reference ref.bed \
  --output-root tracktest \
  --prefix sample
```

The legacy standalone BED mode remains available with `--isoform` and `--read-to-isoform`, but its
unique assignment scope is the supplied merged input rather than individual gene folders.

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
- `sample.isoform_usage.group.tsv` (when at least one sample has a non-empty `group`)
- `sample.unique_assignment.provenance.tsv` (unique assignment mode only)

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

See [`../scripts/run_full_flow_rust.sh`](../scripts/run_full_flow_rust.sh) for
an example end-to-end driver.
