# Changelog

## 0.1.15
- Change `flow` unique counting to select reads directly inside each per-gene output folder, using `{gene}_nano.bed`, `{gene}_simple_coverage[j].bed`, and `{gene}_read_to_isoform.tsv` instead of falling back to manifest or original read inputs.
- Add `flow --count-only` to reuse completed per-gene cluster outputs and rerun only merge, count, multi-sample count, and desc outputs.
- Emit `<prefix>_read_to_isoform.unique.tsv` in unique assignment mode so the exact mapping used for final counts is available for audit/reproduction.
- Skip genes with incomplete per-gene count inputs during unique assignment, rather than failing the whole final count step.
- Fix unique assignment so catalog expansion does not resurrect unused reads or reads excluded by per-gene downsampling.
- Add regression coverage for count-only reruns, unique audit mapping output, missing per-gene inputs, and downsample-scaled final counts.

## 0.1.14
- Derive `count-multi` aggregate counts from the emitted per-sample matrix and write them to `<prefix>.isoform_count.csv`.
- Synchronize `flow --manifest` `<prefix>_isoform_count.csv` from the same multi-sample aggregate count, so the total count table, sample matrix, long usage table, and group usage table share one assignment result.
- Keep `unique` assignment as the default for `count-multi` and `flow`, including catalog-expanded candidate selection for reads that should assign to a closer isoform than the original mapping candidate.
- Document the new aggregate count output and its invariant: each `*.isoform_count.csv` value is exactly the sum of the corresponding row in `*.isoform_counts.matrix.tsv`.
- Add regression coverage that checks `count-multi` aggregate counts and `flow --manifest` main count outputs against the per-sample matrix.

## 0.1.12
- Make `unique` the default assignment mode for `count`, `count-multi`, and `flow`; pass `--assignment-mode fractional` to keep legacy split-count behavior.
- Expand unique counting against the isoform catalog before choosing the closest compatible isoform, so reads are not trapped by an incomplete `read_to_isoform.tsv` candidate set.
- Use structure-aware unique assignment across mapping-backed and embedded-subread counting paths, including multi-sample manifest reads and pooled `flow` outputs.
- Add fuzzy same-junction merging in `clusterj` using the active junction correction offset, with a same-length junction-chain index to keep candidate scans bounded.
- Raise the default weighted junction correction minimum support to `5`; the `rna002` preset now uses junction correction offset `15` with SL 5' offsets `20/25/20`.
- Update CLI and behavior documentation for catalog-aware unique counting, weighted junction support, RNA002 defaults, and the distinction between junction correction and SL/5' terminal offsets.
- Add regression coverage for fuzzy same-junction merging, catalog-expanded unique assignment, default unique counting, and explicit fractional compatibility.

## 0.1.11
- Add configurable junction correction controls for `clusterj`, `flow`, and `clusterj_batch`: `--junction-correction-offset` and `--junction-correction-min-support`.
- Add `--platform-preset generic|rna002|rna004`; `generic` preserves current defaults, `rna002` widens junction correction and SL 5' offsets for RNA002/DEI-style workflows, and `rna004` keeps conservative/default RNA004 cutoffs.
- Let explicit junction correction and SL 5' options override preset values across single-file, flow, and batched clustering.
- Record the active platform preset and junction correction settings in batch summary files for reproducibility.
- Document the distinction between internal junction correction offsets and SL/5' terminal merge/protection offsets.
- Add regression coverage for junction correction offset/min-support behavior, preset expansion, preset override precedence, and summary output.

## 0.1.10
- Speed up SW-aware `clusterj` merging by grouping SL 5' support by junction chain and avoiding redundant merge-target scans.
- Add exact-duplicate representative pruning for non-reference reads with identical corrected transcript structure.
- Add `unc52` regression fixtures for SW-aware junction clustering.

## 0.1.9
- Add SL-aware junction merge controls for `clusterj`, `flow`, and `clusterj_batch`: `--sl-partial-5prime-offset`, `--sl-same-junction-5prime-offset`, `--sl-5prime-cluster-offset`, and `--sl-5prime-min-support`.
- Keep supported alternative SL 5' clusters as candidate isoforms while still merging singleton likely-degradation reads.
- Honor `--sw-score -1` during batched junction merging so truncation collapsing is fully disabled when requested.
- Record SL merge settings in batch summary files and add regression coverage for the new merge boundaries.
- Add an ignored large overlap speed probe for manual 100k+ read performance checks.

## 0.1.8
- Fix overlap `cluster` reference handling so plain BED12 references are protected by input source instead of requiring `ttype=isoform_anno`.
- Report overlap-mode reads with no matching reference chromosome, strand, or locus in `unused.bed` instead of silently dropping them.
- Speed up large overlap loci with sparse exon-overlap candidate generation while preserving exact all-pairs behavior for small or high-cutoff cases.
- Fix `flow --cluster-mode cluster` reruns so per-gene outputs are regenerated unless isoform BED, unused BED, and read-to-isoform TSV are all present.
- Add overlap-mode goldens, unused-read regressions, deterministic threading coverage, and batch-size benchmarks.

## 0.1.7
- Add legacy overlap/intersection clustering to `flow` via `--cluster-mode cluster`.
- Add overlap-mode CLI controls for `cluster`/`flow`: `--batch-size`, `--batch-rounds`, `--sw-score`, cutoff tuning, and `--name2-mode`.
- Keep second-round `SL` reads as their own track when `score == --sw-score`, matching the original TrackCluster boundary behavior.
- Release hygiene: add checked-in license texts, pin CI/release builds to Rust `1.90.0`, and keep release tarballs out of the main branch.

## 0.1.6
- Performance: speed up `clusterj` 5' truncation collapsing on large loci by indexing junction suffixes (avoids quadratic scans).
- Bench: add a large single-locus `clusterj` benchmark to track this workload.

## 0.1.5
- Default `--name2-mode` is now `coverage` (smaller isoform BEDs; rely on mapping TSVs for counting).
- Default `--max-reads-per-gene` is now `50000` for `flow`/`clusterj_batch` (memory-friendly; set `0` to disable). Counts/usage tables are scaled when downsampling occurs.
- Add `--heartbeat-seconds` / `--heartbeat-top` to periodically report progress and in-flight genes during `flow`/`clusterj_batch`.
- `count` and `count-multi` auto-discover `*_read_to_isoform.tsv` next to the isoform BED when present.

## 0.1.4
- Add `--name2-mode` (full/coverage/none) to control isoform `name2` payload size while keeping a read-to-isoform mapping.
- Add `--read-to-isoform` fast path for `count` and `count-multi` to reuse mapping TSVs from `flow`/`clusterj`.
- Add per-gene downsampling for `flow`/`clusterj-batch` (`--max-reads-per-gene`, `--downsample-gene`, `--downsample-seed`) and scale counts/usage tables via `clusterj_batch_downsample.tsv`.
- Change manifest mode: `<prefix>_pooled_reads.bed` is now written only when `flow --emit-pooled-reads` is set.
- Performance: bucketed two-pass `preparedir` for large inputs + faster BED12 reading.
- Default `--sw-score` is now `11` (TrackCluster Python default); use `-1` to disable collapsing.

## 0.1.3
- Add `count-multi` subcommand for per-sample isoform usage from pooled isoforms
- Add `flow --manifest` mode:
  - pool reads from a sample manifest (`<prefix>_pooled_reads.bed`)
  - cluster once into shared isoforms
  - emit multi-sample usage tables (`<prefix>.isoform_usage.long.tsv`, `<prefix>.isoform_counts.matrix.tsv`, optional group table)
- Add manifest TSV parser (`sample`, `reads`, optional `group`) with strict validation and tests
- Add integration fixtures/tests for multi-sample counting and flow manifest mode
- Release: build Linux x86_64 artifact as `x86_64-unknown-linux-musl` to avoid host glibc version mismatch errors

## 0.1.0
- Initial Rust CLI: `validate-bed`, `clusterj`, `cluster`, `count`, `addgene`, `desc`, `preparedir`
- `flow` subcommand: one-command end-to-end pipeline (preparedir + clusterj batch + count + desc)
- Native interval utilities (no runtime shell-out)
- Small fixtures + golden-based integration tests for `clusterj` and `count`
- Pin Rust 1.90.0 via `rust-toolchain.toml`
- CI: lint/test workflow + automated release with pre-built binaries for Linux and macOS
