# Changelog

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
