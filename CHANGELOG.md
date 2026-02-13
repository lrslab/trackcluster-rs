# Changelog

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
