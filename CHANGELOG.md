# Changelog

## 0.1.0
- Initial Rust CLI: `validate-bed`, `clusterj`, `cluster`, `count`, `addgene`, `desc`, `preparedir`
- `flow` subcommand: one-command end-to-end pipeline (preparedir + clusterj batch + count + desc)
- Native interval utilities (no runtime shell-out)
- Small fixtures + golden-based integration tests for `clusterj` and `count`
- Pin Rust 1.90.0 via `rust-toolchain.toml`
- CI: lint/test workflow + automated release with pre-built binaries for Linux and macOS

