# TrackCluster (Rust rewrite)

Pure-Rust rewrite of the TrackCluster long-read isoform clustering/counting pipeline.

Goals:
- No runtime dependency on `bedtools` (native sort/intersect/cluster primitives).
- CLI parity with the legacy Python `trackrun.py` surface (in-progress).

## Toolchain
This repo pins Rust `1.90.0` via `rust-toolchain.toml` to avoid a known `EXDEV`
artifact-write failure seen with newer toolchains in this environment.

## Status
Implemented subcommands:
- `flow`: one-command end-to-end pipeline (recommended)
- `preparedir`: split reads into per-gene folders (and write `<prefix>_gene.txt`, `<prefix>_dedup.bed`, `<prefix>_novel.bed`)
- `clusterj`: junction-chain clustering (fast mode)
- `cluster`: overlap-based clustering (slower, more permissive)
- `count`: isoform expression counting
- `desc`: novel isoform description/classification vs reference
- `addgene`: assign gene names to reads by overlap with reference
- `validate-bed`: basic BED12/bigGenePred input validation

Extra binary:
- `clusterj_batch`: run `clusterj` per gene folder in parallel (can also run `preparedir` first)

## Install

### From source (recommended)
```bash
cargo build --release
./target/release/trackcluster --help
```

### Install binaries into `~/.cargo/bin`
```bash
cargo install --path . --locked --bins
trackcluster --help
clusterj_batch --help
```

## Quickstart (tiny fixtures)
```bash
# One-line flow: prepare per-gene inputs, run clusterj batch, merge outputs, count, and desc
trackcluster flow -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o out --prefix sample

# Validate a BED12/bigGenePred file
trackcluster validate-bed -i tests/fixtures/minimal.bed

# Junction-mode clustering (writes isoform.bed + mapping + unused)
trackcluster clusterj -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o isoform.bed

# Count isoforms
trackcluster count -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -i isoform.bed -o isoform_count.csv

# Describe/classify isoforms vs reference (writes <prefix>_*.txt)
trackcluster desc --isoform isoform.bed --reference tests/fixtures/ref.bed -o desc_out
```

## Docs
- Pipeline tutorial: `docs/PIPELINE.md`
- CLI reference: `docs/CLI.md`
- Formats: `docs/FORMATS.md`
- Design notes: `docs/design/bedtools_audit.md`
- Behavior notes: `docs/behavior/`

## Testing
```bash
cargo test --all --all-features
```

Golden fixtures:
```bash
# Regenerate goldens from the current Rust implementation
bash tests/generate_goldens.sh
```

## License
Licensed under either of:
- MIT license (`LICENSE-MIT`)
- Apache License, Version 2.0 (`LICENSE-APACHE`)

at your option.
