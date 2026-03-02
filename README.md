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
- `clusterj`: junction-chain clustering (fast mode; optimized truncation collapsing for large loci)
- `cluster`: overlap-based clustering (slower, more permissive)
- `count`: isoform expression counting
- `count-multi`: per-sample (and optional per-group) isoform usage from pooled isoforms
- `desc`: novel isoform description/classification vs reference
- `addgene`: assign gene names to reads by overlap with reference
- `validate-bed`: basic BED12/bigGenePred input validation

Extra binary:
- `clusterj_batch`: run `clusterj` per gene folder in parallel (can also run `preparedir` first)

## Install

### Pre-built binaries (recommended)
Download a tarball for your platform from the
[latest GitHub release](https://github.com/lrslab/trackcluster-rs/releases/latest):

```bash
# Example: Linux x86_64
REPO=lrslab/trackcluster-rs
TAG="$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" | sed -n 's/.*"tag_name": "\([^"]*\)".*/\1/p' | head -n1)"
curl -fLO "https://github.com/${REPO}/releases/download/${TAG}/trackcluster-${TAG}-x86_64-unknown-linux-musl.tar.gz"
tar xzf "trackcluster-${TAG}-x86_64-unknown-linux-musl.tar.gz"
./trackcluster --help
```

Available targets: Linux x86_64 (musl static), Linux ARM64 (glibc), macOS Apple Silicon.

### From source
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
# Tip: disable the default per-gene downsampling cap with `--max-reads-per-gene 0` (uses more memory).

# Validate a BED12/bigGenePred file
trackcluster validate-bed -i tests/fixtures/minimal.bed

# Junction-mode clustering (writes isoform.bed + mapping + unused)
trackcluster clusterj -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o isoform.bed

# Count isoforms
trackcluster count -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -i isoform.bed --read-to-isoform isoform.read_to_isoform.tsv -o isoform_count.csv

# Describe/classify isoforms vs reference (writes <prefix>_*.txt)
trackcluster desc --isoform isoform.bed --reference tests/fixtures/ref.bed -o desc_out
```

## Multi-sample pooled usage
Use a manifest TSV to pool reads for clustering once, then quantify per-sample isoform usage.
For a complete real-data walkthrough (including full `samples.tsv` details), see `docs/DEMO_488.md`.

Example manifest (`samples.tsv`):
```tsv
sample	group	reads
S1	control	/path/S1.reads.bed
S2	treated	/path/S2.reads.bed
```

Run full pooled flow:
```bash
trackcluster flow --manifest samples.tsv -r tests/fixtures/ref.bed -o out --prefix pooled
```

Add `--emit-pooled-reads` if you also want `<prefix>_pooled_reads.bed` written.

Or run per-sample quantification from an existing pooled isoform BED:
```bash
trackcluster count-multi --manifest samples.tsv -r tests/fixtures/ref.bed -i out/pooled_isoform.bed -o out/pooled
```

Tip: with default `--name2-mode coverage` (or `none`), use `--read-to-isoform out/pooled_read_to_isoform.tsv` (or keep the TSV next to the isoform BED for auto-discovery).

`count-multi` writes:
- `out/pooled.isoform_usage.long.tsv`
- `out/pooled.isoform_counts.matrix.tsv`
- `out/pooled.isoform_usage.group.tsv` (only when manifest has `group`)

## Docs
- 488 real-data demo (full walkthrough + manifest details): `docs/DEMO_488.md`
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
