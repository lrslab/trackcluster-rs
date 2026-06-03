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
- `clusterj`: junction-chain clustering (fast mode; SL-aware 5' merge controls; optimized truncation collapsing for large loci)
- `cluster`: overlap-based clustering (slower, more permissive)
- `count`: isoform expression counting
- `count-multi`: per-sample (and optional per-group) isoform usage from pooled isoforms
- `desc`: novel isoform description/classification vs reference
- `addgene`: assign gene names to reads by overlap with reference
- `validate-bed`: basic BED12/bigGenePred input validation

Extra binary:
- `clusterj_batch`: run `clusterj` per gene folder in parallel (manual junction-mode batched runner; overlap-mode batching is exposed through `trackcluster flow --cluster-mode cluster`)

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
# One-line flow: prepare per-gene inputs, run per-gene clustering, merge outputs, count, and desc
trackcluster flow -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o out --prefix sample
# Tip: disable the default per-gene downsampling cap with `--max-reads-per-gene 0` (uses more memory).

# If per-gene clustering already finished, rerun only merge/count/desc outputs
trackcluster flow --count-only -r tests/fixtures/ref.bed -o out --prefix sample

# Count from an existing output folder; unique assignment stays inside each gene folder
trackcluster count -r tests/fixtures/ref.bed -o out --prefix sample

# Validate a BED12/bigGenePred file
trackcluster validate-bed -i tests/fixtures/minimal.bed

# Junction-mode clustering (writes isoform.bed + mapping + unused)
trackcluster clusterj -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o isoform.bed
# Platform presets:
#   --platform-preset rna002  # junction offset 15; SL 5' offsets 20/25/20; 3' cluster offset 15
#   --platform-preset rna004  # conservative defaults: junction offset 10; SL 5' offsets 15/25/15; 3' cluster offset 10
# Junction-mode defaults treat reads as no-SL (`--sw-score -1`). Pass
# `--sw-score 11` only when BED score is valid SL/SW 5' evidence.
# SL 5' merge behavior can be tuned with --sl-partial-5prime-offset,
# --sl-same-junction-5prime-offset, --sl-5prime-cluster-offset, and
# --sl-5prime-min-support.
# Same-junction 3' retention can be tuned with --same-junction-3prime-offset,
# --3prime-cluster-offset, and --3prime-min-support.
# SL evidence is optional. Reads without SL information use the normal junction
# correction and 5' truncation collapse path, but are not SL-protected isoforms.
# Supported same-junction 3' terminal clusters are retained as isoforms, including
# minus-strand early-stop clusters where the 3' end is the lower genomic coordinate.

# Overlap-mode clustering (legacy-style two-round exon/intron overlap mode)
trackcluster cluster -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o isoform.bed

# Full flow in overlap mode
trackcluster flow --cluster-mode cluster -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -o out --prefix sample

# Legacy low-level count from a standalone isoform BED
trackcluster count -s tests/fixtures/reads.bed -r tests/fixtures/ref.bed -i isoform.bed --read-to-isoform isoform.read_to_isoform.tsv --out isoform_count.csv

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

If clustering already completed and you only need to regenerate merged count/description outputs, use `--count-only`. Include `--manifest` when you want the multi-sample usage tables regenerated too:
```bash
trackcluster flow --count-only --manifest samples.tsv -r tests/fixtures/ref.bed -o out --prefix pooled
```

Or run per-sample quantification from an existing pooled isoform BED:
```bash
trackcluster count-multi --manifest samples.tsv -r tests/fixtures/ref.bed -i out/pooled_isoform.bed -o out/pooled
```

Tip: with default `--name2-mode coverage` (or `none`), use `--read-to-isoform out/pooled_read_to_isoform.tsv` (or keep the TSV next to the isoform BED for auto-discovery).

For overlap-mode pooled clustering, add `--cluster-mode cluster` to the `flow` command above.

`count-multi` writes:
- `out/pooled.isoform_count.csv`
- `out/pooled.isoform_usage.long.tsv`
- `out/pooled.isoform_counts.matrix.tsv`
- `out/pooled.isoform_usage.group.tsv` (only when manifest has `group`)

In unique assignment mode, `flow` also writes `<prefix>_read_to_isoform.unique.tsv`, the exact read-to-isoform mapping used for final counts. The raw merged `<prefix>_read_to_isoform.tsv` remains the unselected mapping from per-gene clustering.

The aggregate `out/pooled.isoform_count.csv` is derived from the per-sample matrix: each isoform count is the sum of that isoform's sample columns. In `flow --manifest`, the main `<prefix>_isoform_count.csv` is synchronized from the same aggregate count, so total and per-sample counts use the same assignment result.

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
