# TrackCluster (Rust rewrite)

Pure-Rust rewrite of the TrackCluster long-read isoform clustering/counting pipeline.

Goals:
- No runtime dependency on `bedtools` (native sort/intersect/cluster primitives).
- CLI parity with the legacy Python `trackrun.py` surface (in-progress).

## Toolchain
Source checkouts pin Rust `1.90.0` via `rust-toolchain.toml` to avoid a known
`EXDEV` artifact-write failure seen with newer toolchains in this environment.

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
- `validate-bed`: strict BED12/bigGenePred input validation, with explicit lenient repair reports
- `bam2bigg`: convert genome-aligned BAM records to TrackCluster bigGenePred-compatible BED12+8
- `gff2bigg`: convert GFF3 or GTF exon annotations to a TrackCluster reference catalog
- `export`: write transcript catalogs as GTF, GFF3, or a SQANTI3 input-audit table

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
ARCHIVE="trackcluster-${TAG}-x86_64-unknown-linux-musl"
curl -fLO "https://github.com/${REPO}/releases/download/${TAG}/${ARCHIVE}.tar.gz"
curl -fLO "https://github.com/${REPO}/releases/download/${TAG}/SHA256SUMS"
grep -F " ${ARCHIVE}.tar.gz" SHA256SUMS > "${ARCHIVE}.sha256"
test -s "${ARCHIVE}.sha256"
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c "${ARCHIVE}.sha256"
else
  shasum -a 256 -c "${ARCHIVE}.sha256"
fi
# Supply-chain verification when GitHub CLI is installed:
if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  gh attestation verify "${ARCHIVE}.tar.gz" --repo "${REPO}"
fi
tar xzf "${ARCHIVE}.tar.gz"
# Current archives may be flat; newer self-contained archives use one top-level directory.
if [ -d "${ARCHIVE}" ]; then cd "${ARCHIVE}"; fi
./trackcluster --help
# Make the unpacked binaries available to the quickstart commands below.
export PATH="$PWD:$PATH"
```

Available targets: Linux x86_64 (musl static), Linux ARM64 (glibc 2.31+), macOS Apple Silicon.

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

## Quickstart (tiny examples)
```bash
# One-line flow: prepare per-gene inputs, run per-gene clustering, merge outputs, count, and desc
trackcluster flow -s examples/reads.bed -r examples/ref.bed -o out --prefix sample
# Tip: disable the default per-gene downsampling cap with `--max-reads-per-gene 0` (uses more memory).
# Independent per-gene downsampling is rejected when one molecule belongs to multiple genes;
# disable the cap or exclude every affected gene from downsampling in that case.
# Malformed or empty-ID read tracks are skipped individually by default and recorded in
# `<prefix>_rejected_reads.tsv` / `<gene-path-key>/rejected_reads.tsv`.
# Add `--invalid-read-policy fail` to restore strict read-track parsing.
# Gene-local failures are logged and excluded while verified genes continue through merge/count/desc.
# Add `--strict-gene-errors` to stop before downstream outputs when any gene fails.

# If per-gene clustering already finished, rerun only merge/count/desc outputs
trackcluster flow --count-only -r examples/ref.bed -o out --prefix sample

# Count from an existing output folder; unique assignment stays inside each gene folder
trackcluster count -r examples/ref.bed -o out --prefix sample

# Validate a BED12/bigGenePred file
trackcluster validate-bed -i examples/minimal.bed

# Convert your own genome-aligned BAM to TrackCluster BED12+8. The default MAPQ cutoff is 30.
# (A BAM is not bundled with the tiny text examples.)
trackcluster bam2bigg --bamfile alignments.bam --out reads.bed

# Convert the packaged GFF3 model to a deterministic reference BED12+8 catalog.
trackcluster gff2bigg --gff examples/annotation.gff3 --out reference.bed

# Junction-mode clustering (writes isoform.bed + mapping + unused)
trackcluster clusterj -s examples/reads.bed -r examples/ref.bed -o isoform.bed
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
# Supported same-junction 3' terminal clusters are retained as isoforms. On the
# minus strand the 3' end is tx_start; an early stop has a higher tx_start than
# the corresponding full-length isoform.

# Overlap-mode clustering (legacy-style two-round exon/intron overlap mode)
trackcluster cluster -s examples/reads.bed -r examples/ref.bed -o isoform.bed

# Full flow in overlap mode
trackcluster flow --cluster-mode cluster -s examples/reads.bed -r examples/ref.bed -o out --prefix sample
# Flow keeps its shared no-SL default (`--sw-score -1`) in either clustering mode.
# Pass `--sw-score 11` to opt into legacy score-based protection in overlap mode.

# Legacy low-level count from a standalone isoform BED. Default unique mode also
# writes isoform_count.provenance.tsv; fractional mode does not.
trackcluster count -s examples/reads.bed -r examples/ref.bed -i isoform.bed --read-to-isoform isoform.read_to_isoform.tsv --out isoform_count.csv

# Describe/classify isoforms vs reference (writes <prefix>_*.txt)
trackcluster desc --isoform isoform.bed --reference examples/ref.bed -o desc_out
```

## Multi-sample pooled usage
Use a manifest TSV to pool reads for clustering once, then quantify per-sample isoform usage.

Example manifest (`samples.tsv`):
```tsv
sample	group	reads
S1	control	/path/S1.reads.bed
S2	treated	/path/S2.reads.bed
```

Run full pooled flow:
```bash
trackcluster flow --manifest examples/samples.tsv -r examples/ref.bed -o out --prefix pooled
```

Add `--emit-pooled-reads` if you also want `<prefix>_pooled_reads.bed` written.

If clustering already completed and you only need to regenerate merged count/description outputs, use `--count-only`. Include `--manifest` when you want the multi-sample usage tables regenerated too:
```bash
trackcluster flow --count-only --manifest examples/samples.tsv -r examples/ref.bed -o out --prefix pooled
```

Or run per-sample quantification from an existing pooled isoform BED:
```bash
trackcluster count-multi --manifest examples/samples.tsv -r examples/ref.bed -i out/pooled_isoform.bed -o out/pooled
```

Tip: with default `--name2-mode coverage` (or `none`), use `--read-to-isoform out/pooled_read_to_isoform.tsv` (or keep the TSV next to the isoform BED for auto-discovery).

For overlap-mode pooled clustering, add `--cluster-mode cluster` to the `flow` command above.

`count-multi` writes:
- `out/pooled.isoform_count.csv`
- `out/pooled.isoform_usage.long.tsv`
- `out/pooled.isoform_counts.matrix.tsv`
- `out/pooled.isoform_usage.group.tsv` (when at least one sample has a non-empty `group`)
- `out/pooled.unique_assignment.provenance.tsv` (default unique mode)

In unique assignment mode, `flow` also writes `<prefix>_read_to_isoform.unique.tsv`, the exact read-to-isoform mapping used for final counts, plus `<prefix>_unique_assignment.provenance.tsv` with the effective `--unique-assignment-junction-offset` and one-to-one/no-collapse matching policy. The raw merged `<prefix>_read_to_isoform.tsv` remains the unselected mapping from per-gene clustering.

The aggregate `out/pooled.isoform_count.csv` is derived from the per-sample matrix: each isoform count is the sum of that isoform's sample columns. In `flow --manifest`, the main `<prefix>_isoform_count.csv` is synchronized from the same aggregate count, so total and per-sample counts use the same assignment result.

New catalogs use deterministic `tc_novel_v1:` structural IDs for novel
isoforms and a percent-encoded `tc_name2_v1:` payload in `--name2-mode full`.
Count CSVs have columns `gene,isoform_id,count` and use standard CSV escaping.
Repeated read labels are treated as one abundance molecule; conflicting
structures for one label are rejected in unique-assignment mode. See
[`docs/FORMATS.md`](docs/FORMATS.md) for the identity and migration contract.

Within the 0.2.0 format contract, rejected-read reporting does not otherwise
change BED, isoform, count, or description/classification schemas and rules.
Skipped reads do not contribute biological evidence, so result contents can
change when an input contains rejected tracks. I/O, reference, configuration,
and algorithm errors are not downgraded by `--invalid-read-policy skip`.

## Docs
- [Changelog](CHANGELOG.md)
- [Pipeline tutorial](docs/PIPELINE.md)
- [CLI reference](docs/CLI.md)
- [File formats](docs/FORMATS.md)
- [Interchange formats](docs/INTERCHANGE.md)
- [Clustering behavior](docs/behavior/cluster.md)
- [Description/classification behavior](docs/behavior/desc.md)

## Development (source checkout only)

These commands require the repository's source and test fixtures; they are not
included as runnable inputs in pre-built binary archives.

```bash
cargo test --all --all-features
```

Junction-cluster and count golden fixtures:
```bash
# Regenerate the clusterj and count goldens from the current Rust implementation
bash tests/generate_goldens.sh
```

## License
Licensed under either of:
- MIT license (`LICENSE-MIT`)
- Apache License, Version 2.0 (`LICENSE-APACHE`)

at your option.
