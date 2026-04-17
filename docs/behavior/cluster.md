# `cluster` behavior (overlap/intersection mode)

Legacy Python overlap-mode clustering is implemented in `trackcluster/cluster.py` and relies on exon/intron intersections (via `bedtools intersect`) to build distance matrices.

## High-level flow (Python)
1) Build per-read exon/intron BED intervals.
2) Compute pairwise overlap bp for exon and intron intervals.
3) Convert overlaps into distances:
   - `ratio`: `1 - overlap / union_len`
   - `ratio_short`: `1 - overlap / min_len`
4) Combine exon and intron distances with `intronweight`:
   - `(D_exon + intronweight * D_intron) / (1 + intronweight)`
5) Run two passes of filtering/merging:
   - pass 1: `ratio` cutoff (`cutoff1`, default `0.05`)
   - pass 2: `ratio_short` cutoff (`cutoff2`, default `0.01`)
6) Merge reads into the retained representative by exon length, with the special-case SL score cutoff (`scorecutoff`, default `11`) controlling some short-read merges.
   - In the Python behavior, second-pass short-read collapse is gated by `score < scorecutoff`, not `<=`.

## Rust status
Rust provides overlap-mode clustering in `src/cluster/cluster_overlap.rs`, exposed through:
- `trackcluster cluster`
- `trackcluster flow --cluster-mode cluster`

Behavior:
- It performs locus splitting first (span-based), then applies the same two-pass distance/merge idea using native overlap calculations (no external tools).
- For large loci, flow/CLI can optionally pre-merge reads in batches before the final full two-pass overlap clustering (`--batch-size`, `--batch-rounds`).
- Parameters default to the Python flow defaults and are configurable from the CLI/flow surface:
  - `cutoff1=0.05`, `cutoff2=0.01`, `intronweight=0.5`, `scorecutoff=11`
- The second pass matches the Python SL boundary behavior: a short read is collapsed only when `score < scorecutoff`; reads with `score == scorecutoff` are retained as their own track.
- In batched `flow --cluster-mode cluster` runs, per-gene overlap outputs use the `*_simple_coverage.bed` suffix and batch summary files use the `cluster_batch_*` prefix.

Current CLI exposure:
- `trackcluster cluster` runs overlap mode directly on one reads/reference pair.
- `trackcluster flow --cluster-mode cluster` runs the same overlap mode per gene after `preparedir`.
- There is no separate `cluster_batch` binary; the dedicated batch binary remains `clusterj_batch` for junction mode.

This is intended as a starting point; full parity should be validated/adjusted using goldens.
