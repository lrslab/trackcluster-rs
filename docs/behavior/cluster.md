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

## Rust status
Rust currently provides an overlap-mode `cluster` command implemented in `src/cluster/cluster_overlap.rs`.
- It performs locus splitting first (span-based), then applies the same two-pass distance/merge idea using native overlap calculations (no external tools).
- Parameters are currently fixed to Python defaults:
  - `cutoff1=0.05`, `cutoff2=0.01`, `intronweight=0.5`, `scorecutoff=11`

This is intended as a starting point; full parity should be validated/adjusted using goldens.

