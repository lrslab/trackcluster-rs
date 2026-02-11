# Bedtools Audit (Python TrackCluster)

This document lists every `bedtools` operation used by the legacy Python implementation in this repo and maps it to (current or planned) Rust equivalents.

## Operations Used

### `bedtools sort`
- Used to create sorted inputs for downstream intersect/merge.
- Python wrappers:
  - `trackcluster/tracklist.py` (`wrapper_bedtools_intersect2`)
  - `trackcluster/pre.py` (`wrapper_bedtools_merge`, `wrapper_bedtools_subtract`)

**Rust mapping**
- For bigGenePred/BED12 transcripts: `src/interval/sort.rs` (`sort_by_coord`).

### `bedtools intersect`

#### A) Pairwise intersection output (`-wa -wb`)
- Command shape: `bedtools intersect -wa -wb -a A -b B > out`
- Python wrapper:
  - `trackcluster/tracklist.py` (`wrapper_bedtools_intersect2`)
- Used by:
  - `trackcluster/cluster.py` to compute exon/intron overlap bp between all pairs (distance matrix).

**Rust mapping**
- Span-based candidate pairs: `src/interval/intersect_sweep.rs` (`sweep_intersect_pairs`).
- Exact bp overlap for exon/intron intervals: `src/interval/refine.rs` (`exonic_overlap_bp`) plus the same sweep/refinement pattern for introns.
- Current overlap-mode clustering implementation: `src/cluster/cluster_overlap.rs` (does not require writing intermediate BED8 files).

#### B) Reciprocal/thresholded intersection (`-s -f -F -nonamecheck`)
- Command shape: `bedtools intersect -nonamecheck -wa -wb -s -f f1 -F f2 -a reads -b refs > out`
- Python wrapper:
  - `trackcluster/pre.py` (`wrapper_bedtools_intersect2_select`)
- Used by:
  - `trackcluster/flow.py` (`flow_add_gene`) to assign gene names to read tracks.

**Rust mapping**
- Implemented in `src/annotate/addgene.rs`:
  - Uses `sweep_intersect_pairs` (span overlap) and applies the same `-f`/`-F` style fraction thresholds.

### `bedtools merge`
- Command shape (modern): `bedtools merge -nonamecheck -s -c 6 -o distinct,count -i sorted.bed`
- Python wrapper:
  - `trackcluster/pre.py` (`wrapper_bedtools_merge`)
- Used by:
  - Novel-gene region-mark generation in `trackcluster/flow.py`.

**Rust mapping (planned)**
- Sort records (`sort_by_coord`) then merge spans per `(chrom,strand)` while preserving an aggregate count.
- A starting point for the merging primitive exists as locus clustering: `src/interval/cluster_span.rs` (`cluster_by_span`).

### `bedtools subtract`
- Command shape: `bedtools subtract -nonamecheck -s -A -f f1 -F f2 -a A -b B > out`
- Python wrapper:
  - `trackcluster/pre.py` (`wrapper_bedtools_subtract`)
- Used by:
  - Novel-gene region filtering in `trackcluster/flow.py`.

**Rust mapping (planned)**
- Strand-partition both inputs, sweep to find overlaps, and emit residual intervals for `-A` semantics.

## Notes on Semantics
- Coordinate system is BED: 0-based, half-open `[start, end)`.
- `-s` enforces strand matching.
- `-f`/`-F` thresholds are fractions of the feature length overlapped (A and B respectively).
- `-A` (subtract) drops the entire A feature if it overlaps B above thresholds.

## Status Summary
- Implemented in Rust: sort, span intersect pairs, addgene overlap assignment, overlap-mode clustering (no external tools).
- Not yet ported: novel-gene workflow pieces that rely on merge/subtract.

