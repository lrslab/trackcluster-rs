# `desc` behavior (Python parity notes)

This summarizes the legacy Python behavior implemented in `trackcluster/post.py` and `trackcluster/flow.py` (`flow_desc_annotation`).

## Inputs
- Isoform-like BED (often produced by `clusterj`, or `addgene` output run on reads).
- Reference BED containing `isoform_anno` records with gene names.

For meaningful output, each input transcript should have a gene name field populated (Python uses `geneName`; Rust currently reads/writes bigGenePred extra fields).

## Outputs (Rust CLI)
The Rust `desc` command writes four files (prefix controlled by `-o/--out`):
- `<prefix>_desc.txt`: `isoform_id`, `ref_id`, `gene`, `miss_desc`, `extra_desc`
- `<prefix>_class4.txt`: `isoform_id`, `class4`
- `<prefix>_fusion.txt`: `isoform_id`, `gene1;gene2;...` (only for multi-gene entries)
- `<prefix>_class12.txt`: `isoform_id`, `class12`

## Classification Logic

### `class4`
One of:
- `new_junction`: transcript has at least one junction not matching any reference junction (within an offset)
- `all_matched>=_<ref>`: junctions match a reference; transcript span length is >= reference
- `all_matched_<_<ref>`: junctions match a reference; transcript span length is < reference
- `new_combination`: no new junctions, but junction chain does not exactly match any single reference

### `desc` row (miss/extra)
The nearest reference is chosen by minimizing, in order:
1) number of extra junction groups
2) number of extra junctions
3) number of missed junction groups
4) number of missed junctions

Then the miss/extra descriptions are derived by grouping consecutive missed/extra boundary indices into:
- 5′/3′ primer missing/extra
- intron retention
- exon miss/extra

### `class12`
Derived as a bucketization over `class4`, `desc` strings, and fusion gene names.

## Important Parameters
- Junction fuzz/offset is `10bp` in Python; Rust currently uses `10bp`.
- Fusion detection uses overlap-based gene assignment (similar to Python `flow_fusion`):
  - `--fusion-fraction-read` (default `0.1`): minimum fraction of isoform span overlapped by a reference span
  - `--fusion-fraction-ref` (default `0.1`): minimum fraction of reference span overlapped by an isoform span
  - An isoform is labeled as fusion if it overlaps **2+ distinct genes** after applying these thresholds (input gene field is ignored for fusion detection).

## Optional: end-shift (UTR-like) tagging (Rust-only)

Python’s `desc` compares internal exon/intron boundaries (junction chain) and does not explicitly classify isoforms that have the same splice chain but different transcript ends (e.g. longer/shorter first/last exon, “UTR-like” changes).

Rust can optionally add strand-aware end-shift tags when the junction comparison is splice-equal (no missed/extra boundary indices):
- Enable with `trackcluster desc ... --end-shift-bp <N>` (`0` disables; default is `0` for Python parity).
- Adds to the `miss_desc` / `extra_desc` fields:
  - `5 end extension: <bp>` / `5 end truncation: <bp>`
  - `3 end extension: <bp>` / `3 end truncation: <bp>`
- Updates `class12` to bucket these splice-equal end-shift isoforms separately (`5'end_extension`, `3'end_extension`, `5'end_truncation`, `3'end_truncation`) instead of `full_matched>=/<`.
