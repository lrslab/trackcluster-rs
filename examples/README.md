# Packaged examples

These tiny inputs are shipped with the pre-built release archive so the
commands in the top-level README can be run immediately after extraction.

- `minimal.bed`: two plain BED12 transcripts for `validate-bed`.
- `annotation.gff3`: one two-exon GFF3 transcript for `gff2bigg`.
- `reads.bed` and `ref.bed`: one small read/reference pair for clustering,
  flow, counting, and description examples.
- `samples.tsv`, `S1.reads.bed`, and `S2.reads.bed`: a two-sample manifest for
  pooled discovery and `count-multi` examples. Manifest paths are relative to
  the directory containing `samples.tsv`.

All files are synthetic and intentionally small; they demonstrate file layout
and command wiring, not biological performance.
