# Genomics interchange formats

TrackCluster-rs keeps one validated BED12-compatible `Transcript` model inside
the clustering pipeline. Format adapters sit at the boundary so clustering and
counting do not gain format-specific coordinate semantics.

## BAM alignment import

`trackcluster bam2bigg --bamfile alignments.bam --out reads.bed` converts a
genome-aligned BAM directly to TrackCluster bigGenePred-compatible BED12+8
text. `--out` defaults to `bigg.bed`; `--score`/`--min-mapq` defaults to `30`.
Unmapped records are skipped, and secondary (`0x100`) and supplementary
(`0x800`) records are excluded unless `--include-secondary` or
`--include-supplementary` is passed. Records without MAPQ have MAPQ zero.

BAM's one-based alignment start is converted to BED's zero-based coordinate.
Only CIGAR `N` splits exon blocks. Other reference-consuming operations,
including deletions, remain within an exon; insertions and clipping do not
consume reference coordinates. A block containing only deletions is rejected,
as is a span beyond the reference length declared in the BAM header. Flag
`0x10` determines the strand. The emitted BED score is MAPQ, and
forward/reverse records receive item RGB values `250,128,114`/`64,224,208`.

`--group/-g` supplies extra field index `6`; without it, the BAM filename stem
is used. The remaining TrackCluster metadata identifies the row as
`nanopore_read`, leaves gene and `name2` unassigned, and records no CDS or exon
frame. One output row is retained per accepted BAM alignment instance, even
when query names repeat. The completion summary reports each filtering class.
Records are converted as a stream, and a failed conversion leaves any previous
destination untouched.

## GFF3/GTF annotation import

`trackcluster gff2bigg --gff annotation.gff3 --out reference.bed` builds a
reference transcript catalog from exon features. `--out` defaults to
`bigg.bed`. `--input-format auto|gff3|gtf` defaults to `auto`, which detects
GFF3 `key=value` or GTF `key "value"` syntax per row. Blank/comment lines are
ignored and the feature section ends at `##FASTA`.

GFF3 exon rows are grouped by `Parent`; each parent must resolve to a declared
non-gene feature with a unique `ID`. GTF exon-only and full models are grouped
by `transcript_id`. Multiple comma-separated GFF3 parents generate an exon in
each model, with splitting performed before percent decoding so an encoded
comma remains part of an ID. GFF3 gene labels are selected with `--key/-k`
(default `ID`) without changing canonical `ID`/`Parent` graph relationships.
Transcript and exon `gene_id` hints are also accepted; multiple genes are
joined deterministically with `||`, and absent gene annotation becomes `none`.

Annotation coordinates are one-based closed and become zero-based half-open
BED exon intervals. Duplicate exon intervals are collapsed and exons are
sorted, so input order does not affect output. A transcript's outer exon bounds
define its BED span. Cross-contig exons, conflicting known strands, overlapping
blocks, duplicate graph IDs, unresolved or gene-typed exon parents, unsafe
BED fields, invalid reference transcript IDs, and malformed rows fail the
conversion. Declared transcript contig, strand, and containment span are also
validated. A failed conversion leaves any previous destination untouched.

The current adapter is deliberately exon-structure oriented: CDS, UTR, phase,
annotation scores, and declared transcript spans are not transferred. Output
uses score `100`, `itemRgb=0`, no CDS, all exon frames `-1`, `name2=none`, and
`type=isoform_anno`, then sorts models deterministically before writing.

## Transcript export

`trackcluster export --input catalog.bed` accepts any combination of:

- `--gtf catalog.gtf`: GTF 2.2 transcript/exon features with `gene_id` and
  `transcript_id` attributes;
- `--gff3 catalog.gff3`: GFF3 `mRNA`/`exon` features with percent-encoded IDs;
- `--sqanti-input catalog.sqanti.tsv`: an auditable ID/geometry table to retain
  next to the GTF supplied to SQANTI3.

All exports convert BED half-open coordinates to one-based closed coordinates.
TrackCluster does not assign SQANTI structural categories. Use the exported GTF
as input to SQANTI3 when that external classification and QC report is needed;
`--sqanti-input` only writes a compact audit table for the same transcript
catalog.
