# Genomics interchange formats

TrackCluster-rs keeps one validated BED12-compatible `Transcript` model inside
the clustering pipeline. Format adapters sit at the boundary so clustering and
counting do not gain format-specific coordinate semantics.

## Align reads to the genome

Alignment runs upstream of TrackCluster-RS. Start with a genome FASTA and
transcript annotation from the same assembly, with matching sequence names
such as `chr1` versus `1`. The discovery and direct-counting workflows compare
genomic exon coordinates, so the reads BAM must use those genome coordinates.
A BAM aligned to transcript sequences needs external coordinate projection
or realignment to the genome before it can be used for these workflows.
`bam2bigg` converts the coordinates already present in the BAM; it does not
perform that projection or verify the genome assembly against the annotation.

### Minimap2 from FASTQ

For Nanopore direct RNA, a starting command is:

```bash
set -euo pipefail
minimap2 -t 8 -ax splice -uf -k14 genome.fa reads.fastq.gz \
  | samtools sort -@ 4 -o reads.aligned.bam -
samtools index reads.aligned.bam
```

Choose the preset and orientation for the library, following the
[minimap2 RNA alignment guide](https://github.com/lh3/minimap2#map-long-mrnacdna-reads):

| Input | Starting minimap2 options |
| --- | --- |
| Nanopore direct RNA | `-ax splice -uf -k14` |
| Nanopore cDNA with unknown transcript orientation | `-ax splice`; review the strand limitation below before import |
| PacBio Iso-Seq/HiFi transcript reads | `-ax splice:hq -uf` |

`-a` emits SAM, which the pipeline above converts to BAM. `-uf` assumes the
read follows transcript orientation; it does not restrict genomic alignments
to the plus strand. Optional `--junc-bed annotation.bed12` supplies annotated
junctions. When reusing a `.mmi` index, build it with the intended preset and
index parameters: mapping-time `-k14` cannot change an existing index's k-mer
size.

The commands use [coordinate sorting](https://www.htslib.org/doc/samtools-sort.html)
and [BAM indexing](https://www.htslib.org/doc/samtools-index.html) for downstream
tools that need random access. `bam2bigg` itself reads the BAM sequentially
and requires neither sorting nor a `.bai`/`.csi` index. It accepts BAM rather
than the SAM or PAF output from an aligner.

### Dorado from existing BAM basecalls

When basecalls already carry Dorado poly(A) or modification tags, keep the
BAM as the alignment input:

```bash
set -euo pipefail
dorado aligner genome.fa reads.dorado.bam --mm2-opts "-x splice -k 14" \
  | samtools sort -@ 4 -o reads.aligned.bam -
samtools index reads.aligned.bam
```

The [Dorado alignment documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/)
currently describes `lr:hq` as the default preset. Specify `-x splice` for
spliced RNA alignment. Dorado exposes a subset of minimap2 options; check
`dorado aligner genome.fa reads.dorado.bam --mm2-opts "--help"` for the
installed version instead of copying every standalone minimap2 flag.

A plain FASTQ export does not retain BAM tags such as `pt:i` or `MM`/`ML`/`MN`.
Keep the original BAM and its read names and caller `@PG` provenance.
[Dorado's SAM specification](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/)
defines `pt:i`, while its
[modification documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/)
describes modification tags. `bam2bigg` converts exon geometry and does not
copy those tags into BED. Supply the source BAM separately for
[poly(A) aggregation](POLYA.md) or the aligned modBAM for
[`mod-import-dorado`](CLI.md#trackcluster-mod-import-dorado).
Poly(A) aggregation can use the original unaligned Dorado BAM because it
joins by read ID; modification projection requires genome-aligned records.
The modification importer requires `MM`, `ML`, and an integer `MN` matching
the current `SEQ` length. Keep those tags consistent with the sequence through
any trimming or reorientation before import.

For Nanopolish tail estimation, use a sorted/indexed BAM with the matching
reads and genome, plus the raw-signal index required by Nanopolish. Its
[poly(A) tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html)
uses `map-ont` for an unspliced control reference and explicitly advises
splice-aware alignment for native mRNA against a genome. Feed its resulting
TSV into [TrackCluster's Nanopolish poly(A) input](POLYA.md#input-from-nanopolish).

## BAM alignment import

`trackcluster bam2bigg --bamfile alignments.bam --out reads.bed` converts a
genome-aligned BAM directly to TrackCluster bigGenePred-compatible BED12+8
text. `--out` defaults to `bigg.bed`; `--score`/`--min-mapq` defaults to `30`.
Unmapped records are skipped, and secondary (`0x100`) and supplementary
(`0x800`) records are excluded unless `--include-secondary` or
`--include-supplementary` is passed. Missing MAPQ is treated as zero for filtering.

BAM's one-based alignment start is converted to BED's zero-based coordinate.
Only CIGAR `N` splits exon blocks. Other reference-consuming operations,
including deletions, remain within an exon; insertions and clipping do not
consume reference coordinates. A block containing only deletions is rejected,
as is a span beyond the reference length declared in the BAM header. Flag
`0x10` determines the strand. The converter does not interpret transcript-strand
tags such as minimap2's `ts:A` or `XS:A`. It assumes the alignment strand is
the biological transcript strand used for strand-aware gene assignment and
counting. For cDNA containing both transcript orientations, review and normalize
orientation upstream before import; splice-aware alignment alone does not
resolve this limitation. The emitted BED score is `0` (no SL evidence), and
forward/reverse records receive item RGB values `250,128,114`/`64,224,208`.

MAPQ only controls the alignment filter; `--score`/`--min-mapq` does not set
the output BED score. Per-record MAPQ remains in the source BAM.
`bam2bigg` does not compute a Smith-Waterman score or import SLRanger's
`SL_score`. Its `Transcript`/BED output preserves aligned exon boundaries but
does not retain the original CIGAR, softclip lengths, or clipped sequence and
quality evidence. A converted read's score `0` supplies no SL support at the
legacy cutoff of `11`. The new prediction-module
[SL evidence contract](https://github.com/lrslab/trackcluster-rs/blob/main/docs/design/prediction_sl_evidence_contract.md) requires
separate typed MAPQ, SLRanger, and SW values with their evidence provenance.

Older converter output incorrectly used MAPQ as BED score. Regenerate that
output from BAM, or use `--sw-score -1` while reusing those older BED files.

`--group/-g` supplies extra field index `6`; without it, the BAM filename stem
is used. The remaining TrackCluster metadata identifies the row as
`nanopore_read`, leaves gene and `name2` unassigned, and records no CDS or exon
frame. One output row is retained per accepted BAM alignment instance, even
when query names repeat. `--invalid-record-policy skip` is the default: an
independently decoded record with an invalid name, reference/start, CIGAR, span,
or transcript geometry is excluded and counted by a stable reason, while later
records continue. `fail` restores strict record conversion. Header,
BGZF/framing, record-decode, truncation, and write errors remain fatal, as does
an all-invalid candidate set. Records are converted as a stream, and a failed
conversion leaves any previous destination untouched.

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
define its BED span. By default, `--invalid-record-policy skip` quarantines an
entire identifiable transcript for cross-contig exons, conflicting known
strands, overlapping blocks, duplicate graph IDs, unresolved or gene-typed
parents, unsafe BED fields, invalid reference IDs, and attributable malformed
rows. This whole-model boundary prevents one ignored exon from silently
truncating a transcript. Unowned malformed model rows, stream errors, no-exon
inputs, and all-quarantined catalogs remain fatal; `fail` restores strict
all-or-nothing behavior. Rejections are written to `<out>.rejected.tsv`, or to
`--rejected-records PATH`. The report and BED are atomic per file rather than
as one transaction, and the BED is published last as the commit point.

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
