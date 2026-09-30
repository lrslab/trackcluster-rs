# Count against a given isoform catalog

Use `count --assign-against-catalog` for one sample, or
`count-multi --assign-against-catalog` for a sample manifest. These commands
assign reads directly to the supplied isoforms and produce counts without
running preparation, clustering, discovery or classification.

Each distinct input molecule with a candidate contributes one count to its
nearest isoform. The catalog stays fixed, including isoforms with zero counts.
Reads without candidates are reported separately. No discovery mapping,
embedded read memberships, separate reference BED or per-gene output folders
are needed.

## Run the bundled example

Run the commands below from the repository root. For a source checkout, build
and select the current executable first:

```bash
cargo build --release --locked --bin trackcluster
export PATH="$PWD/target/release:$PATH"
trackcluster count --help
```

The help must include `--assign-against-catalog`. If it does not, the executable
predates this feature; rebuild from the updated source and check the executable
selected by `command -v trackcluster`.

### Inputs

The example uses the existing synthetic files in `examples/`:

| Input | Argument | Contents |
| --- | --- | --- |
| [`examples/reads.bed`](../examples/reads.bed) | `--reads` | One aligned read, `read_trunc` |
| [`examples/ref.bed`](../examples/ref.bed) | `--isoform` | Two given isoforms, `ref_a` and `ref_b`, both annotated as `GENEA` |

Both inputs are tab-separated BED12 or bigGenePred-compatible **text**, with
one read alignment or isoform per row. BED column 4 (`name`) supplies the read
ID or isoform ID. Both files must use the same genome assembly, chromosome
names and zero-based, half-open coordinates. Exons are described by
`blockCount`, `blockSizes` and `blockStarts`; each block start is relative to
the row's `chromStart`. See the [BED format](FORMATS.md#required-fields).

In this example, all records are on `chr1`, strand `+`:

| Record | Exon intervals |
| --- | --- |
| `read_trunc` | `[120,130)`, `[140,150)` |
| `ref_a` | `[100,110)`, `[120,130)`, `[140,150)` |
| `ref_b` | `[80,90)`, `[100,110)`, `[120,130)`, `[140,150)` |

### Single-sample command and outputs

```bash
trackcluster count \
  --reads examples/reads.bed \
  --isoform examples/ref.bed \
  --assign-against-catalog \
  --out out/catalog-count/single/count.csv
```

Here `--out` is a CSV filename. Use the full `--out` spelling: for `count`,
`-o` means `--output-root`, which selects the existing-discovery recount mode.
Output parent directories are created automatically.

The command writes these five files under `out/catalog-count/single/`:

| File | Contents |
| --- | --- |
| `count.csv` | Counts for every catalog isoform, including zeros |
| `count.read_to_isoform.tsv` | The exact selected read-to-isoform pairs, without a header |
| `count.unassigned_reads.tsv` | Unassigned read IDs, gene metadata and reasons |
| `count.assignment_stats.tsv` | Input-record and distinct-molecule accounting |
| `count.provenance.tsv` | Effective assignment policy and junction tolerance |

Expected `count.csv`:

```csv
gene,isoform_id,count
GENEA,ref_a,1
GENEA,ref_b,0
```

Expected `count.read_to_isoform.tsv`:

```tsv
read_trunc	ref_a
```

The read's junction matches both candidates. Neither candidate has an extra
intron covered by the read. Their end distances are 20 bp for `ref_a` and
40 bp for `ref_b`, so the read selects `ref_a`. `ref_b` remains in the output
with count zero.

Expected `count.assignment_stats.tsv`:

```tsv
metric	count
input_records	1
total_reads	1
assigned_reads	1
unassigned_reads	0
```

`count.unassigned_reads.tsv` contains only its header in this example:

```tsv
read_id	gene	reason
```

### Multiple samples

[`examples/samples.tsv`](../examples/samples.tsv) points to one reads BED per
sample, using paths relative to the manifest directory:

```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

`sample` and `reads` are required; `group` is optional. Sample names must be
unique and cannot contain `::`. Each sample's reads are genome-aligned BED
records with the same format as the single-sample input.

```bash
trackcluster count-multi \
  --manifest examples/samples.tsv \
  --isoform examples/ref.bed \
  --assign-against-catalog \
  --out out/catalog-count/multi/count
```

For `count-multi`, `--out` (also `-o`) is a **prefix**, not a CSV filename. The
command appends the suffixes below, preserving any dots in the prefix:

| File under `out/catalog-count/multi/` | Contents |
| --- | --- |
| `count.isoform_count.csv` | Aggregate counts, equal to the sum of sample columns |
| `count.isoform_counts.matrix.tsv` | All catalog isoforms, with one column per sample in manifest order |
| `count.isoform_usage.long.tsv` | Nonzero per-sample counts and within-gene usage |
| `count.isoform_usage.group.tsv` | Counts summed within each non-empty group, with within-group gene usage |
| `count.read_to_isoform.tsv` | Selected pairs with `sample::read` IDs |
| `count.unassigned_reads.tsv` | Unassigned `sample::read` IDs and reasons |
| `count.assignment_stats.tsv` | Totals across all samples |
| `count.unique_assignment.provenance.tsv` | Effective assignment policy and tolerance |

The group table is emitted only when at least one sample has a non-empty
`group`. Samples without a group still contribute to sample and aggregate
counts, but not to group summaries.

Expected `count.isoform_counts.matrix.tsv`:

```tsv
gene	isoform_id	S1	S2
GENEA	ref_a	1	1
GENEA	ref_b	0	0
```

Expected `count.isoform_count.csv`:

```csv
gene,isoform_id,count
GENEA,ref_a,2
GENEA,ref_b,0
```

Expected `count.read_to_isoform.tsv`:

```tsv
S1::read_s1	ref_a
S2::read_s2	ref_a
```

The `sample::` prefix is added by the command; supply the original read IDs in
each sample BED. Identical original IDs in different samples remain separate
molecules. Empty read files are allowed and retain zero-valued matrix columns.

For each sample, `proportion = count / gene_total`, where `gene_total` sums
counts for all catalog isoforms with the same gene field. Group proportions
are recomputed from summed counts, not averaged over sample proportions.
Zero-count combinations are omitted from usage tables, but retained in the
count CSV and matrix. These are molecule counts and within-gene fractions;
they are not TPM or library-size-normalized expression values.

## Use your own reads and isoforms

If both inputs are already BED12/bigGenePred text, replace `--reads` and
`--isoform` in the example. The catalog may be a reference annotation, an
external catalog or a filtered discovery catalog. Its IDs and structures are
preserved; the command does not create a new isoform BED.

For genome-aligned BAM and a GTF/GFF3 annotation, convert first:

```bash
trackcluster bam2bigg \
  --bamfile sample.bam --out input/reads.bed
trackcluster gff2bigg \
  --gff annotation.gtf --out input/catalog.bed
trackcluster count \
  --reads input/reads.bed --isoform input/catalog.bed \
  --assign-against-catalog --out out/sample/count.csv
```

Replace `sample.bam` and `annotation.gtf` with your input files. `count` reads
BED text, so FASTQ, BAM and GTF are not direct `--reads`/`--isoform` inputs.
`bam2bigg` uses MAPQ >= 30 by default and excludes secondary and supplementary
alignments; see [conversion rules](INTERCHANGE.md) for other import options.
Fixed-catalog counting considers all reads in the resulting BED without
downsampling.

Gene labels come from the catalog's TrackCluster gene field (column 18,
extra-field index 5). Plain BED12 catalogs are accepted but produce `gene=none`.
For meaningful within-gene usage, use gene-annotated catalogs such as the
output of `gff2bigg`; otherwise all `none` records share one usage denominator.
Groups in the multi-sample outputs come from the manifest's `group` column.

## How the nearest isoform is selected

Candidates must overlap the read on the same chromosome and strand, share at
least one exonic base, and have compatible gene metadata when both sides are
annotated. Unknown strand (`.`) matches only `.`. Multi-gene annotations such
as `GENEA||GENEB` are compatible when their gene sets intersect. If either
side lacks gene metadata, locus, strand and exon overlap determine eligibility.

Among those candidates, minimize these criteria in order:

1. Unmatched read introns.
2. Unmatched catalog introns covered by the read.
3. Absolute start-distance plus absolute end-distance.
4. Exonic symmetric difference within the read's span.
5. Exon-count difference within the read's span.
6. Isoform ID in lexicographic order, for a deterministic tie-break.

Intron matching uses ordered one-to-one matches. The tolerance is
`--unique-assignment-junction-offset 15` by default, applied independently to
each intron boundary. For exact-boundary matching use
`--unique-assignment-junction-offset 0`. This changes how introns are compared;
it does not correct coordinates or impose an endpoint-distance cutoff.

Splice differences are scoring penalties, not rejection criteria. Reads with
extra junctions or single-exon fragments of internal exons still select a
nearest candidate. A selected pair therefore does not imply a perfect
structural match. Reference status, embedded support and previous mappings
give no candidate a preference. Each assigned molecule selects one isoform
globally, even if more than one overlapping gene has eligible isoforms.

The mode supports `--assignment-mode unique` only, which is already the
default. `fractional` is rejected. Existing neighboring mappings and `name2`
read memberships are ignored; do not pass `--read-to-isoform`.

## Unassigned reads and input errors

If a read has no candidate, it receives no count and appears once in
`*.unassigned_reads.tsv`. Reasons identify the first filter that leaves no
candidates: `no_candidate_locus`, `strand_mismatch`, `gene_mismatch` or
`no_exon_overlap`. See [report schemas](FORMATS.md#fixed-catalog-counting-outputs)
for their exact meanings and the accounting fields.

Identical duplicate read alignments count once per molecule. Conflicting
alignments under one read ID, empty read IDs, duplicate/empty isoform IDs,
malformed BED records and empty catalogs fail the run. These input errors are
not unassigned-read records. A valid catalog with no matching reads succeeds
with zero counts. In manifest mode, identity checks apply within each sample.

## Choose the right count mode

| Task | Command |
| --- | --- |
| Assign all input reads to given isoforms, without discovery | `count --assign-against-catalog` or `count-multi --assign-against-catalog` |
| Recount completed per-gene discovery outputs | `count --output-root ... --prefix ... --reference ...` |
| Rerun merge, count and description from completed per-gene outputs | `flow --count-only ...` |
| Count from an existing mapping or embedded memberships | Standalone `count` / `count-multi` without `--assign-against-catalog` |

`flow --count-only` requires completed discovery artifacts; it is not an entry
point for a new external catalog. For direct counting, `--reference` is
optional and unused, and `count --output-root` / `--prefix` cannot be combined
with `--assign-against-catalog`.

See the [CLI reference](CLI.md#trackcluster-count) for all count options and
the [pipeline tutorial](PIPELINE.md) for workflows that discover isoforms.
