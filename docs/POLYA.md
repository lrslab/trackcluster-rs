# Isoform-level poly(A) tail lengths

`polya-aggregate` joins Dorado BAM `pt:i` estimates or Nanopolish `polya` TSV
results to the final read-to-isoform mapping and summarizes tail length for
each sample and isoform.
Use `polya_median_nt` to compare longer and shorter tails, alongside
`polya_reads`, the number of successful estimates contributing to that median.
The output preserves the full distribution as one row per assigned read.
For a runnable example using bundled synthetic data, see the
[packaged poly(A) example](../examples/README.md#join-polya-lengths-to-the-selected-assignments).

## Input from Dorado

Enable tail estimation during basecalling:

```bash
dorado basecaller MODEL reads.pod5 --estimate-poly-a > reads.dorado.bam
```

The BAM must retain the `pt:i` tag and the original query names. It can be
the unaligned Dorado BAM or an aligned BAM with those tags preserved. Tail
length is joined by read ID, so alignment strand, CIGAR and MAPQ are not used
to recalculate the estimate. The final mapping determines which molecules
contribute to each isoform. For cDNA, Dorado's poly(T) estimates are also
included under the poly(A) output column names.

The [Dorado poly(A) documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/)
defines these `pt:i` values:

| Value | Meaning | Contribution to tail-length statistics |
| --- | --- | --- |
| Positive integer | Successful estimated tail length, in nucleotides | Included once per assigned read |
| `-1` | Primer anchor not found | Excluded; counted as `anchor_not_found_reads` |
| `0` | Anchor found, length estimation failed | Excluded; counted as `estimation_failed_reads` |
| Tag absent | No estimate supplied | Excluded; counted as `missing_pt_reads` |

Noninteger tags and integers below `-1` fail validation. A nonempty primary
read stream with no `pt:i` tags fails with instructions to enable estimation.
Missing tags on individual reads are retained as missing data when other
primary reads carry `pt:i`.

## Input from Nanopolish

Use the output of `nanopolish polya`, as described in the
[Nanopolish poly(A) tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html).
The standard TSV has these ten columns:

```tsv
readname	contig	position	leader_start	adapter_start	polya_start	transcript_start	read_rate	polya_length	qc_tag
```

The importer uses `readname`, `polya_length`, and `qc_tag`. A header containing
those names can include other columns in any order. It also accepts the
standard ten-column output without a header, including the tutorial's
`grep 'PASS'` export. An empty export supplies no estimates. Use the original
unfiltered TSV to retain failed-QC records in the audits.

Only `qc_tag=PASS` contributes to length statistics. `polya_length` is a
floating-point estimate in nucleotides; fractional values are retained.
PASS lengths must be finite and non-negative. Zero is valid for Nanopolish:
its [estimator source](https://github.com/jts/nanopolish/blob/master/src/nanopolish_polya_estimator.cpp)
clamps negative estimates to zero. Dorado's zero instead remains a failure
sentinel under its separate caller contract.

Other QC tags, including `SUFFCLIP`, `ADAPTER`, `NOREGION`, and
`READ_FAILED_LOAD`, count as `qc_failed_reads`. Their raw length and QC tag
are retained in the read audit, while `polya_length_nt` is `NA`. An assigned
read absent from the TSV counts as `missing_nanopolish_reads`.

Nanopolish can emit several alignment rows for one read. Agreeing PASS
lengths count once. A PASS row takes precedence over failed-QC rows for the
same read; conflicting PASS lengths fail validation. When several rows have
the same accepted length, or all fail QC, the read audit selects the smallest
`(qc_tag, raw polya_length)` tuple lexically. This makes the result independent
of input row order. The QC table counts every supplied row and inventories
all QC tags, including discarded duplicates and unassigned reads.

## Summarize an existing single-sample assignment

For a discovery run:

```bash
trackcluster polya-aggregate \
  --bam reads.dorado.bam --sample S1 \
  --isoforms out/sample_isoform.bed \
  --read-to-isoform out/sample_read_to_isoform.unique.tsv \
  --out out/sample
```

For Nanopolish, replace `--bam` with `--nanopolish`:

```bash
trackcluster polya-aggregate \
  --nanopolish polya_results.tsv --sample S1 \
  --isoforms out/sample_isoform.bed \
  --read-to-isoform out/sample_read_to_isoform.unique.tsv \
  --out out/sample
```

Single-sample mode joins the exact BAM query name or Nanopolish `readname`
to the mapping read ID. An optional `--group control` supplies a group label.
For fixed-catalog counting,
use the selected `*.read_to_isoform.tsv` from
[`count --assign-against-catalog`](COUNTING.md) with the supplied catalog BED.
No new clustering or assignment is performed.

## Multiple samples

Create a tab-delimited poly(A) manifest:

```tsv
sample	bam	group
S1	S1.dorado.bam	control
S2	S2.dorado.bam	treated
```

For Nanopolish-only inputs, use `sample,nanopolish` and optional `group`
columns. A manifest can also include both source columns:

```tsv
sample	bam	nanopolish	group
S1	S1.dorado.bam	NA	control
S2	NA	S2.polya.tsv	treated
```

`sample` and at least one of `bam` or `nanopolish` are required. Each sample
appears once and supplies exactly one source path; blank or `NA` means no
path. Paths are relative to the manifest. Every output row records its
`caller`, and summaries remain separate per sample. Select one caller per
sample to avoid counting the same molecules twice.

The pooled assignment must use `S1::original_read_id`,
`S2::original_read_id`, and so on. Source files retain original, unprefixed
read names. Samples with the same read name stay separate.

```bash
trackcluster polya-aggregate \
  --polya-manifest polya.tsv \
  --isoforms out/pooled_isoform.bed \
  --read-to-isoform out/pooled_read_to_isoform.unique.tsv \
  --out out/pooled
```

The manifest must cover every sample appearing in the mapping. Unknown
samples, unknown isoforms and a read assigned to multiple distinct isoforms
fail validation. Repeated copies of the same read/isoform pair are counted
once. Flow's gene-local `.unique.tsv` can contain assignments to multiple
genes; poly(A) aggregation requires those assignments to be globally
unambiguous, just like modification aggregation.

## Include summaries in `flow`

For a single sample, add `--polya-bam` and optionally `--polya-sample`:

```bash
trackcluster flow \
  --reads reads.bed --reference ref.bed \
  --output-root out --prefix sample \
  --polya-bam reads.dorado.bam --polya-sample S1
```

Use `--polya-nanopolish polya_results.tsv` for Nanopolish instead. The two
source flags are mutually exclusive. If `--polya-sample` is omitted, the
sample label defaults to the source filename stem.

For pooled discovery, add `--polya-manifest`:

```bash
trackcluster flow \
  --manifest samples.tsv --reference ref.bed \
  --output-root out --prefix pooled \
  --polya-manifest polya.tsv
```

The poly(A) manifest must contain exactly the samples in the reads manifest.
Groups are inherited from the reads manifest; a conflicting declared group
fails validation. Poly(A) processing requires `--assignment-mode unique` and
runs after final assignment. These options also work with `--count-only`.
As with other flow inputs, keep source BAMs, TSVs and manifests outside
`--output-root`.
After argument/input preflight, flow removes previous poly(A) tables before
starting the core pipeline. A failed rerun therefore cannot leave old tail
summaries beside replaced catalogs or assignments. A successful rerun with
poly(A) input regenerates all three tables; a rerun without it leaves them absent.

## Outputs and statistics

| Output | Contents |
| --- | --- |
| `<prefix>.isoform_polya.tsv` | One row per sample and catalog isoform, including zero-read isoforms |
| `<prefix>.read_polya.tsv` | Each assigned read's isoform, caller, raw caller estimate/QC, valid `polya_length_nt`, and status |
| `<prefix>.polya_qc.tsv` | Caller and input path, source record/row counts, read-join rate and failure counts |

The isoform table includes `assigned_reads`, `observed_reads`, `polya_reads`,
the missing/failed categories, and `polya_fraction = polya_reads /
assigned_reads`. `observed_reads` includes failed estimates and missing tags,
but excludes assigned reads not found among primary BAM records or in the
Nanopolish TSV.
`read_join_rate` in the QC table is `observed_reads / assigned_reads`.

Length summaries use accepted estimates under each caller's rules: arithmetic
mean, median, 25th/75th percentiles, minimum, maximum, and sample standard deviation
(`n - 1` denominator). Percentiles use linear interpolation at `(n - 1) * p`.
Length statistics are `NA` with no valid estimates; standard deviation is
also `NA` for a single estimate. A fraction with zero denominator is `NA`.

Dorado secondary and supplementary records are ignored. Repeated assigned
primary records with the same estimate are deduplicated; conflicting estimates fail
instead of depending on BAM record order. Unassigned primary records appear
only in record-level QC and do not enter isoform denominators. Unassigned
Nanopolish rows follow the same rule. A zero read-ID join is visible as
`read_join_rate = 0` and missing source reads.

The original Dorado output columns keep their positions. Nanopolish support
appends `caller` and source-specific audit columns. Dorado-only counters are
zero in Nanopolish molecule summaries and `NA` in its BAM record inventory;
Nanopolish row inventories are `NA` for Dorado. See
[formats](FORMATS.md#isoform-level-polya-formats) for all columns.

Poly(A) counts describe the actual assigned molecules. They are not scaled
by flow's downsampling factors. Sort `polya_median_nt` within a sample to rank
isoforms by tail length, checking `polya_reads` and `polya_fraction` alongside
the ranking. The summaries are descriptive; no differential-tail significance
test or universal long/short threshold is applied.

The existing abundance CSV schemas are preserved. Join the new table to
single-sample counts by `gene,isoform_id`, or to the multi-sample usage table
by `sample,gene,isoform_id`.
