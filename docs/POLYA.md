# Isoform-level poly(A) tail lengths

`polya-aggregate` joins Dorado's per-read `pt:i` estimates to the final
read-to-isoform mapping and summarizes tail length for each sample and isoform.
Use `polya_median_nt` to compare longer and shorter tails, alongside
`polya_reads`, the number of successful estimates contributing to that median.
The output preserves the full distribution as one row per assigned read.

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

## Summarize an existing single-sample assignment

For a discovery run:

```bash
trackcluster polya-aggregate \
  --bam reads.dorado.bam --sample S1 \
  --isoforms out/sample_isoform.bed \
  --read-to-isoform out/sample_read_to_isoform.unique.tsv \
  --out out/sample
```

Single-BAM mode joins the exact BAM query name to the mapping read ID. An
optional `--group control` supplies a group label. For fixed-catalog counting,
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

`sample` and `bam` are required; `group` is optional. BAM paths are relative to
this manifest. Each sample appears exactly once. The pooled assignment must
use `S1::original_read_id`, `S2::original_read_id`, and so on. The BAM files
retain their original, unprefixed query names. Samples with the same query
name stay separate.

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
As with other flow inputs, keep source BAMs and manifests outside `--output-root`.
Rerunning flow without poly(A) input removes previous poly(A) tables after
the replacement count/catalog outputs succeed.

## Outputs and statistics

| Output | Contents |
| --- | --- |
| `<prefix>.isoform_polya.tsv` | One row per sample and catalog isoform, including zero-read isoforms |
| `<prefix>.read_polya.tsv` | Each assigned read's isoform, raw `dorado_pt`, valid `polya_length_nt`, and status |
| `<prefix>.polya_qc.tsv` | Source BAM, available Dorado versions, record filters, read-join rate and failure counts |

The isoform table includes `assigned_reads`, `observed_reads`, `polya_reads`,
the four missing/failed categories, and `polya_fraction = polya_reads /
assigned_reads`. `observed_reads` includes failed estimates and missing tags,
but excludes assigned reads not found among the BAM's primary records.
`read_join_rate` in the QC table is `observed_reads / assigned_reads`.

Length summaries use only positive estimates: arithmetic mean, median,
25th/75th percentiles, minimum, maximum, and sample standard deviation
(`n - 1` denominator). Percentiles use linear interpolation at `(n - 1) * p`.
Length statistics are `NA` with no valid estimates; standard deviation is
also `NA` for a single estimate. A fraction with zero denominator is `NA`.

Secondary and supplementary records are ignored. Repeated assigned primary
records with the same estimate are deduplicated; conflicting estimates fail
instead of depending on BAM record order. Unassigned primary records appear
only in record-level QC and do not enter isoform denominators. A zero read-ID
join is visible as `read_join_rate = 0` and missing BAM reads.

Poly(A) counts describe the actual assigned molecules. They are not scaled
by flow's downsampling factors. Sort `polya_median_nt` within a sample to rank
isoforms by tail length, checking `polya_reads` and `polya_fraction` alongside
the ranking. The summaries are descriptive; no differential-tail significance
test or universal long/short threshold is applied.

The existing abundance CSV schemas are preserved. Join the new table to
single-sample counts by `gene,isoform_id`, or to the multi-sample usage table
by `sample,gene,isoform_id`.
