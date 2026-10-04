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
- `polya.nanopolish.tsv`: a synthetic PASS tail estimate for `read_trunc` in
  `reads.bed`, using the standard Nanopolish ten-column layout.

All files are synthetic and intentionally small; they demonstrate file layout
and command wiring, not biological performance.

## Count given isoforms without discovery

Run these commands from the repository root or unpacked release directory,
using a build whose `count --help` includes `--assign-against-catalog`:

```bash
trackcluster count \
  --reads examples/reads.bed --isoform examples/ref.bed \
  --assign-against-catalog --out out/catalog-count/single/count.csv

trackcluster count-multi \
  --manifest examples/samples.tsv --isoform examples/ref.bed \
  --assign-against-catalog --out out/catalog-count/multi/count
```

The catalog contains `ref_a` and `ref_b` in gene `GENEA`. The single input read
selects `ref_a`, giving `count.csv`:

```csv
gene,isoform_id,count
GENEA,ref_a,1
GENEA,ref_b,0
```

Each manifest sample has one read, giving `count.isoform_counts.matrix.tsv`:

```tsv
gene	isoform_id	S1	S2
GENEA	ref_a	1	1
GENEA	ref_b	0	0
```

Both commands also emit selected mappings, unassigned-read reports and
assignment statistics. Multi-sample counting adds aggregate counts and
sample/group usage. These examples have no unassigned reads. In the source
checkout, see [`docs/COUNTING.md`](../docs/COUNTING.md) for the complete walkthrough
and [`docs/FORMATS.md`](../docs/FORMATS.md#fixed-catalog-counting-outputs) for schemas.

## Join poly(A) lengths to the selected assignments

After the single-sample count above:

```bash
trackcluster polya-aggregate \
  --nanopolish examples/polya.nanopolish.tsv --sample S1 \
  --isoforms examples/ref.bed \
  --read-to-isoform out/catalog-count/single/count.read_to_isoform.tsv \
  --out out/catalog-count/single/count
```

`count.isoform_polya.tsv` contains `polya_reads=1` and
`polya_median_nt=80.250000` for `ref_a`. The zero-read `ref_b` stays in the table
with `NA` length statistics. These are synthetic lengths, not measurements.

To include the same summary after discovery:

```bash
trackcluster flow \
  --reads examples/reads.bed --reference examples/ref.bed \
  --output-root out/polya --prefix sample \
  --polya-nanopolish examples/polya.nanopolish.tsv --polya-sample S1
```

See the [poly(A) guide](../docs/POLYA.md) for Dorado input, multiple samples,
deduplication and failure audits.
