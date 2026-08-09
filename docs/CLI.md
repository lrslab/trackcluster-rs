# CLI reference

This is a reference for the Rust binaries shipped in this repo:
- `trackcluster` (main CLI)
- `clusterj_batch` (gene-batched runner)

## Common conventions
- Core clustering/counting inputs are BED12 / bigGenePred-like. `bam2bigg` and
  `gff2bigg` create that persisted representation from BAM or GFF3/GTF inputs.
- Coordinates are 0-based, half-open.
- Most commands are strand-aware (the legacy TrackCluster behavior).
- **Gene name field**: stored in extra field index `5` (value like `GENE1` or `GENE1||GENE2`; `none` means unassigned).
- **Gene IDs versus path keys**: reports, `<prefix>_gene.txt`,
  `clusterj_batch --gene-list`, and every `--downsample-gene` value use the
  biological gene ID. Per-gene directories and filenames instead use its
  encoded `<gene-path-key>` from `<prefix>_gene_paths.tsv`; do not pass a
  folder name or encoded key where a biological ID is requested. See
  [Biological gene IDs and filesystem keys](FORMATS.md#biological-gene-ids-and-filesystem-keys).
- **Subread payload (`name2`)**: stored in extra field index `0` for isoforms produced by `clusterj`/`cluster`.
  - Default (`--name2-mode coverage`): `|3.5` (no read IDs)
  - `--name2-mode full`: `tc_name2_v1:read1,read2,|3.5` (read IDs are percent encoded)
  - `--name2-mode none`: `none` (no payload)
  - Use the `*_read_to_isoform.tsv` mapping for counting when read IDs are not embedded (and it is auto-discovered by `count` / `count-multi` when present next to the isoform BED).
- Converters, exports, validation reports, `addgene`, `desc`, single-gene
  clustering, and standalone count outputs reject destinations that alias an
  input or one another (including existing hard links) before writing. Each
  destination is replaced atomically after its complete output has been
  written and synced.
- `flow`, `preparedir`, and `count --output-root` own the complete
  `--output-root` tree. Keep the reference, reads, sample manifest, and every
  manifest reads file outside that tree. Before creating or replacing any
  artifact, these commands reject an external input beneath the tree, any
  existing alias to an input, and all symlink or multiply-linked file entries
  in the tree. Use ordinary files/directories in a dedicated output root.

## `trackcluster`

### `trackcluster validate-bed`

Validate a BED12/bigGenePred file and print tab-delimited summary counts.
Strict validation is the default.

Example:

```bash
trackcluster validate-bed --input reads.bed
```

Flags:

- `--input/-i` (required): input BED12/bigGenePred file.
- `--lenient`: accept only the documented legacy repairs below. Every changed
  field is reported; this option does not write a normalized BED.
- `--report <PATH>`: write a tab-delimited
  `trackcluster-bed-validation-v1` summary containing `mode`, `records`,
  `exons`, `repairs`, `normalized_records`, and `errors`. The report is also
  written when record-level validation errors make the command exit nonzero.

Both modes scan the complete file. Standard output contains `records`,
`exons`, `repairs`, `normalized_records`, and `errors`. Record errors are
printed to standard error as `validation_error` rows; any such error makes the
command fail after scanning. In lenient mode, each accepted change is also
printed to standard error as a `repair` row with source path, one-based line,
field, original value, replacement value, and reason.

Lenient mode repairs exactly this subset:

- a non-integer BED score becomes `0`;
- an integer BED score above the BED maximum of `1000` becomes `1000`;
- empty **interior** values in `blockSizes` or `blockStarts` are removed (a
  conventional single trailing comma is already accepted in either mode);
- a block whose start is still before `chromEnd` but whose end exceeds it is
  shortened to end at `chromEnd`.

All other parse and geometry violations remain errors, including malformed
coordinates or strand, invalid thick spans, block-count/list-length mismatch,
out-of-order, overlapping, zero-length, or wholly out-of-span blocks, and a
transcript span that does not match its outer exon bounds.

### `trackcluster bam2bigg`

Convert genome-aligned BAM records to TrackCluster bigGenePred-compatible
BED12+8 text. This command reads BAM, not SAM or CRAM.

Example:

```bash
trackcluster bam2bigg \
  --bamfile alignments.bam \
  --out reads.bed \
  --score 30 \
  --group sample_A
```

Flags and defaults:

- `--bamfile/-b` (required; alias `--input`): input BAM.
- `--out/-o` (alias `--output`, default `bigg.bed`): output BED12+8 text.
- `--score/-s` (alias `--min-mapq`, default `30`): minimum retained MAPQ.
- `--group/-g`: sample/group value written to extra field index `6`. If
  omitted, the BAM filename stem is used; an empty label becomes `none`.
- `--include-secondary`: retain records carrying flag `0x100` (excluded by
  default).
- `--include-supplementary`: retain records carrying flag `0x800` (excluded by
  default).

Unmapped records are always skipped. Missing MAPQ is treated as zero for
filtering and for the BED score. For each retained alignment, only CIGAR `N`
starts a new exon: reference-consuming matches, mismatches, and deletions stay
inside the current block, while insertions and clipping do not consume the
reference. Every emitted block must contain at least one `M`, `=`, or `X`;
deletion-only blocks and alignments beyond the BAM header reference length are
rejected. Reverse-complement flag `0x10` sets the minus strand. Multiple
retained records with the same query name remain separate BED alignment
instances; downstream counting applies the molecule-ID policy.

The BED score is the record MAPQ. The converter writes `name2=none`,
`type=nanopore_read`, `geneName=none`, the selected sample/group, no CDS, and
one `-1` exon frame per block. See [Converter-produced BED12+8](FORMATS.md#converter-produced-bed128)
for the complete field contract.

On success, stderr reports total decoded records, emitted records, and counts
skipped as unmapped, secondary, supplementary, or below MAPQ. Conversion is
streamed and the destination is atomically published only after the complete
BAM succeeds.

### `trackcluster gff2bigg`

Convert GFF3 or GTF transcript annotations to a validated, deterministically
sorted TrackCluster reference BED12+8 catalog.

Example:

```bash
trackcluster gff2bigg \
  --gff annotation.gff3 \
  --out reference.bed \
  --input-format gff3 \
  --key ID
```

Flags and defaults:

- `--gff/-i` (required; alias `--input`): input GFF3 or GTF annotation.
- `--out/-o` (alias `--output`, default `bigg.bed`): output BED12+8 text.
- `--key/-k` (alias `--gene-key`, default `ID`): attribute on GFF3 `gene`
  rows used as the output gene label. Graph relationships continue to use
  canonical `ID`/`Parent`; if the requested attribute is absent, `ID` or
  `gene_id` is used.
- `--input-format` (default `auto`): `auto`, `gff3`, or `gtf`. Auto mode
  detects `key=value` versus `key "value"` attribute syntax per data row.

Every data row must have exactly nine tab-separated columns and valid positive,
one-based inclusive coordinates. Blank/comment lines are ignored, and parsing
stops at `##FASTA`. GFF3 attributes are percent-decoded; GTF quoted attributes
are unescaped. `?` and `.` strands map to unknown.

Only `exon` rows define BED blocks. GFF3 groups them by comma-separated
`Parent`, and every parent must resolve to a declared non-gene feature with a
unique `ID`; GTF groups exon-only or full models by `transcript_id`.
Transcript rows and exon-level `gene_id` attributes provide gene hints.
Multiple distinct hints are sorted and joined with `||`; missing gene
annotation is written as `none`. Duplicate exon intervals are collapsed, while
duplicate graph identities, cross-contig or conflicting-strand parent/child
features, parent spans that do not contain their exons, overlapping blocks, or
otherwise invalid transcript models fail conversion. BED transcript spans are
derived from the outer exon bounds, independent of input order.

CDS, UTR, phase, feature score, and declared transcript span are not projected
into the current output. Each model therefore has score `100`, `itemRgb=0`,
`thickStart=thickEnd=0`, `name2=none`, `cdsStartStat=cdsEndStat=none`,
all exon frames `-1`, `type=isoform_anno`, and no sample/group. Output
transcript IDs are validated as reference IDs and sorted by structure before
writing. The output is atomically published only after the entire annotation
passes. On success, stderr reports the number of transcripts written.

### `trackcluster export`

Export a strict BED12/bigGenePred transcript catalog to one or more
standards-oriented text formats. At least one output flag is required, and all
requested outputs are generated from the same parsed catalog.

```bash
trackcluster export \
  --input catalog.bed \
  --gtf catalog.gtf \
  --gff3 catalog.gff3 \
  --sqanti-input catalog.sqanti.tsv
```

Flags:

- `--input/-i` (required): input transcript catalog.
- `--gtf <PATH>`: GTF 2.2 transcript and exon features with `gene_id` and
  `transcript_id` attributes.
- `--gff3 <PATH>`: GFF3 `mRNA` and exon features with percent-encoded IDs.
- `--sqanti-input <PATH>`: `trackcluster-sqanti-input-v1` audit TSV with
  `isoform_id`, `gene_id`, chromosome, strand, summed exon length, and exon
  count.

GTF and GFF3 coordinates are converted from BED's zero-based, half-open form
to one-based, closed coordinates. `--sqanti-input` is an identity/geometry
audit table, not a SQANTI3 structural-classification result; the exported GTF
is the file to supply to SQANTI3 for that external analysis. See
[`INTERCHANGE.md`](INTERCHANGE.md#transcript-export) for the interchange
contract.

### `trackcluster flow`
Run the full pipeline as a single command:
1) `preparedir` (dedup + gene assignment into per-gene folders)
2) internal per-gene clustering (`clusterj` by default, or overlap-mode `cluster` with `--cluster-mode cluster`)
3) merge per-gene outputs into `<prefix>_isoform.bed` and `<prefix>_unused.bed`
4) `count` and `desc` on the merged isoforms
5) when `--manifest` is used: run per-sample `count-multi` outputs from the pooled isoforms
6) optionally aggregate normalized modification observations using the final unique mapping

Key flags:
- `--cluster-mode`: `clusterj` (default) or `cluster` (overlap/intersection mode)
- `--reads/-s`: reads BED (single-sample mode; mutually exclusive with `--manifest`)
- `--manifest`: sample manifest TSV for pooled clustering + per-sample usage
- `--reference/-r`: reference BED
- `--output-root/-o`: output directory (created if missing)
- `--prefix`: output prefix (used for merged outputs like `<prefix>_isoform.bed`)
- `--threads/-t`: number of worker threads (parallel across genes)
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (flow
  default: `-1`, no SW/SL 5' signal in either cluster mode). Pass a
  non-negative cutoff such as `11` only when BED score is valid SL/SW 5'
  evidence. In overlap mode, the second pass protects a short read only when
  its score is at or above the cutoff; with `-1`, ordinary short-read merging
  still runs. The separate direct `trackcluster cluster` command defaults to
  `11` for legacy compatibility.
- `--batch-size`, `--batch-rounds`: bounds for very large genes; in overlap mode these control iterative pre-merging rounds before the final two-pass overlap clustering
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size; mapping TSVs are used for counting)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls used by junction-mode clustering.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): junction-mode SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): junction-mode same-junction 3' terminal retention controls.
- Junction-mode clustering retains supported same-junction 3' terminal
  clusters as isoforms. The minus-strand 3' end is on the lower-coordinate
  side, but an early stop has a higher `tx_start` than the full-length isoform.
- `--overlap-cutoff1`, `--overlap-cutoff2`, `--overlap-intron-weight`: overlap-mode controls used when `--cluster-mode cluster`
- `--prepare-fraction-read`, `--prepare-fraction-ref`: overlap thresholds for gene assignment
- `--assignment-mode`: final counting mode, `unique` (default) or `fractional`; unique mode expands candidates against the isoform catalog before choosing the closest compatible isoform, including retained 3' early-stop isoforms.
- `--unique-assignment-junction-offset` (default: `15`): maximum per-boundary difference for the ordered one-to-one intron matcher used by unique assignment.
- `--mod-manifest`: optional normalized modification manifest. It requires
  manifest mode and `--assignment-mode unique`; modification aggregation runs
  only after clustering/count artifacts have been published.
- `--mod-analysis-threshold ASSAY=VALUE`: required once per assay when
  `--mod-manifest` is set. `--mod-min-callable`,
  `--mod-min-read-join-rate`, and `--mod-allow-low-join` mirror
  `mod-aggregate`.
- `--mod-contrasts`: optional explicit contrast specification processed after
  the modification design table is written.
- `--emit-pooled-reads`: when using `--manifest`, also write `<prefix>_pooled_reads.bed`
- `--max-reads-per-gene` (default: `50000`; set `0` to disable),
  `--downsample-gene <BIOLOGICAL_GENE_ID>` (repeatable),
  `--downsample-seed`: per-gene downsampling (writes
  `clusterj_batch_downsample.tsv` or `cluster_batch_downsample.tsv` and scales
  counts). Final flow/count-only processing rejects a molecule found in more
  than one selected gene when any of those genes was actually downsampled;
  use `--max-reads-per-gene 0` or exclude every affected gene from
  downsampling rather than publish a biased cross-gene count.
- `--force`: rebuild every per-gene result (otherwise only exact, hash-verified `run.json` results are reused); cannot be combined with `--count-only`
- `--invalid-read-policy`: `skip` (default) excludes only an individual malformed read track or a read track with an empty ID and records it in a rejected-read TSV; `fail` restores strict read parsing and fails the enclosing stage on the first invalid read track.
- `--strict-gene-errors`: restore all-or-nothing behavior. By default, gene-local failures are written to the batch error report, excluded from downstream inputs, and verified genes continue through merge/count/description.
- `--count-only`: reuse existing per-gene outputs and run only merge, unique/fractional count, multi-sample count, and desc outputs
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)

Outputs (under `--output-root`):
- `<prefix>_gene.txt` (biological gene IDs selected during preparation)
- `<prefix>_dedup.bed`
- `<prefix>_novel.bed`
- `<prefix>_gene_paths.tsv` (versioned biological gene ID to filesystem-key mapping)
- `<prefix>_isoform.bed`
- `<prefix>_unused.bed`
- `<prefix>_read_to_isoform.tsv`
- `<prefix>_read_to_isoform.unique.tsv` (unique assignment mode only; exact selected mapping used for final counts)
- `<prefix>_unique_assignment.provenance.tsv` (unique mode; records the effective tolerance and one-to-one/no-collapse policy)
- `<prefix>_isoform_count.csv`
- `<prefix>_desc.txt`, `<prefix>_class4.txt`, `<prefix>_class12.txt`, `<prefix>_fusion.txt`
- `clusterj_batch_summary.txt` / `cluster_batch_summary.txt`; the matching
  `*_errors.txt` exists only when the run records errors
- `clusterj_batch_gene_paths.tsv` / `cluster_batch_gene_paths.tsv`
- `clusterj_batch_downsample.tsv` / `cluster_batch_downsample.tsv` (only when downsampling occurs)
- `<gene-path-key>/run.json` (versioned per-gene request, tool, seed, and output-integrity manifest)
- `<gene-path-key>/rejected_reads.tsv` (read tracks rejected while loading that gene; header-only when none are rejected)
- `<prefix>_rejected_reads.tsv` (read tracks rejected during preparation, before gene assignment)
- `<prefix>_pooled_reads.bed` (manifest mode + `--emit-pooled-reads`)
- `<prefix>.isoform_count.csv` (manifest mode only; aggregate counts derived from the per-sample matrix)
- `<prefix>.isoform_usage.long.tsv` (manifest mode only)
- `<prefix>.isoform_counts.matrix.tsv` (manifest mode only)
- `<prefix>.isoform_usage.group.tsv` (manifest mode only; when at least one sample has a non-empty `group`)
- `<prefix>.mod_join_qc.tsv`, `<prefix>.mod_site_join_qc.tsv`,
  `<prefix>.isoform_mod_sites.tsv`, and `<prefix>.isoform_mod_design.tsv`
  (manifest mode + `--mod-manifest`)
- `<prefix>.isoform_mod_contrasts.tsv` (`--mod-contrasts`)

#### Platform presets and SL/no-SL data

`--platform-preset` changes only the junction correction, SL 5', and same-junction 3' merge/protection defaults. The junction-mode `--sw-score` default remains `-1` for all presets; pass `--sw-score 11` when the BED score should be used as valid SL/SW 5' evidence.
Junction min support is weighted site support: read sites contribute `1`, reference sites contribute `5`.

| Preset | Recommended use | Junction correction offset | Junction min support | SL partial 5' offset | SL same-junction 5' offset | SL 5' cluster offset | SL min support | 3' same-junction offset | 3' cluster offset | 3' min support |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `generic` | Default; conservative setting for mixed or unknown platforms | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |
| `rna002` | RNA002 direct RNA reads; more tolerant of junction and 5' end wobble | `15` | `5` | `20` | `25` | `20` | `2` | `50` | `15` | `5` |
| `rna004` | RNA004 direct RNA reads; use the conservative/default cutoffs | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |

SL information is optional. With the default `--sw-score -1`, all reads are treated as non-SL-supported and ordinary junction correction/truncation merging still runs. When `--sw-score` is non-negative, only reads whose BED score is greater than the cutoff receive SL-cluster protection as alternative 5' isoforms. For no-SL datasets, including human data where BED score is MAPQ or another non-SL value, keep `--sw-score -1`. For SL/SW-scored datasets, pass an explicit cutoff such as `--sw-score 11` together with the appropriate platform preset.

Same-junction 3' terminal clusters with nearby read support are retained as
isoforms independently of SL evidence. By default, protection requires at
least `5` reads within the active junction correction offset and a 3' end more
than `50` bp from the merge target. The rule is strand-aware: on minus-strand
transcripts the 3' end is the lower-coordinate side, but a 3' early stop
truncates that side and therefore has a higher `tx_start` (a higher genomic
terminal boundary) than the full-length isoform.

Example:
```bash
trackcluster flow \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --threads 8
```

Overlap-mode example:
```bash
trackcluster flow \
  --cluster-mode cluster \
  --reads reads.bed \
  --reference ref.bed \
  --output-root out \
  --prefix sample \
  --sw-score 11 \
  --batch-size 500 \
  --batch-rounds 100
```

Manifest example:
```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

Manifest-mode example:
```bash
trackcluster flow \
  --manifest samples.tsv \
  --reference ref.bed \
  --output-root out \
  --prefix pooled
```

Manifest mode with normalized modification observations:
```bash
trackcluster flow \
  --manifest samples.tsv \
  --reference ref.bed \
  --output-root out \
  --prefix pooled \
  --mod-manifest mod_samples.tsv \
  --mod-analysis-threshold dorado_rna004_m6a=0.5
```

If modification aggregation or contrast calculation fails, `flow` returns an
error but leaves the already completed clustering and count artifacts intact.

Count-only rerun example:
```bash
trackcluster flow \
  --count-only \
  --reference ref.bed \
  --output-root out \
  --prefix sample
```

`--count-only` expects completed per-gene output folders under `--output-root`. It uses the prefix-scoped `<prefix>_gene.txt` and `<prefix>_gene_paths.tsv` metadata when present. For a standalone `clusterj_batch` tree that has no prefix-scoped metadata, it uses the versioned `clusterj_batch_gene_paths.tsv` (or `cluster_batch_gene_paths.tsv` for overlap mode); it never infers the active gene set from stale directories. Per-gene directories and filenames use encoded path keys, while reports retain the biological gene ID. Gene IDs are limited to 4096 UTF-8 bytes and reject path separators, absolute paths, `.`/`..`, and control characters. Before publishing any merged output, count-only verifies every selected gene's `run.json`, prepared-input hashes, cluster mode and tool identity, and all output hashes, sizes, and record counts. Missing, legacy, modified, or incomplete per-gene results fail the command; rerun the producing flow or batch command to rebuild them. In unique assignment mode, it selects reads directly from each verified gene folder using the key-based nano BED, per-gene isoform BED, and read-to-isoform TSV. Add `--manifest samples.tsv` to a count-only rerun when you need manifest-mode `*.isoform_usage.*` outputs regenerated.

### `trackcluster mod-import-m6anet`
Normalize m6Anet RNA002 read probabilities and project exact transcript
coordinates through a matching GTF/GFF annotation:

```bash
trackcluster mod-import-m6anet \
  --sample S1 \
  --assay-id m6anet_hct116_rna002 \
  --model-id HCT116_RNA002 \
  --indiv data.indiv_proba.csv.gz \
  --data-info data.info \
  --site-proba data.site_proba.csv.gz \
  --read-map read_index_to_read_id.tsv \
  --reference matching_annotation.gtf \
  --min-reads 20 \
  --out S1.mod
```

`--read-map` has the exact header `read_index<TAB>read_id`. `read_index` is an
opaque source token rather than an integer; this preserves m6Anet multi-input
identifiers such as `966210_0`. Raw read IDs are prefixed as
`<sample>::<read_id>`; already prefixed IDs must match `--sample`.
Transcript IDs, including versions, must match the annotation exactly.
`data.site_proba` is optional QC only and never supplies read observations.
`--candidate-rule` is a 3–31 base odd-length IUPAC motif (default `DRACH`), is
case-normalized, and must contain `A` at its center. When `data.site_proba` is
present, each reported concrete k-mer must match this motif.
Before parsing a genome-wide GTF/GFF, the importer scans the m6Anet tables for
the required transcript IDs and materializes only those transcripts (plus
same-stable-ID versions needed for a precise mismatch diagnostic).
`projection_transcripts_loaded` records the resulting catalog size in import
QC.
Known read-threshold presets are `0.033379376` for `HCT116_RNA002` and
`0.0032978046219796` for `arabidopsis_RNA002`; override the source cross-check
with `--read-probability-threshold`. The final hard-call threshold is still an
explicit `mod-aggregate --analysis-threshold` setting.

### `trackcluster mod-import-dorado`
Decode one modification code from a primary genome-aligned Dorado/modBAM:

```bash
trackcluster mod-import-dorado \
  --sample S1 \
  --assay-id dorado_rna004_m6a \
  --bam calls.aligned.bam \
  --mod-code A+a \
  --model-id rna004_m6a_model \
  --chemistry RNA004 \
  --caller-version 2.0.0 \
  --candidate-rule all-target-canonical-bases \
  --source-emission-threshold 0.05 \
  --out S1.mod
```

The importer validates `MM`/`ML`/`MN`, consumes every ML value before target
filtering, projects `M/=/X/I/D/N/S/H/P`, and excludes secondary,
supplementary, and unmapped records. The candidate rule may be
`all-target-canonical-bases` or an odd-length centered IUPAC DNA motif of 3–31
bases, such as `DRACH`; the motif center must be the target canonical base.
Motifs are matched on the as-sequenced read orientation (including reverse BAM
records), ambiguity in a read base does not create a candidate, and an explicit
target call outside the declared motif is rejected. Motif observations carry
the normalized motif in `context`.

`.` or omitted skip flags become `implicit_below_emission_threshold`; `?`
becomes `unknown`. An emission threshold is required before low-probability
implicit rows can be represented by valid assay metadata.

The safe default `--question-mark-policy unknown` follows the SAM
specification. For a known threshold-sparse Dorado source that uses `?` for
threshold omissions, the user may explicitly pass:

```bash
--candidate-rule all-target-canonical-bases \
--source-emission-threshold 0.05 \
--question-mark-policy below-emission-threshold
```

That override requires a positive source threshold. With a motif rule, only
read bases matching that motif are expanded as omitted candidates; other
canonical bases are outside the declared universe. `--invalid-record-policy
skip` records malformed primary records (including candidate-rule mismatches)
and marks the candidate universe incomplete; fail is default.

Both importers atomically write `<out>.observations.tsv`, `<out>.assay.json`,
and `<out>.import_qc.tsv`. The current modBAM importer materializes and sorts
observations in memory; use a representative aligned subset for initial
genome-wide validation rather than a multi-million-read BAM.
See [MODIFICATION_VALIDATION.md](MODIFICATION_VALIDATION.md) for pinned m6Anet
and ONT public-data checks.

### `trackcluster mod-subsample`
Create synchronized low-coverage pseudo-sample inputs from one high-coverage
sample after final unique assignment:

```bash
trackcluster mod-subsample \
  --manifest samples.tsv \
  --read-to-isoform out/pooled_read_to_isoform.unique.tsv \
  --mod-manifest mod_samples.tsv \
  --source-sample UHRR \
  --sample-prefix UHRR_low \
  --replicates 4 \
  --reads-per-sample 5000 \
  --mode disjoint \
  --seed 17 \
  --out-dir UHRR_low_inputs
```

The sampling unit is a source read molecule from the source sample's reads BED,
not an observation row or an already assigned/read-observation intersection.
Selected reads are synchronously materialized in reads BED, unique assignment,
normalized observations, and coverage BAM outputs. Selected observations
without an assignment are retained, so pseudo-sample join QC is not
artificially raised to 1. BAM records preserve their original order and
auxiliary tags, including `MM`/`ML`/`MN`.

`disjoint` is the default and assigns a source molecule to at most one
pseudo-sample. `independent` draws each pseudo-sample without replacement but
allows overlap between pseudo-samples. Both use versioned SHA-256 ranking and
are independent of BED and assignment input order. Normalized observations
must be in the canonical order emitted by the importers so they can be split
in one streaming pass.

The command publishes a complete new directory atomically. Its `samples.tsv`,
`mod_samples.tsv`, and `read_to_isoform.unique.tsv` can be passed directly to
`mod-aggregate`. Generated sample groups are intentionally empty; the original
group is retained only in `sample_provenance.tsv` and
`subsample_provenance.json`. These samples measure technical coverage
stability. They are not biological replicates, must not be used to increase
condition/interaction sample size, and modification counts are never
multiplied by a downsampling scale factor.

### `trackcluster mod-aggregate`
Join normalized observations to a final unique mapping:

```bash
trackcluster mod-aggregate \
  --manifest samples.tsv \
  --isoforms out/pooled_isoform.bed \
  --read-to-isoform out/pooled_read_to_isoform.unique.tsv \
  --mod-manifest mod_samples.tsv \
  --reference-fasta GRCh38.fa \
  --analysis-threshold dorado_rna004_m6a=0.5 \
  --out out/pooled
```

The command emits sample-level and genomic-site-level join QC, a complete
sample/isoform/site audit table, and a comparison-ready integer-count design
table. It rejects fractional or conflicting assignments, incompatible assay
metadata, low read-ID join rates, and analysis thresholds below a source
emission threshold. `--reference-fasta` requires an adjacent samtools-compatible
`.fai`; strand-oriented canonical-base mismatches remain in audit output as
`reference_base_mismatch` but are not callable. When a modification
manifest row supplies `coverage_bam`, mapped primary alignments are joined by
sample-tagged query name and CIGAR `M`/`=`/`X` blocks provide exact genomic-base
coverage. Every uniquely assigned molecule for that sample must be present;
duplicate primary names or mismatched chromosome/strand assignments fail.
Without a coverage BAM, `n_covering=NA` and `coverage_basis=unavailable`.
The global and site-local join gates use the same
`--min-read-join-rate`. `--allow-low-join` retains failures as ineligible rows;
it never turns them into eligible evidence.

### `trackcluster mod-site-summary`
Reduce one or more complete site tables to a deterministic genomic-site
inventory:

```bash
trackcluster mod-site-summary \
  --sites GAPDH.isoform_mod_sites.tsv \
  --sites ACTB.isoform_mod_sites.tsv \
  --out out/all_genes
```

The command streams each strict 27-column input and writes
`<out>.mod_site_summary.tsv`. It reports how many catalog isoforms are present,
eligible, absent, context dependent, reference-base mismatched, or ineligible
for other reasons at each `(assay, sample, gene, site, mod_code)`. A
`shared_eligible` state requires at least two eligible isoforms. Input row keys
must be unique across all supplied files.

### `trackcluster mod-contrast`
Calculate effect-only shared-site contrasts from an explicit nine-column
contrast specification:

```bash
trackcluster mod-contrast \
  --design out/pooled.isoform_mod_design.tsv \
  --contrasts contrasts.tsv \
  --out out/pooled
```

Supported types are `isoform_effect`, `condition_effect`, and
`isoform_condition_interaction`. V1 does not claim replicate-level inference:
`p_value` and `q_value` are always `NA`. Exact schemas and denominator rules
are in [FORMATS.md](FORMATS.md#isoform-level-modification-formats).

### `trackcluster preparedir`
Split one reads BED + one reference BED into per-gene folders.

`<prefix>_gene.txt` is published last as the prepared-generation commit marker.
If a late output error leaves it empty, rerun preparation; do not cluster or
count the incomplete tree.

When to use:
- For manual gene-batched mode before running `clusterj_batch`.
- Overlap-mode batching is not exposed as a separate `cluster_batch` binary; use `trackcluster flow --cluster-mode cluster` for that path.

Key flags:
- `--reads/-s`: reads BED
- `--reference/-r`: reference BED
- `--output-root/-o`: directory to create (gene folders written here)
- `--prefix`: prefix for summary outputs
- `--fraction-read`, `--fraction-ref`: overlap thresholds for gene assignment
- `--invalid-read-policy`: `skip` (default) omits only malformed/empty-ID read tracks and writes `<prefix>_rejected_reads.tsv`; `fail` stops preparation on the first such track.

Example:
```bash
trackcluster preparedir \
  --reads reads.bed \
  --reference ref.bed \
  --output-root tracktest \
  --prefix sample
```

### `trackcluster addgene`
Assign gene names to reads by overlap with reference and write `reads_gene.bed`.

Notes:
- Uses default overlap thresholds (same as `preparedir` defaults).
- If you need custom thresholds, use `preparedir` (it exposes them via flags).

Example:
```bash
trackcluster addgene \
  --reads reads.bed \
  --reference ref.bed \
  --out reads_gene.bed
```

### `trackcluster clusterj`
Junction-chain clustering (fast mode).

Outputs:
If `--out isoform.bed`, this command writes:
- `isoform.bed`: isoform BED
- `isoform.read_to_isoform.tsv`: read -> isoform mapping
- `isoform.unused.bed`: rare/filtered reads
- `isoform.rejected_reads.tsv`: malformed/empty-ID read tracks excluded before clustering (header-only when none)

Key flags:
- `--reads/-s`, `--reference/-r`, `--out/-o`
- `--threads/-t`: number of worker threads
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (default: `-1`, no SW/SL 5' signal; pass a non-negative cutoff such as `11` only when BED score is valid SL/SW 5' evidence)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls. This offset is separate from the SL/5' and 3' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): same-junction 3' terminal retention controls
- `--invalid-read-policy`: `skip` (default) excludes only the invalid read track; `fail` restores strict parsing

Performance note:
- 5' truncation collapsing is implemented via a junction-suffix index to avoid quadratic scans on large loci; `--batch-size`/`--batch-rounds` remain useful as hard caps for extreme genes.
- SL-supported reads with enough nearby 5' support can be protected as alternative isoforms; singleton likely-degradation reads can still merge into compatible longer/reference tracks.
- Supported same-junction 3' terminal clusters are retained as isoforms and remain compatible with unique counting, which assigns reads to the closest terminal structure.

Example:
```bash
trackcluster clusterj \
  --reads reads.bed \
  --reference ref.bed \
  --out isoform.bed \
  --threads 16
```

### `trackcluster cluster`
Overlap-based clustering (slower, more permissive; intended to mimic the legacy two-round TrackCluster overlap/intersection path).

Outputs are the same as `clusterj`:
- isoform BED, mapping TSV, unused BED, and rejected-read TSV (derived from the `--out` prefix).

Key flags:
- `--reads/-s`, `--reference/-r`, `--out/-o`
- `--threads/-t`: number of worker threads
- `--batch-size`, `--batch-rounds`: optional overlap batching for large loci (`--batch-size 0` disables intermediate batching)
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection in pass 2 (default: `11`; set to `-1` to treat reads as having no SW 5' signal). In pass 2, a short read is protected only when its score is at or above the cutoff; with `-1`, ordinary short-read merging still runs.
- `--cutoff1`, `--cutoff2`: overlap pass 1 / pass 2 cutoffs (default: `0.05`, `0.01`)
- `--intron-weight`: intron contribution to the combined overlap distance (default: `0.5`)
- `--name2-mode`: `coverage` (default), `full`, or `none`
- `--invalid-read-policy`: `skip` (default) excludes only the invalid read track; `fail` restores strict parsing

Behavior summary:
- Pass 1 uses the `ratio` distance (`1 - overlap / union_len`) with `--cutoff1`.
- Pass 2 uses the `ratio_short` distance (`1 - overlap / min_len`) with `--cutoff2`.
- Exon and intron distances are combined with `--intron-weight`.

Example:
```bash
trackcluster cluster \
  --reads reads.bed \
  --reference ref.bed \
  --out isoform.bed \
  --threads 8
```

### `trackcluster count`
Compute isoform counts from an existing `flow`/`clusterj-batch` output directory.
This is the recommended mode because unique assignment is rerun inside each
`<gene-path-key>/` folder using `<gene-path-key>_nano.bed`, the key-named
per-gene isoform BED, and `<gene-path-key>_read_to_isoform.tsv` before merged
counts are written.
When using manual `preparedir` -> `clusterj_batch`, use the same directory for
`clusterj_batch --input-root` and `--output-root` if you want to recount with
`trackcluster count --output-root`; a separate cluster output directory does not
contain `<gene-path-key>_nano.bed` unless you copy the prepared inputs there.

Input:
- `--output-root/-o`: existing output directory containing per-gene folders
- `--prefix`: prefix for merged outputs
- `--reference/-r`: reference BED
- `--cluster-mode`: `clusterj` (default) or `cluster`
- `--assignment-mode`: `unique` (default) or `fractional`
- `--unique-assignment-junction-offset` (default: `15`): unique-assignment intron tolerance in bp

Output:
- `<prefix>_isoform.bed`
- `<prefix>_read_to_isoform.tsv`
- `<prefix>_read_to_isoform.unique.tsv` in unique mode
- `<prefix>_unique_assignment.provenance.tsv` in unique mode
- `<prefix>_isoform_count.csv`
- description files under `<prefix>_*`

Example:
```bash
trackcluster count \
  --reference ref.bed \
  --output-root out \
  --prefix sample
```

Legacy low-level mode can still count a standalone isoform BED. In this mode,
the selector sees the supplied BED/mapping as its full scope, so prefer
`--output-root` when continuing a per-gene cluster run. Unique mode also writes
`<out-stem>.provenance.tsv` next to `--out`; for example,
`--out isoform_count.csv` writes `isoform_count.provenance.tsv`. Fractional mode
does not emit that file and removes a stale unique-mode provenance file on a
successful rerun with the same `--out`.

```bash
trackcluster count \
  --reads reads.bed \
  --reference ref.bed \
  --isoform isoform.bed \
  --read-to-isoform isoform.read_to_isoform.tsv \
  --out isoform_count.csv
```

Single-sample and aggregate count files are standards-compliant CSV with
columns `gene,isoform_id,count`. Novel catalog records use deterministic
`tc_novel_v1:` structural IDs. A repeated read label is one abundance molecule;
duplicate mapping rows are idempotent, while unique assignment rejects
conflicting structures bearing the same label.

### `trackcluster count-multi`
Compute per-sample isoform counts/proportions from pooled isoforms using a sample manifest.

Input:
- `--manifest`: TSV with required columns `sample`, `reads`; optional `group`
- `--reference/-r`: reference BED
- `--isoform/-i`: pooled isoform BED (typically from `flow --manifest` or pooled `clusterj`)
- `--read-to-isoform`: optional mapping TSV (recommended; required when isoform `name2` does not embed read IDs; auto-discovered when next to the isoform BED)
- `--assignment-mode`: `unique` (default; expand candidates against the isoform catalog and assign each read to the closest compatible isoform using read/isoform structure) or `fractional` (split multi-mapped reads across mapped candidates)
- `--unique-assignment-junction-offset` (default: `15`): unique-assignment intron tolerance in bp
- `--out/-o`: output prefix

Outputs (`--out <prefix>`):
- `<prefix>.isoform_count.csv` (aggregate counts derived from the per-sample matrix)
- `<prefix>.isoform_usage.long.tsv`
- `<prefix>.isoform_counts.matrix.tsv`
- `<prefix>.isoform_usage.group.tsv` (when at least one sample has a non-empty `group`)
- `<prefix>.unique_assignment.provenance.tsv` (unique mode)

Aggregate count semantics:
- `count` is exactly the sum of the sample columns in `<prefix>.isoform_counts.matrix.tsv` for the same isoform.
- Columns are `gene,isoform_id,count` and identifiers are CSV escaped.

Long-table semantics:
- One row per `(gene, isoform, sample)` with non-zero count.
- `proportion` is within-gene usage for that `(gene, sample)`:
  `proportion = count / sum(count of all isoforms for the gene+sample)`.

Matrix-table semantics:
- Rows are `(gene, isoform)`.
- Columns are samples in manifest order.
- Missing isoforms in a sample appear as `0`.

Example:
```bash
trackcluster count-multi \
  --manifest samples.tsv \
  --reference ref.bed \
  --isoform pooled_isoform.bed \
  --read-to-isoform pooled_read_to_isoform.tsv \
  --out out/pooled
```

### `trackcluster desc`
Describe/classify isoforms vs reference, and detect fusions.

Outputs (`--out <prefix>`):
- `<prefix>_desc.txt`
- `<prefix>_class4.txt`
- `<prefix>_class12.txt`
- `<prefix>_fusion.txt`

Important:
- `desc` groups isoforms by gene name field. If your isoforms do not have gene names, run `addgene` (or use `flow` / the `preparedir` -> `clusterj_batch` workflow).
- `_class12.txt` follows the paper's Figure 2 contract: 11 novel-event labels plus `reference`, with a 5% summed-exon threshold for Extra/Missing UTR. The exact labels and the “later match overwrites earlier match” priority are documented in [`behavior/desc.md`](behavior/desc.md).
- `--end-shift-bp` controls optional strand-aware text in `_desc.txt` only; it does not add or change Figure 2 `class12` categories.

Example:
```bash
trackcluster desc \
  --isoform isoform.bed \
  --reference ref.bed \
  --out sample
```

## `clusterj_batch`
Run `clusterj` per gene folder in parallel (optionally run `preparedir` first).

Scope:
- This binary is junction-mode only.
- For overlap-mode batched clustering, use `trackcluster flow --cluster-mode cluster` or `trackcluster cluster` directly on a single reads/reference pair.

Useful flags:
- `--gene-list <PATH>`: file containing one biological gene ID per line. It is
  required unless `--prepare-reads`, `--prepare-reference`, and
  `--prepare-prefix` are supplied in the same invocation. Blank lines and
  comments are ignored; entries are not directory names or encoded path keys.
  Batch mode never discovers arbitrary directories, so stale folders from an
  earlier preparation cannot silently enter a run.
- `--sw-score`: Smith-Waterman cutoff for SL-supported 5' protection (default: `-1`, no SW/SL 5' signal; pass a non-negative cutoff such as `11` only when BED score is valid SL/SW 5' evidence)
- `--batch-size`, `--batch-rounds`: bounds for very large genes
- `--name2-mode`: `coverage` (default), `full`, or `none` (controls isoform `name2` payload size)
- `--platform-preset`: `generic` (default), `rna002`, or `rna004`. Presets seed junction correction, SL 5', and same-junction 3' defaults; explicit option values override the preset.
- `--junction-correction-offset` (default: `10`; `rna002`: `15`; `rna004`: `10`), `--junction-correction-min-support` (default: `5`): internal junction-site correction controls. This offset is separate from the SL/5' and 3' terminal merge/protection offsets.
- `--sl-partial-5prime-offset` (default: `15`), `--sl-same-junction-5prime-offset` (default: `25`), `--sl-5prime-cluster-offset` (default: `15`), `--sl-5prime-min-support` (default: `2`): SL 5' merge controls
- `--same-junction-3prime-offset` (default: `50`), `--3prime-cluster-offset` (default: active junction correction offset), `--3prime-min-support` (default: `5`): same-junction 3' terminal retention controls
- `--heartbeat-seconds`, `--heartbeat-top`: periodic status line (and which gene(s) are currently in-flight when progress is not moving)
- `--invalid-read-policy`: `skip` (default) excludes only malformed/empty-ID read tracks within each gene and writes `<gene-path-key>/rejected_reads.tsv`; `fail` restores strict per-gene read parsing.
- `--strict-gene-errors`: return nonzero after all genes finish if any gene failed. Without it, failures are reported and successful/verified genes remain available to downstream callers.
- `--max-reads-per-gene` (default: `50000`; set `0` to disable),
  `--downsample-gene <BIOLOGICAL_GENE_ID>` (repeatable),
  `--downsample-seed`: per-gene downsampling (writes
  `clusterj_batch_downsample.tsv`). A later `flow --count-only` or
  `count --output-root` refuses independently downsampled molecules that occur
  in more than one selected gene; disable downsampling for every affected gene.

Typical usage (after `preparedir`):
```bash
clusterj_batch \
  --input-root tracktest \
  --gene-list tracktest/sample_gene.txt \
  --output-root tracktest \
  --threads 8 \
  --force
```

`clusterj_batch` writes the active platform preset, junction correction settings, SL settings, and same-junction 3' settings into `clusterj_batch_summary.txt` for reproducibility.

Gene-local errors do not stop the default run. The summary uses `status\tpartial`, the detailed
diagnostics remain in `clusterj_batch_errors.txt`, and only genes with complete verified artifacts
are considered mergeable. `--strict-gene-errors` keeps the previous all-or-nothing exit behavior.
If every gene fails, or worker/output infrastructure fails, the command still exits nonzero because
no safe result exists. A downstream `trackcluster flow --count-only` rerun
remains strict and does not accept stale or missing manifests; `clusterj_batch`
itself has no `--count-only` option.

The invalid-read policy is narrower than the gene-error policy. In `skip` mode, only a BED parse
failure attributable to one read record or an empty read ID is downgraded to a rejected-read row.
Read-file I/O failures, malformed references, invalid options, and clustering/counting errors still
fail their enclosing stage; they are never relabeled as rejected reads. In `fail` mode, an invalid
read makes that gene fail, and `--strict-gene-errors` additionally makes any such gene failure stop
batch-level downstream publication.

One-command convenience (prepare + batch):
```bash
clusterj_batch \
  --prepare-reads reads.bed \
  --prepare-reference ref.bed \
  --prepare-prefix sample \
  --input-root tracktest \
  --output-root tracktest \
  --threads 8 \
  --force
```

Use a distinct `--output-root` only when you intend to consume the per-gene
cluster outputs directly or concatenate them manually. For
`trackcluster count --output-root` unique assignment, the count root must also
contain the prepared `<gene-path-key>_nano.bed` files.
