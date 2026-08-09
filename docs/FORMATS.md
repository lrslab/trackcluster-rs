# File formats

See [`INTERCHANGE.md`](INTERCHANGE.md) for BAM/GFF3/GTF conversion and transcript exports.

## 0.2.0 output migration

Version 0.2.0 deliberately changes several
previously ambiguous output contracts. Count CSV is now RFC 4180 with the
header `gene,isoform_id,count`; novel transcript IDs use the `tc_novel_v1:`
structural namespace; and full `name2` payloads use the percent-encoded
`tc_name2_v1:` codec. Description files move from the legacy headerless layout
to `trackcluster-description-v2`. Consumers that assumed two-column count CSV,
representative-read novel IDs, raw comma-split `name2`, or headerless
description files must migrate using the contracts below. Readers retain
support for legacy unescaped `name2` payloads.

Description files use schema `trackcluster-description-v2`. Each starts with a
`#schema` line followed by an explicit column header, replacing the legacy
headerless files. The Figure 2 category values, UTR threshold, and class overlap
priority are defined in [`behavior/desc.md`](behavior/desc.md).

- `_desc.txt`: `isoform_id`, `reference_id`, `gene_id`, `missing_features`,
  `extra_features`
- `_class4.txt`: `isoform_id`, `class`
- `_fusion.txt`: `isoform_id`, `gene_ids`
- `_class12.txt`: `isoform_id`, `class`

The Rust rewrite primarily works with BED12 and bigGenePred-like files.

## BED coordinate conventions
- Coordinates are **0-based, half-open**: `[start, end)`.
- Two intervals overlap iff `a.start < b.end && b.start < a.end`.

## Required fields
At minimum, inputs should include BED12 columns:
- `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`
- `thickStart`, `thickEnd`, `itemRgb`
- `blockCount`, `blockSizes`, `blockStarts`

## bigGenePred extras
Many TrackCluster datasets use bigGenePred extra columns (e.g., gene name / annotation metadata).
The Rust parser preserves extra trailing fields and writes them back out when present.

### Converter-produced BED12+8

`bam2bigg` and `gff2bigg` write tab-separated, bigGenePred-compatible
**BED12+8 text**. They do not create the indexed binary `.bb` representation.
Both use BED's zero-based, half-open coordinates and emit all eight TrackCluster
extension columns in this order:

| Extra index | bigGenePred name | `bam2bigg` | `gff2bigg` |
| ---: | --- | --- | --- |
| 0 | `name2` | `none` | `none` |
| 1 | `cdsStartStat` | `none` | `none` |
| 2 | `cdsEndStat` | `none` | `none` |
| 3 | `exonFrames` | one `-1` per block | one `-1` per block |
| 4 | `type` | `nanopore_read` | `isoform_anno` |
| 5 | `geneName` / TrackCluster gene ID | `none` | resolved annotation gene(s), otherwise `none` |
| 6 | `geneName2` / TrackCluster sample group | `--group` or BAM filename stem | `none` |
| 7 | `geneType` / TrackCluster reserved | `none` | `none` |

The standard BED fields differ as follows:

- `bam2bigg` uses the BAM query name, MAPQ as BED score, flag-derived strand,
  strand color (`250,128,114` for `+`, `64,224,208` for `-`), and CIGAR-derived
  blocks. CIGAR `N` alone splits blocks. `thickStart` and `thickEnd` are both
  zero.
- `gff2bigg` uses the annotation transcript identity, score `100`, resolved
  exon strand, `itemRgb=0`, and exon-derived blocks. GFF/GTF coordinates
  `[start,end]` become BED `[start-1,end)`. CDS/UTR/phase are not transferred,
  so `thickStart` and `thickEnd` are both zero.

`gff2bigg` sorts output deterministically by transcript structure. `bam2bigg`
preserves retained BAM record order and preserves separate alignment instances,
including repeated query names.

### Gene name field (TrackCluster convention)
This Rust rewrite follows the legacy TrackCluster convention of storing a gene name in an extra field:
- Gene name lives at **extra field index 5** (0-based within `extra_fields`).
- Unassigned is `none`.
- Multi-gene values are joined with `||` (example: `GENE1||GENE2`).

### Biological gene IDs and filesystem keys

Per-gene paths never interpolate an unchecked biological gene ID. Gene IDs must be non-empty,
single path components and cannot be absolute paths, `.`/`..`, contain `/` or `\\`, or contain
NUL/control characters. The maximum accepted length is **4096 UTF-8 bytes**; longer IDs fail
validation before any per-gene path is constructed.

TrackCluster derives a deterministic filesystem key for each valid ID. ASCII letters, digits, `_`,
`-`, and `.` keep their historical spelling (for example, `GENE-1`); Unicode and other punctuation
are percent-encoded when the result fits, and oversized path keys use a fixed-size stable hash.
Per-gene directory and artifact filenames use this key. `<prefix>_gene_paths.tsv` and the batch-level
`clusterj_batch_gene_paths.tsv`/`cluster_batch_gene_paths.tsv` are versioned mappings from the
biological ID to its key. Each gene directory also contains `.trackcluster_gene_id`, so metadata
validation can round-trip hashed IDs without treating arbitrary directories as selected genes.
Before any preparation or batch publication, TrackCluster rejects a key that
equals one of its prefix-scoped merged/preparation filenames or fixed batch
report filenames. This keeps gene directories structurally disjoint from every
top-level pipeline artifact; changing the output prefix resolves a
prefix-derived collision, while a collision with a fixed
`clusterj_batch_*`/`cluster_batch_*` name requires a different biological gene
ID. These run-scoped namespace restrictions are checked in addition to the
context-free `GeneId` syntax and path-key rules above.

Gene lists and selectors stay in the biological namespace:
`<prefix>_gene.txt`, `clusterj_batch --gene-list`, and `--downsample-gene`
contain or accept biological gene IDs, not encoded keys or directory names.
The prefix-scoped gene list is also the preparation commit marker: preparation
publishes it last. An empty list means a replacement failed after publication
began, so the prepared generation is incomplete and must not be clustered or
counted until preparation succeeds again.

### Subread list / `name2` field (TrackCluster convention)
Isoforms produced by `clusterj`/`cluster` use extra field index `0` for the TrackCluster-style `name2` payload.

By default (`--name2-mode coverage`), this stores only a coverage value:
- `extra_fields[0] = "|<coverage>"` (no read IDs)
- `--name2-mode none`: `extra_fields[0] = "none"` (no payload)
- `--name2-mode full`: `extra_fields[0] = "tc_name2_v1:<encoded-read1>,<encoded-read2>,...,|<coverage>"`
  - Example: `tc_name2_v1:readA,readB,readC,|2.5`
  - Read IDs are UTF-8 percent encoded. Commas, pipes, percent signs, whitespace,
    and other reserved bytes therefore round-trip without being mistaken for
    payload delimiters.

Readers still accept the pre-v1 unescaped comma-separated payload for backward
compatibility. That legacy form cannot represent a comma inside a read ID;
re-emit the catalog or use the mapping TSV to migrate such data. All newly
written full payloads use `tc_name2_v1`.

When read IDs are omitted from `name2`, use the `*_read_to_isoform.tsv` mapping written by `clusterj`/`cluster`/`flow`. For ordinary recounts, prefer `trackcluster count --output-root <out> --prefix <prefix>`; it reads each gene folder directly, so unique assignment and retained-intron checks stay gene-local. The legacy standalone BED mode can still take `--read-to-isoform` (or auto-discover it next to the isoform BED), but its unique assignment scope is the supplied merged input. Counting defaults to unique best assignment: it expands compatible candidates against the isoform catalog, then selects the closest isoform per read before counting. Pass `--assignment-mode fractional` for compatibility with split multi-mapped counts from the mapping file.

Mapping files contain exactly two raw TSV fields, `read_id` and `isoform_id`,
without a header. Leading and trailing spaces are significant identity bytes and
round-trip unchanged; empty fields, tabs within a field, and embedded line
breaks are invalid.

In `flow` unique assignment mode, `<prefix>_read_to_isoform.tsv` remains the raw merged mapping from per-gene clustering. The selected mapping actually used for final counts is written separately as `<prefix>_read_to_isoform.unique.tsv`; use that file when auditing or reproducing unique-mode counts. `<prefix>_unique_assignment.provenance.tsv` records the effective junction tolerance, ordered one-to-one matcher, and explicit no-collapse policy for microfeatures.

## Stable identity contract

- Reference isoform IDs are preserved. Empty IDs, duplicate reference IDs, and
  references that claim the reserved `tc_novel_v1:` namespace are rejected.
- Novel isoforms use `tc_novel_v1:<gene-hex>:<chromosome-hex>:<strand>:<exons>`.
  This is a lossless structural namespace, not a representative read name or a
  truncated hash, so distinct `(gene, chromosome, strand, exon-chain)` tuples
  cannot silently collide.
- Catalog IDs are validated globally before a merged catalog is atomically
  published and again at counting boundaries.
- A read ID is an abundance **molecule ID**. Multiple BED rows with that ID are
  alignment instances during clustering, but repeated identical mapping rows
  are idempotent during counting. Fractional counting splits one molecule over
  its distinct candidate isoforms. Unique assignment accepts identical duplicate
  structures but rejects conflicting alignments for one molecule ID because
  there is no unambiguous structure to score.
- A read ID equal to a reference transcript ID remains a read; source provenance
  is structural and is never inferred from matching strings.

### Rejected-read diagnostics

`flow`, `preparedir`, `clusterj_batch`, `clusterj`, and `cluster` use `--invalid-read-policy skip` by default. This recovery
is deliberately record-local: a malformed read BED row or a parsed read row whose ID is empty is
excluded, while the remaining valid read tracks continue. Preparation publishes
`<prefix>_rejected_reads.tsv`; per-gene clustering publishes `rejected_reads.tsv` in each gene
directory for read tracks rejected while that gene is loaded. Header-only files mean that no read
was rejected at that boundary.

The direct single-gene commands derive the same diagnostic from `--out`: for example,
`--out isoform.bed` writes `isoform.rejected_reads.tsv`.

Both files have this stable header:

```tsv
source_path	line	read_id	kind	reason
```

- `source_path` and the 1-based physical `line` locate the input record.
- `read_id` is populated when an identifier is available; it is empty when parsing failed before a
  usable ID could be recovered.
- `kind` is `parse` for malformed BED and `identity` for an empty read ID.
- `reason` contains the parser or identity diagnostic. Backslashes, tabs, carriage returns, and
  newlines inside text fields are escaped as `\\`, `\t`, `\r`, and `\n`.

`--invalid-read-policy fail` restores strict read-track parsing: the first malformed or empty-ID
read fails the enclosing preparation/gene stage. Read-file I/O errors, reference parse/identity
errors, invalid configuration, and algorithm/integrity failures are fatal at their existing
boundary under both policies and never appear in these TSVs.

The rejected-read TSV is an additive audit artifact. Existing input BED, isoform BED, mapping,
count, and description/classification schemas are unchanged, as are the clustering and isoform
classification rules. However, rejected reads provide no junction, terminal, coverage, mapping, or
counting evidence, so filtering them can change which isoforms are called and consequently the
contents (but not the schema or category definitions) of classification outputs.

### Per-gene `run.json`

`flow` and `clusterj_batch` publish schema-versioned JSON completion manifests in each gene
directory. Source-aware tool identity uses manifest schema version 3; older per-gene manifests are
rebuilt rather than reused. The manifest has `status: "complete"`, a request fingerprint, input SHA-256 hashes,
effective options (including invalid-read policy, assignment mode, and unique-assignment junction
tolerance), package version, Git commit, deterministic source fingerprint, the per-gene seed, and
output SHA-256 hashes, byte sizes, and record counts. Clean Git checkouts use
`clean`. Cargo package builds and dirty source checkouts use a SHA-256
fingerprint of the actual build-source snapshot; rebuilding after editing an
unpacked package therefore cannot retain the official package's cache
identity. The
manifest is a completion marker rather than a user-editable configuration file: changing it or any
recorded output makes that gene stale on the next normal run. `flow --count-only` and
`count --output-root` require this marker for every selected gene and revalidate its self-fingerprint,
gene/mode/tool identity, current prepared-input contents, and every output hash, size, and record
count before publishing merged results. Legacy folders without a valid manifest must first be
rebuilt by a normal flow run.

Artifact publication is atomic per file, using a synced temporary sibling and rename. The
per-gene manifest is invalidated before replacement and published last, so it is the only reusable
completion marker. A set of flat merged/count/description files is not a cross-file filesystem
transaction; after an interruption during their publish phase, rerun the final stage. Those derived
files are outputs only and are never accepted as evidence that per-gene work completed.

### Batch summary status and recoverable gene errors

Batch summaries use `status=complete` when every gene succeeds, `status=partial` when the default
continue policy excludes one or more failed genes but retains at least one verified result, and
`status=failed` for strict mode, infrastructure failure, or an all-gene failure. The fields
`gene_error_policy`, `mergeable_genes`, `excluded_failed_genes`, and `infrastructure_errors` make
that distinction machine-readable. Detailed gene diagnostics are written to
`*_errors.txt` only when errors are recorded; a clean run removes any stale
error report.

Only processed genes, hash-verified cache hits, and semantic-empty genes with a valid empty
completion manifest are mergeable. Failed genes, their old artifacts, and their old downsampling
scale factors are excluded. This recovery policy applies to normal per-gene execution; count-only
manifest/hash verification remains all-or-nothing.

## Multi-sample manifest TSV
`count-multi` and `flow --manifest` expect a tab-separated manifest with:

- Required columns:
  - `sample`: unique sample name
  - `reads`: BED path (absolute or relative to manifest file)
- Optional columns:
  - `group`: condition/group label for pseudo-bulk summaries

Example:
```tsv
sample	group	reads
S1	control	S1.reads.bed
S2	treated	S2.reads.bed
```

Constraints:
- `sample` must be unique.
- `sample` cannot contain `::` (reserved delimiter for pooled read IDs).
- Missing files fail fast with an error.

## Pooled read IDs
When pooling manifest reads, trackcluster rewrites read IDs as:
`<sample>::<orig_read_id>`

This guarantees sample identity for downstream per-sample counting.

## `count-multi` outputs

### Aggregate count table (`*.isoform_count.csv`)
Columns:
- `gene`
- `isoform_id`
- `count`

Semantics:
- The file is RFC 4180-compatible CSV written with field escaping; commas and
  quotes in identifiers do not create extra columns.
- `count` is exactly the sum of the sample columns in the matching
  `*.isoform_counts.matrix.tsv`.
- In `flow --manifest`, the main `<prefix>_isoform_count.csv` is synchronized
  from the same aggregate table, so total counts and sample counts share one
  assignment result.

### Long table (`*.isoform_usage.long.tsv`)
Columns:
- `gene`
- `isoform_id`
- `sample`
- `group` (present when any sample has a group)
- `count`
- `proportion`
- `gene_total`

Semantics:
- One row per `(gene, isoform_id, sample)` with non-zero count.
- `proportion` is within-gene usage:
  `count / sum(counts for all isoforms in the same gene and sample)`.

### Matrix table (`*.isoform_counts.matrix.tsv`)
Columns:
- `gene`
- `isoform_id`
- one column per sample (manifest order)

Semantics:
- Missing isoforms in a sample are represented as `0`.

### Group table (`*.isoform_usage.group.tsv`)
Only emitted when at least one sample has a non-empty `group` value.

Columns:
- `gene`
- `isoform_id`
- `group`
- `count`
- `proportion`
- `gene_total`

### Unique-assignment provenance (`*.unique_assignment.provenance.tsv`)

Emitted only in unique assignment mode. The table records the effective
`unique_assignment_junction_offset` and the one-molecule/one-isoform matching
policy used to generate the count tables. Fractional mode removes a stale
unique-mode provenance file when the same output prefix is rerun successfully.

## Isoform-level modification formats

The modification boundary is a normalized genomic read/site observation table.
`mod-import-m6anet` and `mod-import-dorado` produce that table;
`mod-aggregate` consumes it and produces three TSVs; `mod-contrast` consumes the
design TSV and an explicit contrast specification.

All TSV schemas below have an exact header and column order. Fields are separated
by one tab. Where a field is nullable, the missing token is the literal `NA`, not
an empty field. A numeric `0` is a measured count, probability, fraction, or
effect and is never a missing-value token. JSON nullable values use JSON `null`,
not the TSV token `NA`.

### Caller import outputs

Both native importers atomically write these three files for output prefix
`<out>`:

- `<out>.observations.tsv`, using the normalized schema below;
- `<out>.assay.json`, using the assay metadata schema below; and
- `<out>.import_qc.tsv`, with exact header `metric<TAB>value` and one
  deterministic caller-specific metric per row.

The m6Anet RNA002 importer requires `data.indiv_proba.csv` and an explicit
two-column `read_index<TAB>read_id` map. `read_index` is compared as an opaque
token so upstream multi-input suffixes such as `_0` are preserved. It projects zero-based
`transcript_position` through an exact-version-matching GTF/GFF transcript.
The importer first derives the required source transcript set and retains only
those annotation records (and same-stable-ID version alternatives), avoiding a
whole-annotation in-memory transcript catalog. QC reports
`projection_transcripts_loaded`.
Optional `data.info` verifies retained sites and read counts; optional
`data.site_proba.csv` cross-checks the same retained sites, `n_reads`, and
`mod_ratio`, but never creates read-level observations. `candidate_rule` is
normalized to an odd-length IUPAC motif of 3–31 bases with a literal `A` at its
center; when `data.site_proba.csv` is supplied, every concrete k-mer must match
that motif. Its QC includes input,
deduplication, read-map, transcript/site, data-info, and site-probability
counters. Known source read-threshold presets are `0.033379376` for
`HCT116_RNA002` and `0.0032978046219796` for `arabidopsis_RNA002`. These are
source `mod_ratio` cross-checks, not site-probability cutoffs and not a hidden
replacement for the aggregation analysis threshold.

The Dorado importer accepts one target modification code from a primary
genome-aligned BAM. `candidate_rule` may be
`all-target-canonical-bases` or an odd-length centered IUPAC DNA motif of 3–31
bases. Motifs are case-normalized (`U` becomes `T`), must have the target
canonical base at the center, and are matched on the as-sequenced read
orientation. Only matching bases belong to the implicit candidate universe;
explicit target calls outside it are invalid. It requires one each of `MM:Z`,
`ML:B:C`, and integer `MN`,
consumes ML in group/position/code order before target filtering, and validates
all CIGAR operations. The SAM canonical base `N` in an MM group matches every
query base rather than a literal `N` sequence character. An ML byte `N`
represents `[N/256,(N+1)/256)` and the
normalized point estimate is `(N+0.5)/256`. Explicit calls on insertions or soft
clips have no genomic observation and are counted in QC. For projectable
canonical bases omitted by `.` or an omitted MM skip marker, it emits
`implicit_below_emission_threshold`; for `?`, it emits `unknown` by default.
The explicit `question_mark_policy=below_emission_threshold` override requires
a positive source emission threshold and records the raw `?` group count while
emitting omitted candidates as below-threshold. Import QC
separately reports record filters, invalid reasons, MM group semantics,
raw target canonical bases, rule-matching candidate bases, explicit/implicit
candidates, unprojectable calls, and emitted observations.
`candidate_observations_complete` describes the retained, projectable genomic
candidate universe; skipped invalid records or a missing target group on a read
that contains at least one rule-matching target base make it false.

### Normalized observation TSV

The required header is:

```tsv
assay_id	sample	read_id	chrom	pos0	strand	mod_code	probability	observation_state	context	source_transcript_id	source_pos0
```

| Column | Type | Contract |
| --- | --- | --- |
| `assay_id` | string | Assay compatibility stratum. Required, not `NA`, and free of control characters. |
| `sample` | string | Biological sample ID. Required and subject to the same identifier rules. |
| `read_id` | string | Final TrackCluster molecule ID in `<sample>::<source_read_id>` form. Both parts must be non-empty and the prefix must equal `sample`. |
| `chrom` | string | Reference sequence name. Required, not `NA`, and free of control characters. |
| `pos0` | unsigned 32-bit integer | Zero-based coordinate of the genomic base, not a half-open interval start for a multi-base feature. |
| `strand` | enum | Exactly `+` or `-`; `.` is rejected. |
| `mod_code` | string | SAM-style canonical-base/strand/modification identity, for example `A+a` or `C+76792`. The canonical base is one of uppercase `A`, `C`, `G`, `T`, `U`, or `N`; the sign is `+` or `-`; the suffix is one or more ASCII alphanumeric characters. |
| `probability` | finite float in `[0,1]` or `NA` | Required for `explicit_probability`; required to be `NA` for the other observation states. `0` and `1` are valid probabilities. |
| `observation_state` | enum | `explicit_probability`, `implicit_below_emission_threshold`, or `unknown`. |
| `context` | string or `NA` | Caller sequence context or candidate-rule label. A present value follows the identifier rules; an empty field is invalid. |
| `source_transcript_id` | string or `NA` | Optional source transcript provenance. A present value follows the identifier rules, and it must be present exactly when `source_pos0` is present. |
| `source_pos0` | unsigned 64-bit integer or `NA` | Zero-based source-transcript offset retained only as provenance. It must be present exactly when `source_transcript_id` is present. |

The observation uniqueness key is
`(assay_id, sample, read_id, chrom, pos0, strand, mod_code)`. Input is sorted by
that key after parsing. Repeated rows with the same key and identical values in
every field are folded and counted as `duplicate_exact`. Reusing a key with a
different probability, state, context, or source-coordinate payload is a fatal
conflicting duplicate; no aggregate output is produced.

`observation_state` records caller semantics, not an aggregate hard call:

| `observation_state` | `probability` | Meaning before aggregation |
| --- | --- | --- |
| `explicit_probability` | numeric, including `0` | The caller emitted a read-level probability. The analysis threshold later determines modified versus unmodified. |
| `implicit_below_emission_threshold` | `NA` | The caller omitted or encoded the candidate implicitly as below its source emission threshold. It is callable only under the metadata and threshold rule below. |
| `unknown` | `NA` | The candidate is represented, but its modification state is not interpretable. It never enters `n_callable`. |

### Assay metadata JSON

Each normalized observation file is paired with one schema-versioned provenance
object. The complete JSON shape is:

```json
{
  "schema_version": 1,
  "assay_id": "rna004_m6a_model_1",
  "caller": "dorado",
  "caller_version": "unknown",
  "model_id": "rna004_model_name",
  "chemistry": "RNA004",
  "candidate_rule": "all-target-canonical-bases",
  "source_emission_threshold": 0.1,
  "source_site_filter": "none",
  "candidate_observations_complete": true,
  "implicit_skip_policy": "low_probability",
  "coordinate_source": "genome_aligned_input",
  "read_id_mapping": "sample_prefixed_source_read_id",
  "source_files": ["normalized-source-identifier"]
}
```

Unknown JSON fields are rejected. On input, `source_emission_threshold` may be
omitted and is then treated as `null`; `source_files` may be omitted and defaults
to an empty array. All other fields are required. The writer emits both optional
fields explicitly and terminates pretty-printed JSON with a newline.

| Field | Type | Contract |
| --- | --- | --- |
| `schema_version` | integer | Must be `1`. |
| `assay_id` | string | Must equal the manifest row and every observation row for that input. It defines a compatibility stratum, not merely a display label. |
| `caller` | string | Caller identity. |
| `caller_version` | string | Caller version; use a meaningful non-`NA` value such as `unknown` if it cannot be recovered. |
| `model_id` | string | Exact model identifier. |
| `chemistry` | string | Sequencing chemistry. |
| `candidate_rule` | string | Candidate-universe rule or context label. |
| `source_emission_threshold` | finite float in `[0,1]` or `null` | Threshold below which the source omitted candidates. It is mandatory when `implicit_skip_policy` is `low_probability`. |
| `source_site_filter` | string | Site/candidate filtering provenance. |
| `candidate_observations_complete` | boolean | Whether every observation in the retained candidate universe is represented in the normalized TSV. This is an assertion; aggregation does not synthesize rows for absent candidates. |
| `implicit_skip_policy` | enum | `low_probability`, `unknown`, or `not_applicable`. |
| `coordinate_source` | string | How source coordinates became genomic coordinates. |
| `read_id_mapping` | string | How source IDs became final TrackCluster read IDs. |
| `source_files` | array of strings | Source paths or immutable source identifiers; each entry must be non-empty, not `NA`, and free of control characters. |

All scalar string fields are required to be non-whitespace, not the literal
`NA`, and free of control characters. Inputs sharing an `assay_id` must match
exactly in caller, caller version, model, chemistry, candidate rule, source
emission threshold, site filter, coordinate source, and read-ID mapping.
`candidate_observations_complete`, `implicit_skip_policy`, and `source_files`
may differ between samples: completeness and the actually observed skip state
are sample-level QC. Each sample's skip policy is still applied independently
when deciding whether implicit observations are callable. An incomplete sample
is retained for audit but made ineligible. Dorado's configured question-mark
policy is embedded in `source_site_filter`, so different overrides remain
incompatible even if no omitted candidate happened to occur in one sample.
Thus caller, model, chemistry, or configured probability semantics cannot be
pooled merely by assigning them the same label.

### Modification manifest TSV

`mod-aggregate --mod-manifest` accepts this exact five-column header:

```tsv
sample	assay_id	observations	assay_metadata	coverage_bam
```

| Column | Contract |
| --- | --- |
| `sample` | Must match a sample in the ordinary TrackCluster sample manifest. It cannot be empty, `NA`, contain a control character, or contain the reserved delimiter `::`. |
| `assay_id` | Must match the paired metadata and observation rows. It cannot be empty, `NA`, or contain a control character. |
| `observations` | Required normalized observation TSV path. |
| `assay_metadata` | Required assay metadata JSON path. |
| `coverage_bam` | Optional genome-aligned BAM for exact assigned-molecule coverage, or `NA`. Raw query names are tagged with the manifest sample; an existing matching `<sample>::` prefix is preserved and a conflicting prefix is rejected. Every uniquely assigned read for the sample must be present, and duplicate primary names are rejected. |

The first non-empty, non-comment line must be the exact header. A comment line
starts with `#` in column one. Every data row has exactly five fields, and at
least one data row is required. Relative paths are resolved relative to the
manifest. Input paths must identify regular files; symbolic links and missing
files are rejected. Each `(sample, assay_id)` pair must be unique.

### Modification pseudo-sample bundle

`mod-subsample` writes a self-contained directory with these stable top-level
artifacts:

- `samples.tsv`: generated samples with empty groups and per-sample reads BEDs.
- `mod_samples.tsv`: generated observation, metadata, and optional BAM inputs.
- `read_to_isoform.unique.tsv`: selected source assignments with rewritten
  pseudo-sample prefixes. Selected unassigned reads have no fabricated row.
- `subsample_read_ids.tsv`: exact selected molecule membership and source
  isoform assignment, or `NA` when the selected source read was unassigned.
- `subsample_qc.tsv`: source/selected/assigned/observation/BAM counts per
  pseudo-sample and assay.
- `sample_provenance.tsv`: `technical_pseudo` status and parent sample/group.
- `overlap_qc.tsv`: pairwise read intersection, union, and Jaccard values.
- `subsample_provenance.json`: algorithm version, mode, seed, depth, source
  paths, and an explicit `biological_replicates=false` declaration.
- `SHA256SUMS`: deterministic checksums for every other published bundle file.
- `samples/*.reads.bed`, `observations/*.observations.tsv`,
  `assays/*.json`, and, when supplied by the source manifest,
  `coverage/*.bam`.

The source reads BED defines the sampling universe. The same read membership is
then applied to all assays without including `assay_id` in the sampling hash.
All read-site rows and modification codes for a molecule therefore move
together. `disjoint` output sets have zero pairwise overlap; `independent`
sets contain no within-sample duplicate but can overlap each other. Coverage
BAMs are complete sequential BAMs for the selected query names and do not
require an index for `mod-aggregate`.

Pseudo-sample groups are empty by contract. Parent group and pseudo status live
in provenance rather than the analysis group column so condition and
interaction contrasts cannot silently treat split molecules as independent
biological samples.

### Assignment, assay, and threshold rules

Modification aggregation requires final unique read-to-isoform assignments.
Read IDs in both the assignment TSV and observation TSV use the same
`<sample>::<source_read_id>` namespace. The assignment sample prefix must exist
in the sample manifest, and every assigned isoform must exist in the catalog.
An identical repeated `(read_id, isoform_id)` mapping is idempotent. Mapping one
read ID to two different isoforms is fatal.

An observation whose read ID has no assignment contributes to join QC as
`unknown_read` but contributes no site universe or modification counts.
Observation sample/assay values that disagree with their modification-manifest
row are fatal. Within an assay and gene, the context for a given
`(chrom, pos0, strand, mod_code)` must be identical, including agreement on
whether it is `NA`; conflicting contexts are fatal.

`--analysis-threshold ASSAY_ID=PROBABILITY` is required exactly once for every
distinct assay in the modification manifest. Missing assays, extra assays, and
duplicate threshold arguments are rejected. Each threshold is finite and in
`[0,1]`. An explicit probability is modified when
`probability >= analysis_threshold` and unmodified when it is lower.
When metadata provides a numeric `source_emission_threshold`, the analysis
threshold must be greater than or equal to it. A lower value is fatal for the
whole aggregation because omitted source candidates would make the requested
hard-call counts incomplete, regardless of the implicit-skip policy.

An `implicit_below_emission_threshold` observation is counted as unmodified only
when all of the following hold:

- the site state is `present` or `context_dependent`;
- `implicit_skip_policy` is `low_probability`;
- `source_emission_threshold` is numeric; and
- `analysis_threshold >= source_emission_threshold`.

Otherwise that observation is counted in `n_unknown`. This prevents a sparse
source encoding from silently turning missing or uninterpretable evidence into
unmodified calls.

The current defaults are `min_callable=1`, `min_read_join_rate=0.9`, and strict
low-join handling. The minimum callable count must be at least one and the join
rate must be finite and in `[0,1]`. If a sample/assay read join rate is below the
minimum, aggregation fails unless `--allow-low-join` is supplied; with that flag,
affected present-site rows are emitted as ineligible.

### `mod-aggregate` outputs

For output prefix `<out>`, the command writes:

- `<out>.mod_join_qc.tsv`
- `<out>.mod_site_join_qc.tsv`
- `<out>.isoform_mod_sites.tsv`
- `<out>.isoform_mod_design.tsv`

The site universe contains joined observations keyed by
`(assay_id, gene, chrom, pos0, strand, mod_code)`. Each universe site is expanded
to every sample participating in that assay and every catalog isoform in that
gene. This expansion makes structural absence and genuine zero-observation rows
explicit. Sites observed only on unassigned reads do not enter the universe.
They are nevertheless retained in `mod_site_join_qc.tsv`, including zero-join
sites.

`site_id` is the string `chrom:pos0:strand`, for example `chr1:1042:+`. It
deliberately excludes `mod_code`; the stable comparison identity is therefore
`(site_id, mod_code)`. The full sites table, design table, contrast spec, and
contrast output all retain `mod_code` as a separate column.

Rows are deterministic. Join QC is ordered by assay then sample. Site-local join
QC is ordered by assay, sample, chromosome, position, strand, and modification
code. Site rows are ordered by assay, sample, gene, isoform, chromosome,
position, strand, and modification code; design rows inherit that order.

#### Join QC (`*.mod_join_qc.tsv`)

The exact header is:

```tsv
assay_id	analysis_threshold	sample	input_rows	valid_rows	projected_rows	joined_rows	joined_reads	read_join_rate	observation_join_rate	unknown_read	unknown_sample	unknown_isoform	duplicate_exact	duplicate_conflict	unprojectable	invalid_probability	candidate_observations_complete
```

| Column | Meaning |
| --- | --- |
| `assay_id` | Assay compatibility stratum. |
| `analysis_threshold` | Hard-call threshold applied to that assay. |
| `sample` | Biological sample ID. |
| `input_rows` | Physical observation data rows, including exact duplicates. |
| `valid_rows` | Unique validated rows after exact deduplication. |
| `projected_rows` | Rows in genomic coordinates. It currently equals `valid_rows` because normalized input is already genomic. |
| `joined_rows` | Unique observation rows whose read ID joined a unique isoform assignment. |
| `joined_reads` | Distinct joined molecule IDs. |
| `read_join_rate` | Distinct joined reads divided by distinct input reads; `0` when there are no input reads. This is the rate used by the low-join gate. |
| `observation_join_rate` | `joined_rows / valid_rows`; `0` when there are no valid rows. |
| `unknown_read` | Unique observation rows with no assignment. |
| `unknown_sample` | `0` in successful current output; sample mismatches fail earlier. |
| `unknown_isoform` | `0` in successful current output; unknown assignment isoforms fail earlier. |
| `duplicate_exact` | Physical rows folded because every field matched an earlier row with the same observation key. |
| `duplicate_conflict` | `0` in successful current output; conflicting duplicates are fatal. |
| `unprojectable` | `0` for the current normalized-genomic input boundary. |
| `invalid_probability` | `0` in successful current output; invalid probabilities are fatal. |
| `candidate_observations_complete` | Metadata assertion copied for audit. |

#### Site-local join QC (`*.mod_site_join_qc.tsv`)

The exact header is:

```tsv
assay_id	sample	site_id	chrom	pos0	strand	mod_code	input_rows	joined_rows	observation_join_rate	passes_min_join_rate
```

Each normalized observation key contains one read and one genomic site, so
`input_rows` and `joined_rows` are distinct read-site counts after exact
deduplication. The rate is evaluated independently for each source genomic site.
Sites with zero joined rows remain present with rate `0`. A present isoform/site
row is ineligible with `site_join_rate_low` when its sample/site rate is below
`--min-read-join-rate`, even if the sample/assay global join rate passes.

#### Complete site table (`*.isoform_mod_sites.tsv`)

The exact header is:

```tsv
assay_id	analysis_threshold	sample	group	gene	isoform_id	site_id	chrom	pos0	strand	mod_code	context	site_state	coverage_basis	n_assigned	n_covering	n_candidate	n_callable	n_modified	n_unmodified	n_unknown	mod_fraction	mean_probability	ci_low	ci_high	eligibility	eligibility_reason
```

| Column | Meaning |
| --- | --- |
| `assay_id` | Assay compatibility stratum. |
| `analysis_threshold` | Applied hard-call threshold. |
| `sample` | Biological sample ID. |
| `group` | Sample-manifest group, or `NA`. |
| `gene` | Single gene ID from isoform metadata; missing or multi-gene (`||`) metadata is rejected. |
| `isoform_id` | Catalog isoform ID. |
| `site_id` | Genomic `chrom:pos0:strand` identity, excluding modification code. |
| `chrom` | Reference sequence name. |
| `pos0` | Zero-based genomic base coordinate. |
| `strand` | `+` or `-`. |
| `mod_code` | SAM-style modification identity retained separately from `site_id`. |
| `context` | Shared caller context/rule label, or `NA`. |
| `site_state` | `present`, `structurally_absent`, `context_dependent`, `reference_base_mismatch`, or `unprojectable`. |
| `coverage_basis` | Coverage provenance token: `bam_exact` when a coverage BAM was supplied, otherwise `unavailable`. `bed_approximate` remains reserved. |
| `n_assigned` | Distinct molecules uniquely assigned to this sample/isoform, independent of whether they have a candidate observation at this site. |
| `n_covering` | Assigned molecules whose primary alignment covers the genomic base with a CIGAR `M`, `=`, or `X` operation. Deletions, reference skips, insertions, and clips do not cover a genomic base. It is numeric, including zero, for `bam_exact` and `NA` for `unavailable`. |
| `n_candidate` | Joined candidate observations represented for this sample/isoform/site. |
| `n_callable` | Candidate observations with an interpretable thresholded state. |
| `n_modified` | Callable observations with explicit probability at or above the analysis threshold. |
| `n_unmodified` | Explicit probabilities below threshold plus any implicit-below-threshold observations that satisfy the implicit-call rule. |
| `n_unknown` | Represented candidates excluded from the callable denominator. Always `n_candidate - n_callable`. |
| `mod_fraction` | `n_modified / n_callable` when the denominator is defined; otherwise `NA`. |
| `mean_probability` | Arithmetic mean of explicit probabilities only. Implicit unmodified observations are not inserted as zeroes. |
| `ci_low` | Lower bound of the Wilson 95% interval when `mod_fraction` is defined; otherwise `NA`. |
| `ci_high` | Upper bound of the Wilson 95% interval when `mod_fraction` is defined; otherwise `NA`. |
| `eligibility` | Derived token `eligible` only when `eligibility_reason=ok`; otherwise `ineligible`. |
| `eligibility_reason` | Stable reason described under eligibility below. |

Site state is calculated against each catalog isoform. The chromosome and strand
must match and `pos0` must lie inside an exon. Flank inference uses `context`, or
falls back to metadata `candidate_rule` when `context=NA`. The required flank is
two bases for `DRACH` (case-insensitive); for an odd-length string longer than one
containing only IUPAC DNA/RNA symbols (case-insensitive), it is half the string
length; other labels require no flank. A base inside an exon but too close
to that exon's boundary for the flank is `context_dependent`. A
chromosome/strand mismatch or a base outside all exons is `structurally_absent`.
An isoform with unknown strand is `unprojectable`. When `--reference-fasta` is
supplied, the genomic base is oriented to the site strand and compared with the
canonical base at the start of `mod_code` (`U` is compared as genomic `T`).
A mismatch is `reference_base_mismatch` and is not callable.

#### Site summary (`*.mod_site_summary.tsv`)

`mod-site-summary` streams one or more exact 27-column complete site tables and
writes:

```tsv
assay_id	analysis_threshold	sample	group	gene	site_id	chrom	pos0	strand	mod_code	context	coverage_basis	n_isoforms_total	n_isoforms_assigned	n_isoforms_present	n_isoforms_eligible	n_isoforms_site_absent	n_isoforms_context_dependent	n_isoforms_reference_base_mismatch	n_isoforms_unprojectable	n_isoforms_incomplete_candidate_universe	n_isoforms_join_rate_low	n_isoforms_low_callable	n_isoforms_other_ineligible	min_eligible_n_covering	min_eligible_n_callable	summary_state
```

`summary_state` is `shared_eligible`, `single_eligible`, or
`no_eligible_isoform` according to whether 2+, 1, or 0 isoforms are eligible.
`n_isoforms_other_ineligible` retains stable or future present-site reasons not
assigned their own summary column, including `site_join_rate_low` and
`unknown_denominator`. Eligible minima are `NA` when no isoform is eligible;
the covering minimum is also `NA` when coverage is unavailable.

#### Comparison design (`*.isoform_mod_design.tsv`)

The exact header is:

```tsv
assay_id	analysis_threshold	sample	group	gene	site_id	mod_code	isoform_id	n_modified	n_unmodified	mod_fraction	eligibility	eligibility_reason
```

This is a lossless projection of the comparison fields from the complete site
table; it includes eligible and ineligible rows. `group` and `mod_fraction` may
be `NA`. Counts are always non-negative integers, including zero. The other
columns have the meanings above. In particular, `mod_code` remains part of the
comparison key even though it is not embedded in `site_id`.

`mod-contrast` accepts a design file only when this exact header is present and
at least one row exists. Its row uniqueness key is
`(assay_id, sample, gene, site_id, mod_code, isoform_id)`. Thresholds must be
finite, in `[0,1]`, and identical for every row in an assay. A numeric
`mod_fraction` must equal `n_modified / (n_modified + n_unmodified)` within
floating-point tolerance and therefore requires a nonzero denominator.
`eligibility=eligible` requires `eligibility_reason=ok` and a numeric fraction;
`eligibility=ineligible` requires a reason other than `ok`.

### Denominator and eligibility truth tables

For a represented candidate, accounting follows this table. Here “callable site
state” means `present` or `context_dependent`, `T` is the analysis threshold, and
`E` is the source emission threshold.

| Input condition | `n_candidate` | `n_modified` | `n_unmodified` | `n_unknown` | Explicit-probability mean |
| --- | ---: | ---: | ---: | ---: | --- |
| Explicit probability, callable site state, `p >= T` | +1 | +1 | +0 | +0 | Include `p` |
| Explicit probability, callable site state, `p < T` | +1 | +0 | +1 | +0 | Include `p` |
| Implicit below threshold, callable site state, policy `low_probability`, and `T >= E` | +1 | +0 | +1 | +0 | No contribution |
| Implicit below threshold at a callable site state with a policy other than `low_probability` | +1 | +0 | +0 | +1 | No contribution |
| `unknown` observation | +1 | +0 | +0 | +1 | No contribution |
| Any represented observation at a `structurally_absent` or `unprojectable` site | +1 | +0 | +0 | +1 | Explicit probabilities are not summarized |
| No represented observation in an expanded sample/isoform/site row | 0 | 0 | 0 | 0 | `NA` |

For every row, `n_callable = n_modified + n_unmodified` and
`n_candidate = n_callable + n_unknown`. The output meaning of zero versus `NA`
is:

| Condition | Counts | `mod_fraction` and Wilson CI | Interpretation |
| --- | --- | --- | --- |
| Complete candidate universe, callable site state, `n_callable > 0` | Numeric, possibly zero | Numeric; `mod_fraction` may be exactly `0` | A denominator exists. Zero modified calls is an observed zero. |
| Same as above with `n_unknown > 0` | Numeric | `NA` | The represented denominator is partly unknown and is fail-closed. |
| `n_callable = 0` | Counts are `0` or represented as unknown | `NA` | No callable denominator exists. |
| `candidate_observations_complete=false` | Thresholded counts and explicit mean may still be numeric | `NA` | The retained candidate universe is incomplete, so a fraction is not estimated. |
| `structurally_absent` or `unprojectable` | Counts remain explicit audit integers | `NA` | Modification fraction is not biologically defined for that isoform/site. |
| `context_dependent`, complete universe, `n_callable > 0` | Numeric | Numeric | The descriptive denominator exists, but the row remains ineligible for contrast. |
| No coverage BAM supplied | `n_covering=NA` | Not applicable | `coverage_basis=unavailable`; zero must not be substituted for missing coverage. |

Eligibility is independent of whether a descriptive fraction happens to be
numeric. The first applicable rule below determines `eligibility_reason`:

| Priority | Condition | `eligibility_reason` | `eligibility` |
| ---: | --- | --- | --- |
| 1 | `site_state=structurally_absent` | `site_absent` | `ineligible` |
| 2 | `site_state=context_dependent` | `context_dependent` | `ineligible` |
| 3 | `site_state=reference_base_mismatch` | `reference_base_mismatch` | `ineligible` |
| 4 | `site_state=unprojectable` | `unprojectable` | `ineligible` |
| 5 | Present site with incomplete candidate universe | `incomplete_candidate_universe` | `ineligible` |
| 6 | Present site in a sample/assay retained with `--allow-low-join` after failing the global read-join threshold | `join_rate_low` | `ineligible` |
| 7 | Present site whose site-local observation join rate fails the same threshold | `site_join_rate_low` | `ineligible` |
| 8 | Present site with an unknown observation or unknown implicit denominator policy | `unknown_denominator` | `ineligible` |
| 9 | Present site with `n_callable < min_callable` | `low_callable` | `ineligible` |
| 10 | All preceding conditions pass | `ok` | `eligible` |

Consequently, an isoform or interaction contrast uses a site only when both
isoform rows are independently `eligible` in the same sample. Structural exon
differences and splice-dependent contexts are not converted into modification
differences.

### Explicit contrast specification

`mod-contrast --contrasts` requires this exact header:

```tsv
contrast_type	assay_id	gene	site_id	mod_code	isoform_a	isoform_b	group_a	group_b
```

Required strings cannot be empty, `NA`, or contain control characters. Optional
values use literal `NA`, not an empty field. `site_id` and `mod_code` are exact,
separate selectors. The three contrast forms are:

| `contrast_type` | `isoform_a` | `isoform_b` | `group_a` | `group_b` |
| --- | --- | --- | --- | --- |
| `isoform_effect` | required | required | optional filter or `NA` | must be `NA` |
| `condition_effect` | required | must be `NA` | required baseline group | required comparator group |
| `isoform_condition_interaction` | required | required | required baseline group | required comparator group |

At least one specification row is required. Each specification row generates
one result; duplicate specification rows are not folded. For isoform and
interaction effects, `isoform_a` and `isoform_b` must differ. For condition and
interaction effects, `group_a` and `group_b` must differ. Referencing an assay
absent from the design is fatal. Other selectors with no eligible match produce
an ineligible output row rather than an error.

### Effect-only contrast output

For output prefix `<out>`, `mod-contrast` writes
`<out>.isoform_mod_contrasts.tsv` with this exact header:

```tsv
contrast_type	assay_id	analysis_threshold	gene	site_id	mod_code	isoform_a	isoform_b	group_a	group_b	n_eligible_samples	delta_fraction	odds_ratio	interaction_delta	p_value	q_value	method	eligibility_reason
```

Only design rows with `eligibility=eligible` participate. Matching is exact on
`assay_id`, `gene`, `site_id`, and `mod_code`. No contrast pools different assay
strata. Results are sorted by contrast type, then assay, gene, site, modification
code, isoforms, and groups; specification-file order is not preserved. The
output columns are:

| Column | Meaning |
| --- | --- |
| `contrast_type` | Requested contrast family. |
| `assay_id` | Single assay stratum used for the calculation. |
| `analysis_threshold` | Assay threshold copied from the validated design. |
| `gene`, `site_id`, `mod_code` | Exact requested genomic/modification identity. `mod_code` is retained and never collapsed into or dropped from the comparison key. |
| `isoform_a`, `isoform_b` | Requested isoform identities; non-applicable `isoform_b` is `NA`. |
| `group_a`, `group_b` | Requested filters/comparison groups; non-applicable values are `NA`. |
| `n_eligible_samples` | Number of contributing eligible sample rows or within-sample isoform pairs, as defined below. |
| `delta_fraction` | Mean descriptive fraction difference for isoform or condition effects; otherwise `NA`. |
| `odds_ratio` | Descriptive pooled-count odds ratio with a `0.5` correction in all four cells; available only for estimable isoform and condition effects. |
| `interaction_delta` | Difference-in-differences for an interaction; otherwise `NA`. |
| `p_value` | Always `NA` in the current effect-only implementation. |
| `q_value` | Always `NA` in the current effect-only implementation. |
| `method` | Always `effect_only`. |
| `eligibility_reason` | `ok`, `no_shared_eligible_samples`, `missing_eligible_group`, or `missing_paired_group`. |

The implemented effect definitions are:

| Contrast | Effect | `n_eligible_samples` | Eligibility requirement |
| --- | --- | ---: | --- |
| `isoform_effect` | Mean across samples of `fraction(isoform_a) - fraction(isoform_b)` | Number of same-sample eligible isoform pairs | At least one pair. If `group_a` is supplied, it filters the pairs; both rows in a pair must have the same group. |
| `condition_effect` | `mean(group_b) - mean(group_a)` for `isoform_a` | Number of eligible rows in both groups combined | Both groups must contain at least one eligible row. |
| `isoform_condition_interaction` | `mean_group_b[fraction(a)-fraction(b)] - mean_group_a[fraction(a)-fraction(b)]` | Number of same-sample eligible isoform pairs in both groups combined | Both groups must contain at least one eligible pair, and paired rows must have the same group. |

For isoform effects, the odds ratio compares pooled counts for isoform A against
isoform B. For condition effects, it compares group B against group A. The
formula for left versus right is
`((modified_left + 0.5) * (unmodified_right + 0.5)) /
((unmodified_left + 0.5) * (modified_right + 0.5))`. It is a descriptive effect
size only: reads are not treated as biological replicates, and no inferential
test or multiple-testing adjustment is performed.
