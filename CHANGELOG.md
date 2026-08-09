# Changelog

## Unreleased

## 0.3.0

## 0.3.0

- Add the isoform-level RNA-modification V1 commands: `mod-import-dorado` for strict MM/ML/MN decoding from primary genome-aligned modBAM, `mod-import-m6anet` for RNA002 read-probability projection, `mod-aggregate` for unique read-to-isoform aggregation, `mod-site-summary` for deterministic site QC inventories, `mod-contrast` for explicit effect-only comparisons, and `mod-subsample` for synchronized technical coverage partitions. `flow --mod-manifest` can run aggregation, and optional contrasts, after final unique assignment.
- Standardize caller-import sidecars as `<prefix>.observations.tsv`, `<prefix>.assay.json`, and `<prefix>.import_qc.tsv`. Aggregation publishes `<prefix>.mod_join_qc.tsv`, `<prefix>.mod_site_join_qc.tsv`, `<prefix>.isoform_mod_sites.tsv`, and `<prefix>.isoform_mod_design.tsv`; site summaries and contrasts publish `<prefix>.mod_site_summary.tsv` and `<prefix>.isoform_mod_contrasts.tsv`. Subsampling emits ready-to-run sample/modification manifests, synchronized reads, assignments, observations, optional coverage BAMs, QC, checksums, and provenance in one output directory.
- Keep V1 modification fractions fail-closed: missing or unknown calls are not treated as unmodified, source and analysis thresholds remain distinct, incomplete candidate universes and low global/site join rates make rows ineligible, and optional exact BAM coverage plus strand-oriented FASTA checks audit denominators and canonical bases.
- Limit V1 to exact unique read-to-isoform assignments and compatible caller/model/chemistry assay strata. Contrasts are descriptive only (`p_value` and `q_value` are `NA`), technical pseudo-samples are not biological replicates, splice-junction motif contexts remain ineligible, the Dorado importer and normalized observation reader currently materialize data in memory, and related output sets are atomic per file rather than one multi-file transaction.

## 0.2.0

This is a breaking pre-1.0 release. Before upgrading an existing workflow,
review the [0.2.0 output migration](docs/FORMATS.md#020-output-migration),
especially the versioned novel-isoform/name2 identities and count/description
schemas.

- Add `bam2bigg` (`bam-to-bigg` alias) for streaming pure-Rust conversion of genome-aligned BAM records to TrackCluster bigGenePred-compatible BED12+8, with a default MAPQ cutoff of `30`, BAM-stem sample grouping, and opt-in secondary/supplementary retention.
- Add `gff2bigg` (`gff-to-bigg` alias) for deterministic conversion of GFF3 or GTF exon annotations to validated reference BED12+8 catalogs, with quote-aware syntax auto-detection, strict GFF3 parent identity checks, and configurable GFF3 gene-label attributes.
- Add `export` for GTF 2.2, GFF3, and SQANTI3 input-audit TSV generation from a BED transcript catalog.
- Expand `validate-bed` with an explicit `--lenient` legacy-repair subset and a schema-versioned `--report` summary while retaining strict validation by default.
- Define the converter output contract: BAM CIGAR `N` and annotation exons form BED blocks; converter outputs carry all eight TrackCluster metadata fields and intentionally do not project CDS, UTR, or annotation phase.
- Publish converter outputs atomically so malformed BAM/GFF/GTF input cannot replace a previous successful destination.
- Version ambiguous output contracts: count CSV is RFC 4180 with `gene,isoform_id,count`, description tables use `trackcluster-description-v2` headers, novel isoforms use lossless `tc_novel_v1:` structural IDs, and full `name2` payloads use percent-encoded `tc_name2_v1:` values.
- Enforce catalog and molecule identity at merge/count boundaries, including duplicate/reserved reference-ID rejection, globally unique structural novel IDs, idempotent duplicate mapping rows, and rejection of conflicting alignments for one molecule ID during unique assignment.
- Encode validated biological gene IDs into deterministic filesystem keys, publish versioned biological-ID/key maps and per-directory identity markers, keep `--gene-list`/`--downsample-gene` in the biological-ID namespace, and reject keys reserved by top-level preparation, flow, or batch artifacts before mutation.
- Add versioned per-gene completion manifests with request, input, tool, seed, and output-integrity fingerprints; reuse only exact verified results, publish the manifest last, and make count-only verification all-or-nothing.
- Add record-local invalid-read recovery with auditable rejected-read TSVs, plus recoverable gene-local batch failures that exclude stale artifacts under the default continue policy; strict read and gene-error modes remain available, and batch error reports are emitted only when errors occur.
- Make binary release archives self-contained with the changelog, core offline documentation, and synthetic runnable examples; release builds now depend on the reusable full CI workflow and smoke-test only unpacked archive contents.
- Keep mode-specific count artifacts exact across reruns: unique modes publish auditable mapping/provenance files, while successful fractional reruns remove stale unique-only outputs.
- Make count-only reuse fail closed on empty or conflicting prefix-scoped gene metadata, use a versioned cluster-batch gene map instead of directory discovery when prefix metadata is absent, rebuild scale factors from manifest-verified per-gene downsampling records, and refuse independently downsampled multi-gene molecules rather than publish biased abundance estimates.
- Preserve structurally distinct alignments that share a read ID across junction-batching boundaries, reject junction snaps that would create an empty exon, and make standalone multi-sample counting tag delimiter-bearing raw read IDs exactly as manifest flow does.
- Preserve leading/trailing whitespace in read-to-isoform TSV identities and use overflow-safe terminal-distance comparisons across the full BED coordinate domain.
- Reject input/output and pairwise output aliases before standalone commands publish files, stage prepared generations behind an authoritative gene-list commit marker, require explicit batch gene selection outside inline preparation, and reject symlink-substituted cache manifests, artifacts, or atomic temporary files.
- Treat `flow`, `preparedir`, and output-root recount directories as pipeline-owned trees: external inputs must stay outside them, and pre-existing symlink or hard-link aliases are rejected before any mutation; bucket scratch directories are uniquely reserved and cleaned on failure.
- Preserve the source commit from Cargo's packaged VCS metadata, keep clean Git builds identifiable, and hash the actual tree for Cargo-package or dirty-checkout builds so edited sources cannot share cached per-gene results.
- Remove the non-portable internal 488 walkthrough from public documentation and keep maintainer/design documents out of the crates.io source package.

## 0.1.18
- Change the junction-mode default `--sw-score` for `clusterj`, `flow`, and `clusterj_batch` from `11` to `-1`, so no-SL/no-valid-5'-score datasets use ordinary truncation merging by default instead of treating BED scores as SL/SW 5' evidence.
- Keep SL/SW-supported 5' protection available for scored datasets by passing an explicit non-negative cutoff such as `--sw-score 11`; platform presets still control junction, SL 5', and 3' offsets but no longer opt into score-based 5' protection implicitly.
- Fix terminal single-exon containment in junction-mode merging, including exact single-exon matches, terminal exon boundary cases, and minus-strand genomic containment checks.
- Add regression coverage showing that high-score terminal single-exon reads merge under the no-SL default, while explicit `--sw-score 11` can still retain supported terminal SL clusters.
- Update README, CLI, pipeline, and behavior docs for the new no-SL default and explicit `--sw-score 11` guidance.
- Refine the m6A/modification analysis plan with call-state/provenance fields, negative-call completeness handling, read-id audit requirements, site eligibility/QC outputs, and transcript-to-genomic projection constraints.

## 0.1.17
- Interpret `--sw-score -1` as "no SW-supported 5' signal" for junction and overlap clustering instead of disabling ordinary truncation merging; non-SL-supported reads still merge through the normal cluster rules, including batched runs.
- Speed up catalog-aware unique assignment by indexing isoforms by gene/chromosome/strand span bins and scanning retained-intron regions with sorted intron bounds.
- Keep unique assignment inside the inferred gene domain: reads with an available gene no longer fall back to same-locus isoforms from other genes when no overlapping in-gene catalog candidate exists.
- Prefer an original read-to-isoform mapping over a catalog-expanded isoform when unique assignment scores are exactly tied, avoiding lexical tie switches.
- Clarify manual `preparedir` -> `clusterj_batch` docs: use the prepared gene root as the batch output root when the result will be consumed by `trackcluster count --output-root`.
- Add regression coverage for negative `--sw-score` behavior, catalog tie handling, gene-domain fallback, and `-1` CLI parsing.

## 0.1.16
- Preserve supported same-junction 3' terminal clusters in junction-mode clustering, so high-support 3' early-stop isoforms are not collapsed into longer/reference isoforms solely because the splice chain matches.
- Expose 3' terminal protection controls in `clusterj`, `flow`, and `clusterj_batch`: `--same-junction-3prime-offset`, `--3prime-cluster-offset`, and `--3prime-min-support`.
- Change same-junction 3' terminal defaults to require at least `5` nearby supported reads and a biological 3' end more than `50` bp from the merge target; the nearby 3' support window defaults to the active junction correction offset.
- Apply the 3' terminal support rule strand-aware, including minus-strand genes where the biological 3' end is the lower genomic coordinate.
- Keep unique counting aligned with retained 3' early-stop isoforms: catalog-aware unique assignment chooses the closest compatible terminal structure instead of assigning those reads back to a longer/reference same-junction isoform.
- Add regression coverage for retained same-junction 3' clusters and minus-strand 3' early-stop unique counting.

## 0.1.15
- Change `flow` unique counting to select reads directly inside each per-gene output folder, using `{gene}_nano.bed`, `{gene}_simple_coverage[j].bed`, and `{gene}_read_to_isoform.tsv` instead of falling back to manifest or original read inputs.
- Add `flow --count-only` to reuse completed per-gene cluster outputs and rerun only merge, count, multi-sample count, and desc outputs.
- Emit `<prefix>_read_to_isoform.unique.tsv` in unique assignment mode so the exact mapping used for final counts is available for audit/reproduction.
- Skip genes with incomplete per-gene count inputs during unique assignment, rather than failing the whole final count step.
- Fix unique assignment so catalog expansion does not resurrect unused reads or reads excluded by per-gene downsampling.
- Add regression coverage for count-only reruns, unique audit mapping output, missing per-gene inputs, and downsample-scaled final counts.

## 0.1.14
- Derive `count-multi` aggregate counts from the emitted per-sample matrix and write them to `<prefix>.isoform_count.csv`.
- Synchronize `flow --manifest` `<prefix>_isoform_count.csv` from the same multi-sample aggregate count, so the total count table, sample matrix, long usage table, and group usage table share one assignment result.
- Keep `unique` assignment as the default for `count-multi` and `flow`, including catalog-expanded candidate selection for reads that should assign to a closer isoform than the original mapping candidate.
- Document the new aggregate count output and its invariant: each `*.isoform_count.csv` value is exactly the sum of the corresponding row in `*.isoform_counts.matrix.tsv`.
- Add regression coverage that checks `count-multi` aggregate counts and `flow --manifest` main count outputs against the per-sample matrix.

## 0.1.12
- Make `unique` the default assignment mode for `count`, `count-multi`, and `flow`; pass `--assignment-mode fractional` to keep legacy split-count behavior.
- Expand unique counting against the isoform catalog before choosing the closest compatible isoform, so reads are not trapped by an incomplete `read_to_isoform.tsv` candidate set.
- Use structure-aware unique assignment across mapping-backed and embedded-subread counting paths, including multi-sample manifest reads and pooled `flow` outputs.
- Add fuzzy same-junction merging in `clusterj` using the active junction correction offset, with a same-length junction-chain index to keep candidate scans bounded.
- Raise the default weighted junction correction minimum support to `5`; the `rna002` preset now uses junction correction offset `15` with SL 5' offsets `20/25/20`.
- Update CLI and behavior documentation for catalog-aware unique counting, weighted junction support, RNA002 defaults, and the distinction between junction correction and SL/5' terminal offsets.
- Add regression coverage for fuzzy same-junction merging, catalog-expanded unique assignment, default unique counting, and explicit fractional compatibility.

## 0.1.11
- Add configurable junction correction controls for `clusterj`, `flow`, and `clusterj_batch`: `--junction-correction-offset` and `--junction-correction-min-support`.
- Add `--platform-preset generic|rna002|rna004`; `generic` preserves current defaults, `rna002` widens junction correction and SL 5' offsets for RNA002/DEI-style workflows, and `rna004` keeps conservative/default RNA004 cutoffs.
- Let explicit junction correction and SL 5' options override preset values across single-file, flow, and batched clustering.
- Record the active platform preset and junction correction settings in batch summary files for reproducibility.
- Document the distinction between internal junction correction offsets and SL/5' terminal merge/protection offsets.
- Add regression coverage for junction correction offset/min-support behavior, preset expansion, preset override precedence, and summary output.

## 0.1.10
- Speed up SW-aware `clusterj` merging by grouping SL 5' support by junction chain and avoiding redundant merge-target scans.
- Add exact-duplicate representative pruning for non-reference reads with identical corrected transcript structure.
- Add `unc52` regression fixtures for SW-aware junction clustering.

## 0.1.9
- Add SL-aware junction merge controls for `clusterj`, `flow`, and `clusterj_batch`: `--sl-partial-5prime-offset`, `--sl-same-junction-5prime-offset`, `--sl-5prime-cluster-offset`, and `--sl-5prime-min-support`.
- Keep supported alternative SL 5' clusters as candidate isoforms while still merging singleton likely-degradation reads.
- Honor `--sw-score -1` during batched junction merging so truncation collapsing is fully disabled when requested.
- Record SL merge settings in batch summary files and add regression coverage for the new merge boundaries.
- Add an ignored large overlap speed probe for manual 100k+ read performance checks.

## 0.1.8
- Fix overlap `cluster` reference handling so plain BED12 references are protected by input source instead of requiring `ttype=isoform_anno`.
- Report overlap-mode reads with no matching reference chromosome, strand, or locus in `unused.bed` instead of silently dropping them.
- Speed up large overlap loci with sparse exon-overlap candidate generation while preserving exact all-pairs behavior for small or high-cutoff cases.
- Fix `flow --cluster-mode cluster` reruns so per-gene outputs are regenerated unless isoform BED, unused BED, and read-to-isoform TSV are all present.
- Add overlap-mode goldens, unused-read regressions, deterministic threading coverage, and batch-size benchmarks.

## 0.1.7
- Add legacy overlap/intersection clustering to `flow` via `--cluster-mode cluster`.
- Add overlap-mode CLI controls for `cluster`/`flow`: `--batch-size`, `--batch-rounds`, `--sw-score`, cutoff tuning, and `--name2-mode`.
- Keep second-round `SL` reads as their own track when `score == --sw-score`, matching the original TrackCluster boundary behavior.
- Release hygiene: add checked-in license texts, pin CI/release builds to Rust `1.90.0`, and keep release tarballs out of the main branch.

## 0.1.6
- Performance: speed up `clusterj` 5' truncation collapsing on large loci by indexing junction suffixes (avoids quadratic scans).
- Bench: add a large single-locus `clusterj` benchmark to track this workload.

## 0.1.5
- Default `--name2-mode` is now `coverage` (smaller isoform BEDs; rely on mapping TSVs for counting).
- Default `--max-reads-per-gene` is now `50000` for `flow`/`clusterj_batch` (memory-friendly; set `0` to disable). Counts/usage tables are scaled when downsampling occurs.
- Add `--heartbeat-seconds` / `--heartbeat-top` to periodically report progress and in-flight genes during `flow`/`clusterj_batch`.
- `count` and `count-multi` auto-discover `*_read_to_isoform.tsv` next to the isoform BED when present.

## 0.1.4
- Add `--name2-mode` (full/coverage/none) to control isoform `name2` payload size while keeping a read-to-isoform mapping.
- Add `--read-to-isoform` fast path for `count` and `count-multi` to reuse mapping TSVs from `flow`/`clusterj`.
- Add per-gene downsampling for `flow`/`clusterj-batch` (`--max-reads-per-gene`, `--downsample-gene`, `--downsample-seed`) and scale counts/usage tables via `clusterj_batch_downsample.tsv`.
- Change manifest mode: `<prefix>_pooled_reads.bed` is now written only when `flow --emit-pooled-reads` is set.
- Performance: bucketed two-pass `preparedir` for large inputs + faster BED12 reading.
- Default `--sw-score` is now `11` (TrackCluster Python default); use `-1` to disable collapsing.

## 0.1.3
- Add `count-multi` subcommand for per-sample isoform usage from pooled isoforms
- Add `flow --manifest` mode:
  - pool reads from a sample manifest (`<prefix>_pooled_reads.bed`)
  - cluster once into shared isoforms
  - emit multi-sample usage tables (`<prefix>.isoform_usage.long.tsv`, `<prefix>.isoform_counts.matrix.tsv`, optional group table)
- Add manifest TSV parser (`sample`, `reads`, optional `group`) with strict validation and tests
- Add integration fixtures/tests for multi-sample counting and flow manifest mode
- Release: build Linux x86_64 artifact as `x86_64-unknown-linux-musl` to avoid host glibc version mismatch errors

## 0.1.0
- Initial Rust CLI: `validate-bed`, `clusterj`, `cluster`, `count`, `addgene`, `desc`, `preparedir`
- `flow` subcommand: one-command end-to-end pipeline (preparedir + clusterj batch + count + desc)
- Native interval utilities (no runtime shell-out)
- Small fixtures + golden-based integration tests for `clusterj` and `count`
- Pin Rust 1.90.0 via `rust-toolchain.toml`
- CI: lint/test workflow + automated release with pre-built binaries for Linux and macOS
