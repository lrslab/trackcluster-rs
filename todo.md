# TrackCluster Rust Revision TODO

Review date: 2026-04-25

Scope reviewed:
- Rust cluster paths: `src/cluster/clusterj.rs`, `src/cluster/cluster_overlap.rs`, `src/flow/full.rs`, `src/flow/preparedir.rs`, CLI wrappers, tests, and behavior docs.
- Original package reference: local Python package at `../trackclusterRust/trackcluster` and PyPI `trackcluster` 0.1.7 metadata.

## Current Decisions

- [x] Keep the current Rust `clusterj` many-to-many assignment behavior.
  - A truncated read may map to multiple compatible isoforms.
  - This differs from Python TrackCluster 0.1.7, which breaks after the first compatible merge.
  - Do not change `src/cluster/clusterj.rs` to first-compatible assignment.
  - Follow-up: document this as intentional Rust behavior and keep mapping-driven counts as the source of truth.

- [x] Prioritize overlap `cluster` mode before further broad parity work.
  - Immediate target: `trackcluster cluster` and `trackcluster flow --cluster-mode cluster`.
  - Main risks: slow O(n^2) merge/filter behavior, weak overlap-mode goldens, and direct cluster input loss around unmatched references.

## Current Flow Summary

Rust exposes three clustering routes:
- `trackcluster clusterj`: parse reads/reference, run junction correction and junction truncation collapse, write isoforms, `*.read_to_isoform.tsv`, and unused BED.
- `trackcluster cluster`: parse reads/reference, run two-pass exon/intron overlap clustering, write the same output shape.
- `trackcluster flow`: prepare per-gene folders, run per-gene `clusterj` or `cluster`, merge gene outputs, then run count and desc.

The high-level shape matches the Python package: `flow.py` prepares per-gene inputs, `batch.py` calls `clusterj.py`/`clustercj.py` for junction mode or `cluster.py::flow_cluster` for overlap mode, then count/description steps consume the merged isoform BED.

## Milestone 1 - Fix Overlap `cluster` Correctness First

- [x] Make reference identity in overlap `cluster` independent of `extra_fields[4] == "isoform_anno"`.
  - Current `cluster_overlap` treats a record as protected annotation only when the bigGenePred tail has `ttype=isoform_anno`.
  - CLI/docs accept ordinary BED12 references, so plain BED12 references can be treated as reads or dropped during filtering.
  - Add explicit source state to internal overlap tracks, e.g. `Track { tx, source: Reference|Read, ... }`.
  - Keep reference protection and read mapping based on source, not string metadata.
  - Files: `src/cluster/cluster_overlap.rs`, `src/cli/cluster.rs`, `src/flow/full.rs`.

- [x] Define and implement unmatched-read behavior for overlap `cluster`.
  - Current reads whose `(chrom, strand)` has no reference partition are silently skipped.
  - Proposed behavior for review: emit direct-CLI unmatched reads to `unused.bed`; in `flow`, preserve per-gene prepared-input behavior and report any impossible/unmatched records in the batch summary.
  - Add tests for wrong chromosome, wrong strand, and empty/no matching reference.
  - Files: `src/cluster/cluster_overlap.rs`, `src/cli/cluster.rs`, `src/flow/full.rs`, `tests/`.

- [ ] Add overlap-mode goldens before changing the algorithm.
  - Current `flow_overlap_mode_runs_end_to_end` mostly checks existence and total count.
  - Completed: direct `trackcluster cluster` golden for plain BED12 references plus wrong-chromosome, wrong-strand, and disjoint same-strand reads.
  - Add direct `trackcluster cluster` and `flow --cluster-mode cluster` golden comparisons for:
    - isoform BED
    - read-to-isoform TSV
    - unused BED
    - count CSV
  - Include fixtures with plain BED12 references, multi-gene input, wrong-strand reads, duplicate read names, single-exon reads, and truncation around `score == sw_score`.
  - Files: `tests/integration_flow.rs`, `tests/integration_cluster*.rs`, new `tests/golden/cluster/*`.

- [ ] Fix or document `preparedir` gene fallback for plain BED12 reads before adding plain-BED flow goldens.
  - `preparedir` currently falls back to the transcript name when a record lacks a gene-name extra field.
  - That is useful for plain BED12 references, but plain BED12 reads with no assignment can be treated as their own gene instead of novel/unassigned.
  - Direct `trackcluster cluster` now handles plain BED12 references correctly; flow-level plain-BED fixtures should use read records with explicit `none` gene fields until this is resolved.
  - Files: `src/flow/preparedir.rs`, `src/annotate/addgene.rs`, `tests/integration_flow.rs`.

- [x] Preserve overlap score-cutoff semantics.
  - Python/Rust overlap second pass collapses short reads only when `score < scorecutoff`; `score == scorecutoff` stays separate.
  - Keep this behavior while optimizing.
  - Added regression tests around `score <`, `score ==`, and `score >`.
  - Files: `src/cluster/cluster_overlap.rs`.

## Milestone 2 - Speed Up Overlap `cluster`

- [x] Add a Criterion performance baseline for overlap mode.
  - Add Criterion cases for representative cluster sizes:
    - small: `50 refs / 400 reads`
    - medium: `200 refs / 2000 reads`
    - large: `500 refs / 5000 reads`
  - Completed in `benches/perf.rs`.
  - Smoke command run on 2026-04-26: `cargo bench --bench perf -- cluster_overlap --warm-up-time 1 --measurement-time 1`.
  - Current smoke timings:
    - `50_refs_400_reads`: about `3.88 ms`
    - `200_refs_2000_reads`: about `56.4 ms`
    - `500_refs_5000_reads`: about `326.8 ms` on a focused rerun
  - Longer rerun after sparse aggregation tuning (`--warm-up-time 2 --measurement-time 3`):
    - clean HEAD baseline with same sizes: `3.34 ms`, `43.7 ms`, `232.6 ms`
    - current optimized path: `4.03 ms`, `38.1 ms`, `196.0 ms`
    - interpretation: small loci are roughly neutral/slightly slower, while medium/large dense loci improved by about 13-16%.
  - Remaining: add one `flow --cluster-mode cluster` fixture/timing smoke test or document a manual flow benchmark command.
  - Files: `benches/perf.rs`, `docs/behavior/cluster.md`, `tests/integration_flow.rs`.

- [x] Replace all-pairs overlap filtering with sparse candidate generation for large loci.
  - Current `filter_pass` loops over every pair in each locus and computes exon/intron overlaps directly.
  - Implemented a local exon-interval sweep that aggregates exon overlap bp per track pair and evaluates only exon-overlap candidates for normal low cutoffs.
  - The code falls back to exact all-pairs scanning for small loci and for high cutoffs where no-exon-overlap pairs can still pass through intron overlap.
  - Intron overlap is now computed only for retained sparse candidates on large low-cutoff loci.
  - Added regression coverage for sparse ordering and high-cutoff fallback behavior.
  - Files: `src/cluster/cluster_overlap.rs`.

- [ ] Add cheap pruning before distance calculation.
  - Use length-ratio bounds to skip pairs that cannot satisfy `cutoff1` or `cutoff2` even with perfect containment.
  - Skip impossible intron comparisons when either track has no introns.
  - Note: the sparse-candidate pass already avoids intron work for pairs with no exon-overlap candidate under normal cutoffs.
  - Keep exact current distance formulas for retained candidate pairs.
  - Files: `src/cluster/cluster_overlap.rs`.

- [ ] Make batching in overlap mode reduce worst-case work.
  - Current overlap batching still ends with a full final `cluster_once`, which can restore O(n^2) behavior on large loci.
  - Batch-size timing on the current 500-ref/5000-read synthetic locus:
    - `batch_size=0`: about `197 ms`
    - `batch_size=250`: about `1.19 s`
    - `batch_size=500`: about `1.77 s`
    - `batch_size=1000`: about `3.03 s`
  - 50k-read release-mode probe data, 500 refs, generated single-chromosome overlap workload:
    - `locus_span=500000`: `batch_size=0` `1.76 s`; `500` `4.26 s`; `2000` `8.63 s`.
    - `locus_span=100000`: `batch_size=0` `7.67 s`; `500` `10.80 s`; `2000` `22.66 s`.
    - Batching changed output counts slightly in the 50k runs, so it is not only a performance knob.
    - No large real BED fixture is available in this repo or the adjacent local TrackCluster checkout; these are synthetic stress measurements.
  - Oversized batch follow-up after fixing `read_count <= batch_size` to use the no-batch path:
    - `locus_span=500000`: `batch_size=0` `1.78 s`; `100000` `1.84 s`; output counts matched.
    - `locus_span=100000`: `batch_size=0` `7.40 s`; `100000` `7.55 s`; output counts matched.
    - Criterion 5k-read case for `batch_size=100000`: about `188.8 ms`.
  - 120k-read ignored release-mode speed test:
    - Command: `cargo test --release --test perf_overlap_large -- --ignored --nocapture`.
    - `locus_span=500000`: `batch_size=0` `9.05 s`; `500` `14.64 s`; `100000` `12.26 s`.
    - `batch_size=100000` is slower than no batching once read count exceeds 100k because it performs one large intermediate merge before the final merge.
  - Recommendation before changing defaults: keep direct overlap `cluster` at `batch_size=0`; consider changing `flow --cluster-mode cluster` to use an overlap-specific default of `0` while keeping junction `clusterj` at `500`.
  - After sparse candidate generation is in place, verify whether final full merge is acceptable.
  - If still slow, implement a bounded final reconciliation pass using the same sparse candidate index rather than all-pairs scanning.
  - Files: `src/cluster/cluster_overlap.rs`, `src/flow/full.rs`.

- [ ] Keep output semantics stable while optimizing.
  - `read_to_isoform.tsv` must remain deterministic.
  - `name2-mode coverage|full|none` must keep current meaning.
  - Read assignment can remain many-to-many where overlap semantics currently produce it.
  - Completed: added single-thread vs multi-thread determinism tests for overlap mode.
  - Files: `src/cluster/cluster_overlap.rs`, `tests/`.

## Milestone 3 - Fix `flow --cluster-mode cluster`

- [x] Ensure `flow --cluster-mode cluster` uses the optimized overlap core.
  - Per-gene clustering should call the same corrected/optimized code as direct `trackcluster cluster`.
  - Preserve merged output naming: `*_simple_coverage.bed`, `*_unused.bed`, `*_read_to_isoform.tsv`, `<prefix>_isoform.bed`, `<prefix>_unused.bed`, `<prefix>_read_to_isoform.tsv`.
  - Fixed per-gene skip/reuse logic so a rerun without `--force` only skips when isoform, unused, and read-to-isoform outputs all exist.
  - Added an overlap-flow rerun regression test that removes a per-gene unused BED and confirms it is regenerated.
  - Files: `src/flow/full.rs`.

- [ ] Review overlap-mode flow defaults.
  - Upstream `trackrun.py cluster` defaults `--batchsize 2000` and exposes `--cutoff2 0.05`.
  - Rust direct `cluster` defaults to no batching and `cutoff2=0.01`; Rust `flow` has its own batching defaults.
  - Current recommendation from benchmarks: keep `cutoff1=0.05`, keep `cutoff2=0.01` for now, and set overlap-mode batch size to `0` by default unless memory pressure on extreme genes requires an explicit cap.
  - If an explicit cap is needed for a very large overlap gene, prefer `--batch-size 500` over `2000`; it was consistently faster in the 50k-read probe. Treat that as a compatibility/memory fallback, not the speed default.
  - Rationale: current overlap batching is slower on the tested synthetic large locus because it repeats intermediate merges and still performs the final full reconciliation.
  - Decide after benchmarks whether to change defaults, add a Python-compatible preset, or keep current defaults with clearer docs.
  - Files: `src/cli/cluster.rs`, `src/cli/flow.rs`, `docs/CLI.md`, `docs/PIPELINE.md`.

- [ ] Add flow-level reporting for overlap-mode large genes.
  - Summary should make it clear when batching, downsampling, or unmatched-read handling affected a gene.
  - Keep heartbeat/progress behavior useful for long overlap runs.
  - Files: `src/flow/full.rs`.

## Milestone 4 - Document Accepted Junction-Mode Differences

- [ ] Document `clusterj` many-to-many assignment as intentional.
  - Update CLI/behavior docs to say a read can map to multiple compatible isoforms and counts are split by occurrence.
  - Make the existing two-reference truncation fixture explicitly cover this contract.
  - Files: `docs/CLI.md`, `docs/behavior/cluster.md`, `tests/golden/clusterj/*`.

- [ ] Move remaining `clusterj` correctness issues after overlap mode.
  - Reference identity should eventually become source-based in `clusterj` too.
  - Direct `clusterj` unmatched-read behavior still needs a decision.
  - Minus-strand and single-exon semantic tests should still be broadened.
  - Files: `src/cluster/clusterj.rs`, `tests/`.

## Later Backlog

- [ ] Reconcile upstream flow output set.
  - Upstream flow writes `<prefix>_cov5_isoform.bed` and optionally novel-gene isoform outputs.
  - Rust currently omits these.

- [ ] Reconcile count output compatibility.
  - Upstream count writes `<prefix>_exp.csv` with `gene,isoform,coverageall,<groups...>`.
  - Rust single-sample count writes `isoform_id,count`.

- [ ] Reconcile fusion classification evidence.
  - Upstream `desc` identifies fusion from existing isoform `geneName` values split on `||`.
  - Rust recomputes fusion by overlap against references.

- [ ] Reduce string-heavy subread propagation.
  - Current `HashSet<String>` cloning is expensive.
  - Prefer interned read IDs or `Arc<str>` first, then consider numeric IDs or deferred membership materialization.

## Verification Checklist

Last run on 2026-04-26 for the overlap correctness + sparse-candidate speed pass:
- `cargo fmt --all`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all --all-features`
- `cargo bench --bench perf -- cluster_overlap --warm-up-time 1 --measurement-time 1`

Run after each behavior or optimization change:
- [ ] `cargo fmt --all`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all --all-features`
- [ ] `cargo bench --bench perf` when overlap clustering hot paths change
  - Last overlap smoke run: `cargo bench --bench perf -- cluster_overlap --warm-up-time 1 --measurement-time 1`
- [ ] Regenerate or review goldens only after deciding whether each difference is parity or intentional Rust behavior.
