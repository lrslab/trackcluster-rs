# `cluster` behavior (overlap/intersection mode)

Legacy Python overlap-mode clustering is implemented in `trackcluster/cluster.py` and relies on exon/intron intersections (via `bedtools intersect`) to build distance matrices.

## High-level flow (Python)
1) Build per-read exon/intron BED intervals.
2) Compute pairwise overlap bp for exon and intron intervals.
3) Convert overlaps into distances:
   - `ratio`: `1 - overlap / union_len`
   - `ratio_short`: `1 - overlap / min_len`
4) Combine exon and intron distances with `intronweight`:
   - `(D_exon + intronweight * D_intron) / (1 + intronweight)`
5) Run two passes of filtering/merging:
   - pass 1: `ratio` cutoff (`cutoff1`, default `0.05`)
   - pass 2: `ratio_short` cutoff (`cutoff2`, default `0.01`)
6) Merge reads into the retained representative by exon length, with the special-case SL score cutoff (`scorecutoff`, default `11`) controlling some short-read merges.
   - In the Python behavior, second-pass short-read collapse is gated by `score < scorecutoff`, not `<=`.

## Rust status
Rust provides overlap-mode clustering in `src/cluster/cluster_overlap.rs`, exposed through:
- `trackcluster cluster`
- `trackcluster flow --cluster-mode cluster`

Behavior:
- It performs locus splitting first (span-based), then applies the same two-pass distance/merge idea using native overlap calculations (no external tools).
- For large loci, flow/CLI can optionally pre-merge reads in batches before the final full two-pass overlap clustering (`--batch-size`, `--batch-rounds`).
- The distance defaults are `cutoff1=0.05`, `cutoff2=0.01`, and
  `intronweight=0.5` on both surfaces. The score cutoff differs by entry point:
  direct `trackcluster cluster` defaults `--sw-score` to `11`, while
  `trackcluster flow --cluster-mode cluster` defaults it to `-1`, consistent
  with the flow-wide no-SL default. Pass `--sw-score 11` to flow explicitly
  when BED score is valid SL/SW 5' evidence.
- With a non-negative `scorecutoff`, the second pass treats `score >= scorecutoff` as SL-supported. After the pair passes the exon/intron distance cutoff, an SL-supported read of a different length may collapse into a longer read or reference when their biological 5' ends differ by at most `--sl-partial-5prime-offset` (default `15` bp). A read container must itself meet the SL score cutoff and must not already be dropped; an ordinary read cannot replace the SL representative through this exception. This compares `tx_start` on plus and `tx_end` on minus; unknown strand does not enable this exception. Similar equal-length read pairs still choose one representative under the existing rule, even outside the 5' window. First-pass behavior is unchanged. With `--sw-score -1`, ordinary short-read merging runs without SL protection.
- In batched `flow --cluster-mode cluster` runs, per-gene overlap outputs use the `*_simple_coverage.bed` suffix and batch summary files use the `cluster_batch_*` prefix.

Current CLI exposure:
- `trackcluster cluster` runs overlap mode directly on one reads/reference pair.
- `trackcluster flow --cluster-mode cluster` runs the same overlap mode per gene after `preparedir`.
- There is no separate `cluster_batch` binary; the dedicated batch binary remains `clusterj_batch` for junction mode.

Compatibility is guarded by frozen legacy and scientific-truth goldens. Changes
to these merge boundaries or defaults require corresponding behavior-document
and regression-fixture updates.

## Junction-mode correction controls

Junction-mode clustering (`clusterj`, `flow` default mode, and `clusterj_batch`) first corrects low-support splice-junction sites before the SL-aware merge pass:
- `--junction-correction-offset` (default `10`; `rna002` preset `15`; `rna004` preset `10`)
- `--junction-correction-min-support` (default `5`)

The minimum support is weighted site support: a read junction site contributes `1`, and a reference junction site contributes `5`. The correction offset controls internal splice-junction coordinate snapping. It is distinct from the SL/5' and 3' terminal offsets below, which protect or merge transcript ends after junction correction. Widening junction correction can reduce rare/unused reads, but it can also erase real nearby splice sites.

`--platform-preset generic|rna002|rna004` seeds junction correction, SL, and same-junction 3' defaults. `rna002` sets junction correction offset to `15`, SL partial 5' offset to `20`, SL same-junction 5' offset to `25`, SL 5' cluster offset to `20`, SL 5' minimum support to `2`, and 3' cluster offset to `15`. `rna004` intentionally uses the conservative default cutoffs: junction correction offset `10`, SL partial 5' offset `15`, SL same-junction 5' offset `25`, SL 5' cluster offset `15`, SL 5' minimum support `2`, and 3' cluster offset `10`. Explicit CLI values override the preset.

| Preset | Junction correction offset | Junction min support | SL partial 5' offset | SL same-junction 5' offset | SL 5' cluster offset | SL min support | 3' same-junction offset | 3' cluster offset | 3' min support |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `generic` | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |
| `rna002` | `15` | `5` | `20` | `25` | `20` | `2` | `50` | `15` | `5` |
| `rna004` | `10` | `5` | `15` | `25` | `15` | `2` | `50` | `10` | `5` |

## Junction-mode SL 5' merge controls

Junction-mode clustering (`clusterj`, `flow` default mode, and `clusterj_batch`) exposes SL-aware 5' merge controls:
- `--sl-partial-5prime-offset` (default `15`)
- `--sl-same-junction-5prime-offset` (default `25`)
- `--sl-5prime-cluster-offset` (default `15`)
- `--sl-5prime-min-support` (default `2`)

When `--sw-score` is non-negative, reads with score greater than `--sw-score` are treated as SL-supported. A supported candidate with enough nearby same-junction 5' support is protected from merging when its biological 5' end is outside the relevant offset from the longer/reference track. Singleton supported reads can still merge as likely degradation.

SL information is optional and many datasets do not have it for every read. With the junction-mode default `--sw-score -1`, all reads are handled as non-SL-supported reads: they still participate in junction correction and normal 5' truncation collapsing, but they do not receive SL-cluster protection as alternative 5' isoforms. Use a non-negative cutoff such as `--sw-score 11` only when the BED score should be used as valid SL/SW 5' evidence.

## Junction-mode 3' terminal support

Same-junction reads with enough nearby 3' terminal support are retained as isoforms when their biological 3' end is outside the same-junction terminal tolerance from a compatible longer/reference track. This protects high-expression 3' early-stop isoforms that share the same splice chain as a longer isoform.

SL 5' and same-junction 3' support are computed once over the complete corrected
locus before any read batching. Each representative retains support for its
original endpoints through every batch and the final merge. Splitting a
supported terminal cluster across batches cannot erase its support, and reads
absorbed from other endpoints do not become new terminal evidence.

After support is frozen, exact same-structure reads are coalesced before batching.
The representative prefers an independently protected SL read and then stable
original metadata, not input position. Reference identity and gene/strand/structure
boundaries are preserved by this coalescing step, and every read instance remains
in its membership. Ordinary duplicate members do not become SL observations.

A protected 3' endpoint cannot be replaced by a read container without its own
frozen minimum 3' support, even if the endpoints are close. When a merge is legal,
original protected endpoint coordinates travel with the membership: every later
target must be within tolerance of all those coordinates. Several individually
small shifts cannot add up to a large displacement. SL 5' protection follows
the relevant merge-kind offset; 3' protection applies to same-junction merging
and is not transferred from a truncated chain as support for a different chain.

The CLI controls are:

- `--same-junction-3prime-offset` (default `50`): maximum distance from every retained protected original 3' endpoint to a same-junction merge target; read targets must also have independent minimum 3' support.
- `--3prime-cluster-offset` (default: active `--junction-correction-offset`): window used to sum nearby same-junction 3' support.
- `--3prime-min-support` (default `5`): minimum nearby non-reference read support required for protection.

The rule is strand-aware. On plus-strand transcripts, the biological 3' end is
`tx_end`; on minus-strand transcripts, it is `tx_start`. The minus-strand 3'
end lies on the lower-coordinate side, but an early stop truncates that side
and therefore appears as a **higher** `tx_start` (a higher genomic terminal
boundary) than the full-length isoform.

Unique counting remains catalog-aware after these isoforms are retained: reads are assigned to the closest compatible isoform by junction compatibility and terminal distance, so retained 3' early-stop reads count to the shorter terminal isoform instead of being reassigned to a longer same-junction reference.

## Standalone `clusterj` locus cap

`trackcluster clusterj` (the single-file command) reservoir-samples each overlapping locus to `--max-reads-per-locus` (default `5000`). The seed is `--downsample-seed` mixed with chrom, strand, locus span, and locus size. Dropped reads are unused; unlike `flow`/`clusterj_batch` per-gene sampling, standalone clustering does not emit a downsample TSV or scale later counts. Library callers such as `flow` pass a zero locus cap because they already downsample per gene. Setting `--max-reads-per-locus 0` is a time warning: high-diversity loci can still be expensive after the candidate-index and final-merge fixes.
