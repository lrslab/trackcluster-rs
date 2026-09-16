# Performance policy

The targeted single-locus `clusterj` benchmarks increase molecule abundance
while keeping the splice-site catalog fixed. They use one thread, batches of up
to 500 read representatives, no downsampling, and the production defaults for
junction correction, terminal retention, and output (`sw_score=-1`, coverage-only `name2`). Both
start from the same six-exon locus with explicitly specified exon coordinates.

| Benchmark | Input structures | Read counts | Required output |
| --- | --- | --- | --- |
| `clusterj_abundance/four_fixed_isoforms` | Full chain, exon 3 skipping, exon 4 skipping, and exons 3 + 4 skipping; abundance ratio 60:25:10:5 | 2,000 / 10,000 / 20,000 | The four catalog structures, every read mapped once, zero unused |
| `clusterj_terminal_variation/fixed_junctions_two_three_prime_clusters` | One unchanged internal splice chain; two supported 3' clusters 160 bp apart in an 80:20 mixture; 0–14 bp terminal trimming within each cluster | 2,000 / 10,000 / 20,000 | Two terminal isoforms, every read mapped once, zero unused |

Only the full chain is supplied as a reference. The skipping isoforms must
survive as novel structures. The terminal fixture changes only `tx_start` and
`tx_end`, with more than 100 distinct endpoint pairs; assertions verify that
every internal donor and acceptor stays at its catalog coordinate. Thus one
case exercises repeated reads over a small isoform catalog and the other
exercises merging endpoint variants and retaining supported terminal clusters.
Increasing read count never introduces additional internal splice sites.

Junction clustering freezes terminal support before coalescing exact read
structures. During merging, endpoint evidence is updated immediately, but read
memberships are expanded only for retained representatives by traversing the
recorded merge sources. This avoids copying large sets through every intermediate
container without changing the pair decisions or multi-target mapping semantics.
The source graph costs O(N + E) memory per batch; dense unbatched loci can still
have quadratic edge counts, so this is not a general linear-time guarantee.

These are controlled synthetic workloads. The selected abundance mixtures
and terminal offsets are not fitted to an experimental dataset. Real-locus
measurements are still needed to determine which mechanism dominates a
particular dataset. The older random-coordinate `clusterj` grouping and overlap
cases remain synthetic algorithm probes, not biological performance models.

The former `clusterj_high_diversity` and derived rare-filter benchmark have
been removed. Their middle exon moved by 7 or 40 bp for each new structure,
creating thousands of distinct splice positions in one locus. In the 7-bp
version, per-coordinate support was one (default threshold five), so only
the reference-matching read and its 7-bp neighbor survived correction:
1,998 of 2,000 reads were filtered. Repeating each structure five times in
the 40-bp version fixed survival but did not fix the unsuitable splice model.
Neither case supports a claim about performance on high-abundance genes.

Validate the targeted fixtures and their output assertions with:

```bash
cargo bench --locked --all-features --bench perf -- 'clusterj_(abundance|terminal_variation)' --test
```

Remove `--test` to collect Criterion timing samples; a fixture-validation run
alone is not a performance measurement.

The scheduled `performance.yml` workflow builds with `--locked`, runs the full
Criterion suite, and runs the 120k-read overlap probe under `/usr/bin/time -v`.
Criterion output and peak-RSS logs are retained as run artifacts for 90 days.

Overlap clustering switches loci of at least 512 tracks to a sweep-based
exon-overlap candidate index when the active cutoff makes zero-exon-overlap
pairs mathematically impossible to merge. Smaller loci retain the all-pairs
path. Large loci also retain the quadratic fallback for parameter combinations
where an intron-weighted distance can legitimately merge a pair without exon
overlap; replacing that path with the same index would change scientific
results. Introducing a broader hierarchical approximation in a separately
versioned mode requires real-locus measurements as well as the synthetic
regression benchmarks.
