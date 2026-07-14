# Performance policy

The default and high-diversity `clusterj` Criterion cases use the production
defaults (`sw_score=-1`, coverage-only `name2`) and retain a separately labeled
SL/full-payload scenario. Low-level overlap cases use `ClusterOptions` defaults
(`sw_score=11`, full `name2`) unless the individual benchmark explicitly
overrides `name2` to coverage mode. The `clusterj_high_diversity` case keeps
thousands of mostly nonmergeable reads in one overlapping locus to expose
repeated batching and reduction costs.

The scheduled `performance.yml` workflow builds with `--locked`, runs the full
Criterion suite, and runs the 120k-read overlap probe under `/usr/bin/time -v`.
Criterion output and peak-RSS logs are retained as run artifacts for 90 days.

Overlap clustering switches loci of at least 512 tracks to a sweep-based
exon-overlap candidate index when the active cutoff makes zero-exon-overlap
pairs mathematically impossible to merge. Smaller loci retain the all-pairs
path. Large loci also retain the quadratic fallback for parameter combinations
where an intron-weighted distance can legitimately merge a pair without exon
overlap; replacing that path with the same index would change scientific
results. The scheduled high-diversity measurement is the gate for introducing a
broader hierarchical approximation in a separately versioned mode.
