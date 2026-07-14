# Rust API policy

TrackCluster-rs is primarily a command-line application. The supported,
semver-stable library surface consists of the core model types and their
validating constructors, strict BED/manifest I/O, and interval algorithms
documented by the crate. Public clustering, counting, annotation, and flow modules are packaged
for compatibility with the binaries, but may evolve between minor releases.

`Interval::new`, `Transcript::new`, and the strict BED readers enforce their
documented invariants at construction time. Model fields remain public and
mutable for source compatibility, so direct struct literals or later field
mutation can bypass that validation. Library callers should prefer the
constructors and readers and must preserve the invariants when changing public
fields directly; model-consuming algorithms assume they remain valid.

The Clap command tree, sample-name encoding helpers, worker implementation, and
artifact plumbing are internal. Packaged binaries enter through the narrow
`trackcluster_rs::run_cli_from_env` wrapper, so CLI implementation details do not
need to remain public.

The crate enables `warn(missing_docs)`. The core model/I/O/interval surface is
fully documented and includes a compiling BED-loading example. Callers that
choose to use the compatibility pipeline modules should compose validated
`PrepareConfig`, `ClusteringConfig` (including junction and overlap settings),
`CountingConfig`, `DownsampleConfig`, and `RuntimeConfig` instead of passing
long positional option lists, while treating those types as a minor-version
compatibility surface rather than the stable core API.

`RuntimeConfig::gene_error_policy` defaults to `GeneErrorPolicy::Continue`: an isolated gene
failure is reported and excluded while verified genes remain available for downstream work.
Library callers that require all-or-nothing batches can select `GeneErrorPolicy::Strict`. Global
configuration, preparation, infrastructure, integrity, and output-publication errors remain fatal.

BED extension columns remain round-trippable, but supported code should use
`Transcript::metadata()`, `Transcript::metadata_mut()`, or `BigGenePredAttrs`
instead of numeric `extra_fields` indices. `Transcript::geometry()` exposes an
immutable geometry-only view when an algorithm must not depend on annotations.

## Feature flags

- `index-binned`: enables a reusable fixed-bin transcript-span index. This is
  useful when many query collections are intersected against one catalog. The
  default sweep backend remains available and accepts unsorted input.
