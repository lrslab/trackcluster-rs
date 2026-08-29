# Isoform-modification validation strategy

The modification workflow has three distinct correctness boundaries:

1. caller-format decoding and coordinate projection;
2. read-to-isoform joining, denominator/QC handling, and effect calculation; and
3. biological modified-base accuracy.

No single fixture validates all three. The project therefore uses a layered
strategy and keeps large public data out of normal CI.

## Validation tiers

| Tier | Data | What it validates | When to run |
| --- | --- | --- | --- |
| Specification/unit | Hand-built MM/ML/MN, CIGAR, GTF, FASTA, and normalized rows | Strict parsing, reverse orientation, all CIGAR operations, duplicate handling, skip semantics, schema invariants | Every CI run |
| Realistic simulation | Six samples, two conditions, two splice isoforms, exact coverage BAMs, unmatched observations, and known integer effects | Sample namespaces, joins, coverage, isoform/site expansion, eligibility, fractions, and all three descriptive contrast types | Every CI run |
| Pinned m6Anet fixture | Official m6Anet 2.1.0 test output at commit `590ec277...` | Real CSV/gzip dialect, 5,595 read-site rows, opaque read indexes, `data.info` retention audit, site-probability cross-check, and deterministic projection arithmetic | Before release or importer changes |
| m6Anet reference audit | The fixture's 82 candidate-universe transcript IDs versus pinned Ensembl 91 or GENCODE 39 | Exact transcript-version provenance; prevents synthetic/mixed fixture IDs from being assigned biological coordinates | Before claiming genomic validation |
| ONT all-5-mer truth | ONT `rna-mod-validation-all5mer-2026.07` control-A and m6A BAM prefixes plus all 256 truth sites | Real Dorado 2.0 multi-mod MM/ML/MN, clipping/projection, sparse calls, and expected modified/unmodified direction | Before release or Dorado importer changes |
| Biological locus panel | A bounded HCT116/SG-NEx locus with final TrackCluster isoforms and matched m6Anet calls | End-to-end behavior on natural splice variation and heterogeneous coverage | Manual milestone validation |

The integration contrast is deliberately descriptive. V1 emits no inferential
`p_value` or `q_value`; pseudo-samples created by `mod-subsample` are technical
coverage partitions and are not biological replicates.

For analysis-facing validation, run aggregation with
`--eligibility-profile strict`. A row then needs an indexed reference FASTA,
exact primary-alignment coverage, configured covering/callable minima,
candidate/covering and callable/covering rate gates, complete/known
denominators, and passing global assigned-read join QC. Site-local observation
assignment rates remain audit-only. Dorado rows additionally require version,
explicit model, and source-threshold provenance recovered consistently from one
BAM `@PG` record. The shipped defaults (20 covering, 10 callable, and 0.8 for both
rates) are conservative engineering guardrails that must be calibrated against
the selected caller/model/chemistry controls.

## Reproduce the public-data checks

The opt-in script caches downloads and writes reports under `target/`, which is
not committed:

```bash
# Small official m6Anet fixture; uses a deterministic synthetic projection
# catalog to validate formats and projection arithmetic.
scripts/validate_modification_real_data.sh m6anet

# Audit the fixture against real annotations. These commands intentionally
# exit nonzero because the fixture IDs do not all match one annotation; they
# still write reference_compatibility*.tsv with the exact reasons.
M6ANET_REFERENCE_MODE=ensembl91 \
  scripts/validate_modification_real_data.sh m6anet
M6ANET_REFERENCE_MODE=gencode39 \
  scripts/validate_modification_real_data.sh m6anet

# Stream only the first 1,000 records from each 108–123 MB public ONT BAM.
# Requires samtools; the 1.7 TB full POD5 collection is not downloaded.
scripts/validate_modification_real_data.sh ont

# Increase the bounded BAM prefix for a release candidate (minimum: 1,000).
ONT_RECORD_LIMIT=10000 scripts/validate_modification_real_data.sh ont
```

The m6Anet inputs and GENCODE 39 annotation are SHA-256 pinned. Ensembl release
91 is checked against its published BSD checksum. The ONT reference/truth files
are SHA-256 pinned; the source BAM size and multipart ETag are recorded in the
generated report.

The default m6Anet synthetic catalog includes every transcript in both
`data.indiv_proba` and `data.info`. It proves schema, exact transcript-version
matching, retention auditing, and projection arithmetic, but its chromosomes
are not biological.

The fixture has 52 retained transcripts in `data.indiv_proba` and 82 across the
complete `data.info` candidate universe. The real-annotation audit found that
it is not tied to a single usable annotation: Ensembl 91 has 0/82 exact IDs, 76
version mismatches, and 6 absent stable IDs; GENCODE 39 has 17/82 exact IDs, 58
version mismatches, and 7 absent stable IDs. The importer therefore rejects
biological projection instead of silently stripping transcript versions. Use a
naturally generated m6Anet run with its exact reference annotation for
biological validation.

The full Ensembl 91 negative-projection scan also served as a loader stress
test on the same development machine. Before transcript prefiltering it took
168.61 s and reached 6,555,877,376 bytes maximum RSS; after prefiltering it took
47.16 s and reached 13,189,120 bytes. These are local engineering measurements,
not cross-platform performance guarantees.

At `ONT_RECORD_LIMIT=1000`, the current smoke test covers all 256 control and
all 256 m6A truth sites. It requires hard-call accuracy of at least 0.90 for the
unmodified control, at least 0.80 for m6A, and a m6A-minus-control modified
fraction of at least 0.70. These are regression gates for this importer, not a
claim about model performance. The public BAM header does not record an
explicit modified-base emission threshold, so the script records Dorado's
0.05 default. This validates decoding but does not satisfy strict
`verified_from_pg` provenance. For publication-grade benchmarking, regenerate
BAMs with the release recipe's explicit `--modified-bases-threshold 0` and an
auditable model command in `@PG`.

Sources:

- [m6Anet repository and pinned fixture](https://github.com/GoekeLab/m6anet/tree/590ec277cb48d61774f0872395099e466022e810/m6anet/tests/data)
- [m6Anet documentation](https://m6anet.readthedocs.io/)
- [SG-NEx data and Ensembl release-91 reference instructions](https://github.com/GoekeLab/sg-nex-data)
- [GENCODE 39 annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/)
- [ONT RNA all-5-mer ground-truth release](https://epi2me.nanoporetech.com/rna-mod-validation-all5mer-2026.07/)
- [SAM MM/ML/MN specification](https://samtools.github.io/hts-specs/SAMtags.pdf)
- [Dorado modified-base encoding](https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/)

## Scientific limitations to preserve in interpretation

- The Dorado importer supports all-context candidates and centered odd-length
  IUPAC motifs such as `DRACH`. Motifs are reconstructed from the as-sequenced
  read sequence, so basecalling errors in the motif flank can change whether a
  read contributes a candidate. An explicit call outside the declared motif is
  rejected rather than contaminating the denominator.
- A motif flank that crosses an exon boundary is reported as
  `context_dependent`. The current aggregator does not reconstruct a spliced
  transcript motif from FASTA, so genuine splice-junction candidates remain
  ineligible.
- The m6Anet annotation loader now prefilters a genome-wide GTF/GFF to the
  source transcript set. `mod-import-dorado` and `mod-aggregate` still
  materialize observations and sample×isoform×site rows in memory. Validate
  bounded locus/read subsets before attempting whole-genome all-context inputs.
- The effect-only contrast is suitable for QC and effect recovery, not
  biological inference. A later inferential layer should consume the integer
  design table and model biological samples, coverage/overdispersion, and
  multiple testing explicitly.

## Recommended next biological panel

Build one bounded HCT116 or A549 panel rather than downloading a whole study:

1. choose 5–10 genes with at least two well-covered splice isoforms and m6Anet
   DRACH sites;
2. pin the exact reads and the same versioned annotation used to generate the
   m6Anet source, plus reference FASTA, model, caller version, and thresholds;
3. retain at least three genuine biological samples per condition when testing
   condition effects;
4. include high, medium, and near-threshold sites, structurally absent sites,
   splice-boundary contexts, reverse-strand loci, unmatched reads, and one
   deliberately incomplete sample;
5. record wall time, peak RSS, join rates, callable coverage, recovered effect
   direction, and all exclusion reasons in a versioned report.

This panel should complement, not replace, the synthetic ground truth: natural
RNA supplies realistic isoforms, while the ONT constructs supply known
modified/unmodified labels.
