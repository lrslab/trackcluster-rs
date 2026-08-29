use std::ffi::OsString;
use std::path::{Path, PathBuf};

use anyhow::Context;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Clustering algorithm used inside flow: `clusterj` (junction mode) or `cluster` (overlap mode)
    #[arg(long = "cluster-mode", default_value_t = crate::flow::full::ClusterMode::Clusterj)]
    pub cluster_mode: crate::flow::full::ClusterMode,

    /// Reads BED (single-sample mode; mutually exclusive with --manifest)
    #[arg(short = 's', long = "reads")]
    pub reads: Option<PathBuf>,

    /// Multi-sample manifest TSV with columns: sample, reads, optional group
    #[arg(long = "manifest")]
    pub manifest: Option<PathBuf>,

    /// Reference BED
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Output directory (created if missing)
    #[arg(short = 'o', long = "output-root")]
    pub output_root: PathBuf,

    /// Prefix for merged outputs like `<prefix>_isoform.bed`
    #[arg(long = "prefix")]
    pub prefix: String,

    /// Number of worker threads (parallel across genes)
    #[arg(short = 't', long = "threads", default_value_t = 1, allow_hyphen_values = true, value_parser = crate::config::parse_worker_threads)]
    pub threads: usize,

    /// Smith-Waterman score cutoff for SL-supported 5' protection; default -1 treats reads as having no SW 5' signal
    #[arg(long = "sw-score", default_value_t = crate::cluster::clusterj::DEFAULT_SW_SCORE, allow_hyphen_values = true)]
    pub sw_score: i64,

    /// Batch size to bound O(n^2) merge cost for very large genes (TrackCluster Python default: 500)
    #[arg(long = "batch-size", default_value_t = 500)]
    pub batch_size: usize,

    /// Maximum number of batching rounds before a final merge (TrackCluster Python default: 100)
    #[arg(long = "batch-rounds", default_value_t = 100, value_parser = crate::config::parse_batch_rounds)]
    pub batch_rounds: usize,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Coverage)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,

    /// Platform preset used to seed junction correction, SL 5', and 3' defaults: generic, rna002, or rna004
    #[arg(long = "platform-preset", default_value_t = crate::cluster::clusterj::PlatformPreset::Generic)]
    pub platform_preset: crate::cluster::clusterj::PlatformPreset,

    /// Internal junction correction offset in bp; distinct from SL/5' and 3' terminal offsets
    #[arg(long = "junction-correction-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub junction_correction_offset: Option<u32>,

    /// Minimum weighted support for a junction site to avoid correction/filtering
    #[arg(long = "junction-correction-min-support", allow_hyphen_values = true, value_parser = crate::config::parse_weighted_minimum_support)]
    pub junction_correction_min_support: Option<u32>,

    /// SL-supported partial/5' truncation biological 5' offset tolerated for merging
    #[arg(long = "sl-partial-5prime-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub sl_partial_5prime_offset: Option<u32>,

    /// SL-supported same-junction biological 5' offset tolerated for merging
    #[arg(long = "sl-same-junction-5prime-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub sl_same_junction_5prime_offset: Option<u32>,

    /// Offset used to group SL-supported reads into the same biological 5' cluster
    #[arg(long = "sl-5prime-cluster-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub sl_5prime_cluster_offset: Option<u32>,

    /// Minimum read support required for an SL 5' cluster to protect a candidate isoform
    #[arg(long = "sl-5prime-min-support", allow_hyphen_values = true, value_parser = crate::config::parse_minimum_support)]
    pub sl_5prime_min_support: Option<usize>,

    /// Same-junction biological 3' offset tolerated for merging
    #[arg(long = "same-junction-3prime-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub same_junction_3prime_offset: Option<u32>,

    /// Offset used to group reads into the same biological 3' cluster; defaults to the active junction correction offset
    #[arg(long = "3prime-cluster-offset", allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub three_prime_cluster_offset: Option<u32>,

    /// Minimum read support required for a 3' cluster to protect a candidate isoform
    #[arg(long = "3prime-min-support", allow_hyphen_values = true, value_parser = crate::config::parse_minimum_support)]
    pub three_prime_min_support: Option<usize>,

    /// Overlap-mode pass 1 cutoff (`cluster-mode=cluster` only)
    #[arg(long = "overlap-cutoff1", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF1, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub overlap_cutoff1: f64,

    /// Overlap-mode pass 2 cutoff (`cluster-mode=cluster` only)
    #[arg(long = "overlap-cutoff2", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF2, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub overlap_cutoff2: f64,

    /// Overlap-mode intron weighting (`cluster-mode=cluster` only)
    #[arg(long = "overlap-intron-weight", default_value_t = crate::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT, allow_hyphen_values = true, value_parser = crate::config::parse_nonnegative_weight)]
    pub overlap_intron_weight: f64,

    /// Prepare step: minimum fraction of read span overlapping a reference span
    #[arg(long = "prepare-fraction-read", default_value_t = 0.01, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub prepare_fraction_read: f64,

    /// Prepare step: minimum fraction of reference span overlapping a read span
    #[arg(long = "prepare-fraction-ref", default_value_t = 0.05, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub prepare_fraction_ref: f64,

    /// Counting mode for reads with multiple isoform candidates: fractional or unique
    #[arg(long = "assignment-mode", default_value_t = crate::count::AssignmentMode::Unique)]
    pub assignment_mode: crate::count::AssignmentMode,

    /// Junction tolerance in bp used only for unique read-to-isoform assignment
    #[arg(
        long = "unique-assignment-junction-offset",
        default_value_t = crate::count::DEFAULT_UNIQUE_ASSIGNMENT_JUNCTION_OFFSET
    )]
    pub unique_assignment_junction_offset: u32,

    /// Optional normalized modification manifest processed after unique assignment.
    #[arg(long = "mod-manifest")]
    pub mod_manifest: Option<PathBuf>,

    /// Indexed genomic FASTA used to reject modification canonical-base mismatches.
    #[arg(long = "mod-reference-fasta", requires = "mod_manifest")]
    pub mod_reference_fasta: Option<PathBuf>,

    /// Per-assay modification hard-call threshold as ASSAY_ID=PROBABILITY.
    #[arg(
        long = "mod-analysis-threshold",
        requires = "mod_manifest",
        value_parser = crate::cli::mod_aggregate::parse_analysis_threshold
    )]
    mod_analysis_thresholds: Vec<crate::cli::mod_aggregate::AnalysisThreshold>,

    /// Modification eligibility policy: exploratory or strict.
    #[arg(
        long = "mod-eligibility-profile",
        requires = "mod_manifest",
        default_value_t = crate::modification::EligibilityProfile::Exploratory
    )]
    pub mod_eligibility_profile: crate::modification::EligibilityProfile,

    /// Minimum independently covering molecules for strict modification eligibility.
    #[arg(
        long = "mod-min-covering",
        requires = "mod_manifest",
        default_value_t = 20,
        value_parser = clap::value_parser!(u64).range(1..)
    )]
    pub mod_min_covering: u64,

    /// Minimum callable molecules for an eligible modification site row.
    #[arg(
        long = "mod-min-callable",
        requires = "mod_manifest",
        value_parser = clap::value_parser!(u64).range(1..)
    )]
    pub mod_min_callable: Option<u64>,

    /// Minimum candidate/covering molecule fraction for strict eligibility.
    #[arg(
        long = "mod-min-candidate-rate",
        requires = "mod_manifest",
        default_value = "0.8",
        value_parser = crate::cli::mod_aggregate::parse_rate
    )]
    pub mod_min_candidate_rate: f64,

    /// Minimum callable/covering molecule fraction for strict eligibility.
    #[arg(
        long = "mod-min-callable-rate",
        requires = "mod_manifest",
        default_value = "0.8",
        value_parser = crate::cli::mod_aggregate::parse_rate
    )]
    pub mod_min_callable_rate: f64,

    /// Minimum assigned-read join rate for each sample/assay.
    #[arg(
        long = "mod-min-read-join-rate",
        requires = "mod_manifest",
        default_value = "0.9",
        value_parser = crate::cli::mod_aggregate::parse_rate
    )]
    pub mod_min_read_join_rate: f64,

    /// Emit low-join modification rows as ineligible instead of failing.
    #[arg(
        long = "mod-allow-low-global-join",
        visible_alias = "mod-allow-low-join",
        requires = "mod_manifest"
    )]
    pub mod_allow_low_join: bool,

    /// Optional explicit contrast specification run after modification aggregation.
    #[arg(long = "mod-contrasts", requires = "mod_manifest")]
    pub mod_contrasts: Option<PathBuf>,

    /// Manifest mode: also write `<prefix>_pooled_reads.bed` for debugging/compatibility
    #[arg(long = "emit-pooled-reads")]
    pub emit_pooled_reads: bool,

    /// Rebuild every gene instead of reusing an exact, hash-verified completion manifest
    #[arg(long = "force")]
    pub force: bool,

    /// Stop before merge/count/description if any gene fails; by default failed genes are logged and excluded
    #[arg(long = "strict-gene-errors")]
    pub strict_gene_errors: bool,

    /// Handling of malformed read tracks: skip only that track, or fail the enclosing stage
    #[arg(long = "invalid-read-policy", default_value_t = crate::flow::config::InvalidReadPolicy::Skip)]
    pub invalid_read_policy: crate::flow::config::InvalidReadPolicy,

    /// Reuse exact, hash-verified per-gene completion manifests and run only merge/count/desc outputs
    #[arg(long = "count-only")]
    pub count_only: bool,

    /// Print a progress line every N genes
    #[arg(long = "progress-every", default_value_t = 1000)]
    pub progress_every: usize,

    /// Emit a heartbeat status line every N seconds during per-gene clustering (0 disables).
    /// Useful when large genes make progress appear "stuck" because no gene completes for a while.
    #[arg(long = "heartbeat-seconds", default_value_t = 60)]
    pub heartbeat_seconds: u64,

    /// When a heartbeat sees no progress, print up to this many in-flight genes (0 => 1).
    #[arg(long = "heartbeat-top", default_value_t = 5)]
    pub heartbeat_top: usize,

    /// Restrict downsampling to these biological gene ID(s) only (repeatable).
    /// If omitted and `--max-reads-per-gene > 0`, downsampling applies to all genes.
    #[arg(long = "downsample-gene")]
    pub downsample_genes: Vec<String>,

    /// Per-gene cap: reservoir-sample reads down to this count (set to 0 to disable downsampling).
    #[arg(
        long = "max-reads-per-gene",
        default_value_t = crate::flow::config::DEFAULT_MAX_READS_PER_GENE
    )]
    pub max_reads_per_gene: usize,

    /// Deterministic RNG seed used for per-gene downsampling.
    #[arg(long = "downsample-seed", default_value_t = 1)]
    pub downsample_seed: u64,
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

const MOD_AGGREGATE_OUTPUT_SUFFIXES: [&str; 4] = [
    ".mod_join_qc.tsv",
    ".mod_site_join_qc.tsv",
    ".isoform_mod_sites.tsv",
    ".isoform_mod_design.tsv",
];
const MOD_CONTRAST_OUTPUT_SUFFIX: &str = ".isoform_mod_contrasts.tsv";

fn remove_stale_optional_output(prefix: &Path, suffix: &str) -> anyhow::Result<()> {
    let path = append_suffix(prefix, suffix);
    let metadata = match std::fs::symlink_metadata(&path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(error) => {
            return Err(error).with_context(|| format!("inspect stale flow output {path:?}"));
        }
    };
    if !metadata.file_type().is_file() {
        anyhow::bail!("refusing to remove stale flow output {path:?}: expected a regular file");
    }
    std::fs::remove_file(&path).with_context(|| format!("remove stale flow output {path:?}"))
}

fn remove_stale_modification_outputs(prefix: &Path) -> anyhow::Result<()> {
    for suffix in MOD_AGGREGATE_OUTPUT_SUFFIXES {
        remove_stale_optional_output(prefix, suffix)?;
    }
    remove_stale_optional_output(prefix, MOD_CONTRAST_OUTPUT_SUFFIX)
}

impl Args {
    pub fn junction_config(&self) -> crate::flow::config::JunctionConfig {
        crate::flow::config::JunctionConfig::resolve(
            self.platform_preset,
            crate::flow::config::JunctionOverrides {
                correction_offset: self.junction_correction_offset,
                correction_min_support: self.junction_correction_min_support,
                sl_partial_five_prime_offset: self.sl_partial_5prime_offset,
                sl_same_junction_five_prime_offset: self.sl_same_junction_5prime_offset,
                sl_five_prime_cluster_offset: self.sl_5prime_cluster_offset,
                sl_five_prime_min_support: self.sl_5prime_min_support,
                same_junction_three_prime_offset: self.same_junction_3prime_offset,
                three_prime_cluster_offset: self.three_prime_cluster_offset,
                three_prime_min_support: self.three_prime_min_support,
            },
        )
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.reads.is_some() && args.manifest.is_some() {
        anyhow::bail!("flow: use either --reads or --manifest, not both");
    }
    if !args.count_only && args.reads.is_none() && args.manifest.is_none() {
        anyhow::bail!("flow: one of --reads or --manifest is required");
    }
    if args.count_only && args.strict_gene_errors {
        anyhow::bail!(
            "flow: --strict-gene-errors applies to per-gene execution and cannot be used with --count-only"
        );
    }
    if args.count_only && args.force {
        anyhow::bail!(
            "flow: --force rebuilds per-gene results and cannot be used with --count-only"
        );
    }
    if args.mod_manifest.is_some() {
        if args.manifest.is_none() {
            anyhow::bail!("flow: --mod-manifest requires multi-sample --manifest input");
        }
        if args.assignment_mode != crate::count::AssignmentMode::Unique {
            anyhow::bail!("flow: modification aggregation requires --assignment-mode unique");
        }
        if args.mod_analysis_thresholds.is_empty() {
            anyhow::bail!("flow: --mod-manifest requires one --mod-analysis-threshold per assay");
        }
    } else if !args.mod_analysis_thresholds.is_empty()
        || args.mod_contrasts.is_some()
        || args.mod_reference_fasta.is_some()
        || args.mod_min_callable.is_some()
        || args.mod_eligibility_profile != crate::modification::EligibilityProfile::Exploratory
        || args.mod_min_covering != 20
        || args.mod_min_candidate_rate != 0.8
        || args.mod_min_callable_rate != 0.8
        || args.mod_min_read_join_rate != 0.9
        || args.mod_allow_low_join
    {
        anyhow::bail!("flow: modification analysis options require --mod-manifest");
    }
    let junction = args.junction_config();
    let flow_manifest = args.manifest.clone();
    let mod_manifest = args.mod_manifest.clone();
    let mod_contrasts = args.mod_contrasts.clone();
    let mod_thresholds =
        crate::cli::mod_aggregate::collect_thresholds(args.mod_analysis_thresholds.clone())?;
    let mod_min_callable = args
        .mod_min_callable
        .unwrap_or_else(|| args.mod_eligibility_profile.default_min_callable());
    let mod_output_prefix = args.output_root.join(&args.prefix);
    if mod_manifest.is_some() {
        crate::modification::generation::ensure_managed(&mod_output_prefix)?;
    }
    crate::modification::generation::invalidate_current(&mod_output_prefix)?;

    let result = crate::flow::full::run_full_flow(crate::flow::full::FullFlowOptions {
        cluster_mode: args.cluster_mode,
        reads: args.reads,
        manifest: args.manifest,
        reference: args.reference,
        output_root: args.output_root,
        prefix: args.prefix,
        prepare: crate::flow::config::PrepareConfig {
            fraction_read: args.prepare_fraction_read,
            fraction_ref: args.prepare_fraction_ref,
        },
        clustering: crate::flow::config::ClusteringConfig {
            sw_score: args.sw_score,
            batch_size: args.batch_size,
            batch_rounds: args.batch_rounds,
            name2_mode: args.name2_mode,
            junction,
            overlap: crate::flow::config::OverlapConfig {
                cutoff1: args.overlap_cutoff1,
                cutoff2: args.overlap_cutoff2,
                intron_weight: args.overlap_intron_weight,
            },
        },
        counting: crate::flow::config::CountingConfig {
            assignment_mode: args.assignment_mode,
            unique_assignment: crate::count::UniqueAssignmentOptions {
                junction_offset: args.unique_assignment_junction_offset,
            },
        },
        runtime: crate::flow::config::RuntimeConfig {
            threads: args.threads,
            force: args.force,
            progress_every: args.progress_every,
            heartbeat_seconds: args.heartbeat_seconds,
            heartbeat_top: args.heartbeat_top,
            gene_error_policy: if args.strict_gene_errors {
                crate::flow::config::GeneErrorPolicy::Strict
            } else {
                crate::flow::config::GeneErrorPolicy::Continue
            },
            invalid_read_policy: args.invalid_read_policy,
        },
        downsample: crate::flow::config::DownsampleConfig {
            genes: args.downsample_genes,
            max_reads_per_gene: args.max_reads_per_gene,
            seed: args.downsample_seed,
        },
        emit_pooled_reads: args.emit_pooled_reads,
        count_only: args.count_only,
    })?;

    if let Some(mod_manifest) = mod_manifest.as_deref() {
        let sample_manifest = flow_manifest
            .as_deref()
            .expect("validated modification flow manifest");
        let read_to_isoform = result
            .unique_read_to_isoform_tsv
            .as_deref()
            .expect("validated unique assignment mode");
        let generation = crate::modification::generation::begin(&mod_output_prefix).with_context(
            || {
                "flow modification generation failed; no current pointer was published and flat modification outputs must not be used"
            },
        )?;
        let generation_outputs = crate::cli::mod_aggregate::AggregateOutputPaths {
            join_qc: generation.join_qc.clone(),
            site_join_qc: generation.site_join_qc.clone(),
            sites: generation.sites.clone(),
            design: generation.design.clone(),
        };
        let generation_options = crate::modification::generation::GenerationOptions::new(
            mod_thresholds.clone(),
            args.mod_eligibility_profile,
            args.mod_min_covering,
            mod_min_callable,
            args.mod_min_candidate_rate,
            args.mod_min_callable_rate,
            args.mod_min_read_join_rate,
            args.mod_allow_low_join,
            args.mod_reference_fasta.is_some(),
            mod_contrasts.is_some(),
        );
        let inputs_before = crate::modification::generation::collect_inputs(
            sample_manifest,
            &result.isoform_bed,
            read_to_isoform,
            mod_manifest,
            args.mod_reference_fasta.as_deref(),
            mod_contrasts.as_deref(),
        )
        .with_context(|| {
            "flow modification generation failed; no current pointer was published and flat modification outputs must not be used"
        })?;
        let generation_result = (|| -> anyhow::Result<_> {
            let mod_result = crate::cli::mod_aggregate::run_with_output_paths(
                sample_manifest,
                &result.isoform_bed,
                read_to_isoform,
                mod_manifest,
                args.mod_reference_fasta.as_deref(),
                mod_thresholds,
                args.mod_eligibility_profile,
                args.mod_min_covering,
                mod_min_callable,
                args.mod_min_candidate_rate,
                args.mod_min_callable_rate,
                args.mod_min_read_join_rate,
                args.mod_allow_low_join,
                &generation_outputs,
            )?;
            if let Some(contrasts) = mod_contrasts.as_deref() {
                crate::cli::mod_contrast::run_with_output_path(
                    &generation.design,
                    contrasts,
                    &generation.contrasts,
                )?;
            }
            let inputs_after = crate::modification::generation::collect_inputs(
                sample_manifest,
                &result.isoform_bed,
                read_to_isoform,
                mod_manifest,
                args.mod_reference_fasta.as_deref(),
                mod_contrasts.as_deref(),
            )?;
            if inputs_before != inputs_after {
                anyhow::bail!("modification inputs changed while the generation was running");
            }
            crate::modification::generation::publish(
                &mod_output_prefix,
                &generation,
                inputs_after,
                generation_options,
                mod_contrasts.is_some(),
            )?;
            Ok(mod_result)
        })()
        .with_context(|| {
            "flow modification generation failed; no current pointer was published and flat modification outputs must not be used"
        })?;
        eprintln!(
            "flow: modification join_qc_rows={} site_rows={} design_rows={}",
            generation_result.join_qc.len(),
            generation_result.sites.len(),
            generation_result.design.len()
        );
    } else {
        remove_stale_modification_outputs(&mod_output_prefix)?;
    }

    eprintln!(
        "flow: isoform={:?} unused={:?} count={:?} desc_prefix={:?}",
        result.isoform_bed, result.unused_bed, result.count_csv, result.desc_prefix
    );
    eprintln!(
        "flow: batch total={} processed={} skipped={} errors={} elapsed_s={}",
        result.batch.total_genes,
        result.batch.processed,
        result.batch.skipped,
        result.batch.errors,
        result.batch.elapsed_seconds
    );
    if let Some(multi) = result.multi_sample {
        eprintln!(
            "flow: multi-sample count={:?} long={:?} matrix={:?} group={:?}",
            multi.count_csv, multi.long_tsv, multi.matrix_tsv, multi.group_tsv
        );
    }

    Ok(())
}
