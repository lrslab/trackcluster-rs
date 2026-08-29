use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED (sorted recommended)
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Output isoform BED
    #[arg(short = 'o', long = "out", default_value = "isoform.bed")]
    pub out: PathBuf,

    /// Number of threads
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

    /// Handling of malformed read tracks: skip only that track, or fail the command
    #[arg(long = "invalid-read-policy", default_value_t = crate::flow::config::InvalidReadPolicy::Skip)]
    pub invalid_read_policy: crate::flow::config::InvalidReadPolicy,

    /// Per-locus cap: reservoir-sample reads down to this count (set to 0 to disable).
    /// Dropped reads are written to unused.bed and are not scaled back into counts.
    #[arg(
        long = "max-reads-per-locus",
        default_value_t = crate::flow::config::DEFAULT_MAX_READS_PER_GENE
    )]
    pub max_reads_per_locus: usize,

    /// Deterministic RNG seed mixed with chrom, strand, and locus span for per-locus sampling.
    #[arg(long = "downsample-seed", default_value_t = 1)]
    pub downsample_seed: u64,

    /// Emit a heartbeat status line every N seconds during chrom/strand partitioning (0 disables).
    #[arg(long = "heartbeat-seconds", default_value_t = 60)]
    pub heartbeat_seconds: u64,

    /// When a heartbeat sees no progress, print up to this many in-flight partitions (0 => 1).
    #[arg(long = "heartbeat-top", default_value_t = 5)]
    pub heartbeat_top: usize,
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

    #[cfg(test)]
    pub fn resolved_platform_options(&self) -> crate::cluster::clusterj::ResolvedPlatformOptions {
        self.junction_config().resolved()
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    let unused_path = args.out.with_extension("unused.bed");
    let rejected_path = args.out.with_extension("rejected_reads.tsv");
    super::ensure_distinct_inputs_and_outputs(
        &[
            ("reads input", args.reads.as_path()),
            ("reference input", args.reference.as_path()),
        ],
        &[
            ("isoform BED output", args.out.as_path()),
            ("read-to-isoform output", mapping_path.as_path()),
            ("unused-read BED output", unused_path.as_path()),
            ("rejected-read output", rejected_path.as_path()),
        ],
    )?;
    let (reads, rejected_reads) = super::read_read_tracks(&args.reads, args.invalid_read_policy)?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    if args.max_reads_per_locus == 0 {
        eprintln!(
            "clusterj: note: --max-reads-per-locus=0 disables the per-locus read cap; \
large overlapping loci can take a long time"
        );
    }
    let junction = args.junction_config();
    let (result, summary) =
        crate::cluster::clusterj::try_clusterj_with_runtime_options_and_summary(
            &reads,
            Some(&refs),
            args.threads,
            args.sw_score,
            args.batch_size,
            args.batch_rounds,
            args.name2_mode,
            junction.sl,
            junction.three_prime,
            junction.correction,
            crate::cluster::clusterj::ClusterjRuntimeOptions {
                max_reads_per_locus: args.max_reads_per_locus,
                downsample_seed: args.downsample_seed,
                heartbeat_seconds: args.heartbeat_seconds,
                heartbeat_top: args.heartbeat_top,
            },
        )?;
    summary.emit();

    crate::flow::artifact_manifest::atomic_write_with(&args.out, |temporary| {
        crate::cluster::output::write_isoforms_bed_to_writer(temporary, &result.isoforms)
            .map_err(Into::into)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&mapping_path, |temporary| {
        crate::cluster::output::write_read_to_isoform_tsv_writer(temporary, &result.read_to_isoform)
            .map_err(Into::into)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&unused_path, |temporary| {
        crate::cluster::output::write_isoforms_bed_to_writer(temporary, &result.unused)
            .map_err(Into::into)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&rejected_path, |temporary| {
        crate::io::bed::write_rejected_reads_tsv_to_writer(temporary, &rejected_reads)
            .map_err(Into::into)
    })?;
    if !rejected_reads.is_empty() {
        eprintln!(
            "clusterj: warning: excluded {} malformed read track(s); details: {:?}",
            rejected_reads.len(),
            rejected_path
        );
    }

    Ok(())
}
