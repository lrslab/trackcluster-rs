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
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,

    /// Smith-Waterman score cutoff for collapsing 5' truncations; set to -1 to disable collapsing
    #[arg(long = "sw-score", default_value_t = crate::cluster::clusterj::DEFAULT_SW_SCORE, allow_hyphen_values = true)]
    pub sw_score: i64,

    /// Batch size to bound O(n^2) merge cost for very large genes (TrackCluster Python default: 500)
    #[arg(long = "batch-size", default_value_t = 500)]
    pub batch_size: usize,

    /// Maximum number of batching rounds before a final merge (TrackCluster Python default: 100)
    #[arg(long = "batch-rounds", default_value_t = 100)]
    pub batch_rounds: usize,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Coverage)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,

    /// SL-supported partial/5' truncation biological 5' offset tolerated for merging
    #[arg(long = "sl-partial-5prime-offset", default_value_t = crate::cluster::clusterj::DEFAULT_SL_PARTIAL_FIVE_PRIME_END_OFFSET)]
    pub sl_partial_5prime_offset: u32,

    /// SL-supported same-junction biological 5' offset tolerated for merging
    #[arg(long = "sl-same-junction-5prime-offset", default_value_t = crate::cluster::clusterj::DEFAULT_SL_SAME_JUNCTION_FIVE_PRIME_END_OFFSET)]
    pub sl_same_junction_5prime_offset: u32,

    /// Offset used to group SL-supported reads into the same biological 5' cluster
    #[arg(long = "sl-5prime-cluster-offset", default_value_t = crate::cluster::clusterj::DEFAULT_SL_FIVE_PRIME_CLUSTER_OFFSET)]
    pub sl_5prime_cluster_offset: u32,

    /// Minimum read support required for an SL 5' cluster to protect a candidate isoform
    #[arg(long = "sl-5prime-min-support", default_value_t = crate::cluster::clusterj::DEFAULT_MIN_SL_FIVE_PRIME_CLUSTER_SUPPORT)]
    pub sl_5prime_min_support: usize,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::cluster::clusterj::clusterj_with_options(
        &reads,
        Some(&refs),
        args.threads,
        args.sw_score,
        args.batch_size,
        args.batch_rounds,
        args.name2_mode,
        crate::cluster::clusterj::SlMergeOptions {
            partial_five_prime_end_offset: args.sl_partial_5prime_offset,
            same_junction_five_prime_end_offset: args.sl_same_junction_5prime_offset,
            five_prime_cluster_offset: args.sl_5prime_cluster_offset,
            min_five_prime_cluster_support: args.sl_5prime_min_support,
        },
    );

    crate::cluster::output::write_isoforms_bed(&args.out, &result.isoforms)?;

    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    crate::cluster::output::write_read_to_isoform_tsv(mapping_path, &result.read_to_isoform)?;

    let unused_path = args.out.with_extension("unused.bed");
    crate::cluster::output::write_isoforms_bed(unused_path, &result.unused)?;

    Ok(())
}
