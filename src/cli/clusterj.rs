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

    /// Smith-Waterman score cutoff for SL-supported 5' protection; set to -1 to treat reads as having no SW 5' signal
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

    /// Platform preset used to seed junction correction, SL 5', and 3' defaults: generic, rna002, or rna004
    #[arg(long = "platform-preset", default_value_t = crate::cluster::clusterj::PlatformPreset::Generic)]
    pub platform_preset: crate::cluster::clusterj::PlatformPreset,

    /// Internal junction correction offset in bp; distinct from SL/5' and 3' terminal offsets
    #[arg(long = "junction-correction-offset")]
    pub junction_correction_offset: Option<u32>,

    /// Minimum weighted support for a junction site to avoid correction/filtering
    #[arg(long = "junction-correction-min-support")]
    pub junction_correction_min_support: Option<u32>,

    /// SL-supported partial/5' truncation biological 5' offset tolerated for merging
    #[arg(long = "sl-partial-5prime-offset")]
    pub sl_partial_5prime_offset: Option<u32>,

    /// SL-supported same-junction biological 5' offset tolerated for merging
    #[arg(long = "sl-same-junction-5prime-offset")]
    pub sl_same_junction_5prime_offset: Option<u32>,

    /// Offset used to group SL-supported reads into the same biological 5' cluster
    #[arg(long = "sl-5prime-cluster-offset")]
    pub sl_5prime_cluster_offset: Option<u32>,

    /// Minimum read support required for an SL 5' cluster to protect a candidate isoform
    #[arg(long = "sl-5prime-min-support")]
    pub sl_5prime_min_support: Option<usize>,

    /// Same-junction biological 3' offset tolerated for merging
    #[arg(long = "same-junction-3prime-offset")]
    pub same_junction_3prime_offset: Option<u32>,

    /// Offset used to group reads into the same biological 3' cluster; defaults to the active junction correction offset
    #[arg(long = "3prime-cluster-offset")]
    pub three_prime_cluster_offset: Option<u32>,

    /// Minimum read support required for a 3' cluster to protect a candidate isoform
    #[arg(long = "3prime-min-support")]
    pub three_prime_min_support: Option<usize>,
}

impl Args {
    pub fn resolved_platform_options(&self) -> crate::cluster::clusterj::ResolvedPlatformOptions {
        crate::cluster::clusterj::resolve_platform_options(
            self.platform_preset,
            self.junction_correction_offset,
            self.junction_correction_min_support,
            self.sl_partial_5prime_offset,
            self.sl_same_junction_5prime_offset,
            self.sl_5prime_cluster_offset,
            self.sl_5prime_min_support,
            self.same_junction_3prime_offset,
            self.three_prime_cluster_offset,
            self.three_prime_min_support,
        )
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let resolved_options = args.resolved_platform_options();

    let result = crate::cluster::clusterj::clusterj_with_options(
        &reads,
        Some(&refs),
        args.threads,
        args.sw_score,
        args.batch_size,
        args.batch_rounds,
        args.name2_mode,
        resolved_options.sl_options,
        resolved_options.three_prime_options,
        resolved_options.junction_correction,
    );

    crate::cluster::output::write_isoforms_bed(&args.out, &result.isoforms)?;

    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    crate::cluster::output::write_read_to_isoform_tsv(mapping_path, &result.read_to_isoform)?;

    let unused_path = args.out.with_extension("unused.bed");
    crate::cluster::output::write_isoforms_bed(unused_path, &result.unused)?;

    Ok(())
}
