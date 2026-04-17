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

    /// Batch size for iterative overlap clustering on large loci (0 disables batching)
    #[arg(long = "batch-size", default_value_t = 0)]
    pub batch_size: usize,

    /// Maximum number of intermediate batching rounds before the final full overlap merge
    #[arg(long = "batch-rounds", default_value_t = 100)]
    pub batch_rounds: usize,

    /// Smith-Waterman score cutoff for collapsing 5' truncations; set to -1 to disable collapsing
    #[arg(long = "sw-score", default_value_t = crate::cluster::cluster_overlap::DEFAULT_SW_SCORE, allow_hyphen_values = true)]
    pub sw_score: i64,

    /// Overlap-mode pass 1 cutoff
    #[arg(long = "cutoff1", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF1)]
    pub cutoff1: f64,

    /// Overlap-mode pass 2 cutoff
    #[arg(long = "cutoff2", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF2)]
    pub cutoff2: f64,

    /// Overlap-mode intron weighting
    #[arg(long = "intron-weight", default_value_t = crate::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT)]
    pub intron_weight: f64,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Coverage)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::cluster::cluster_overlap::cluster_with_options(
        &reads,
        Some(&refs),
        args.threads,
        crate::cluster::cluster_overlap::ClusterOptions {
            cutoff1: args.cutoff1,
            cutoff2: args.cutoff2,
            intron_weight: args.intron_weight,
            sw_score: args.sw_score,
            name2_mode: args.name2_mode,
            batch_size: args.batch_size,
            batch_rounds: args.batch_rounds,
        },
    );

    crate::cluster::output::write_isoforms_bed(&args.out, &result.isoforms)?;

    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    crate::cluster::output::write_read_to_isoform_tsv(mapping_path, &result.read_to_isoform)?;

    let unused_path = args.out.with_extension("unused.bed");
    crate::cluster::output::write_isoforms_bed(unused_path, &result.unused)?;

    Ok(())
}
