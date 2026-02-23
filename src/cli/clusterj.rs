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
    #[arg(long = "sw-score", default_value_t = 11, allow_hyphen_values = true)]
    pub sw_score: i64,

    /// Batch size to bound O(n^2) merge cost for very large genes (TrackCluster Python default: 500)
    #[arg(long = "batch-size", default_value_t = 500)]
    pub batch_size: usize,

    /// Maximum number of batching rounds before a final merge (TrackCluster Python default: 100)
    #[arg(long = "batch-rounds", default_value_t = 100)]
    pub batch_rounds: usize,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Full)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::cluster::clusterj::clusterj_with_name2_mode(
        &reads,
        Some(&refs),
        args.threads,
        args.sw_score,
        args.batch_size,
        args.batch_rounds,
        args.name2_mode,
    );

    crate::cluster::output::write_isoforms_bed(&args.out, &result.isoforms)?;

    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    crate::cluster::output::write_read_to_isoform_tsv(mapping_path, &result.read_to_isoform)?;

    let unused_path = args.out.with_extension("unused.bed");
    crate::cluster::output::write_isoforms_bed(unused_path, &result.unused)?;

    Ok(())
}
