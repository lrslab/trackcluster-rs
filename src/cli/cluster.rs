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

    /// Batch size for iterative overlap clustering on large loci (0 disables batching)
    #[arg(long = "batch-size", default_value_t = 0)]
    pub batch_size: usize,

    /// Maximum number of intermediate batching rounds before the final full overlap merge
    #[arg(long = "batch-rounds", default_value_t = 100, value_parser = crate::config::parse_batch_rounds)]
    pub batch_rounds: usize,

    /// Smith-Waterman score cutoff for SL-supported 5' protection; set to -1 to treat reads as having no SW 5' signal
    #[arg(long = "sw-score", default_value_t = crate::cluster::cluster_overlap::DEFAULT_SW_SCORE, allow_hyphen_values = true)]
    pub sw_score: i64,

    /// Overlap-mode pass 1 cutoff
    #[arg(long = "cutoff1", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF1, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub cutoff1: f64,

    /// Overlap-mode pass 2 cutoff
    #[arg(long = "cutoff2", default_value_t = crate::cluster::cluster_overlap::DEFAULT_CUTOFF2, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub cutoff2: f64,

    /// Overlap-mode intron weighting
    #[arg(long = "intron-weight", default_value_t = crate::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT, allow_hyphen_values = true, value_parser = crate::config::parse_nonnegative_weight)]
    pub intron_weight: f64,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Coverage)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,

    /// Handling of malformed read tracks: skip only that track, or fail the command
    #[arg(long = "invalid-read-policy", default_value_t = crate::flow::config::InvalidReadPolicy::Skip)]
    pub invalid_read_policy: crate::flow::config::InvalidReadPolicy,
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

    let result = crate::flow::config::ClusteringConfig {
        sw_score: args.sw_score,
        batch_size: args.batch_size,
        batch_rounds: args.batch_rounds,
        name2_mode: args.name2_mode,
        junction: crate::flow::config::JunctionConfig::default(),
        overlap: crate::flow::config::OverlapConfig {
            cutoff1: args.cutoff1,
            cutoff2: args.cutoff2,
            intron_weight: args.intron_weight,
        },
    }
    .cluster_gene(
        crate::flow::full::ClusterMode::Cluster,
        &reads,
        &refs,
        args.threads,
    )?;

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
            "cluster: warning: excluded {} malformed read track(s); details: {:?}",
            rejected_reads.len(),
            rejected_path
        );
    }

    Ok(())
}
