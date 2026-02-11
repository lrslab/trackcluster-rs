use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

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

    /// Prepare step: minimum fraction of read span overlapping a reference span
    #[arg(long = "prepare-fraction-read", default_value_t = 0.01)]
    pub prepare_fraction_read: f64,

    /// Prepare step: minimum fraction of reference span overlapping a read span
    #[arg(long = "prepare-fraction-ref", default_value_t = 0.05)]
    pub prepare_fraction_ref: f64,

    /// Overwrite existing outputs (default: skip genes whose output file already exists)
    #[arg(long = "force")]
    pub force: bool,

    /// Print a progress line every N genes
    #[arg(long = "progress-every", default_value_t = 1000)]
    pub progress_every: usize,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let result = crate::flow::full::run_full_flow(crate::flow::full::FullFlowOptions {
        reads: args.reads,
        reference: args.reference,
        output_root: args.output_root,
        prefix: args.prefix,
        threads: args.threads,
        sw_score: args.sw_score,
        batch_size: args.batch_size,
        batch_rounds: args.batch_rounds,
        prepare_fraction_read: args.prepare_fraction_read,
        prepare_fraction_ref: args.prepare_fraction_ref,
        force: args.force,
        progress_every: args.progress_every,
    })?;

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

    Ok(())
}
