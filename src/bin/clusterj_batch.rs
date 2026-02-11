use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "clusterj-batch",
    about = "Batch runner for per-gene TrackCluster junction clustering"
)]
struct Args {
    /// Optional: prepare per-gene folders from a single reads BED before running clustering.
    /// Requires also setting --prepare-reference and --prepare-prefix.
    #[arg(long = "prepare-reads")]
    prepare_reads: Option<PathBuf>,

    /// Optional: reference BED for prepare step (used with --prepare-reads).
    #[arg(long = "prepare-reference")]
    prepare_reference: Option<PathBuf>,

    /// Optional: prefix for prepare step outputs (e.g. 488_aba_1). Required when using --prepare-reads.
    #[arg(long = "prepare-prefix")]
    prepare_prefix: Option<String>,

    /// Prepare step: minimum fraction of read span overlapping a reference span
    #[arg(long = "prepare-fraction-read", default_value_t = 0.01)]
    prepare_fraction_read: f64,

    /// Prepare step: minimum fraction of reference span overlapping a read span
    #[arg(long = "prepare-fraction-ref", default_value_t = 0.05)]
    prepare_fraction_ref: f64,

    /// Directory containing per-gene folders (e.g. `/.../tracktest/<gene>/`)
    #[arg(long = "input-root")]
    input_root: PathBuf,

    /// Optional file containing gene folder names (one per line)
    #[arg(long = "gene-list")]
    gene_list: Option<PathBuf>,

    /// Output directory to write per-gene results
    #[arg(long = "output-root")]
    output_root: PathBuf,

    /// Number of worker threads (parallel across genes)
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Smith-Waterman score cutoff for collapsing 5' truncations; set to -1 to disable collapsing
    #[arg(long = "sw-score", default_value_t = 11, allow_hyphen_values = true)]
    sw_score: i64,

    /// Batch size to bound O(n^2) merge cost for very large genes (TrackCluster Python default: 500)
    #[arg(long = "batch-size", default_value_t = 500)]
    batch_size: usize,

    /// Maximum number of batching rounds before a final merge (TrackCluster Python default: 100)
    #[arg(long = "batch-rounds", default_value_t = 100)]
    batch_rounds: usize,

    /// Overwrite existing outputs (default: skip genes whose output file already exists)
    #[arg(long = "force")]
    force: bool,

    /// Print a progress line every N genes
    #[arg(long = "progress-every", default_value_t = 1000)]
    progress_every: usize,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    trackcluster_rs::flow::full::run_clusterj_batch(
        trackcluster_rs::flow::full::BatchRunOptions {
            prepare_reads: args.prepare_reads,
            prepare_reference: args.prepare_reference,
            prepare_prefix: args.prepare_prefix,
            prepare_fraction_read: args.prepare_fraction_read,
            prepare_fraction_ref: args.prepare_fraction_ref,
            input_root: args.input_root,
            gene_list: args.gene_list,
            output_root: args.output_root,
            threads: args.threads,
            sw_score: args.sw_score,
            batch_size: args.batch_size,
            batch_rounds: args.batch_rounds,
            force: args.force,
            progress_every: args.progress_every,
        },
    )?;
    Ok(())
}
