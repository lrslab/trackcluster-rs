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
    #[arg(long = "sw-score", default_value_t = trackcluster_rs::cluster::clusterj::DEFAULT_SW_SCORE, allow_hyphen_values = true)]
    sw_score: i64,

    /// Batch size to bound O(n^2) merge cost for very large genes (TrackCluster Python default: 500)
    #[arg(long = "batch-size", default_value_t = 500)]
    batch_size: usize,

    /// Maximum number of batching rounds before a final merge (TrackCluster Python default: 100)
    #[arg(long = "batch-rounds", default_value_t = 100)]
    batch_rounds: usize,

    /// name2 output mode: full read list + coverage, coverage only, or none
    #[arg(long = "name2-mode", default_value_t = trackcluster_rs::cluster::clusterj::Name2Mode::Coverage)]
    name2_mode: trackcluster_rs::cluster::clusterj::Name2Mode,

    /// SL-supported partial/5' truncation biological 5' offset tolerated for merging.
    #[arg(long = "sl-partial-5prime-offset", default_value_t = trackcluster_rs::cluster::clusterj::DEFAULT_SL_PARTIAL_FIVE_PRIME_END_OFFSET)]
    sl_partial_5prime_offset: u32,

    /// SL-supported same-junction biological 5' offset tolerated for merging.
    #[arg(long = "sl-same-junction-5prime-offset", default_value_t = trackcluster_rs::cluster::clusterj::DEFAULT_SL_SAME_JUNCTION_FIVE_PRIME_END_OFFSET)]
    sl_same_junction_5prime_offset: u32,

    /// Offset used to group SL-supported reads into the same biological 5' cluster.
    #[arg(long = "sl-5prime-cluster-offset", default_value_t = trackcluster_rs::cluster::clusterj::DEFAULT_SL_FIVE_PRIME_CLUSTER_OFFSET)]
    sl_5prime_cluster_offset: u32,

    /// Minimum read support required for an SL 5' cluster to protect a candidate isoform.
    #[arg(long = "sl-5prime-min-support", default_value_t = trackcluster_rs::cluster::clusterj::DEFAULT_MIN_SL_FIVE_PRIME_CLUSTER_SUPPORT)]
    sl_5prime_min_support: usize,

    /// Overwrite existing outputs (default: skip genes whose output file already exists)
    #[arg(long = "force")]
    force: bool,

    /// Print a progress line every N genes
    #[arg(long = "progress-every", default_value_t = 1000)]
    progress_every: usize,

    /// Emit a heartbeat status line every N seconds (0 disables).
    /// Useful when large genes make progress appear "stuck" because no gene completes for a while.
    #[arg(long = "heartbeat-seconds", default_value_t = 60)]
    heartbeat_seconds: u64,

    /// When a heartbeat sees no progress, print up to this many in-flight genes (0 => 1).
    #[arg(long = "heartbeat-top", default_value_t = 5)]
    heartbeat_top: usize,

    /// Restrict downsampling to these gene(s) only (repeatable; exact folder names under --input-root).
    /// If omitted and `--max-reads-per-gene > 0`, downsampling applies to all genes.
    #[arg(long = "downsample-gene")]
    downsample_genes: Vec<String>,

    /// Per-gene cap: reservoir-sample reads down to this count (set to 0 to disable downsampling).
    #[arg(long = "max-reads-per-gene", default_value_t = 50000)]
    max_reads_per_gene: usize,

    /// Deterministic RNG seed used for per-gene downsampling.
    #[arg(long = "downsample-seed", default_value_t = 1)]
    downsample_seed: u64,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    trackcluster_rs::flow::full::run_clusterj_batch(
        trackcluster_rs::flow::full::BatchRunOptions {
            cluster_mode: trackcluster_rs::flow::full::ClusterMode::Clusterj,
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
            name2_mode: args.name2_mode,
            sl_options: trackcluster_rs::cluster::clusterj::SlMergeOptions {
                partial_five_prime_end_offset: args.sl_partial_5prime_offset,
                same_junction_five_prime_end_offset: args.sl_same_junction_5prime_offset,
                five_prime_cluster_offset: args.sl_5prime_cluster_offset,
                min_five_prime_cluster_support: args.sl_5prime_min_support,
            },
            overlap_cutoff1: trackcluster_rs::cluster::cluster_overlap::DEFAULT_CUTOFF1,
            overlap_cutoff2: trackcluster_rs::cluster::cluster_overlap::DEFAULT_CUTOFF2,
            overlap_intron_weight: trackcluster_rs::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT,
            force: args.force,
            progress_every: args.progress_every,
            heartbeat_seconds: args.heartbeat_seconds,
            heartbeat_top: args.heartbeat_top,
            downsample_genes: args.downsample_genes,
            max_reads_per_gene: args.max_reads_per_gene,
            downsample_seed: args.downsample_seed,
        },
    )?;
    Ok(())
}
