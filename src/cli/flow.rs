use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
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
    #[arg(long = "name2-mode", default_value_t = crate::cluster::clusterj::Name2Mode::Coverage)]
    pub name2_mode: crate::cluster::clusterj::Name2Mode,

    /// Prepare step: minimum fraction of read span overlapping a reference span
    #[arg(long = "prepare-fraction-read", default_value_t = 0.01)]
    pub prepare_fraction_read: f64,

    /// Prepare step: minimum fraction of reference span overlapping a read span
    #[arg(long = "prepare-fraction-ref", default_value_t = 0.05)]
    pub prepare_fraction_ref: f64,

    /// Manifest mode: also write `<prefix>_pooled_reads.bed` for debugging/compatibility
    #[arg(long = "emit-pooled-reads")]
    pub emit_pooled_reads: bool,

    /// Overwrite existing outputs (default: skip genes whose output file already exists)
    #[arg(long = "force")]
    pub force: bool,

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

    /// Restrict downsampling to these gene(s) only (repeatable; exact gene folder names).
    /// If omitted and `--max-reads-per-gene > 0`, downsampling applies to all genes.
    #[arg(long = "downsample-gene")]
    pub downsample_genes: Vec<String>,

    /// Per-gene cap: reservoir-sample reads down to this count (set to 0 to disable downsampling).
    #[arg(long = "max-reads-per-gene", default_value_t = 50000)]
    pub max_reads_per_gene: usize,

    /// Deterministic RNG seed used for per-gene downsampling.
    #[arg(long = "downsample-seed", default_value_t = 1)]
    pub downsample_seed: u64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match (&args.reads, &args.manifest) {
        (Some(_), Some(_)) => anyhow::bail!("flow: use either --reads or --manifest, not both"),
        (None, None) => anyhow::bail!("flow: one of --reads or --manifest is required"),
        _ => {}
    }

    let result = crate::flow::full::run_full_flow(crate::flow::full::FullFlowOptions {
        reads: args.reads,
        manifest: args.manifest,
        reference: args.reference,
        output_root: args.output_root,
        prefix: args.prefix,
        threads: args.threads,
        sw_score: args.sw_score,
        batch_size: args.batch_size,
        batch_rounds: args.batch_rounds,
        name2_mode: args.name2_mode,
        prepare_fraction_read: args.prepare_fraction_read,
        prepare_fraction_ref: args.prepare_fraction_ref,
        emit_pooled_reads: args.emit_pooled_reads,
        force: args.force,
        progress_every: args.progress_every,
        heartbeat_seconds: args.heartbeat_seconds,
        heartbeat_top: args.heartbeat_top,
        downsample_genes: args.downsample_genes,
        max_reads_per_gene: args.max_reads_per_gene,
        downsample_seed: args.downsample_seed,
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
    if let Some(multi) = result.multi_sample {
        eprintln!(
            "flow: multi-sample long={:?} matrix={:?} group={:?}",
            multi.long_tsv, multi.matrix_tsv, multi.group_tsv
        );
    }

    Ok(())
}
