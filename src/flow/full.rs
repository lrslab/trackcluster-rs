use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use anyhow::Context;

use crate::annotate::addgene::AddGeneOpts;
use crate::annotate::desc::{describe, DescOpts};
use crate::count::multi::{run_count_multi_from_read_to_isoform, MultiSampleOutputPaths};
use crate::count::{count_by_read_to_isoform, write_counts_csv};
use crate::flow::preparedir::{
    prepare_dir_from_manifest_rows, prepare_dir_from_paths, PrepareDirResult,
};
use crate::io::bed::{read_bed12, BedError};
use crate::io::manifest::{read_manifest_tsv, SampleRow};
use crate::model::Transcript;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ClusterMode {
    #[default]
    Clusterj,
    Cluster,
}

impl ClusterMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Clusterj => "clusterj",
            Self::Cluster => "cluster",
        }
    }

    fn batch_log_label(self) -> &'static str {
        match self {
            Self::Clusterj => "clusterj-batch",
            Self::Cluster => "cluster-batch",
        }
    }

    fn batch_file_prefix(self) -> &'static str {
        match self {
            Self::Clusterj => "clusterj_batch",
            Self::Cluster => "cluster_batch",
        }
    }

    fn per_gene_isoform_suffix(self) -> &'static str {
        match self {
            Self::Clusterj => "_simple_coveragej.bed",
            Self::Cluster => "_simple_coverage.bed",
        }
    }
}

impl std::fmt::Display for ClusterMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for ClusterMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("clusterj") {
            return Ok(Self::Clusterj);
        }
        if s.eq_ignore_ascii_case("cluster") {
            return Ok(Self::Cluster);
        }
        Err(format!(
            "invalid cluster mode {s:?}; expected one of: clusterj, cluster"
        ))
    }
}

#[derive(Clone, Debug)]
pub struct BatchRunOptions {
    pub cluster_mode: ClusterMode,
    pub prepare_reads: Option<PathBuf>,
    pub prepare_reference: Option<PathBuf>,
    pub prepare_prefix: Option<String>,
    pub prepare_fraction_read: f64,
    pub prepare_fraction_ref: f64,
    pub input_root: PathBuf,
    pub gene_list: Option<PathBuf>,
    pub output_root: PathBuf,
    pub threads: usize,
    pub sw_score: i64,
    pub batch_size: usize,
    pub batch_rounds: usize,
    pub name2_mode: crate::cluster::clusterj::Name2Mode,
    pub platform_preset: crate::cluster::clusterj::PlatformPreset,
    pub junction_correction_options: crate::cluster::clusterj::JunctionCorrectionOptions,
    pub sl_options: crate::cluster::clusterj::SlMergeOptions,
    pub overlap_cutoff1: f64,
    pub overlap_cutoff2: f64,
    pub overlap_intron_weight: f64,
    pub force: bool,
    pub progress_every: usize,
    /// Emit a heartbeat status line every N seconds (0 disables).
    pub heartbeat_seconds: u64,
    /// When a heartbeat sees no progress, print up to this many in-flight genes (0 => 1).
    pub heartbeat_top: usize,
    /// Per-gene downsampling: if non-empty, only downsample these gene folders (exact names).
    /// If empty and `max_reads_per_gene > 0`, downsampling applies to all genes.
    pub downsample_genes: Vec<String>,
    /// Per-gene downsampling: cap reads per selected gene to this many (0 disables).
    pub max_reads_per_gene: usize,
    /// Per-gene downsampling: deterministic RNG seed.
    pub downsample_seed: u64,
}

#[derive(Clone, Debug)]
pub struct BatchRunResult {
    pub prepared: Option<PrepareDirResult>,
    pub total_genes: usize,
    pub processed: usize,
    pub skipped: usize,
    pub errors: usize,
    pub elapsed_seconds: f64,
    pub summary_path: PathBuf,
    pub error_path: PathBuf,
    pub downsample_path: PathBuf,
}

#[derive(Clone, Debug)]
pub struct FullFlowOptions {
    pub cluster_mode: ClusterMode,
    pub reads: Option<PathBuf>,
    pub manifest: Option<PathBuf>,
    pub reference: PathBuf,
    pub output_root: PathBuf,
    pub prefix: String,
    pub threads: usize,
    pub sw_score: i64,
    pub batch_size: usize,
    pub batch_rounds: usize,
    pub name2_mode: crate::cluster::clusterj::Name2Mode,
    pub platform_preset: crate::cluster::clusterj::PlatformPreset,
    pub junction_correction_options: crate::cluster::clusterj::JunctionCorrectionOptions,
    pub sl_options: crate::cluster::clusterj::SlMergeOptions,
    pub overlap_cutoff1: f64,
    pub overlap_cutoff2: f64,
    pub overlap_intron_weight: f64,
    pub prepare_fraction_read: f64,
    pub prepare_fraction_ref: f64,
    pub assignment_mode: crate::count::AssignmentMode,
    pub emit_pooled_reads: bool,
    pub force: bool,
    pub progress_every: usize,
    /// Emit a heartbeat status line every N seconds during per-gene clustering (0 disables).
    pub heartbeat_seconds: u64,
    /// When a heartbeat sees no progress, print up to this many in-flight genes (0 => 1).
    pub heartbeat_top: usize,
    /// Per-gene downsampling: if non-empty, only downsample these gene folders (exact names).
    /// If empty and `max_reads_per_gene > 0`, downsampling applies to all genes.
    pub downsample_genes: Vec<String>,
    /// Per-gene downsampling: cap reads per selected gene to this many (0 disables).
    pub max_reads_per_gene: usize,
    /// Per-gene downsampling: deterministic RNG seed.
    pub downsample_seed: u64,
}

#[derive(Clone, Debug)]
pub struct FullFlowResult {
    pub batch: BatchRunResult,
    pub isoform_bed: PathBuf,
    pub unused_bed: PathBuf,
    pub count_csv: PathBuf,
    pub desc_prefix: PathBuf,
    pub multi_sample: Option<MultiSampleOutputPaths>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneOutcome {
    Processed,
    Skipped,
}

#[derive(Clone, Copy, Debug)]
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    fn gen_below(&mut self, upper: usize) -> usize {
        debug_assert!(upper > 0);
        (self.next_u64() % (upper as u64)) as usize
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 14695981039346656037;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(1099511628211);
    }
    hash
}

fn seed_for_gene(base_seed: u64, gene: &str) -> u64 {
    base_seed ^ fnv1a64(gene.as_bytes())
}

#[derive(Clone, Debug)]
struct DownsampleRecord {
    gene: String,
    original_reads: usize,
    sampled_reads: usize,
    target_reads: usize,
    seed: u64,
    scale_factor: f64,
}

#[derive(Clone, Debug)]
struct ProcessGeneResult {
    outcome: GeneOutcome,
    downsample: Option<DownsampleRecord>,
}

fn should_downsample_gene(gene: &str, args: &BatchRunOptions) -> bool {
    if args.max_reads_per_gene == 0 {
        return false;
    }
    if args.downsample_genes.is_empty() {
        return true;
    }
    args.downsample_genes.iter().any(|g| g == gene)
}

fn reservoir_sample_reads(
    path: &Path,
    target_reads: usize,
    seed: u64,
) -> anyhow::Result<(Vec<Transcript>, usize)> {
    let mut rng = Lcg64::new(seed);
    // Avoid huge preallocation when downsampling is enabled by default. The vector will grow to
    // `target_reads` as needed, but most genes are far smaller than the cap.
    let mut sampled: Vec<Transcript> = Vec::with_capacity(target_reads.min(4096));
    let mut total_reads = 0usize;

    let reader = read_bed12(path).with_context(|| format!("open reads {path:?}"))?;
    for record in reader {
        let tx = record.with_context(|| format!("parse reads {path:?}"))?;
        total_reads += 1;

        if sampled.len() < target_reads {
            sampled.push(tx);
            continue;
        }

        let idx = rng.gen_below(total_reads);
        if idx < target_reads {
            sampled[idx] = tx;
        }
    }

    Ok((sampled, total_reads))
}

fn write_downsample_tsv(path: &Path, record: &DownsampleRecord) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(
        fs::File::create(path).with_context(|| format!("write downsample info {path:?}"))?,
    );
    writeln!(
        writer,
        "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}",
        record.gene,
        record.original_reads,
        record.sampled_reads,
        record.scale_factor,
        record.seed,
        record.target_reads
    )?;
    Ok(())
}

fn read_downsample_scales(path: &Path) -> anyhow::Result<HashMap<String, f64>> {
    if !path.exists() {
        return Ok(HashMap::new());
    }

    let file = fs::File::open(path).with_context(|| format!("open downsample file {path:?}"))?;
    let reader = BufReader::new(file);
    let mut out: HashMap<String, f64> = HashMap::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("read downsample file {path:?}"))?;
        if line_no == 0 {
            continue;
        }
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields = line.split('\t');
        let Some(gene) = fields.next() else {
            continue;
        };
        let _original_reads = fields.next();
        let _sampled_reads = fields.next();
        let Some(scale_factor) = fields.next() else {
            continue;
        };

        let gene = gene.trim();
        if gene.is_empty() || gene == "none" {
            continue;
        }
        let scale_factor: f64 = scale_factor.trim().parse().with_context(|| {
            format!("parse scale_factor {scale_factor:?} for gene {gene:?} in {path:?}")
        })?;
        if scale_factor > 0.0 {
            out.insert(gene.to_owned(), scale_factor);
        }
    }
    Ok(out)
}

fn scale_for_gene_field(gene_field: &str, scales: &HashMap<String, f64>) -> Option<f64> {
    let gene_field = gene_field.trim();
    if gene_field.is_empty() || gene_field == "none" {
        return None;
    }

    let mut scale: Option<f64> = None;
    for gene in gene_field
        .split("||")
        .map(str::trim)
        .filter(|g| !g.is_empty() && *g != "none")
    {
        let Some(candidate) = scales.get(gene).copied() else {
            continue;
        };
        if scale.replace(candidate).is_some() {
            return None;
        }
    }
    scale
}

fn read_gene_list(path: &Path) -> anyhow::Result<Vec<String>> {
    let file = fs::File::open(path).with_context(|| format!("open gene list {path:?}"))?;
    let reader = BufReader::new(file);
    let mut genes: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("read gene list {path:?}"))?;
        let gene = line.trim();
        if gene.is_empty() || gene.starts_with('#') {
            continue;
        }
        genes.push(gene.to_owned());
    }
    genes.sort();
    genes.dedup();
    Ok(genes)
}

fn discover_genes(root: &Path) -> anyhow::Result<Vec<String>> {
    let mut genes: Vec<String> = Vec::new();
    for entry in fs::read_dir(root).with_context(|| format!("read_dir {root:?}"))? {
        let entry = entry.with_context(|| format!("read_dir entry in {root:?}"))?;
        let ty = entry
            .file_type()
            .with_context(|| format!("stat {:?}", entry.path()))?;
        if !ty.is_dir() {
            continue;
        }
        genes.push(entry.file_name().to_string_lossy().to_string());
    }
    genes.sort();
    genes.dedup();
    Ok(genes)
}

fn per_gene_outputs_complete(out_isoforms: &Path, out_unused: &Path, out_mapping: &Path) -> bool {
    out_isoforms.exists() && out_unused.exists() && out_mapping.exists()
}

fn process_gene(gene: &str, args: &BatchRunOptions) -> anyhow::Result<ProcessGeneResult> {
    let gene_dir = args.input_root.join(gene);
    let reads = gene_dir.join(format!("{gene}_nano.bed"));
    let reference = gene_dir.join(format!("{gene}_gff.bed"));
    if !reads.exists() || !reference.exists() {
        return Ok(ProcessGeneResult {
            outcome: GeneOutcome::Skipped,
            downsample: None,
        });
    }

    let reads_len = fs::metadata(&reads)
        .with_context(|| format!("stat reads {reads:?}"))?
        .len();
    if reads_len == 0 {
        return Ok(ProcessGeneResult {
            outcome: GeneOutcome::Skipped,
            downsample: None,
        });
    }

    let out_dir = args.output_root.join(gene);
    fs::create_dir_all(&out_dir).with_context(|| format!("create {out_dir:?}"))?;

    let out_isoforms = out_dir.join(format!(
        "{gene}{}",
        args.cluster_mode.per_gene_isoform_suffix()
    ));
    let out_unused = out_dir.join(format!("{gene}_unused.bed"));
    let out_mapping = out_dir.join(format!("{gene}_read_to_isoform.tsv"));

    if !args.force && per_gene_outputs_complete(&out_isoforms, &out_unused, &out_mapping) {
        return Ok(ProcessGeneResult {
            outcome: GeneOutcome::Skipped,
            downsample: None,
        });
    }

    let (reads, downsample) = if should_downsample_gene(gene, args) {
        let target_reads = args.max_reads_per_gene.max(1);
        let seed = seed_for_gene(args.downsample_seed, gene);
        let (sampled, total_reads) = reservoir_sample_reads(&reads, target_reads, seed)?;
        let sampled_reads = sampled.len();
        let downsample = if total_reads > sampled_reads && sampled_reads > 0 {
            Some(DownsampleRecord {
                gene: gene.to_owned(),
                original_reads: total_reads,
                sampled_reads,
                target_reads,
                seed,
                scale_factor: total_reads as f64 / sampled_reads as f64,
            })
        } else {
            None
        };
        (sampled, downsample)
    } else {
        let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&reads)
            .with_context(|| format!("open reads {reads:?}"))?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
            .with_context(|| format!("parse reads {reads:?}"))?;
        (reads, None)
    };

    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&reference)
        .with_context(|| format!("open reference {reference:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reference {reference:?}"))?;

    let result = match args.cluster_mode {
        ClusterMode::Clusterj => crate::cluster::clusterj::clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            args.sw_score,
            args.batch_size,
            args.batch_rounds,
            args.name2_mode,
            args.sl_options,
            args.junction_correction_options,
        ),
        ClusterMode::Cluster => crate::cluster::cluster_overlap::cluster_with_options(
            &reads,
            Some(&refs),
            1,
            crate::cluster::cluster_overlap::ClusterOptions {
                cutoff1: args.overlap_cutoff1,
                cutoff2: args.overlap_cutoff2,
                intron_weight: args.overlap_intron_weight,
                sw_score: args.sw_score,
                name2_mode: args.name2_mode,
                batch_size: args.batch_size,
                batch_rounds: args.batch_rounds,
            },
        ),
    };

    crate::cluster::output::write_isoforms_bed(&out_isoforms, &result.isoforms)
        .with_context(|| format!("write {out_isoforms:?}"))?;
    crate::cluster::output::write_isoforms_bed(&out_unused, &result.unused)
        .with_context(|| format!("write {out_unused:?}"))?;
    crate::cluster::output::write_read_to_isoform_tsv(&out_mapping, &result.read_to_isoform)
        .with_context(|| format!("write {out_mapping:?}"))?;

    if let Some(record) = downsample.as_ref() {
        let path = out_dir.join("downsample.tsv");
        write_downsample_tsv(&path, record)
            .with_context(|| format!("write per-gene downsample info {path:?}"))?;
    }

    Ok(ProcessGeneResult {
        outcome: GeneOutcome::Processed,
        downsample,
    })
}

pub fn run_clusterj_batch(args: BatchRunOptions) -> anyhow::Result<BatchRunResult> {
    fs::create_dir_all(&args.output_root)
        .with_context(|| format!("create {:?}", args.output_root))?;

    let batch_log_label = args.cluster_mode.batch_log_label();
    let batch_file_prefix = args.cluster_mode.batch_file_prefix();

    let prepared = if args.prepare_reads.is_some()
        || args.prepare_reference.is_some()
        || args.prepare_prefix.is_some()
    {
        let reads = args
            .prepare_reads
            .as_ref()
            .context("--prepare-reads is required when using prepare")?;
        let reference = args
            .prepare_reference
            .as_ref()
            .context("--prepare-reference is required when using prepare")?;
        let prefix = args
            .prepare_prefix
            .as_deref()
            .context("--prepare-prefix is required when using prepare")?;

        eprintln!(
            "prepare: reads={:?} reference={:?} input_root={:?} prefix={}",
            reads, reference, args.input_root, prefix
        );
        let res = prepare_dir_from_paths(
            reads,
            reference,
            &args.input_root,
            prefix,
            AddGeneOpts {
                fraction_read: args.prepare_fraction_read,
                fraction_ref: args.prepare_fraction_ref,
            },
        )?;
        eprintln!(
            "prepare: genes={}, dedup_reads={}, novel_reads={}",
            res.genes.len(),
            res.dedup_reads,
            res.novel_reads
        );
        Some(res)
    } else {
        None
    };

    let mut genes = match &args.gene_list {
        Some(list) => read_gene_list(list)?,
        None => discover_genes(&args.input_root)?,
    };
    genes.sort();
    genes.dedup();

    let total = genes.len();
    if total == 0 {
        anyhow::bail!("no genes found (input-root={:?})", args.input_root);
    }

    let started = Instant::now();
    eprintln!(
        "{batch_log_label}: {} genes, {} worker threads",
        total,
        args.threads.max(1)
    );

    if args.threads > 1 && args.max_reads_per_gene == 0 {
        eprintln!(
            "{batch_log_label}: note: --threads > 1 with --max-reads-per-gene=0 can use a lot of memory on large genes; \
consider --max-reads-per-gene and/or --name2-mode coverage"
        );
    }
    if args.cluster_mode == ClusterMode::Cluster && args.batch_size == 0 {
        eprintln!(
            "{batch_log_label}: note: overlap mode batching is disabled because --batch-size=0; \
each gene will run one full two-pass overlap merge"
        );
    }

    let processed = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));
    let errors = Arc::new(AtomicUsize::new(0));
    let done = Arc::new(AtomicUsize::new(0));
    let error_lines: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let downsample_records: Arc<Mutex<Vec<DownsampleRecord>>> = Arc::new(Mutex::new(Vec::new()));

    let queue = Arc::new(Mutex::new(genes));
    let worker_count = args.threads.max(1).min(total);

    #[derive(Debug, Default)]
    struct WorkerState {
        gene: Option<String>,
        started_at: Option<Instant>,
    }

    let worker_states: Arc<Vec<Mutex<WorkerState>>> = Arc::new(
        (0..worker_count)
            .map(|_| Mutex::new(WorkerState::default()))
            .collect(),
    );

    let (heartbeat_stop_tx, heartbeat_handle) = if args.heartbeat_seconds > 0 {
        use std::sync::mpsc::{self, RecvTimeoutError};
        let (stop_tx, stop_rx) = mpsc::channel::<()>();

        let processed = Arc::clone(&processed);
        let skipped = Arc::clone(&skipped);
        let errors = Arc::clone(&errors);
        let done = Arc::clone(&done);
        let queue = Arc::clone(&queue);
        let worker_states = Arc::clone(&worker_states);
        let heartbeat_seconds = args.heartbeat_seconds;
        let heartbeat_top = args.heartbeat_top.max(1);

        let handle = std::thread::spawn(move || {
            let mut last_done = done.load(Ordering::Relaxed);
            loop {
                match stop_rx.recv_timeout(Duration::from_secs(heartbeat_seconds)) {
                    Ok(()) => break,
                    Err(RecvTimeoutError::Disconnected) => break,
                    Err(RecvTimeoutError::Timeout) => {}
                }

                let done_now = done.load(Ordering::Relaxed);
                let processed_now = processed.load(Ordering::Relaxed);
                let skipped_now = skipped.load(Ordering::Relaxed);
                let errors_now = errors.load(Ordering::Relaxed);
                let queue_remaining = queue.lock().map(|guard| guard.len()).unwrap_or_default();

                eprintln!(
                    "heartbeat: {done_now}/{total} (processed={processed_now}, skipped={skipped_now}, errors={errors_now}) queue_remaining={queue_remaining} elapsed={:?}",
                    started.elapsed()
                );

                if done_now == last_done && done_now < total {
                    let mut inflight: Vec<(String, Duration)> = Vec::new();
                    for state_lock in worker_states.iter() {
                        let Ok(state) = state_lock.lock() else {
                            continue;
                        };
                        let (Some(gene), Some(started_at)) =
                            (state.gene.as_ref(), state.started_at)
                        else {
                            continue;
                        };
                        inflight.push((gene.clone(), started_at.elapsed()));
                    }

                    if inflight.is_empty() {
                        eprintln!("heartbeat: no in-flight genes (all workers idle?)");
                    } else {
                        inflight.sort_by(|a, b| b.1.cmp(&a.1));
                        let top = heartbeat_top.min(inflight.len());
                        let mut line = String::from("in_flight(top):");
                        for (gene, dur) in inflight.into_iter().take(top) {
                            line.push(' ');
                            line.push_str(&format!("{gene}={:.1}s", dur.as_secs_f64()));
                        }
                        eprintln!("{line}");
                    }
                }

                last_done = done_now;
                if done_now >= total {
                    break;
                }
            }
        });

        (Some(stop_tx), Some(handle))
    } else {
        (None, None)
    };

    let mut handles = Vec::with_capacity(worker_count);
    for worker_idx in 0..worker_count {
        let queue = Arc::clone(&queue);
        let processed = Arc::clone(&processed);
        let skipped = Arc::clone(&skipped);
        let errors = Arc::clone(&errors);
        let done = Arc::clone(&done);
        let error_lines = Arc::clone(&error_lines);
        let downsample_records = Arc::clone(&downsample_records);
        let worker_states = Arc::clone(&worker_states);
        let args = args.clone();

        handles.push(std::thread::spawn(move || loop {
            let gene = {
                let mut guard = queue.lock().expect("work queue poisoned");
                guard.pop()
            };
            let Some(gene) = gene else {
                break;
            };

            if let Ok(mut state) = worker_states[worker_idx].lock() {
                state.gene = Some(gene.clone());
                state.started_at = Some(Instant::now());
            }

            let result = match panic::catch_unwind(AssertUnwindSafe(|| process_gene(&gene, &args)))
            {
                Ok(Ok(result)) => result,
                Ok(Err(err)) => {
                    errors.fetch_add(1, Ordering::Relaxed);
                    if let Ok(mut guard) = error_lines.lock() {
                        guard.push(format!("{gene}\t{err}"));
                    }
                    ProcessGeneResult {
                        outcome: GeneOutcome::Skipped,
                        downsample: None,
                    }
                }
                Err(payload) => {
                    errors.fetch_add(1, Ordering::Relaxed);
                    let msg = payload
                        .downcast_ref::<&str>()
                        .copied()
                        .map(str::to_owned)
                        .or_else(|| payload.downcast_ref::<String>().cloned())
                        .unwrap_or_else(|| "<non-string panic>".to_owned());
                    if let Ok(mut guard) = error_lines.lock() {
                        guard.push(format!("{gene}\tpanic\t{msg}"));
                    }
                    ProcessGeneResult {
                        outcome: GeneOutcome::Skipped,
                        downsample: None,
                    }
                }
            };

            if let Ok(mut state) = worker_states[worker_idx].lock() {
                state.gene = None;
                state.started_at = None;
            }

            if let Some(record) = result.downsample {
                if let Ok(mut guard) = downsample_records.lock() {
                    guard.push(record);
                }
            }

            match result.outcome {
                GeneOutcome::Processed => {
                    processed.fetch_add(1, Ordering::Relaxed);
                }
                GeneOutcome::Skipped => {
                    skipped.fetch_add(1, Ordering::Relaxed);
                }
            }

            let done_now = done.fetch_add(1, Ordering::Relaxed) + 1;
            let every = args.progress_every.max(1);
            if done_now.is_multiple_of(every) || done_now == total {
                eprintln!(
                    "progress: {done_now}/{total} (processed={}, skipped={}, errors={})",
                    processed.load(Ordering::Relaxed),
                    skipped.load(Ordering::Relaxed),
                    errors.load(Ordering::Relaxed),
                );
            }
        }));
    }

    for handle in handles {
        if handle.join().is_err() {
            errors.fetch_add(1, Ordering::Relaxed);
        }
    }

    if let Some(tx) = heartbeat_stop_tx {
        let _ = tx.send(());
    }
    if let Some(handle) = heartbeat_handle {
        let _ = handle.join();
    }

    let elapsed = started.elapsed();
    let processed = processed.load(Ordering::Relaxed);
    let skipped = skipped.load(Ordering::Relaxed);
    let errors = errors.load(Ordering::Relaxed);

    eprintln!(
        "done: processed={processed}, skipped={skipped}, errors={errors}, elapsed={:?}",
        elapsed
    );

    let summary_path = args
        .output_root
        .join(format!("{batch_file_prefix}_summary.txt"));
    let mut summary =
        fs::File::create(&summary_path).with_context(|| format!("write {summary_path:?}"))?;
    writeln!(summary, "input_root\t{:?}", args.input_root)?;
    writeln!(summary, "gene_list\t{:?}", args.gene_list)?;
    writeln!(summary, "output_root\t{:?}", args.output_root)?;
    writeln!(summary, "cluster_mode\t{}", args.cluster_mode)?;
    writeln!(summary, "threads\t{}", args.threads)?;
    writeln!(summary, "sw_score\t{}", args.sw_score)?;
    writeln!(summary, "batch_size\t{}", args.batch_size)?;
    writeln!(summary, "batch_rounds\t{}", args.batch_rounds)?;
    writeln!(summary, "name2_mode\t{}", args.name2_mode)?;
    writeln!(summary, "platform_preset\t{}", args.platform_preset)?;
    writeln!(
        summary,
        "junction_correction_offset\t{}",
        args.junction_correction_options.offset
    )?;
    writeln!(
        summary,
        "junction_correction_min_support\t{}",
        args.junction_correction_options.min_support
    )?;
    writeln!(
        summary,
        "sl_partial_5prime_offset\t{}",
        args.sl_options.partial_five_prime_end_offset
    )?;
    writeln!(
        summary,
        "sl_same_junction_5prime_offset\t{}",
        args.sl_options.same_junction_five_prime_end_offset
    )?;
    writeln!(
        summary,
        "sl_5prime_cluster_offset\t{}",
        args.sl_options.five_prime_cluster_offset
    )?;
    writeln!(
        summary,
        "sl_5prime_min_support\t{}",
        args.sl_options.min_five_prime_cluster_support
    )?;
    writeln!(summary, "overlap_cutoff1\t{}", args.overlap_cutoff1)?;
    writeln!(summary, "overlap_cutoff2\t{}", args.overlap_cutoff2)?;
    writeln!(
        summary,
        "overlap_intron_weight\t{}",
        args.overlap_intron_weight
    )?;
    writeln!(summary, "force\t{}", args.force)?;
    writeln!(summary, "progress_every\t{}", args.progress_every)?;
    writeln!(summary, "heartbeat_seconds\t{}", args.heartbeat_seconds)?;
    writeln!(summary, "heartbeat_top\t{}", args.heartbeat_top)?;
    writeln!(summary, "max_reads_per_gene\t{}", args.max_reads_per_gene)?;
    writeln!(summary, "downsample_seed\t{}", args.downsample_seed)?;
    if args.downsample_genes.is_empty() {
        writeln!(summary, "downsample_genes\t[]")?;
    } else {
        writeln!(
            summary,
            "downsample_genes\t{}",
            args.downsample_genes.join(",")
        )?;
    }
    writeln!(summary, "total_genes\t{}", total)?;
    writeln!(summary, "processed\t{}", processed)?;
    writeln!(summary, "skipped\t{}", skipped)?;
    writeln!(summary, "errors\t{}", errors)?;
    writeln!(summary, "elapsed_seconds\t{}", elapsed.as_secs_f64())?;

    let error_path = args
        .output_root
        .join(format!("{batch_file_prefix}_errors.txt"));
    let error_lines = error_lines
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if !error_lines.is_empty() {
        let mut errors_file =
            fs::File::create(&error_path).with_context(|| format!("write {error_path:?}"))?;
        for line in error_lines.iter() {
            writeln!(errors_file, "{line}")?;
        }
    } else if error_path.exists() {
        let _ = fs::remove_file(&error_path);
    }

    let downsample_path = args
        .output_root
        .join(format!("{batch_file_prefix}_downsample.tsv"));
    let mut downsample_records = downsample_records
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if !downsample_records.is_empty() {
        downsample_records.sort_by(|a, b| a.gene.cmp(&b.gene));
        let mut writer = std::io::BufWriter::new(
            fs::File::create(&downsample_path)
                .with_context(|| format!("write {downsample_path:?}"))?,
        );
        writeln!(
            writer,
            "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads"
        )?;
        for record in downsample_records.iter() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                record.gene,
                record.original_reads,
                record.sampled_reads,
                record.scale_factor,
                record.seed,
                record.target_reads
            )?;
        }
    } else if downsample_path.exists() {
        let _ = fs::remove_file(&downsample_path);
    }

    Ok(BatchRunResult {
        prepared,
        total_genes: total,
        processed,
        skipped,
        errors,
        elapsed_seconds: elapsed.as_secs_f64(),
        summary_path,
        error_path,
        downsample_path,
    })
}

fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

fn merge_files(inputs: &[PathBuf], out: &Path) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(
        fs::File::create(out).with_context(|| format!("create merged output {out:?}"))?,
    );
    let mut buffer = vec![0u8; 1024 * 1024];
    for input in inputs {
        let mut reader = std::io::BufReader::new(
            fs::File::open(input).with_context(|| format!("open {input:?}"))?,
        );
        let mut saw_bytes = false;
        let mut last_byte: u8 = b'\n';
        loop {
            let read_len = reader
                .read(&mut buffer)
                .with_context(|| format!("read {input:?}"))?;
            if read_len == 0 {
                break;
            }
            saw_bytes = true;
            last_byte = buffer[read_len - 1];
            writer
                .write_all(&buffer[..read_len])
                .with_context(|| format!("append {input:?} into {out:?}"))?;
        }
        if saw_bytes && last_byte != b'\n' {
            writer
                .write_all(b"\n")
                .with_context(|| format!("final newline after {input:?}"))?;
        }
    }
    Ok(())
}

fn merge_gene_outputs(
    output_root: &Path,
    genes: &[String],
    per_gene_suffix: &str,
    merged_out: &Path,
) -> anyhow::Result<()> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for gene in genes {
        let path = output_root
            .join(gene)
            .join(format!("{gene}{per_gene_suffix}"));
        if path.exists() {
            paths.push(path);
        }
    }
    paths.sort();
    merge_files(&paths, merged_out)?;
    Ok(())
}

fn read_bed12_records(path: &Path, kind: &str) -> anyhow::Result<Vec<Transcript>> {
    read_bed12(path)
        .with_context(|| format!("open {kind} {path:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse {kind} {path:?}"))
}

fn run_count_and_desc(
    isoforms: &[Transcript],
    refs: &[Transcript],
    read_to_isoform: &[(String, String)],
    count_csv: &Path,
    desc_prefix: &Path,
    downsample_scales: Option<&HashMap<String, f64>>,
) -> anyhow::Result<()> {
    let mut counts = count_by_read_to_isoform(isoforms, read_to_isoform);
    if let Some(scales) = downsample_scales {
        const GENE_NAME_COL: usize = 5;
        for (record, isoform) in counts.iter_mut().zip(isoforms.iter()) {
            let gene_field = isoform
                .extra_fields
                .get(GENE_NAME_COL)
                .map(|value| value.as_str())
                .unwrap_or("")
                .trim();
            if let Some(scale) = scale_for_gene_field(gene_field, scales) {
                record.count *= scale;
            }
        }
    }
    write_counts_csv(count_csv, &counts)
        .with_context(|| format!("write count csv {count_csv:?}"))?;

    let desc = describe(isoforms, refs, DescOpts::default());

    let desc_path = append_suffix(desc_prefix, "_desc.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&desc_path).with_context(|| format!("write {desc_path:?}"))?,
    );
    for row in &desc.desc_rows {
        writeln!(
            &mut writer,
            "{}\t{}\t{}\t{}\t{}",
            row.isoform_id, row.ref_id, row.gene, row.miss, row.extra
        )?;
    }

    let class4_path = append_suffix(desc_prefix, "_class4.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&class4_path).with_context(|| format!("write {class4_path:?}"))?,
    );
    for row in &desc.class4_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    let fusion_path = append_suffix(desc_prefix, "_fusion.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&fusion_path).with_context(|| format!("write {fusion_path:?}"))?,
    );
    for row in &desc.fusion_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.genes.join(";"))?;
    }

    let class12_path = append_suffix(desc_prefix, "_class12.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&class12_path).with_context(|| format!("write {class12_path:?}"))?,
    );
    for row in &desc.class12_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    Ok(())
}

fn run_count_multi_scaled(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    out_prefix: &Path,
    downsample_scales: &HashMap<String, f64>,
) -> anyhow::Result<MultiSampleOutputPaths> {
    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }

    let mut result =
        crate::count::multi::count_multi_by_read_to_isoform(isoforms, read_to_isoform, sample_rows)
            .with_context(|| format!("count-multi {:?}", out_prefix))?;
    let include_group = sample_rows.iter().any(|sample| sample.group.is_some());

    for row in &mut result.matrix_rows {
        let Some(scale) = scale_for_gene_field(&row.gene, downsample_scales) else {
            continue;
        };
        for count in &mut row.counts {
            *count *= scale;
        }
    }

    for row in &mut result.long_rows {
        let Some(scale) = scale_for_gene_field(&row.gene, downsample_scales) else {
            continue;
        };
        row.count *= scale;
        row.gene_total *= scale;
    }

    for row in &mut result.group_rows {
        let Some(scale) = scale_for_gene_field(&row.gene, downsample_scales) else {
            continue;
        };
        row.count *= scale;
        row.gene_total *= scale;
    }

    let long_tsv = append_suffix(out_prefix, ".isoform_usage.long.tsv");
    let matrix_tsv = append_suffix(out_prefix, ".isoform_counts.matrix.tsv");
    crate::count::multi::write_usage_long_tsv(&long_tsv, &result.long_rows, include_group)
        .with_context(|| format!("write long output {long_tsv:?}"))?;
    crate::count::multi::write_counts_matrix_tsv(&matrix_tsv, &result.matrix_rows, sample_rows)
        .with_context(|| format!("write matrix output {matrix_tsv:?}"))?;

    let group_tsv = if result.group_rows.is_empty() {
        None
    } else {
        let path = append_suffix(out_prefix, ".isoform_usage.group.tsv");
        crate::count::multi::write_group_usage_tsv(&path, &result.group_rows)
            .with_context(|| format!("write group output {path:?}"))?;
        Some(path)
    };

    Ok(MultiSampleOutputPaths {
        long_tsv,
        matrix_tsv,
        group_tsv,
    })
}

pub fn run_full_flow(opts: FullFlowOptions) -> anyhow::Result<FullFlowResult> {
    fs::create_dir_all(&opts.output_root)
        .with_context(|| format!("create {:?}", opts.output_root))?;

    let gene_list = opts.output_root.join(format!("{}_gene.txt", opts.prefix));
    let mut sample_rows: Option<Vec<SampleRow>> = None;

    let batch = match (&opts.reads, &opts.manifest) {
        (Some(reads), None) => run_clusterj_batch(BatchRunOptions {
            cluster_mode: opts.cluster_mode,
            prepare_reads: Some(reads.clone()),
            prepare_reference: Some(opts.reference.clone()),
            prepare_prefix: Some(opts.prefix.clone()),
            prepare_fraction_read: opts.prepare_fraction_read,
            prepare_fraction_ref: opts.prepare_fraction_ref,
            input_root: opts.output_root.clone(),
            gene_list: Some(gene_list.clone()),
            output_root: opts.output_root.clone(),
            threads: opts.threads,
            sw_score: opts.sw_score,
            batch_size: opts.batch_size,
            batch_rounds: opts.batch_rounds,
            name2_mode: opts.name2_mode,
            platform_preset: opts.platform_preset,
            junction_correction_options: opts.junction_correction_options,
            sl_options: opts.sl_options,
            overlap_cutoff1: opts.overlap_cutoff1,
            overlap_cutoff2: opts.overlap_cutoff2,
            overlap_intron_weight: opts.overlap_intron_weight,
            force: opts.force,
            progress_every: opts.progress_every,
            heartbeat_seconds: opts.heartbeat_seconds,
            heartbeat_top: opts.heartbeat_top,
            downsample_genes: opts.downsample_genes.clone(),
            max_reads_per_gene: opts.max_reads_per_gene,
            downsample_seed: opts.downsample_seed,
        })?,
        (None, Some(manifest_path)) => {
            let rows = read_manifest_tsv(manifest_path)?;
            let pooled_reads_out = if opts.emit_pooled_reads {
                Some(
                    opts.output_root
                        .join(format!("{}_pooled_reads.bed", opts.prefix)),
                )
            } else {
                None
            };
            let prepared = prepare_dir_from_manifest_rows(
                &rows,
                &opts.reference,
                &opts.output_root,
                &opts.prefix,
                AddGeneOpts {
                    fraction_read: opts.prepare_fraction_read,
                    fraction_ref: opts.prepare_fraction_ref,
                },
                pooled_reads_out.as_deref(),
            )
            .with_context(|| format!("prepare from manifest {manifest_path:?}"))?;
            if let Some(path) = pooled_reads_out {
                eprintln!(
                    "flow: emitted pooled reads from manifest {:?} (samples={}) -> {:?}",
                    manifest_path,
                    rows.len(),
                    path
                );
            }
            eprintln!(
                "prepare: genes={}, dedup_reads={}, novel_reads={}",
                prepared.genes.len(),
                prepared.dedup_reads,
                prepared.novel_reads
            );

            let mut batch = run_clusterj_batch(BatchRunOptions {
                cluster_mode: opts.cluster_mode,
                prepare_reads: None,
                prepare_reference: None,
                prepare_prefix: None,
                prepare_fraction_read: opts.prepare_fraction_read,
                prepare_fraction_ref: opts.prepare_fraction_ref,
                input_root: opts.output_root.clone(),
                gene_list: Some(gene_list.clone()),
                output_root: opts.output_root.clone(),
                threads: opts.threads,
                sw_score: opts.sw_score,
                batch_size: opts.batch_size,
                batch_rounds: opts.batch_rounds,
                name2_mode: opts.name2_mode,
                platform_preset: opts.platform_preset,
                junction_correction_options: opts.junction_correction_options,
                sl_options: opts.sl_options,
                overlap_cutoff1: opts.overlap_cutoff1,
                overlap_cutoff2: opts.overlap_cutoff2,
                overlap_intron_weight: opts.overlap_intron_weight,
                force: opts.force,
                progress_every: opts.progress_every,
                heartbeat_seconds: opts.heartbeat_seconds,
                heartbeat_top: opts.heartbeat_top,
                downsample_genes: opts.downsample_genes.clone(),
                max_reads_per_gene: opts.max_reads_per_gene,
                downsample_seed: opts.downsample_seed,
            })?;
            batch.prepared = Some(prepared);
            sample_rows = Some(rows);
            batch
        }
        (Some(_), Some(_)) => anyhow::bail!("flow: use either reads or manifest, not both"),
        (None, None) => anyhow::bail!("flow: reads or manifest is required"),
    };

    let genes = read_gene_list(&gene_list).with_context(|| format!("read {:?}", gene_list))?;
    let isoform_bed = opts
        .output_root
        .join(format!("{}_isoform.bed", opts.prefix));
    let unused_bed = opts.output_root.join(format!("{}_unused.bed", opts.prefix));
    let read_to_isoform_tsv = opts
        .output_root
        .join(format!("{}_read_to_isoform.tsv", opts.prefix));
    let count_csv = opts
        .output_root
        .join(format!("{}_isoform_count.csv", opts.prefix));
    let desc_prefix = opts.output_root.join(&opts.prefix);

    eprintln!("flow: merge isoforms -> {:?}", isoform_bed);
    merge_gene_outputs(
        &opts.output_root,
        &genes,
        opts.cluster_mode.per_gene_isoform_suffix(),
        &isoform_bed,
    )?;

    eprintln!("flow: merge unused -> {:?}", unused_bed);
    merge_gene_outputs(&opts.output_root, &genes, "_unused.bed", &unused_bed)?;

    eprintln!("flow: merge read-to-isoform -> {:?}", read_to_isoform_tsv);
    merge_gene_outputs(
        &opts.output_root,
        &genes,
        "_read_to_isoform.tsv",
        &read_to_isoform_tsv,
    )?;

    let downsample_path = batch.downsample_path.clone();
    let downsample_scales = read_downsample_scales(&downsample_path)
        .with_context(|| format!("read downsample scales {downsample_path:?}"))?;
    let downsample_scales = if downsample_scales.is_empty() {
        None
    } else {
        eprintln!(
            "flow: applying downsample scale factors (genes={}) from {:?}",
            downsample_scales.len(),
            downsample_path
        );
        Some(downsample_scales)
    };

    let isoforms = read_bed12_records(&isoform_bed, "isoform")?;
    let refs = read_bed12_records(&opts.reference, "reference")?;
    let read_to_isoform = crate::count::read_read_to_isoform_tsv(&read_to_isoform_tsv)
        .with_context(|| format!("read merged read_to_isoform {read_to_isoform_tsv:?}"))?;

    let selected_read_to_isoform;
    let count_read_to_isoform = if opts.assignment_mode == crate::count::AssignmentMode::Unique {
        let reads_for_assignment = if let Some(rows) = sample_rows.as_ref() {
            crate::count::multi::read_tagged_sample_reads(rows)?
        } else if let Some(reads_path) = opts.reads.as_ref() {
            read_bed12_records(reads_path, "reads")?
        } else {
            Vec::new()
        };
        selected_read_to_isoform = crate::count::select_unique_best_read_to_isoform(
            &reads_for_assignment,
            &isoforms,
            &read_to_isoform,
        )?;
        &selected_read_to_isoform
    } else {
        &read_to_isoform
    };

    eprintln!("flow: count + desc");
    run_count_and_desc(
        &isoforms,
        &refs,
        count_read_to_isoform,
        &count_csv,
        &desc_prefix,
        downsample_scales.as_ref(),
    )?;

    let multi_sample = if let Some(rows) = sample_rows.as_ref() {
        let output_prefix = opts.output_root.join(&opts.prefix);
        if let Some(scales) = downsample_scales.as_ref() {
            Some(run_count_multi_scaled(
                rows,
                &isoforms,
                count_read_to_isoform,
                &output_prefix,
                scales,
            )?)
        } else {
            Some(
                run_count_multi_from_read_to_isoform(
                    rows,
                    &isoforms,
                    count_read_to_isoform,
                    &output_prefix,
                )
                .with_context(|| {
                    format!("run multi-sample counting for flow {:?}", output_prefix)
                })?,
            )
        }
    } else {
        None
    };

    Ok(FullFlowResult {
        batch,
        isoform_bed,
        unused_bed,
        count_csv,
        desc_prefix,
        multi_sample,
    })
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn fresh_temp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "trackcluster_rs_full_{}_{}_{}",
            prefix,
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn merge_files_adds_newline_only_when_missing() {
        let dir = fresh_temp_dir("merge_files_newline");
        let in1 = dir.join("in1.txt");
        let in2 = dir.join("in2.txt");
        let in3 = dir.join("in3.txt");
        let out = dir.join("out.txt");

        fs::write(&in1, "first\n").unwrap();
        fs::write(&in2, "second").unwrap();
        fs::write(&in3, "").unwrap();

        merge_files(&[in1, in2, in3], &out).unwrap();
        let merged = fs::read_to_string(out).unwrap();
        assert_eq!(merged, "first\nsecond\n");
    }
}
