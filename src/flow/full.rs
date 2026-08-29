use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::annotate::addgene::AddGeneOpts;
use crate::count::multi::MultiSampleOutputPaths;
use crate::flow::artifact_layout::{
    preflight_gene_artifacts, validate_pipeline_gene_namespace, FlowArtifactLayout,
    GeneArtifactRequirements,
};
use crate::flow::artifact_manifest::{
    assess_cache, atomic_copy, atomic_write_with, invalidate_completion_manifest,
    validate_recorded_completion, write_completion_manifest, CacheDecision, EffectiveOptions,
    InputArtifact, OutputSpec, RecordCountKind, RunManifest, RunRequest, ToolIdentity,
    MANIFEST_FILE_NAME,
};
use crate::flow::config::{
    ClusteringConfig, CountingConfig, DownsampleConfig, InvalidReadPolicy, PrepareConfig,
    RuntimeConfig,
};
use crate::flow::counting::{
    read_bed12_records, run_count_and_desc, run_count_multi_atomic,
    select_unique_read_to_isoform_by_gene,
};
#[cfg(test)]
use crate::flow::merge::{merge_files, merge_isoform_files};
use crate::flow::merge::{merge_gene_isoform_outputs, merge_gene_outputs};
use crate::flow::path_key::{
    ensure_destination_within, gene_artifact_path, gene_dir_path, read_gene_path_map,
    reject_external_inputs_in_output_root, validate_gene_ids, write_gene_id_marker,
    write_gene_path_map, GeneId, SafePathComponent,
};
use crate::flow::preparedir::{
    prepare_dir_from_manifest_rows_with_policy, prepare_dir_from_paths_with_policy,
    PrepareDirResult,
};
use crate::io::bed::{
    count_rejected_reads_tsv, read_bed12, write_rejected_reads_tsv_to_writer, RejectedReadRecord,
};
use crate::io::manifest::read_manifest_tsv;
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

    pub(super) fn per_gene_isoform_suffix(self) -> &'static str {
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
    pub prepare: PrepareConfig,
    /// Rejections already recorded by an upstream manifest preparation step.
    pub prepare_rejected_read_tracks: usize,
    pub input_root: PathBuf,
    pub gene_list: Option<PathBuf>,
    pub output_root: PathBuf,
    pub clustering: ClusteringConfig,
    pub counting: CountingConfig,
    pub runtime: RuntimeConfig,
    pub downsample: DownsampleConfig,
}

impl BatchRunOptions {
    /// Validate all scientific and runtime parameters before any output is created.
    pub fn validate(&self) -> Result<(), crate::config::ParameterError> {
        self.runtime.validate()?;
        self.prepare.validate()?;
        self.clustering.validate()?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct BatchRunResult {
    pub prepared: Option<PrepareDirResult>,
    pub total_genes: usize,
    /// Genes whose complete outputs are safe to merge in this run.
    pub mergeable_genes: Vec<String>,
    pub processed: usize,
    pub skipped: usize,
    pub skipped_completed_outputs: usize,
    pub skipped_empty_reads: usize,
    pub skipped_no_usable_reads: usize,
    pub rejected_read_tracks: usize,
    pub genes_with_rejected_reads: usize,
    pub errors: usize,
    pub failed_missing_inputs: usize,
    pub failed_processing: usize,
    pub failed_panics: usize,
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
    pub prepare: PrepareConfig,
    pub clustering: ClusteringConfig,
    pub counting: CountingConfig,
    pub runtime: RuntimeConfig,
    pub downsample: DownsampleConfig,
    pub emit_pooled_reads: bool,
    pub count_only: bool,
}

impl FullFlowOptions {
    /// Validate all scientific and runtime parameters before any output is created.
    pub fn validate(&self) -> Result<(), crate::config::ParameterError> {
        self.runtime.validate()?;
        self.prepare.validate()?;
        self.clustering.validate()?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct FullFlowResult {
    pub batch: BatchRunResult,
    pub isoform_bed: PathBuf,
    pub unused_bed: PathBuf,
    pub count_csv: PathBuf,
    pub desc_prefix: PathBuf,
    pub multi_sample: Option<MultiSampleOutputPaths>,
    /// Final unique mapping, present only when unique assignment was selected.
    pub unique_read_to_isoform_tsv: Option<PathBuf>,
}

fn remove_optional_output(path: &Path, kind: &str) -> anyhow::Result<()> {
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("remove stale {kind} output {path:?}")),
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum GeneOutcome {
    Processed,
    Skipped(GeneSkipReason),
    Failed(GeneFailure),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum GeneSkipReason {
    CompletedOutputs,
    EmptyReads,
    NoUsableReads,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum GeneFailureKind {
    MissingInputs,
    Processing,
    Panic,
}

impl GeneFailureKind {
    pub(super) fn as_str(self) -> &'static str {
        match self {
            Self::MissingInputs => "missing_inputs",
            Self::Processing => "processing",
            Self::Panic => "panic",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct GeneFailure {
    pub(super) kind: GeneFailureKind,
    pub(super) message: String,
}

fn seed_for_gene(base_seed: u64, gene: &str) -> u64 {
    base_seed ^ crate::rng::fnv1a64(gene.as_bytes())
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct DownsampleRecord {
    pub(super) gene: String,
    original_reads: usize,
    sampled_reads: usize,
    target_reads: usize,
    seed: u64,
    scale_factor: f64,
    input_fingerprint: Option<String>,
}

#[derive(Clone, Debug)]
pub(super) struct ProcessGeneResult {
    pub(super) outcome: GeneOutcome,
    pub(super) downsample: Option<DownsampleRecord>,
    pub(super) resume_reason: String,
    pub(super) rejected_read_tracks: usize,
    pub(super) all_reads_rejected: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ResumeDecision {
    pub(super) gene: String,
    pub(super) action: &'static str,
    pub(super) reason: String,
}

impl ProcessGeneResult {
    pub(super) fn failed(kind: GeneFailureKind, message: impl Into<String>) -> Self {
        Self {
            outcome: GeneOutcome::Failed(GeneFailure {
                kind,
                message: message.into(),
            }),
            downsample: None,
            resume_reason: "not_checked_due_to_failure".to_owned(),
            rejected_read_tracks: 0,
            all_reads_rejected: false,
        }
    }
}

fn should_downsample_gene(gene: &str, args: &BatchRunOptions) -> bool {
    args.downsample.selects(gene)
}

fn reservoir_sample_reads(
    path: &Path,
    target_reads: usize,
    seed: u64,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<(Vec<Transcript>, usize, Vec<RejectedReadRecord>)> {
    let mut rng = crate::rng::Lcg64::new(seed);
    // Avoid huge preallocation when downsampling is enabled by default. The vector will grow to
    // `target_reads` as needed, but most genes are far smaller than the cap.
    let mut sampled: Vec<Transcript> = Vec::with_capacity(target_reads.min(4096));
    let mut total_reads = 0usize;

    let mut reader = read_bed12(path).with_context(|| format!("open reads {path:?}"))?;
    loop {
        let next = match invalid_read_policy {
            InvalidReadPolicy::Skip => reader.next_recovering_read(),
            InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse reads {path:?}"))?;
        let Some(tx) = next else {
            break;
        };
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

    Ok((sampled, total_reads, reader.take_rejected_reads()))
}

fn read_gene_read_records(
    path: &Path,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<(Vec<Transcript>, Vec<RejectedReadRecord>)> {
    let mut reader = read_bed12(path).with_context(|| format!("open reads {path:?}"))?;
    let mut records = Vec::new();
    loop {
        let next = match invalid_read_policy {
            InvalidReadPolicy::Skip => reader.next_recovering_read(),
            InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse reads {path:?}"))?;
        let Some(record) = next else {
            break;
        };
        records.push(record);
    }
    Ok((records, reader.take_rejected_reads()))
}

fn fingerprint_file(path: &Path) -> anyhow::Result<String> {
    let file =
        fs::File::open(path).with_context(|| format!("open input for fingerprint {path:?}"))?;
    let mut reader = BufReader::new(file);
    let mut buffer = vec![0u8; 1024 * 1024];
    let mut hash = 14695981039346656037;
    loop {
        let read_len = reader
            .read(&mut buffer)
            .with_context(|| format!("fingerprint input {path:?}"))?;
        if read_len == 0 {
            break;
        }
        crate::rng::update_fnv1a64(&mut hash, &buffer[..read_len]);
    }
    Ok(format!("fnv1a64:{hash:016x}"))
}

impl DownsampleRecord {
    fn validate(&self, path: &Path) -> anyhow::Result<()> {
        if self.gene.trim().is_empty() || self.gene == "none" {
            anyhow::bail!("invalid empty gene in downsample metadata {path:?}");
        }
        if self.sampled_reads == 0 || self.original_reads < self.sampled_reads {
            anyhow::bail!(
                "invalid read counts for gene {:?} in {path:?}: expected original_reads >= sampled_reads > 0, got original_reads={} sampled_reads={}",
                self.gene,
                self.original_reads,
                self.sampled_reads
            );
        }
        if self.target_reads == 0 || self.sampled_reads > self.target_reads {
            anyhow::bail!(
                "invalid target_reads for gene {:?} in {path:?}: expected target_reads >= sampled_reads > 0, got target_reads={} sampled_reads={}",
                self.gene,
                self.target_reads,
                self.sampled_reads
            );
        }
        if !self.scale_factor.is_finite() || self.scale_factor <= 0.0 {
            anyhow::bail!(
                "invalid scale_factor for gene {:?} in {path:?}: expected a finite positive value, got {}",
                self.gene,
                self.scale_factor
            );
        }
        let expected_scale = self.original_reads as f64 / self.sampled_reads as f64;
        if (self.scale_factor - expected_scale).abs() > f64::EPSILON * expected_scale.max(1.0) * 4.0
        {
            anyhow::bail!(
                "inconsistent scale_factor for gene {:?} in {path:?}: expected {}, got {}",
                self.gene,
                expected_scale,
                self.scale_factor
            );
        }
        if let Some(fingerprint) = self.input_fingerprint.as_deref() {
            let Some(hex) = fingerprint.strip_prefix("fnv1a64:") else {
                anyhow::bail!(
                    "invalid input_fingerprint for gene {:?} in {path:?}: {fingerprint:?}",
                    self.gene
                );
            };
            if hex.len() != 16 || !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
                anyhow::bail!(
                    "invalid input_fingerprint for gene {:?} in {path:?}: {fingerprint:?}",
                    self.gene
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
fn write_downsample_records(path: &Path, records: &[DownsampleRecord]) -> anyhow::Result<()> {
    let mut file =
        fs::File::create(path).with_context(|| format!("write downsample info {path:?}"))?;
    write_downsample_records_to_writer(&mut file, path, records)
}

pub(super) fn write_downsample_records_to_writer<W: Write>(
    writer: &mut W,
    diagnostic_path: &Path,
    records: &[DownsampleRecord],
) -> anyhow::Result<()> {
    for record in records {
        record.validate(diagnostic_path)?;
    }
    writeln!(
        writer,
        "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads\tinput_fingerprint"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.gene,
            record.original_reads,
            record.sampled_reads,
            record.scale_factor,
            record.seed,
            record.target_reads,
            record.input_fingerprint.as_deref().unwrap_or("")
        )?;
    }
    writer
        .flush()
        .with_context(|| format!("flush downsample info {diagnostic_path:?}"))
}

fn read_downsample_records(path: &Path) -> anyhow::Result<Vec<DownsampleRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }

    let file = fs::File::open(path).with_context(|| format!("open downsample file {path:?}"))?;
    let reader = BufReader::new(file);
    let mut header: Option<Vec<String>> = None;
    let mut records = Vec::new();
    for (zero_based_line_no, line) in reader.lines().enumerate() {
        let line_no = zero_based_line_no + 1;
        let line = line.with_context(|| format!("read downsample file {path:?}"))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if header.is_none() {
            let columns: Vec<String> = line.split('\t').map(str::to_owned).collect();
            for required in [
                "gene",
                "original_reads",
                "sampled_reads",
                "scale_factor",
                "seed",
                "target_reads",
            ] {
                if !columns.iter().any(|column| column == required) {
                    anyhow::bail!(
                        "downsample metadata {path:?}:{line_no} is missing required column {required:?}"
                    );
                }
            }
            header = Some(columns);
            continue;
        }

        let columns = header.as_ref().expect("header initialized above");
        let fields: Vec<&str> = line.split('\t').collect();
        let field = |name: &str| -> anyhow::Result<&str> {
            let index = columns
                .iter()
                .position(|column| column == name)
                .expect("required header was validated");
            fields.get(index).copied().ok_or_else(|| {
                anyhow::anyhow!(
                    "downsample metadata {path:?}:{line_no} has no value for column {name:?}"
                )
            })
        };
        let parse_usize = |name: &str| -> anyhow::Result<usize> {
            let value = field(name)?;
            value.parse().with_context(|| {
                format!("parse {name} value {value:?} in downsample metadata {path:?}:{line_no}")
            })
        };
        let parse_u64 = |name: &str| -> anyhow::Result<u64> {
            let value = field(name)?;
            value.parse().with_context(|| {
                format!("parse {name} value {value:?} in downsample metadata {path:?}:{line_no}")
            })
        };
        let scale_factor_value = field("scale_factor")?;
        let input_fingerprint = columns
            .iter()
            .position(|column| column == "input_fingerprint")
            .and_then(|index| fields.get(index))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        let record = DownsampleRecord {
            gene: field("gene")?.trim().to_owned(),
            original_reads: parse_usize("original_reads")?,
            sampled_reads: parse_usize("sampled_reads")?,
            scale_factor: scale_factor_value.parse().with_context(|| {
                format!(
                    "parse scale_factor value {scale_factor_value:?} in downsample metadata {path:?}:{line_no}"
                )
            })?,
            seed: parse_u64("seed")?,
            target_reads: parse_usize("target_reads")?,
            input_fingerprint,
        };
        record.validate(path)?;
        records.push(record);
    }
    if header.is_none() {
        anyhow::bail!("downsample metadata {path:?} has no header");
    }
    Ok(records)
}

fn read_per_gene_downsample(
    path: &Path,
    expected_gene: &str,
) -> anyhow::Result<Option<DownsampleRecord>> {
    let mut records = read_downsample_records(path)?;
    if records.is_empty() {
        return Ok(None);
    }
    if records.len() != 1 {
        anyhow::bail!(
            "per-gene downsample metadata {path:?} must contain exactly one record, found {}",
            records.len()
        );
    }
    let record = records.pop().expect("one record checked above");
    if record.gene != expected_gene {
        anyhow::bail!(
            "per-gene downsample metadata {path:?} belongs to gene {:?}, expected {expected_gene:?}",
            record.gene
        );
    }
    Ok(Some(record))
}

fn read_per_gene_downsample_records(
    output_root: &Path,
    genes: &[String],
) -> anyhow::Result<Vec<DownsampleRecord>> {
    let mut records = Vec::new();
    for gene in genes {
        let gene_id = GeneId::parse(gene)?;
        let path = gene_dir_path(output_root, &gene_id)?.join("downsample.tsv");
        ensure_destination_within(output_root, &path)?;
        let metadata = fs::symlink_metadata(&path)
            .with_context(|| format!("inspect verified per-gene downsample metadata {path:?}"))?;
        if !metadata.is_file() || metadata.file_type().is_symlink() {
            anyhow::bail!("verified per-gene downsample metadata is not a regular file: {path:?}");
        }
        if let Some(record) = read_per_gene_downsample(&path, gene_id.as_str())
            .with_context(|| format!("read verified per-gene downsample metadata {path:?}"))?
        {
            records.push(record);
        }
    }
    records.sort_by(|left, right| left.gene.cmp(&right.gene));
    for pair in records.windows(2) {
        if pair[0].gene == pair[1].gene {
            anyhow::bail!(
                "duplicate verified per-gene downsample metadata for gene {:?}",
                pair[0].gene
            );
        }
    }
    Ok(records)
}

fn reject_cross_gene_downsampling(
    output_root: &Path,
    genes: &[String],
    downsample_records: &[DownsampleRecord],
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<()> {
    if downsample_records.is_empty() {
        return Ok(());
    }
    let downsampled_genes: HashSet<&str> = downsample_records
        .iter()
        .map(|record| record.gene.as_str())
        .collect();
    let mut first_gene_by_read: HashMap<String, String> = HashMap::new();

    for gene in genes {
        let gene_id = GeneId::parse(gene)?;
        let reads_path = gene_artifact_path(output_root, &gene_id, "_nano.bed")?;
        let mut reader = read_bed12(&reads_path)
            .with_context(|| format!("open prepared reads for gene {gene:?}: {reads_path:?}"))?;
        loop {
            let next = match invalid_read_policy {
                InvalidReadPolicy::Skip => reader.next_recovering_read(),
                InvalidReadPolicy::Fail => reader.next_strict_read(),
            }
            .with_context(|| format!("parse prepared reads for gene {gene:?}: {reads_path:?}"))?;
            let Some(read) = next else {
                break;
            };
            match first_gene_by_read.entry(read.name.clone()) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(gene.clone());
                }
                std::collections::hash_map::Entry::Occupied(entry)
                    if entry.get() != gene
                        && (downsampled_genes.contains(entry.get().as_str())
                            || downsampled_genes.contains(gene.as_str())) =>
                {
                    anyhow::bail!(
                        "independent per-gene downsampling is unsafe for molecule {:?}, which is assigned to genes {:?} and {:?}; rerun with --max-reads-per-gene 0 or exclude every multi-gene candidate from downsampling",
                        read.name,
                        entry.get(),
                        gene
                    );
                }
                std::collections::hash_map::Entry::Occupied(_) => {}
            }
        }
    }
    Ok(())
}

fn read_gene_list(path: &Path) -> anyhow::Result<Vec<String>> {
    let file = fs::File::open(path).with_context(|| format!("open gene list {path:?}"))?;
    let reader = BufReader::new(file);
    let mut genes: Vec<String> = Vec::new();
    for (zero_based, line) in reader.lines().enumerate() {
        let line_no = zero_based + 1;
        let line = line.with_context(|| format!("read gene list {path:?}:{line_no}"))?;
        let gene = line.trim();
        if gene.is_empty() || gene.starts_with('#') {
            continue;
        }
        GeneId::parse(gene)
            .with_context(|| format!("invalid gene-list entry {path:?}:{line_no}"))?;
        genes.push(gene.to_owned());
    }
    genes.sort();
    genes.dedup();
    validate_gene_ids(genes.iter().map(String::as_str))?;
    Ok(genes)
}

fn per_gene_run_request(
    gene: &GeneId,
    args: &BatchRunOptions,
    reads: &Path,
    reference: &Path,
) -> anyhow::Result<RunRequest> {
    let downsample_selected = should_downsample_gene(gene.as_str(), args);
    Ok(RunRequest {
        gene: gene.as_str().to_owned(),
        inputs: vec![
            InputArtifact::from_file("reads", reads)?,
            InputArtifact::from_file("reference", reference)?,
        ],
        options: EffectiveOptions {
            cluster_mode: args.cluster_mode.as_str().to_owned(),
            prepare_fraction_read: args.prepare.fraction_read,
            prepare_fraction_ref: args.prepare.fraction_ref,
            sw_score: args.clustering.sw_score,
            batch_size: args.clustering.batch_size,
            batch_rounds: args.clustering.batch_rounds,
            name2_mode: args.clustering.name2_mode.as_str().to_owned(),
            platform_preset: args.clustering.junction.platform_preset.as_str().to_owned(),
            junction_correction_offset: args.clustering.junction.correction.offset,
            junction_correction_min_support: args.clustering.junction.correction.min_support,
            sl_partial_five_prime_offset: args.clustering.junction.sl.partial_five_prime_end_offset,
            sl_same_junction_five_prime_offset: args
                .clustering
                .junction
                .sl
                .same_junction_five_prime_end_offset,
            sl_five_prime_cluster_offset: args.clustering.junction.sl.five_prime_cluster_offset,
            sl_min_five_prime_cluster_support: args
                .clustering
                .junction
                .sl
                .min_five_prime_cluster_support,
            same_junction_three_prime_offset: args
                .clustering
                .junction
                .three_prime
                .same_junction_three_prime_end_offset,
            three_prime_cluster_offset: args
                .clustering
                .junction
                .three_prime
                .three_prime_cluster_offset,
            min_three_prime_cluster_support: args
                .clustering
                .junction
                .three_prime
                .min_three_prime_cluster_support,
            overlap_cutoff1: args.clustering.overlap.cutoff1,
            overlap_cutoff2: args.clustering.overlap.cutoff2,
            overlap_intron_weight: args.clustering.overlap.intron_weight,
            assignment_mode: args.counting.assignment_mode.to_string(),
            unique_assignment_junction_offset: args.counting.unique_assignment.junction_offset,
            invalid_read_policy: args.runtime.invalid_read_policy.to_string(),
            downsample_selected,
            max_reads_per_gene: if downsample_selected {
                args.downsample.max_reads_per_gene
            } else {
                0
            },
        },
        tool: ToolIdentity::current(),
        seed: if downsample_selected {
            seed_for_gene(args.downsample.seed, gene.as_str())
        } else {
            0
        },
    })
}

fn per_gene_output_specs(
    out_isoforms: &Path,
    out_unused: &Path,
    out_mapping: &Path,
    downsample: &Path,
    rejected_reads: &Path,
) -> Vec<OutputSpec> {
    vec![
        OutputSpec::new("isoforms", out_isoforms, RecordCountKind::NonEmptyLines),
        OutputSpec::new("unused", out_unused, RecordCountKind::NonEmptyLines),
        OutputSpec::new(
            "read_to_isoform",
            out_mapping,
            RecordCountKind::NonEmptyLines,
        ),
        OutputSpec::new(
            "downsample",
            downsample,
            RecordCountKind::HeaderThenNonEmptyLines,
        ),
        OutputSpec::new(
            "rejected_reads",
            rejected_reads,
            RecordCountKind::HeaderThenNonEmptyLines,
        ),
    ]
}

fn validate_count_only_gene_manifests(
    output_root: &Path,
    genes: &[String],
    cluster_mode: ClusterMode,
) -> anyhow::Result<()> {
    for gene in genes {
        let gene_id = GeneId::parse(gene)?;
        let reads = gene_artifact_path(output_root, &gene_id, "_nano.bed")?;
        let reference = gene_artifact_path(output_root, &gene_id, "_gff.bed")?;
        let isoforms = gene_artifact_path(
            output_root,
            &gene_id,
            cluster_mode.per_gene_isoform_suffix(),
        )?;
        let unused = gene_artifact_path(output_root, &gene_id, "_unused.bed")?;
        let mapping = gene_artifact_path(output_root, &gene_id, "_read_to_isoform.tsv")?;
        let gene_dir = gene_dir_path(output_root, &gene_id)?;
        let downsample = gene_dir.join("downsample.tsv");
        let rejected_reads = gene_dir.join("rejected_reads.tsv");
        let manifest = gene_dir.join(MANIFEST_FILE_NAME);
        let expected_inputs = vec![
            InputArtifact::from_file("reads", &reads),
            InputArtifact::from_file("reference", &reference),
        ]
        .into_iter()
        .collect::<anyhow::Result<Vec<_>>>()
        .with_context(|| {
            format!(
                "validate current prepared inputs for count-only gene {:?}",
                gene_id.as_str()
            )
        })?;
        let output_specs =
            per_gene_output_specs(&isoforms, &unused, &mapping, &downsample, &rejected_reads);
        validate_recorded_completion(
            &manifest,
            gene_id.as_str(),
            cluster_mode.as_str(),
            &expected_inputs,
            &output_specs,
        )
        .with_context(|| {
            format!(
                "count-only refuses stale or unverified results for gene {:?}; rerun flow without --count-only",
                gene_id.as_str()
            )
        })?;
    }
    Ok(())
}

fn process_gene(gene: &str, args: &BatchRunOptions) -> ProcessGeneResult {
    let gene_id = match GeneId::parse(gene) {
        Ok(gene) => gene,
        Err(error) => {
            return ProcessGeneResult::failed(
                GeneFailureKind::Processing,
                format!("invalid biological gene ID {gene:?}: {error:#}"),
            );
        }
    };
    let reads = match gene_artifact_path(&args.input_root, &gene_id, "_nano.bed") {
        Ok(path) => path,
        Err(error) => {
            return ProcessGeneResult::failed(GeneFailureKind::Processing, format!("{error:#}"));
        }
    };
    let reference = match gene_artifact_path(&args.input_root, &gene_id, "_gff.bed") {
        Ok(path) => path,
        Err(error) => {
            return ProcessGeneResult::failed(GeneFailureKind::Processing, format!("{error:#}"));
        }
    };
    let missing: Vec<PathBuf> = [&reads, &reference]
        .into_iter()
        .filter(|path| !path.exists())
        .cloned()
        .collect();
    if !missing.is_empty() {
        return ProcessGeneResult::failed(
            GeneFailureKind::MissingInputs,
            format!("missing required per-gene input(s): {missing:?}"),
        );
    }

    match process_gene_inputs(&gene_id, args, &reads, &reference) {
        Ok(result) => result,
        Err(error) => ProcessGeneResult::failed(GeneFailureKind::Processing, format!("{error:#}")),
    }
}

fn publish_empty_gene_result(
    manifest_path: &Path,
    request: RunRequest,
    output_specs: &[OutputSpec],
    rejected_reads: &[RejectedReadRecord],
    skip_reason: GeneSkipReason,
    resume_reason: String,
) -> anyhow::Result<ProcessGeneResult> {
    invalidate_completion_manifest(manifest_path)?;
    for spec in output_specs {
        match fs::symlink_metadata(&spec.path) {
            Ok(metadata) if metadata.is_file() || metadata.file_type().is_symlink() => {
                fs::remove_file(&spec.path)
                    .with_context(|| format!("remove stale empty-gene artifact {:?}", spec.path))?;
            }
            Ok(_) => {
                anyhow::bail!(
                    "refusing to remove non-file stale empty-gene artifact {:?}",
                    spec.path
                )
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("inspect stale empty-gene artifact {:?}", spec.path));
            }
        }
    }

    for spec in output_specs {
        atomic_write_with(&spec.path, |file| match spec.role {
            "downsample" => write_downsample_records_to_writer(file, &spec.path, &[]),
            "rejected_reads" => write_rejected_reads_tsv_to_writer(file, rejected_reads)
                .map_err(anyhow::Error::from),
            _ => Ok(()),
        })?;
    }
    let manifest = RunManifest::complete(request, output_specs)
        .context("build empty-gene completion manifest")?;
    write_completion_manifest(manifest_path, &manifest)
        .context("publish empty-gene completion manifest")?;

    Ok(ProcessGeneResult {
        outcome: GeneOutcome::Skipped(skip_reason),
        downsample: None,
        resume_reason,
        rejected_read_tracks: rejected_reads.len(),
        all_reads_rejected: skip_reason == GeneSkipReason::NoUsableReads,
    })
}

fn process_gene_inputs(
    gene: &GeneId,
    args: &BatchRunOptions,
    reads: &Path,
    reference: &Path,
) -> anyhow::Result<ProcessGeneResult> {
    let out_dir = gene_dir_path(&args.output_root, gene)?;
    fs::create_dir_all(&out_dir).with_context(|| format!("create {out_dir:?}"))?;
    ensure_destination_within(&args.output_root, &out_dir)?;
    write_gene_id_marker(&out_dir, gene)?;

    let out_isoforms = gene_artifact_path(
        &args.output_root,
        gene,
        args.cluster_mode.per_gene_isoform_suffix(),
    )?;
    let out_unused = gene_artifact_path(&args.output_root, gene, "_unused.bed")?;
    let out_mapping = gene_artifact_path(&args.output_root, gene, "_read_to_isoform.tsv")?;
    let per_gene_downsample_path = out_dir.join("downsample.tsv");
    let rejected_reads_path = out_dir.join("rejected_reads.tsv");
    ensure_destination_within(&args.output_root, &per_gene_downsample_path)?;
    ensure_destination_within(&args.output_root, &rejected_reads_path)?;
    let manifest_path = out_dir.join(MANIFEST_FILE_NAME);
    ensure_destination_within(&args.output_root, &manifest_path)?;

    let request = per_gene_run_request(gene, args, reads, reference)?;
    let output_specs = per_gene_output_specs(
        &out_isoforms,
        &out_unused,
        &out_mapping,
        &per_gene_downsample_path,
        &rejected_reads_path,
    );
    let cache_decision = if args.runtime.force {
        CacheDecision::Rebuild("forced".to_owned())
    } else {
        assess_cache(&manifest_path, &request, &output_specs)
    };

    if cache_decision == CacheDecision::Reuse {
        let downsample = read_per_gene_downsample(&per_gene_downsample_path, gene.as_str())
            .with_context(|| {
                format!(
                    "load persisted downsampling state for gene {:?}",
                    gene.as_str()
                )
            })?;
        let rejected_read_tracks = count_rejected_reads_tsv(&rejected_reads_path)
            .with_context(|| format!("load rejected-read state {rejected_reads_path:?}"))?;
        let all_reads_rejected = rejected_read_tracks > 0
            && [&out_isoforms, &out_unused, &out_mapping]
                .into_iter()
                .all(|path| fs::metadata(path).is_ok_and(|metadata| metadata.len() == 0));
        return Ok(ProcessGeneResult {
            outcome: GeneOutcome::Skipped(GeneSkipReason::CompletedOutputs),
            downsample,
            resume_reason: cache_decision.reason().to_owned(),
            rejected_read_tracks,
            all_reads_rejected,
        });
    }
    let resume_reason = cache_decision.reason().to_owned();
    invalidate_completion_manifest(&manifest_path)?;

    let (read_records, downsample, rejected_reads) = if should_downsample_gene(gene.as_str(), args)
    {
        let target_reads = args.downsample.max_reads_per_gene.max(1);
        let seed = seed_for_gene(args.downsample.seed, gene.as_str());
        let input_fingerprint = fingerprint_file(reads)?;
        let (sampled, total_reads, rejected_reads) =
            reservoir_sample_reads(reads, target_reads, seed, args.runtime.invalid_read_policy)?;
        let sampled_reads = sampled.len();
        let downsample = if total_reads > sampled_reads && sampled_reads > 0 {
            Some(DownsampleRecord {
                gene: gene.as_str().to_owned(),
                original_reads: total_reads,
                sampled_reads,
                target_reads,
                seed,
                scale_factor: total_reads as f64 / sampled_reads as f64,
                input_fingerprint: Some(input_fingerprint),
            })
        } else {
            None
        };
        if let Some(record) = downsample.as_ref() {
            eprintln!(
                "{}: subsample gene={} original_reads={} sampled_reads={} scale_factor={:.6} seed={}",
                args.cluster_mode.batch_log_label(),
                record.gene,
                record.original_reads,
                record.sampled_reads,
                record.scale_factor,
                record.seed
            );
        }
        (sampled, downsample, rejected_reads)
    } else {
        let (records, rejected_reads) =
            read_gene_read_records(reads, args.runtime.invalid_read_policy)?;
        (records, None, rejected_reads)
    };

    if read_records.is_empty() {
        let skip_reason = if rejected_reads.is_empty() {
            GeneSkipReason::EmptyReads
        } else {
            GeneSkipReason::NoUsableReads
        };
        let reason = if rejected_reads.is_empty() {
            "empty_reads_stale_outputs_removed"
        } else {
            "all_read_tracks_rejected"
        };
        return publish_empty_gene_result(
            &manifest_path,
            request,
            &output_specs,
            &rejected_reads,
            skip_reason,
            reason.to_owned(),
        );
    }

    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(reference)
        .with_context(|| format!("open reference {reference:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reference {reference:?}"))?;

    let result = args
        .clustering
        .cluster_gene(args.cluster_mode, &read_records, &refs, 1)?;

    atomic_write_with(&out_isoforms, |file| {
        crate::cluster::output::write_isoforms_bed_to_writer(file, &result.isoforms)
            .map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write {out_isoforms:?}"))?;
    atomic_write_with(&out_unused, |file| {
        crate::cluster::output::write_isoforms_bed_to_writer(file, &result.unused)
            .map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write {out_unused:?}"))?;
    atomic_write_with(&out_mapping, |file| {
        crate::cluster::output::write_read_to_isoform_tsv_writer(file, &result.read_to_isoform)
            .map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write {out_mapping:?}"))?;
    atomic_write_with(&per_gene_downsample_path, |file| {
        write_downsample_records_to_writer(file, &per_gene_downsample_path, downsample.as_slice())
    })
    .with_context(|| {
        format!("atomically write per-gene downsample info {per_gene_downsample_path:?}")
    })?;
    atomic_write_with(&rejected_reads_path, |file| {
        write_rejected_reads_tsv_to_writer(file, &rejected_reads).map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write rejected reads {rejected_reads_path:?}"))?;

    let request_after_processing = per_gene_run_request(gene, args, reads, reference)?;
    if request_after_processing != request {
        anyhow::bail!(
            "per-gene inputs or effective options changed while processing gene {:?}; outputs were left without a completion manifest and will be rebuilt",
            gene.as_str()
        );
    }
    let manifest = RunManifest::complete(request, &output_specs)
        .with_context(|| format!("build completion manifest {manifest_path:?}"))?;
    write_completion_manifest(&manifest_path, &manifest)
        .with_context(|| format!("publish completion manifest {manifest_path:?}"))?;

    Ok(ProcessGeneResult {
        outcome: GeneOutcome::Processed,
        downsample,
        resume_reason,
        rejected_read_tracks: rejected_reads.len(),
        all_reads_rejected: false,
    })
}

pub fn run_clusterj_batch(args: BatchRunOptions) -> anyhow::Result<BatchRunResult> {
    args.validate().context("invalid batch options")?;
    if let Some(prefix) = args.prepare_prefix.as_deref() {
        SafePathComponent::parse("output prefix", prefix)?;
    }
    args.downsample.validate()?;
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
        let res = prepare_dir_from_paths_with_policy(
            reads,
            reference,
            &args.input_root,
            prefix,
            AddGeneOpts {
                fraction_read: args.prepare.fraction_read,
                fraction_ref: args.prepare.fraction_ref,
            },
            args.runtime.invalid_read_policy,
        )?;
        eprintln!(
            "prepare: genes={}, dedup_reads={}, novel_reads={}, rejected_read_tracks={}",
            res.genes.len(),
            res.dedup_reads,
            res.novel_reads,
            res.rejected_read_tracks
        );
        Some(res)
    } else {
        None
    };

    let mut genes = match (&args.gene_list, prepared.as_ref()) {
        (Some(list), Some(prepared)) => {
            let selected = read_gene_list(list)?;
            let prepared_genes: HashSet<&str> = prepared.genes.iter().map(String::as_str).collect();
            let stale: Vec<&str> = selected
                .iter()
                .map(String::as_str)
                .filter(|gene| !prepared_genes.contains(gene))
                .collect();
            if !stale.is_empty() {
                anyhow::bail!(
                    "--gene-list selects gene(s) absent from the just-prepared generation: {stale:?}"
                );
            }
            selected
        }
        (None, Some(prepared)) => prepared.genes.clone(),
        (Some(list), None) => read_gene_list(list)?,
        (None, None) => {
            anyhow::bail!(
                "--gene-list is required when clusterj_batch is run without inline --prepare-reads/--prepare-reference/--prepare-prefix; directory discovery is intentionally disabled to prevent stale-gene reuse"
            )
        }
    };
    genes.sort();
    genes.dedup();
    let gene_ids = validate_gene_ids(genes.iter().map(String::as_str))
        .context("validate selected biological gene IDs")?;
    validate_pipeline_gene_namespace(&gene_ids, None)
        .context("validate selected gene paths against batch-level artifacts")?;

    let total = genes.len();
    if total == 0 {
        anyhow::bail!("no genes found (input-root={:?})", args.input_root);
    }

    let gene_path_map = args
        .output_root
        .join(format!("{batch_file_prefix}_gene_paths.tsv"));
    ensure_destination_within(&args.output_root, &gene_path_map)?;
    write_gene_path_map(&gene_path_map, &gene_ids)?;

    eprintln!(
        "{batch_log_label}: {} genes, {} worker threads",
        total, args.runtime.threads
    );

    if args.downsample.max_reads_per_gene > 0 {
        let scope = if args.downsample.genes.is_empty() {
            "all genes".to_owned()
        } else {
            format!("{} selected gene(s)", args.downsample.genes.len())
        };
        eprintln!(
            "{batch_log_label}: per-gene subsampling enabled: cap={} scope={scope} seed={}",
            args.downsample.max_reads_per_gene, args.downsample.seed
        );
    }

    if args.runtime.threads > 1 && args.downsample.max_reads_per_gene == 0 {
        eprintln!(
            "{batch_log_label}: note: --threads > 1 with --max-reads-per-gene=0 can use a lot of memory on large genes; \
consider --max-reads-per-gene and/or --name2-mode coverage"
        );
    }
    if args.cluster_mode == ClusterMode::Cluster && args.clustering.batch_size == 0 {
        eprintln!(
            "{batch_log_label}: note: overlap mode batching is disabled because --batch-size=0; \
each gene will run one full two-pass overlap merge"
        );
    }

    let mut execution = crate::flow::executor::execute_genes(genes, &args, process_gene);
    execution.prepare_rejected_read_tracks = args.prepare_rejected_read_tracks.saturating_add(
        prepared
            .as_ref()
            .map_or(0, |result| result.rejected_read_tracks),
    );
    let elapsed = execution.elapsed;
    let processed = execution.processed;
    let skipped = execution.skipped;
    let skipped_completed_outputs = execution.skipped_completed_outputs;
    let skipped_empty_reads = execution.skipped_empty_reads;
    let skipped_no_usable_reads = execution.skipped_no_usable_reads;
    let rejected_read_tracks = execution
        .prepare_rejected_read_tracks
        .saturating_add(execution.rejected_read_tracks);
    let genes_with_rejected_reads = execution.genes_with_rejected_reads;
    let errors = execution.errors;
    let mergeable_genes = execution.mergeable_genes.clone();
    let infrastructure_errors = execution.infrastructure_error_count();
    let failed_missing_inputs = execution.failed_missing_inputs;
    let failed_processing = execution.failed_processing;
    let failed_panics = execution.failed_panics;
    eprintln!(
        "done: processed={processed}, skipped={skipped}, errors={errors}, elapsed={:?}",
        elapsed
    );
    if rejected_read_tracks > 0 {
        eprintln!(
            "{batch_log_label}: warning: excluded {rejected_read_tracks} malformed read track(s); see *_rejected_reads.tsv diagnostics"
        );
    }

    let published = crate::flow::reporting::publish_batch_report(
        crate::flow::reporting::BatchReportContext {
            args: &args,
            batch_file_prefix,
            gene_path_map: &gene_path_map,
            total,
        },
        &mut execution,
    )?;
    if infrastructure_errors > 0 {
        let first_error = published
            .first_error
            .as_deref()
            .unwrap_or("<error details unavailable>");
        anyhow::bail!(
            "{batch_log_label} failed: {infrastructure_errors} worker/reporting infrastructure error(s); downstream stages are unsafe; details: {:?}; first failure: {first_error}",
            published.error_path
        );
    }
    if errors > 0 && mergeable_genes.is_empty() {
        let first_error = published
            .first_error
            .as_deref()
            .unwrap_or("<error details unavailable>");
        anyhow::bail!(
            "{batch_log_label} failed: all {errors} attempted gene(s) failed; no verified gene result is available for downstream stages; details: {:?}; first failure: {first_error}",
            published.error_path
        );
    }
    if errors > 0 && args.runtime.gene_error_policy == crate::flow::config::GeneErrorPolicy::Strict
    {
        let first_error = published
            .first_error
            .as_deref()
            .unwrap_or("<error details unavailable>");
        anyhow::bail!(
            "{batch_log_label} failed in strict mode: {errors} required gene(s) failed; no merge/count/description stages were run; details: {:?}; first failure: {first_error}",
            published.error_path
        );
    }
    if errors > 0 {
        eprintln!(
            "{batch_log_label}: warning: {errors} gene(s) failed and were excluded; continuing with {} verified gene(s); details: {:?}",
            mergeable_genes.len(),
            published.error_path
        );
    }
    let summary_path = published.summary_path;
    let error_path = published.error_path;
    let downsample_path = published.downsample_path;
    Ok(BatchRunResult {
        prepared,
        total_genes: total,
        mergeable_genes,
        processed,
        skipped,
        skipped_completed_outputs,
        skipped_empty_reads,
        skipped_no_usable_reads,
        rejected_read_tracks,
        genes_with_rejected_reads,
        errors,
        failed_missing_inputs,
        failed_processing,
        failed_panics,
        elapsed_seconds: elapsed.as_secs_f64(),
        summary_path,
        error_path,
        downsample_path,
    })
}

pub fn run_full_flow(opts: FullFlowOptions) -> anyhow::Result<FullFlowResult> {
    opts.validate().context("invalid flow options")?;
    if opts.reads.is_some() && opts.manifest.is_some() {
        anyhow::bail!("flow: use either reads or manifest, not both");
    }
    if !opts.count_only && opts.reads.is_none() && opts.manifest.is_none() {
        anyhow::bail!("flow: reads or manifest is required");
    }
    if opts.count_only && opts.runtime.force {
        anyhow::bail!(
            "flow: --force rebuilds per-gene results and cannot be used with --count-only"
        );
    }
    if opts.count_only
        && opts.runtime.gene_error_policy == crate::flow::config::GeneErrorPolicy::Strict
    {
        anyhow::bail!(
            "flow: strict gene-error handling applies to per-gene execution and cannot be used with count-only mode"
        );
    }
    let prefix = SafePathComponent::parse("output prefix", &opts.prefix)?;
    opts.downsample.validate()?;
    let mut sample_rows = opts
        .manifest
        .as_deref()
        .map(read_manifest_tsv)
        .transpose()?;
    {
        let mut external_inputs: Vec<(&str, &Path)> = vec![("reference input", &opts.reference)];
        if let Some(reads) = opts.reads.as_deref() {
            external_inputs.push(("reads input", reads));
        }
        if let Some(manifest) = opts.manifest.as_deref() {
            external_inputs.push(("manifest input", manifest));
        }
        if let Some(rows) = sample_rows.as_ref() {
            external_inputs.extend(
                rows.iter()
                    .map(|row| ("sample reads input", row.reads.as_path())),
            );
        }
        reject_external_inputs_in_output_root(&opts.output_root, external_inputs)?;
    }
    fs::create_dir_all(&opts.output_root)
        .with_context(|| format!("create {:?}", opts.output_root))?;

    let layout = FlowArtifactLayout::new(&opts.output_root, &prefix)?;
    let gene_list = layout.gene_list.clone();
    let mut count_only_genes: Option<Vec<String>> = None;

    let batch = if opts.count_only {
        let batch_file_prefix = opts.cluster_mode.batch_file_prefix();
        let batch_gene_path_map = opts
            .output_root
            .join(format!("{batch_file_prefix}_gene_paths.tsv"));
        let list_genes = gene_list
            .exists()
            .then(|| read_gene_list(&gene_list).with_context(|| format!("read {gene_list:?}")))
            .transpose()?;
        if list_genes.as_ref().is_some_and(Vec::is_empty) {
            anyhow::bail!(
                "flow: count-only selected no genes from the prefix-scoped gene list {:?}; rerun flow without --count-only",
                gene_list
            );
        }
        let map_genes = layout
            .gene_path_map
            .exists()
            .then(|| -> anyhow::Result<Vec<String>> {
                let genes = read_gene_path_map(&layout.gene_path_map).with_context(|| {
                    format!(
                        "read prefix-scoped gene path mapping {:?}",
                        layout.gene_path_map
                    )
                })?;
                Ok(genes
                    .into_iter()
                    .map(|gene| gene.as_str().to_owned())
                    .collect())
            })
            .transpose()?;
        let mut genes = match (list_genes, map_genes) {
            (Some(list), Some(mut map)) => {
                map.sort();
                map.dedup();
                if list != map {
                    anyhow::bail!(
                        "flow: count-only prefix-scoped gene metadata disagrees between {:?} and {:?}; rerun flow without --count-only",
                        gene_list,
                        layout.gene_path_map
                    );
                }
                list
            }
            (Some(list), None) => list,
            (None, Some(map)) => {
                eprintln!(
                    "flow: count-only gene list {:?} missing; recovering the prefix-scoped gene set from {:?}",
                    gene_list, layout.gene_path_map
                );
                map
            }
            (None, None) => {
                if !batch_gene_path_map.exists() {
                    anyhow::bail!(
                        "flow: count-only requires the prefix-scoped gene list {:?}, prefix-scoped gene path mapping {:?}, or cluster-batch gene path mapping {:?}; rerun flow or clusterj_batch without --count-only",
                        gene_list,
                        layout.gene_path_map,
                        batch_gene_path_map
                    );
                }
                let map = read_gene_path_map(&batch_gene_path_map)
                    .with_context(|| {
                        format!(
                            "read cluster-batch gene path mapping {:?}",
                            batch_gene_path_map
                        )
                    })?
                    .into_iter()
                    .map(|gene| gene.as_str().to_owned())
                    .collect();
                eprintln!(
                    "flow: count-only prefix-scoped gene metadata missing; recovering the gene set from {:?}",
                    batch_gene_path_map
                );
                map
            }
        };
        genes.sort();
        genes.dedup();
        if genes.is_empty() {
            anyhow::bail!(
                "flow: count-only selected no genes from the prefix-scoped metadata for {:?}; rerun flow without --count-only",
                opts.prefix
            );
        }
        let batch = BatchRunResult {
            prepared: None,
            total_genes: genes.len(),
            mergeable_genes: genes.clone(),
            processed: 0,
            skipped: genes.len(),
            skipped_completed_outputs: genes.len(),
            skipped_empty_reads: 0,
            skipped_no_usable_reads: 0,
            rejected_read_tracks: 0,
            genes_with_rejected_reads: 0,
            errors: 0,
            failed_missing_inputs: 0,
            failed_processing: 0,
            failed_panics: 0,
            elapsed_seconds: 0.0,
            summary_path: opts
                .output_root
                .join(format!("{batch_file_prefix}_summary.txt")),
            error_path: opts
                .output_root
                .join(format!("{batch_file_prefix}_errors.txt")),
            downsample_path: opts
                .output_root
                .join(format!("{batch_file_prefix}_downsample.tsv")),
        };
        eprintln!(
            "flow: count-only: reusing existing per-gene outputs (genes={})",
            genes.len()
        );
        count_only_genes = Some(genes);
        batch
    } else {
        match (&opts.reads, &opts.manifest) {
            (Some(reads), None) => run_clusterj_batch(BatchRunOptions {
                cluster_mode: opts.cluster_mode,
                prepare_reads: Some(reads.clone()),
                prepare_reference: Some(opts.reference.clone()),
                prepare_prefix: Some(opts.prefix.clone()),
                prepare: opts.prepare,
                prepare_rejected_read_tracks: 0,
                input_root: opts.output_root.clone(),
                gene_list: Some(gene_list.clone()),
                output_root: opts.output_root.clone(),
                clustering: opts.clustering,
                counting: opts.counting,
                runtime: opts.runtime,
                downsample: opts.downsample.clone(),
            })?,
            (None, Some(manifest_path)) => {
                let rows = sample_rows
                    .take()
                    .context("flow manifest rows were not preloaded")?;
                let pooled_reads_out = if opts.emit_pooled_reads {
                    Some(
                        opts.output_root
                            .join(format!("{}_pooled_reads.bed", opts.prefix)),
                    )
                } else {
                    None
                };
                let prepared = prepare_dir_from_manifest_rows_with_policy(
                    &rows,
                    &opts.reference,
                    &opts.output_root,
                    &opts.prefix,
                    AddGeneOpts {
                        fraction_read: opts.prepare.fraction_read,
                        fraction_ref: opts.prepare.fraction_ref,
                    },
                    pooled_reads_out.as_deref(),
                    opts.runtime.invalid_read_policy,
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
                    "prepare: genes={}, dedup_reads={}, novel_reads={}, rejected_read_tracks={}",
                    prepared.genes.len(),
                    prepared.dedup_reads,
                    prepared.novel_reads,
                    prepared.rejected_read_tracks
                );

                let mut batch = run_clusterj_batch(BatchRunOptions {
                    cluster_mode: opts.cluster_mode,
                    prepare_reads: None,
                    prepare_reference: None,
                    prepare_prefix: None,
                    prepare: opts.prepare,
                    prepare_rejected_read_tracks: prepared.rejected_read_tracks,
                    input_root: opts.output_root.clone(),
                    gene_list: Some(gene_list.clone()),
                    output_root: opts.output_root.clone(),
                    clustering: opts.clustering,
                    counting: opts.counting,
                    runtime: opts.runtime,
                    downsample: opts.downsample.clone(),
                })?;
                batch.prepared = Some(prepared);
                sample_rows = Some(rows);
                batch
            }
            (Some(_), Some(_)) | (None, None) => unreachable!("validated flow input mode"),
        }
    };

    let selected_genes = if let Some(genes) = count_only_genes {
        genes
    } else {
        read_gene_list(&gene_list).with_context(|| format!("read {:?}", gene_list))?
    };
    if opts.count_only {
        preflight_gene_artifacts(
            &opts.output_root,
            &selected_genes,
            GeneArtifactRequirements {
                isoform_suffix: opts.cluster_mode.per_gene_isoform_suffix(),
                require_reads: opts.counting.assignment_mode
                    == crate::count::AssignmentMode::Unique,
            },
        )?;
        validate_count_only_gene_manifests(&opts.output_root, &selected_genes, opts.cluster_mode)?;
    }
    let gene_ids = validate_gene_ids(selected_genes.iter().map(String::as_str))?;
    validate_pipeline_gene_namespace(&gene_ids, Some(&prefix))
        .context("validate selected gene paths against flow-level artifacts")?;
    write_gene_path_map(&layout.gene_path_map, &gene_ids)?;
    let merge_genes = if opts.count_only {
        selected_genes.clone()
    } else {
        batch.mergeable_genes.clone()
    };
    let downsample_path = batch.downsample_path.clone();
    ensure_destination_within(&opts.output_root, &downsample_path)?;
    let downsample_records = read_per_gene_downsample_records(&opts.output_root, &merge_genes)
        .context("read scale factors from verified per-gene downsample metadata")?;
    reject_cross_gene_downsampling(
        &opts.output_root,
        &merge_genes,
        &downsample_records,
        opts.runtime.invalid_read_policy,
    )?;
    if opts.count_only {
        if downsample_records.is_empty() {
            remove_optional_output(&downsample_path, "aggregate downsample state")?;
        } else {
            atomic_write_with(&downsample_path, |file| {
                write_downsample_records_to_writer(file, &downsample_path, &downsample_records)
            })
            .with_context(|| format!("rebuild aggregate downsample state {downsample_path:?}"))?;
        }
    }
    let downsample_scales = if downsample_records.is_empty() {
        None
    } else {
        eprintln!(
            "flow: applying downsample scale factors (genes={}) from verified per-gene metadata",
            downsample_records.len()
        );
        Some(
            downsample_records
                .into_iter()
                .map(|record| (record.gene, record.scale_factor))
                .collect::<HashMap<_, _>>(),
        )
    };
    let isoform_bed = layout.isoform_bed.clone();
    let unused_bed = layout.unused_bed.clone();
    let read_to_isoform_tsv = layout.read_to_isoform_tsv.clone();
    let unique_read_to_isoform_tsv = layout.unique_read_to_isoform_tsv.clone();
    let unique_assignment_provenance_tsv = layout.unique_assignment_provenance_tsv.clone();
    let count_csv = layout.count_csv.clone();
    let desc_prefix = layout.desc_prefix.clone();

    eprintln!("flow: merge isoforms -> {:?}", isoform_bed);
    merge_gene_isoform_outputs(
        &opts.output_root,
        &merge_genes,
        opts.cluster_mode.per_gene_isoform_suffix(),
        &isoform_bed,
    )?;

    eprintln!("flow: merge unused -> {:?}", unused_bed);
    merge_gene_outputs(&opts.output_root, &merge_genes, "_unused.bed", &unused_bed)?;

    eprintln!("flow: merge read-to-isoform -> {:?}", read_to_isoform_tsv);
    merge_gene_outputs(
        &opts.output_root,
        &merge_genes,
        "_read_to_isoform.tsv",
        &read_to_isoform_tsv,
    )?;

    let isoforms = read_bed12_records(&isoform_bed, "isoform")?;
    let refs = read_bed12_records(&opts.reference, "reference")?;

    let merged_read_to_isoform;
    let selected_read_to_isoform;
    let count_read_to_isoform = if opts.counting.assignment_mode
        == crate::count::AssignmentMode::Unique
    {
        eprintln!("flow: unique assignment from per-gene folders");
        selected_read_to_isoform = select_unique_read_to_isoform_by_gene(
            &opts.output_root,
            &merge_genes,
            opts.cluster_mode,
            opts.counting.unique_assignment,
            opts.runtime.invalid_read_policy,
        )?;
        eprintln!(
            "flow: write unique read-to-isoform -> {:?}",
            unique_read_to_isoform_tsv
        );
        atomic_write_with(&unique_read_to_isoform_tsv, |file| {
            crate::cluster::output::write_read_to_isoform_tsv_writer(
                file,
                &selected_read_to_isoform,
            )
            .map_err(anyhow::Error::from)
        })
        .with_context(|| {
            format!("atomically write unique read-to-isoform {unique_read_to_isoform_tsv:?}")
        })?;
        atomic_write_with(&unique_assignment_provenance_tsv, |file| {
            crate::count::write_unique_assignment_provenance_to_writer(
                file,
                opts.counting.unique_assignment,
            )
            .map_err(anyhow::Error::from)
        })
        .with_context(|| {
            format!(
                "atomically write unique-assignment provenance {unique_assignment_provenance_tsv:?}"
            )
        })?;
        &selected_read_to_isoform
    } else {
        merged_read_to_isoform = crate::count::read_read_to_isoform_tsv(&read_to_isoform_tsv)
            .with_context(|| format!("read merged read_to_isoform {read_to_isoform_tsv:?}"))?;
        &merged_read_to_isoform
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
        Some(run_count_multi_atomic(
            rows,
            &isoforms,
            count_read_to_isoform,
            &output_prefix,
            downsample_scales.as_ref(),
        )?)
    } else {
        None
    };

    if let Some(multi) = multi_sample.as_ref() {
        atomic_copy(&multi.count_csv, &count_csv).with_context(|| {
            format!(
                "atomically sync aggregate count csv from multi-sample matrix {:?} -> {:?}",
                multi.count_csv, count_csv
            )
        })?;
    }

    if opts.counting.assignment_mode != crate::count::AssignmentMode::Unique {
        remove_optional_output(
            &unique_read_to_isoform_tsv,
            "unique read-to-isoform mapping",
        )?;
        remove_optional_output(
            &unique_assignment_provenance_tsv,
            "unique-assignment provenance",
        )?;
    }

    Ok(FullFlowResult {
        batch,
        isoform_bed,
        unused_bed,
        count_csv,
        desc_prefix,
        multi_sample,
        unique_read_to_isoform_tsv: (opts.counting.assignment_mode
            == crate::count::AssignmentMode::Unique)
            .then_some(unique_read_to_isoform_tsv),
    })
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;
    use crate::flow::artifact_manifest::read_run_manifest;

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

    fn batch_options(root: &Path) -> BatchRunOptions {
        BatchRunOptions {
            cluster_mode: ClusterMode::Clusterj,
            prepare_reads: None,
            prepare_reference: None,
            prepare_prefix: None,
            prepare: PrepareConfig::default(),
            prepare_rejected_read_tracks: 0,
            input_root: root.to_path_buf(),
            gene_list: None,
            output_root: root.to_path_buf(),
            clustering: ClusteringConfig::default(),
            counting: CountingConfig::default(),
            runtime: RuntimeConfig {
                heartbeat_seconds: 0,
                ..RuntimeConfig::default()
            },
            downsample: DownsampleConfig {
                max_reads_per_gene: 0,
                ..DownsampleConfig::default()
            },
        }
    }

    #[test]
    fn empty_reads_invalidate_stale_manifest_and_outputs() {
        let root = fresh_temp_dir("empty_reads_stale_outputs");
        let gene = GeneId::parse("GENEA").unwrap();
        let gene_dir = gene_dir_path(&root, &gene).unwrap();
        fs::create_dir_all(&gene_dir).unwrap();
        let reads = gene_artifact_path(&root, &gene, "_nano.bed").unwrap();
        let reference = gene_artifact_path(&root, &gene, "_gff.bed").unwrap();
        fs::write(&reads, "").unwrap();
        fs::write(
            &reference,
            "chr1\t0\t10\tref\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tGENEA\n",
        )
        .unwrap();
        let isoform = gene_artifact_path(&root, &gene, "_simple_coveragej.bed").unwrap();
        let unused = gene_artifact_path(&root, &gene, "_unused.bed").unwrap();
        let mapping = gene_artifact_path(&root, &gene, "_read_to_isoform.tsv").unwrap();
        let downsample = gene_dir.join("downsample.tsv");
        let manifest = gene_dir.join(MANIFEST_FILE_NAME);
        for path in [&isoform, &unused, &mapping, &downsample] {
            fs::write(path, "stale biological result\n").unwrap();
        }
        fs::write(&manifest, "{\"stale\":true}\n").unwrap();

        let result = process_gene_inputs(&gene, &batch_options(&root), &reads, &reference).unwrap();
        assert_eq!(
            result.outcome,
            GeneOutcome::Skipped(GeneSkipReason::EmptyReads)
        );
        assert_eq!(result.resume_reason, "empty_reads_stale_outputs_removed");
        let empty_manifest = read_run_manifest(&manifest).unwrap();
        assert_eq!(empty_manifest.status, "complete");
        assert!(empty_manifest
            .outputs
            .iter()
            .all(|output| output.records == 0));
        for path in [&isoform, &unused, &mapping] {
            assert_eq!(fs::read(path).unwrap(), Vec::<u8>::new());
        }
        let downsample_text = fs::read_to_string(&downsample).unwrap();
        assert!(downsample_text.starts_with("gene\toriginal_reads\tsampled_reads"));
        assert!(!downsample_text.contains("stale biological result"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn comment_only_reads_are_a_distinct_empty_outcome() {
        let root = fresh_temp_dir("comment_only_reads");
        let gene = GeneId::parse("GENEA").unwrap();
        let gene_dir = gene_dir_path(&root, &gene).unwrap();
        fs::create_dir_all(&gene_dir).unwrap();
        let reads = gene_artifact_path(&root, &gene, "_nano.bed").unwrap();
        let reference = gene_artifact_path(&root, &gene, "_gff.bed").unwrap();
        fs::write(&reads, "# comment only\n\n   \n# another comment\n").unwrap();
        fs::write(
            &reference,
            "chr1\t0\t10\tref\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tGENEA\n",
        )
        .unwrap();

        let result = process_gene_inputs(&gene, &batch_options(&root), &reads, &reference).unwrap();
        assert_eq!(
            result.outcome,
            GeneOutcome::Skipped(GeneSkipReason::EmptyReads)
        );
        assert!(gene_dir.join(MANIFEST_FILE_NAME).is_file());
        let isoform = gene_artifact_path(&root, &gene, "_simple_coveragej.bed").unwrap();
        assert_eq!(fs::metadata(isoform).unwrap().len(), 0);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn count_only_preflight_preserves_merged_outputs_when_gene_artifact_is_missing() {
        let root = fresh_temp_dir("count_only_missing_artifact");
        let input_root = fresh_temp_dir("count_only_missing_artifact_input");
        let reference = input_root.join("reference.bed");
        fs::write(
            &reference,
            "chr1\t0\t10\tref\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tGENEA\n",
        )
        .unwrap();
        fs::write(root.join("sample_gene.txt"), "GENEA\n").unwrap();
        let gene = GeneId::parse("GENEA").unwrap();
        let gene_dir = gene_dir_path(&root, &gene).unwrap();
        fs::create_dir_all(&gene_dir).unwrap();
        fs::write(
            gene_artifact_path(&root, &gene, "_simple_coveragej.bed").unwrap(),
            "",
        )
        .unwrap();
        fs::write(gene_artifact_path(&root, &gene, "_unused.bed").unwrap(), "").unwrap();
        fs::write(
            gene_artifact_path(&root, &gene, "_read_to_isoform.tsv").unwrap(),
            "",
        )
        .unwrap();
        let merged_isoform = root.join("sample_isoform.bed");
        let merged_unused = root.join("sample_unused.bed");
        let merged_mapping = root.join("sample_read_to_isoform.tsv");
        for path in [&merged_isoform, &merged_unused, &merged_mapping] {
            fs::write(path, "old complete generation\n").unwrap();
        }

        let error = run_full_flow(FullFlowOptions {
            cluster_mode: ClusterMode::Clusterj,
            reads: None,
            manifest: None,
            reference,
            output_root: root.clone(),
            prefix: "sample".to_owned(),
            prepare: PrepareConfig::default(),
            clustering: ClusteringConfig::default(),
            counting: CountingConfig::default(),
            runtime: RuntimeConfig {
                heartbeat_seconds: 0,
                ..RuntimeConfig::default()
            },
            downsample: DownsampleConfig {
                max_reads_per_gene: 0,
                ..DownsampleConfig::default()
            },
            emit_pooled_reads: false,
            count_only: true,
        })
        .expect_err("missing prepared reads must fail before merging");
        assert!(format!("{error:#}").contains("prepared_reads"));
        for path in [&merged_isoform, &merged_unused, &merged_mapping] {
            assert_eq!(
                fs::read_to_string(path).unwrap(),
                "old complete generation\n"
            );
        }
        fs::remove_dir_all(root).unwrap();
        fs::remove_dir_all(input_root).unwrap();
    }

    #[test]
    fn downsample_metadata_round_trips_and_accepts_legacy_records() {
        let dir = fresh_temp_dir("downsample_metadata_round_trip");
        let path = dir.join("downsample.tsv");
        let record = DownsampleRecord {
            gene: "GENEA".to_owned(),
            original_reads: 10,
            sampled_reads: 4,
            target_reads: 4,
            seed: 17,
            scale_factor: 2.5,
            input_fingerprint: Some("fnv1a64:0123456789abcdef".to_owned()),
        };
        write_downsample_records(&path, std::slice::from_ref(&record)).expect("write metadata");
        assert_eq!(
            read_per_gene_downsample(&path, "GENEA").expect("read metadata"),
            Some(record)
        );

        fs::write(
            &path,
            "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads\nGENEA\t10\t4\t2.5\t17\t4\n",
        )
        .expect("write legacy metadata");
        let legacy = read_per_gene_downsample(&path, "GENEA")
            .expect("read legacy metadata")
            .expect("legacy record");
        assert_eq!(legacy.input_fingerprint, None);
        assert_eq!(legacy.scale_factor, 2.5);
    }

    #[test]
    fn downsample_metadata_rejects_invalid_counts_and_scales() {
        let dir = fresh_temp_dir("downsample_metadata_invalid");
        let path = dir.join("downsample.tsv");
        let header = "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads\tinput_fingerprint\n";
        for (record, expected_error) in [
            (
                "GENEA\t10\t0\t1\t1\t4\tfnv1a64:0123456789abcdef\n",
                "original_reads >= sampled_reads > 0",
            ),
            (
                "GENEA\t3\t4\t0.75\t1\t4\tfnv1a64:0123456789abcdef\n",
                "original_reads >= sampled_reads > 0",
            ),
            (
                "GENEA\t10\t4\tNaN\t1\t4\tfnv1a64:0123456789abcdef\n",
                "finite positive value",
            ),
            (
                "GENEA\t10\t4\t2.5\t1\t3\tfnv1a64:0123456789abcdef\n",
                "target_reads >= sampled_reads > 0",
            ),
            (
                "GENEA\t10\t4\t2.5\t1\t4\tnot-a-fingerprint\n",
                "invalid input_fingerprint",
            ),
        ] {
            fs::write(&path, format!("{header}{record}")).expect("write invalid metadata");
            let error = read_downsample_records(&path).expect_err("metadata must be rejected");
            assert!(format!("{error:#}").contains(expected_error), "{error:#}");
        }
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

    #[test]
    fn failed_merge_preserves_previously_published_output() {
        let dir = fresh_temp_dir("merge_files_atomic_failure");
        let input = dir.join("input.txt");
        let missing = dir.join("missing.txt");
        let out = dir.join("out.txt");
        fs::write(&input, "new partial content\n").unwrap();
        fs::write(&out, "old complete content\n").unwrap();

        let error = merge_files(&[input, missing], &out).expect_err("merge must fail");
        assert!(format!("{error:#}").contains("missing.txt"));
        assert_eq!(fs::read_to_string(&out).unwrap(), "old complete content\n");
        assert!(fs::read_dir(&dir).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains(".tmp.")));
    }

    #[test]
    fn duplicate_merged_isoform_ids_are_rejected_before_atomic_publish() {
        let dir = fresh_temp_dir("duplicate_merged_isoform_ids");
        let first = dir.join("first.bed");
        let second = dir.join("second.bed");
        let out = dir.join("catalog.bed");
        fs::write(
            &first,
            "chr1\t0\t10\tduplicate\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tnone\tGENEA\n",
        )
        .unwrap();
        fs::write(
            &second,
            "chr2\t20\t30\tduplicate\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tnone\tGENEB\n",
        )
        .unwrap();
        fs::write(&out, "old valid catalog\n").unwrap();

        let error = merge_isoform_files(&[first, second], &out).unwrap_err();
        assert!(format!("{error:#}").contains("duplicate isoform id"));
        assert_eq!(fs::read_to_string(&out).unwrap(), "old valid catalog\n");
    }

    #[test]
    fn identical_multi_gene_reference_isoforms_are_merged_once() {
        let dir = fresh_temp_dir("identical_multi_gene_reference_isoform");
        let first = dir.join("first.bed");
        let second = dir.join("second.bed");
        let out = dir.join("catalog.bed");
        let shared = "chr1\t0\t10\tshared_ref\t0\t+\t0\t0\t0\t1\t10,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tGENEA||GENEB\n";
        fs::write(&first, shared).unwrap();
        fs::write(&second, shared).unwrap();

        merge_isoform_files(&[first, second], &out).unwrap();
        let records = read_bed12(&out)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "shared_ref");
    }

    #[test]
    fn unique_assignment_rejects_missing_gene_folder_reads() {
        let dir = fresh_temp_dir("unique_requires_gene_reads");
        let gene = "GENEA".to_owned();
        let gene_dir = dir.join(&gene);
        fs::create_dir_all(&gene_dir).unwrap();
        fs::write(
            gene_dir.join(format!("{gene}_simple_coveragej.bed")),
            "chr1\t100\t150\tiso1\t0\t+\t0\t0\t0\t1\t50,\t0,\tnone\tnone\tnone\tnone\tnone\tGENEA\n",
        )
        .unwrap();
        fs::write(
            gene_dir.join(format!("{gene}_read_to_isoform.tsv")),
            "read1\tiso1\n",
        )
        .unwrap();

        let error = select_unique_read_to_isoform_by_gene(
            &dir,
            &[gene],
            ClusterMode::Clusterj,
            crate::count::UniqueAssignmentOptions::default(),
            InvalidReadPolicy::Skip,
        )
        .expect_err("missing per-gene reads must fail strict assignment");
        assert!(format!("{error:#}").contains("requires every selected gene artifact"));
    }

    #[test]
    fn unique_assignment_isolated_to_each_gene_folder() {
        let dir = fresh_temp_dir("unique_gene_folder_isolation");
        let gene_a = "GENEA".to_owned();
        let gene_b = "GENEB".to_owned();
        let gene_a_dir = dir.join(&gene_a);
        let gene_b_dir = dir.join(&gene_b);
        fs::create_dir_all(&gene_a_dir).unwrap();
        fs::create_dir_all(&gene_b_dir).unwrap();

        fs::write(
            gene_a_dir.join(format!("{gene_a}_nano.bed")),
            "chr1\t140\t160\tread_a\t0\t+\t0\t0\t0\t1\t20,\t0,\tnone\tnone\tnone\tnone\tnone\tGENEA\n",
        )
        .unwrap();
        fs::write(
            gene_a_dir.join(format!("{gene_a}_simple_coveragej.bed")),
            "chr1\t100\t200\tiso_a_retained_like\t0\t+\t0\t0\t0\t1\t100,\t0,\tnone\tnone\tnone\tnone\tnone\tGENEA\n",
        )
        .unwrap();
        fs::write(
            gene_a_dir.join(format!("{gene_a}_read_to_isoform.tsv")),
            "read_a\tiso_a_retained_like\n",
        )
        .unwrap();

        fs::write(
            gene_b_dir.join(format!("{gene_b}_nano.bed")),
            "chr1\t100\t200\tread_b\t0\t+\t0\t0\t0\t2\t30,20,\t0,80,\tnone\tnone\tnone\tnone\tnone\tGENEB\n",
        )
        .unwrap();
        fs::write(
            gene_b_dir.join(format!("{gene_b}_simple_coveragej.bed")),
            "chr1\t100\t200\tiso_b_spliced\t0\t+\t0\t0\t0\t2\t30,20,\t0,80,\tnone\tnone\tnone\tnone\tnone\tGENEB\n",
        )
        .unwrap();
        fs::write(
            gene_b_dir.join(format!("{gene_b}_read_to_isoform.tsv")),
            "read_b\tiso_b_spliced\n",
        )
        .unwrap();

        let selected = select_unique_read_to_isoform_by_gene(
            &dir,
            &[gene_a, gene_b],
            ClusterMode::Clusterj,
            crate::count::UniqueAssignmentOptions::default(),
            InvalidReadPolicy::Skip,
        )
        .expect("unique assignment should run per gene folder");
        assert_eq!(
            selected,
            vec![
                ("read_a".to_owned(), "iso_a_retained_like".to_owned()),
                ("read_b".to_owned(), "iso_b_spliced".to_owned()),
            ]
        );
    }
}
