use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Context;

use crate::annotate::addgene::{add_gene_no_dedup, dedup_longest_by_name, AddGeneOpts};
use crate::flow::artifact_layout::validate_pipeline_gene_namespace;
use crate::flow::artifact_manifest::{atomic_copy, atomic_write_with};
use crate::flow::config::InvalidReadPolicy;
use crate::flow::path_key::{
    ensure_destination_within, gene_artifact_path, gene_dir_path,
    reject_external_inputs_in_output_root, validate_gene_ids, write_gene_id_marker,
    write_gene_path_map, GeneId, SafePathComponent,
};
use crate::io::bed::{read_bed12, write_bed12_to_writer, BedError, RejectedReadRecord};
use crate::io::manifest::SampleRow;
use crate::model::{Strand, Transcript};
use crate::sample::tagged_read_name;
// For large inputs, `preparedir` is implemented as a bucketed two-pass pipeline to bound memory.
// We keep the original in-memory implementation for smaller inputs to preserve deterministic
// per-gene sorting and the existing golden outputs.
const PREPARE_IN_MEMORY_MAX_BYTES: u64 = 64 * 1024 * 1024;

// Must be a power of two for `bucket = hash & (BUCKET_COUNT - 1)`.
const BUCKET_COUNT: usize = 256;

// Reference span index bin size.
const BIN_SIZE: u32 = 16_384;

// Max number of open per-gene read writers.
const GENE_WRITER_CACHE_CAPACITY: usize = 128;

static PREPARE_TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

struct PrepareTempDir(PathBuf);

impl PrepareTempDir {
    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for PrepareTempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[derive(Clone, Debug)]
pub struct PrepareDirResult {
    pub genes: Vec<String>,
    pub dedup_reads: usize,
    pub novel_reads: usize,
    /// Number of malformed/identity-invalid input read tracks excluded during preparation.
    pub rejected_read_tracks: usize,
    /// Auditable TSV containing the excluded read-track diagnostics.
    pub rejected_reads_path: PathBuf,
}

fn read_input_reads(
    path: &Path,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<(Vec<Transcript>, Vec<RejectedReadRecord>)> {
    let mut reader = read_bed12(path).with_context(|| format!("open reads {path:?}"))?;
    let mut reads = Vec::new();
    loop {
        let next = match invalid_read_policy {
            InvalidReadPolicy::Skip => reader.next_recovering_read(),
            InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse reads {path:?}"))?;
        let Some(read) = next else {
            break;
        };
        reads.push(read);
    }
    Ok((reads, reader.take_rejected_reads()))
}

fn publish_prepare_rejections(
    mut result: PrepareDirResult,
    output_root: &Path,
    prefix: &str,
    rejected_reads: &[RejectedReadRecord],
) -> anyhow::Result<PrepareDirResult> {
    let prefix = SafePathComponent::parse("output prefix", prefix)?;
    let path = prefixed_output_path(output_root, &prefix, "_rejected_reads.tsv")?;
    atomic_write_with(&path, |file| {
        crate::io::bed::write_rejected_reads_tsv_to_writer(file, rejected_reads)
            .map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write rejected reads {path:?}"))?;
    result.rejected_read_tracks = rejected_reads.len();
    result.rejected_reads_path = path;
    publish_prepared_gene_list(output_root, &prefix, &result.genes)?;
    Ok(result)
}

fn reference_gene_field(tx: &Transcript) -> &str {
    tx.metadata().gene_id_field().unwrap_or(tx.name.as_str())
}

fn assigned_gene_field(tx: &Transcript) -> &str {
    tx.metadata().gene_id_field().unwrap_or("none")
}

fn genes(field: &str) -> impl Iterator<Item = &str> {
    field
        .split("||")
        .map(str::trim)
        .filter(|g| !g.is_empty() && *g != "none")
}

fn prefixed_output_path(
    output_root: &Path,
    prefix: &SafePathComponent,
    suffix: &str,
) -> anyhow::Result<PathBuf> {
    let path = output_root.join(format!("{}{suffix}", prefix.as_str()));
    ensure_destination_within(output_root, &path)?;
    Ok(path)
}

fn validate_pooled_reads_output(
    output_root: &Path,
    prefix: &SafePathComponent,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<()> {
    let Some(path) = pooled_reads_out else {
        return Ok(());
    };
    ensure_destination_within(output_root, path).with_context(|| {
        format!("pooled reads output must remain beneath output_root: {path:?}")
    })?;
    let expected = prefixed_output_path(output_root, prefix, "_pooled_reads.bed")?;
    if path != expected {
        anyhow::bail!(
            "pooled reads output must use the reserved path {expected:?}; found {path:?}"
        );
    }
    Ok(())
}

fn create_gene_directory(output_root: &Path, gene: &GeneId) -> anyhow::Result<PathBuf> {
    let gene_dir = gene_dir_path(output_root, gene)?;
    fs::create_dir_all(&gene_dir).with_context(|| format!("create {gene_dir:?}"))?;
    ensure_destination_within(output_root, &gene_dir)?;
    write_gene_id_marker(&gene_dir, gene)?;
    Ok(gene_dir)
}

fn prepared_gene_list_path(
    output_root: &Path,
    prefix: &SafePathComponent,
) -> anyhow::Result<PathBuf> {
    prefixed_output_path(output_root, prefix, "_gene.txt")
}

fn invalidate_prepared_gene_list(
    output_root: &Path,
    prefix: &SafePathComponent,
) -> anyhow::Result<()> {
    let gene_list_path = prepared_gene_list_path(output_root, prefix)?;
    atomic_write_with(&gene_list_path, |_file| Ok(()))
        .with_context(|| format!("invalidate active prepared gene list {gene_list_path:?}"))
}

fn publish_prepared_gene_list(
    output_root: &Path,
    prefix: &SafePathComponent,
    genes: &[String],
) -> anyhow::Result<()> {
    validate_gene_ids(genes.iter().map(String::as_str))?;
    let gene_list_path = prepared_gene_list_path(output_root, prefix)?;
    atomic_write_with(&gene_list_path, |file| {
        let mut gene_list = BufWriter::new(file);
        for gene in genes {
            writeln!(gene_list, "{gene}")?;
        }
        gene_list
            .flush()
            .with_context(|| format!("flush {gene_list_path:?}"))
    })
    .with_context(|| format!("publish active prepared gene list {gene_list_path:?}"))
}

fn write_prepared_gene_path_map(
    output_root: &Path,
    prefix: &SafePathComponent,
    genes: &[String],
) -> anyhow::Result<Vec<GeneId>> {
    let gene_ids = validate_gene_ids(genes.iter().map(String::as_str))?;

    let mapping_path = prefixed_output_path(output_root, prefix, "_gene_paths.tsv")?;
    write_gene_path_map(&mapping_path, &gene_ids)?;
    Ok(gene_ids)
}

fn write_bed12_indices(
    path: &Path,
    transcripts: &[Transcript],
    indices: &[usize],
) -> Result<(), BedError> {
    let file = File::create(path).map_err(|source| BedError::IoWrite {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = std::io::BufWriter::new(file);
    write_bed12_indices_to_writer(&mut writer, transcripts, indices).map_err(|source| {
        BedError::IoWrite {
            path: path.to_path_buf(),
            source,
        }
    })
}

fn write_bed12_indices_to_writer<W: Write>(
    writer: &mut W,
    transcripts: &[Transcript],
    indices: &[usize],
) -> std::io::Result<()> {
    write_bed12_to_writer(writer, indices.iter().map(|&idx| &transcripts[idx]))?;
    writer.flush()
}

fn write_bed12_indices_for_gene_to_writer<W: Write>(
    writer: &mut W,
    transcripts: &[Transcript],
    indices: &[usize],
    gene: &GeneId,
) -> std::io::Result<()> {
    for &idx in indices {
        let mut transcript = transcripts[idx].clone();
        transcript.metadata_mut().set_gene_id(gene.as_str());
        write_bed12_to_writer(writer, std::iter::once(&transcript))?;
    }
    writer.flush()
}

fn sort_indices_by_coord(transcripts: &[Transcript], indices: &mut [usize]) {
    indices.sort_by(|&a, &b| {
        let left = &transcripts[a];
        let right = &transcripts[b];
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.tx_start.cmp(&right.tx_start))
            .then_with(|| left.tx_end.cmp(&right.tx_end))
            .then_with(|| left.name.cmp(&right.name))
    });
}

fn prepare_dir(
    reads_raw: &[Transcript],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    let prefix = SafePathComponent::parse("output prefix", prefix)?;
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    // Build and validate the complete generation before replacing any active output.
    let reads_dedup = dedup_longest_by_name(reads_raw);
    let reads_annotated = add_gene_no_dedup(&reads_dedup, refs, addgene_opts);

    let mut ref_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tx) in refs.iter().enumerate() {
        for gene in genes(reference_gene_field(tx)) {
            let gene = GeneId::parse(gene)?;
            ref_by_gene
                .entry(gene.as_str().to_owned())
                .or_default()
                .push(idx);
        }
    }

    let mut reads_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    let mut novel_indices: Vec<usize> = Vec::new();
    for (idx, tx) in reads_annotated.iter().enumerate() {
        let field = assigned_gene_field(tx);
        if field.trim() == "none" || field.trim().is_empty() {
            novel_indices.push(idx);
            continue;
        }
        for gene in genes(field) {
            let gene = GeneId::parse(gene)?;
            reads_by_gene
                .entry(gene.as_str().to_owned())
                .or_default()
                .push(idx);
        }
    }

    let mut genes: Vec<String> = reads_by_gene.keys().cloned().collect();
    genes.sort();
    let gene_ids = validate_gene_ids(genes.iter().map(String::as_str))?;
    validate_pipeline_gene_namespace(&gene_ids, Some(&prefix))?;

    // The gene list is the commit marker for a prepared generation. Empty it before
    // publishing any replacement and restore it only after every other artifact succeeds.
    invalidate_prepared_gene_list(output_root, &prefix)?;

    let dedup_path = prefixed_output_path(output_root, &prefix, "_dedup.bed")?;
    atomic_write_with(&dedup_path, |file| {
        write_bed12_to_writer(file, reads_dedup.iter()).map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write {dedup_path:?}"))?;

    let novel_path = prefixed_output_path(output_root, &prefix, "_novel.bed")?;
    atomic_write_with(&novel_path, |file| {
        write_bed12_indices_to_writer(file, &reads_annotated, &novel_indices)
            .map_err(anyhow::Error::from)
    })
    .with_context(|| format!("atomically write {novel_path:?}"))?;

    for gene in &gene_ids {
        create_gene_directory(output_root, gene)?;

        let mut ref_indices = ref_by_gene.get(gene.as_str()).cloned().unwrap_or_default();
        sort_indices_by_coord(refs, &mut ref_indices);
        let ref_path = gene_artifact_path(output_root, gene, "_gff.bed")?;
        atomic_write_with(&ref_path, |file| {
            write_bed12_indices_to_writer(file, refs, &ref_indices).map_err(anyhow::Error::from)
        })
        .with_context(|| format!("atomically write {ref_path:?}"))?;

        let mut read_indices = reads_by_gene
            .get(gene.as_str())
            .cloned()
            .unwrap_or_default();
        sort_indices_by_coord(&reads_annotated, &mut read_indices);
        let reads_path = gene_artifact_path(output_root, gene, "_nano.bed")?;
        atomic_write_with(&reads_path, |file| {
            write_bed12_indices_for_gene_to_writer(file, &reads_annotated, &read_indices, gene)
                .map_err(anyhow::Error::from)
        })
        .with_context(|| format!("atomically write {reads_path:?}"))?;
    }
    write_prepared_gene_path_map(output_root, &prefix, &genes)?;

    Ok(PrepareDirResult {
        genes,
        dedup_reads: reads_dedup.len(),
        novel_reads: novel_indices.len(),
        rejected_read_tracks: 0,
        rejected_reads_path: PathBuf::new(),
    })
}

#[cfg(test)]
fn prepare_dir_for_test(
    reads_raw: &[Transcript],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    let result = prepare_dir(reads_raw, refs, output_root, prefix, addgene_opts)?;
    publish_prepare_rejections(result, output_root, prefix, &[])
}

pub fn prepare_dir_from_paths(
    reads_bed: &Path,
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    prepare_dir_from_paths_with_policy(
        reads_bed,
        reference_bed,
        output_root,
        prefix,
        addgene_opts,
        InvalidReadPolicy::Skip,
    )
}

pub fn prepare_dir_from_paths_with_policy(
    reads_bed: &Path,
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<PrepareDirResult> {
    SafePathComponent::parse("output prefix", prefix)?;
    addgene_opts
        .validate()
        .context("invalid gene-assignment options")?;
    reject_external_inputs_in_output_root(
        output_root,
        [
            ("reads input", reads_bed),
            ("reference input", reference_bed),
        ],
    )?;
    let refs: Vec<Transcript> = read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    let reads_len = fs::metadata(reads_bed)
        .with_context(|| format!("stat reads {reads_bed:?}"))?
        .len();

    if reads_len <= PREPARE_IN_MEMORY_MAX_BYTES {
        let (reads_raw, rejected_reads) = read_input_reads(reads_bed, invalid_read_policy)?;
        let result = prepare_dir(&reads_raw, &refs, output_root, prefix, addgene_opts)?;
        return publish_prepare_rejections(result, output_root, prefix, &rejected_reads);
    }

    prepare_dir_from_paths_bucketed_with_policy(
        reads_bed,
        &refs,
        output_root,
        prefix,
        addgene_opts,
        invalid_read_policy,
    )
}

/// Prepare tagged manifest reads, optionally publishing the pooled BED at the reserved
/// `<output_root>/<prefix>_pooled_reads.bed` path.
pub fn prepare_dir_from_manifest_rows(
    sample_rows: &[SampleRow],
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<PrepareDirResult> {
    prepare_dir_from_manifest_rows_with_policy(
        sample_rows,
        reference_bed,
        output_root,
        prefix,
        addgene_opts,
        pooled_reads_out,
        InvalidReadPolicy::Skip,
    )
}

/// Policy-selectable form of [`prepare_dir_from_manifest_rows`].
///
/// When present, `pooled_reads_out` must be exactly
/// `<output_root>/<prefix>_pooled_reads.bed`; arbitrary destinations are rejected before mutation
/// so the pooled file cannot alias preparation or per-gene artifacts.
pub fn prepare_dir_from_manifest_rows_with_policy(
    sample_rows: &[SampleRow],
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<PrepareDirResult> {
    let prefix_component = SafePathComponent::parse("output prefix", prefix)?;
    addgene_opts
        .validate()
        .context("invalid gene-assignment options")?;
    reject_external_inputs_in_output_root(
        output_root,
        std::iter::once(("reference input", reference_bed)).chain(
            sample_rows
                .iter()
                .map(|row| ("sample reads input", row.reads.as_path())),
        ),
    )?;
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;
    validate_pooled_reads_output(output_root, &prefix_component, pooled_reads_out)?;

    let refs: Vec<Transcript> = read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    let total_reads_len: u64 = sample_rows
        .iter()
        .map(|row| {
            fs::metadata(&row.reads)
                .with_context(|| format!("stat reads {:?}", row.reads))
                .map(|meta| meta.len())
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .sum();

    if total_reads_len <= PREPARE_IN_MEMORY_MAX_BYTES {
        return prepare_dir_from_manifest_rows_in_memory(
            sample_rows,
            &refs,
            output_root,
            prefix,
            addgene_opts,
            pooled_reads_out,
            invalid_read_policy,
        );
    }

    prepare_dir_from_manifest_rows_bucketed_with_policy(
        sample_rows,
        &refs,
        output_root,
        prefix,
        addgene_opts,
        pooled_reads_out,
        invalid_read_policy,
    )
}

fn prepare_dir_from_manifest_rows_in_memory(
    sample_rows: &[SampleRow],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<PrepareDirResult> {
    let mut reads_raw: Vec<Transcript> = Vec::new();
    let mut rejected_reads = Vec::new();
    for row in sample_rows {
        let (reads, mut row_rejected_reads) = read_input_reads(&row.reads, invalid_read_policy)?;
        rejected_reads.append(&mut row_rejected_reads);
        for mut tx in reads {
            tx.name = tagged_read_name(&row.sample, &tx.name);
            reads_raw.push(tx);
        }
    }

    let result = prepare_dir(&reads_raw, refs, output_root, prefix, addgene_opts)?;
    if let Some(path) = pooled_reads_out {
        atomic_write_with(path, |file| {
            write_bed12_to_writer(file, reads_raw.iter()).map_err(anyhow::Error::from)
        })
        .with_context(|| format!("atomically write pooled reads output {path:?}"))?;
    }
    publish_prepare_rejections(result, output_root, prefix, &rejected_reads)
}

#[cfg(test)]
fn prepare_dir_from_paths_bucketed(
    reads_bed: &Path,
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    prepare_dir_from_paths_bucketed_with_policy(
        reads_bed,
        refs,
        output_root,
        prefix,
        addgene_opts,
        InvalidReadPolicy::Skip,
    )
}

fn prepare_dir_from_paths_bucketed_with_policy(
    reads_bed: &Path,
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<PrepareDirResult> {
    SafePathComponent::parse("output prefix", prefix)?;
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let (_bucket_dir, bucket_paths) = create_bucket_paths(output_root, prefix)?;
    let rejected_reads = partition_reads_file_into_buckets_with_policy(
        reads_bed,
        &bucket_paths,
        invalid_read_policy,
    )
    .with_context(|| format!("bucket reads {reads_bed:?}"))?;
    let res = prepare_dir_from_buckets(&bucket_paths, refs, output_root, prefix, addgene_opts)?;
    publish_prepare_rejections(res, output_root, prefix, &rejected_reads)
}

#[cfg(test)]
fn prepare_dir_from_manifest_rows_bucketed(
    sample_rows: &[SampleRow],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<PrepareDirResult> {
    prepare_dir_from_manifest_rows_bucketed_with_policy(
        sample_rows,
        refs,
        output_root,
        prefix,
        addgene_opts,
        pooled_reads_out,
        InvalidReadPolicy::Skip,
    )
}

fn prepare_dir_from_manifest_rows_bucketed_with_policy(
    sample_rows: &[SampleRow],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<PrepareDirResult> {
    let prefix_component = SafePathComponent::parse("output prefix", prefix)?;
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;
    validate_pooled_reads_output(output_root, &prefix_component, pooled_reads_out)?;

    let (bucket_dir, bucket_paths) = create_bucket_paths(output_root, prefix)?;
    let staged_pooled_reads = pooled_reads_out.map(|_| bucket_dir.path().join("pooled_reads.bed"));
    let rejected_reads = partition_manifest_into_buckets_with_policy(
        sample_rows,
        &bucket_paths,
        staged_pooled_reads.as_deref(),
        invalid_read_policy,
    )
    .with_context(|| "bucket reads from manifest".to_owned())?;
    let res = prepare_dir_from_buckets(&bucket_paths, refs, output_root, prefix, addgene_opts)?;
    if let (Some(staged), Some(final_path)) = (staged_pooled_reads.as_deref(), pooled_reads_out) {
        atomic_copy(staged, final_path)
            .with_context(|| format!("atomically publish pooled reads {final_path:?}"))?;
    }
    publish_prepare_rejections(res, output_root, prefix, &rejected_reads)
}

fn create_bucket_paths(
    output_root: &Path,
    prefix: &str,
) -> anyhow::Result<(PrepareTempDir, Vec<PathBuf>)> {
    let _prefix = SafePathComponent::parse("output prefix", prefix)?;
    let mut reserved = None;
    for _ in 0..1_000 {
        let sequence = PREPARE_TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        // Gene path keys are capped at 180 bytes. A longer scratch component makes
        // collision with any valid gene directory impossible.
        let candidate = output_root.join(format!(
            ".tc_tmp_{}_{}_{}",
            "x".repeat(181),
            std::process::id(),
            sequence
        ));
        ensure_destination_within(output_root, &candidate)?;
        match fs::create_dir(&candidate) {
            Ok(()) => {
                reserved = Some(candidate);
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("reserve preparedir temp dir {candidate:?}"));
            }
        }
    }
    let dir =
        PrepareTempDir(reserved.context("could not reserve a unique preparedir temp directory")?);

    let mut paths: Vec<PathBuf> = Vec::with_capacity(BUCKET_COUNT);
    for bucket in 0..BUCKET_COUNT {
        let path = dir.path().join(format!("reads_bucket_{bucket:03}.bed"));
        ensure_destination_within(output_root, &path)?;
        paths.push(path);
    }
    Ok((dir, paths))
}

fn partition_reads_file_into_buckets_with_policy(
    reads_bed: &Path,
    bucket_paths: &[PathBuf],
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<Vec<RejectedReadRecord>> {
    let mut writers: Vec<BufWriter<File>> = Vec::with_capacity(bucket_paths.len());
    for path in bucket_paths {
        let file = File::create(path).with_context(|| format!("create bucket {path:?}"))?;
        writers.push(BufWriter::new(file));
    }

    let mut reader = read_bed12(reads_bed).with_context(|| format!("open reads {reads_bed:?}"))?;
    loop {
        let next = match invalid_read_policy {
            InvalidReadPolicy::Skip => reader.next_recovering_read(),
            InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse reads {reads_bed:?}"))?;
        let Some(read) = next else {
            break;
        };
        let bucket = bucket_for_hash(fnv1a_hash(read.name.as_bytes()));
        write_bed12_to_writer(&mut writers[bucket], std::iter::once(&read))
            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
    }

    for (bucket, writer) in writers.iter_mut().enumerate() {
        writer
            .flush()
            .with_context(|| format!("flush bucket {:?}", bucket_paths[bucket]))?;
    }

    Ok(reader.take_rejected_reads())
}

fn partition_manifest_into_buckets_with_policy(
    sample_rows: &[SampleRow],
    bucket_paths: &[PathBuf],
    pooled_reads_out: Option<&Path>,
    invalid_read_policy: InvalidReadPolicy,
) -> anyhow::Result<Vec<RejectedReadRecord>> {
    let mut writers: Vec<BufWriter<File>> = Vec::with_capacity(bucket_paths.len());
    for path in bucket_paths {
        let file = File::create(path).with_context(|| format!("create bucket {path:?}"))?;
        writers.push(BufWriter::new(file));
    }

    let mut pooled_writer = if let Some(path) = pooled_reads_out {
        Some((
            path.to_path_buf(),
            BufWriter::new(
                File::create(path)
                    .with_context(|| format!("create pooled reads output {path:?}"))?,
            ),
        ))
    } else {
        None
    };

    let mut rejected_reads = Vec::new();
    for row in sample_rows {
        let mut reader =
            read_bed12(&row.reads).with_context(|| format!("open reads {:?}", row.reads))?;
        loop {
            let next = match invalid_read_policy {
                InvalidReadPolicy::Skip => reader.next_recovering_read(),
                InvalidReadPolicy::Fail => reader.next_strict_read(),
            }
            .with_context(|| format!("parse reads {:?}", row.reads))?;
            let Some(mut read) = next else {
                break;
            };
            read.name = tagged_read_name(&row.sample, &read.name);
            let bucket = bucket_for_hash(fnv1a_hash(read.name.as_bytes()));
            write_bed12_to_writer(&mut writers[bucket], std::iter::once(&read))
                .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
            if let Some((path, writer)) = pooled_writer.as_mut() {
                write_bed12_to_writer(writer, std::iter::once(&read))
                    .with_context(|| format!("write pooled reads {path:?}"))?;
            }
        }
        rejected_reads.extend(reader.take_rejected_reads());
    }

    for (bucket, writer) in writers.iter_mut().enumerate() {
        writer
            .flush()
            .with_context(|| format!("flush bucket {:?}", bucket_paths[bucket]))?;
    }
    if let Some((path, writer)) = pooled_writer.as_mut() {
        writer
            .flush()
            .with_context(|| format!("flush pooled reads output {path:?}"))?;
    }

    Ok(rejected_reads)
}

fn prepare_dir_from_buckets(
    bucket_paths: &[PathBuf],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    let prefix = SafePathComponent::parse("output prefix", prefix)?;
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let bucket_root = bucket_paths
        .first()
        .and_then(|path| path.parent())
        .context("bucketed preparation requires at least one bucket path")?;
    if bucket_paths
        .iter()
        .any(|path| path.parent() != Some(bucket_root))
    {
        anyhow::bail!("bucketed preparation paths do not share one staging directory");
    }
    let stage_root = bucket_root.join("prepared_outputs");
    let stage_global_root = stage_root.join("global");
    let stage_gene_root = stage_root.join("genes");
    ensure_destination_within(output_root, &stage_root)?;
    fs::create_dir_all(&stage_global_root)
        .with_context(|| format!("create prepared global staging {stage_global_root:?}"))?;
    fs::create_dir_all(&stage_gene_root)
        .with_context(|| format!("create prepared gene staging {stage_gene_root:?}"))?;

    let staged_dedup_path = stage_global_root.join("dedup.bed");
    let staged_novel_path = stage_global_root.join("novel.bed");
    let mut dedup_writer = BufWriter::new(
        File::create(&staged_dedup_path)
            .with_context(|| format!("write staged dedup {staged_dedup_path:?}"))?,
    );
    let mut novel_writer = BufWriter::new(
        File::create(&staged_novel_path)
            .with_context(|| format!("write staged novel {staged_novel_path:?}"))?,
    );

    let mut ref_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tx) in refs.iter().enumerate() {
        for gene in genes(reference_gene_field(tx)) {
            let gene = GeneId::parse(gene)?;
            ref_by_gene
                .entry(gene.as_str().to_owned())
                .or_default()
                .push(idx);
        }
    }

    let ref_index = RefOverlapIndex::build(refs);
    let mut overlap_scratch = RefOverlapScratch::new(refs.len());
    let mut gene_name_buf: Vec<&str> = Vec::new();

    let mut genes_seen: HashSet<String> = HashSet::new();
    let mut path_keys_seen: HashMap<String, String> = HashMap::new();
    let mut gene_writer_cache = GeneWriterCache::new(GENE_WRITER_CACHE_CAPACITY);

    let mut dedup_reads_total: usize = 0;
    let mut novel_reads_total: usize = 0;

    for bucket_path in bucket_paths {
        let bucket_len = fs::metadata(bucket_path)
            .with_context(|| format!("stat bucket {bucket_path:?}"))?
            .len();
        if bucket_len == 0 {
            continue;
        }

        let mut dedup_reads: Vec<Transcript> = Vec::new();
        let mut pos: HashMap<String, usize> = HashMap::new();

        let reader = read_bed12(bucket_path)
            .with_context(|| format!("open reads bucket {bucket_path:?}"))?;
        for record in reader {
            let tx = record.with_context(|| format!("parse reads bucket {bucket_path:?}"))?;
            match pos.get(&tx.name).copied() {
                None => {
                    dedup_reads.push(tx);
                    let idx = dedup_reads.len() - 1;
                    pos.insert(dedup_reads[idx].name.clone(), idx);
                }
                Some(existing_idx) => {
                    let existing_len = exon_len(&dedup_reads[existing_idx]);
                    let candidate_len = exon_len(&tx);
                    if candidate_len > existing_len {
                        dedup_reads[existing_idx] = tx;
                    }
                }
            }
        }

        dedup_reads_total += dedup_reads.len();
        write_bed12_to_writer(&mut dedup_writer, dedup_reads.iter())
            .with_context(|| format!("append staged dedup {staged_dedup_path:?}"))?;

        for tx in &mut dedup_reads {
            tx.metadata_mut().set_gene_id("none");
            ref_index.collect_gene_names(
                tx,
                refs,
                addgene_opts,
                &mut overlap_scratch,
                &mut gene_name_buf,
            );
            if !gene_name_buf.is_empty() {
                let joined = gene_name_buf.join("||");
                tx.metadata_mut().set_gene_id(joined);
            }

            let field = assigned_gene_field(tx);
            if field.trim() == "none" || field.trim().is_empty() {
                novel_reads_total += 1;
                write_bed12_to_writer(&mut novel_writer, std::iter::once(&*tx))
                    .with_context(|| format!("append staged novel {staged_novel_path:?}"))?;
                continue;
            }

            let assigned_genes = genes(field)
                .map(GeneId::parse)
                .collect::<anyhow::Result<Vec<_>>>()?;
            let original_gene_field = tx.metadata().gene_id_field().unwrap_or("none").to_owned();
            for gene in assigned_genes {
                let writer = ensure_gene_writer(
                    &gene,
                    &mut genes_seen,
                    &mut path_keys_seen,
                    &mut gene_writer_cache,
                    &ref_by_gene,
                    refs,
                    &stage_gene_root,
                )?;
                tx.metadata_mut().set_gene_id(gene.as_str());
                write_bed12_to_writer(writer, std::iter::once(&*tx)).with_context(|| {
                    format!(
                        "write staged reads for gene {:?} (stage_root={stage_gene_root:?})",
                        gene.as_str()
                    )
                })?;
            }
            tx.metadata_mut().set_gene_id(original_gene_field);
        }
    }

    dedup_writer
        .flush()
        .with_context(|| format!("flush staged dedup {staged_dedup_path:?}"))?;
    novel_writer
        .flush()
        .with_context(|| format!("flush staged novel {staged_novel_path:?}"))?;
    gene_writer_cache.flush_all()?;
    drop(dedup_writer);
    drop(novel_writer);
    drop(gene_writer_cache);

    let mut genes: Vec<String> = genes_seen.into_iter().collect();
    genes.sort();
    let gene_ids = validate_gene_ids(genes.iter().map(String::as_str))?;
    validate_pipeline_gene_namespace(&gene_ids, Some(&prefix))?;

    // Detect deterministic namespace conflicts before invalidating an earlier generation.
    for gene in &gene_ids {
        let gene_dir = gene_dir_path(output_root, gene)?;
        match fs::symlink_metadata(&gene_dir) {
            Ok(metadata) if metadata.file_type().is_symlink() || !metadata.is_dir() => {
                anyhow::bail!(
                    "prepared gene directory collides with a non-directory entry: {gene_dir:?}"
                );
            }
            Ok(_) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                return Err(error).with_context(|| format!("inspect {gene_dir:?}"));
            }
        }
    }

    invalidate_prepared_gene_list(output_root, &prefix)?;
    let dedup_path = prefixed_output_path(output_root, &prefix, "_dedup.bed")?;
    atomic_copy(&staged_dedup_path, &dedup_path)
        .with_context(|| format!("publish staged dedup {dedup_path:?}"))?;
    let novel_path = prefixed_output_path(output_root, &prefix, "_novel.bed")?;
    atomic_copy(&staged_novel_path, &novel_path)
        .with_context(|| format!("publish staged novel {novel_path:?}"))?;

    for gene in &gene_ids {
        create_gene_directory(output_root, gene)?;
        for suffix in ["_gff.bed", "_nano.bed"] {
            let staged = gene_artifact_path(&stage_gene_root, gene, suffix)?;
            let published = gene_artifact_path(output_root, gene, suffix)?;
            atomic_copy(&staged, &published).with_context(|| {
                format!("publish staged gene artifact {staged:?} -> {published:?}")
            })?;
        }
    }
    write_prepared_gene_path_map(output_root, &prefix, &genes)?;

    Ok(PrepareDirResult {
        genes,
        dedup_reads: dedup_reads_total,
        novel_reads: novel_reads_total,
        rejected_read_tracks: 0,
        rejected_reads_path: PathBuf::new(),
    })
}

fn ensure_gene_writer<'a>(
    gene: &GeneId,
    genes_seen: &mut HashSet<String>,
    path_keys_seen: &mut HashMap<String, String>,
    writer_cache: &'a mut GeneWriterCache,
    ref_by_gene: &HashMap<String, Vec<usize>>,
    refs: &[Transcript],
    stage_gene_root: &Path,
) -> anyhow::Result<&'a mut BufWriter<File>> {
    let path_key = gene.path_key();
    if let Some(previous_gene) = path_keys_seen.get(path_key.as_str()) {
        if previous_gene != gene.as_str() {
            anyhow::bail!(
                "gene path-key collision: biological IDs {previous_gene:?} and {:?} both map to {:?}",
                gene.as_str(),
                path_key.as_str()
            );
        }
    } else {
        path_keys_seen.insert(path_key.as_str().to_owned(), gene.as_str().to_owned());
    }
    let is_new = if genes_seen.contains(gene.as_str()) {
        false
    } else {
        genes_seen.insert(gene.as_str().to_owned());
        true
    };

    if is_new {
        create_gene_directory(stage_gene_root, gene)?;

        let mut ref_indices = ref_by_gene.get(gene.as_str()).cloned().unwrap_or_default();
        sort_indices_by_coord(refs, &mut ref_indices);
        let ref_path = gene_artifact_path(stage_gene_root, gene, "_gff.bed")?;
        write_bed12_indices(&ref_path, refs, &ref_indices)
            .with_context(|| format!("write {ref_path:?}"))?;
    }

    let reads_path = gene_artifact_path(stage_gene_root, gene, "_nano.bed")?;
    writer_cache.get_or_open(gene.as_str(), &reads_path, is_new)
}

fn exon_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|exon| u64::from(exon.len())).sum()
}

fn span_len(tx: &Transcript) -> u32 {
    tx.tx_end.get().saturating_sub(tx.tx_start.get())
}

fn span_overlap_len(a: &Transcript, b: &Transcript) -> u32 {
    let start = a.tx_start.get().max(b.tx_start.get());
    let end = a.tx_end.get().min(b.tx_end.get());
    end.saturating_sub(start)
}

const FNV_OFFSET_BASIS: u64 = 14695981039346656037;
const FNV_PRIME: u64 = 1099511628211;

fn fnv1a_hash(bytes: &[u8]) -> u64 {
    let mut hash: u64 = FNV_OFFSET_BASIS;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

fn bucket_for_hash(hash: u64) -> usize {
    debug_assert!(BUCKET_COUNT.is_power_of_two());
    (hash as usize) & (BUCKET_COUNT - 1)
}

struct BinnedPartitionIndex {
    max_bin: u32,
    bins: Vec<Vec<usize>>,
}

impl BinnedPartitionIndex {
    fn build(refs: &[Transcript], indices: &[usize]) -> Self {
        let max_end = indices
            .iter()
            .map(|&idx| refs[idx].tx_end.get())
            .max()
            .unwrap_or(0);
        let max_bin = max_end / BIN_SIZE;

        let mut bins: Vec<Vec<usize>> = vec![Vec::new(); (max_bin + 1) as usize];
        for &idx in indices {
            let start = refs[idx].tx_start.get();
            let end = refs[idx].tx_end.get();
            if start >= end {
                continue;
            }

            let start_bin = start / BIN_SIZE;
            let end_bin = end.saturating_sub(1) / BIN_SIZE;
            for bin in start_bin..=end_bin {
                bins[bin as usize].push(idx);
            }
        }

        Self { max_bin, bins }
    }

    fn collect_candidates(
        &self,
        query_start: u32,
        query_end: u32,
        seen: &mut [u32],
        stamp: u32,
        out: &mut Vec<usize>,
    ) {
        if query_start >= query_end {
            return;
        }

        let start_bin = query_start / BIN_SIZE;
        if start_bin > self.max_bin {
            return;
        }

        let end_bin = (query_end.saturating_sub(1) / BIN_SIZE).min(self.max_bin);
        for bin in start_bin..=end_bin {
            for &idx in &self.bins[bin as usize] {
                if seen[idx] != stamp {
                    seen[idx] = stamp;
                    out.push(idx);
                }
            }
        }
    }
}

#[derive(Default)]
struct ChromBinsBuilder {
    plus: Vec<usize>,
    minus: Vec<usize>,
    unknown: Vec<usize>,
}

struct ChromBins {
    plus: Option<BinnedPartitionIndex>,
    minus: Option<BinnedPartitionIndex>,
    unknown: Option<BinnedPartitionIndex>,
}

impl ChromBins {
    fn for_strand(&self, strand: Strand) -> Option<&BinnedPartitionIndex> {
        match strand {
            Strand::Plus => self.plus.as_ref(),
            Strand::Minus => self.minus.as_ref(),
            Strand::Unknown => self.unknown.as_ref(),
        }
    }
}

struct RefOverlapIndex {
    by_chrom: HashMap<String, ChromBins>,
}

impl RefOverlapIndex {
    fn build(refs: &[Transcript]) -> Self {
        let mut builders: HashMap<String, ChromBinsBuilder> = HashMap::new();
        for (idx, tx) in refs.iter().enumerate() {
            let builder = builders.entry(tx.chrom.clone()).or_default();
            match tx.strand {
                Strand::Plus => builder.plus.push(idx),
                Strand::Minus => builder.minus.push(idx),
                Strand::Unknown => builder.unknown.push(idx),
            }
        }

        let mut by_chrom: HashMap<String, ChromBins> = HashMap::with_capacity(builders.len());
        for (chrom, builder) in builders {
            let plus = if builder.plus.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.plus))
            };
            let minus = if builder.minus.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.minus))
            };
            let unknown = if builder.unknown.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.unknown))
            };
            by_chrom.insert(
                chrom,
                ChromBins {
                    plus,
                    minus,
                    unknown,
                },
            );
        }

        Self { by_chrom }
    }

    fn collect_gene_names<'a>(
        &self,
        read: &Transcript,
        refs: &'a [Transcript],
        opts: AddGeneOpts,
        scratch: &mut RefOverlapScratch,
        out: &mut Vec<&'a str>,
    ) {
        out.clear();

        let Some(chrom_bins) = self.by_chrom.get(read.chrom.as_str()) else {
            return;
        };
        let Some(index) = chrom_bins.for_strand(read.strand) else {
            return;
        };

        let start = read.tx_start.get();
        let end = read.tx_end.get();
        if start >= end {
            return;
        }

        let stamp = scratch.next_stamp();
        scratch.candidates.clear();
        index.collect_candidates(
            start,
            end,
            &mut scratch.seen,
            stamp,
            &mut scratch.candidates,
        );

        let read_len = span_len(read);
        if read_len == 0 {
            return;
        }

        for &ref_idx in &scratch.candidates {
            let reference = &refs[ref_idx];
            let overlap = span_overlap_len(read, reference);
            if overlap == 0 {
                continue;
            }

            let ref_len = span_len(reference);
            if ref_len == 0 {
                continue;
            }

            let overlap_f = overlap as f64;
            if overlap_f / (read_len as f64) < opts.fraction_read {
                continue;
            }
            if overlap_f / (ref_len as f64) < opts.fraction_ref {
                continue;
            }

            out.push(reference_gene_field(reference));
        }

        out.sort_unstable();
        out.dedup();
    }
}

struct RefOverlapScratch {
    seen: Vec<u32>,
    stamp: u32,
    candidates: Vec<usize>,
}

impl RefOverlapScratch {
    fn new(ref_len: usize) -> Self {
        Self {
            seen: vec![0; ref_len],
            stamp: 0,
            candidates: Vec::new(),
        }
    }

    fn next_stamp(&mut self) -> u32 {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.seen.fill(0);
            self.stamp = 1;
        }
        self.stamp
    }
}

struct GeneWriterEntry {
    writer: BufWriter<File>,
    last_used: u64,
}

struct GeneWriterCache {
    capacity: usize,
    tick: u64,
    entries: HashMap<String, GeneWriterEntry>,
}

impl GeneWriterCache {
    fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            tick: 0,
            entries: HashMap::new(),
        }
    }

    fn get_or_open<'a>(
        &'a mut self,
        gene: &str,
        path: &Path,
        truncate: bool,
    ) -> anyhow::Result<&'a mut BufWriter<File>> {
        self.tick = self.tick.wrapping_add(1);
        if self.entries.contains_key(gene) {
            let tick = self.tick;
            let entry = self
                .entries
                .get_mut(gene)
                .expect("contains_key checked in branch");
            entry.last_used = tick;
            return Ok(&mut entry.writer);
        }

        if self.entries.len() >= self.capacity {
            let lru_key = self
                .entries
                .iter()
                .min_by_key(|(_, entry)| entry.last_used)
                .map(|(key, _)| key.clone());
            if let Some(key) = lru_key {
                if let Some(mut entry) = self.entries.remove(&key) {
                    entry
                        .writer
                        .flush()
                        .with_context(|| format!("flush evicted writer {key:?} -> {path:?}"))?;
                }
            }
        }

        let file = if truncate {
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(path)
        } else {
            OpenOptions::new().create(true).append(true).open(path)
        }
        .with_context(|| format!("open gene reads {gene:?} at {path:?}"))?;

        self.entries.insert(
            gene.to_owned(),
            GeneWriterEntry {
                writer: BufWriter::new(file),
                last_used: self.tick,
            },
        );

        Ok(&mut self
            .entries
            .get_mut(gene)
            .expect("inserted entry missing")
            .writer)
    }

    fn flush_all(&mut self) -> anyhow::Result<()> {
        for (gene, entry) in self.entries.iter_mut() {
            entry
                .writer
                .flush()
                .with_context(|| format!("flush gene writer {gene:?}"))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::model::{Coord, Interval};
    use proptest::prelude::*;

    use super::*;

    #[derive(Default)]
    struct FlushFails;

    impl Write for FlushFails {
        fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
            Ok(buffer.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::other("injected flush failure"))
        }
    }

    #[test]
    fn indexed_bed_writer_propagates_final_flush_errors() {
        let mut writer = FlushFails;
        assert!(write_bed12_indices_to_writer(&mut writer, &[], &[]).is_err());
    }

    fn fresh_temp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "trackcluster_rs_{}_{}_{}",
            prefix,
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn bucket_temp_reservation_never_removes_a_preexisting_directory() {
        let root = fresh_temp_dir("bucket_temp_reservation");
        let (first, _) = create_bucket_paths(&root, "sample").unwrap();
        let scratch_name = first.path().file_name().unwrap().to_str().unwrap();
        assert!(scratch_name.len() > 180);
        let same_spelling_gene = GeneId::parse(scratch_name).unwrap();
        assert_ne!(same_spelling_gene.path_key().as_str(), scratch_name);
        let sentinel = first.path().join("keep.txt");
        fs::write(&sentinel, "keep\n").unwrap();

        let (second, _) = create_bucket_paths(&root, "sample").unwrap();
        assert_ne!(first.path(), second.path());
        assert_eq!(fs::read_to_string(&sentinel).unwrap(), "keep\n");

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preparedir_bucketed_creates_gene_folders_and_files() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads = repo_root.join("tests/fixtures/reads.bed");
        let reference = repo_root.join("tests/fixtures/ref.bed");

        let out_root = fresh_temp_dir("preparedir_bucketed");
        let prefix = "sample";

        let refs: Vec<Transcript> = read_bed12(&reference)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();

        let res = prepare_dir_from_paths_bucketed(
            &reads,
            &refs,
            &out_root,
            prefix,
            AddGeneOpts::default(),
        )
        .expect("preparedir bucketed run");

        assert_eq!(res.genes, vec!["GENEA".to_owned()]);

        let gene_list = out_root.join(format!("{prefix}_gene.txt"));
        assert!(gene_list.exists());
        assert_eq!(fs::read_to_string(gene_list).unwrap(), "GENEA\n");

        let dedup = out_root.join(format!("{prefix}_dedup.bed"));
        assert!(dedup.exists());
        let dedup_reads: Vec<Transcript> = read_bed12(&dedup)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(dedup_reads.len(), 1);
        assert_eq!(dedup_reads[0].name, "read_trunc");
        assert_eq!(
            dedup_reads[0].extra_fields.get(5).map(|s| s.as_str()),
            Some("none")
        );

        let novel = out_root.join(format!("{prefix}_novel.bed"));
        assert!(novel.exists());
        let novel_reads: Vec<Transcript> = read_bed12(&novel)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(novel_reads.len(), 0);

        let gene_dir = out_root.join("GENEA");
        let gene_reads_path = gene_dir.join("GENEA_nano.bed");
        let gene_ref_path = gene_dir.join("GENEA_gff.bed");
        assert!(gene_reads_path.exists());
        assert!(gene_ref_path.exists());

        let gene_reads: Vec<Transcript> = read_bed12(&gene_reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(gene_reads.len(), 1);
        assert_eq!(
            gene_reads[0].extra_fields.get(5).map(|s| s.as_str()),
            Some("GENEA")
        );

        let gene_refs: Vec<Transcript> = read_bed12(&gene_ref_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(gene_refs.len(), 2);
        for tx in &gene_refs {
            assert_eq!(tx.extra_fields.get(5).map(|s| s.as_str()), Some("GENEA"));
        }
    }

    #[test]
    fn preparedir_manifest_bucketed_writes_pooled_reads_with_sample_tags() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let manifest = repo_root.join("tests/fixtures/samples.tsv");
        let reference = repo_root.join("tests/fixtures/ref.bed");

        let out_root = fresh_temp_dir("preparedir_manifest_bucketed");
        let prefix = "pooled";

        let sample_rows = crate::io::manifest::read_manifest_tsv(&manifest).unwrap();
        let refs: Vec<Transcript> = read_bed12(&reference)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let pooled_out = out_root.join(format!("{prefix}_pooled_reads.bed"));

        let res = prepare_dir_from_manifest_rows_bucketed(
            &sample_rows,
            &refs,
            &out_root,
            prefix,
            AddGeneOpts::default(),
            Some(&pooled_out),
        )
        .expect("preparedir manifest bucketed run");

        assert_eq!(res.genes, vec!["GENEA".to_owned()]);
        assert!(pooled_out.exists());
        let pooled = fs::read_to_string(pooled_out).expect("read pooled reads");
        assert!(pooled.contains("\tS1::read_s1\t"));
        assert!(pooled.contains("\tS2::read_s2\t"));
    }

    #[test]
    fn manifest_paths_reject_empty_read_id_before_sample_tagging() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reference = repo_root.join("tests/fixtures/ref.bed");
        let out_root = fresh_temp_dir("preparedir_manifest_empty_read_id");
        let input_root = fresh_temp_dir("preparedir_manifest_empty_read_id_input");
        let reads = input_root.join("reads.bed");
        fs::write(
            &reads,
            concat!(
                "chr1\t120\t150\t\t0\t+\t0\t0\t0\t2\t10,10,\t0,20,\n",
                "chr1\t120\t150\tgood\t0\t+\t0\t0\t0\t2\t10,10,\t0,20,\n",
            ),
        )
        .unwrap();
        let rows = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: reads.clone(),
        }];
        let pooled = out_root.join("sample_pooled_reads.bed");
        let result = prepare_dir_from_manifest_rows_with_policy(
            &rows,
            &reference,
            &out_root,
            "sample",
            AddGeneOpts::default(),
            Some(&pooled),
            InvalidReadPolicy::Skip,
        )
        .unwrap();
        assert_eq!(result.rejected_read_tracks, 1);
        let pooled_text = fs::read_to_string(&pooled).unwrap();
        assert!(pooled_text.contains("\tS1::good\t"), "{pooled_text}");
        assert!(!pooled_text.contains("\tS1::\t"), "{pooled_text}");

        let bucket_root = fresh_temp_dir("preparedir_manifest_empty_read_id_bucket");
        let (_bucket_dir, bucket_paths) = create_bucket_paths(&bucket_root, "sample").unwrap();
        let bucket_pooled = bucket_root.join("sample_pooled_reads.bed");
        let rejected = partition_manifest_into_buckets_with_policy(
            &rows,
            &bucket_paths,
            Some(&bucket_pooled),
            InvalidReadPolicy::Skip,
        )
        .unwrap();
        assert_eq!(rejected.len(), 1);
        let bucket_pooled_text = fs::read_to_string(bucket_pooled).unwrap();
        assert!(bucket_pooled_text.contains("\tS1::good\t"));
        assert!(!bucket_pooled_text.contains("\tS1::\t"));
        fs::remove_dir_all(input_root).unwrap();
    }

    #[test]
    fn manifest_preparation_rejects_noncanonical_pooled_output_before_mutation() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let manifest = repo_root.join("tests/fixtures/samples.tsv");
        let reference = repo_root.join("tests/fixtures/ref.bed");
        let rows = crate::io::manifest::read_manifest_tsv(&manifest).unwrap();
        let output_root = fresh_temp_dir("preparedir_manifest_bad_pooled_path");
        let alias = output_root.join("sample_gene.txt");
        fs::write(&alias, "previous committed genes\n").unwrap();

        let error = prepare_dir_from_manifest_rows_with_policy(
            &rows,
            &reference,
            &output_root,
            "sample",
            AddGeneOpts::default(),
            Some(&alias),
            InvalidReadPolicy::Skip,
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("must use the reserved path"),
            "{error:#}"
        );
        assert_eq!(
            fs::read_to_string(&alias).unwrap(),
            "previous committed genes\n"
        );
        assert!(!output_root.join("sample_dedup.bed").exists());
        fs::remove_dir_all(output_root).unwrap();
    }

    fn assert_plain_bed_novel_result(out_root: &Path, prefix: &str, result: &PrepareDirResult) {
        assert_eq!(result.genes, vec!["GENEA".to_owned()]);
        assert_eq!(result.dedup_reads, 2);
        assert_eq!(result.novel_reads, 1);

        let novel_path = out_root.join(format!("{prefix}_novel.bed"));
        let novel_reads: Vec<Transcript> = read_bed12(&novel_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(novel_reads.len(), 1);
        assert_eq!(novel_reads[0].name, "tx2");
        assert_eq!(novel_reads[0].metadata().gene_id_field(), Some("none"));
        assert!(!out_root.join("tx2").exists());

        let gene_reads_path = out_root.join("GENEA").join("GENEA_nano.bed");
        let gene_reads: Vec<Transcript> = read_bed12(&gene_reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(gene_reads.len(), 1);
        assert_eq!(gene_reads[0].name, "tx1");
        assert_eq!(gene_reads[0].metadata().gene_id(), Some("GENEA"));
    }

    #[test]
    fn preparedir_plain_bed_unmatched_read_is_novel_in_memory_and_bucketed() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/minimal.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let prefix = "plain";

        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();

        let in_memory_root = fresh_temp_dir("preparedir_plain_in_memory");
        let in_memory = prepare_dir_for_test(
            &reads,
            &refs,
            &in_memory_root,
            prefix,
            AddGeneOpts::default(),
        )
        .expect("in-memory plain BED preparation");
        assert_plain_bed_novel_result(&in_memory_root, prefix, &in_memory);

        let bucketed_root = fresh_temp_dir("preparedir_plain_bucketed");
        let bucketed = prepare_dir_from_paths_bucketed(
            &reads_path,
            &refs,
            &bucketed_root,
            prefix,
            AddGeneOpts::default(),
        )
        .expect("bucketed plain BED preparation");
        assert_plain_bed_novel_result(&bucketed_root, prefix, &bucketed);

        fs::remove_dir_all(in_memory_root).unwrap();
        fs::remove_dir_all(bucketed_root).unwrap();
    }

    #[test]
    fn preparedir_rejects_unsafe_gene_ids_and_prefixes() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/reads.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let mut refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        for reference in &mut refs {
            reference.metadata_mut().set_gene_id("../escape");
        }

        let out_root = fresh_temp_dir("preparedir_unsafe_gene");
        let escaped = out_root.parent().unwrap().join("escape");
        let err = prepare_dir(&reads, &refs, &out_root, "sample", AddGeneOpts::default())
            .unwrap_err()
            .to_string();
        assert!(err.contains("gene id"));
        assert!(!escaped.exists());

        let bucketed_root = fresh_temp_dir("preparedir_unsafe_gene_bucketed");
        let err = prepare_dir_from_paths_bucketed(
            &reads_path,
            &refs,
            &bucketed_root,
            "sample",
            AddGeneOpts::default(),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("gene id"), "{err}");
        assert!(!escaped.exists());

        let err = prepare_dir(
            &reads,
            &refs,
            &out_root,
            "../unsafe-prefix",
            AddGeneOpts::default(),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("output prefix"));

        fs::remove_dir_all(out_root).unwrap();
        fs::remove_dir_all(bucketed_root).unwrap();
    }

    #[test]
    fn failed_reprepare_preserves_the_previous_committed_generation() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/reads.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let mut invalid_refs = refs.clone();
        for reference in &mut invalid_refs {
            reference.metadata_mut().set_gene_id("../invalid");
        }

        let relative_paths = [
            "sample_dedup.bed",
            "sample_novel.bed",
            "sample_gene.txt",
            "sample_gene_paths.tsv",
            "sample_rejected_reads.tsv",
            "GENEA/.trackcluster_gene_id",
            "GENEA/GENEA_gff.bed",
            "GENEA/GENEA_nano.bed",
        ];

        for (label, bucketed) in [("in_memory", false), ("bucketed", true)] {
            let root = fresh_temp_dir(&format!("preparedir_failed_reprepare_{label}"));
            if bucketed {
                prepare_dir_from_paths_bucketed(
                    &reads_path,
                    &refs,
                    &root,
                    "sample",
                    AddGeneOpts::default(),
                )
                .unwrap();
            } else {
                prepare_dir_for_test(&reads, &refs, &root, "sample", AddGeneOpts::default())
                    .unwrap();
            }
            let before = relative_paths
                .iter()
                .map(|relative| fs::read(root.join(relative)).unwrap())
                .collect::<Vec<_>>();

            let error = if bucketed {
                prepare_dir_from_paths_bucketed(
                    &reads_path,
                    &invalid_refs,
                    &root,
                    "sample",
                    AddGeneOpts::default(),
                )
                .unwrap_err()
            } else {
                prepare_dir(
                    &reads,
                    &invalid_refs,
                    &root,
                    "sample",
                    AddGeneOpts::default(),
                )
                .unwrap_err()
            };
            assert!(error.to_string().contains("gene id"), "{error:#}");
            for (relative, expected) in relative_paths.iter().zip(before) {
                assert_eq!(
                    fs::read(root.join(relative)).unwrap(),
                    expected,
                    "failed {label} rerun changed {relative}"
                );
            }
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn late_publish_failure_leaves_the_generation_commit_marker_empty() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/reads.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();

        for (label, bucketed) in [("in_memory", false), ("bucketed", true)] {
            let root = fresh_temp_dir(&format!("preparedir_late_failure_{label}"));
            if bucketed {
                prepare_dir_from_paths_bucketed(
                    &reads_path,
                    &refs,
                    &root,
                    "sample",
                    AddGeneOpts::default(),
                )
                .unwrap();
            } else {
                prepare_dir_for_test(&reads, &refs, &root, "sample", AddGeneOpts::default())
                    .unwrap();
            }
            assert_eq!(
                fs::read_to_string(root.join("sample_gene.txt")).unwrap(),
                "GENEA\n"
            );

            let mapping = root.join("sample_gene_paths.tsv");
            fs::remove_file(&mapping).unwrap();
            fs::create_dir(&mapping).unwrap();
            let result = if bucketed {
                prepare_dir_from_paths_bucketed(
                    &reads_path,
                    &refs,
                    &root,
                    "sample",
                    AddGeneOpts::default(),
                )
            } else {
                prepare_dir_for_test(&reads, &refs, &root, "sample", AddGeneOpts::default())
            };
            assert!(
                result.is_err(),
                "late {label} publication unexpectedly succeeded"
            );
            assert_eq!(
                fs::read(root.join("sample_gene.txt")).unwrap(),
                b"",
                "late {label} failure left a mixed generation active"
            );
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn preparedir_rejects_gene_keys_that_collide_with_global_artifacts() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/reads.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let original_refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        for reserved_name in [
            "sample_dedup.bed",
            "sample_isoform.bed",
            "clusterj_batch_gene_paths.tsv",
        ] {
            let mut refs = original_refs.clone();
            for reference in &mut refs {
                reference.metadata_mut().set_gene_id(reserved_name);
            }
            let root = fresh_temp_dir("preparedir_reserved_global_name");
            let error =
                prepare_dir(&reads, &refs, &root, "sample", AddGeneOpts::default()).unwrap_err();
            assert!(
                error.to_string().contains("reserved top-level"),
                "reserved={reserved_name:?}: {error:#}"
            );
            assert!(!root.join(reserved_name).exists());
            assert!(!root.join("sample_gene.txt").exists());
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn preparedir_round_trips_unicode_and_long_gene_ids_in_both_paths() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads_path = repo_root.join("tests/fixtures/reads.bed");
        let reference_path = repo_root.join("tests/fixtures/ref.bed");
        let reads: Vec<Transcript> = read_bed12(&reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let original_refs: Vec<Transcript> = read_bed12(&reference_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();

        for (label, biological_gene) in [
            ("unicode", "基因-α.1".to_owned()),
            ("long", "very-long-biological-gene-".repeat(30)),
        ] {
            let mut refs = original_refs.clone();
            for reference in &mut refs {
                reference
                    .metadata_mut()
                    .set_gene_id(biological_gene.clone());
            }
            let gene = GeneId::parse(&biological_gene).unwrap();
            let key = gene.path_key();
            assert!(key.as_str().len() < 200);
            if label == "unicode" {
                assert_ne!(key.as_str(), gene.as_str());
            }

            let in_memory_root = fresh_temp_dir(&format!("preparedir_{label}_in_memory"));
            let in_memory = prepare_dir_for_test(
                &reads,
                &refs,
                &in_memory_root,
                "sample",
                AddGeneOpts::default(),
            )
            .expect("in-memory special gene preparation");

            let bucketed_root = fresh_temp_dir(&format!("preparedir_{label}_bucketed"));
            let bucketed = prepare_dir_from_paths_bucketed(
                &reads_path,
                &refs,
                &bucketed_root,
                "sample",
                AddGeneOpts::default(),
            )
            .expect("bucketed special gene preparation");

            for (root, result) in [(&in_memory_root, &in_memory), (&bucketed_root, &bucketed)] {
                assert_eq!(result.genes, vec![biological_gene.clone()]);
                assert_eq!(
                    fs::read_to_string(root.join("sample_gene.txt")).unwrap(),
                    format!("{biological_gene}\n")
                );
                let mapping =
                    crate::flow::path_key::read_gene_path_map(&root.join("sample_gene_paths.tsv"))
                        .unwrap();
                assert_eq!(mapping, vec![gene.clone()]);

                let gene_dir = root.join(key.as_str());
                assert!(gene_dir.is_dir());
                assert_eq!(
                    fs::read_to_string(gene_dir.join(crate::flow::path_key::GENE_ID_MARKER_FILE))
                        .unwrap(),
                    format!("{biological_gene}\n")
                );
                assert!(gene_dir.join(format!("{}_gff.bed", key.as_str())).exists());
                assert!(gene_dir.join(format!("{}_nano.bed", key.as_str())).exists());
            }

            fs::remove_dir_all(in_memory_root).unwrap();
            fs::remove_dir_all(bucketed_root).unwrap();
        }
    }

    fn property_transcript(
        name: String,
        start: u32,
        end: u32,
        strand: Strand,
        gene: &str,
        is_reference: bool,
    ) -> Transcript {
        let mut extra_fields = vec!["none".to_owned(); 8];
        extra_fields[4] = if is_reference {
            "isoform_anno".to_owned()
        } else {
            "nanopore_read".to_owned()
        };
        extra_fields[5] = gene.to_owned();
        Transcript::new(
            "chrP".to_owned(),
            strand,
            Coord::new(start),
            Coord::new(end),
            name,
            vec![Interval::new(Coord::new(start), Coord::new(end)).unwrap()],
            crate::model::Bed12Attrs {
                score: 0,
                thick_start: Coord::new(0),
                thick_end: Coord::new(0),
                item_rgb: "0".to_owned(),
                extra_fields,
            },
        )
        .unwrap()
    }

    fn assert_per_gene_read_metadata(root: &Path) {
        for gene in ["GENEA", "GENEB"] {
            let reads_path = root.join(gene).join(format!("{gene}_nano.bed"));
            let records = read_bed12(&reads_path)
                .unwrap()
                .collect::<Result<Vec<_>, BedError>>()
                .unwrap();
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].name, "multi_gene_read");
            assert_eq!(records[0].metadata().gene_id(), Some(gene));
        }

        let dedup = read_bed12(root.join("sample_dedup.bed"))
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(dedup.len(), 1);
        assert_eq!(dedup[0].metadata().gene_id_field(), Some("none"));
    }

    #[test]
    fn preparedir_localizes_multi_gene_reads_in_memory_and_bucketed_outputs() {
        let references = vec![
            property_transcript("ref_a".to_owned(), 100, 200, Strand::Plus, "GENEA", true),
            property_transcript("ref_b".to_owned(), 100, 200, Strand::Plus, "GENEB", true),
        ];
        let read = property_transcript(
            "multi_gene_read".to_owned(),
            120,
            180,
            Strand::Plus,
            "none",
            false,
        );

        let in_memory_root = fresh_temp_dir("preparedir_multi_gene_in_memory");
        let in_memory = prepare_dir_for_test(
            std::slice::from_ref(&read),
            &references,
            &in_memory_root,
            "sample",
            AddGeneOpts::default(),
        )
        .unwrap();
        assert_eq!(in_memory.genes, ["GENEA", "GENEB"]);
        assert_per_gene_read_metadata(&in_memory_root);

        let bucketed_root = fresh_temp_dir("preparedir_multi_gene_bucketed");
        let reads_path = bucketed_root.join("reads.bed");
        crate::io::bed::write_bed12(&reads_path, std::iter::once(&read)).unwrap();
        let bucketed = prepare_dir_from_paths_bucketed(
            &reads_path,
            &references,
            &bucketed_root,
            "sample",
            AddGeneOpts::default(),
        )
        .unwrap();
        assert_eq!(bucketed.genes, ["GENEA", "GENEB"]);
        assert_per_gene_read_metadata(&bucketed_root);

        fs::remove_dir_all(in_memory_root).unwrap();
        fs::remove_dir_all(bucketed_root).unwrap();
    }

    fn records_by_name(path: &Path) -> HashMap<String, Transcript> {
        read_bed12(path)
            .unwrap()
            .map(|record| {
                let record = record.unwrap();
                (record.name.clone(), record)
            })
            .collect()
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        #[test]
        fn in_memory_and_bucketed_preparation_are_equivalent(
            raw_reads in proptest::collection::vec((0u32..280, 1u32..80, 0u8..3), 0..8)
        ) {
            let reference = property_transcript(
                "ref".to_owned(),
                100,
                200,
                Strand::Plus,
                "GENE-P",
                true,
            );
            let reads: Vec<_> = raw_reads
                .into_iter()
                .enumerate()
                .map(|(index, (start, length, strand_code))| {
                    let strand = match strand_code {
                        0 => Strand::Plus,
                        1 => Strand::Minus,
                        _ => Strand::Unknown,
                    };
                    property_transcript(
                        format!("read-{index}"),
                        start,
                        start + length,
                        strand,
                        "none",
                        false,
                    )
                })
                .collect();

            let in_memory_root = fresh_temp_dir("property_in_memory");
            let bucketed_root = fresh_temp_dir("property_bucketed");
            let reads_path = bucketed_root.join("input.bed");
            crate::io::bed::write_bed12(&reads_path, reads.iter()).unwrap();

            let in_memory = prepare_dir_for_test(
                &reads,
                std::slice::from_ref(&reference),
                &in_memory_root,
                "sample",
                AddGeneOpts::default(),
            ).unwrap();
            let bucketed = prepare_dir_from_paths_bucketed(
                &reads_path,
                &[reference],
                &bucketed_root,
                "sample",
                AddGeneOpts::default(),
            ).unwrap();

            prop_assert_eq!(&in_memory.genes, &bucketed.genes);
            prop_assert_eq!(in_memory.dedup_reads, bucketed.dedup_reads);
            prop_assert_eq!(in_memory.novel_reads, bucketed.novel_reads);
            for relative in ["sample_dedup.bed", "sample_novel.bed"] {
                prop_assert_eq!(
                    records_by_name(&in_memory_root.join(relative)),
                    records_by_name(&bucketed_root.join(relative)),
                );
            }
            for relative in ["sample_gene.txt", "sample_gene_paths.tsv"] {
                prop_assert_eq!(
                    fs::read(in_memory_root.join(relative)).unwrap(),
                    fs::read(bucketed_root.join(relative)).unwrap(),
                );
            }
            for gene in &in_memory.genes {
                let gene = GeneId::parse(gene).unwrap();
                let key = gene.path_key();
                for suffix in ["_gff.bed", "_nano.bed"] {
                    let file = format!("{}{suffix}", key.as_str());
                    prop_assert_eq!(
                        records_by_name(&in_memory_root.join(key.as_str()).join(&file)),
                        records_by_name(&bucketed_root.join(key.as_str()).join(&file)),
                    );
                }
            }

            fs::remove_dir_all(in_memory_root).unwrap();
            fs::remove_dir_all(bucketed_root).unwrap();
        }
    }
}
