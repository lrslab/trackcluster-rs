//! Deterministic molecule-level pseudo-sample bundle generation.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::io::bam_subset::{split_bam_by_query_name, BamSubsetTarget};
use crate::io::mod_calls::{
    new_observation_writer, scan_canonical_observations_tsv, write_assay_metadata_to_writer,
    write_observation_record,
};
use crate::sample::{split_tagged_read_name, tagged_read_name, SAMPLE_DELIM};

/// Relationship among generated pseudo-samples.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SubsampleMode {
    /// Assign each selected source molecule to at most one pseudo-sample.
    #[default]
    Disjoint,
    /// Draw each pseudo-sample independently without replacement within that sample.
    Independent,
}

impl SubsampleMode {
    /// Stable CLI and provenance token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Disjoint => "disjoint",
            Self::Independent => "independent",
        }
    }
}

impl std::fmt::Display for SubsampleMode {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for SubsampleMode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "disjoint" => Ok(Self::Disjoint),
            "independent" => Ok(Self::Independent),
            _ => Err(format!(
                "invalid subsample mode {value:?}; expected disjoint or independent"
            )),
        }
    }
}

/// Inputs and deterministic sampling settings for one pseudo-sample bundle.
#[derive(Clone, Debug)]
pub struct SubsampleOptions {
    /// TrackCluster sample manifest containing the high-coverage source sample.
    pub manifest: PathBuf,
    /// Final unique read-to-isoform assignment TSV.
    pub read_to_isoform: PathBuf,
    /// Modification manifest containing one or more assays for the source sample.
    pub mod_manifest: PathBuf,
    /// Source sample identifier.
    pub source_sample: String,
    /// Prefix used to generate `<prefix>_001`, `<prefix>_002`, and so on.
    pub sample_prefix: String,
    /// Number of pseudo-samples.
    pub replicates: usize,
    /// Source read molecules selected per pseudo-sample.
    pub reads_per_sample: usize,
    /// Sampling relationship among pseudo-samples.
    pub mode: SubsampleMode,
    /// Deterministic base seed.
    pub seed: u64,
    /// New output directory; existing paths are rejected.
    pub out_dir: PathBuf,
}

/// Summary of a published pseudo-sample bundle.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubsampleResult {
    /// Published output directory.
    pub out_dir: PathBuf,
    /// Number of source read molecules available for sampling.
    pub available_reads: usize,
    /// Generated pseudo-sample identifiers.
    pub samples: Vec<String>,
    /// Number of source assay rows split for every pseudo-sample.
    pub assays: usize,
}

#[derive(Serialize)]
struct BundleProvenance {
    format_version: u32,
    tool: &'static str,
    tool_version: &'static str,
    source_sample: String,
    parent_group: Option<String>,
    source_read_set_sha256: String,
    sample_prefix: String,
    mode: String,
    seed: u64,
    replicates: usize,
    reads_per_sample: usize,
    sampling_unit: &'static str,
    biological_replicates: bool,
    intended_use: &'static str,
    source_manifest: String,
    source_read_to_isoform: String,
    source_mod_manifest: String,
}

struct StagedDirectory {
    path: PathBuf,
    armed: bool,
}

impl Drop for StagedDirectory {
    fn drop(&mut self) {
        if self.armed {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}

#[derive(Clone, Debug, Default)]
struct AssaySampleQc {
    observation_rows: usize,
    observation_reads: usize,
    bam_records: Option<usize>,
    bam_primary_reads: Option<usize>,
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn staged_output_directory(out_dir: &Path) -> anyhow::Result<StagedDirectory> {
    if out_dir.file_name().is_none() {
        anyhow::bail!("subsample output directory must have a final path component");
    }
    match fs::symlink_metadata(out_dir) {
        Ok(_) => anyhow::bail!("refusing to overwrite existing output directory {out_dir:?}"),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error).with_context(|| format!("inspect output directory {out_dir:?}"));
        }
    }
    let parent = output_parent(out_dir);
    fs::create_dir_all(parent)
        .with_context(|| format!("create output parent directory {parent:?}"))?;
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before Unix epoch")?
        .as_nanos();
    let staged = parent.join(format!(
        ".trackcluster-mod-subsample-{}-{nonce}",
        std::process::id()
    ));
    fs::create_dir(&staged)
        .with_context(|| format!("create staged output directory {staged:?}"))?;
    Ok(StagedDirectory {
        path: staged,
        armed: true,
    })
}

fn validate_options(options: &SubsampleOptions) -> anyhow::Result<()> {
    if options.replicates == 0 {
        anyhow::bail!("replicates must be at least 1");
    }
    if options.reads_per_sample == 0 {
        anyhow::bail!("reads-per-sample must be at least 1");
    }
    if options.source_sample.trim().is_empty()
        || options.source_sample == "NA"
        || options.source_sample.contains(SAMPLE_DELIM)
        || options.source_sample.chars().any(char::is_control)
    {
        anyhow::bail!("source sample is empty, unsafe, or contains {SAMPLE_DELIM:?}");
    }
    crate::flow::path_key::SafePathComponent::parse(
        "subsample sample prefix",
        &options.sample_prefix,
    )?;
    if options.sample_prefix.contains(SAMPLE_DELIM) {
        anyhow::bail!("sample prefix must not contain {SAMPLE_DELIM:?}");
    }
    Ok(())
}

fn score(seed: u64, replicate: Option<usize>, read_id: &str) -> [u8; 32] {
    let mut digest = Sha256::new();
    digest.update(b"trackcluster-mod-subsample-v1\0");
    digest.update(seed.to_le_bytes());
    if let Some(replicate) = replicate {
        digest.update((replicate as u64).to_le_bytes());
    }
    digest.update(b"\0");
    digest.update(read_id.as_bytes());
    digest.finalize().into()
}

fn read_set_sha256(read_ids: &BTreeSet<String>) -> String {
    let mut digest = Sha256::new();
    for read_id in read_ids {
        digest.update(read_id.as_bytes());
        digest.update(b"\n");
    }
    format!("{:x}", digest.finalize())
}

fn select_reads(
    source_reads: &BTreeSet<String>,
    replicates: usize,
    reads_per_sample: usize,
    seed: u64,
    mode: SubsampleMode,
) -> anyhow::Result<Vec<BTreeSet<String>>> {
    if reads_per_sample > source_reads.len() {
        anyhow::bail!(
            "requested {reads_per_sample} reads per pseudo-sample, but only {} source reads are \
             available",
            source_reads.len()
        );
    }
    if mode == SubsampleMode::Disjoint {
        let total = replicates
            .checked_mul(reads_per_sample)
            .context("replicates * reads-per-sample overflows usize")?;
        if total > source_reads.len() {
            anyhow::bail!(
                "disjoint mode requests {total} reads across {replicates} pseudo-samples, but only \
                 {} source reads are available",
                source_reads.len()
            );
        }
        let mut ranked = source_reads
            .iter()
            .map(|read_id| (score(seed, None, read_id), read_id.as_str()))
            .collect::<Vec<_>>();
        ranked.sort_unstable();
        return Ok((0..replicates)
            .map(|replicate| {
                let start = replicate * reads_per_sample;
                ranked[start..start + reads_per_sample]
                    .iter()
                    .map(|(_, read_id)| (*read_id).to_owned())
                    .collect()
            })
            .collect());
    }

    Ok((0..replicates)
        .map(|replicate| {
            let mut ranked = source_reads
                .iter()
                .map(|read_id| (score(seed, Some(replicate), read_id), read_id.as_str()))
                .collect::<Vec<_>>();
            ranked.sort_unstable();
            ranked
                .into_iter()
                .take(reads_per_sample)
                .map(|(_, read_id)| read_id.to_owned())
                .collect()
        })
        .collect())
}

fn read_source_assignments(
    path: &Path,
    source_sample: &str,
) -> anyhow::Result<BTreeMap<String, String>> {
    let pairs = crate::count::read_read_to_isoform_tsv(path)
        .with_context(|| format!("read unique assignment {path:?}"))?;
    let mut assignments = BTreeMap::new();
    for (read_id, isoform_id) in pairs {
        let (sample, raw_read_id) = split_tagged_read_name(&read_id).with_context(|| {
            format!("read-to-isoform ID {read_id:?} must use the <sample>::<read_id> form")
        })?;
        if sample != source_sample {
            continue;
        }
        match assignments.get(raw_read_id) {
            Some(existing) if existing == &isoform_id => {}
            Some(existing) => {
                anyhow::bail!(
                    "source read {read_id:?} is assigned to both {existing:?} and {isoform_id:?}; \
                     mod-subsample requires unique assignment"
                );
            }
            None => {
                assignments.insert(raw_read_id.to_owned(), isoform_id);
            }
        }
    }
    if assignments.is_empty() {
        anyhow::bail!(
            "read-to-isoform {path:?} has no assignments for source sample {source_sample:?}"
        );
    }
    Ok(assignments)
}

fn read_source_tracks(
    path: &Path,
) -> anyhow::Result<(Vec<crate::model::Transcript>, BTreeSet<String>)> {
    let tracks = crate::io::bed::read_bed12(path)
        .with_context(|| format!("open source reads {path:?}"))?
        .collect::<Result<Vec<_>, _>>()
        .with_context(|| format!("parse source reads {path:?}"))?;
    let mut names = BTreeSet::new();
    for track in &tracks {
        if !names.insert(track.name.clone()) {
            anyhow::bail!(
                "source reads {path:?} contain duplicate read name {:?}",
                track.name
            );
        }
    }
    Ok((tracks, names))
}

fn selected_membership(selected: &[BTreeSet<String>]) -> HashMap<&str, Vec<usize>> {
    let mut membership = HashMap::new();
    for (sample_index, read_ids) in selected.iter().enumerate() {
        for read_id in read_ids {
            membership
                .entry(read_id.as_str())
                .or_insert_with(Vec::new)
                .push(sample_index);
        }
    }
    membership
}

fn write_sample_reads(
    directory: &Path,
    tracks: &[crate::model::Transcript],
    selected: &[BTreeSet<String>],
    sample_names: &[String],
) -> anyhow::Result<Vec<PathBuf>> {
    let paths = sample_names
        .iter()
        .map(|sample| directory.join(format!("{sample}.reads.bed")))
        .collect::<Vec<_>>();
    let mut writers = paths
        .iter()
        .map(|path| {
            File::create(path)
                .map(BufWriter::new)
                .with_context(|| format!("create pseudo-sample reads {path:?}"))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let membership = selected_membership(selected);
    let mut counts = vec![0usize; selected.len()];
    for track in tracks {
        let Some(sample_indices) = membership.get(track.name.as_str()) else {
            continue;
        };
        for &sample_index in sample_indices {
            crate::io::bed::write_bed12_to_writer(
                &mut writers[sample_index],
                std::iter::once(track),
            )
            .with_context(|| format!("write pseudo-sample read {:?}", track.name))?;
            counts[sample_index] += 1;
        }
    }
    for (sample_index, writer) in writers.iter_mut().enumerate() {
        writer.flush().with_context(|| {
            format!(
                "flush pseudo-sample reads for {:?}",
                sample_names[sample_index]
            )
        })?;
        if counts[sample_index] != selected[sample_index].len() {
            anyhow::bail!(
                "pseudo-sample {:?} wrote {} reads, expected {}",
                sample_names[sample_index],
                counts[sample_index],
                selected[sample_index].len()
            );
        }
    }
    Ok(paths)
}

fn write_assignments(
    path: &Path,
    assignments: &BTreeMap<String, String>,
    selected: &[BTreeSet<String>],
    sample_names: &[String],
) -> anyhow::Result<()> {
    let mut pairs = Vec::new();
    for (sample_index, read_ids) in selected.iter().enumerate() {
        for raw_read_id in read_ids {
            if let Some(isoform_id) = assignments.get(raw_read_id) {
                pairs.push((
                    tagged_read_name(&sample_names[sample_index], raw_read_id),
                    isoform_id.clone(),
                ));
            }
        }
    }
    crate::cluster::output::write_read_to_isoform_tsv(path, &pairs)
        .with_context(|| format!("write pseudo-sample assignments {path:?}"))
}

fn write_sample_manifest(path: &Path, sample_names: &[String]) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("create sample manifest {path:?}"))?,
    );
    writeln!(writer, "sample\tgroup\treads")?;
    for sample in sample_names {
        writeln!(writer, "{sample}\t\tsamples/{sample}.reads.bed")?;
    }
    writer
        .flush()
        .with_context(|| format!("flush sample manifest {path:?}"))
}

fn write_sample_provenance(
    path: &Path,
    source_sample: &str,
    parent_group: Option<&str>,
    sample_names: &[String],
    mode: SubsampleMode,
    seed: u64,
) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("create sample provenance {path:?}"))?,
    );
    writeln!(
        writer,
        "sample\tparent_sample\tparent_group\tsample_kind\tmode\tseed"
    )?;
    for sample in sample_names {
        writeln!(
            writer,
            "{sample}\t{source_sample}\t{}\ttechnical_pseudo\t{mode}\t{seed}",
            parent_group.unwrap_or("NA")
        )?;
    }
    writer
        .flush()
        .with_context(|| format!("flush sample provenance {path:?}"))
}

fn write_selected_read_ids(
    path: &Path,
    assignments: &BTreeMap<String, String>,
    selected: &[BTreeSet<String>],
    sample_names: &[String],
) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("create selected-read audit {path:?}"))?,
    );
    writeln!(
        writer,
        "output_sample\tsource_read_id\tisoform_id\tselection_index"
    )?;
    for (sample_index, read_ids) in selected.iter().enumerate() {
        for (selection_index, read_id) in read_ids.iter().enumerate() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                sample_names[sample_index],
                read_id,
                assignments.get(read_id).map(String::as_str).unwrap_or("NA"),
                selection_index + 1
            )?;
        }
    }
    writer
        .flush()
        .with_context(|| format!("flush selected-read audit {path:?}"))
}

fn write_overlap_qc(
    path: &Path,
    selected: &[BTreeSet<String>],
    sample_names: &[String],
) -> anyhow::Result<()> {
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("create overlap QC {path:?}"))?);
    writeln!(
        writer,
        "sample_a\tsample_b\tintersection_reads\tunion_reads\tjaccard"
    )?;
    for left in 0..selected.len() {
        for right in (left + 1)..selected.len() {
            let intersection = selected[left].intersection(&selected[right]).count();
            let union = selected[left].len() + selected[right].len() - intersection;
            let jaccard = if union == 0 {
                0.0
            } else {
                intersection as f64 / union as f64
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                sample_names[left], sample_names[right], intersection, union, jaccard
            )?;
        }
    }
    writer
        .flush()
        .with_context(|| format!("flush overlap QC {path:?}"))
}

fn optional_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

fn collect_bundle_files(
    root: &Path,
    directory: &Path,
    files: &mut Vec<PathBuf>,
) -> anyhow::Result<()> {
    for entry in fs::read_dir(directory)
        .with_context(|| format!("list staged bundle directory {directory:?}"))?
    {
        let entry = entry.with_context(|| format!("read staged bundle entry in {directory:?}"))?;
        let path = entry.path();
        let file_type = entry
            .file_type()
            .with_context(|| format!("inspect staged bundle entry {path:?}"))?;
        if file_type.is_dir() {
            collect_bundle_files(root, &path, files)?;
        } else if file_type.is_file() {
            let relative = path.strip_prefix(root).expect("bundle entry is below root");
            if relative != Path::new("SHA256SUMS") {
                files.push(relative.to_path_buf());
            }
        } else {
            anyhow::bail!("staged bundle contains a non-file, non-directory entry {path:?}");
        }
    }
    Ok(())
}

fn write_bundle_checksums(root: &Path) -> anyhow::Result<()> {
    let mut files = Vec::new();
    collect_bundle_files(root, root, &mut files)?;
    files.sort();
    let checksum_path = root.join("SHA256SUMS");
    let mut output = BufWriter::new(
        File::create(&checksum_path)
            .with_context(|| format!("create bundle checksums {checksum_path:?}"))?,
    );
    let mut buffer = vec![0u8; 1024 * 1024];
    for relative in files {
        let path = root.join(&relative);
        let mut input = BufReader::new(
            File::open(&path).with_context(|| format!("open bundle artifact {path:?}"))?,
        );
        let mut digest = Sha256::new();
        loop {
            let read = input
                .read(&mut buffer)
                .with_context(|| format!("hash bundle artifact {path:?}"))?;
            if read == 0 {
                break;
            }
            digest.update(&buffer[..read]);
        }
        writeln!(output, "{:x}  {}", digest.finalize(), relative.display())?;
    }
    output
        .flush()
        .with_context(|| format!("flush bundle checksums {checksum_path:?}"))
}

/// Build and atomically publish a deterministic pseudo-sample input bundle.
pub fn create_subsample_bundle(options: &SubsampleOptions) -> anyhow::Result<SubsampleResult> {
    validate_options(options)?;

    let samples = crate::io::manifest::read_manifest_tsv(&options.manifest)?;
    let source_row = samples
        .iter()
        .find(|row| row.sample == options.source_sample)
        .with_context(|| {
            format!(
                "sample manifest {:?} has no source sample {:?}",
                options.manifest, options.source_sample
            )
        })?;
    let sample_names = (1..=options.replicates)
        .map(|index| format!("{}_{index:03}", options.sample_prefix))
        .collect::<Vec<_>>();
    let existing_samples = samples
        .iter()
        .map(|row| row.sample.as_str())
        .collect::<HashSet<_>>();
    if let Some(collision) = sample_names
        .iter()
        .find(|sample| existing_samples.contains(sample.as_str()))
    {
        anyhow::bail!("generated pseudo-sample {collision:?} collides with an input sample");
    }

    let assignments = read_source_assignments(&options.read_to_isoform, &options.source_sample)?;
    let (source_tracks, source_track_names) = read_source_tracks(&source_row.reads)?;
    if let Some(missing) = assignments
        .keys()
        .find(|read_id| !source_track_names.contains(read_id.as_str()))
    {
        anyhow::bail!(
            "source assignment {:?} is missing from source reads {:?}",
            tagged_read_name(&options.source_sample, missing),
            source_row.reads
        );
    }
    let selected = select_reads(
        &source_track_names,
        options.replicates,
        options.reads_per_sample,
        options.seed,
        options.mode,
    )?;
    let membership = selected_membership(&selected);

    let mod_rows = crate::io::mod_manifest::read_mod_manifest_tsv(&options.mod_manifest)?
        .into_iter()
        .filter(|row| row.sample == options.source_sample)
        .collect::<Vec<_>>();
    if mod_rows.is_empty() {
        anyhow::bail!(
            "modification manifest {:?} has no rows for source sample {:?}",
            options.mod_manifest,
            options.source_sample
        );
    }

    let mut stage = staged_output_directory(&options.out_dir)?;
    let samples_dir = stage.path.join("samples");
    let observations_dir = stage.path.join("observations");
    let assays_dir = stage.path.join("assays");
    let coverage_dir = stage.path.join("coverage");
    for directory in [
        samples_dir.as_path(),
        observations_dir.as_path(),
        assays_dir.as_path(),
        coverage_dir.as_path(),
    ] {
        fs::create_dir(directory)
            .with_context(|| format!("create staged bundle directory {directory:?}"))?;
    }

    write_sample_reads(&samples_dir, &source_tracks, &selected, &sample_names)?;
    write_assignments(
        &stage.path.join("read_to_isoform.unique.tsv"),
        &assignments,
        &selected,
        &sample_names,
    )?;
    write_sample_manifest(&stage.path.join("samples.tsv"), &sample_names)?;
    write_sample_provenance(
        &stage.path.join("sample_provenance.tsv"),
        &options.source_sample,
        source_row.group.as_deref(),
        &sample_names,
        options.mode,
        options.seed,
    )?;
    write_selected_read_ids(
        &stage.path.join("subsample_read_ids.tsv"),
        &assignments,
        &selected,
        &sample_names,
    )?;
    write_overlap_qc(&stage.path.join("overlap_qc.tsv"), &selected, &sample_names)?;

    let mut all_qc = vec![vec![AssaySampleQc::default(); sample_names.len()]; mod_rows.len()];
    let mut assay_metadata_paths = Vec::with_capacity(mod_rows.len());
    let mut observation_paths = vec![Vec::new(); mod_rows.len()];
    let mut coverage_paths = vec![Vec::new(); mod_rows.len()];

    for (assay_index, mod_row) in mod_rows.iter().enumerate() {
        let metadata = crate::io::mod_calls::read_assay_metadata(&mod_row.assay_metadata)?;
        if metadata.assay_id != mod_row.assay_id {
            anyhow::bail!(
                "assay metadata {:?} declares assay {:?}, expected {:?}",
                mod_row.assay_metadata,
                metadata.assay_id,
                mod_row.assay_id
            );
        }
        let mut derived_metadata = metadata;
        derived_metadata.read_id_mapping = format!(
            "{};trackcluster_mod_subsample_v1",
            derived_metadata.read_id_mapping
        );
        derived_metadata
            .source_files
            .push(mod_row.observations.display().to_string());
        if let Some(path) = mod_row.coverage_bam.as_ref() {
            derived_metadata
                .source_files
                .push(path.display().to_string());
        }
        let metadata_relative = PathBuf::from(format!("assays/assay_{:03}.json", assay_index + 1));
        let metadata_path = stage.path.join(&metadata_relative);
        write_assay_metadata_to_writer(
            BufWriter::new(
                File::create(&metadata_path)
                    .with_context(|| format!("create assay metadata {metadata_path:?}"))?,
            ),
            &derived_metadata,
        )?;
        assay_metadata_paths.push(metadata_relative);

        let assay_observation_paths = sample_names
            .iter()
            .map(|sample| {
                PathBuf::from(format!(
                    "observations/{sample}.assay_{:03}.observations.tsv",
                    assay_index + 1
                ))
            })
            .collect::<Vec<_>>();
        let mut observation_writers = assay_observation_paths
            .iter()
            .map(|relative| {
                let path = stage.path.join(relative);
                let file = File::create(&path)
                    .with_context(|| format!("create pseudo-sample observations {path:?}"))?;
                new_observation_writer(BufWriter::new(file))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        let mut observation_reads = vec![HashSet::new(); sample_names.len()];
        scan_canonical_observations_tsv(&mod_row.observations, |observation| {
            if observation.key.assay_id != mod_row.assay_id
                || observation.key.sample != options.source_sample
            {
                anyhow::bail!(
                    "observation key ({:?}, {:?}) does not match mod manifest ({:?}, {:?})",
                    observation.key.sample,
                    observation.key.assay_id,
                    options.source_sample,
                    mod_row.assay_id
                );
            }
            let (sample, raw_read_id) = split_tagged_read_name(&observation.key.read_id)
                .with_context(|| {
                    format!(
                        "observation read {:?} must use <sample>::<read_id>",
                        observation.key.read_id
                    )
                })?;
            if sample != options.source_sample {
                anyhow::bail!(
                    "observation read {:?} has source prefix {sample:?}, expected {:?}",
                    observation.key.read_id,
                    options.source_sample
                );
            }
            let Some(sample_indices) = membership.get(raw_read_id) else {
                return Ok(());
            };
            for &sample_index in sample_indices {
                let mut derived = observation.clone();
                derived.key.sample = sample_names[sample_index].clone();
                derived.key.read_id = tagged_read_name(&sample_names[sample_index], raw_read_id);
                write_observation_record(&mut observation_writers[sample_index], &derived)?;
                all_qc[assay_index][sample_index].observation_rows += 1;
                observation_reads[sample_index].insert(raw_read_id.to_owned());
            }
            Ok(())
        })?;
        for (sample_index, writer) in observation_writers.iter_mut().enumerate() {
            writer.flush().with_context(|| {
                format!(
                    "flush observations for assay {:?}, sample {:?}",
                    mod_row.assay_id, sample_names[sample_index]
                )
            })?;
            all_qc[assay_index][sample_index].observation_reads =
                observation_reads[sample_index].len();
        }
        observation_paths[assay_index] = assay_observation_paths;

        if let Some(source_bam) = mod_row.coverage_bam.as_ref() {
            let assay_coverage_paths = sample_names
                .iter()
                .map(|sample| {
                    PathBuf::from(format!(
                        "coverage/{sample}.assay_{:03}.bam",
                        assay_index + 1
                    ))
                })
                .collect::<Vec<_>>();
            let absolute_paths = assay_coverage_paths
                .iter()
                .map(|relative| stage.path.join(relative))
                .collect::<Vec<_>>();
            let targets = absolute_paths
                .iter()
                .zip(selected.iter())
                .map(|(path, read_names)| BamSubsetTarget {
                    path: path.as_path(),
                    read_names,
                })
                .collect::<Vec<_>>();
            let bam_qc = split_bam_by_query_name(source_bam, &targets)?;
            for (sample_index, target_qc) in bam_qc.targets.iter().enumerate() {
                all_qc[assay_index][sample_index].bam_records = Some(target_qc.written_records);
                all_qc[assay_index][sample_index].bam_primary_reads = Some(target_qc.primary_reads);
            }
            coverage_paths[assay_index] = assay_coverage_paths.into_iter().map(Some).collect();
        } else {
            coverage_paths[assay_index] = vec![None; sample_names.len()];
        }
    }

    let mut mod_manifest = BufWriter::new(
        File::create(stage.path.join("mod_samples.tsv"))
            .context("create pseudo-sample modification manifest")?,
    );
    writeln!(
        mod_manifest,
        "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam"
    )?;
    for (assay_index, mod_row) in mod_rows.iter().enumerate() {
        for (sample_index, sample) in sample_names.iter().enumerate() {
            let coverage = coverage_paths[assay_index][sample_index]
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| "NA".to_owned());
            writeln!(
                mod_manifest,
                "{sample}\t{}\t{}\t{}\t{coverage}",
                mod_row.assay_id,
                observation_paths[assay_index][sample_index].display(),
                assay_metadata_paths[assay_index].display(),
            )?;
        }
    }
    mod_manifest
        .flush()
        .context("flush pseudo-sample modification manifest")?;

    let mut qc_writer = BufWriter::new(
        File::create(stage.path.join("subsample_qc.tsv"))
            .context("create pseudo-sample QC table")?,
    );
    writeln!(
        qc_writer,
        "source_sample\toutput_sample\tassay_id\tmode\tseed\tavailable_assigned_reads\t\
         available_reads\trequested_reads\tselected_reads\tselected_assigned_reads\t\
         observation_reads\tobservation_rows\tbam_records\tbam_primary_reads"
    )?;
    for (assay_index, mod_row) in mod_rows.iter().enumerate() {
        for (sample_index, sample) in sample_names.iter().enumerate() {
            let qc = &all_qc[assay_index][sample_index];
            writeln!(
                qc_writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                options.source_sample,
                sample,
                mod_row.assay_id,
                options.mode,
                options.seed,
                assignments.len(),
                source_track_names.len(),
                options.reads_per_sample,
                selected[sample_index].len(),
                selected[sample_index]
                    .iter()
                    .filter(|read_id| assignments.contains_key(read_id.as_str()))
                    .count(),
                qc.observation_reads,
                qc.observation_rows,
                optional_usize(qc.bam_records),
                optional_usize(qc.bam_primary_reads),
            )?;
        }
    }
    qc_writer.flush().context("flush pseudo-sample QC table")?;

    let provenance = BundleProvenance {
        format_version: 1,
        tool: "trackcluster mod-subsample",
        tool_version: env!("CARGO_PKG_VERSION"),
        source_sample: options.source_sample.clone(),
        parent_group: source_row.group.clone(),
        source_read_set_sha256: read_set_sha256(&source_track_names),
        sample_prefix: options.sample_prefix.clone(),
        mode: options.mode.to_string(),
        seed: options.seed,
        replicates: options.replicates,
        reads_per_sample: options.reads_per_sample,
        sampling_unit: "primary_gene_read_molecule",
        biological_replicates: false,
        intended_use: "technical_coverage_stability_only",
        source_manifest: options.manifest.display().to_string(),
        source_read_to_isoform: options.read_to_isoform.display().to_string(),
        source_mod_manifest: options.mod_manifest.display().to_string(),
    };
    let provenance_path = stage.path.join("subsample_provenance.json");
    let mut provenance_writer = BufWriter::new(
        File::create(&provenance_path)
            .with_context(|| format!("create subsample provenance {provenance_path:?}"))?,
    );
    serde_json::to_writer_pretty(&mut provenance_writer, &provenance)
        .context("serialize subsample provenance")?;
    provenance_writer.write_all(b"\n")?;
    provenance_writer.flush()?;
    write_bundle_checksums(&stage.path)?;

    match fs::symlink_metadata(&options.out_dir) {
        Ok(_) => anyhow::bail!(
            "output directory {:?} appeared while the bundle was being built",
            options.out_dir
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error)
                .with_context(|| format!("recheck output directory {:?}", options.out_dir));
        }
    }
    fs::rename(&stage.path, &options.out_dir).with_context(|| {
        format!(
            "publish staged pseudo-sample bundle {:?} to {:?}",
            stage.path, options.out_dir
        )
    })?;
    stage.armed = false;

    Ok(SubsampleResult {
        out_dir: options.out_dir.clone(),
        available_reads: source_track_names.len(),
        samples: sample_names,
        assays: mod_rows.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn universe(size: usize) -> BTreeSet<String> {
        (0..size).map(|index| format!("read-{index:04}")).collect()
    }

    #[test]
    fn hash_selection_is_deterministic_disjoint_and_order_independent() {
        let reads = universe(100);
        let selected =
            select_reads(&reads, 4, 20, 17, SubsampleMode::Disjoint).expect("select reads");
        let repeated =
            select_reads(&reads, 4, 20, 17, SubsampleMode::Disjoint).expect("repeat selection");
        assert_eq!(selected, repeated);
        assert!(selected
            .iter()
            .all(|sample| sample.len() == 20 && sample.is_subset(&reads)));
        for left in 0..selected.len() {
            for right in (left + 1)..selected.len() {
                assert!(selected[left].is_disjoint(&selected[right]));
            }
        }
        let different_seed =
            select_reads(&reads, 4, 20, 18, SubsampleMode::Disjoint).expect("different seed");
        assert_ne!(selected, different_seed);
    }

    #[test]
    fn independent_selection_allows_overlap_but_never_duplicates_within_sample() {
        let reads = universe(40);
        let selected =
            select_reads(&reads, 4, 20, 9, SubsampleMode::Independent).expect("select reads");
        assert!(selected.iter().all(|sample| sample.len() == 20));
        assert!(selected.iter().all(|sample| sample.is_subset(&reads)));
        assert!(selected.iter().enumerate().any(|(left, sample)| selected
            .iter()
            .skip(left + 1)
            .any(|other| !sample.is_disjoint(other))));
    }

    #[test]
    fn disjoint_selection_rejects_an_oversubscribed_source_universe() {
        let error = select_reads(&universe(10), 3, 4, 1, SubsampleMode::Disjoint)
            .expect_err("oversubscribed disjoint selection must fail");
        assert!(error.to_string().contains("requests 12 reads"));
    }
}
