//! Transactional publication of flow-integrated modification output sets.

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use serde::{Deserialize, Serialize};

use crate::flow::artifact_manifest::{
    atomic_copy, atomic_write_with, invalidate_completion_manifest, sha256_file, InputArtifact,
    ToolIdentity,
};
use crate::modification::EligibilityProfile;

const GENERATION_SCHEMA_VERSION: u32 = 1;
static RUN_SEQUENCE: AtomicU64 = AtomicU64::new(0);

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

/// Effective modification settings recorded in a generation manifest.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct GenerationOptions {
    pub(crate) analysis_thresholds: BTreeMap<String, f64>,
    pub(crate) eligibility_profile: String,
    pub(crate) min_covering: u64,
    pub(crate) min_callable: u64,
    pub(crate) min_candidate_rate: f64,
    pub(crate) min_callable_rate: f64,
    pub(crate) min_read_join_rate: f64,
    pub(crate) allow_low_global_join: bool,
    pub(crate) reference_validation: bool,
    pub(crate) contrast_requested: bool,
}

impl GenerationOptions {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        analysis_thresholds: BTreeMap<String, f64>,
        eligibility_profile: EligibilityProfile,
        min_covering: u64,
        min_callable: u64,
        min_candidate_rate: f64,
        min_callable_rate: f64,
        min_read_join_rate: f64,
        allow_low_global_join: bool,
        reference_validation: bool,
        contrast_requested: bool,
    ) -> Self {
        Self {
            analysis_thresholds,
            eligibility_profile: eligibility_profile.to_string(),
            min_covering,
            min_callable,
            min_candidate_rate,
            min_callable_rate,
            min_read_join_rate,
            allow_low_global_join,
            reference_validation,
            contrast_requested,
        }
    }
}

/// Paths reserved for one unpublished modification generation.
#[derive(Clone, Debug)]
pub(crate) struct GenerationPaths {
    pub(crate) run_id: String,
    pub(crate) directory: PathBuf,
    pub(crate) join_qc: PathBuf,
    pub(crate) site_join_qc: PathBuf,
    pub(crate) sites: PathBuf,
    pub(crate) design: PathBuf,
    pub(crate) contrasts: PathBuf,
    pub(crate) manifest: PathBuf,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
struct GenerationOutput {
    role: String,
    path: String,
    sha256: String,
    bytes: u64,
}

impl GenerationOutput {
    fn from_file(role: &str, path: &Path) -> anyhow::Result<Self> {
        let metadata = fs::symlink_metadata(path)
            .with_context(|| format!("inspect modification generation output {path:?}"))?;
        if metadata.file_type().is_symlink() || !metadata.is_file() {
            anyhow::bail!("modification generation output is not a regular file: {path:?}");
        }
        let name = path
            .file_name()
            .and_then(|value| value.to_str())
            .with_context(|| format!("generation output has no UTF-8 file name: {path:?}"))?;
        Ok(Self {
            role: role.to_owned(),
            path: name.to_owned(),
            sha256: sha256_file(path)?,
            bytes: metadata.len(),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
struct GenerationManifest {
    schema_version: u32,
    complete: bool,
    run_id: String,
    tool: ToolIdentity,
    inputs: Vec<InputArtifact>,
    options: GenerationOptions,
    outputs: Vec<GenerationOutput>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
struct CurrentPointer {
    schema_version: u32,
    complete: bool,
    run_id: String,
    generation: String,
    manifest_sha256: String,
}

pub(crate) fn generations_root(prefix: &Path) -> PathBuf {
    append_suffix(prefix, ".mod.generations")
}

pub(crate) fn current_pointer_path(prefix: &Path) -> PathBuf {
    append_suffix(prefix, ".mod.current.json")
}

/// Mark an output prefix as flow-managed, even before a concrete generation is reserved.
pub(crate) fn ensure_managed(prefix: &Path) -> anyhow::Result<()> {
    let root = generations_root(prefix);
    fs::create_dir_all(&root)
        .with_context(|| format!("create modification generations directory {root:?}"))?;
    let root_metadata = fs::symlink_metadata(&root)
        .with_context(|| format!("inspect modification generations directory {root:?}"))?;
    if root_metadata.file_type().is_symlink() || !root_metadata.is_dir() {
        anyhow::bail!("modification generations path is not a regular directory: {root:?}");
    }
    Ok(())
}

/// Keep standalone aggregation from silently replacing flow-managed compatibility files.
pub(crate) fn ensure_standalone_prefix(prefix: &Path) -> anyhow::Result<()> {
    for path in [generations_root(prefix), current_pointer_path(prefix)] {
        match fs::symlink_metadata(&path) {
            Ok(_) => anyhow::bail!(
                "modification output prefix {prefix:?} is flow-managed; rerun flow or choose a different standalone --out prefix"
            ),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                return Err(error).with_context(|| format!("inspect modification state {path:?}"));
            }
        }
    }
    Ok(())
}

/// Invalidate the authoritative pointer before any core or modification output is replaced.
pub(crate) fn invalidate_current(prefix: &Path) -> anyhow::Result<()> {
    let path = current_pointer_path(prefix);
    match fs::symlink_metadata(&path) {
        Ok(metadata) => {
            if !metadata.file_type().is_file() {
                anyhow::bail!(
                    "refusing to invalidate modification current pointer {path:?}: expected a regular file"
                );
            }
            invalidate_completion_manifest(&path)?;
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error)
                .with_context(|| format!("inspect modification current pointer {path:?}"));
        }
    }
    Ok(())
}

/// Reserve a unique, unpublished generation directory.
pub(crate) fn begin(prefix: &Path) -> anyhow::Result<GenerationPaths> {
    let root = generations_root(prefix);
    ensure_managed(prefix)?;
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system clock is before the Unix epoch")?
        .as_nanos();
    for _ in 0..1000 {
        let sequence = RUN_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let run_id = format!("{nanos}-{}-{sequence}", std::process::id());
        let directory = root.join(&run_id);
        match fs::create_dir(&directory) {
            Ok(()) => {
                return Ok(GenerationPaths {
                    run_id,
                    join_qc: directory.join("mod_join_qc.tsv"),
                    site_join_qc: directory.join("mod_site_join_qc.tsv"),
                    sites: directory.join("isoform_mod_sites.tsv"),
                    design: directory.join("isoform_mod_design.tsv"),
                    contrasts: directory.join("isoform_mod_contrasts.tsv"),
                    manifest: directory.join("manifest.json"),
                    directory,
                });
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("reserve modification generation {directory:?}"));
            }
        }
    }
    anyhow::bail!("could not reserve a unique modification generation directory")
}

fn fai_path(reference: &Path) -> PathBuf {
    append_suffix(reference, ".fai")
}

/// Hash every direct and nested input that can affect a modification generation.
pub(crate) fn collect_inputs(
    sample_manifest: &Path,
    isoforms: &Path,
    read_to_isoform: &Path,
    mod_manifest: &Path,
    reference_fasta: Option<&Path>,
    contrasts: Option<&Path>,
) -> anyhow::Result<Vec<InputArtifact>> {
    let mut inputs = vec![
        InputArtifact::from_file("sample_manifest", sample_manifest)?,
        InputArtifact::from_file("isoform_bed", isoforms)?,
        InputArtifact::from_file("read_to_isoform", read_to_isoform)?,
        InputArtifact::from_file("mod_manifest", mod_manifest)?,
    ];
    if let Some(reference) = reference_fasta {
        inputs.push(InputArtifact::from_file("reference_fasta", reference)?);
        inputs.push(InputArtifact::from_file(
            "reference_fai",
            &fai_path(reference),
        )?);
    }
    if let Some(contrasts) = contrasts {
        inputs.push(InputArtifact::from_file(
            "contrast_specification",
            contrasts,
        )?);
    }
    for row in crate::io::mod_manifest::read_mod_manifest_tsv(mod_manifest)? {
        let prefix = format!("{}:{}", row.sample, row.assay_id);
        inputs.push(InputArtifact::from_file(
            &format!("mod_observations:{prefix}"),
            &row.observations,
        )?);
        inputs.push(InputArtifact::from_file(
            &format!("assay_metadata:{prefix}"),
            &row.assay_metadata,
        )?);
        if let Some(coverage) = row.coverage_bam {
            inputs.push(InputArtifact::from_file(
                &format!("coverage_bam:{prefix}"),
                &coverage,
            )?);
        }
    }
    inputs.sort_by(|left, right| {
        left.role
            .cmp(&right.role)
            .then_with(|| left.path.cmp(&right.path))
    });
    Ok(inputs)
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> anyhow::Result<()> {
    atomic_write_with(path, |writer| {
        let mut buffered = std::io::BufWriter::new(writer);
        serde_json::to_writer_pretty(&mut buffered, value)
            .with_context(|| format!("serialize modification generation JSON {path:?}"))?;
        use std::io::Write as _;
        buffered.write_all(b"\n")?;
        buffered.flush()?;
        Ok(())
    })
}

fn flat_outputs(prefix: &Path) -> [(PathBuf, &'static str); 5] {
    [
        (append_suffix(prefix, ".mod_join_qc.tsv"), "mod_join_qc"),
        (
            append_suffix(prefix, ".mod_site_join_qc.tsv"),
            "mod_site_join_qc",
        ),
        (
            append_suffix(prefix, ".isoform_mod_sites.tsv"),
            "isoform_mod_sites",
        ),
        (
            append_suffix(prefix, ".isoform_mod_design.tsv"),
            "isoform_mod_design",
        ),
        (
            append_suffix(prefix, ".isoform_mod_contrasts.tsv"),
            "isoform_mod_contrasts",
        ),
    ]
}

fn remove_optional_regular_file(path: &Path) -> anyhow::Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            if !metadata.file_type().is_file() {
                anyhow::bail!("refusing to remove stale output {path:?}: expected a regular file");
            }
            fs::remove_file(path).with_context(|| format!("remove stale output {path:?}"))
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("inspect stale output {path:?}")),
    }
}

/// Validate, record, and publish a generation; the current pointer is written last.
pub(crate) fn publish(
    prefix: &Path,
    paths: &GenerationPaths,
    inputs: Vec<InputArtifact>,
    options: GenerationOptions,
    has_contrasts: bool,
) -> anyhow::Result<()> {
    let mut outputs = vec![
        GenerationOutput::from_file("mod_join_qc", &paths.join_qc)?,
        GenerationOutput::from_file("mod_site_join_qc", &paths.site_join_qc)?,
        GenerationOutput::from_file("isoform_mod_sites", &paths.sites)?,
        GenerationOutput::from_file("isoform_mod_design", &paths.design)?,
    ];
    if has_contrasts {
        outputs.push(GenerationOutput::from_file(
            "isoform_mod_contrasts",
            &paths.contrasts,
        )?);
    }
    outputs.sort_by(|left, right| left.role.cmp(&right.role));
    let manifest = GenerationManifest {
        schema_version: GENERATION_SCHEMA_VERSION,
        complete: true,
        run_id: paths.run_id.clone(),
        tool: ToolIdentity::current(),
        inputs,
        options,
        outputs,
    };
    write_json(&paths.manifest, &manifest)?;

    let generation_by_role = BTreeMap::from([
        ("mod_join_qc", paths.join_qc.as_path()),
        ("mod_site_join_qc", paths.site_join_qc.as_path()),
        ("isoform_mod_sites", paths.sites.as_path()),
        ("isoform_mod_design", paths.design.as_path()),
        ("isoform_mod_contrasts", paths.contrasts.as_path()),
    ]);
    for (flat, role) in flat_outputs(prefix) {
        if role == "isoform_mod_contrasts" && !has_contrasts {
            remove_optional_regular_file(&flat)?;
        } else {
            atomic_copy(generation_by_role[role], &flat)?;
        }
    }

    let pointer = CurrentPointer {
        schema_version: GENERATION_SCHEMA_VERSION,
        complete: true,
        run_id: paths.run_id.clone(),
        generation: paths
            .directory
            .strip_prefix(generations_root(prefix))
            .unwrap_or(paths.directory.as_path())
            .to_string_lossy()
            .into_owned(),
        manifest_sha256: sha256_file(&paths.manifest)?,
    };
    write_json(&current_pointer_path(prefix), &pointer)
}

fn managed_prefix(path: &Path, suffix: &str) -> Option<PathBuf> {
    let name = path.file_name()?.to_str()?;
    let prefix_name = name.strip_suffix(suffix)?;
    if prefix_name.is_empty() {
        return None;
    }
    Some(
        path.parent()
            .unwrap_or_else(|| Path::new("."))
            .join(prefix_name),
    )
}

/// Reject a stale or tampered flat compatibility file when it belongs to a flow generation.
/// Standalone `mod-aggregate` outputs, which have no generations directory, remain supported.
pub(crate) fn validate_current_flat(path: &Path, suffix: &str, role: &str) -> anyhow::Result<()> {
    let Some(prefix) = managed_prefix(path, suffix) else {
        return Ok(());
    };
    let root = generations_root(&prefix);
    let pointer_path = current_pointer_path(&prefix);
    let managed = root.exists() || pointer_path.exists();
    if !managed {
        return Ok(());
    }
    let pointer_metadata = fs::symlink_metadata(&pointer_path).with_context(|| {
        format!(
            "flow-managed modification output {path:?} has no valid current pointer; it is stale and must not be used"
        )
    })?;
    if pointer_metadata.file_type().is_symlink() || !pointer_metadata.is_file() {
        anyhow::bail!("modification current pointer is not a regular file: {pointer_path:?}");
    }
    let pointer: CurrentPointer = serde_json::from_slice(
        &fs::read(&pointer_path)
            .with_context(|| format!("read modification current pointer {pointer_path:?}"))?,
    )
    .with_context(|| format!("parse modification current pointer {pointer_path:?}"))?;
    if pointer.schema_version != GENERATION_SCHEMA_VERSION || !pointer.complete {
        anyhow::bail!("modification current pointer is incomplete or uses an unsupported schema");
    }
    let generation_component = Path::new(&pointer.generation);
    if generation_component.components().count() != 1
        || generation_component
            .file_name()
            .and_then(|value| value.to_str())
            != Some(pointer.run_id.as_str())
    {
        anyhow::bail!("modification current pointer contains an invalid generation path");
    }
    let directory = root.join(&pointer.generation);
    let manifest_path = directory.join("manifest.json");
    if sha256_file(&manifest_path)? != pointer.manifest_sha256 {
        anyhow::bail!("modification generation manifest hash does not match the current pointer");
    }
    let manifest: GenerationManifest = serde_json::from_slice(
        &fs::read(&manifest_path)
            .with_context(|| format!("read modification generation manifest {manifest_path:?}"))?,
    )
    .with_context(|| format!("parse modification generation manifest {manifest_path:?}"))?;
    if manifest.schema_version != GENERATION_SCHEMA_VERSION
        || !manifest.complete
        || manifest.run_id != pointer.run_id
    {
        anyhow::bail!("modification generation manifest is incomplete or inconsistent");
    }
    let output = manifest
        .outputs
        .iter()
        .find(|output| output.role == role)
        .with_context(|| format!("current modification generation has no {role} output"))?;
    let generation_file = directory.join(&output.path);
    let generation_hash = sha256_file(&generation_file)?;
    if generation_hash != output.sha256 {
        anyhow::bail!("current modification generation {role} output failed hash validation");
    }
    if sha256_file(path)? != generation_hash {
        anyhow::bail!(
            "flat modification output {path:?} does not match the current generation and must not be used"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_root(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("test clock is before Unix epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-mod-generation-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir(&path).expect("create generation test directory");
        path
    }

    #[test]
    fn suffixes_preserve_dotted_prefixes() {
        let prefix = Path::new("out/result.v1");
        assert_eq!(
            current_pointer_path(prefix),
            PathBuf::from("out/result.v1.mod.current.json")
        );
        assert_eq!(
            generations_root(prefix),
            PathBuf::from("out/result.v1.mod.generations")
        );
    }

    #[test]
    fn standalone_prefix_cannot_overlap_flow_managed_state() {
        let root = test_root("standalone-prefix");
        let prefix = root.join("result");
        ensure_standalone_prefix(&prefix).expect("unused prefix is standalone-safe");
        ensure_managed(&prefix).expect("mark prefix flow-managed");
        let error = ensure_standalone_prefix(&prefix)
            .expect_err("flow-managed prefix must reject standalone aggregation");
        assert!(error.to_string().contains("flow-managed"), "{error:#}");
        fs::remove_dir_all(&root).expect("remove generation test directory");
    }

    #[test]
    fn failure_after_first_flat_output_never_publishes_a_current_pointer() {
        let root = test_root("partial-publish");
        let prefix = root.join("result");
        let paths = begin(&prefix).expect("reserve generation");
        for output in [
            &paths.join_qc,
            &paths.site_join_qc,
            &paths.sites,
            &paths.design,
        ] {
            fs::write(output, b"header\n").expect("write generation output");
        }

        let blocked_output = append_suffix(&prefix, ".mod_site_join_qc.tsv");
        fs::create_dir(&blocked_output).expect("block second flat output");
        let options = GenerationOptions::new(
            BTreeMap::new(),
            EligibilityProfile::Exploratory,
            20,
            1,
            0.8,
            0.8,
            0.9,
            false,
            false,
            false,
        );
        let error = publish(&prefix, &paths, Vec::new(), options, false)
            .expect_err("second flat output must fail");
        assert!(!error.to_string().is_empty());

        let first_flat = append_suffix(&prefix, ".mod_join_qc.tsv");
        assert!(first_flat.is_file(), "first compatibility file was written");
        assert!(!current_pointer_path(&prefix).exists());
        let stale = validate_current_flat(&first_flat, ".mod_join_qc.tsv", "mod_join_qc")
            .expect_err("partial compatibility output must be stale");
        assert!(stale.to_string().contains("stale"), "{stale:#}");

        fs::remove_dir_all(&root).expect("remove generation test directory");
    }
}
