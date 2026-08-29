use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Context;
#[cfg(not(unix))]
use same_file::Handle;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub(crate) const MANIFEST_SCHEMA_VERSION: u32 = 3;
pub(crate) const MANIFEST_FILE_NAME: &str = "run.json";

static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct ToolIdentity {
    pub package_name: String,
    pub package_version: String,
    pub git_commit: String,
    #[serde(default = "legacy_source_fingerprint")]
    pub source_fingerprint: String,
}

fn legacy_source_fingerprint() -> String {
    "unknown".to_owned()
}

impl ToolIdentity {
    pub(crate) fn current() -> Self {
        Self {
            package_name: env!("CARGO_PKG_NAME").to_owned(),
            package_version: env!("CARGO_PKG_VERSION").to_owned(),
            git_commit: option_env!("TRACKCLUSTER_GIT_COMMIT")
                .unwrap_or("unknown")
                .to_owned(),
            source_fingerprint: option_env!("TRACKCLUSTER_SOURCE_FINGERPRINT")
                .unwrap_or("unknown")
                .to_owned(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct InputArtifact {
    pub role: String,
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
}

impl InputArtifact {
    pub(crate) fn from_file(role: &str, path: &Path) -> anyhow::Result<Self> {
        let metadata = fs::symlink_metadata(path)
            .with_context(|| format!("inspect manifest input {path:?}"))?;
        if metadata.file_type().is_symlink() || !metadata.is_file() {
            anyhow::bail!("manifest input is not a regular file: {path:?}");
        }
        Ok(Self {
            role: role.to_owned(),
            path: path.to_string_lossy().into_owned(),
            sha256: sha256_file(path)?,
            bytes: metadata.len(),
        })
    }
}

/// Every effective option that can affect clustering or later assignment.
///
/// Floats are safe here because callers validate them as finite before building
/// a request. Keeping this as a concrete structure makes newly added options a
/// deliberate schema change instead of silently falling out of the cache key.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct EffectiveOptions {
    pub cluster_mode: String,
    pub prepare_fraction_read: f64,
    pub prepare_fraction_ref: f64,
    pub sw_score: i64,
    pub batch_size: usize,
    pub batch_rounds: usize,
    pub name2_mode: String,
    pub platform_preset: String,
    pub junction_correction_offset: u32,
    pub junction_correction_min_support: u32,
    pub sl_partial_five_prime_offset: u32,
    pub sl_same_junction_five_prime_offset: u32,
    pub sl_five_prime_cluster_offset: u32,
    pub sl_min_five_prime_cluster_support: usize,
    pub same_junction_three_prime_offset: u32,
    pub three_prime_cluster_offset: u32,
    pub min_three_prime_cluster_support: usize,
    pub overlap_cutoff1: f64,
    pub overlap_cutoff2: f64,
    pub overlap_intron_weight: f64,
    pub assignment_mode: String,
    pub unique_assignment_junction_offset: u32,
    #[serde(default = "legacy_invalid_read_policy")]
    pub invalid_read_policy: String,
    pub downsample_selected: bool,
    pub max_reads_per_gene: usize,
}

fn legacy_invalid_read_policy() -> String {
    "fail".to_owned()
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct RunRequest {
    pub gene: String,
    pub inputs: Vec<InputArtifact>,
    pub options: EffectiveOptions,
    pub tool: ToolIdentity,
    pub seed: u64,
}

impl RunRequest {
    pub(crate) fn fingerprint(&self) -> anyhow::Result<String> {
        let encoded = serde_json::to_vec(self).context("serialize run request fingerprint")?;
        Ok(sha256_bytes(&encoded))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RecordCountKind {
    NonEmptyLines,
    HeaderThenNonEmptyLines,
}

#[derive(Clone, Debug)]
pub(crate) struct OutputSpec {
    pub role: &'static str,
    pub path: PathBuf,
    pub record_count_kind: RecordCountKind,
}

impl OutputSpec {
    pub(crate) fn new(
        role: &'static str,
        path: impl Into<PathBuf>,
        record_count_kind: RecordCountKind,
    ) -> Self {
        Self {
            role,
            path: path.into(),
            record_count_kind,
        }
    }

    fn stored_path(&self) -> anyhow::Result<String> {
        self.path
            .file_name()
            .and_then(|name| name.to_str())
            .map(str::to_owned)
            .with_context(|| format!("manifest output has no UTF-8 file name: {:?}", self.path))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct OutputArtifact {
    pub role: String,
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
    pub records: u64,
}

impl OutputArtifact {
    fn from_spec(spec: &OutputSpec) -> anyhow::Result<Self> {
        Ok(Self {
            role: spec.role.to_owned(),
            path: spec.stored_path()?,
            sha256: sha256_file(&spec.path)?,
            bytes: fs::metadata(&spec.path)
                .with_context(|| format!("stat manifest output {:?}", spec.path))?
                .len(),
            records: count_records(&spec.path, spec.record_count_kind)?,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct RunManifest {
    pub schema_version: u32,
    pub status: String,
    pub request_fingerprint: String,
    pub request: RunRequest,
    pub outputs: Vec<OutputArtifact>,
}

impl RunManifest {
    pub(crate) fn complete(
        request: RunRequest,
        output_specs: &[OutputSpec],
    ) -> anyhow::Result<Self> {
        let request_fingerprint = request.fingerprint()?;
        let mut outputs = output_specs
            .iter()
            .map(OutputArtifact::from_spec)
            .collect::<anyhow::Result<Vec<_>>>()?;
        outputs.sort_by(|left, right| left.role.cmp(&right.role));
        Ok(Self {
            schema_version: MANIFEST_SCHEMA_VERSION,
            status: "complete".to_owned(),
            request_fingerprint,
            request,
            outputs,
        })
    }
}

pub(crate) fn read_run_manifest(path: &Path) -> anyhow::Result<RunManifest> {
    let bytes = fs::read(path).with_context(|| format!("read run manifest {path:?}"))?;
    serde_json::from_slice(&bytes).with_context(|| format!("parse run manifest {path:?}"))
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum CacheDecision {
    Reuse,
    Rebuild(String),
}

impl CacheDecision {
    pub(crate) fn reason(&self) -> &str {
        match self {
            Self::Reuse => "exact_manifest_match",
            Self::Rebuild(reason) => reason,
        }
    }
}

pub(crate) fn assess_cache(
    manifest_path: &Path,
    expected_request: &RunRequest,
    expected_outputs: &[OutputSpec],
) -> CacheDecision {
    match assess_cache_inner(manifest_path, expected_request, expected_outputs) {
        Ok(()) => CacheDecision::Reuse,
        Err(reason) => CacheDecision::Rebuild(reason),
    }
}

/// Validate an already-recorded per-gene completion without requiring the
/// caller to reproduce unrelated clustering options from the original run.
///
/// Count-only execution intentionally consumes a prior completed run. Its
/// trust boundary is therefore the manifest itself: the manifest must be
/// internally self-consistent, belong to the requested gene/mode and current
/// tool, and still describe the exact input contents and output artifacts on
/// disk. Paths of prepared inputs may change when a complete output directory
/// is relocated, so input identity is bound by role, byte length, and SHA-256.
pub(crate) fn validate_recorded_completion(
    manifest_path: &Path,
    expected_gene: &str,
    expected_cluster_mode: &str,
    expected_inputs: &[InputArtifact],
    expected_outputs: &[OutputSpec],
) -> anyhow::Result<()> {
    let metadata = fs::symlink_metadata(manifest_path)
        .with_context(|| format!("inspect completion manifest {manifest_path:?}"))?;
    if metadata.file_type().is_symlink() || !metadata.is_file() {
        anyhow::bail!("completion manifest is not a regular file: {manifest_path:?}");
    }

    let manifest = read_run_manifest(manifest_path)?;
    if manifest.schema_version != MANIFEST_SCHEMA_VERSION {
        anyhow::bail!(
            "completion manifest schema changed: expected {}, found {}",
            MANIFEST_SCHEMA_VERSION,
            manifest.schema_version
        );
    }
    if manifest.status != "complete" {
        anyhow::bail!(
            "completion manifest status is {:?}, not complete",
            manifest.status
        );
    }
    let actual_fingerprint = manifest.request.fingerprint()?;
    if manifest.request_fingerprint != actual_fingerprint {
        anyhow::bail!("completion manifest request fingerprint is inconsistent");
    }
    if manifest.request.gene != expected_gene {
        anyhow::bail!(
            "completion manifest gene mismatch: expected {expected_gene:?}, found {:?}",
            manifest.request.gene
        );
    }
    if manifest.request.options.cluster_mode != expected_cluster_mode {
        anyhow::bail!(
            "completion manifest cluster mode mismatch: expected {expected_cluster_mode:?}, found {:?}",
            manifest.request.options.cluster_mode
        );
    }
    if manifest.request.tool != ToolIdentity::current() {
        anyhow::bail!("completion manifest was produced by a different tool build");
    }

    if manifest.request.inputs.len() != expected_inputs.len() {
        anyhow::bail!("completion manifest input set does not match the required input set");
    }
    for expected in expected_inputs {
        let matches = manifest
            .request
            .inputs
            .iter()
            .filter(|recorded| recorded.role == expected.role)
            .collect::<Vec<_>>();
        if matches.len() != 1 {
            anyhow::bail!(
                "completion manifest must contain exactly one {:?} input",
                expected.role
            );
        }
        let recorded = matches[0];
        if recorded.sha256 != expected.sha256 || recorded.bytes != expected.bytes {
            anyhow::bail!(
                "completion manifest input content changed for role {:?}",
                expected.role
            );
        }
    }

    if manifest.outputs.len() != expected_outputs.len() {
        anyhow::bail!("completion manifest output set does not match the required output set");
    }
    for spec in expected_outputs {
        let stored_path = spec.stored_path()?;
        let matches = manifest
            .outputs
            .iter()
            .filter(|output| output.role == spec.role && output.path == stored_path)
            .collect::<Vec<_>>();
        if matches.len() != 1 {
            anyhow::bail!(
                "completion manifest must contain exactly one {:?} output at {:?}",
                spec.role,
                stored_path
            );
        }
        let output_metadata = fs::symlink_metadata(&spec.path)
            .with_context(|| format!("inspect recorded output {:?}", spec.path))?;
        if output_metadata.file_type().is_symlink() || !output_metadata.is_file() {
            anyhow::bail!("recorded output is not a regular file: {:?}", spec.path);
        }
        let actual = OutputArtifact::from_spec(spec)?;
        if *matches[0] != actual {
            anyhow::bail!(
                "recorded output content or record count changed for role {:?}",
                spec.role
            );
        }
    }
    Ok(())
}

fn assess_cache_inner(
    manifest_path: &Path,
    expected_request: &RunRequest,
    expected_outputs: &[OutputSpec],
) -> Result<(), String> {
    let manifest_metadata = match fs::symlink_metadata(manifest_path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            return Err("manifest_missing".to_owned());
        }
        Err(_) => return Err("manifest_unreadable".to_owned()),
    };
    if manifest_metadata.file_type().is_symlink() || !manifest_metadata.is_file() {
        return Err("manifest_not_regular".to_owned());
    }
    let bytes = fs::read(manifest_path).map_err(|_| "manifest_unreadable".to_owned())?;
    let manifest: RunManifest =
        serde_json::from_slice(&bytes).map_err(|_| "manifest_invalid_json".to_owned())?;
    if manifest.schema_version != MANIFEST_SCHEMA_VERSION {
        return Err("manifest_schema_changed".to_owned());
    }
    if manifest.status != "complete" {
        return Err("manifest_incomplete".to_owned());
    }

    let expected_fingerprint = expected_request
        .fingerprint()
        .map_err(|_| "request_fingerprint_failed".to_owned())?;
    if manifest.request_fingerprint != expected_fingerprint || manifest.request != *expected_request
    {
        return Err(request_change_reason(&manifest.request, expected_request).to_owned());
    }

    if manifest.outputs.len() != expected_outputs.len() {
        return Err("output_set_changed".to_owned());
    }
    for spec in expected_outputs {
        let stored_path = spec
            .stored_path()
            .map_err(|_| "output_path_invalid".to_owned())?;
        let Some(recorded) = manifest
            .outputs
            .iter()
            .find(|output| output.role == spec.role && output.path == stored_path)
        else {
            return Err("output_set_changed".to_owned());
        };
        let output_metadata = match fs::symlink_metadata(&spec.path) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return Err(format!("output_missing:{}", spec.role));
            }
            Err(_) => return Err(format!("output_unreadable:{}", spec.role)),
        };
        if output_metadata.file_type().is_symlink() || !output_metadata.is_file() {
            return Err(format!("output_not_regular:{}", spec.role));
        }
        let actual_hash =
            sha256_file(&spec.path).map_err(|_| format!("output_unreadable:{}", spec.role))?;
        if recorded.sha256 != actual_hash {
            return Err(format!("output_hash_mismatch:{}", spec.role));
        }
        let actual_bytes = fs::metadata(&spec.path)
            .map_err(|_| format!("output_unreadable:{}", spec.role))?
            .len();
        if recorded.bytes != actual_bytes {
            return Err(format!("output_size_mismatch:{}", spec.role));
        }
        let actual_records = count_records(&spec.path, spec.record_count_kind)
            .map_err(|_| format!("output_unreadable:{}", spec.role))?;
        if recorded.records != actual_records {
            return Err(format!("output_record_count_mismatch:{}", spec.role));
        }
    }
    Ok(())
}

fn request_change_reason(previous: &RunRequest, current: &RunRequest) -> &'static str {
    if previous.gene != current.gene {
        return "gene_changed";
    }
    if previous.tool != current.tool {
        return "tool_version_changed";
    }
    if previous.seed != current.seed {
        return "seed_changed";
    }
    if previous.inputs != current.inputs {
        let reads_changed = input_changed(&previous.inputs, &current.inputs, "reads");
        let reference_changed = input_changed(&previous.inputs, &current.inputs, "reference");
        return match (reads_changed, reference_changed) {
            (true, false) => "reads_changed",
            (false, true) => "reference_changed",
            _ => "inputs_changed",
        };
    }
    if previous.options.cluster_mode != current.options.cluster_mode {
        return "cluster_mode_changed";
    }
    if previous.options != current.options {
        return "effective_options_changed";
    }
    "request_fingerprint_mismatch"
}

fn input_changed(previous: &[InputArtifact], current: &[InputArtifact], role: &str) -> bool {
    previous.iter().find(|input| input.role == role)
        != current.iter().find(|input| input.role == role)
}

pub(crate) fn write_completion_manifest(path: &Path, manifest: &RunManifest) -> anyhow::Result<()> {
    atomic_write_with(path, |file| {
        let mut writer = BufWriter::new(file);
        serde_json::to_writer_pretty(&mut writer, manifest)
            .with_context(|| format!("serialize manifest {path:?}"))?;
        writer.write_all(b"\n")?;
        writer.flush()?;
        Ok(())
    })
}

/// Remove the completion marker before publishing any replacement output.
pub(crate) fn invalidate_completion_manifest(path: &Path) -> anyhow::Result<()> {
    if path.exists() {
        fs::remove_file(path)
            .with_context(|| format!("invalidate completion manifest {path:?}"))?;
        sync_parent(path)?;
    }
    Ok(())
}

pub(crate) fn atomic_write_with<T, F>(path: &Path, write: F) -> anyhow::Result<T>
where
    F: FnOnce(&mut fs::File) -> anyhow::Result<T>,
{
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent).with_context(|| format!("create output directory {parent:?}"))?;
    let (temporary, reserved_file) = reserve_temporary_path(path)?;
    let guard = TemporaryGuard(temporary.clone());
    // Declare the open descriptor after the cleanup guard so every early return
    // closes it before the guard unlinks the temporary path (required on Windows).
    let mut reservation = reserved_file;
    let reserved_metadata = reservation
        .metadata()
        .with_context(|| format!("stat reserved temporary output {temporary:?}"))?;
    verify_regular_file(&reserved_metadata, &temporary)?;
    verify_single_link(&reserved_metadata, &temporary)?;
    let reserved_identity = snapshot_file_identity(&reservation, &temporary)?;

    // The create_new descriptor is deliberately the only write target exposed
    // to the callback. Reopening the predictable temporary pathname would let
    // a concurrent replacement redirect writes into another file.
    let value = write(&mut reservation)?;
    reservation
        .flush()
        .with_context(|| format!("flush temporary output {temporary:?}"))?;
    reservation
        .sync_all()
        .with_context(|| format!("sync temporary output {temporary:?}"))?;

    let descriptor_metadata = reservation
        .metadata()
        .with_context(|| format!("restat temporary output descriptor {temporary:?}"))?;
    verify_regular_file(&descriptor_metadata, &temporary)?;
    verify_single_link(&descriptor_metadata, &temporary)?;
    let linked_metadata = fs::symlink_metadata(&temporary)
        .with_context(|| format!("restat linked temporary output before publish {temporary:?}"))?;
    verify_regular_file(&linked_metadata, &temporary)?;
    verify_single_link(&linked_metadata, &temporary)?;
    verify_path_identity(&reserved_identity, &temporary)?;

    fs::rename(&temporary, path)
        .with_context(|| format!("publish temporary output {temporary:?} -> {path:?}"))?;
    let published_metadata =
        fs::symlink_metadata(path).with_context(|| format!("stat published output {path:?}"))?;
    verify_regular_file(&published_metadata, path)?;
    verify_single_link(&published_metadata, path)?;
    verify_path_identity(&reserved_identity, path)?;
    sync_parent(path)?;
    drop(reservation);
    std::mem::forget(guard);
    Ok(value)
}

pub(crate) fn atomic_copy(from: &Path, to: &Path) -> anyhow::Result<u64> {
    atomic_write_with(to, |output| {
        let source_metadata =
            fs::symlink_metadata(from).with_context(|| format!("inspect copy source {from:?}"))?;
        if source_metadata.file_type().is_symlink() || !source_metadata.is_file() {
            anyhow::bail!("copy source is not a regular file: {from:?}");
        }
        let mut input =
            fs::File::open(from).with_context(|| format!("open copy source {from:?}"))?;
        let source_identity = snapshot_file_identity(&input, from)?;
        verify_path_identity(&source_identity, from)?;
        let copied = std::io::copy(&mut input, output)
            .with_context(|| format!("copy {from:?} into atomic output {to:?}"))?;
        output
            .set_permissions(source_metadata.permissions())
            .with_context(|| format!("copy permissions from {from:?} to {to:?}"))?;
        Ok(copied)
    })
}

fn reserve_temporary_path(path: &Path) -> anyhow::Result<(PathBuf, fs::File)> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .with_context(|| format!("atomic output has no UTF-8 file name: {path:?}"))?;
    for _ in 0..1000 {
        let sequence = TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let candidate = parent.join(format!(
            ".{file_name}.tmp.{}.{}",
            std::process::id(),
            sequence
        ));
        match fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => return Ok((candidate, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("reserve temporary output {candidate:?}"));
            }
        }
    }
    anyhow::bail!("could not allocate temporary path next to {path:?}")
}

fn verify_regular_file(metadata: &fs::Metadata, path: &Path) -> anyhow::Result<()> {
    if !metadata.is_file() {
        anyhow::bail!("atomic output is not a regular file: {path:?}");
    }
    Ok(())
}

struct FileIdentity {
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
    #[cfg(not(unix))]
    handle: Handle,
}

#[cfg(unix)]
fn snapshot_file_identity(file: &fs::File, path: &Path) -> anyhow::Result<FileIdentity> {
    use std::os::unix::fs::MetadataExt;

    let metadata = file
        .metadata()
        .with_context(|| format!("snapshot atomic output identity {path:?}"))?;
    Ok(FileIdentity {
        device: metadata.dev(),
        inode: metadata.ino(),
    })
}

#[cfg(not(unix))]
fn snapshot_file_identity(file: &fs::File, path: &Path) -> anyhow::Result<FileIdentity> {
    let duplicate = file
        .try_clone()
        .with_context(|| format!("duplicate atomic output descriptor {path:?}"))?;
    let handle = Handle::from_file(duplicate)
        .with_context(|| format!("snapshot atomic output identity {path:?}"))?;
    Ok(FileIdentity { handle })
}

fn verify_path_identity(expected: &FileIdentity, path: &Path) -> anyhow::Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;

        let actual = fs::symlink_metadata(path)
            .with_context(|| format!("stat linked atomic output identity {path:?}"))?;
        if expected.device != actual.dev() || expected.inode != actual.ino() {
            anyhow::bail!("temporary output changed identity while publishing: {path:?}");
        }
    }
    #[cfg(not(unix))]
    {
        let actual = Handle::from_path(path)
            .with_context(|| format!("open linked atomic output for identity check {path:?}"))?;
        if expected.handle != actual {
            anyhow::bail!("temporary output changed identity while publishing: {path:?}");
        }
    }
    Ok(())
}

#[cfg(unix)]
fn verify_single_link(metadata: &fs::Metadata, path: &Path) -> anyhow::Result<()> {
    use std::os::unix::fs::MetadataExt;

    if metadata.nlink() != 1 {
        anyhow::bail!(
            "atomic output must retain exactly one directory link while publishing: {path:?}"
        );
    }
    Ok(())
}

#[cfg(not(unix))]
fn verify_single_link(_metadata: &fs::Metadata, _path: &Path) -> anyhow::Result<()> {
    Ok(())
}

fn sync_parent(path: &Path) -> anyhow::Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let directory =
        fs::File::open(parent).with_context(|| format!("open output dir {parent:?}"))?;
    directory
        .sync_all()
        .with_context(|| format!("sync output dir {parent:?}"))
}

struct TemporaryGuard(PathBuf);

impl Drop for TemporaryGuard {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

pub(crate) fn sha256_file(path: &Path) -> anyhow::Result<String> {
    let file = fs::File::open(path).with_context(|| format!("open file for SHA-256 {path:?}"))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 1024 * 1024];
    loop {
        let read = reader
            .read(&mut buffer)
            .with_context(|| format!("hash file {path:?}"))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("sha256:{:x}", hasher.finalize()))
}

fn sha256_bytes(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("sha256:{:x}", hasher.finalize())
}

fn count_records(path: &Path, kind: RecordCountKind) -> anyhow::Result<u64> {
    let file =
        fs::File::open(path).with_context(|| format!("open output for record count {path:?}"))?;
    let reader = BufReader::new(file);
    let nonempty = reader.lines().try_fold(0u64, |count, line| {
        let line = line?;
        Ok::<_, std::io::Error>(count + u64::from(!line.trim().is_empty()))
    })?;
    Ok(match kind {
        RecordCountKind::NonEmptyLines => nonempty,
        RecordCountKind::HeaderThenNonEmptyLines => nonempty.saturating_sub(1),
    })
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn temp_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-artifact-manifest-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn request(root: &Path) -> RunRequest {
        let reads = root.join("reads.bed");
        let reference = root.join("reference.bed");
        fs::write(&reads, "read\n").unwrap();
        fs::write(&reference, "reference\n").unwrap();
        RunRequest {
            gene: "GENEA".to_owned(),
            inputs: vec![
                InputArtifact::from_file("reads", &reads).unwrap(),
                InputArtifact::from_file("reference", &reference).unwrap(),
            ],
            options: EffectiveOptions {
                cluster_mode: "clusterj".to_owned(),
                prepare_fraction_read: 0.01,
                prepare_fraction_ref: 0.05,
                sw_score: -1,
                batch_size: 500,
                batch_rounds: 100,
                name2_mode: "coverage".to_owned(),
                platform_preset: "generic".to_owned(),
                junction_correction_offset: 10,
                junction_correction_min_support: 1,
                sl_partial_five_prime_offset: 10,
                sl_same_junction_five_prime_offset: 10,
                sl_five_prime_cluster_offset: 10,
                sl_min_five_prime_cluster_support: 1,
                same_junction_three_prime_offset: 10,
                three_prime_cluster_offset: 10,
                min_three_prime_cluster_support: 1,
                overlap_cutoff1: 0.05,
                overlap_cutoff2: 0.01,
                overlap_intron_weight: 0.5,
                assignment_mode: "unique".to_owned(),
                unique_assignment_junction_offset: 15,
                invalid_read_policy: "skip".to_owned(),
                downsample_selected: true,
                max_reads_per_gene: 50_000,
            },
            tool: ToolIdentity::current(),
            seed: 7,
        }
    }

    fn published(root: &Path, request: RunRequest) -> (PathBuf, Vec<OutputSpec>) {
        let output = root.join("isoforms.bed");
        atomic_write_with(&output, |temporary| {
            temporary.write_all(b"isoform-1\nisoform-2\n")?;
            Ok(())
        })
        .unwrap();
        let specs = vec![OutputSpec::new(
            "isoforms",
            output,
            RecordCountKind::NonEmptyLines,
        )];
        let manifest_path = root.join(MANIFEST_FILE_NAME);
        let manifest = RunManifest::complete(request, &specs).unwrap();
        write_completion_manifest(&manifest_path, &manifest).unwrap();
        (manifest_path, specs)
    }

    #[test]
    fn exact_request_and_verified_output_are_reused() {
        let root = temp_dir("reuse");
        let request = request(&root);
        let (manifest, specs) = published(&root, request.clone());
        assert_eq!(
            assess_cache(&manifest, &request, &specs),
            CacheDecision::Reuse
        );
    }

    #[test]
    fn version_two_manifest_without_a_source_fingerprint_is_explicitly_stale() {
        let root = temp_dir("v2-tool-identity");
        let request = request(&root);
        let (manifest, specs) = published(&root, request.clone());
        let mut legacy: serde_json::Value =
            serde_json::from_slice(&fs::read(&manifest).unwrap()).unwrap();
        legacy["schema_version"] = serde_json::json!(2);
        legacy["request"]["tool"]
            .as_object_mut()
            .unwrap()
            .remove("source_fingerprint");
        fs::write(&manifest, serde_json::to_vec_pretty(&legacy).unwrap()).unwrap();

        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "manifest_schema_changed"
        );
    }

    #[cfg(unix)]
    #[test]
    fn normal_cache_rejects_symlink_manifests_and_outputs() {
        use std::os::unix::fs::symlink;

        let root = temp_dir("symlinks");
        let request = request(&root);
        let (manifest, specs) = published(&root, request.clone());

        let manifest_target = root.join("real-run.json");
        fs::rename(&manifest, &manifest_target).unwrap();
        symlink(&manifest_target, &manifest).unwrap();
        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "manifest_not_regular"
        );

        fs::remove_file(&manifest).unwrap();
        fs::rename(&manifest_target, &manifest).unwrap();
        let output_target = root.join("real-isoforms.bed");
        fs::rename(&specs[0].path, &output_target).unwrap();
        symlink(&output_target, &specs[0].path).unwrap();
        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "output_not_regular:isoforms"
        );
    }

    #[test]
    fn recorded_completion_revalidates_inputs_and_outputs_for_count_only_use() {
        let root = temp_dir("recorded-completion");
        let request = request(&root);
        let expected_inputs = request.inputs.clone();
        let (manifest, specs) = published(&root, request);

        validate_recorded_completion(&manifest, "GENEA", "clusterj", &expected_inputs, &specs)
            .expect("valid recorded completion");

        fs::write(&specs[0].path, "tampered\n").unwrap();
        let error =
            validate_recorded_completion(&manifest, "GENEA", "clusterj", &expected_inputs, &specs)
                .expect_err("tampered output must not be trusted");
        assert!(
            error
                .to_string()
                .contains("recorded output content or record count changed"),
            "{error:#}"
        );

        fs::write(root.join("reads.bed"), "changed input\n").unwrap();
        let changed_inputs = vec![
            InputArtifact::from_file("reads", &root.join("reads.bed")).unwrap(),
            InputArtifact::from_file("reference", &root.join("reference.bed")).unwrap(),
        ];
        let error =
            validate_recorded_completion(&manifest, "GENEA", "clusterj", &changed_inputs, &specs)
                .expect_err("changed prepared input must not be trusted");
        assert!(
            error
                .to_string()
                .contains("input content changed for role \"reads\""),
            "{error:#}"
        );
    }

    #[test]
    fn every_required_invalidation_class_is_reported() {
        let root = temp_dir("invalidation");
        let request = request(&root);
        let (manifest, specs) = published(&root, request.clone());

        let mut changed = request.clone();
        changed.options.cluster_mode = "cluster".to_owned();
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "cluster_mode_changed"
        );

        let mut changed = request.clone();
        changed.options.overlap_cutoff1 = 0.25;
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "effective_options_changed"
        );

        let mut changed = request.clone();
        changed.options.unique_assignment_junction_offset += 1;
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "effective_options_changed"
        );

        let mut changed = request.clone();
        changed.seed += 1;
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "seed_changed"
        );

        let mut changed = request.clone();
        changed.tool.package_version = "999.0.0".to_owned();
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "tool_version_changed"
        );

        let mut changed = request.clone();
        changed.tool.source_fingerprint = "sha256:dirty-source".to_owned();
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "tool_version_changed"
        );

        fs::write(root.join("reads.bed"), "changed reads\n").unwrap();
        let mut changed = request.clone();
        changed.inputs[0] = InputArtifact::from_file("reads", &root.join("reads.bed")).unwrap();
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "reads_changed"
        );

        fs::write(root.join("reference.bed"), "changed reference\n").unwrap();
        let mut changed = request;
        changed.inputs[1] =
            InputArtifact::from_file("reference", &root.join("reference.bed")).unwrap();
        assert_eq!(
            assess_cache(&manifest, &changed, &specs).reason(),
            "reference_changed"
        );
    }

    #[test]
    fn missing_legacy_or_corrupt_outputs_are_never_reused() {
        let root = temp_dir("legacy-corrupt");
        let request = request(&root);
        let output = root.join("isoforms.bed");
        fs::write(&output, "legacy\n").unwrap();
        let specs = vec![OutputSpec::new(
            "isoforms",
            &output,
            RecordCountKind::NonEmptyLines,
        )];
        let manifest = root.join(MANIFEST_FILE_NAME);
        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "manifest_missing"
        );

        let (manifest, specs) = published(&root, request.clone());
        fs::write(&output, "corrupt\n").unwrap();
        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "output_hash_mismatch:isoforms"
        );
    }

    #[test]
    fn failed_or_interrupted_publish_cannot_replace_a_valid_result() {
        let root = temp_dir("interrupted");
        let output = root.join("result.txt");
        fs::write(&output, "old complete output\n").unwrap();

        let error = atomic_write_with::<(), _>(&output, |temporary| {
            temporary.write_all(b"partial")?;
            anyhow::bail!("simulated termination")
        })
        .unwrap_err();
        assert!(error.to_string().contains("simulated termination"));
        assert_eq!(
            fs::read_to_string(&output).unwrap(),
            "old complete output\n"
        );
        assert!(fs::read_dir(&root).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains(".tmp.")));

        let request = request(&root);
        let (manifest, specs) = published(&root, request.clone());
        invalidate_completion_manifest(&manifest).unwrap();
        fs::write(&specs[0].path, "partially replaced\n").unwrap();
        assert_eq!(
            assess_cache(&manifest, &request, &specs).reason(),
            "manifest_missing"
        );
    }

    #[cfg(unix)]
    #[test]
    fn replacing_the_reserved_temp_path_cannot_clobber_a_valid_output() {
        use std::os::unix::fs::symlink;

        let root = temp_dir("temp-path-swap");
        let output = root.join("result.txt");
        let canary = root.join("canary.txt");
        fs::write(&output, "old complete output\n").unwrap();
        fs::write(&canary, "must remain unchanged\n").unwrap();

        let error = atomic_write_with(&output, |temporary| {
            temporary.write_all(b"intended replacement\n")?;
            let reserved_path = fs::read_dir(&root)?
                .map(|entry| entry.map(|entry| entry.path()))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .find(|path| {
                    path.file_name()
                        .is_some_and(|name| name.to_string_lossy().contains(".result.txt.tmp."))
                })
                .context("find reserved atomic temp path")?;
            fs::remove_file(&reserved_path)?;
            symlink(&canary, &reserved_path)?;
            temporary.write_all(b"descriptor-only follow-up\n")?;
            Ok(())
        })
        .unwrap_err();

        assert!(
            error.to_string().contains("exactly one directory link")
                || error.to_string().contains("changed identity"),
            "unexpected error: {error:#}"
        );
        assert_eq!(
            fs::read_to_string(&output).unwrap(),
            "old complete output\n"
        );
        assert_eq!(
            fs::read_to_string(&canary).unwrap(),
            "must remain unchanged\n"
        );
        assert!(fs::read_dir(&root).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains(".tmp.")));
    }
}
