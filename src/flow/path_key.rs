use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Component, Path, PathBuf};

use anyhow::Context;
use sha2::{Digest, Sha256};

use crate::flow::artifact_manifest::atomic_write_with;

/// Gene identifiers are accepted up to this many UTF-8 bytes. Longer values are rejected before
/// any output path is constructed. This limit bounds metadata and diagnostic memory while still
/// allowing substantially longer identifiers than a filesystem component can hold.
pub(crate) const MAX_GENE_ID_BYTES: usize = 4_096;

const MAX_PREFIX_BYTES: usize = 128;
const MAX_REVERSIBLE_PATH_KEY_BYTES: usize = 180;
const ENCODED_KEY_PREFIX: &str = "g~";
const HASHED_KEY_PREFIX: &str = "g~h~";
const PATH_MAP_VERSION: u32 = 1;

pub(crate) const GENE_ID_MARKER_FILE: &str = ".trackcluster_gene_id";

#[derive(Debug)]
struct ExternalInput {
    label: String,
    original: PathBuf,
    resolved: PathBuf,
    #[cfg(unix)]
    identity: (u64, u64),
}

fn lexical_absolute(path: &Path) -> anyhow::Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .context("resolve current working directory")?
            .join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            Component::Prefix(_) | Component::RootDir => {
                normalized.push(component.as_os_str());
            }
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            Component::Normal(value) => normalized.push(value),
        }
    }
    Ok(normalized)
}

fn resolve_allow_missing(path: &Path) -> anyhow::Result<PathBuf> {
    let absolute = lexical_absolute(path)?;
    let mut ancestor = absolute.as_path();
    loop {
        match fs::canonicalize(ancestor) {
            Ok(resolved_ancestor) => {
                let suffix = absolute
                    .strip_prefix(ancestor)
                    .expect("ancestor is derived from the absolute path");
                return Ok(resolved_ancestor.join(suffix));
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                ancestor = ancestor
                    .parent()
                    .with_context(|| format!("resolve nearest existing ancestor of {path:?}"))?;
            }
            Err(error) => {
                return Err(error).with_context(|| format!("resolve path {path:?}"));
            }
        }
    }
}

fn inspect_output_entry_for_alias(entry: &Path, inputs: &[ExternalInput]) -> anyhow::Result<()> {
    let link_metadata = fs::symlink_metadata(entry)
        .with_context(|| format!("inspect output-root entry {entry:?}"))?;
    if link_metadata.file_type().is_symlink() {
        anyhow::bail!(
            "pipeline-owned output_root contains a symlink entry {entry:?}; remove symlinks before running the pipeline"
        );
    }
    #[cfg(unix)]
    if link_metadata.is_file() {
        use std::os::unix::fs::MetadataExt as _;
        if link_metadata.nlink() > 1 {
            anyhow::bail!(
                "pipeline-owned output_root contains a multiply-linked file {entry:?}; remove hard links before running the pipeline"
            );
        }
    }

    let resolved = match fs::canonicalize(entry) {
        Ok(path) => Some(path),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => {
            return Err(error).with_context(|| format!("resolve output-root entry {entry:?}"));
        }
    };
    let metadata = match fs::metadata(entry) {
        Ok(metadata) => Some(metadata),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => {
            return Err(error).with_context(|| format!("inspect output-root entry {entry:?}"));
        }
    };

    for input in inputs {
        let aliases_target = resolved
            .as_ref()
            .is_some_and(|path| path == &input.resolved);
        #[cfg(unix)]
        let aliases_inode = metadata.as_ref().is_some_and(|metadata| {
            use std::os::unix::fs::MetadataExt as _;
            (metadata.dev(), metadata.ino()) == input.identity
        });
        #[cfg(not(unix))]
        let aliases_inode = false;

        if aliases_target || aliases_inode {
            anyhow::bail!(
                "external {} {:?} aliases existing entry {entry:?} in pipeline-owned output_root; move the input outside output_root",
                input.label,
                input.original
            );
        }
    }
    Ok(())
}

fn scan_output_tree_for_input_aliases(
    directory: &Path,
    canonical_root: &Path,
    inputs: &[ExternalInput],
    visited_directories: &mut HashSet<PathBuf>,
) -> anyhow::Result<()> {
    for entry in
        fs::read_dir(directory).with_context(|| format!("read output_root {directory:?}"))?
    {
        let entry = entry.with_context(|| format!("read entry in output_root {directory:?}"))?;
        let path = entry.path();
        inspect_output_entry_for_alias(&path, inputs)?;

        let metadata = match fs::metadata(&path) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(error) => {
                return Err(error).with_context(|| format!("inspect output-root entry {path:?}"));
            }
        };
        if !metadata.is_dir() {
            continue;
        }
        let resolved = fs::canonicalize(&path)
            .with_context(|| format!("resolve output-root directory {path:?}"))?;
        if resolved.starts_with(canonical_root) && visited_directories.insert(resolved.clone()) {
            scan_output_tree_for_input_aliases(
                &resolved,
                canonical_root,
                inputs,
                visited_directories,
            )?;
        }
    }
    Ok(())
}

/// Treat `output_root` as a pipeline-owned tree and reject external source files that could be
/// overwritten by any current or future artifact path beneath it.
///
/// Inputs are rejected when their resolved target is beneath the managed tree, or when an
/// existing entry anywhere in that tree resolves to or hard-links the same file. The latter check
/// prevents an apparently external source from being clobbered through a pre-existing alias.
pub(crate) fn reject_external_inputs_in_output_root<'a>(
    output_root: &Path,
    inputs: impl IntoIterator<Item = (&'a str, &'a Path)>,
) -> anyhow::Result<()> {
    let resolved_root = resolve_allow_missing(output_root)
        .with_context(|| format!("resolve pipeline-owned output_root {output_root:?}"))?;
    let mut protected = Vec::new();
    for (label, path) in inputs {
        let resolved =
            fs::canonicalize(path).with_context(|| format!("resolve external {label} {path:?}"))?;
        if resolved.starts_with(&resolved_root) {
            anyhow::bail!(
                "external {label} {path:?} resolves beneath pipeline-owned output_root {output_root:?}; move the input outside output_root"
            );
        }
        let metadata =
            fs::metadata(path).with_context(|| format!("inspect external {label} {path:?}"))?;
        #[cfg(not(unix))]
        let _ = &metadata;
        #[cfg(unix)]
        let identity = {
            use std::os::unix::fs::MetadataExt as _;
            (metadata.dev(), metadata.ino())
        };
        protected.push(ExternalInput {
            label: label.to_owned(),
            original: path.to_path_buf(),
            resolved,
            #[cfg(unix)]
            identity,
        });
    }

    match fs::metadata(output_root) {
        Ok(metadata) if metadata.is_dir() => {
            let canonical_root = fs::canonicalize(output_root)
                .with_context(|| format!("resolve existing output_root {output_root:?}"))?;
            let mut visited_directories = HashSet::from([canonical_root.clone()]);
            scan_output_tree_for_input_aliases(
                &canonical_root,
                &canonical_root,
                &protected,
                &mut visited_directories,
            )
        }
        Ok(_) => anyhow::bail!("pipeline-owned output_root is not a directory: {output_root:?}"),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("inspect output_root {output_root:?}")),
    }
}

/// A biological gene identifier. It is never used directly as a path component without first
/// deriving a [`GenePathKey`].
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct GeneId(String);

impl GeneId {
    pub(crate) fn parse(value: &str) -> anyhow::Result<Self> {
        validate_single_component("gene id", value, MAX_GENE_ID_BYTES)?;
        Ok(Self(value.to_owned()))
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }

    pub(crate) fn path_key(&self) -> GenePathKey {
        GenePathKey::from_gene_id(self)
    }
}

/// A deterministic, bounded filesystem representation of a biological [`GeneId`]. Common safe
/// ASCII IDs retain their historical spelling. Other IDs are percent encoded when that fits, or
/// represented by a stable SHA-256 hash when it does not. The mapping file is authoritative for
/// hashed keys and collisions are rejected before output is written.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct GenePathKey(String);

impl GenePathKey {
    fn from_gene_id(gene: &GeneId) -> Self {
        let bytes = gene.as_str().as_bytes();
        if bytes.len() <= MAX_REVERSIBLE_PATH_KEY_BYTES
            && bytes.iter().copied().all(is_unencoded_path_byte)
        {
            return Self(gene.as_str().to_owned());
        }

        let mut encoded = String::with_capacity(ENCODED_KEY_PREFIX.len() + bytes.len() * 3);
        encoded.push_str(ENCODED_KEY_PREFIX);
        for &byte in bytes {
            if is_unencoded_path_byte(byte) {
                encoded.push(char::from(byte));
            } else {
                encoded.push('%');
                encoded.push(hex_digit(byte >> 4));
                encoded.push(hex_digit(byte & 0x0f));
            }
        }
        if encoded.len() <= MAX_REVERSIBLE_PATH_KEY_BYTES {
            return Self(encoded);
        }

        let mut hasher = Sha256::new();
        hasher.update(b"trackcluster-gene-path-v1\0");
        hasher.update(bytes);
        Self(format!("{HASHED_KEY_PREFIX}{:x}", hasher.finalize()))
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }

    /// Decode keys that retain the biological ID. Hashed keys require the persisted mapping.
    #[cfg(test)]
    pub(crate) fn decode_reversible(value: &str) -> anyhow::Result<Option<GeneId>> {
        if value.starts_with(HASHED_KEY_PREFIX) {
            return Ok(None);
        }
        if !value.starts_with(ENCODED_KEY_PREFIX) {
            return GeneId::parse(value).map(Some);
        }

        let encoded = &value[ENCODED_KEY_PREFIX.len()..];
        let bytes = encoded.as_bytes();
        let mut decoded = Vec::with_capacity(bytes.len());
        let mut index = 0;
        while index < bytes.len() {
            if bytes[index] != b'%' {
                decoded.push(bytes[index]);
                index += 1;
                continue;
            }
            if index + 2 >= bytes.len() {
                anyhow::bail!("invalid truncated percent escape in gene path key {value:?}");
            }
            let high = parse_hex_digit(bytes[index + 1])
                .with_context(|| format!("invalid percent escape in gene path key {value:?}"))?;
            let low = parse_hex_digit(bytes[index + 2])
                .with_context(|| format!("invalid percent escape in gene path key {value:?}"))?;
            decoded.push((high << 4) | low);
            index += 3;
        }
        let decoded = String::from_utf8(decoded)
            .with_context(|| format!("gene path key {value:?} does not decode to UTF-8"))?;
        let gene = GeneId::parse(&decoded)?;
        if gene.path_key().as_str() != value {
            anyhow::bail!("non-canonical gene path key {value:?}");
        }
        Ok(Some(gene))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SafePathComponent(String);

impl SafePathComponent {
    pub(crate) fn parse(kind: &str, value: &str) -> anyhow::Result<Self> {
        validate_single_component(kind, value, MAX_PREFIX_BYTES)?;
        Ok(Self(value.to_owned()))
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

fn validate_single_component(kind: &str, value: &str, max_bytes: usize) -> anyhow::Result<()> {
    if value.is_empty() {
        anyhow::bail!("{kind} must not be empty");
    }
    if value.len() > max_bytes {
        anyhow::bail!(
            "{kind} is too long: {} UTF-8 bytes (maximum {max_bytes})",
            value.len()
        );
    }
    if value.chars().any(char::is_control) {
        anyhow::bail!("{kind} {value:?} contains a NUL or control character");
    }
    if value.contains('/') || value.contains('\\') {
        anyhow::bail!("{kind} {value:?} must not contain path separators");
    }
    if value == "." || value == ".." {
        anyhow::bail!("{kind} {value:?} is not allowed");
    }

    let path = Path::new(value);
    if path.is_absolute() {
        anyhow::bail!("{kind} {value:?} must not be an absolute path");
    }
    let mut components = path.components();
    if !matches!(components.next(), Some(Component::Normal(_))) || components.next().is_some() {
        anyhow::bail!("{kind} {value:?} must be one normal path component");
    }
    Ok(())
}

fn is_unencoded_path_byte(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.')
}

fn hex_digit(value: u8) -> char {
    char::from(if value < 10 {
        b'0' + value
    } else {
        b'A' + value - 10
    })
}

#[cfg(test)]
fn parse_hex_digit(value: u8) -> anyhow::Result<u8> {
    match value {
        b'0'..=b'9' => Ok(value - b'0'),
        b'a'..=b'f' => Ok(value - b'a' + 10),
        b'A'..=b'F' => Ok(value - b'A' + 10),
        _ => anyhow::bail!("expected hexadecimal digit"),
    }
}

pub(crate) fn validate_gene_ids<'a>(
    values: impl IntoIterator<Item = &'a str>,
) -> anyhow::Result<Vec<GeneId>> {
    let mut genes = Vec::new();
    for value in values {
        genes.push(GeneId::parse(value)?);
    }
    validate_unique_path_keys(&genes)?;
    Ok(genes)
}

fn validate_unique_path_keys(genes: &[GeneId]) -> anyhow::Result<()> {
    let mut by_key: HashMap<GenePathKey, &GeneId> = HashMap::new();
    for gene in genes {
        let key = gene.path_key();
        if let Some(previous) = by_key.insert(key.clone(), gene) {
            if previous != gene {
                anyhow::bail!(
                    "gene path-key collision: biological IDs {:?} and {:?} both map to {:?}",
                    previous.as_str(),
                    gene.as_str(),
                    key.as_str()
                );
            }
        }
    }
    Ok(())
}

/// Persist the biological-to-filesystem mapping as versioned run metadata.
pub(crate) fn write_gene_path_map(path: &Path, genes: &[GeneId]) -> anyhow::Result<()> {
    validate_unique_path_keys(genes)?;
    let mut expected = genes.to_vec();
    expected.sort();
    expected.dedup();

    atomic_write_with(path, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(
            writer,
            "# trackcluster_gene_path_map_version={PATH_MAP_VERSION}"
        )?;
        writeln!(writer, "gene_id\tpath_key")?;
        for gene in &expected {
            writeln!(writer, "{}\t{}", gene.as_str(), gene.path_key().as_str())?;
        }
        writer
            .flush()
            .with_context(|| format!("flush temporary gene path mapping {path:?}"))
    })?;

    let persisted = read_gene_path_map(path)?;
    if persisted != expected {
        anyhow::bail!("gene path mapping {path:?} did not round-trip exactly");
    }
    Ok(())
}

pub(crate) fn read_gene_path_map(path: &Path) -> anyhow::Result<Vec<GeneId>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("open gene path mapping {path:?}"))?,
    );
    let mut saw_version = false;
    let mut saw_header = false;
    let mut genes = Vec::new();
    for (zero_based, line) in reader.lines().enumerate() {
        let line_no = zero_based + 1;
        let line = line.with_context(|| format!("read gene path mapping {path:?}:{line_no}"))?;
        if line == format!("# trackcluster_gene_path_map_version={PATH_MAP_VERSION}") {
            saw_version = true;
            continue;
        }
        if line == "gene_id\tpath_key" {
            saw_header = true;
            continue;
        }
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if !saw_header {
            anyhow::bail!("gene path mapping {path:?}:{line_no} has data before its header");
        }
        let Some((gene_value, key_value)) = line.split_once('\t') else {
            anyhow::bail!("gene path mapping {path:?}:{line_no} must have two TSV fields");
        };
        let gene = GeneId::parse(gene_value)
            .with_context(|| format!("invalid gene in path mapping {path:?}:{line_no}"))?;
        let expected = gene.path_key();
        if expected.as_str() != key_value {
            anyhow::bail!(
                "gene path mapping {path:?}:{line_no} has key {key_value:?}, expected {:?}",
                expected.as_str()
            );
        }
        genes.push(gene);
    }
    if !saw_version {
        anyhow::bail!("gene path mapping {path:?} has no supported version marker");
    }
    if !saw_header {
        anyhow::bail!("gene path mapping {path:?} has no gene_id/path_key header");
    }
    let mut unique_genes = genes.clone();
    unique_genes.sort();
    unique_genes.dedup();
    if unique_genes.len() != genes.len() {
        anyhow::bail!("gene path mapping {path:?} contains a duplicate gene row");
    }
    validate_unique_path_keys(&genes)?;
    Ok(genes)
}

pub(crate) fn write_gene_id_marker(gene_dir: &Path, gene: &GeneId) -> anyhow::Result<()> {
    let marker = gene_dir.join(GENE_ID_MARKER_FILE);
    let root = gene_dir
        .parent()
        .context("gene directory has no output-root parent")?;
    ensure_destination_within(root, &marker)?;
    atomic_write_with(&marker, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(writer, "{}", gene.as_str())?;
        writer
            .flush()
            .with_context(|| format!("flush temporary gene identity marker {marker:?}"))
    })
}

/// Verify that a candidate (including its resolved existing ancestor) stays beneath `root`.
/// This catches lexical traversal and symlink redirection before a file is opened.
pub(crate) fn ensure_destination_within(root: &Path, candidate: &Path) -> anyhow::Result<()> {
    let relative = candidate.strip_prefix(root).with_context(|| {
        format!("destination {candidate:?} is not lexically beneath output root {root:?}")
    })?;
    if relative
        .components()
        .any(|component| !matches!(component, Component::Normal(_) | Component::CurDir))
    {
        anyhow::bail!(
            "destination {candidate:?} contains a non-normal component beneath root {root:?}"
        );
    }

    let canonical_root = fs::canonicalize(root).with_context(|| {
        format!("resolve output root before validating {candidate:?}: {root:?}")
    })?;
    let mut ancestor = candidate;
    loop {
        match fs::symlink_metadata(ancestor) {
            Ok(_) => break,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                ancestor = ancestor.parent().with_context(|| {
                    format!("destination {candidate:?} has no existing ancestor")
                })?;
            }
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("inspect destination ancestor {ancestor:?}"));
            }
        }
    }
    let resolved_ancestor = fs::canonicalize(ancestor)
        .with_context(|| format!("resolve destination ancestor {ancestor:?}"))?;
    if !resolved_ancestor.starts_with(&canonical_root) {
        anyhow::bail!(
            "resolved destination {candidate:?} escapes output root {root:?}: existing ancestor resolves to {resolved_ancestor:?}"
        );
    }
    Ok(())
}

pub(crate) fn gene_dir_path(root: &Path, gene: &GeneId) -> anyhow::Result<PathBuf> {
    let path = root.join(gene.path_key().as_str());
    ensure_destination_within(root, &path)?;
    Ok(path)
}

pub(crate) fn gene_artifact_path(
    root: &Path,
    gene: &GeneId,
    suffix: &str,
) -> anyhow::Result<PathBuf> {
    if suffix.is_empty()
        || suffix.contains('/')
        || suffix.contains('\\')
        || suffix.chars().any(char::is_control)
    {
        anyhow::bail!("invalid per-gene artifact suffix {suffix:?}");
    }
    let dir = gene_dir_path(root, gene)?;
    let path = dir.join(format!("{}{suffix}", gene.path_key().as_str()));
    ensure_destination_within(root, &path)?;
    Ok(path)
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
        let path = std::env::temp_dir().join(format!(
            "trackcluster_path_key_{prefix}_{}_{}",
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&path).expect("create temp dir");
        path
    }

    #[test]
    fn separates_biological_ids_from_bounded_path_keys() {
        let ascii = GeneId::parse("GENEA.1-alpha").unwrap();
        assert_eq!(ascii.path_key().as_str(), "GENEA.1-alpha");

        let unicode = GeneId::parse("基因-α.1").unwrap();
        let unicode_key = unicode.path_key();
        assert_ne!(unicode_key.as_str(), unicode.as_str());
        assert!(unicode_key.as_str().len() <= MAX_REVERSIBLE_PATH_KEY_BYTES);
        assert_eq!(
            GenePathKey::decode_reversible(unicode_key.as_str()).unwrap(),
            Some(unicode)
        );

        let long = GeneId::parse(&"长".repeat(400)).unwrap();
        let long_key = long.path_key();
        assert!(long_key.as_str().starts_with(HASHED_KEY_PREFIX));
        assert!(long_key.as_str().len() < 80);
        assert_eq!(
            GenePathKey::decode_reversible(long_key.as_str()).unwrap(),
            None
        );
    }

    #[test]
    fn rejects_traversal_absolute_separators_controls_and_oversize_ids() {
        for value in [
            "",
            ".",
            "..",
            "../escape",
            "/absolute",
            "a/b",
            "a\\b",
            "a\0b",
            "a\nb",
        ] {
            assert!(
                GeneId::parse(value).is_err(),
                "unexpectedly accepted {value:?}"
            );
        }
        let error = GeneId::parse(&"x".repeat(MAX_GENE_ID_BYTES + 1))
            .expect_err("oversized gene must fail")
            .to_string();
        assert!(error.contains("maximum 4096"), "{error}");
    }

    #[test]
    fn mapping_round_trips_unicode_and_hashed_ids() {
        let dir = fresh_temp_dir("mapping");
        let path = dir.join("sample_gene_paths.tsv");
        let genes = vec![
            GeneId::parse("GENEA").unwrap(),
            GeneId::parse(&"very-long-gene-".repeat(40)).unwrap(),
            GeneId::parse("基因-α.1").unwrap(),
        ];
        write_gene_path_map(&path, &genes).unwrap();
        assert_eq!(read_gene_path_map(&path).unwrap(), genes);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn pipeline_owned_tree_rejects_external_inputs_beneath_it() {
        let parent = fresh_temp_dir("owned_input");
        let root = parent.join("output");
        fs::create_dir_all(&root).unwrap();
        let reads = root.join("sample_novel.bed");
        fs::write(&reads, "source reads\n").unwrap();

        let error =
            reject_external_inputs_in_output_root(&root, [("reads input", reads.as_path())])
                .expect_err("an input beneath output_root must be rejected");
        assert!(
            format!("{error:#}").contains("resolves beneath pipeline-owned output_root"),
            "{error:#}"
        );
        assert_eq!(fs::read_to_string(&reads).unwrap(), "source reads\n");
        fs::remove_dir_all(parent).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn pipeline_owned_tree_rejects_symlinks_and_hard_links() {
        use std::os::unix::fs::symlink;

        let parent = fresh_temp_dir("owned_aliases");
        let root = parent.join("output");
        fs::create_dir_all(&root).unwrap();
        let input = parent.join("reads.bed");
        let target = parent.join("target.bed");
        fs::write(&input, "external reads\n").unwrap();
        fs::write(&target, "target\n").unwrap();

        let alias = root.join("generated.bed");
        symlink(&target, &alias).unwrap();
        let error =
            reject_external_inputs_in_output_root(&root, [("reads input", input.as_path())])
                .expect_err("a managed-tree symlink must be rejected");
        assert!(format!("{error:#}").contains("contains a symlink entry"));

        fs::remove_file(&alias).unwrap();
        fs::write(&alias, "generated\n").unwrap();
        fs::hard_link(&alias, parent.join("generated-hardlink.bed")).unwrap();
        let error =
            reject_external_inputs_in_output_root(&root, [("reads input", input.as_path())])
                .expect_err("a multiply-linked managed file must be rejected");
        assert!(format!("{error:#}").contains("contains a multiply-linked file"));

        assert_eq!(fs::read_to_string(&input).unwrap(), "external reads\n");
        fs::remove_dir_all(parent).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn rejects_a_destination_redirected_outside_root_by_symlink() {
        use std::os::unix::fs::symlink;

        let parent = fresh_temp_dir("symlink");
        let root = parent.join("root");
        let outside = parent.join("outside");
        fs::create_dir_all(&root).unwrap();
        fs::create_dir_all(&outside).unwrap();
        symlink(&outside, root.join("GENEA")).unwrap();

        let gene = GeneId::parse("GENEA").unwrap();
        let error = gene_artifact_path(&root, &gene, "_nano.bed")
            .expect_err("outside symlink must fail")
            .to_string();
        assert!(error.contains("escapes output root"), "{error}");
        fs::remove_dir_all(parent).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn identity_marker_cannot_overwrite_an_external_symlink_target() {
        use std::os::unix::fs::symlink;

        let parent = fresh_temp_dir("marker_symlink");
        let root = parent.join("root");
        let outside = parent.join("outside.txt");
        let gene_dir = root.join("GENEA");
        fs::create_dir_all(&gene_dir).unwrap();
        fs::write(&outside, "do-not-overwrite\n").unwrap();
        symlink(&outside, gene_dir.join(GENE_ID_MARKER_FILE)).unwrap();

        let gene = GeneId::parse("GENEA").unwrap();
        let error = write_gene_id_marker(&gene_dir, &gene)
            .expect_err("external marker symlink must fail")
            .to_string();
        assert!(error.contains("escapes output root"), "{error}");
        assert_eq!(fs::read_to_string(&outside).unwrap(), "do-not-overwrite\n");
        fs::remove_dir_all(parent).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn metadata_publication_replaces_internal_hard_links_without_clobbering_targets() {
        let parent = fresh_temp_dir("metadata_hard_links");
        let root = parent.join("root");
        let gene_dir = root.join("GENEA");
        fs::create_dir_all(&gene_dir).unwrap();
        let protected = gene_dir.join("GENEA_nano.bed");
        fs::write(&protected, "protected reads\n").unwrap();

        let marker = gene_dir.join(GENE_ID_MARKER_FILE);
        fs::hard_link(&protected, &marker).unwrap();
        let gene = GeneId::parse("GENEA").unwrap();
        write_gene_id_marker(&gene_dir, &gene).unwrap();
        assert_eq!(fs::read_to_string(&protected).unwrap(), "protected reads\n");
        assert_eq!(fs::read_to_string(&marker).unwrap(), "GENEA\n");

        let mapping = root.join("clusterj_batch_gene_paths.tsv");
        fs::hard_link(&protected, &mapping).unwrap();
        write_gene_path_map(&mapping, std::slice::from_ref(&gene)).unwrap();
        assert_eq!(fs::read_to_string(&protected).unwrap(), "protected reads\n");
        assert_eq!(read_gene_path_map(&mapping).unwrap(), vec![gene]);
        fs::remove_dir_all(parent).unwrap();
    }
}
