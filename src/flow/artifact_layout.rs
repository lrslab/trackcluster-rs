//! Validated paths for the merged artifacts produced by a full flow run.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::flow::path_key::{ensure_destination_within, SafePathComponent};
use crate::flow::path_key::{gene_artifact_path, GeneId};

const BATCH_TOP_LEVEL_NAMES: &[&str] = &[
    "clusterj_batch_gene_paths.tsv",
    "clusterj_batch_summary.txt",
    "clusterj_batch_errors.txt",
    "clusterj_batch_downsample.tsv",
    "cluster_batch_gene_paths.tsv",
    "cluster_batch_summary.txt",
    "cluster_batch_errors.txt",
    "cluster_batch_downsample.tsv",
];

const PREFIX_TOP_LEVEL_SUFFIXES: &[&str] = &[
    "_dedup.bed",
    "_novel.bed",
    "_gene.txt",
    "_gene_paths.tsv",
    "_rejected_reads.tsv",
    "_pooled_reads.bed",
    "_isoform.bed",
    "_unused.bed",
    "_read_to_isoform.tsv",
    "_read_to_isoform.unique.tsv",
    "_unique_assignment.provenance.tsv",
    "_isoform_count.csv",
    "_desc.txt",
    "_class4.txt",
    "_fusion.txt",
    "_class12.txt",
    "_sqanti_structural_category.tsv",
    ".isoform_count.csv",
    ".isoform_usage.long.tsv",
    ".isoform_counts.matrix.tsv",
    ".isoform_usage.group.tsv",
    ".unique_assignment.provenance.tsv",
];

/// Keep biological gene directories disjoint from every top-level pipeline artifact.
///
/// Common ASCII gene IDs retain their spelling as a path key, so a valid ID can otherwise equal
/// a later summary or merged-output filename. Rejecting the complete reserved set up front avoids
/// a clean run failing only after expensive per-gene work has already completed.
pub(super) fn validate_pipeline_gene_namespace(
    genes: &[GeneId],
    prefix: Option<&SafePathComponent>,
) -> anyhow::Result<()> {
    let mut reserved = BATCH_TOP_LEVEL_NAMES
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<std::collections::HashSet<_>>();
    if let Some(prefix) = prefix {
        reserved.extend(
            PREFIX_TOP_LEVEL_SUFFIXES
                .iter()
                .map(|suffix| format!("{}{suffix}", prefix.as_str())),
        );
    }
    for gene in genes {
        let key = gene.path_key();
        if reserved.contains(key.as_str()) {
            if BATCH_TOP_LEVEL_NAMES.contains(&key.as_str()) {
                anyhow::bail!(
                    "gene {:?} maps to reserved top-level pipeline artifact name {:?} (fixed name); choose a different biological gene ID",
                    gene.as_str(),
                    key.as_str()
                );
            }
            anyhow::bail!(
                "gene {:?} maps to reserved top-level pipeline artifact name {:?} (prefix-scoped name); choose a different output prefix or biological gene ID",
                gene.as_str(),
                key.as_str()
            );
        }
    }
    Ok(())
}

/// Add a suffix without assuming that the base path has a UTF-8 extension.
pub(super) fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

/// Complete, containment-checked layout for batch-level flow artifacts.
#[derive(Clone, Debug)]
pub(super) struct FlowArtifactLayout {
    pub(super) gene_list: PathBuf,
    pub(super) gene_path_map: PathBuf,
    pub(super) isoform_bed: PathBuf,
    pub(super) unused_bed: PathBuf,
    pub(super) read_to_isoform_tsv: PathBuf,
    pub(super) unique_read_to_isoform_tsv: PathBuf,
    pub(super) unique_assignment_provenance_tsv: PathBuf,
    pub(super) count_csv: PathBuf,
    pub(super) desc_prefix: PathBuf,
}

impl FlowArtifactLayout {
    /// Build and validate every published path before any flow stage runs.
    pub(super) fn new(output_root: &Path, prefix: &SafePathComponent) -> anyhow::Result<Self> {
        let named = |suffix: &str| output_root.join(format!("{}{suffix}", prefix.as_str()));
        let layout = Self {
            gene_list: named("_gene.txt"),
            gene_path_map: named("_gene_paths.tsv"),
            isoform_bed: named("_isoform.bed"),
            unused_bed: named("_unused.bed"),
            read_to_isoform_tsv: named("_read_to_isoform.tsv"),
            unique_read_to_isoform_tsv: named("_read_to_isoform.unique.tsv"),
            unique_assignment_provenance_tsv: named("_unique_assignment.provenance.tsv"),
            count_csv: named("_isoform_count.csv"),
            desc_prefix: output_root.join(prefix.as_str()),
        };

        for path in [
            &layout.gene_list,
            &layout.gene_path_map,
            &layout.isoform_bed,
            &layout.unused_bed,
            &layout.read_to_isoform_tsv,
            &layout.unique_read_to_isoform_tsv,
            &layout.unique_assignment_provenance_tsv,
            &layout.count_csv,
        ] {
            ensure_destination_within(output_root, path)?;
        }
        for suffix in [
            "_desc.txt",
            "_class4.txt",
            "_fusion.txt",
            "_class12.txt",
            ".isoform_count.csv",
            ".isoform_usage.long.tsv",
            ".isoform_counts.matrix.tsv",
            ".isoform_usage.group.tsv",
        ] {
            ensure_destination_within(output_root, &append_suffix(&layout.desc_prefix, suffix))?;
        }
        Ok(layout)
    }
}

/// Required per-gene inputs for a count-only merge.
#[derive(Clone, Copy, Debug)]
pub(super) struct GeneArtifactRequirements<'a> {
    /// Mode-specific per-gene isoform suffix.
    pub(super) isoform_suffix: &'a str,
    /// Whether unique assignment requires the prepared read BED.
    pub(super) require_reads: bool,
}

/// Require a complete per-gene artifact set before any merged output is published.
pub(super) fn preflight_gene_artifacts(
    output_root: &Path,
    genes: &[String],
    requirements: GeneArtifactRequirements<'_>,
) -> anyhow::Result<()> {
    let mut missing = Vec::new();
    for gene in genes {
        let gene_id = GeneId::parse(gene)?;
        let mut required = vec![
            (
                "isoform",
                gene_artifact_path(output_root, &gene_id, requirements.isoform_suffix)?,
            ),
            (
                "unused",
                gene_artifact_path(output_root, &gene_id, "_unused.bed")?,
            ),
            (
                "read_to_isoform",
                gene_artifact_path(output_root, &gene_id, "_read_to_isoform.tsv")?,
            ),
        ];
        if requirements.require_reads {
            required.push((
                "prepared_reads",
                gene_artifact_path(output_root, &gene_id, "_nano.bed")?,
            ));
        }
        for (role, path) in required {
            match std::fs::symlink_metadata(&path) {
                Ok(metadata) if metadata.is_file() && !metadata.file_type().is_symlink() => {}
                Ok(_) => missing.push(format!(
                    "gene={gene:?} role={role} path={path:?} (not a regular file)"
                )),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                    missing.push(format!("gene={gene:?} role={role} path={path:?}"));
                }
                Err(error) => {
                    return Err(error)
                        .with_context(|| format!("inspect required gene artifact {path:?}"));
                }
            }
        }
    }
    if !missing.is_empty() {
        anyhow::bail!(
            "count-only preflight found {} missing or invalid required per-gene artifact(s): {}",
            missing.len(),
            missing.join("; ")
        );
    }
    Ok(())
}
