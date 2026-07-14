//! Counting, assignment, and description stages for merged flow artifacts.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Context;

use crate::annotate::desc::{describe, DescOpts};
use crate::count::count_by_read_to_isoform;
use crate::count::multi::MultiSampleOutputPaths;
use crate::flow::artifact_layout::append_suffix;
use crate::flow::artifact_manifest::atomic_copy;
use crate::flow::full::ClusterMode;
use crate::flow::path_key::{gene_artifact_path, GeneId};
use crate::io::bed::{read_bed12, BedError};
use crate::io::manifest::SampleRow;
use crate::model::Transcript;

static STAGING_SEQUENCE: AtomicU64 = AtomicU64::new(0);

struct StagingDir(PathBuf);

impl StagingDir {
    fn create(destination: &Path) -> anyhow::Result<Self> {
        let parent = destination
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        std::fs::create_dir_all(parent)
            .with_context(|| format!("create derived-output parent {parent:?}"))?;
        for _ in 0..1_000 {
            let sequence = STAGING_SEQUENCE.fetch_add(1, Ordering::Relaxed);
            let path = parent.join(format!(
                ".trackcluster-stage-{}-{sequence}",
                std::process::id()
            ));
            match std::fs::create_dir(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    return Err(error)
                        .with_context(|| format!("create derived-output staging dir {path:?}"));
                }
            }
        }
        anyhow::bail!("could not reserve a unique derived-output staging directory")
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for StagingDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Validate the complete staged set before publishing any file. Publication then uses one
/// sync-and-rename per destination. Portable filesystems cannot provide a transaction spanning
/// multiple names, but generation and staging failures preserve the complete previous set.
fn publish_staged_files(files: &[(PathBuf, PathBuf)]) -> anyhow::Result<()> {
    for (staged, _) in files {
        let metadata = std::fs::symlink_metadata(staged)
            .with_context(|| format!("validate staged derived output {staged:?}"))?;
        if metadata.file_type().is_symlink() || !metadata.is_file() {
            anyhow::bail!("staged derived output is not a regular file: {staged:?}");
        }
    }
    for (staged, destination) in files {
        atomic_copy(staged, destination).with_context(|| {
            format!("atomically publish derived output {staged:?} -> {destination:?}")
        })?;
    }
    Ok(())
}

fn scale_for_gene_field(gene_field: &str, scales: &HashMap<String, f64>) -> Option<f64> {
    let gene_field = gene_field.trim();
    if gene_field.is_empty() || gene_field == "none" {
        return None;
    }
    let mut scale = None;
    for gene in gene_field
        .split("||")
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
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

pub(super) fn read_bed12_records(path: &Path, kind: &str) -> anyhow::Result<Vec<Transcript>> {
    read_bed12(path)
        .with_context(|| format!("open {kind} {path:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse {kind} {path:?}"))
}

fn read_read_records(
    path: &Path,
    invalid_read_policy: crate::flow::config::InvalidReadPolicy,
) -> anyhow::Result<Vec<Transcript>> {
    let mut reader =
        read_bed12(path).with_context(|| format!("open prepared gene reads {path:?}"))?;
    let mut records = Vec::new();
    loop {
        let next = match invalid_read_policy {
            crate::flow::config::InvalidReadPolicy::Skip => reader.next_recovering_read(),
            crate::flow::config::InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse prepared gene reads {path:?}"))?;
        let Some(record) = next else {
            break;
        };
        records.push(record);
    }
    Ok(records)
}

pub(super) fn select_unique_read_to_isoform_by_gene(
    output_root: &Path,
    genes: &[String],
    cluster_mode: ClusterMode,
    options: crate::count::UniqueAssignmentOptions,
    invalid_read_policy: crate::flow::config::InvalidReadPolicy,
) -> anyhow::Result<Vec<(String, String)>> {
    let mut selected = Vec::new();
    for gene in genes {
        let gene_id = GeneId::parse(gene)?;
        let reads_path = gene_artifact_path(output_root, &gene_id, "_nano.bed")?;
        let isoform_path = gene_artifact_path(
            output_root,
            &gene_id,
            cluster_mode.per_gene_isoform_suffix(),
        )?;
        let mapping_path = gene_artifact_path(output_root, &gene_id, "_read_to_isoform.tsv")?;
        let required = [&reads_path, &isoform_path, &mapping_path];
        if let Some(missing) = required.iter().find(|path| !path.is_file()) {
            anyhow::bail!(
                "unique assignment requires every selected gene artifact; gene {gene:?} is missing {missing:?}"
            );
        }

        let reads = read_read_records(&reads_path, invalid_read_policy)?;
        let isoforms = read_bed12_records(&isoform_path, "per-gene isoforms")?;
        let read_to_isoform = crate::count::read_read_to_isoform_tsv(&mapping_path)
            .with_context(|| format!("read per-gene read_to_isoform {mapping_path:?}"))?;
        selected.extend(
            crate::count::select_unique_best_read_to_isoform_with_options(
                &reads,
                &isoforms,
                &read_to_isoform,
                options,
            )?,
        );
    }
    Ok(selected)
}

pub(super) fn run_count_and_desc(
    isoforms: &[Transcript],
    refs: &[Transcript],
    read_to_isoform: &[(String, String)],
    count_csv: &Path,
    desc_prefix: &Path,
    downsample_scales: Option<&HashMap<String, f64>>,
) -> anyhow::Result<()> {
    let mut counts = count_by_read_to_isoform(isoforms, read_to_isoform)?;
    if let Some(scales) = downsample_scales {
        for (record, isoform) in counts.iter_mut().zip(isoforms.iter()) {
            let gene_field = isoform.metadata().gene_id().unwrap_or("").trim();
            if let Some(scale) = scale_for_gene_field(gene_field, scales) {
                record.count *= scale;
            }
        }
    }
    let descriptions = describe(isoforms, refs, DescOpts::default());
    let staging = StagingDir::create(count_csv)?;
    let staged_count = staging.path().join("isoform_count.csv");
    crate::count::write_counts_csv(&staged_count, &counts)
        .with_context(|| format!("stage count csv {staged_count:?}"))?;
    let staged_desc = crate::annotate::desc_output::write_desc_outputs(
        &staging.path().join("description"),
        &descriptions,
    )?;
    crate::annotate::desc_output::remove_retired_description_outputs(desc_prefix)?;
    publish_staged_files(&[
        (staged_count, count_csv.to_path_buf()),
        (staged_desc.desc, append_suffix(desc_prefix, "_desc.txt")),
        (
            staged_desc.class4,
            append_suffix(desc_prefix, "_class4.txt"),
        ),
        (
            staged_desc.fusion,
            append_suffix(desc_prefix, "_fusion.txt"),
        ),
        (
            staged_desc.class12,
            append_suffix(desc_prefix, "_class12.txt"),
        ),
    ])
}

pub(super) fn run_count_multi_atomic(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    out_prefix: &Path,
    downsample_scales: Option<&HashMap<String, f64>>,
) -> anyhow::Result<MultiSampleOutputPaths> {
    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }
    let mut result =
        crate::count::multi::count_multi_by_read_to_isoform(isoforms, read_to_isoform, sample_rows)
            .with_context(|| format!("count-multi {out_prefix:?}"))?;
    if let Some(scales) = downsample_scales {
        for row in &mut result.matrix_rows {
            let Some(scale) = scale_for_gene_field(&row.gene, scales) else {
                continue;
            };
            for count in &mut row.counts {
                *count *= scale;
            }
        }
        for row in &mut result.long_rows {
            let Some(scale) = scale_for_gene_field(&row.gene, scales) else {
                continue;
            };
            row.count *= scale;
            row.gene_total *= scale;
        }
        for row in &mut result.group_rows {
            let Some(scale) = scale_for_gene_field(&row.gene, scales) else {
                continue;
            };
            row.count *= scale;
            row.gene_total *= scale;
        }
    }

    let count_csv = append_suffix(out_prefix, ".isoform_count.csv");
    let long_tsv = append_suffix(out_prefix, ".isoform_usage.long.tsv");
    let matrix_tsv = append_suffix(out_prefix, ".isoform_counts.matrix.tsv");
    let staging = StagingDir::create(&count_csv)?;
    let staged_count = staging.path().join("isoform_count.csv");
    let staged_long = staging.path().join("isoform_usage.long.tsv");
    let staged_matrix = staging.path().join("isoform_counts.matrix.tsv");
    let count_records =
        crate::count::multi::total_count_records_from_matrix_rows(&result.matrix_rows);
    crate::count::write_counts_csv(&staged_count, &count_records)
        .with_context(|| format!("stage aggregate count output {staged_count:?}"))?;
    let include_group = sample_rows.iter().any(|sample| sample.group.is_some());
    crate::count::multi::write_usage_long_tsv(&staged_long, &result.long_rows, include_group)
        .with_context(|| format!("stage long output {staged_long:?}"))?;
    crate::count::multi::write_counts_matrix_tsv(&staged_matrix, &result.matrix_rows, sample_rows)
        .with_context(|| format!("stage matrix output {staged_matrix:?}"))?;
    let mut files = vec![
        (staged_count, count_csv.clone()),
        (staged_long, long_tsv.clone()),
        (staged_matrix, matrix_tsv.clone()),
    ];
    let group_tsv = if include_group {
        let destination = append_suffix(out_prefix, ".isoform_usage.group.tsv");
        let staged = staging.path().join("isoform_usage.group.tsv");
        crate::count::multi::write_group_usage_tsv(&staged, &result.group_rows)
            .with_context(|| format!("stage group output {staged:?}"))?;
        files.push((staged, destination.clone()));
        Some(destination)
    } else {
        None
    };
    publish_staged_files(&files)?;

    if group_tsv.is_none() {
        let stale_group = append_suffix(out_prefix, ".isoform_usage.group.tsv");
        if stale_group.exists() {
            std::fs::remove_file(&stale_group)
                .with_context(|| format!("remove stale group output {stale_group:?}"))?;
        }
    }
    Ok(MultiSampleOutputPaths {
        count_csv,
        long_tsv,
        matrix_tsv,
        group_tsv,
        unique_assignment_provenance_tsv: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fused_gene_field_is_not_scaled_ambiguously() {
        let scales = HashMap::from([("A".to_owned(), 2.0), ("B".to_owned(), 3.0)]);
        assert_eq!(scale_for_gene_field("A", &scales), Some(2.0));
        assert_eq!(scale_for_gene_field("A||B", &scales), None);
        assert_eq!(scale_for_gene_field("none", &scales), None);
    }

    #[test]
    fn staging_validation_failure_preserves_every_existing_destination() {
        let dir = std::env::temp_dir().join(format!(
            "trackcluster-counting-stage-{}-{}",
            std::process::id(),
            STAGING_SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let staged = dir.join("staged-one");
        let missing = dir.join("missing-two");
        let first = dir.join("first");
        let second = dir.join("second");
        std::fs::write(&staged, "new-one").unwrap();
        std::fs::write(&first, "old-one").unwrap();
        std::fs::write(&second, "old-two").unwrap();
        let error = publish_staged_files(&[(staged, first.clone()), (missing, second.clone())])
            .expect_err("missing staged file must fail preflight");
        assert!(format!("{error:#}").contains("missing-two"));
        assert_eq!(std::fs::read_to_string(&first).unwrap(), "old-one");
        assert_eq!(std::fs::read_to_string(&second).unwrap(), "old-two");
        std::fs::remove_dir_all(dir).unwrap();
    }
}
