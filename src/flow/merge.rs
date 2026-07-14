//! Atomic, deterministic merging of per-gene flow artifacts.

use std::collections::HashMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::flow::artifact_manifest::atomic_write_with;
use crate::flow::path_key::{ensure_destination_within, gene_artifact_path, GeneId};
use crate::io::bed::{read_bed12, BedError};
use crate::model::Transcript;

pub(super) fn merge_files(inputs: &[PathBuf], out: &Path) -> anyhow::Result<()> {
    atomic_write_with(out, |writer| write_merged_file(inputs, writer, out))
        .with_context(|| format!("atomically publish merged output {out:?}"))
}

pub(super) fn merge_isoform_files(inputs: &[PathBuf], out: &Path) -> anyhow::Result<()> {
    atomic_write_with(out, |temporary| {
        let mut isoforms = Vec::new();
        let mut index_by_id = HashMap::<String, usize>::new();
        for input in inputs {
            for isoform in read_records(input, "per-gene isoform")? {
                if let Some(&existing_index) = index_by_id.get(&isoform.name) {
                    if isoforms[existing_index] == isoform {
                        // A reference annotated to multiple genes is copied to
                        // each gene folder. Collapse only the byte-equivalent
                        // biological record; conflicting reuse of an ID is an
                        // integrity error.
                        continue;
                    }
                    anyhow::bail!(
                        "duplicate isoform id {:?} describes conflicting records in per-gene outputs",
                        isoform.name
                    );
                }
                index_by_id.insert(isoform.name.clone(), isoforms.len());
                isoforms.push(isoform);
            }
        }
        crate::identity::validate_isoform_ids(&isoforms).with_context(|| {
            format!("validate globally unique isoform IDs before publishing {out:?}")
        })?;
        crate::io::bed::write_bed12_to_writer(temporary, isoforms.iter())
            .with_context(|| format!("write merged isoform catalog {out:?}"))?;
        Ok(())
    })
    .with_context(|| format!("atomically publish validated merged isoform catalog {out:?}"))
}

fn read_records(path: &Path, kind: &str) -> anyhow::Result<Vec<Transcript>> {
    read_bed12(path)
        .with_context(|| format!("open {kind} {path:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse {kind} {path:?}"))
}

fn write_merged_file<W: Write>(
    inputs: &[PathBuf],
    output: &mut W,
    out: &Path,
) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(output);
    let mut buffer = vec![0u8; 1024 * 1024];
    for input in inputs {
        let mut reader = std::io::BufReader::new(
            fs::File::open(input).with_context(|| format!("open {input:?}"))?,
        );
        let mut saw_bytes = false;
        let mut last_byte = b'\n';
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
    writer
        .flush()
        .with_context(|| format!("flush merged output {out:?}"))?;
    Ok(())
}

fn existing_gene_artifacts(
    output_root: &Path,
    genes: &[String],
    per_gene_suffix: &str,
) -> anyhow::Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for gene in genes {
        let gene = GeneId::parse(gene)?;
        let path = gene_artifact_path(output_root, &gene, per_gene_suffix)?;
        let metadata = fs::symlink_metadata(&path)
            .with_context(|| format!("required per-gene merge artifact is missing {path:?}"))?;
        if metadata.file_type().is_symlink() || !metadata.is_file() {
            anyhow::bail!("required per-gene merge artifact is not a regular file: {path:?}");
        }
        paths.push(path);
    }
    paths.sort();
    Ok(paths)
}

pub(super) fn merge_gene_outputs(
    output_root: &Path,
    genes: &[String],
    per_gene_suffix: &str,
    merged_out: &Path,
) -> anyhow::Result<()> {
    ensure_destination_within(output_root, merged_out)?;
    merge_files(
        &existing_gene_artifacts(output_root, genes, per_gene_suffix)?,
        merged_out,
    )
}

pub(super) fn merge_gene_isoform_outputs(
    output_root: &Path,
    genes: &[String],
    per_gene_suffix: &str,
    merged_out: &Path,
) -> anyhow::Result<()> {
    ensure_destination_within(output_root, merged_out)?;
    merge_isoform_files(
        &existing_gene_artifacts(output_root, genes, per_gene_suffix)?,
        merged_out,
    )
}
