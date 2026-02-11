use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Context;

use crate::annotate::addgene::{add_gene, AddGeneOpts};
use crate::io::bed::{read_bed12, write_bed12_to_writer, BedError};
use crate::model::Transcript;

const GENE_NAME_COL: usize = 5;

#[derive(Clone, Debug)]
pub struct PrepareDirResult {
    pub genes: Vec<String>,
    pub dedup_reads: usize,
    pub novel_reads: usize,
}

fn exon_len(tx: &Transcript) -> u32 {
    tx.exons.iter().map(|exon| exon.len()).sum()
}

fn dedup_longest_by_name(reads: &[Transcript]) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    let mut pos: HashMap<String, usize> = HashMap::new();

    for read in reads {
        match pos.get(&read.name).copied() {
            None => {
                out.push(read.clone());
                pos.insert(read.name.clone(), out.len() - 1);
            }
            Some(existing_idx) => {
                let existing_len = exon_len(&out[existing_idx]);
                let candidate_len = exon_len(read);
                if candidate_len > existing_len {
                    out[existing_idx] = read.clone();
                }
            }
        }
    }

    out
}

fn gene_field(tx: &Transcript) -> &str {
    tx.extra_fields
        .get(GENE_NAME_COL)
        .map(|value| value.as_str())
        .unwrap_or(tx.name.as_str())
}

fn genes(tx: &Transcript) -> impl Iterator<Item = &str> {
    gene_field(tx)
        .split("||")
        .map(str::trim)
        .filter(|g| !g.is_empty() && *g != "none")
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
    write_bed12_to_writer(&mut writer, indices.iter().map(|&idx| &transcripts[idx])).map_err(
        |source| BedError::IoWrite {
            path: path.to_path_buf(),
            source,
        },
    )?;
    Ok(())
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

pub fn prepare_dir_from_paths(
    reads_bed: &Path,
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let reads_raw: Vec<Transcript> = read_bed12(reads_bed)
        .with_context(|| format!("open reads {reads_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reads {reads_bed:?}"))?;

    let refs: Vec<Transcript> = read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    // Step 1 (Python: list_to_dic + write prefix_dedup.bed)
    let reads_dedup = dedup_longest_by_name(&reads_raw);
    let dedup_path = output_root.join(format!("{prefix}_dedup.bed"));
    crate::io::bed::write_bed12(&dedup_path, reads_dedup.iter())
        .with_context(|| format!("write {dedup_path:?}"))?;

    // Step 2 (Python: bedtools intersect + tracklist_add_gene)
    let reads_annotated = add_gene(&reads_dedup, &refs, addgene_opts);

    // Group references by gene (Python: group_bigg_by_gene(bigg_ref))
    let mut ref_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tx) in refs.iter().enumerate() {
        for gene in genes(tx) {
            ref_by_gene.entry(gene.to_owned()).or_default().push(idx);
        }
    }

    // Group reads by gene (Python: group_bigg_by_gene(bigg_new))
    let mut reads_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    let mut novel_indices: Vec<usize> = Vec::new();
    for (idx, tx) in reads_annotated.iter().enumerate() {
        let field = gene_field(tx);
        if field.trim() == "none" || field.trim().is_empty() {
            novel_indices.push(idx);
            continue;
        }
        for gene in genes(tx) {
            reads_by_gene.entry(gene.to_owned()).or_default().push(idx);
        }
    }

    // Write novel reads file (Python: prefix_novel.bed)
    let novel_path = output_root.join(format!("{prefix}_novel.bed"));
    write_bed12_indices(&novel_path, &reads_annotated, &novel_indices)
        .with_context(|| format!("write {novel_path:?}"))?;

    // Materialize gene list (Python: name2file(genename_l, prefix_gene.txt))
    let mut genes: Vec<String> = reads_by_gene.keys().cloned().collect();
    genes.sort();

    let gene_list_path = output_root.join(format!("{prefix}_gene.txt"));
    let mut gene_list = BufWriter::new(
        File::create(&gene_list_path).with_context(|| format!("write {gene_list_path:?}"))?,
    );
    for gene in &genes {
        writeln!(gene_list, "{gene}")?;
    }

    // Create per-gene folders and write inputs (Python: write {gene}_gff.bed and {gene}_nano.bed)
    for gene in &genes {
        let gene_dir = output_root.join(gene);
        fs::create_dir_all(&gene_dir).with_context(|| format!("create {gene_dir:?}"))?;

        let mut ref_indices = ref_by_gene.get(gene).cloned().unwrap_or_default();
        sort_indices_by_coord(&refs, &mut ref_indices);
        let ref_path = gene_dir.join(format!("{gene}_gff.bed"));
        write_bed12_indices(&ref_path, &refs, &ref_indices)
            .with_context(|| format!("write {ref_path:?}"))?;

        let mut read_indices = reads_by_gene.get(gene).cloned().unwrap_or_default();
        sort_indices_by_coord(&reads_annotated, &mut read_indices);
        let reads_path = gene_dir.join(format!("{gene}_nano.bed"));
        write_bed12_indices(&reads_path, &reads_annotated, &read_indices)
            .with_context(|| format!("write {reads_path:?}"))?;
    }

    Ok(PrepareDirResult {
        genes,
        dedup_reads: reads_dedup.len(),
        novel_reads: novel_indices.len(),
    })
}
