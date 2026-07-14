use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::count::parse_subreads;
use crate::count::CountRecord;
use crate::io::manifest::{read_manifest_tsv, SampleRow};
use crate::model::Transcript;
use crate::sample::{split_tagged_read_name, tagged_read_name};

const EPSILON: f64 = 1e-12;

#[derive(Clone, Debug, PartialEq)]
pub struct UsageLongRow {
    pub gene: String,
    pub isoform_id: String,
    pub sample: String,
    pub group: Option<String>,
    pub count: f64,
    pub proportion: f64,
    pub gene_total: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct UsageMatrixRow {
    pub gene: String,
    pub isoform_id: String,
    pub counts: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GroupUsageRow {
    pub gene: String,
    pub isoform_id: String,
    pub group: String,
    pub count: f64,
    pub proportion: f64,
    pub gene_total: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultiSampleCountResult {
    pub matrix_rows: Vec<UsageMatrixRow>,
    pub long_rows: Vec<UsageLongRow>,
    pub group_rows: Vec<GroupUsageRow>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiSampleOutputPaths {
    pub count_csv: PathBuf,
    pub long_tsv: PathBuf,
    pub matrix_tsv: PathBuf,
    pub group_tsv: Option<PathBuf>,
    pub unique_assignment_provenance_tsv: Option<PathBuf>,
}

fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

fn unique_assignment_provenance_path(out_prefix: &Path) -> PathBuf {
    append_suffix(out_prefix, ".unique_assignment.provenance.tsv")
}

fn remove_optional_output(path: &Path, kind: &str) -> anyhow::Result<()> {
    match std::fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("remove stale {kind} output {path:?}")),
    }
}

fn isoform_gene(tx: &Transcript) -> String {
    tx.metadata()
        .gene_id()
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
        .unwrap_or("none")
        .to_owned()
}

fn group_order(samples: &[SampleRow]) -> Vec<String> {
    let mut seen: HashSet<&str> = HashSet::new();
    let mut groups: Vec<String> = Vec::new();
    for sample in samples {
        let Some(group) = sample.group.as_deref() else {
            continue;
        };
        if seen.insert(group) {
            groups.push(group.to_owned());
        }
    }
    groups
}

fn usage_from_matrix_rows(
    mut matrix_rows: Vec<UsageMatrixRow>,
    samples: &[SampleRow],
) -> MultiSampleCountResult {
    matrix_rows.sort_by(|a, b| {
        a.gene
            .cmp(&b.gene)
            .then_with(|| a.isoform_id.cmp(&b.isoform_id))
    });

    let mut gene_totals: HashMap<String, Vec<f64>> = HashMap::new();
    for row in &matrix_rows {
        let totals = gene_totals
            .entry(row.gene.clone())
            .or_insert_with(|| vec![0.0f64; samples.len()]);
        for (idx, count) in row.counts.iter().copied().enumerate() {
            totals[idx] += count;
        }
    }

    let mut long_rows: Vec<UsageLongRow> = Vec::new();
    for row in &matrix_rows {
        let totals = gene_totals
            .get(&row.gene)
            .expect("gene totals generated from same row set");
        for (sample_idx, count) in row.counts.iter().copied().enumerate() {
            if count <= EPSILON {
                continue;
            }
            let gene_total = totals[sample_idx];
            let proportion = if gene_total > 0.0 {
                count / gene_total
            } else {
                0.0
            };
            long_rows.push(UsageLongRow {
                gene: row.gene.clone(),
                isoform_id: row.isoform_id.clone(),
                sample: samples[sample_idx].sample.clone(),
                group: samples[sample_idx].group.clone(),
                count,
                proportion,
                gene_total,
            });
        }
    }

    let groups = group_order(samples);
    let group_rows = if groups.is_empty() {
        Vec::new()
    } else {
        let group_to_idx: HashMap<&str, usize> = groups
            .iter()
            .enumerate()
            .map(|(idx, group)| (group.as_str(), idx))
            .collect();
        let mut group_totals: HashMap<String, Vec<f64>> = HashMap::new();
        let mut grouped_rows: Vec<(String, String, Vec<f64>)> = Vec::new();

        for row in &matrix_rows {
            let mut counts = vec![0.0f64; groups.len()];
            for (sample_idx, sample_count) in row.counts.iter().copied().enumerate() {
                let Some(group) = samples[sample_idx].group.as_deref() else {
                    continue;
                };
                let group_idx = *group_to_idx
                    .get(group)
                    .expect("group index map derived from sample groups");
                counts[group_idx] += sample_count;
            }

            let totals = group_totals
                .entry(row.gene.clone())
                .or_insert_with(|| vec![0.0f64; groups.len()]);
            for (idx, count) in counts.iter().copied().enumerate() {
                totals[idx] += count;
            }
            grouped_rows.push((row.gene.clone(), row.isoform_id.clone(), counts));
        }

        let mut out: Vec<GroupUsageRow> = Vec::new();
        for (gene, isoform_id, counts) in grouped_rows {
            let totals = group_totals
                .get(&gene)
                .expect("group totals generated from same row set");
            for (group_idx, count) in counts.into_iter().enumerate() {
                if count <= EPSILON {
                    continue;
                }
                let gene_total = totals[group_idx];
                let proportion = if gene_total > 0.0 {
                    count / gene_total
                } else {
                    0.0
                };
                out.push(GroupUsageRow {
                    gene: gene.clone(),
                    isoform_id: isoform_id.clone(),
                    group: groups[group_idx].clone(),
                    count,
                    proportion,
                    gene_total,
                });
            }
        }
        out
    };

    MultiSampleCountResult {
        matrix_rows,
        long_rows,
        group_rows,
    }
}

pub fn count_multi_by_subreads(
    isoforms: &[Transcript],
    _references: &[Transcript],
    samples: &[SampleRow],
) -> anyhow::Result<MultiSampleCountResult> {
    if samples.is_empty() {
        anyhow::bail!("count-multi requires at least one sample");
    }

    let mut mapping = Vec::new();
    for isoform in isoforms {
        for read_name in parse_subreads(isoform)? {
            mapping.push((read_name, isoform.name.clone()));
        }
    }
    count_multi_by_read_to_isoform(isoforms, &mapping, samples)
}

pub fn count_multi_by_read_to_isoform(
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    samples: &[SampleRow],
) -> anyhow::Result<MultiSampleCountResult> {
    if samples.is_empty() {
        anyhow::bail!("count-multi requires at least one sample");
    }
    crate::identity::validate_isoform_ids(isoforms)?;

    let isoform_ids: HashSet<&str> = isoforms.iter().map(|tx| tx.name.as_str()).collect();
    let mut unique_pairs: std::collections::BTreeSet<(&str, &str)> =
        std::collections::BTreeSet::new();
    for (read_name, isoform_id) in read_to_isoform {
        if read_name.is_empty() {
            anyhow::bail!("read_to_isoform contains an empty read id");
        }
        if !isoform_ids.contains(isoform_id.as_str()) {
            anyhow::bail!(
                "read_to_isoform references isoform id {isoform_id:?} that is missing from isoform BED"
            );
        }
        unique_pairs.insert((read_name.as_str(), isoform_id.as_str()));
    }

    let sample_to_idx: HashMap<&str, usize> = samples
        .iter()
        .enumerate()
        .map(|(idx, row)| (row.sample.as_str(), idx))
        .collect();

    let mut read_occurrence: HashMap<&str, u32> = HashMap::new();
    let mut read_to_sample_idx: HashMap<&str, usize> = HashMap::new();
    for &(read_name, _isoform_id) in &unique_pairs {
        let (sample_name, _raw_read_id) = split_tagged_read_name(read_name).ok_or_else(|| {
            anyhow::anyhow!(
                "read id {read_name:?} is missing sample prefix; expected format '<sample>::<read_id>'"
            )
        })?;
        let sample_idx = sample_to_idx.get(sample_name).copied().ok_or_else(|| {
            anyhow::anyhow!(
                "read id {read_name:?} references unknown sample {sample_name:?} (not present in manifest)"
            )
        })?;

        match read_to_sample_idx.get(read_name).copied() {
            Some(existing) if existing != sample_idx => {
                anyhow::bail!(
                    "read id {read_name:?} maps to multiple samples ({:?} and {:?})",
                    samples[existing].sample,
                    samples[sample_idx].sample
                );
            }
            Some(_) => {}
            None => {
                read_to_sample_idx.insert(read_name, sample_idx);
            }
        }

        *read_occurrence.entry(read_name).or_insert(0) += 1;
    }

    let isoform_to_idx: HashMap<&str, usize> = isoforms
        .iter()
        .enumerate()
        .map(|(idx, isoform)| (isoform.name.as_str(), idx))
        .collect();
    let mut counts_by_isoform: Vec<Vec<f64>> = vec![vec![0.0f64; samples.len()]; isoforms.len()];

    for &(read_name, isoform_id) in &unique_pairs {
        let Some(sample_idx) = read_to_sample_idx.get(read_name).copied() else {
            anyhow::bail!("internal error: sample lookup missing for read id {read_name:?}");
        };
        let Some(isoform_idx) = isoform_to_idx.get(isoform_id).copied() else {
            anyhow::bail!(
                "read_to_isoform references isoform id {isoform_id:?} that is missing from isoform BED"
            );
        };

        let denom = read_occurrence.get(read_name).copied().unwrap_or(0);
        if denom > 0 {
            counts_by_isoform[isoform_idx][sample_idx] += 1.0f64 / denom as f64;
        }
    }

    let matrix_rows: Vec<UsageMatrixRow> = isoforms
        .iter()
        .zip(counts_by_isoform)
        .map(|(isoform, counts)| UsageMatrixRow {
            gene: isoform_gene(isoform),
            isoform_id: isoform.name.clone(),
            counts,
        })
        .collect();

    Ok(usage_from_matrix_rows(matrix_rows, samples))
}

pub fn write_usage_long_tsv(
    path: &Path,
    rows: &[UsageLongRow],
    include_group: bool,
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    write_usage_long_to_writer(&mut writer, rows, include_group)
}

fn write_usage_long_to_writer<W: Write>(
    writer: &mut W,
    rows: &[UsageLongRow],
    include_group: bool,
) -> Result<(), std::io::Error> {
    if include_group {
        writeln!(
            writer,
            "gene\tisoform_id\tsample\tgroup\tcount\tproportion\tgene_total"
        )?;
        for row in rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.gene,
                row.isoform_id,
                row.sample,
                row.group.as_deref().unwrap_or(""),
                row.count,
                row.proportion,
                row.gene_total
            )?;
        }
    } else {
        writeln!(
            writer,
            "gene\tisoform_id\tsample\tcount\tproportion\tgene_total"
        )?;
        for row in rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                row.gene, row.isoform_id, row.sample, row.count, row.proportion, row.gene_total
            )?;
        }
    }
    writer.flush()
}

pub fn write_counts_matrix_tsv(
    path: &Path,
    rows: &[UsageMatrixRow],
    samples: &[SampleRow],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    write_counts_matrix_to_writer(&mut writer, rows, samples)
}

fn write_counts_matrix_to_writer<W: Write>(
    writer: &mut W,
    rows: &[UsageMatrixRow],
    samples: &[SampleRow],
) -> Result<(), std::io::Error> {
    write!(writer, "gene\tisoform_id")?;
    for sample in samples {
        write!(writer, "\t{}", sample.sample)?;
    }
    writeln!(writer)?;

    for row in rows {
        write!(writer, "{}\t{}", row.gene, row.isoform_id)?;
        for count in &row.counts {
            write!(writer, "\t{}", count)?;
        }
        writeln!(writer)?;
    }
    writer.flush()
}

pub fn total_count_records_from_matrix_rows(rows: &[UsageMatrixRow]) -> Vec<CountRecord> {
    rows.iter()
        .map(|row| CountRecord {
            gene: row.gene.clone(),
            isoform_id: row.isoform_id.clone(),
            count: row.counts.iter().sum(),
        })
        .collect()
}

pub fn write_group_usage_tsv(path: &Path, rows: &[GroupUsageRow]) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    write_group_usage_to_writer(&mut writer, rows)
}

fn write_group_usage_to_writer<W: Write>(
    writer: &mut W,
    rows: &[GroupUsageRow],
) -> Result<(), std::io::Error> {
    writeln!(
        writer,
        "gene\tisoform_id\tgroup\tcount\tproportion\tgene_total"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.gene, row.isoform_id, row.group, row.count, row.proportion, row.gene_total
        )?;
    }
    writer.flush()
}

pub fn run_count_multi_from_paths(
    manifest: &Path,
    reference: &Path,
    isoform: &Path,
    out_prefix: &Path,
) -> anyhow::Result<MultiSampleOutputPaths> {
    let sample_rows = read_manifest_tsv(manifest)?;
    let isoforms: Vec<Transcript> = crate::io::bed::read_bed12(isoform)
        .with_context(|| format!("open isoform {isoform:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse isoform {isoform:?}"))?;
    let refs: Vec<Transcript> = crate::io::bed::read_bed12(reference)
        .with_context(|| format!("open reference {reference:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reference {reference:?}"))?;

    run_count_multi(&sample_rows, &isoforms, &refs, out_prefix)
}

pub fn write_count_multi_outputs(
    sample_rows: &[SampleRow],
    result: &MultiSampleCountResult,
    out_prefix: &Path,
) -> anyhow::Result<MultiSampleOutputPaths> {
    let include_group = sample_rows.iter().any(|sample| sample.group.is_some());

    let count_csv = append_suffix(out_prefix, ".isoform_count.csv");
    let long_tsv = append_suffix(out_prefix, ".isoform_usage.long.tsv");
    let matrix_tsv = append_suffix(out_prefix, ".isoform_counts.matrix.tsv");
    let count_records = total_count_records_from_matrix_rows(&result.matrix_rows);
    crate::flow::artifact_manifest::atomic_write_with(&count_csv, |writer| {
        crate::count::write_counts_csv_to_writer(writer, &count_records).map_err(Into::into)
    })
    .with_context(|| format!("write aggregate count output {count_csv:?}"))?;
    crate::flow::artifact_manifest::atomic_write_with(&long_tsv, |writer| {
        write_usage_long_to_writer(writer, &result.long_rows, include_group).map_err(Into::into)
    })
    .with_context(|| format!("write long output {long_tsv:?}"))?;
    crate::flow::artifact_manifest::atomic_write_with(&matrix_tsv, |writer| {
        write_counts_matrix_to_writer(writer, &result.matrix_rows, sample_rows).map_err(Into::into)
    })
    .with_context(|| format!("write matrix output {matrix_tsv:?}"))?;

    let group_tsv = if include_group {
        let path = append_suffix(out_prefix, ".isoform_usage.group.tsv");
        crate::flow::artifact_manifest::atomic_write_with(&path, |writer| {
            write_group_usage_to_writer(writer, &result.group_rows).map_err(Into::into)
        })
        .with_context(|| format!("write group output {path:?}"))?;
        Some(path)
    } else {
        remove_optional_output(
            &append_suffix(out_prefix, ".isoform_usage.group.tsv"),
            "group-usage",
        )?;
        None
    };

    // Unique-mode callers publish fresh provenance after this common output set succeeds.
    remove_optional_output(
        &unique_assignment_provenance_path(out_prefix),
        "unique-assignment provenance",
    )?;

    Ok(MultiSampleOutputPaths {
        count_csv,
        long_tsv,
        matrix_tsv,
        group_tsv,
        unique_assignment_provenance_tsv: None,
    })
}

pub fn run_count_multi(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    references: &[Transcript],
    out_prefix: &Path,
) -> anyhow::Result<MultiSampleOutputPaths> {
    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }

    let result = count_multi_by_subreads(isoforms, references, sample_rows)?;
    write_count_multi_outputs(sample_rows, &result, out_prefix)
}

pub fn run_count_multi_from_read_to_isoform(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    out_prefix: &Path,
) -> anyhow::Result<MultiSampleOutputPaths> {
    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }

    let result = count_multi_by_read_to_isoform(isoforms, read_to_isoform, sample_rows)?;
    write_count_multi_outputs(sample_rows, &result, out_prefix)
}

pub fn read_tagged_sample_reads(sample_rows: &[SampleRow]) -> anyhow::Result<Vec<Transcript>> {
    let mut reads = Vec::new();
    for row in sample_rows {
        let mut sample_reads: Vec<Transcript> = crate::io::bed::read_bed12(&row.reads)
            .with_context(|| format!("open reads {:?}", row.reads))?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
            .with_context(|| format!("parse reads {:?}", row.reads))?;
        for read in &mut sample_reads {
            read.name = tagged_read_name(&row.sample, &read.name);
        }
        reads.extend(sample_reads);
    }
    Ok(reads)
}

pub fn run_count_multi_from_read_to_isoform_unique(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    reads: &[Transcript],
    read_to_isoform: &[(String, String)],
    out_prefix: &Path,
) -> anyhow::Result<MultiSampleOutputPaths> {
    run_count_multi_from_read_to_isoform_unique_with_options(
        sample_rows,
        isoforms,
        reads,
        read_to_isoform,
        out_prefix,
        crate::count::UniqueAssignmentOptions::default(),
    )
}

pub fn run_count_multi_from_read_to_isoform_unique_with_options(
    sample_rows: &[SampleRow],
    isoforms: &[Transcript],
    reads: &[Transcript],
    read_to_isoform: &[(String, String)],
    out_prefix: &Path,
    options: crate::count::UniqueAssignmentOptions,
) -> anyhow::Result<MultiSampleOutputPaths> {
    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }

    let unique_pairs = crate::count::select_unique_best_read_to_isoform_with_options(
        reads,
        isoforms,
        read_to_isoform,
        options,
    )?;
    let result = count_multi_by_read_to_isoform(isoforms, &unique_pairs, sample_rows)?;
    let mut outputs = write_count_multi_outputs(sample_rows, &result, out_prefix)?;
    let provenance = unique_assignment_provenance_path(out_prefix);
    crate::flow::artifact_manifest::atomic_write_with(&provenance, |writer| {
        crate::count::write_unique_assignment_provenance_to_writer(writer, options)
            .map_err(Into::into)
    })
    .with_context(|| format!("write unique-assignment provenance {provenance:?}"))?;
    outputs.unique_assignment_provenance_tsv = Some(provenance);
    Ok(outputs)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

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
    fn all_tsv_writers_propagate_final_flush_errors() {
        let mut writer = FlushFails;
        assert!(write_usage_long_to_writer(&mut writer, &[], false).is_err());

        let mut writer = FlushFails;
        assert!(write_counts_matrix_to_writer(&mut writer, &[], &[]).is_err());

        let mut writer = FlushFails;
        assert!(write_group_usage_to_writer(&mut writer, &[]).is_err());
    }

    #[test]
    fn provenance_suffix_preserves_dotted_prefix() {
        assert_eq!(
            unique_assignment_provenance_path(Path::new("out/result.v1")),
            PathBuf::from("out/result.v1.unique_assignment.provenance.tsv")
        );
    }

    #[test]
    fn grouped_outputs_include_an_empty_group_table() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time before Unix epoch")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "trackcluster-count-multi-empty-group-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("result");
        let result = MultiSampleCountResult {
            matrix_rows: Vec::new(),
            long_rows: Vec::new(),
            group_rows: Vec::new(),
        };

        let outputs = write_count_multi_outputs(&sample_rows(), &result, &prefix).unwrap();
        let group = outputs.group_tsv.expect("grouped manifest has group table");
        assert_eq!(
            fs::read_to_string(group).unwrap(),
            "gene\tisoform_id\tgroup\tcount\tproportion\tgene_total\n"
        );

        fs::write(unique_assignment_provenance_path(&prefix), "stale\n").unwrap();
        let ungrouped_samples = sample_rows()
            .into_iter()
            .map(|mut sample| {
                sample.group = None;
                sample
            })
            .collect::<Vec<_>>();
        let outputs = write_count_multi_outputs(&ungrouped_samples, &result, &prefix).unwrap();
        assert!(outputs.group_tsv.is_none());
        assert!(!append_suffix(&prefix, ".isoform_usage.group.tsv").exists());
        assert!(!unique_assignment_provenance_path(&prefix).exists());

        fs::remove_dir_all(dir).unwrap();
    }

    fn make_tx(name: &str, name2: &str, gene: &str) -> Transcript {
        make_tx_with_exons(name, name2, gene, &[(100, 200)])
    }

    #[test]
    fn reading_sample_reads_always_tags_the_raw_identifier() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time before Unix epoch")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "trackcluster-count-multi-tagged-raw-id-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let reads_path = dir.join("reads.bed");
        let read = make_tx("S1::r1", "none", "GENEA");
        crate::io::bed::write_bed12(&reads_path, std::iter::once(&read)).unwrap();
        let rows = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: reads_path,
        }];

        let reads = read_tagged_sample_reads(&rows).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].name, "S1::S1::r1");

        fs::remove_dir_all(dir).unwrap();
    }

    fn make_tx_with_exons(name: &str, name2: &str, gene: &str, exons: &[(u32, u32)]) -> Transcript {
        let tx_start = exons.iter().map(|(start, _)| *start).min().unwrap();
        let tx_end = exons.iter().map(|(_, end)| *end).max().unwrap();
        let exons = exons
            .iter()
            .map(|(start, end)| Interval::new(Coord::new(*start), Coord::new(*end)).unwrap())
            .collect::<Vec<_>>();
        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    name2.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "-1,".to_owned(),
                    "isoform_anno".to_owned(),
                    gene.to_owned(),
                ],
            },
        )
        .unwrap()
    }

    fn sample_rows() -> Vec<SampleRow> {
        vec![
            SampleRow {
                sample: "S1".to_owned(),
                group: Some("control".to_owned()),
                reads: PathBuf::from("S1.reads.bed"),
            },
            SampleRow {
                sample: "S2".to_owned(),
                group: Some("treated".to_owned()),
                reads: PathBuf::from("S2.reads.bed"),
            },
        ]
    }

    #[test]
    fn counts_by_sample_and_gene_proportions() {
        let references = vec![make_tx("ref_a", "ref_a", "GENEA")];
        let isoforms = vec![
            make_tx("iso1", "S1::r1,S1::r2,S2::r3,|0", "GENEA"),
            make_tx("iso2", "S1::r2,S2::r4,|0", "GENEA"),
        ];

        let result = count_multi_by_subreads(&isoforms, &references, &sample_rows()).unwrap();

        assert_eq!(result.matrix_rows.len(), 2);
        assert_eq!(result.matrix_rows[0].isoform_id, "iso1");
        assert!((result.matrix_rows[0].counts[0] - 1.5).abs() < 1e-9);
        assert!((result.matrix_rows[0].counts[1] - 1.0).abs() < 1e-9);
        assert_eq!(result.matrix_rows[1].isoform_id, "iso2");
        assert!((result.matrix_rows[1].counts[0] - 0.5).abs() < 1e-9);
        assert!((result.matrix_rows[1].counts[1] - 1.0).abs() < 1e-9);

        assert_eq!(result.long_rows.len(), 4);
        let s1_total: f64 = result
            .long_rows
            .iter()
            .filter(|row| row.sample == "S1")
            .map(|row| row.proportion)
            .sum();
        let s2_total: f64 = result
            .long_rows
            .iter()
            .filter(|row| row.sample == "S2")
            .map(|row| row.proportion)
            .sum();
        assert!((s1_total - 1.0).abs() < 1e-9);
        assert!((s2_total - 1.0).abs() < 1e-9);

        assert_eq!(result.group_rows.len(), 4);
    }

    #[test]
    fn mapping_counts_match_subread_counts() {
        let references = vec![make_tx("ref_a", "ref_a", "GENEA")];
        let isoforms = vec![
            make_tx("iso1", "S1::r1,S1::r2,S2::r3,|0", "GENEA"),
            make_tx("iso2", "S1::r2,S2::r4,|0", "GENEA"),
        ];
        let mapping = vec![
            ("S1::r1".to_owned(), "iso1".to_owned()),
            ("S1::r2".to_owned(), "iso1".to_owned()),
            ("S1::r2".to_owned(), "iso2".to_owned()),
            ("S2::r3".to_owned(), "iso1".to_owned()),
            ("S2::r4".to_owned(), "iso2".to_owned()),
        ];

        let by_subreads = count_multi_by_subreads(&isoforms, &references, &sample_rows()).unwrap();
        let by_mapping =
            count_multi_by_read_to_isoform(&isoforms, &mapping, &sample_rows()).unwrap();

        assert_eq!(by_subreads.matrix_rows.len(), by_mapping.matrix_rows.len());
        for (left, right) in by_subreads
            .matrix_rows
            .iter()
            .zip(by_mapping.matrix_rows.iter())
        {
            assert_eq!(left.gene, right.gene);
            assert_eq!(left.isoform_id, right.isoform_id);
            assert_eq!(left.counts.len(), right.counts.len());
            for (l, r) in left.counts.iter().zip(right.counts.iter()) {
                assert!((l - r).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn duplicate_mapping_rows_do_not_change_multi_sample_molecule_counts() {
        let isoforms = vec![
            make_tx("iso1", "none", "GENEA"),
            make_tx("iso2", "none", "GENEA"),
        ];
        let mapping = vec![
            ("S1::molecule".to_owned(), "iso1".to_owned()),
            ("S1::molecule".to_owned(), "iso1".to_owned()),
            ("S1::molecule".to_owned(), "iso2".to_owned()),
        ];

        let result = count_multi_by_read_to_isoform(&isoforms, &mapping, &sample_rows()).unwrap();
        assert_eq!(result.matrix_rows[0].counts, vec![0.5, 0.0]);
        assert_eq!(result.matrix_rows[1].counts, vec![0.5, 0.0]);
    }

    #[test]
    fn embedded_read_matching_reference_id_is_still_counted() {
        let references = vec![make_tx("S1::same", "none", "GENEA")];
        let isoforms = vec![make_tx("iso1", "S1::same,|0", "GENEA")];

        let result = count_multi_by_subreads(&isoforms, &references, &sample_rows()).unwrap();
        assert_eq!(result.matrix_rows[0].counts, vec![1.0, 0.0]);
    }

    #[test]
    fn multi_sample_count_rejects_duplicate_catalog_ids() {
        let isoforms = vec![
            make_tx_with_exons("duplicate", "none", "GENEA", &[(100, 200)]),
            make_tx_with_exons("duplicate", "none", "GENEB", &[(300, 400)]),
        ];
        let error = count_multi_by_read_to_isoform(&isoforms, &[], &sample_rows()).unwrap_err();
        assert!(error.to_string().contains("duplicate isoform id"));
    }

    #[test]
    fn unique_mapping_counts_each_read_once() {
        let isoforms = vec![
            make_tx_with_exons("long_ref", "S1::r1,|0", "GENEA", &[(50, 60), (100, 200)]),
            make_tx_with_exons("closest_novel", "S1::r1,|0", "GENEA", &[(100, 200)]),
        ];
        let read = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "S1::r1".to_owned(),
            vec![Interval::new(Coord::new(100), Coord::new(200)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(100),
                thick_end: Coord::new(200),
                item_rgb: "0".to_owned(),
                extra_fields: vec![],
            },
        )
        .unwrap();
        let mapping = vec![
            ("S1::r1".to_owned(), "long_ref".to_owned()),
            ("S1::r1".to_owned(), "closest_novel".to_owned()),
        ];

        let unique =
            crate::count::select_unique_best_read_to_isoform(&[read], &isoforms, &mapping).unwrap();
        let result = count_multi_by_read_to_isoform(&isoforms, &unique, &sample_rows()).unwrap();
        let long_ref = result
            .matrix_rows
            .iter()
            .find(|row| row.isoform_id == "long_ref")
            .unwrap();
        let closest = result
            .matrix_rows
            .iter()
            .find(|row| row.isoform_id == "closest_novel")
            .unwrap();

        assert_eq!(long_ref.counts, vec![0.0, 0.0]);
        assert_eq!(closest.counts, vec![1.0, 0.0]);
    }

    #[test]
    fn errors_when_read_name_has_no_sample_prefix() {
        let references = vec![make_tx("ref_a", "ref_a", "GENEA")];
        let isoforms = vec![make_tx("iso1", "r1,|0", "GENEA")];

        let err = count_multi_by_subreads(&isoforms, &references, &sample_rows())
            .unwrap_err()
            .to_string();
        assert!(err.contains("missing sample prefix"));
    }

    #[test]
    fn errors_on_unknown_sample_prefix() {
        let references = vec![make_tx("ref_a", "ref_a", "GENEA")];
        let isoforms = vec![make_tx("iso1", "S3::r1,|0", "GENEA")];

        let err = count_multi_by_subreads(&isoforms, &references, &sample_rows())
            .unwrap_err()
            .to_string();
        assert!(err.contains("unknown sample"));
    }
}
