use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::count::parse_subreads;
use crate::io::manifest::{read_manifest_tsv, SampleRow};
use crate::model::Transcript;
use crate::sample::split_tagged_read_name;

const GENE_NAME_COL: usize = 5;
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
    pub long_tsv: PathBuf,
    pub matrix_tsv: PathBuf,
    pub group_tsv: Option<PathBuf>,
}

fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

fn isoform_gene(tx: &Transcript) -> String {
    let gene = tx
        .extra_fields
        .get(GENE_NAME_COL)
        .map(|value| value.trim())
        .unwrap_or("");
    if gene.is_empty() || gene == "none" {
        return "none".to_owned();
    }
    gene.to_owned()
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

pub fn count_multi_by_subreads(
    isoforms: &[Transcript],
    references: &[Transcript],
    samples: &[SampleRow],
) -> anyhow::Result<MultiSampleCountResult> {
    if samples.is_empty() {
        anyhow::bail!("count-multi requires at least one sample");
    }

    let ref_names: HashSet<&str> = references.iter().map(|tx| tx.name.as_str()).collect();
    let sample_to_idx: HashMap<&str, usize> = samples
        .iter()
        .enumerate()
        .map(|(idx, row)| (row.sample.as_str(), idx))
        .collect();

    let mut read_occurrence: HashMap<String, u32> = HashMap::new();
    let mut read_to_sample_idx: HashMap<String, usize> = HashMap::new();
    for isoform in isoforms {
        for read_name in parse_subreads(isoform) {
            if ref_names.contains(read_name) {
                continue;
            }

            let (sample_name, _raw_read_id) = split_tagged_read_name(read_name).ok_or_else(|| {
                anyhow::anyhow!(
                    "read id {read_name:?} in isoform {:?} is missing sample prefix; expected format '<sample>::<read_id>'",
                    isoform.name
                )
            })?;
            let sample_idx = sample_to_idx.get(sample_name).copied().ok_or_else(|| {
                anyhow::anyhow!(
                    "read id {read_name:?} in isoform {:?} references unknown sample {sample_name:?} (not present in manifest)",
                    isoform.name
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
                    read_to_sample_idx.insert(read_name.to_owned(), sample_idx);
                }
            }

            *read_occurrence.entry(read_name.to_owned()).or_insert(0) += 1;
        }
    }

    let mut matrix_rows: Vec<UsageMatrixRow> = isoforms
        .iter()
        .map(|isoform| {
            let mut counts = vec![0.0f64; samples.len()];
            for read_name in parse_subreads(isoform) {
                if ref_names.contains(read_name) {
                    continue;
                }
                let sample_idx = *read_to_sample_idx.get(read_name).ok_or_else(|| {
                    anyhow::anyhow!(
                        "internal error: sample lookup missing for read id {read_name:?}"
                    )
                })?;
                let denom = read_occurrence.get(read_name).copied().unwrap_or(0);
                if denom > 0 {
                    counts[sample_idx] += 1.0f64 / (denom as f64);
                }
            }

            Ok(UsageMatrixRow {
                gene: isoform_gene(isoform),
                isoform_id: isoform.name.clone(),
                counts,
            })
        })
        .collect::<Result<Vec<_>, anyhow::Error>>()?;
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

    Ok(MultiSampleCountResult {
        matrix_rows,
        long_rows,
        group_rows,
    })
}

pub fn write_usage_long_tsv(
    path: &Path,
    rows: &[UsageLongRow],
    include_group: bool,
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    if include_group {
        writeln!(
            &mut writer,
            "gene\tisoform_id\tsample\tgroup\tcount\tproportion\tgene_total"
        )?;
        for row in rows {
            writeln!(
                &mut writer,
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
            &mut writer,
            "gene\tisoform_id\tsample\tcount\tproportion\tgene_total"
        )?;
        for row in rows {
            writeln!(
                &mut writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                row.gene, row.isoform_id, row.sample, row.count, row.proportion, row.gene_total
            )?;
        }
    }
    Ok(())
}

pub fn write_counts_matrix_tsv(
    path: &Path,
    rows: &[UsageMatrixRow],
    samples: &[SampleRow],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    write!(&mut writer, "gene\tisoform_id")?;
    for sample in samples {
        write!(&mut writer, "\t{}", sample.sample)?;
    }
    writeln!(&mut writer)?;

    for row in rows {
        write!(&mut writer, "{}\t{}", row.gene, row.isoform_id)?;
        for count in &row.counts {
            write!(&mut writer, "\t{}", count)?;
        }
        writeln!(&mut writer)?;
    }
    Ok(())
}

pub fn write_group_usage_tsv(path: &Path, rows: &[GroupUsageRow]) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        &mut writer,
        "gene\tisoform_id\tgroup\tcount\tproportion\tgene_total"
    )?;
    for row in rows {
        writeln!(
            &mut writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.gene, row.isoform_id, row.group, row.count, row.proportion, row.gene_total
        )?;
    }
    Ok(())
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

    if let Some(parent) = out_prefix
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).with_context(|| format!("create output dir {parent:?}"))?;
    }

    let result = count_multi_by_subreads(&isoforms, &refs, &sample_rows)?;
    let include_group = sample_rows.iter().any(|sample| sample.group.is_some());

    let long_tsv = append_suffix(out_prefix, ".isoform_usage.long.tsv");
    let matrix_tsv = append_suffix(out_prefix, ".isoform_counts.matrix.tsv");
    write_usage_long_tsv(&long_tsv, &result.long_rows, include_group)
        .with_context(|| format!("write long output {long_tsv:?}"))?;
    write_counts_matrix_tsv(&matrix_tsv, &result.matrix_rows, &sample_rows)
        .with_context(|| format!("write matrix output {matrix_tsv:?}"))?;

    let group_tsv = if result.group_rows.is_empty() {
        None
    } else {
        let path = append_suffix(out_prefix, ".isoform_usage.group.tsv");
        write_group_usage_tsv(&path, &result.group_rows)
            .with_context(|| format!("write group output {path:?}"))?;
        Some(path)
    };

    Ok(MultiSampleOutputPaths {
        long_tsv,
        matrix_tsv,
        group_tsv,
    })
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(name: &str, name2: &str, gene: &str) -> Transcript {
        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            name.to_owned(),
            vec![Interval::new(Coord::new(100), Coord::new(200)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(100),
                thick_end: Coord::new(200),
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
