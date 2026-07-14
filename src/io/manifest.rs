use std::collections::HashSet;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::sample::SAMPLE_DELIM;

#[derive(Clone, Debug, PartialEq, Eq)]
/// One validated row from a multi-sample reads manifest.
pub struct SampleRow {
    /// Unique sample identifier.
    pub sample: String,
    /// Optional experimental group.
    pub group: Option<String>,
    /// Reads BED path, resolved relative to the manifest.
    pub reads: PathBuf,
}

fn split_tsv_line(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn column_index(columns: &[&str], name: &str) -> Option<usize> {
    columns
        .iter()
        .position(|column| column.trim().eq_ignore_ascii_case(name))
}

fn resolve_reads_path(manifest_path: &Path, reads_field: &str) -> PathBuf {
    let reads_path = PathBuf::from(reads_field);
    if reads_path.is_absolute() {
        return reads_path;
    }
    manifest_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(reads_path)
}

/// Read a TSV manifest with required `sample` and `reads` columns.
pub fn read_manifest_tsv(path: &Path) -> anyhow::Result<Vec<SampleRow>> {
    let file = std::fs::File::open(path).with_context(|| format!("open manifest {path:?}"))?;
    let reader = std::io::BufReader::new(file);

    let mut header: Option<Vec<String>> = None;
    let mut sample_rows: Vec<SampleRow> = Vec::new();
    let mut seen_samples: HashSet<String> = HashSet::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line_no = line_no + 1;
        let line = line.with_context(|| format!("read manifest {path:?}:{line_no}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if header.is_none() {
            let columns = split_tsv_line(trimmed)
                .into_iter()
                .map(|column| column.trim().to_owned())
                .collect::<Vec<_>>();

            if !columns.iter().any(|column| !column.is_empty()) {
                continue;
            }

            let refs = columns.iter().map(String::as_str).collect::<Vec<_>>();
            if column_index(&refs, "sample").is_none() || column_index(&refs, "reads").is_none() {
                anyhow::bail!(
                    "manifest {path:?}:{line_no}: header must include 'sample' and 'reads' columns"
                );
            }
            header = Some(columns);
            continue;
        }

        let columns = header
            .as_ref()
            .expect("header checked in previous branch")
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let values = split_tsv_line(trimmed);

        let sample_idx = column_index(&columns, "sample").expect("validated at header parse");
        let reads_idx = column_index(&columns, "reads").expect("validated at header parse");
        let group_idx = column_index(&columns, "group");

        let sample = values
            .get(sample_idx)
            .map(|value| value.trim())
            .unwrap_or_default();
        if sample.is_empty() {
            anyhow::bail!("manifest {path:?}:{line_no}: sample is empty");
        }
        if sample.contains(SAMPLE_DELIM) {
            anyhow::bail!(
                "manifest {path:?}:{line_no}: sample {sample:?} must not contain {SAMPLE_DELIM:?}"
            );
        }
        if !seen_samples.insert(sample.to_owned()) {
            anyhow::bail!("manifest {path:?}:{line_no}: duplicate sample {sample:?}");
        }

        let reads_field = values
            .get(reads_idx)
            .map(|value| value.trim())
            .unwrap_or_default();
        if reads_field.is_empty() {
            anyhow::bail!("manifest {path:?}:{line_no}: reads path is empty");
        }
        let reads = resolve_reads_path(path, reads_field);
        if !reads.exists() {
            anyhow::bail!(
                "manifest {path:?}:{line_no}: reads file not found for sample {sample:?}: {reads:?}"
            );
        }

        let group = group_idx
            .and_then(|idx| values.get(idx))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(str::to_owned);

        sample_rows.push(SampleRow {
            sample: sample.to_owned(),
            group,
            reads,
        });
    }

    if header.is_none() {
        anyhow::bail!("manifest {path:?}: missing header row");
    }
    if sample_rows.is_empty() {
        anyhow::bail!("manifest {path:?}: no sample rows found");
    }

    Ok(sample_rows)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn fresh_temp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "trackcluster_rs_manifest_{}_{}_{}",
            prefix,
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn parses_manifest_with_optional_group() {
        let dir = fresh_temp_dir("ok");
        let reads1 = dir.join("S1.reads.bed");
        let reads2 = dir.join("S2.reads.bed");
        fs::write(&reads1, "chr1\t0\t10\tr1\t0\t+\t0\t10\t0\t1\t10,\t0,\n").unwrap();
        fs::write(&reads2, "chr1\t0\t10\tr2\t0\t+\t0\t10\t0\t1\t10,\t0,\n").unwrap();

        let manifest = dir.join("samples.tsv");
        fs::write(
            &manifest,
            "sample\tgroup\treads\nS1\tcontrol\tS1.reads.bed\nS2\ttreated\tS2.reads.bed\n",
        )
        .unwrap();

        let rows = read_manifest_tsv(&manifest).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].sample, "S1");
        assert_eq!(rows[0].group.as_deref(), Some("control"));
        assert!(rows[0].reads.is_absolute());
        assert_eq!(rows[1].sample, "S2");
        assert_eq!(rows[1].group.as_deref(), Some("treated"));
    }

    #[test]
    fn errors_on_missing_required_columns() {
        let dir = fresh_temp_dir("missing_cols");
        let manifest = dir.join("samples.tsv");
        fs::write(&manifest, "sample\tgroup\nS1\tcontrol\n").unwrap();

        let err = read_manifest_tsv(&manifest).unwrap_err().to_string();
        assert!(err.contains("header must include 'sample' and 'reads' columns"));
    }

    #[test]
    fn errors_on_duplicate_samples() {
        let dir = fresh_temp_dir("dup");
        let reads = dir.join("reads.bed");
        fs::write(&reads, "chr1\t0\t10\tr1\t0\t+\t0\t10\t0\t1\t10,\t0,\n").unwrap();
        let manifest = dir.join("samples.tsv");
        fs::write(&manifest, "sample\treads\nS1\treads.bed\nS1\treads.bed\n").unwrap();

        let err = read_manifest_tsv(&manifest).unwrap_err().to_string();
        assert!(err.contains("duplicate sample"));
    }
}
