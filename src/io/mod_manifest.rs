//! Modification sample manifest parsing.

use std::collections::HashSet;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::sample::SAMPLE_DELIM;

/// One validated `(sample, assay)` modification input row.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModSampleRow {
    /// Sample identifier matching the TrackCluster sample manifest.
    pub sample: String,
    /// Assay compatibility stratum matching the assay metadata.
    pub assay_id: String,
    /// Normalized observation TSV path.
    pub observations: PathBuf,
    /// Assay provenance JSON path.
    pub assay_metadata: PathBuf,
    /// Optional genome-aligned BAM used only for exact coverage.
    pub coverage_bam: Option<PathBuf>,
}

const COLUMNS: [&str; 5] = [
    "sample",
    "assay_id",
    "observations",
    "assay_metadata",
    "coverage_bam",
];

fn resolve_path(manifest: &Path, value: &str) -> PathBuf {
    let value = PathBuf::from(value);
    if value.is_absolute() {
        value
    } else {
        manifest
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(value)
    }
}

fn require_input(path: PathBuf, role: &str, line: usize) -> anyhow::Result<PathBuf> {
    let metadata = std::fs::symlink_metadata(&path)
        .with_context(|| format!("mod manifest line {line}: inspect {role} {path:?}"))?;
    if metadata.file_type().is_symlink() || !metadata.is_file() {
        anyhow::bail!("mod manifest line {line}: {role} is not a regular file: {path:?}");
    }
    Ok(path)
}

/// Read a strict five-column modification sample manifest.
pub fn read_mod_manifest_tsv(path: &Path) -> anyhow::Result<Vec<ModSampleRow>> {
    let file = std::fs::File::open(path).with_context(|| format!("open mod manifest {path:?}"))?;
    let reader = std::io::BufReader::new(file);
    let mut header_seen = false;
    let mut rows = Vec::new();
    let mut seen = HashSet::new();

    for (line_index, result) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = result.with_context(|| format!("read mod manifest {path:?}:{line_number}"))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if !header_seen {
            if fields != COLUMNS {
                anyhow::bail!(
                    "mod manifest {path:?}:{line_number} header mismatch: expected {:?}, found {:?}",
                    COLUMNS,
                    fields
                );
            }
            header_seen = true;
            continue;
        }
        if fields.len() != COLUMNS.len() {
            anyhow::bail!(
                "mod manifest {path:?}:{line_number} must have exactly {} fields, found {}",
                COLUMNS.len(),
                fields.len()
            );
        }
        let sample = fields[0];
        let assay_id = fields[1];
        for (name, value) in [("sample", sample), ("assay_id", assay_id)] {
            if value.trim().is_empty() || value == "NA" || value.chars().any(char::is_control) {
                anyhow::bail!("mod manifest {path:?}:{line_number} has invalid {name} {value:?}");
            }
        }
        if sample.contains(SAMPLE_DELIM) {
            anyhow::bail!(
                "mod manifest {path:?}:{line_number} sample {sample:?} must not contain {SAMPLE_DELIM:?}"
            );
        }
        if !seen.insert((sample.to_owned(), assay_id.to_owned())) {
            anyhow::bail!(
                "mod manifest {path:?}:{line_number} duplicates sample/assay pair ({sample:?}, {assay_id:?})"
            );
        }
        if fields[2].is_empty() || fields[2] == "NA" {
            anyhow::bail!("mod manifest {path:?}:{line_number} observations path is missing");
        }
        if fields[3].is_empty() || fields[3] == "NA" {
            anyhow::bail!("mod manifest {path:?}:{line_number} assay_metadata path is missing");
        }
        let observations =
            require_input(resolve_path(path, fields[2]), "observations", line_number)?;
        let assay_metadata =
            require_input(resolve_path(path, fields[3]), "assay metadata", line_number)?;
        let coverage_bam = match fields[4] {
            "NA" => None,
            "" => {
                anyhow::bail!(
                    "mod manifest {path:?}:{line_number} coverage_bam must be a path or NA"
                );
            }
            value => Some(require_input(
                resolve_path(path, value),
                "coverage BAM",
                line_number,
            )?),
        };
        rows.push(ModSampleRow {
            sample: sample.to_owned(),
            assay_id: assay_id.to_owned(),
            observations,
            assay_metadata,
            coverage_bam,
        });
    }

    if !header_seen {
        anyhow::bail!("mod manifest {path:?} is missing its header");
    }
    if rows.is_empty() {
        anyhow::bail!("mod manifest {path:?} has no data rows");
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn temp_dir(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-mod-manifest-{label}-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn resolves_paths_and_rejects_duplicate_sample_assays() {
        let dir = temp_dir("valid");
        fs::write(dir.join("obs.tsv"), "header\n").unwrap();
        fs::write(dir.join("assay.json"), "{}\n").unwrap();
        let manifest = dir.join("mods.tsv");
        fs::write(
            &manifest,
            concat!(
                "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n",
                "S1\ta1\tobs.tsv\tassay.json\tNA\n",
            ),
        )
        .unwrap();
        let rows = read_mod_manifest_tsv(&manifest).unwrap();
        assert_eq!(rows.len(), 1);
        assert!(rows[0].observations.is_absolute());
        assert_eq!(rows[0].coverage_bam, None);

        fs::write(
            &manifest,
            concat!(
                "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n",
                "S1\ta1\tobs.tsv\tassay.json\tNA\n",
                "S1\ta1\tobs.tsv\tassay.json\tNA\n",
            ),
        )
        .unwrap();
        assert!(read_mod_manifest_tsv(&manifest)
            .unwrap_err()
            .to_string()
            .contains("duplicates sample/assay pair"));
        let _ = fs::remove_dir_all(dir);
    }
}
