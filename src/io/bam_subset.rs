//! Deterministic BAM splitting by query-name membership.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use anyhow::Context;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;

/// One requested BAM subset.
#[derive(Clone, Copy, Debug)]
pub struct BamSubsetTarget<'a> {
    /// Destination BAM path.
    pub path: &'a Path,
    /// Raw BAM query names to retain.
    pub read_names: &'a BTreeSet<String>,
}

/// Audit counts for one written BAM subset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BamSubsetQc {
    /// Destination BAM path.
    pub path: PathBuf,
    /// BAM records written, including selected secondary and supplementary records.
    pub written_records: usize,
    /// Distinct mapped primary query names written.
    pub primary_reads: usize,
}

/// Audit counts for one source BAM split.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BamSplitResult {
    /// All decoded source BAM records.
    pub input_records: usize,
    /// Per-destination counts in the same order as the requested targets.
    pub targets: Vec<BamSubsetQc>,
}

struct TargetWriter<'a> {
    path: PathBuf,
    selected: &'a BTreeSet<String>,
    writer: bam::io::Writer<noodles_bgzf::io::Writer<BufWriter<File>>>,
    seen_primary: HashSet<String>,
    written_records: usize,
}

/// Split one BAM into query-name subsets in a single source scan.
///
/// Records retain their original order and auxiliary tags. A selected query
/// name must have exactly one mapped primary alignment; secondary and
/// supplementary records for that query name are copied as well.
pub fn split_bam_by_query_name(
    input: &Path,
    targets: &[BamSubsetTarget<'_>],
) -> anyhow::Result<BamSplitResult> {
    if targets.is_empty() {
        anyhow::bail!("BAM split requires at least one output target");
    }

    let mut seen_paths = HashSet::new();
    let mut membership: HashMap<&str, Vec<usize>> = HashMap::new();
    for (target_index, target) in targets.iter().enumerate() {
        if target.read_names.is_empty() {
            anyhow::bail!("BAM subset target {:?} has no selected reads", target.path);
        }
        let path = target.path.to_path_buf();
        if !seen_paths.insert(path.clone()) {
            anyhow::bail!("BAM subset output path is repeated: {path:?}");
        }
        for read_name in target.read_names {
            membership
                .entry(read_name.as_str())
                .or_default()
                .push(target_index);
        }
    }

    let file = File::open(input).with_context(|| format!("open source BAM {input:?}"))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .with_context(|| format!("read source BAM header {input:?}"))?;

    let mut writers = Vec::with_capacity(targets.len());
    for target in targets {
        let file = File::create(target.path)
            .with_context(|| format!("create BAM subset {:?}", target.path))?;
        let mut writer = bam::io::Writer::new(BufWriter::new(file));
        writer
            .write_header(&header)
            .with_context(|| format!("write BAM subset header {:?}", target.path))?;
        writers.push(TargetWriter {
            path: target.path.to_path_buf(),
            selected: target.read_names,
            writer,
            seen_primary: HashSet::new(),
            written_records: 0,
        });
    }

    let mut input_records = 0usize;
    for (record_index, record) in reader.records().enumerate() {
        let ordinal = record_index + 1;
        let record = record.with_context(|| format!("read BAM record {ordinal} from {input:?}"))?;
        input_records += 1;
        let Some(read_name) = record.name().map(ToString::to_string) else {
            continue;
        };
        let Some(target_indices) = membership.get(read_name.as_str()) else {
            continue;
        };
        let flags = record.flags();
        let mapped_primary =
            !flags.is_unmapped() && !flags.is_secondary() && !flags.is_supplementary();
        for &target_index in target_indices {
            let target = &mut writers[target_index];
            if mapped_primary && !target.seen_primary.insert(read_name.clone()) {
                anyhow::bail!(
                    "source BAM {input:?} contains duplicate mapped primary query name \
                     {read_name:?} selected for {:?}",
                    target.path
                );
            }
            target
                .writer
                .write_alignment_record(&header, &record)
                .with_context(|| {
                    format!(
                        "write selected BAM record {ordinal} ({read_name:?}) to {:?}",
                        target.path
                    )
                })?;
            target.written_records += 1;
        }
    }

    let mut qc = Vec::with_capacity(writers.len());
    for mut target in writers {
        if target.seen_primary.len() != target.selected.len() {
            let first_missing = target
                .selected
                .iter()
                .find(|name| !target.seen_primary.contains(name.as_str()))
                .expect("different set lengths imply a missing name");
            anyhow::bail!(
                "source BAM {input:?} is missing a mapped primary alignment for {} of {} \
                 selected reads in {:?}; first missing query name is {first_missing:?}",
                target.selected.len() - target.seen_primary.len(),
                target.selected.len(),
                target.path
            );
        }
        target
            .writer
            .try_finish()
            .with_context(|| format!("finish BAM subset {:?}", target.path))?;
        qc.push(BamSubsetQc {
            path: target.path,
            written_records: target.written_records,
            primary_reads: target.seen_primary.len(),
        });
    }

    Ok(BamSplitResult {
        input_records,
        targets: qc,
    })
}
