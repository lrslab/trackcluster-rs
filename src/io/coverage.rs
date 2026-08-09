//! Exact genomic-base coverage from primary BAM alignments.

use std::collections::HashSet;
use std::fs::File;
use std::path::Path;

use anyhow::Context;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;

use crate::model::{Coord, Interval, Strand};

/// Match-bearing reference blocks for one primary alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadCoverage {
    /// TrackCluster molecule ID in `<sample>::<source_read_id>` form.
    pub read_id: String,
    /// Reference sequence name.
    pub chrom: String,
    /// Alignment strand.
    pub strand: Strand,
    /// Half-open blocks covered by `M`, `=`, or `X` CIGAR operations.
    pub match_blocks: Vec<Interval>,
}

/// Deterministic primary-alignment coverage scan.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BamCoverageResult {
    /// Retained primary mapped alignments.
    pub reads: Vec<ReadCoverage>,
    /// All decoded BAM records.
    pub total_records: usize,
    /// Unmapped records omitted.
    pub skipped_unmapped: usize,
    /// Secondary records omitted.
    pub skipped_secondary: usize,
    /// Supplementary records omitted.
    pub skipped_supplementary: usize,
}

fn checked_advance(position: u32, length: usize) -> anyhow::Result<u32> {
    let length = u32::try_from(length).context("CIGAR operation length exceeds u32")?;
    position
        .checked_add(length)
        .context("CIGAR reference span exceeds u32")
}

fn push_match_block(blocks: &mut Vec<Interval>, start: u32, end: u32) -> anyhow::Result<()> {
    if let Some(previous) = blocks.last_mut() {
        if previous.end.get() == start {
            previous.end = Coord::new(end);
            return Ok(());
        }
    }
    blocks.push(Interval::new(Coord::new(start), Coord::new(end))?);
    Ok(())
}

fn normalized_read_id(sample: &str, source_name: &str) -> anyhow::Result<String> {
    if source_name.trim().is_empty()
        || source_name == "*"
        || source_name.chars().any(char::is_control)
    {
        anyhow::bail!("coverage BAM query name is empty, '*', or contains a control character");
    }
    if source_name.contains(crate::sample::SAMPLE_DELIM) {
        let (source_sample, _) =
            crate::sample::split_tagged_read_name(source_name).with_context(|| {
                format!("coverage BAM query name {source_name:?} has a malformed sample prefix")
            })?;
        if source_sample != sample {
            anyhow::bail!(
                "coverage BAM query name {source_name:?} has sample prefix {source_sample:?}, expected {sample:?}"
            );
        }
        Ok(source_name.to_owned())
    } else {
        Ok(crate::sample::tagged_read_name(sample, source_name))
    }
}

/// Read mapped primary BAM records and retain exact base-covering CIGAR blocks.
///
/// Read names are tagged exactly as TrackCluster manifest reads are tagged.
/// Duplicate primary query names are rejected because molecule coverage would
/// otherwise be counted more than once.
pub fn read_primary_bam_coverage(path: &Path, sample: &str) -> anyhow::Result<BamCoverageResult> {
    if sample.trim().is_empty()
        || sample.chars().any(char::is_control)
        || sample.contains(crate::sample::SAMPLE_DELIM)
    {
        anyhow::bail!(
            "coverage BAM sample must not be empty, contain control characters, or contain {:?}",
            crate::sample::SAMPLE_DELIM
        );
    }
    let file = File::open(path).with_context(|| format!("open coverage BAM {path:?}"))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .with_context(|| format!("read coverage BAM header {path:?}"))?;
    let mut result = BamCoverageResult::default();
    let mut seen_names = HashSet::new();

    for (record_index, record) in reader.records().enumerate() {
        let ordinal = record_index + 1;
        let record =
            record.with_context(|| format!("read coverage BAM record {ordinal} from {path:?}"))?;
        result.total_records += 1;
        let flags = record.flags();
        if flags.is_unmapped() {
            result.skipped_unmapped += 1;
            continue;
        }
        if flags.is_secondary() {
            result.skipped_secondary += 1;
            continue;
        }
        if flags.is_supplementary() {
            result.skipped_supplementary += 1;
            continue;
        }

        let source_name = record
            .name()
            .map(ToString::to_string)
            .filter(|name| !name.trim().is_empty() && name != "*")
            .with_context(|| {
                format!("mapped primary coverage BAM record {ordinal} has no query name")
            })?;
        let read_id = normalized_read_id(sample, &source_name).with_context(|| {
            format!("validate coverage BAM record {ordinal} query name {source_name:?}")
        })?;
        if !seen_names.insert(read_id.clone()) {
            anyhow::bail!(
                "coverage BAM {path:?} contains duplicate primary query name {source_name:?}"
            );
        }

        let reference_id = record
            .reference_sequence_id()
            .transpose()
            .with_context(|| {
                format!("decode reference ID for coverage BAM record {ordinal} ({source_name:?})")
            })?
            .with_context(|| {
                format!(
                    "mapped coverage BAM record {ordinal} ({source_name:?}) has no reference sequence"
                )
            })?;
        let (reference_name, reference_sequence) = header
            .reference_sequences()
            .get_index(reference_id)
            .with_context(|| {
                format!(
                    "coverage BAM record {ordinal} references missing header sequence ID {reference_id}"
                )
            })?;
        let chrom = reference_name.to_string();
        let reference_length = u32::try_from(usize::from(reference_sequence.length()))
            .context("coverage BAM reference length exceeds u32")?;
        let alignment_start = record
            .alignment_start()
            .transpose()
            .with_context(|| {
                format!(
                    "decode alignment start for coverage BAM record {ordinal} ({source_name:?})"
                )
            })?
            .with_context(|| {
                format!(
                    "mapped coverage BAM record {ordinal} ({source_name:?}) has no alignment start"
                )
            })?;
        let mut reference_pos0 = u32::try_from(usize::from(alignment_start) - 1)
            .context("coverage BAM alignment start exceeds u32")?;
        let mut match_blocks = Vec::new();
        for operation in record.cigar().iter() {
            let operation = operation.with_context(|| {
                format!("parse CIGAR for coverage BAM record {ordinal} ({source_name:?})")
            })?;
            let start = reference_pos0;
            if operation.kind().consumes_reference() {
                reference_pos0 = checked_advance(reference_pos0, operation.len())?;
            }
            if matches!(
                operation.kind(),
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
            ) {
                push_match_block(&mut match_blocks, start, reference_pos0)?;
            }
        }
        if match_blocks.is_empty() {
            anyhow::bail!(
                "mapped primary coverage BAM record {ordinal} ({source_name:?}) has no M, =, or X CIGAR bases"
            );
        }
        if reference_pos0 > reference_length {
            anyhow::bail!(
                "coverage BAM record {ordinal} ({source_name:?}) ends at {reference_pos0} beyond reference {chrom:?} length {reference_length}"
            );
        }
        result.reads.push(ReadCoverage {
            read_id,
            chrom,
            strand: if flags.is_reverse_complemented() {
                Strand::Minus
            } else {
                Strand::Plus
            },
            match_blocks,
        });
    }
    result
        .reads
        .sort_by(|left, right| left.read_id.cmp(&right.read_id));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::num::NonZero;
    use std::time::{SystemTime, UNIX_EPOCH};

    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record_buf::{Cigar, Sequence};
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};

    use super::*;

    fn temp_bam() -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "trackcluster-coverage-{}-{nonce}.bam",
            std::process::id()
        ))
    }

    #[test]
    fn exact_coverage_excludes_deletions_splices_and_query_only_operations() {
        let path = temp_bam();
        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
            )
            .build();
        let cigar: Cigar = [
            Op::new(Kind::SoftClip, 2),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::SequenceMatch, 2),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::SequenceMismatch, 2),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 1),
        ]
        .into_iter()
        .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_name("read1")
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start("101".parse().unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(b"AAAAAAAAAAA".to_vec()))
            .build();
        let mut writer = bam::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_header(&header).unwrap();
        writer.write_alignment_record(&header, &record).unwrap();
        writer.try_finish().unwrap();

        let result = read_primary_bam_coverage(&path, "S1").unwrap();
        assert_eq!(result.total_records, 1);
        assert_eq!(result.reads.len(), 1);
        assert_eq!(result.reads[0].read_id, "S1::read1");
        assert_eq!(
            result.reads[0]
                .match_blocks
                .iter()
                .map(|block| (block.start.get(), block.end.get()))
                .collect::<Vec<_>>(),
            vec![(100, 105), (107, 109), (114, 115)]
        );
        let _ = fs::remove_file(path);
    }

    #[test]
    fn coverage_read_ids_preserve_a_matching_existing_sample_prefix() {
        assert_eq!(normalized_read_id("S1", "read1").unwrap(), "S1::read1");
        assert_eq!(normalized_read_id("S1", "S1::read1").unwrap(), "S1::read1");
        assert!(normalized_read_id("S1", "S2::read1").is_err());
        assert!(normalized_read_id("S1", "S1::").is_err());
    }
}
