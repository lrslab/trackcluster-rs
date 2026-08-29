//! Pure-Rust conversion of spliced BAM alignments into TrackCluster transcripts.

use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use anyhow::Context;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

/// Policy for a decoded BAM record that cannot be converted to a transcript.
///
/// BAM header, framing, decompression, and record-decode errors always fail:
/// after those errors, continuing at the next record boundary is not guaranteed
/// to be safe.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum InvalidRecordPolicy {
    /// Exclude only the invalid decoded record and count a stable reason.
    #[default]
    Skip,
    /// Stop at the first decoded record that cannot be converted.
    Fail,
}

/// Stable reason for excluding a decoded BAM record under skip policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidRecordReason {
    /// A mapped record has no usable query name.
    MissingQueryName,
    /// A query name is not valid UTF-8 or contains a control character.
    InvalidQueryName,
    /// A mapped record has no valid entry in the BAM reference dictionary.
    InvalidReference,
    /// A mapped record has no valid BED-representable alignment start.
    InvalidAlignmentStart,
    /// A CIGAR operation could not be decoded.
    InvalidCigar,
    /// A decoded CIGAR cannot form nonempty aligned exon blocks.
    InvalidCigarStructure,
    /// CIGAR coordinate arithmetic exceeded the BED coordinate range.
    CoordinateOverflow,
    /// The alignment extends beyond the reference length declared in the header.
    ReferenceOutOfBounds,
    /// The converted fields violate TrackCluster transcript geometry.
    InvalidTranscript,
}

impl InvalidRecordReason {
    /// Stable token used in command diagnostics.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::MissingQueryName => "missing_query_name",
            Self::InvalidQueryName => "invalid_query_name",
            Self::InvalidReference => "invalid_reference",
            Self::InvalidAlignmentStart => "invalid_alignment_start",
            Self::InvalidCigar => "invalid_cigar",
            Self::InvalidCigarStructure => "invalid_cigar_structure",
            Self::CoordinateOverflow => "coordinate_overflow",
            Self::ReferenceOutOfBounds => "reference_out_of_bounds",
            Self::InvalidTranscript => "invalid_transcript",
        }
    }
}

/// Aggregated invalid-record count with one bounded diagnostic example.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InvalidRecordStats {
    /// Records excluded for this reason.
    pub count: usize,
    /// One-based ordinal of the first excluded record with this reason.
    pub first_record_ordinal: usize,
    /// Diagnostic for the first excluded record with this reason.
    pub first_error: String,
}

#[derive(Debug, thiserror::Error)]
#[error("{message}")]
struct InvalidRecordError {
    reason: InvalidRecordReason,
    message: String,
}

impl InvalidRecordError {
    fn new(reason: InvalidRecordReason, message: impl Into<String>) -> Self {
        Self {
            reason,
            message: message.into(),
        }
    }
}

/// Controls which mapped BAM alignment records are converted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BamToBiggOptions {
    /// Minimum mapping quality retained.
    pub min_mapq: u8,
    /// Whether records marked secondary (`0x100`) are retained.
    pub include_secondary: bool,
    /// Whether records marked supplementary (`0x800`) are retained.
    pub include_supplementary: bool,
    /// Handling of decoded records that cannot be converted to BED.
    pub invalid_record_policy: InvalidRecordPolicy,
    /// Sample/group label written to the TrackCluster sample metadata field.
    pub group: String,
}

impl Default for BamToBiggOptions {
    fn default() -> Self {
        Self {
            min_mapq: 30,
            include_secondary: false,
            include_supplementary: false,
            invalid_record_policy: InvalidRecordPolicy::Skip,
            group: "none".to_owned(),
        }
    }
}

/// Converted transcripts and deterministic BAM filtering counts.
#[derive(Clone, Debug, Default)]
pub struct BamToBiggResult {
    /// BAM records converted to BED12+8 transcripts.
    pub transcripts: Vec<Transcript>,
    /// Total decoded BAM records inspected.
    pub total_records: usize,
    /// Unmapped records omitted.
    pub skipped_unmapped: usize,
    /// Secondary records omitted by policy.
    pub skipped_secondary: usize,
    /// Supplementary records omitted by policy.
    pub skipped_supplementary: usize,
    /// Records omitted because MAPQ was below the configured threshold.
    pub skipped_below_mapq: usize,
    /// Invalid decoded records grouped by a stable exclusion reason.
    pub invalid_records: BTreeMap<InvalidRecordReason, InvalidRecordStats>,
}

impl BamToBiggResult {
    /// Total records excluded because a decoded record could not be converted.
    pub fn skipped_invalid_records(&self) -> usize {
        self.invalid_records.values().map(|stats| stats.count).sum()
    }
}

/// Deterministic counts from a streaming BAM-to-BED conversion.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BamToBiggSummary {
    /// Total decoded BAM records inspected.
    pub total_records: usize,
    /// BAM records converted to BED12+8 transcripts.
    pub converted_records: usize,
    /// Unmapped records omitted.
    pub skipped_unmapped: usize,
    /// Secondary records omitted by policy.
    pub skipped_secondary: usize,
    /// Supplementary records omitted by policy.
    pub skipped_supplementary: usize,
    /// Records omitted because MAPQ was below the configured threshold.
    pub skipped_below_mapq: usize,
    /// Invalid decoded records grouped by a stable exclusion reason.
    pub invalid_records: BTreeMap<InvalidRecordReason, InvalidRecordStats>,
}

impl BamToBiggSummary {
    /// Total records excluded because a decoded record could not be converted.
    pub fn skipped_invalid_records(&self) -> usize {
        self.invalid_records.values().map(|stats| stats.count).sum()
    }

    fn record_invalid(&mut self, record_ordinal: usize, error: &InvalidRecordError) {
        self.invalid_records
            .entry(error.reason)
            .and_modify(|stats| stats.count += 1)
            .or_insert_with(|| InvalidRecordStats {
                count: 1,
                first_record_ordinal: record_ordinal,
                first_error: error.to_string(),
            });
    }
}

fn checked_advance(position: u32, length: usize) -> Result<u32, InvalidRecordError> {
    let length = u32::try_from(length).map_err(|_| {
        InvalidRecordError::new(
            InvalidRecordReason::CoordinateOverflow,
            "CIGAR operation length exceeds u32",
        )
    })?;
    position.checked_add(length).ok_or_else(|| {
        InvalidRecordError::new(
            InvalidRecordReason::CoordinateOverflow,
            "CIGAR reference span exceeds u32",
        )
    })
}

fn validate_group(group: &str) -> anyhow::Result<&str> {
    let group = group.trim();
    if group.is_empty() {
        return Ok("none");
    }
    if group.chars().any(char::is_control) {
        anyhow::bail!("BAM group label must not contain control characters");
    }
    Ok(group)
}

fn record_to_transcript(
    record: &bam::Record,
    header: &sam::Header,
    group: &str,
    record_ordinal: usize,
) -> Result<Transcript, InvalidRecordError> {
    let raw_name = record.name().ok_or_else(|| {
        InvalidRecordError::new(
            InvalidRecordReason::MissingQueryName,
            format!("mapped BAM record {record_ordinal} has no query name"),
        )
    })?;
    let name = std::str::from_utf8(raw_name.as_ref())
        .map_err(|_| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidQueryName,
                format!("mapped BAM record {record_ordinal} has a non-UTF-8 query name"),
            )
        })?
        .to_owned();
    if name.trim().is_empty() || name == "*" {
        return Err(InvalidRecordError::new(
            InvalidRecordReason::MissingQueryName,
            format!("mapped BAM record {record_ordinal} has no query name"),
        ));
    }
    if name.chars().any(char::is_control) {
        return Err(InvalidRecordError::new(
            InvalidRecordReason::InvalidQueryName,
            format!("mapped BAM record {record_ordinal} has an unsafe query name {name:?}"),
        ));
    }
    let reference_id = record
        .reference_sequence_id()
        .transpose()
        .map_err(|error| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidReference,
                format!("decode reference ID for BAM record {record_ordinal}: {error}"),
            )
        })?
        .ok_or_else(|| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidReference,
                format!("mapped BAM record {record_ordinal} has no reference sequence"),
            )
        })?;
    let (reference_name, reference_sequence) = header
        .reference_sequences()
        .get_index(reference_id)
        .ok_or_else(|| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidReference,
                format!(
                    "BAM record {record_ordinal} references missing header sequence ID {reference_id}"
                ),
            )
        })?;
    let chrom = reference_name.to_string();
    let reference_length =
        u64::try_from(usize::from(reference_sequence.length())).map_err(|_| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidReference,
                "BAM reference length exceeds u64",
            )
        })?;
    let alignment_start = record
        .alignment_start()
        .transpose()
        .map_err(|error| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidAlignmentStart,
                format!("decode alignment start for BAM record {record_ordinal}: {error}"),
            )
        })?
        .ok_or_else(|| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidAlignmentStart,
                format!("mapped BAM record {record_ordinal} has no alignment start"),
            )
        })?;
    let start0 = u32::try_from(usize::from(alignment_start) - 1).map_err(|_| {
        InvalidRecordError::new(
            InvalidRecordReason::InvalidAlignmentStart,
            format!("BAM record {record_ordinal} alignment start exceeds BED coordinate range"),
        )
    })?;

    let mut cursor = start0;
    let mut exon_start = start0;
    let mut exon_has_alignment = false;
    let mut exons = Vec::new();
    for (operation_index, result) in record.cigar().iter().enumerate() {
        let operation = result.map_err(|error| {
            InvalidRecordError::new(
                InvalidRecordReason::InvalidCigar,
                format!(
                    "parse CIGAR operation {} for BAM record {record_ordinal} ({name:?}): {error}",
                    operation_index + 1
                ),
            )
        })?;
        if operation.is_empty() {
            return Err(InvalidRecordError::new(
                InvalidRecordReason::InvalidCigarStructure,
                format!(
                    "BAM record {record_ordinal} ({name:?}) has a zero-length CIGAR operation at position {}",
                    operation_index + 1
                ),
            ));
        }
        match operation.kind() {
            Kind::Skip => {
                if exon_start == cursor {
                    return Err(InvalidRecordError::new(
                        InvalidRecordReason::InvalidCigarStructure,
                        format!(
                            "BAM record {record_ordinal} ({name:?}) has a CIGAR N without a preceding exon"
                        ),
                    ));
                }
                if !exon_has_alignment {
                    return Err(InvalidRecordError::new(
                        InvalidRecordReason::InvalidCigarStructure,
                        format!(
                            "BAM record {record_ordinal} ({name:?}) has an exon block containing only deletions"
                        ),
                    ));
                }
                exons.push(
                    Interval::new(Coord::new(exon_start), Coord::new(cursor)).map_err(|error| {
                        InvalidRecordError::new(
                            InvalidRecordReason::InvalidCigarStructure,
                            format!(
                                "BAM record {record_ordinal} ({name:?}) has invalid exon geometry: {error}"
                            ),
                        )
                    })?,
                );
                cursor = checked_advance(cursor, operation.len())?;
                exon_start = cursor;
                exon_has_alignment = false;
            }
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                cursor = checked_advance(cursor, operation.len())?;
                exon_has_alignment = true;
            }
            kind if kind.consumes_reference() => {
                cursor = checked_advance(cursor, operation.len())?;
            }
            _ => {}
        }
    }
    if exon_start < cursor {
        if !exon_has_alignment {
            return Err(InvalidRecordError::new(
                InvalidRecordReason::InvalidCigarStructure,
                format!(
                    "BAM record {record_ordinal} ({name:?}) has an exon block containing only deletions"
                ),
            ));
        }
        exons.push(
            Interval::new(Coord::new(exon_start), Coord::new(cursor)).map_err(|error| {
                InvalidRecordError::new(
                    InvalidRecordReason::InvalidCigarStructure,
                    format!(
                        "BAM record {record_ordinal} ({name:?}) has invalid exon geometry: {error}"
                    ),
                )
            })?,
        );
    } else if !exons.is_empty() {
        return Err(InvalidRecordError::new(
            InvalidRecordReason::InvalidCigarStructure,
            format!(
                "BAM record {record_ordinal} ({name:?}) has a CIGAR ending in N without a following exon"
            ),
        ));
    }
    if exons.is_empty() {
        return Err(InvalidRecordError::new(
            InvalidRecordReason::InvalidCigarStructure,
            format!(
                "mapped BAM record {record_ordinal} ({name:?}) has no reference-consuming CIGAR operation"
            ),
        ));
    }

    let tx_start = exons[0].start;
    let tx_end = exons[exons.len() - 1].end;
    if u64::from(tx_end.get()) > reference_length {
        return Err(InvalidRecordError::new(
            InvalidRecordReason::ReferenceOutOfBounds,
            format!(
                "BAM record {record_ordinal} ({name:?}) ends at {} beyond reference {chrom:?} length {reference_length}",
                tx_end.get()
            ),
        ));
    }
    let strand = if record.flags().is_reverse_complemented() {
        Strand::Minus
    } else {
        Strand::Plus
    };
    let score = record
        .mapping_quality()
        .map(u8::from)
        .map(u32::from)
        .unwrap_or(0);
    let exon_frames = format!(
        "{},",
        std::iter::repeat_n("-1", exons.len())
            .collect::<Vec<_>>()
            .join(",")
    );
    let item_rgb = match strand {
        Strand::Plus => "250,128,114",
        Strand::Minus => "64,224,208",
        Strand::Unknown => "0",
    };

    Transcript::new(
        chrom,
        strand,
        tx_start,
        tx_end,
        name.clone(),
        exons,
        Bed12Attrs {
            score,
            thick_start: Coord::new(0),
            thick_end: Coord::new(0),
            item_rgb: item_rgb.to_owned(),
            extra_fields: vec![
                "none".to_owned(),
                "none".to_owned(),
                "none".to_owned(),
                exon_frames,
                "nanopore_read".to_owned(),
                "none".to_owned(),
                group.to_owned(),
                "none".to_owned(),
            ],
        },
    )
    .map_err(|error| {
        InvalidRecordError::new(
            InvalidRecordReason::InvalidTranscript,
            format!("BAM record {record_ordinal} ({name:?}) is not a valid transcript: {error}"),
        )
    })
}

/// Read a BAM file and convert mapped CIGAR `N` operations into transcript introns.
///
/// The converter preserves one BED record per retained BAM alignment instance.
/// Insertions and clipping do not consume reference coordinates, deletions stay
/// within exons, and only CIGAR `N` splits exon blocks.
fn scan_bam<F>(
    path: &Path,
    options: &BamToBiggOptions,
    mut emit: F,
) -> anyhow::Result<BamToBiggSummary>
where
    F: FnMut(Transcript) -> anyhow::Result<()>,
{
    let group = validate_group(&options.group)?;
    let file = File::open(path).with_context(|| format!("open BAM input {path:?}"))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .with_context(|| format!("read BAM header {path:?}"))?;
    let mut summary = BamToBiggSummary::default();

    for (record_index, result) in reader.records().enumerate() {
        let record_ordinal = record_index + 1;
        let record =
            result.with_context(|| format!("read BAM record {record_ordinal} from {path:?}"))?;
        summary.total_records += 1;
        let flags = record.flags();
        if flags.is_unmapped() {
            summary.skipped_unmapped += 1;
            continue;
        }
        if flags.is_secondary() && !options.include_secondary {
            summary.skipped_secondary += 1;
            continue;
        }
        if flags.is_supplementary() && !options.include_supplementary {
            summary.skipped_supplementary += 1;
            continue;
        }
        let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
        if mapq < options.min_mapq {
            summary.skipped_below_mapq += 1;
            continue;
        }
        let transcript = match record_to_transcript(&record, &header, group, record_ordinal) {
            Ok(transcript) => transcript,
            Err(error) if options.invalid_record_policy == InvalidRecordPolicy::Skip => {
                summary.record_invalid(record_ordinal, &error);
                continue;
            }
            Err(error) => {
                return Err(anyhow::Error::new(error)
                    .context(format!("convert BAM record {record_ordinal} from {path:?}")));
            }
        };
        emit(transcript).with_context(|| {
            format!("write BAM record {record_ordinal} converted from {path:?}")
        })?;
        summary.converted_records += 1;
    }

    let skipped_invalid = summary.skipped_invalid_records();
    if options.invalid_record_policy == InvalidRecordPolicy::Skip
        && summary.converted_records == 0
        && skipped_invalid > 0
    {
        let reasons = summary
            .invalid_records
            .iter()
            .map(|(reason, stats)| {
                format!(
                    "{}={} (first record {}: {})",
                    reason.as_str(),
                    stats.count,
                    stats.first_record_ordinal,
                    stats.first_error
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        anyhow::bail!(
            "BAM conversion produced no valid records after skipping {skipped_invalid} invalid decoded record(s): {reasons}"
        );
    }

    Ok(summary)
}

/// Stream converted BAM records to a BED12+8 writer without retaining the
/// complete alignment file in memory.
pub fn write_bam_to_bed_writer<W: std::io::Write>(
    path: &Path,
    options: &BamToBiggOptions,
    writer: &mut W,
) -> anyhow::Result<BamToBiggSummary> {
    scan_bam(path, options, |transcript| {
        crate::io::bed::write_bed12_to_writer(writer, std::iter::once(&transcript))
            .map_err(Into::into)
    })
}

/// Read a BAM file and retain all converted records in memory.
///
/// Prefer [`write_bam_to_bed_writer`] in command-line workflows processing
/// large BAM files.
pub fn read_bam(path: &Path, options: &BamToBiggOptions) -> anyhow::Result<BamToBiggResult> {
    let mut transcripts = Vec::new();
    let summary = scan_bam(path, options, |transcript| {
        transcripts.push(transcript);
        Ok(())
    })?;
    Ok(BamToBiggResult {
        transcripts,
        total_records: summary.total_records,
        skipped_unmapped: summary.skipped_unmapped,
        skipped_secondary: summary.skipped_secondary,
        skipped_supplementary: summary.skipped_supplementary,
        skipped_below_mapq: summary.skipped_below_mapq,
        invalid_records: summary.invalid_records,
    })
}
