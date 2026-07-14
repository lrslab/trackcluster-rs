//! Pure-Rust conversion of spliced BAM alignments into TrackCluster transcripts.

use std::fs::File;
use std::path::Path;

use anyhow::Context;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

/// Controls which mapped BAM alignment records are converted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BamToBiggOptions {
    /// Minimum mapping quality retained.
    pub min_mapq: u8,
    /// Whether records marked secondary (`0x100`) are retained.
    pub include_secondary: bool,
    /// Whether records marked supplementary (`0x800`) are retained.
    pub include_supplementary: bool,
    /// Sample/group label written to the TrackCluster sample metadata field.
    pub group: String,
}

impl Default for BamToBiggOptions {
    fn default() -> Self {
        Self {
            min_mapq: 30,
            include_secondary: false,
            include_supplementary: false,
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
}

fn checked_advance(position: u32, length: usize) -> anyhow::Result<u32> {
    let length = u32::try_from(length).context("CIGAR operation length exceeds u32")?;
    position
        .checked_add(length)
        .context("CIGAR reference span exceeds u32")
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
) -> anyhow::Result<Transcript> {
    let name = record
        .name()
        .map(ToString::to_string)
        .filter(|name| !name.trim().is_empty() && name != "*")
        .with_context(|| format!("mapped BAM record {record_ordinal} has no query name"))?;
    if name.chars().any(char::is_control) {
        anyhow::bail!("mapped BAM record {record_ordinal} has an unsafe query name {name:?}");
    }
    let reference_id = record
        .reference_sequence_id()
        .transpose()
        .with_context(|| format!("decode reference ID for BAM record {record_ordinal}"))?
        .with_context(|| format!("mapped BAM record {record_ordinal} has no reference sequence"))?;
    let (reference_name, reference_sequence) = header
        .reference_sequences()
        .get_index(reference_id)
        .with_context(|| {
            format!(
                "BAM record {record_ordinal} references missing header sequence ID {reference_id}"
            )
        })?;
    let chrom = reference_name.to_string();
    let reference_length = u64::try_from(usize::from(reference_sequence.length()))
        .context("BAM reference length exceeds u64")?;
    let alignment_start = record
        .alignment_start()
        .transpose()
        .with_context(|| format!("decode alignment start for BAM record {record_ordinal}"))?
        .with_context(|| format!("mapped BAM record {record_ordinal} has no alignment start"))?;
    let start0 = u32::try_from(usize::from(alignment_start) - 1)
        .context("BAM alignment start exceeds BED coordinate range")?;

    let mut cursor = start0;
    let mut exon_start = start0;
    let mut exon_has_alignment = false;
    let mut exons = Vec::new();
    for result in record.cigar().iter() {
        let operation = result
            .with_context(|| format!("parse CIGAR for BAM record {record_ordinal} ({name:?})"))?;
        match operation.kind() {
            Kind::Skip => {
                if exon_start == cursor {
                    anyhow::bail!(
                        "BAM record {record_ordinal} ({name:?}) has a CIGAR N without a preceding exon"
                    );
                }
                if !exon_has_alignment {
                    anyhow::bail!(
                        "BAM record {record_ordinal} ({name:?}) has an exon block containing only deletions"
                    );
                }
                exons.push(Interval::new(Coord::new(exon_start), Coord::new(cursor))?);
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
            anyhow::bail!(
                "BAM record {record_ordinal} ({name:?}) has an exon block containing only deletions"
            );
        }
        exons.push(Interval::new(Coord::new(exon_start), Coord::new(cursor))?);
    } else if !exons.is_empty() {
        anyhow::bail!(
            "BAM record {record_ordinal} ({name:?}) has a CIGAR ending in N without a following exon"
        );
    }
    if exons.is_empty() {
        anyhow::bail!(
            "mapped BAM record {record_ordinal} ({name:?}) has no reference-consuming CIGAR operation"
        );
    }

    let tx_start = exons[0].start;
    let tx_end = exons[exons.len() - 1].end;
    if u64::from(tx_end.get()) > reference_length {
        anyhow::bail!(
            "BAM record {record_ordinal} ({name:?}) ends at {} beyond reference {chrom:?} length {reference_length}",
            tx_end.get()
        );
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
        name,
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
    .map_err(Into::into)
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
        let transcript = record_to_transcript(&record, &header, group, record_ordinal)
            .with_context(|| format!("convert BAM record {record_ordinal} from {path:?}"))?;
        emit(transcript).with_context(|| {
            format!("write BAM record {record_ordinal} converted from {path:?}")
        })?;
        summary.converted_records += 1;
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
    })
}
