use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript, TranscriptError};

#[derive(Error, Debug)]
/// File-level BED read, parse, and write errors.
#[allow(missing_docs)]
pub enum BedError {
    /// Opening or reading the input failed.
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        path: PathBuf,
        source: std::io::Error,
    },

    /// A record failed validation at the reported line.
    #[error("{path:?}:{line}: {source}")]
    Parse {
        path: PathBuf,
        line: usize,
        source: BedParseError,
    },

    /// A parsed read record has no usable molecule identifier.
    #[error("{path:?}:{line}: read id must not be empty")]
    InvalidReadIdentity { path: PathBuf, line: usize },

    /// Creating or writing the output failed.
    #[error("I/O error writing {path:?}: {source}")]
    IoWrite {
        path: PathBuf,
        source: std::io::Error,
    },
}

#[derive(Error, Debug)]
/// Structural error in one BED12 record.
#[allow(missing_docs)]
pub enum BedParseError {
    /// The record has fewer than twelve columns.
    #[error("expected at least 12 columns, got {got}")]
    TooFewColumns { got: usize },

    /// A required unsigned integer field is malformed.
    #[error("invalid integer for {field}: {value:?}")]
    InvalidInt { field: &'static str, value: String },

    /// BED score is outside the standard range.
    #[error("BED score must be between 0 and 1000, got {value}")]
    ScoreOutOfRange { value: u32 },

    /// `blockCount` is malformed.
    #[error("invalid blockCount: {value:?}")]
    InvalidBlockCount { value: String },

    /// A block list length does not match `blockCount`.
    #[error("blockCount {block_count} does not match {field_name} length {list_len}")]
    BlockListLengthMismatch {
        block_count: usize,
        field_name: &'static str,
        list_len: usize,
    },

    /// A block list contains an empty value before its optional trailing comma.
    #[error("empty value inside {field_name} at list position {index}")]
    EmptyBlockListToken {
        field_name: &'static str,
        index: usize,
    },

    /// Adding a block offset or size overflowed `u32`.
    #[error("block start+size overflows u32")]
    BlockOverflow,

    /// A block lies outside the declared transcript span.
    #[error(
        "block {block_index} [{block_start}, {block_end}) lies outside transcript span [{tx_start}, {tx_end})"
    )]
    BlockOutsideSpan {
        block_index: usize,
        block_start: u32,
        block_end: u32,
        tx_start: u32,
        tx_end: u32,
    },

    /// Relative block starts are not nondecreasing.
    #[error(
        "blockStarts must be nondecreasing, but block {block_index} starts at {current} after {previous}"
    )]
    BlocksOutOfOrder {
        block_index: usize,
        previous: u32,
        current: u32,
    },

    /// Coding/thick coordinates are ordered incorrectly or outside the transcript.
    #[error(
        "invalid thick span [{thick_start}, {thick_end}) for transcript [{tx_start}, {tx_end})"
    )]
    InvalidThickSpan {
        thick_start: u32,
        thick_end: u32,
        tx_start: u32,
        tx_end: u32,
    },

    /// The strand token is invalid.
    #[error(transparent)]
    Strand(#[from] crate::model::strand::StrandParseError),

    /// The resulting transcript violates geometry invariants.
    #[error(transparent)]
    Transcript(#[from] TranscriptError),
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// Controls whether selected legacy BED defects are rejected or explicitly repaired.
pub enum BedParseMode {
    /// Reject every structural defect.
    #[default]
    Strict,
    /// Repair the small documented legacy subset and emit [`BedRepairWarning`] values.
    Lenient,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// Structured description of one lenient-mode field repair.
pub struct BedRepairWarning {
    /// Source file path.
    pub path: PathBuf,
    /// One-based source line number.
    pub line: usize,
    /// BED field name.
    pub field: &'static str,
    /// Original serialized field value.
    pub original: String,
    /// Replacement serialized field value.
    pub replacement: String,
    /// Stable human-readable repair reason.
    pub reason: &'static str,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// One read track rejected by explicit recovering-read iteration.
pub struct RejectedReadRecord {
    /// Source BED path.
    pub path: PathBuf,
    /// One-based physical source line.
    pub line: usize,
    /// Read identifier when the BED record parsed far enough to recover it.
    pub read_id: Option<String>,
    /// Stable rejection class such as `parse` or `identity`.
    pub kind: &'static str,
    /// Human-readable reason for exclusion.
    pub reason: String,
}

#[derive(Clone, Debug)]
struct LineRepair {
    field: &'static str,
    original: String,
    replacement: String,
    reason: &'static str,
}

/// Streaming BED12 reader that attaches path and line context to errors.
pub struct Bed12Reader<R: BufRead> {
    path: PathBuf,
    line_number: usize,
    reader: R,
    line_buf: String,
    mode: BedParseMode,
    warnings: Vec<BedRepairWarning>,
    rejected_reads: Vec<RejectedReadRecord>,
}

const BED_READ_BUFFER_BYTES: usize = 1024 * 1024;

/// Open a BED12 file in strict mode.
pub fn read_bed12<P: AsRef<Path>>(path: P) -> Result<Bed12Reader<BufReader<File>>, BedError> {
    read_bed12_with_mode(path, BedParseMode::Strict)
}

/// Open a BED12 file using an explicit parse mode.
pub fn read_bed12_with_mode<P: AsRef<Path>>(
    path: P,
    mode: BedParseMode,
) -> Result<Bed12Reader<BufReader<File>>, BedError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| BedError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::with_capacity(BED_READ_BUFFER_BYTES, file);
    Ok(Bed12Reader {
        path,
        line_number: 0,
        reader,
        line_buf: String::new(),
        mode,
        warnings: Vec::new(),
        rejected_reads: Vec::new(),
    })
}

impl<R: BufRead> Bed12Reader<R> {
    /// Return all repairs emitted for records already consumed from this reader.
    pub fn warnings(&self) -> &[BedRepairWarning] {
        &self.warnings
    }

    /// Read the next valid read track while rejecting only record-local parse
    /// and empty-read-ID failures. File I/O errors remain fatal.
    ///
    /// This method is intentionally separate from [`Iterator::next`], whose
    /// strict behavior is unchanged for references and general BED consumers.
    pub fn next_recovering_read(&mut self) -> Result<Option<Transcript>, BedError> {
        loop {
            match self.next() {
                None => return Ok(None),
                Some(Ok(transcript)) if transcript.name.trim().is_empty() => {
                    self.rejected_reads.push(RejectedReadRecord {
                        path: self.path.clone(),
                        line: self.line_number,
                        read_id: None,
                        kind: "identity",
                        reason: "read id must not be empty".to_owned(),
                    });
                }
                Some(Ok(transcript)) => return Ok(Some(transcript)),
                Some(Err(BedError::Parse { path, line, source })) => {
                    self.rejected_reads.push(RejectedReadRecord {
                        path,
                        line,
                        read_id: read_id_hint(&self.line_buf),
                        kind: "parse",
                        reason: source.to_string(),
                    });
                }
                Some(Err(error)) => return Err(error),
            }
        }
    }

    /// Read the next read track without recovery, including strict validation
    /// of the read identifier.
    pub fn next_strict_read(&mut self) -> Result<Option<Transcript>, BedError> {
        match self.next() {
            None => Ok(None),
            Some(Ok(transcript)) if transcript.name.trim().is_empty() => {
                Err(BedError::InvalidReadIdentity {
                    path: self.path.clone(),
                    line: self.line_number,
                })
            }
            Some(Ok(transcript)) => Ok(Some(transcript)),
            Some(Err(error)) => Err(error),
        }
    }

    /// Remove and return all read-track rejections collected so far.
    pub fn take_rejected_reads(&mut self) -> Vec<RejectedReadRecord> {
        std::mem::take(&mut self.rejected_reads)
    }
}

fn read_id_hint(line: &str) -> Option<String> {
    let value = if line.contains('\t') {
        line.trim_end_matches(['\r', '\n']).split('\t').nth(3)
    } else {
        line.split_whitespace().nth(3)
    }?;
    let value = value.trim();
    (!value.is_empty()).then(|| value.to_owned())
}

/// Header used by per-gene rejected-read diagnostics.
pub const REJECTED_READS_TSV_HEADER: &str = "source_path\tline\tread_id\tkind\treason";

fn escape_rejected_tsv_field(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

/// Write auditable rejected-read diagnostics as a stable TSV.
pub fn write_rejected_reads_tsv<P: AsRef<Path>>(
    path: P,
    records: &[RejectedReadRecord],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    write_rejected_reads_tsv_to_writer(&mut writer, records)
}

/// Write auditable rejected-read diagnostics to an existing writer.
pub fn write_rejected_reads_tsv_to_writer<W: Write>(
    writer: &mut W,
    records: &[RejectedReadRecord],
) -> Result<(), std::io::Error> {
    writeln!(writer, "{REJECTED_READS_TSV_HEADER}")?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            escape_rejected_tsv_field(&record.path.to_string_lossy()),
            record.line,
            escape_rejected_tsv_field(record.read_id.as_deref().unwrap_or("")),
            record.kind,
            escape_rejected_tsv_field(&record.reason)
        )?;
    }
    writer.flush()
}

/// Validate a rejected-read TSV header and return its record count.
pub fn count_rejected_reads_tsv<P: AsRef<Path>>(path: P) -> Result<usize, std::io::Error> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    if header != REJECTED_READS_TSV_HEADER {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "invalid rejected-read header: expected {REJECTED_READS_TSV_HEADER:?}, found {header:?}"
            ),
        ));
    }
    lines.try_fold(0usize, |count, line| {
        let line = line?;
        Ok(count + usize::from(!line.trim().is_empty()))
    })
}

impl<R: BufRead> Iterator for Bed12Reader<R> {
    type Item = Result<Transcript, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();
            let read_len = match self.reader.read_line(&mut self.line_buf) {
                Ok(read_len) => read_len,
                Err(source) => {
                    return Some(Err(BedError::IoRead {
                        path: self.path.clone(),
                        source,
                    }))
                }
            };
            if read_len == 0 {
                return None;
            }
            self.line_number += 1;

            let line = self.line_buf.trim_end_matches(['\r', '\n']);
            let control = line.trim();
            if control.is_empty() || control.starts_with('#') {
                continue;
            }

            let mut repairs = Vec::new();
            let record = match parse_bed12_line_with_mode(line, self.mode, &mut repairs) {
                Ok(record) => record,
                Err(source) => {
                    return Some(Err(BedError::Parse {
                        path: self.path.clone(),
                        line: self.line_number,
                        source,
                    }))
                }
            };

            self.warnings
                .extend(repairs.into_iter().map(|repair| BedRepairWarning {
                    path: self.path.clone(),
                    line: self.line_number,
                    field: repair.field,
                    original: repair.original,
                    replacement: repair.replacement,
                    reason: repair.reason,
                }));

            return Some(Ok(record));
        }
    }
}

fn parse_u32(field: &'static str, value: &str) -> Result<u32, BedParseError> {
    value.parse::<u32>().map_err(|_| BedParseError::InvalidInt {
        field,
        value: value.to_owned(),
    })
}

fn parse_usize(field: &'static str, value: &str) -> Result<usize, BedParseError> {
    value
        .parse::<usize>()
        .map_err(|_| BedParseError::InvalidInt {
            field,
            value: value.to_owned(),
        })
}

fn parse_comma_u32_list(
    field_name: &'static str,
    value: &str,
    mode: BedParseMode,
    repairs: &mut Vec<LineRepair>,
) -> Result<Vec<u32>, BedParseError> {
    let tokens: Vec<&str> = value.split(',').collect();
    let mut parsed = Vec::with_capacity(tokens.len());
    let normalized = || {
        tokens
            .iter()
            .filter(|token| !token.is_empty())
            .copied()
            .collect::<Vec<_>>()
            .join(",")
    };

    for (index, token) in tokens.iter().enumerate() {
        if token.is_empty() {
            if index + 1 == tokens.len() {
                continue;
            }
            if mode == BedParseMode::Strict {
                return Err(BedParseError::EmptyBlockListToken { field_name, index });
            }
            repairs.push(LineRepair {
                field: field_name,
                original: value.to_owned(),
                replacement: normalized(),
                reason: "removed empty interior list value",
            });
            continue;
        }
        parsed.push(parse_u32(field_name, token)?);
    }
    Ok(parsed)
}

#[cfg(test)]
fn parse_bed12_line(line: &str) -> Result<Transcript, BedParseError> {
    parse_bed12_line_with_mode(line, BedParseMode::Strict, &mut Vec::new())
}

fn parse_bed12_line_with_mode(
    line: &str,
    mode: BedParseMode,
    repairs: &mut Vec<LineRepair>,
) -> Result<Transcript, BedParseError> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() == 1 {
        fields = line.split_whitespace().collect();
    }

    if fields.len() < 12 {
        return Err(BedParseError::TooFewColumns { got: fields.len() });
    }

    let chrom = fields[0].to_owned();
    let tx_start = Coord::new(parse_u32("chromStart", fields[1])?);
    let tx_end = Coord::new(parse_u32("chromEnd", fields[2])?);
    let name = fields[3].to_owned();
    let score = match parse_u32("score", fields[4]) {
        Ok(value @ 0..=1000) => value,
        Ok(value) if mode == BedParseMode::Strict => {
            return Err(BedParseError::ScoreOutOfRange { value });
        }
        Ok(value) => {
            repairs.push(LineRepair {
                field: "score",
                original: value.to_string(),
                replacement: "1000".to_owned(),
                reason: "clamped score to BED range",
            });
            1000
        }
        Err(error) if mode == BedParseMode::Strict => return Err(error),
        Err(_) => {
            repairs.push(LineRepair {
                field: "score",
                original: fields[4].to_owned(),
                replacement: "0".to_owned(),
                reason: "replaced invalid score",
            });
            0
        }
    };
    let strand = Strand::try_from(fields[5])?;
    let thick_start = Coord::new(parse_u32("thickStart", fields[6])?);
    let thick_end = Coord::new(parse_u32("thickEnd", fields[7])?);
    let thick_is_legacy_sentinel = thick_start.get() == 0 && thick_end.get() == 0;
    if !thick_is_legacy_sentinel
        && (thick_start > thick_end || thick_start < tx_start || thick_end > tx_end)
    {
        return Err(BedParseError::InvalidThickSpan {
            thick_start: thick_start.get(),
            thick_end: thick_end.get(),
            tx_start: tx_start.get(),
            tx_end: tx_end.get(),
        });
    }
    let item_rgb = fields[8].to_owned();
    let block_count =
        parse_usize("blockCount", fields[9]).map_err(|_| BedParseError::InvalidBlockCount {
            value: fields[9].to_owned(),
        })?;

    let block_sizes = parse_comma_u32_list("blockSizes", fields[10], mode, repairs)?;
    if block_sizes.len() != block_count {
        return Err(BedParseError::BlockListLengthMismatch {
            block_count,
            field_name: "blockSizes",
            list_len: block_sizes.len(),
        });
    }

    let block_starts = parse_comma_u32_list("blockStarts", fields[11], mode, repairs)?;
    if block_starts.len() != block_count {
        return Err(BedParseError::BlockListLengthMismatch {
            block_count,
            field_name: "blockStarts",
            list_len: block_starts.len(),
        });
    }

    let mut exons = Vec::with_capacity(block_count);
    for i in 0..block_count {
        let rel_start = block_starts[i];
        let block_size = block_sizes[i];

        if i > 0 && rel_start < block_starts[i - 1] {
            return Err(BedParseError::BlocksOutOfOrder {
                block_index: i,
                previous: block_starts[i - 1],
                current: rel_start,
            });
        }

        let exon_start_u32 = tx_start
            .get()
            .checked_add(rel_start)
            .ok_or(BedParseError::BlockOverflow)?;
        let exon_end_u32 = exon_start_u32
            .checked_add(block_size)
            .ok_or(BedParseError::BlockOverflow)?;

        let exon_end_u32 = if exon_start_u32 < tx_start.get() || exon_end_u32 > tx_end.get() {
            if mode == BedParseMode::Strict {
                return Err(BedParseError::BlockOutsideSpan {
                    block_index: i,
                    block_start: exon_start_u32,
                    block_end: exon_end_u32,
                    tx_start: tx_start.get(),
                    tx_end: tx_end.get(),
                });
            }
            if exon_start_u32 >= tx_end.get() {
                return Err(BedParseError::BlockOutsideSpan {
                    block_index: i,
                    block_start: exon_start_u32,
                    block_end: exon_end_u32,
                    tx_start: tx_start.get(),
                    tx_end: tx_end.get(),
                });
            }
            let clamped_end = exon_end_u32.min(tx_end.get());
            let clamped_size = clamped_end - exon_start_u32;
            repairs.push(LineRepair {
                field: "blockSizes",
                original: block_size.to_string(),
                replacement: clamped_size.to_string(),
                reason: "clamped block end to chromEnd",
            });
            clamped_end
        } else {
            exon_end_u32
        };

        let exon = Interval::new(Coord::new(exon_start_u32), Coord::new(exon_end_u32))
            .map_err(|_| BedParseError::BlockOverflow)?;
        exons.push(exon);
    }

    let extra_fields = fields[12..]
        .iter()
        .map(|value| (*value).to_owned())
        .collect();

    Ok(Transcript::new(
        chrom,
        strand,
        tx_start,
        tx_end,
        name,
        exons,
        Bed12Attrs {
            score,
            thick_start,
            thick_end,
            item_rgb,
            extra_fields,
        },
    )?)
}

/// Write transcripts to a BED12/bigGenePred-compatible file.
pub fn write_bed12<'a, P, I>(path: P, transcripts: I) -> Result<(), BedError>
where
    P: AsRef<Path>,
    I: IntoIterator<Item = &'a Transcript>,
{
    let path = path.as_ref().to_path_buf();
    let file = File::create(&path).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    let mut writer = std::io::BufWriter::new(file);
    write_bed12_to_writer(&mut writer, transcripts).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    writer
        .flush()
        .map_err(|source| BedError::IoWrite { path, source })?;
    Ok(())
}

/// Serialize transcripts to an existing writer.
pub fn write_bed12_to_writer<'a, W, I>(writer: &mut W, transcripts: I) -> Result<(), std::io::Error>
where
    W: Write,
    I: IntoIterator<Item = &'a Transcript>,
{
    for transcript in transcripts {
        let block_count = transcript.exons.len();
        let mut block_sizes = String::new();
        let mut block_starts = String::new();

        for exon in &transcript.exons {
            let size = exon.len();
            let rel_start = exon.start.get().saturating_sub(transcript.tx_start.get());
            block_sizes.push_str(&format!("{size},"));
            block_starts.push_str(&format!("{rel_start},"));
        }

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            transcript.chrom,
            transcript.tx_start.get(),
            transcript.tx_end.get(),
            transcript.name,
            transcript.score,
            transcript.strand.as_char(),
            transcript.thick_start.get(),
            transcript.thick_end.get(),
            transcript.item_rgb,
            block_count,
            block_sizes,
            block_starts
        )?;

        for extra in &transcript.extra_fields {
            write!(writer, "\t{extra}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn test_reader(input: &str) -> Bed12Reader<Cursor<Vec<u8>>> {
        Bed12Reader {
            path: PathBuf::from("reads.bed"),
            line_number: 0,
            reader: Cursor::new(input.as_bytes().to_vec()),
            line_buf: String::new(),
            mode: BedParseMode::Strict,
            warnings: Vec::new(),
            rejected_reads: Vec::new(),
        }
    }

    #[test]
    fn bed12_roundtrip_normalized() {
        let input = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t50,30,\t0,70,\n";
        let transcript = parse_bed12_line(input.trim()).unwrap();

        let mut buffer = Vec::new();
        write_bed12_to_writer(&mut buffer, [&transcript]).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        let reparsed = parse_bed12_line(output.trim()).unwrap();
        assert_eq!(transcript, reparsed);
    }

    #[test]
    fn reader_preserves_trailing_empty_extension_columns() {
        let input = concat!(
            "  # comment\r\n",
            " \t \r\n",
            "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t100,\t0,\tfuture\t\t\r\n",
        );
        let mut reader = test_reader(input);
        let transcript = reader.next().unwrap().unwrap();

        assert_eq!(transcript.extra_fields, ["future", "", ""]);
        assert!(reader.next().is_none());

        let mut buffer = Vec::new();
        write_bed12_to_writer(&mut buffer, [&transcript]).unwrap();
        assert_eq!(
            String::from_utf8(buffer).unwrap(),
            "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t100,\t0,\tfuture\t\t\n"
        );
    }

    #[test]
    fn recovering_reader_skips_only_bad_tracks_and_keeps_physical_lines() {
        let input = concat!(
            "# comment\n",
            "\n",
            "not-a-bed-record\n",
            "chr1\t0\t10\tgood1\t0\t+\t0\t0\t0\t1\t10,\t0,\n",
            "chr1\t10\t20\t\t0\t+\t0\t0\t0\t1\t10,\t0,\n",
            "chr1\t20\t30\tbad_score\tNaN\t+\t0\t0\t0\t1\t10,\t0,\n",
            "chr1\t30\t40\tgood2\t0\t+\t0\t0\t0\t1\t10,\t0,\n",
        );
        let mut reader = test_reader(input);
        let mut names = Vec::new();
        while let Some(read) = reader.next_recovering_read().unwrap() {
            names.push(read.name);
        }
        assert_eq!(names, ["good1", "good2"]);
        let rejected = reader.take_rejected_reads();
        assert_eq!(rejected.len(), 3);
        assert_eq!(
            rejected
                .iter()
                .map(|record| (record.line, record.kind))
                .collect::<Vec<_>>(),
            [(3, "parse"), (5, "identity"), (6, "parse")]
        );
        assert_eq!(rejected[2].read_id.as_deref(), Some("bad_score"));
    }

    #[test]
    fn strict_read_iteration_rejects_empty_identity() {
        let mut reader = test_reader("chr1\t0\t10\t\t0\t+\t0\t0\t0\t1\t10,\t0,\n");
        assert!(matches!(
            reader.next_strict_read(),
            Err(BedError::InvalidReadIdentity { line: 1, .. })
        ));
    }

    #[test]
    fn recovering_reader_never_swallows_io_failures() {
        struct FailingRead {
            bytes: Cursor<Vec<u8>>,
        }

        impl Read for FailingRead {
            fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
                let read = self.bytes.read(buffer)?;
                if read == 0 {
                    return Err(std::io::Error::other("synthetic read failure"));
                }
                Ok(read)
            }
        }

        let input = b"chr1\t0\t10\tgood\t0\t+\t0\t0\t0\t1\t10,\t0,\n".to_vec();
        let mut reader = Bed12Reader {
            path: PathBuf::from("failing.bed"),
            line_number: 0,
            reader: BufReader::with_capacity(
                1,
                FailingRead {
                    bytes: Cursor::new(input),
                },
            ),
            line_buf: String::new(),
            mode: BedParseMode::Strict,
            warnings: Vec::new(),
            rejected_reads: Vec::new(),
        };
        assert_eq!(reader.next_recovering_read().unwrap().unwrap().name, "good");
        assert!(matches!(
            reader.next_recovering_read(),
            Err(BedError::IoRead { .. })
        ));
        assert!(reader.take_rejected_reads().is_empty());
    }

    #[test]
    fn rejected_read_tsv_round_trips_count_and_escapes_fields() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-rejected-{}-{nanos}.tsv",
            std::process::id()
        ));
        let records = vec![RejectedReadRecord {
            path: PathBuf::from("source\treads.bed"),
            line: 7,
            read_id: Some("read\\id".to_owned()),
            kind: "parse",
            reason: "bad\tfield\ncontinued".to_owned(),
        }];
        write_rejected_reads_tsv(&path, &records).unwrap();
        assert_eq!(count_rejected_reads_tsv(&path).unwrap(), 1);
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with(&format!("{REJECTED_READS_TSV_HEADER}\n")));
        assert!(
            content.contains("source\\treads.bed\t7\tread\\\\id\tparse\tbad\\tfield\\ncontinued")
        );

        write_rejected_reads_tsv(&path, &[]).unwrap();
        assert_eq!(count_rejected_reads_tsv(&path).unwrap(), 0);
        std::fs::write(&path, "wrong\n").unwrap();
        assert_eq!(
            count_rejected_reads_tsv(&path).unwrap_err().kind(),
            std::io::ErrorKind::InvalidData
        );
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn strict_parser_rejects_invalid_score_and_out_of_range_score() {
        let invalid = "chr1\t100\t200\ttx1\tnot-a-score\t+\t100\t200\t0\t1\t100,\t0,";
        assert!(matches!(
            parse_bed12_line(invalid),
            Err(BedParseError::InvalidInt { field: "score", .. })
        ));

        let out_of_range = "chr1\t100\t200\ttx1\t1001\t+\t100\t200\t0\t1\t100,\t0,";
        assert!(matches!(
            parse_bed12_line(out_of_range),
            Err(BedParseError::ScoreOutOfRange { value: 1001 })
        ));
    }

    #[test]
    fn lenient_parser_reports_every_legacy_repair() {
        let input = "chr1\t100\t200\ttx1\tbad\t+\t100\t200\t0\t2\t50,,30,\t0,70,";
        let mut repairs = Vec::new();
        let transcript =
            parse_bed12_line_with_mode(input, BedParseMode::Lenient, &mut repairs).unwrap();

        assert_eq!(transcript.score, 0);
        assert_eq!(transcript.exons.len(), 2);
        assert_eq!(repairs.len(), 2);
        assert!(repairs.iter().any(|repair| repair.field == "score"));
        assert!(repairs.iter().any(|repair| repair.field == "blockSizes"));
    }

    #[test]
    fn strict_parser_rejects_empty_tokens_and_blocks_outside_span() {
        let empty_token = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t50,,30,\t0,70,";
        assert!(matches!(
            parse_bed12_line(empty_token),
            Err(BedParseError::EmptyBlockListToken {
                field_name: "blockSizes",
                index: 1
            })
        ));

        let outside = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t150,\t0,";
        assert!(matches!(
            parse_bed12_line(outside),
            Err(BedParseError::BlockOutsideSpan {
                block_index: 0,
                block_start: 100,
                block_end: 250,
                tx_start: 100,
                tx_end: 200
            })
        ));
    }

    #[test]
    fn strict_parser_rejects_zero_length_and_overlapping_exons() {
        let zero = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t0,100,\t0,0,";
        assert!(matches!(
            parse_bed12_line(zero),
            Err(BedParseError::Transcript(TranscriptError::EmptyExon { .. }))
        ));

        let overlap = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t70,60,\t0,40,";
        assert!(matches!(
            parse_bed12_line(overlap),
            Err(BedParseError::Transcript(
                TranscriptError::OverlappingExons { .. }
            ))
        ));

        let out_of_order = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t40,40,\t60,0,";
        assert!(matches!(
            parse_bed12_line(out_of_order),
            Err(BedParseError::BlocksOutOfOrder {
                block_index: 1,
                previous: 60,
                current: 0
            })
        ));
    }

    #[test]
    fn strict_parser_rejects_invalid_thick_span_and_exon_bounds() {
        let invalid_thick = "chr1\t100\t200\ttx1\t0\t+\t90\t200\t0\t1\t100,\t0,";
        assert!(matches!(
            parse_bed12_line(invalid_thick),
            Err(BedParseError::InvalidThickSpan { .. })
        ));

        let mismatched_bounds = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t80,\t10,";
        assert!(matches!(
            parse_bed12_line(mismatched_bounds),
            Err(BedParseError::Transcript(
                TranscriptError::ExonBoundsDoNotMatchSpan { .. }
            ))
        ));
    }
}
