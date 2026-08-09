//! Indexed FASTA access for genomic reference-base validation.

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::model::Strand;

#[derive(Clone, Debug)]
struct FaiRecord {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// Random-access reader backed by a samtools-compatible `.fai` index.
#[derive(Debug)]
pub struct IndexedFasta {
    reader: BufReader<File>,
    records: BTreeMap<String, FaiRecord>,
}

fn fai_path(path: &Path) -> PathBuf {
    let mut value: OsString = path.as_os_str().to_os_string();
    value.push(".fai");
    PathBuf::from(value)
}

fn parse_u64(field: &str, value: &str, line_number: usize) -> anyhow::Result<u64> {
    value
        .parse::<u64>()
        .with_context(|| format!("parse FASTA index {field} {value:?} at line {line_number}"))
}

impl IndexedFasta {
    /// Open an uncompressed FASTA and its adjacent `<fasta>.fai` index.
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        let index_path = fai_path(path);
        let index =
            File::open(&index_path).with_context(|| format!("open FASTA index {index_path:?}"))?;
        let mut records = BTreeMap::new();
        for (line_index, line) in BufReader::new(index).lines().enumerate() {
            let line_number = line_index + 1;
            let line = line
                .with_context(|| format!("read FASTA index {index_path:?} line {line_number}"))?;
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() < 5 {
                anyhow::bail!(
                    "FASTA index {index_path:?} line {line_number} has {} fields; expected at least 5",
                    fields.len()
                );
            }
            let name = fields[0];
            if name.is_empty() || name.chars().any(char::is_control) {
                anyhow::bail!(
                    "FASTA index {index_path:?} line {line_number} has invalid sequence name {name:?}"
                );
            }
            let record = FaiRecord {
                length: parse_u64("length", fields[1], line_number)?,
                offset: parse_u64("offset", fields[2], line_number)?,
                line_bases: parse_u64("line_bases", fields[3], line_number)?,
                line_width: parse_u64("line_width", fields[4], line_number)?,
            };
            if record.line_bases == 0 || record.line_width < record.line_bases {
                anyhow::bail!(
                    "FASTA index {index_path:?} line {line_number} has invalid line geometry {} bases in {} bytes",
                    record.line_bases,
                    record.line_width
                );
            }
            if records.insert(name.to_owned(), record).is_some() {
                anyhow::bail!("FASTA index {index_path:?} contains duplicate sequence {name:?}");
            }
        }
        if records.is_empty() {
            anyhow::bail!("FASTA index {index_path:?} contains no sequences");
        }

        let fasta = File::open(path).with_context(|| format!("open FASTA {path:?}"))?;
        Ok(Self {
            reader: BufReader::new(fasta),
            records,
        })
    }

    /// Read one uppercase genomic reference base at a zero-based coordinate.
    pub fn base(&mut self, chrom: &str, pos0: u32) -> anyhow::Result<u8> {
        let record = self
            .records
            .get(chrom)
            .with_context(|| format!("reference FASTA has no sequence {chrom:?}"))?;
        let pos0 = u64::from(pos0);
        if pos0 >= record.length {
            anyhow::bail!(
                "reference coordinate {chrom}:{pos0} is outside sequence length {}",
                record.length
            );
        }
        let byte_offset = record.offset
            + (pos0 / record.line_bases) * record.line_width
            + pos0 % record.line_bases;
        self.reader
            .seek(SeekFrom::Start(byte_offset))
            .with_context(|| format!("seek FASTA to {chrom}:{pos0}"))?;
        let mut base = [0u8; 1];
        self.reader
            .read_exact(&mut base)
            .with_context(|| format!("read FASTA base at {chrom}:{pos0}"))?;
        let base = base[0].to_ascii_uppercase();
        if !base.is_ascii_alphabetic() {
            anyhow::bail!(
                "reference FASTA base at {chrom}:{pos0} is not alphabetic: byte {}",
                base
            );
        }
        Ok(base)
    }

    /// Read one base in the orientation of the genomic strand.
    pub fn oriented_base(&mut self, chrom: &str, pos0: u32, strand: Strand) -> anyhow::Result<u8> {
        let base = self.base(chrom, pos0)?;
        match strand {
            Strand::Plus => Ok(base),
            Strand::Minus => complement(base)
                .with_context(|| format!("complement reference base at {chrom}:{pos0}")),
            Strand::Unknown => {
                anyhow::bail!("cannot orient reference base at {chrom}:{pos0} on unknown strand")
            }
        }
    }
}

fn complement(base: u8) -> anyhow::Result<u8> {
    match base {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' | b'U' => Ok(b'A'),
        b'R' => Ok(b'Y'),
        b'Y' => Ok(b'R'),
        b'M' => Ok(b'K'),
        b'K' => Ok(b'M'),
        b'S' => Ok(b'S'),
        b'W' => Ok(b'W'),
        b'B' => Ok(b'V'),
        b'V' => Ok(b'B'),
        b'D' => Ok(b'H'),
        b'H' => Ok(b'D'),
        b'N' => Ok(b'N'),
        _ => anyhow::bail!("unsupported reference base byte {base}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn reads_wrapped_fasta_bases_in_both_orientations() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root =
            std::env::temp_dir().join(format!("trackcluster-fasta-{}-{nonce}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let fasta = root.join("reference.fa");
        fs::write(&fasta, b">chr1\nACGT\nNATC\n").unwrap();
        fs::write(fai_path(&fasta), b"chr1\t8\t6\t4\t5\n").unwrap();

        let mut reader = IndexedFasta::open(&fasta).unwrap();
        assert_eq!(reader.base("chr1", 0).unwrap(), b'A');
        assert_eq!(reader.base("chr1", 4).unwrap(), b'N');
        assert_eq!(reader.oriented_base("chr1", 6, Strand::Plus).unwrap(), b'T');
        assert_eq!(
            reader.oriented_base("chr1", 6, Strand::Minus).unwrap(),
            b'A'
        );
        assert!(reader.base("chr1", 8).is_err());
        assert!(reader.base("missing", 0).is_err());

        fs::remove_dir_all(root).unwrap();
    }
}
