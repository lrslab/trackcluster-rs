use std::io::Write;
use std::path::Path;

use crate::io::bed::BedError;
use crate::model::Transcript;

pub fn write_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
    pairs: &[(String, String)],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    write_read_to_isoform_tsv_writer(&mut writer, pairs)
}

pub fn write_read_to_isoform_tsv_writer<W: Write>(
    writer: &mut W,
    pairs: &[(String, String)],
) -> Result<(), std::io::Error> {
    for (read, isoform) in pairs {
        if read.is_empty()
            || isoform.is_empty()
            || read.contains(['\t', '\r', '\n'])
            || isoform.contains(['\t', '\r', '\n'])
        {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "read-to-isoform TSV fields must be non-empty and contain no tabs or line breaks: read={read:?}, isoform={isoform:?}"
                ),
            ));
        }
        writeln!(writer, "{read}\t{isoform}")?;
    }
    writer.flush()
}

pub fn write_isoforms_bed<P: AsRef<Path>>(
    path: P,
    isoforms: &[Transcript],
) -> Result<(), BedError> {
    crate::io::bed::write_bed12(path, isoforms.iter())
}

pub fn write_isoforms_bed_to_writer<W: Write>(
    writer: &mut W,
    isoforms: &[Transcript],
) -> Result<(), std::io::Error> {
    crate::io::bed::write_bed12_to_writer(writer, isoforms.iter())
}
