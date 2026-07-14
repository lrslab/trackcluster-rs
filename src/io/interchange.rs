//! Standards-oriented transcript catalog exports.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Context;

use crate::model::Transcript;

fn gene_id(transcript: &Transcript) -> &str {
    transcript.metadata().gene_id().unwrap_or("unassigned")
}

fn gtf_escape(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}

fn gff3_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-' | b':') {
            out.push(char::from(byte));
        } else {
            out.push_str(&format!("%{byte:02X}"));
        }
    }
    out
}

fn synced_writer(path: &Path) -> anyhow::Result<BufWriter<File>> {
    let file = File::create(path).with_context(|| format!("create {path:?}"))?;
    Ok(BufWriter::new(file))
}

fn finish_writer(mut writer: BufWriter<File>) -> anyhow::Result<()> {
    writer.flush()?;
    writer.get_ref().sync_all()?;
    Ok(())
}

/// Write a GTF 2.2 transcript catalog with standardized `gene_id` and
/// `transcript_id` attributes.
///
/// BED's zero-based, half-open coordinates are converted to GTF's one-based,
/// closed coordinates. Each transcript is followed by its exon features.
pub fn write_gtf(path: &Path, transcripts: &[Transcript]) -> anyhow::Result<()> {
    let mut writer = synced_writer(path)?;
    write_gtf_to_writer(&mut writer, transcripts)?;
    finish_writer(writer)
}

/// Write a GTF 2.2 transcript catalog to an existing writer.
pub fn write_gtf_to_writer<W: Write>(
    writer: &mut W,
    transcripts: &[Transcript],
) -> anyhow::Result<()> {
    writeln!(writer, "# trackcluster-rs GTF 2.2")?;
    for transcript in transcripts {
        let gene = gtf_escape(gene_id(transcript));
        let transcript_id = gtf_escape(&transcript.name);
        let attributes = format!("gene_id \"{gene}\"; transcript_id \"{transcript_id}\";");
        writeln!(
            writer,
            "{}\ttrackcluster-rs\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
            transcript.chrom,
            transcript.tx_start.get() + 1,
            transcript.tx_end.get(),
            transcript.strand.as_char(),
            attributes
        )?;
        for (index, exon) in transcript.exons.iter().enumerate() {
            writeln!(
                writer,
                "{}\ttrackcluster-rs\texon\t{}\t{}\t.\t{}\t.\t{} exon_number \"{}\";",
                transcript.chrom,
                exon.start.get() + 1,
                exon.end.get(),
                transcript.strand.as_char(),
                attributes,
                index + 1
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write a GFF3 transcript catalog using `mRNA` and `exon` features.
pub fn write_gff3(path: &Path, transcripts: &[Transcript]) -> anyhow::Result<()> {
    let mut writer = synced_writer(path)?;
    write_gff3_to_writer(&mut writer, transcripts)?;
    finish_writer(writer)
}

/// Write a GFF3 transcript catalog to an existing writer.
pub fn write_gff3_to_writer<W: Write>(
    writer: &mut W,
    transcripts: &[Transcript],
) -> anyhow::Result<()> {
    writeln!(writer, "##gff-version 3")?;
    for transcript in transcripts {
        let transcript_id = gff3_escape(&transcript.name);
        let gene = gff3_escape(gene_id(transcript));
        writeln!(
            writer,
            "{}\ttrackcluster-rs\tmRNA\t{}\t{}\t.\t{}\t.\tID={};gene_id={}",
            transcript.chrom,
            transcript.tx_start.get() + 1,
            transcript.tx_end.get(),
            transcript.strand.as_char(),
            transcript_id,
            gene
        )?;
        for (index, exon) in transcript.exons.iter().enumerate() {
            writeln!(
                writer,
                "{}\ttrackcluster-rs\texon\t{}\t{}\t.\t{}\t.\tID={}.exon{};Parent={}",
                transcript.chrom,
                exon.start.get() + 1,
                exon.end.get(),
                transcript.strand.as_char(),
                transcript_id,
                index + 1,
                transcript_id
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write the compact transcript table used to audit inputs passed to SQANTI3.
///
/// The companion GTF from [`write_gtf`] is the actual SQANTI3 classification
/// input; this table makes identifiers and basic geometry easy to validate.
pub fn write_sqanti3_input_table(path: &Path, transcripts: &[Transcript]) -> anyhow::Result<()> {
    let mut writer = synced_writer(path)?;
    write_sqanti3_input_table_to_writer(&mut writer, transcripts)?;
    finish_writer(writer)
}

/// Write the SQANTI3 audit table to an existing writer.
pub fn write_sqanti3_input_table_to_writer<W: Write>(
    writer: &mut W,
    transcripts: &[Transcript],
) -> anyhow::Result<()> {
    writeln!(writer, "#schema\ttrackcluster-sqanti-input-v1")?;
    writeln!(
        writer,
        "isoform_id\tgene_id\tchrom\tstrand\tlength\texon_count"
    )?;
    for transcript in transcripts {
        let length: u64 = transcript
            .exons
            .iter()
            .map(|exon| u64::from(exon.len()))
            .sum();
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            transcript.name,
            gene_id(transcript),
            transcript.chrom,
            transcript.strand.as_char(),
            length,
            transcript.exons.len()
        )?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::model::{Bed12Attrs, Coord, Interval, Strand};

    use super::*;

    fn fixture() -> Transcript {
        Transcript::new(
            "chr1".to_owned(),
            Strand::Minus,
            Coord::new(99),
            Coord::new(220),
            "novel transcript/1".to_owned(),
            vec![
                Interval::new(Coord::new(99), Coord::new(120)).unwrap(),
                Interval::new(Coord::new(199), Coord::new(220)).unwrap(),
            ],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(0),
                thick_end: Coord::new(0),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "nanopore_read".to_owned(),
                    "GENE\"A".to_owned(),
                ],
            },
        )
        .unwrap()
    }

    fn temp_path(extension: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "trackcluster-interchange-{}-{nonce}.{extension}",
            std::process::id()
        ))
    }

    #[test]
    fn gtf_and_gff3_convert_coordinates_and_escape_ids() {
        let transcript = fixture();
        let gtf = temp_path("gtf");
        let gff3 = temp_path("gff3");
        write_gtf(&gtf, std::slice::from_ref(&transcript)).unwrap();
        write_gff3(&gff3, std::slice::from_ref(&transcript)).unwrap();

        let gtf_text = fs::read_to_string(&gtf).unwrap();
        assert!(gtf_text.contains("transcript\t100\t220\t.\t-"));
        assert!(gtf_text.contains("gene_id \"GENE\\\"A\""));
        let gff3_text = fs::read_to_string(&gff3).unwrap();
        assert!(gff3_text.contains("ID=novel%20transcript%2F1"));
        assert!(gff3_text.contains("exon\t100\t120"));

        fs::remove_file(gtf).unwrap();
        fs::remove_file(gff3).unwrap();
    }
}
