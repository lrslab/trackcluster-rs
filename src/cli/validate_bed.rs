use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Input BED12 / bigGenePred file
    #[arg(short, long)]
    pub input: PathBuf,

    /// Accept a limited set of legacy repairs and report every changed field
    #[arg(long)]
    pub lenient: bool,

    /// Write a tab-delimited validation report, including counts for errors and repaired records
    #[arg(long)]
    pub report: Option<PathBuf>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if let Some(report) = args.report.as_deref() {
        super::ensure_distinct_inputs_and_outputs(
            &[("BED input", args.input.as_path())],
            &[("validation-report output", report)],
        )?;
    }
    let mode = if args.lenient {
        crate::io::bed::BedParseMode::Lenient
    } else {
        crate::io::bed::BedParseMode::Strict
    };
    let mut reader = crate::io::bed::read_bed12_with_mode(&args.input, mode)?;

    let mut transcript_count: u64 = 0;
    let mut exon_count: u64 = 0;
    let mut error_count: u64 = 0;
    let mut normalized_record_count: u64 = 0;
    let mut first_error = None;

    loop {
        let repairs_before = reader.warnings().len();
        let Some(result) = reader.next() else {
            break;
        };
        match result {
            Ok(transcript) => {
                transcript_count += 1;
                exon_count += transcript.exons.len() as u64;
                if reader.warnings().len() > repairs_before {
                    normalized_record_count += 1;
                }
            }
            Err(error) => {
                error_count += 1;
                let message = error.to_string();
                if first_error.is_none() {
                    first_error = Some(message.clone());
                }
                eprintln!("validation_error\t{}", escape_tsv(&message));
            }
        }
    }

    println!("records\t{transcript_count}");
    println!("exons\t{exon_count}");
    println!("repairs\t{}", reader.warnings().len());
    println!("normalized_records\t{normalized_record_count}");
    println!("errors\t{error_count}");
    for warning in reader.warnings() {
        eprintln!(
            "repair\t{}\t{}\t{}\t{}\t{}\t{}",
            escape_tsv(&warning.path.display().to_string()),
            warning.line,
            warning.field,
            escape_tsv(&warning.original),
            escape_tsv(&warning.replacement),
            escape_tsv(warning.reason)
        );
    }

    if let Some(report) = &args.report {
        crate::flow::artifact_manifest::atomic_write_with(report, |temporary| {
            write_report(
                temporary,
                mode,
                transcript_count,
                exon_count,
                reader.warnings().len() as u64,
                normalized_record_count,
                error_count,
            )
        })?;
    }

    if error_count > 0 {
        anyhow::bail!(
            "BED validation failed with {error_count} error(s); first error: {}",
            first_error.unwrap_or_else(|| "unknown validation error".to_owned())
        );
    }

    Ok(())
}

fn escape_tsv(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

fn write_report<W: Write>(
    output: &mut W,
    mode: crate::io::bed::BedParseMode,
    records: u64,
    exons: u64,
    repairs: u64,
    normalized_records: u64,
    errors: u64,
) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(output);
    writeln!(writer, "schema\ttrackcluster-bed-validation-v1")?;
    writeln!(
        writer,
        "mode\t{}",
        match mode {
            crate::io::bed::BedParseMode::Strict => "strict",
            crate::io::bed::BedParseMode::Lenient => "lenient",
        }
    )?;
    writeln!(writer, "records\t{records}")?;
    writeln!(writer, "exons\t{exons}")?;
    writeln!(writer, "repairs\t{repairs}")?;
    writeln!(writer, "normalized_records\t{normalized_records}")?;
    writeln!(writer, "errors\t{errors}")?;
    writer.flush()?;
    Ok(())
}
