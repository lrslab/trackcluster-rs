use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Input BED12/bigGenePred transcript catalog
    #[arg(short, long)]
    pub input: PathBuf,

    /// Output GTF 2.2 path
    #[arg(long)]
    pub gtf: Option<PathBuf>,

    /// Output GFF3 path
    #[arg(long)]
    pub gff3: Option<PathBuf>,

    /// Output SQANTI3 input-audit TSV path
    #[arg(long)]
    pub sqanti_input: Option<PathBuf>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.gtf.is_none() && args.gff3.is_none() && args.sqanti_input.is_none() {
        anyhow::bail!("export requires at least one of --gtf, --gff3, or --sqanti-input");
    }
    let mut outputs = Vec::new();
    if let Some(path) = args.gtf.as_deref() {
        outputs.push(("GTF output", path));
    }
    if let Some(path) = args.gff3.as_deref() {
        outputs.push(("GFF3 output", path));
    }
    if let Some(path) = args.sqanti_input.as_deref() {
        outputs.push(("SQANTI input-audit output", path));
    }
    super::ensure_distinct_inputs_and_outputs(&[("BED input", args.input.as_path())], &outputs)?;
    let transcripts = crate::io::bed::read_bed12(&args.input)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    if let Some(path) = args.gtf.as_deref() {
        crate::flow::artifact_manifest::atomic_write_with(path, |temporary| {
            crate::io::interchange::write_gtf_to_writer(temporary, &transcripts)
        })?;
    }
    if let Some(path) = args.gff3.as_deref() {
        crate::flow::artifact_manifest::atomic_write_with(path, |temporary| {
            crate::io::interchange::write_gff3_to_writer(temporary, &transcripts)
        })?;
    }
    if let Some(path) = args.sqanti_input.as_deref() {
        crate::flow::artifact_manifest::atomic_write_with(path, |temporary| {
            crate::io::interchange::write_sqanti3_input_table_to_writer(temporary, &transcripts)
        })?;
    }
    Ok(())
}
