use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reference BED
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Reads BED
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Output BED
    #[arg(short = 'o', long = "out", default_value = "reads_gene.bed")]
    pub out: PathBuf,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    super::ensure_distinct_inputs_and_outputs(
        &[
            ("reads input", args.reads.as_path()),
            ("reference input", args.reference.as_path()),
        ],
        &[("annotated BED output", args.out.as_path())],
    )?;
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let annotated = crate::annotate::addgene::add_gene(
        &reads,
        &refs,
        crate::annotate::addgene::AddGeneOpts::default(),
    );
    crate::flow::artifact_manifest::atomic_write_with(&args.out, |temporary| {
        crate::io::bed::write_bed12_to_writer(temporary, annotated.iter()).map_err(Into::into)
    })?;
    Ok(())
}
