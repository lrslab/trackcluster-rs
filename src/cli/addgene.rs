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
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let annotated = crate::annotate::addgene::add_gene(
        &reads,
        &refs,
        crate::annotate::addgene::AddGeneOpts::default(),
    );
    crate::io::bed::write_bed12(&args.out, annotated.iter())?;
    Ok(())
}
