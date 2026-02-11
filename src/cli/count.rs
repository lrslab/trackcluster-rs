use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED (sorted recommended)
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Isoform BED (from cluster/clusterj)
    #[arg(short = 'i', long = "isoform")]
    pub isoform: PathBuf,

    /// Output CSV
    #[arg(short = 'o', long = "out", default_value = "isoform_count.csv")]
    pub out: PathBuf,
}

pub fn run(_args: Args) -> anyhow::Result<()> {
    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.isoform)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;

    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;

    let records = crate::count::count_by_subreads(&isoforms, &refs);
    crate::count::write_counts_csv(&_args.out, &records)?;

    Ok(())
}
