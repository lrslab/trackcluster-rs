use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED (sorted recommended)
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Output isoform BED
    #[arg(short = 'o', long = "out", default_value = "isoform.bed")]
    pub out: PathBuf,

    /// Number of threads
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reads)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::cluster::cluster_overlap::cluster(&reads, Some(&refs), args.threads);

    crate::cluster::output::write_isoforms_bed(&args.out, &result.isoforms)?;

    let mapping_path = args.out.with_extension("read_to_isoform.tsv");
    crate::cluster::output::write_read_to_isoform_tsv(mapping_path, &result.read_to_isoform)?;

    let unused_path = args.out.with_extension("unused.bed");
    crate::cluster::output::write_isoforms_bed(unused_path, &result.unused)?;

    Ok(())
}
