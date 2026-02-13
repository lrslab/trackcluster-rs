use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Sample manifest TSV with columns: sample, reads, optional group
    #[arg(long = "manifest")]
    pub manifest: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Isoform BED (from pooled cluster/flow)
    #[arg(short = 'i', long = "isoform")]
    pub isoform: PathBuf,

    /// Output file prefix (writes .isoform_usage.long.tsv and .isoform_counts.matrix.tsv)
    #[arg(short = 'o', long = "out")]
    pub out_prefix: PathBuf,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let outputs = crate::count::multi::run_count_multi_from_paths(
        &args.manifest,
        &args.reference,
        &args.isoform,
        &args.out_prefix,
    )?;

    eprintln!(
        "count-multi: long={:?} matrix={:?} group={:?}",
        outputs.long_tsv, outputs.matrix_tsv, outputs.group_tsv
    );

    Ok(())
}
