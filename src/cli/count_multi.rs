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

    /// Optional read-to-isoform TSV mapping (fast path; from pooled flow output)
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: Option<PathBuf>,

    /// Output file prefix (writes .isoform_usage.long.tsv and .isoform_counts.matrix.tsv)
    #[arg(short = 'o', long = "out")]
    pub out_prefix: PathBuf,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let outputs = if let Some(mapping_path) = args.read_to_isoform.as_ref() {
        let sample_rows = crate::io::manifest::read_manifest_tsv(&args.manifest)?;
        let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.isoform)?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
        )?;
        let pairs = crate::count::read_read_to_isoform_tsv(mapping_path)?;

        crate::count::multi::run_count_multi_from_read_to_isoform(
            &sample_rows,
            &isoforms,
            &pairs,
            &args.out_prefix,
        )?
    } else {
        crate::count::multi::run_count_multi_from_paths(
            &args.manifest,
            &args.reference,
            &args.isoform,
            &args.out_prefix,
        )?
    };

    eprintln!(
        "count-multi: long={:?} matrix={:?} group={:?}",
        outputs.long_tsv, outputs.matrix_tsv, outputs.group_tsv
    );

    Ok(())
}
