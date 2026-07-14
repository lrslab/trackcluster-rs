use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Input GFF3 or GTF annotation file
    #[arg(short = 'i', long = "gff", visible_alias = "input")]
    pub gff: PathBuf,

    /// Output TrackCluster bigGenePred-compatible BED path
    #[arg(
        short = 'o',
        long = "out",
        visible_alias = "output",
        default_value = "bigg.bed"
    )]
    pub out: PathBuf,

    /// GFF3 gene-feature attribute written as the BED gene ID
    #[arg(
        short = 'k',
        long = "key",
        visible_alias = "gene-key",
        default_value = "ID"
    )]
    pub key: String,

    /// Annotation attribute syntax
    #[arg(long = "input-format", value_enum, default_value_t = crate::io::gff::AnnotationFormat::Auto)]
    pub input_format: crate::io::gff::AnnotationFormat,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    super::ensure_distinct_input_output(&args.gff, &args.out, "annotation")?;
    let options = crate::io::gff::GffToBiggOptions {
        format: args.input_format,
        gene_key: args.key,
    };
    let transcripts = crate::io::gff::read_annotation_transcripts(&args.gff, &options)?;
    crate::flow::artifact_manifest::atomic_write_with(&args.out, |temporary| {
        crate::io::bed::write_bed12_to_writer(temporary, transcripts.iter()).map_err(Into::into)
    })?;
    eprintln!("gff2bigg: transcripts={}", transcripts.len());
    Ok(())
}
