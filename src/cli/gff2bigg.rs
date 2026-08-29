use std::ffi::OsString;
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

    /// Handling of malformed annotation records: quarantine safe model-local failures, or fail
    #[arg(long = "invalid-record-policy", value_enum, default_value_t = crate::io::gff::InvalidAnnotationPolicy::Skip)]
    pub invalid_record_policy: crate::io::gff::InvalidAnnotationPolicy,

    /// Rejected-record audit TSV; defaults to <out>.rejected.tsv
    #[arg(long = "rejected-records")]
    pub rejected_records: Option<PathBuf>,
}

fn default_rejected_records_path(output: &std::path::Path) -> PathBuf {
    let mut value: OsString = output.as_os_str().to_os_string();
    value.push(".rejected.tsv");
    PathBuf::from(value)
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let rejected_path = args
        .rejected_records
        .unwrap_or_else(|| default_rejected_records_path(&args.out));
    super::ensure_distinct_inputs_and_outputs(
        &[("annotation input", args.gff.as_path())],
        &[
            ("BED output", args.out.as_path()),
            ("rejected-record output", rejected_path.as_path()),
        ],
    )?;
    let options = crate::io::gff::GffToBiggOptions {
        format: args.input_format,
        gene_key: args.key,
    };
    let result = crate::io::gff::read_annotation_transcripts_with_policy(
        &args.gff,
        &options,
        args.invalid_record_policy,
    )?;
    crate::flow::artifact_manifest::atomic_write_with(&rejected_path, |temporary| {
        crate::io::gff::write_rejected_annotations_tsv_to_writer(
            temporary,
            &result.rejected_records,
        )
        .map_err(Into::into)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&args.out, |temporary| {
        crate::io::bed::write_bed12_to_writer(temporary, result.transcripts.iter())
            .map_err(Into::into)
    })?;
    if !result.rejected_records.is_empty() {
        eprintln!(
            "gff2bigg: warning: excluded {} invalid record(s) affecting {} transcript model(s); details: {:?}",
            result.rejected_records.len(),
            result.rejected_transcripts,
            rejected_path
        );
    }
    eprintln!(
        "gff2bigg: transcripts={} rejected_records={} rejected_transcripts={}",
        result.transcripts.len(),
        result.rejected_records.len(),
        result.rejected_transcripts
    );
    Ok(())
}
