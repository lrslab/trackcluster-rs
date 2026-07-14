use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Isoform BED (or reads_gene BED)
    #[arg(long = "isoform")]
    pub isoform: PathBuf,

    /// Reference BED
    #[arg(long = "reference")]
    pub reference: PathBuf,

    /// Output prefix
    #[arg(short = 'o', long = "out", default_value = "desc")]
    pub out: String,

    /// Junction fuzz/offset in bp (Python default: 10)
    #[arg(long = "offset-bp", default_value_t = 10, allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub offset_bp: u32,

    /// Optional strand-aware end-shift tags in _desc.txt; does not change class12; 0 disables
    #[arg(long = "end-shift-bp", default_value_t = 0, allow_hyphen_values = true, value_parser = crate::config::parse_base_pair_offset)]
    pub end_shift_bp: u32,

    /// Minimum fraction of isoform span overlapping a reference span for fusion detection (Python flow_fusion default: 0.1)
    #[arg(long = "fusion-fraction-read", default_value_t = 0.1, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub fusion_fraction_read: f64,

    /// Minimum fraction of reference span overlapping an isoform span for fusion detection (Python flow_fusion default: 0.1)
    #[arg(long = "fusion-fraction-ref", default_value_t = 0.1, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub fusion_fraction_ref: f64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let prefix = PathBuf::from(&args.out);
    let output_paths = crate::annotate::desc_output::DescOutputPaths::for_prefix(&prefix);
    let retired = crate::annotate::desc_output::retired_description_output_path(&prefix);
    super::ensure_distinct_inputs_and_outputs(
        &[
            ("isoform input", args.isoform.as_path()),
            ("reference input", args.reference.as_path()),
        ],
        &[
            ("description output", output_paths.desc.as_path()),
            ("class4 output", output_paths.class4.as_path()),
            ("fusion output", output_paths.fusion.as_path()),
            ("class12 output", output_paths.class12.as_path()),
            ("retired SQANTI output", retired.as_path()),
        ],
    )?;
    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.isoform)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::annotate::desc::try_describe(
        &isoforms,
        &refs,
        crate::annotate::desc::DescOpts {
            offset_bp: args.offset_bp,
            end_shift_bp: args.end_shift_bp,
            fusion_fraction_read: args.fusion_fraction_read,
            fusion_fraction_ref: args.fusion_fraction_ref,
        },
    )?;

    crate::annotate::desc_output::write_desc_outputs(&prefix, &result)?;

    Ok(())
}
