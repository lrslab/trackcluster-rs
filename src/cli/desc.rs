use std::path::PathBuf;
use std::{fs::File, io::Write};

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
    #[arg(long = "offset-bp", default_value_t = 10)]
    pub offset_bp: u32,

    /// Strand-aware 5'/3' end-shift threshold in bp; 0 disables end-shift tagging (Python parity)
    #[arg(long = "end-shift-bp", default_value_t = 0)]
    pub end_shift_bp: u32,

    /// Minimum fraction of isoform span overlapping a reference span for fusion detection (Python flow_fusion default: 0.1)
    #[arg(long = "fusion-fraction-read", default_value_t = 0.1)]
    pub fusion_fraction_read: f64,

    /// Minimum fraction of reference span overlapping an isoform span for fusion detection (Python flow_fusion default: 0.1)
    #[arg(long = "fusion-fraction-ref", default_value_t = 0.1)]
    pub fusion_fraction_ref: f64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.isoform)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;

    let result = crate::annotate::desc::describe(
        &isoforms,
        &refs,
        crate::annotate::desc::DescOpts {
            offset_bp: args.offset_bp,
            end_shift_bp: args.end_shift_bp,
            fusion_fraction_read: args.fusion_fraction_read,
            fusion_fraction_ref: args.fusion_fraction_ref,
        },
    );

    let desc_path = PathBuf::from(format!("{}_desc.txt", args.out));
    let mut writer = std::io::BufWriter::new(File::create(desc_path)?);
    for row in &result.desc_rows {
        writeln!(
            &mut writer,
            "{}\t{}\t{}\t{}\t{}",
            row.isoform_id, row.ref_id, row.gene, row.miss, row.extra
        )?;
    }

    let class4_path = PathBuf::from(format!("{}_class4.txt", args.out));
    let mut writer = std::io::BufWriter::new(File::create(class4_path)?);
    for row in &result.class4_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    let fusion_path = PathBuf::from(format!("{}_fusion.txt", args.out));
    let mut writer = std::io::BufWriter::new(File::create(fusion_path)?);
    for row in &result.fusion_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.genes.join(";"))?;
    }

    let class12_path = PathBuf::from(format!("{}_class12.txt", args.out));
    let mut writer = std::io::BufWriter::new(File::create(class12_path)?);
    for row in &result.class12_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    Ok(())
}
