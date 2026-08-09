use std::ffi::OsString;
use std::path::{Path, PathBuf};

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Comparison-ready isoform modification design TSV.
    #[arg(long)]
    pub design: PathBuf,
    /// Explicit contrast specification TSV.
    #[arg(long)]
    pub contrasts: PathBuf,
    /// Output prefix.
    #[arg(short, long)]
    pub out: PathBuf,
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let output = append_suffix(&args.out, ".isoform_mod_contrasts.tsv");
    crate::cli::ensure_distinct_inputs_and_outputs(
        &[
            ("modification design", args.design.as_path()),
            ("contrast specification", args.contrasts.as_path()),
        ],
        &[("contrast output", output.as_path())],
    )?;
    let specs = crate::modification::contrast::read_contrast_specs(&args.contrasts)?;
    let rows = crate::modification::contrast::calculate_contrasts(&args.design, &specs)?;
    crate::flow::artifact_manifest::atomic_write_with(&output, |writer| {
        crate::modification::contrast::write_contrasts_tsv(writer, &rows)
    })?;
    eprintln!("mod-contrast: contrasts={}", rows.len());
    Ok(())
}
