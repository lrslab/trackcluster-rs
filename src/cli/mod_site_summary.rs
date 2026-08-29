use std::ffi::OsString;
use std::path::{Path, PathBuf};

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Complete isoform modification site TSV; repeat to combine files.
    #[arg(long, required = true)]
    pub sites: Vec<PathBuf>,
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
    let output = append_suffix(&args.out, ".mod_site_summary.tsv");
    let labels = (1..=args.sites.len())
        .map(|index| format!("site table {index}"))
        .collect::<Vec<_>>();
    let inputs = labels
        .iter()
        .zip(&args.sites)
        .map(|(label, path)| (label.as_str(), path.as_path()))
        .collect::<Vec<_>>();
    crate::cli::ensure_distinct_inputs_and_outputs(
        &inputs,
        &[("site summary output", output.as_path())],
    )?;
    for path in &args.sites {
        crate::modification::generation::validate_current_flat(
            path,
            ".isoform_mod_sites.tsv",
            "isoform_mod_sites",
        )?;
    }

    let result = crate::modification::site_summary::summarize_site_files(&args.sites)?;
    crate::flow::artifact_manifest::atomic_write_with(&output, |writer| {
        crate::modification::site_summary::write_site_summary_tsv(writer, &result)
    })?;
    eprintln!(
        "mod-site-summary: inputs={} input_rows={} sites={}",
        args.sites.len(),
        result.input_rows(),
        result.site_count()
    );
    Ok(())
}
