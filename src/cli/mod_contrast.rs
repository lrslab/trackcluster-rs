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
    crate::modification::generation::ensure_standalone_prefix(&args.out)?;
    let output = append_suffix(&args.out, ".isoform_mod_contrasts.tsv");
    let count = run_with_output_path(&args.design, &args.contrasts, &output)?;
    eprintln!("mod-contrast: contrasts={count}");
    Ok(())
}

pub(crate) fn run_with_output_path(
    design: &Path,
    contrasts: &Path,
    output: &Path,
) -> anyhow::Result<usize> {
    crate::modification::generation::validate_current_flat(
        design,
        ".isoform_mod_design.tsv",
        "isoform_mod_design",
    )?;
    crate::cli::ensure_distinct_inputs_and_outputs(
        &[
            ("modification design", design),
            ("contrast specification", contrasts),
        ],
        &[("contrast output", output)],
    )?;
    let specs = crate::modification::contrast::read_contrast_specs(contrasts)?;
    let rows = crate::modification::contrast::calculate_contrasts(design, &specs)?;
    crate::flow::artifact_manifest::atomic_write_with(output, |writer| {
        crate::modification::contrast::write_contrasts_tsv(writer, &rows)
    })?;
    Ok(rows.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn standalone_command_rejects_a_flow_managed_prefix_before_reading_inputs() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("test clock is before Unix epoch")
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "trackcluster-mod-contrast-managed-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir(&root).expect("create temporary directory");
        let prefix = root.join("result");
        crate::modification::generation::ensure_managed(&prefix)
            .expect("mark output prefix as flow-managed");

        let error = run(Args {
            design: root.join("missing-design.tsv"),
            contrasts: root.join("missing-contrasts.tsv"),
            out: prefix,
        })
        .expect_err("standalone contrast must not replace a managed flat output");

        assert!(error.to_string().contains("flow-managed"), "{error:#}");
        fs::remove_dir_all(root).expect("remove temporary directory");
    }
}
