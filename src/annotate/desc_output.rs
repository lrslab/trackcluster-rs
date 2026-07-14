use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;

use super::desc::DescResult;

pub const DESC_SCHEMA_VERSION: &str = "trackcluster-description-v2";
const RETIRED_SQANTI_CLASSIFICATION_SUFFIX: &str = "_sqanti_structural_category.tsv";

#[derive(Clone, Debug)]
pub struct DescOutputPaths {
    pub desc: PathBuf,
    pub class4: PathBuf,
    pub fusion: PathBuf,
    pub class12: PathBuf,
}

impl DescOutputPaths {
    pub fn for_prefix(prefix: &Path) -> Self {
        Self {
            desc: append_suffix(prefix, "_desc.txt"),
            class4: append_suffix(prefix, "_class4.txt"),
            fusion: append_suffix(prefix, "_fusion.txt"),
            class12: append_suffix(prefix, "_class12.txt"),
        }
    }
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

pub fn retired_description_output_path(prefix: &Path) -> PathBuf {
    append_suffix(prefix, RETIRED_SQANTI_CLASSIFICATION_SUFFIX)
}

fn finish<W: Write>(mut writer: BufWriter<W>) -> anyhow::Result<()> {
    writer.flush()?;
    Ok(())
}

pub(crate) fn remove_retired_description_outputs(prefix: &Path) -> anyhow::Result<()> {
    let path = retired_description_output_path(prefix);
    match std::fs::remove_file(&path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("remove retired output {path:?}")),
    }
}

pub fn write_desc_outputs(prefix: &Path, result: &DescResult) -> anyhow::Result<DescOutputPaths> {
    remove_retired_description_outputs(prefix)?;
    let paths = DescOutputPaths::for_prefix(prefix);

    crate::flow::artifact_manifest::atomic_write_with(&paths.desc, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(writer, "#schema\t{DESC_SCHEMA_VERSION}\tdesc")?;
        writeln!(
            writer,
            "isoform_id\treference_id\tgene_id\tmissing_features\textra_features"
        )?;
        for row in &result.desc_rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                row.isoform_id, row.ref_id, row.gene, row.miss, row.extra
            )?;
        }
        finish(writer)
    })?;

    crate::flow::artifact_manifest::atomic_write_with(&paths.class4, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(writer, "#schema\t{DESC_SCHEMA_VERSION}\tclass4")?;
        writeln!(writer, "isoform_id\tclass")?;
        for row in &result.class4_rows {
            writeln!(writer, "{}\t{}", row.isoform_id, row.class)?;
        }
        finish(writer)
    })?;

    crate::flow::artifact_manifest::atomic_write_with(&paths.fusion, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(writer, "#schema\t{DESC_SCHEMA_VERSION}\tfusion")?;
        writeln!(writer, "isoform_id\tgene_ids")?;
        for row in &result.fusion_rows {
            writeln!(writer, "{}\t{}", row.isoform_id, row.genes.join(";"))?;
        }
        finish(writer)
    })?;

    crate::flow::artifact_manifest::atomic_write_with(&paths.class12, |temporary| {
        let mut writer = BufWriter::new(temporary);
        writeln!(writer, "#schema\t{DESC_SCHEMA_VERSION}\tclass12")?;
        writeln!(writer, "isoform_id\tclass")?;
        for row in &result.class12_rows {
            writeln!(writer, "{}\t{}", row.isoform_id, row.class)?;
        }
        finish(writer)
    })?;

    Ok(paths)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    #[test]
    fn outputs_have_schema_and_column_headers() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let prefix = std::env::temp_dir().join(format!(
            "trackcluster-desc-output-{}-{nonce}",
            std::process::id()
        ));
        let retired = append_suffix(&prefix, RETIRED_SQANTI_CLASSIFICATION_SUFFIX);
        fs::write(&retired, "stale classification\n").unwrap();
        let paths = write_desc_outputs(&prefix, &DescResult::default()).unwrap();
        assert!(!retired.exists());
        for path in [&paths.desc, &paths.class4, &paths.fusion, &paths.class12] {
            let text = fs::read_to_string(path).unwrap();
            assert!(text.starts_with("#schema\ttrackcluster-description-v2"));
            assert!(text.lines().count() >= 2);
            fs::remove_file(path).unwrap();
        }
    }
}
