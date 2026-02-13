use std::io::Write;
use std::path::Path;

use anyhow::Context;

use crate::io::bed::{read_bed12, write_bed12_to_writer, BedError};
use crate::io::manifest::SampleRow;
use crate::sample::tagged_read_name;

pub fn write_pooled_reads(sample_rows: &[SampleRow], out_path: &Path) -> anyhow::Result<usize> {
    let file = std::fs::File::create(out_path)
        .with_context(|| format!("create pooled reads output {out_path:?}"))?;
    let mut writer = std::io::BufWriter::new(file);

    let mut written = 0usize;
    for row in sample_rows {
        let reader =
            read_bed12(&row.reads).with_context(|| format!("open reads {:?}", row.reads))?;
        for record in reader {
            let mut tx = record.with_context(|| format!("parse reads {:?}", row.reads))?;
            tx.name = tagged_read_name(&row.sample, &tx.name);
            write_bed12_to_writer(&mut writer, std::iter::once(&tx)).map_err(|source| {
                BedError::IoWrite {
                    path: out_path.to_path_buf(),
                    source,
                }
            })?;
            written += 1;
        }
    }
    writer
        .flush()
        .with_context(|| format!("flush pooled reads output {out_path:?}"))?;

    Ok(written)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::io::manifest::SampleRow;

    use super::*;

    fn fresh_temp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "trackcluster_rs_pool_{}_{}_{}",
            prefix,
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn writes_pooled_reads_with_sample_prefixes() {
        let dir = fresh_temp_dir("prefix");
        let s1 = dir.join("S1.reads.bed");
        let s2 = dir.join("S2.reads.bed");
        fs::write(&s1, "chr1\t0\t10\tr1\t0\t+\t0\t10\t0\t1\t10,\t0,\n").unwrap();
        fs::write(&s2, "chr1\t0\t10\tr2\t0\t+\t0\t10\t0\t1\t10,\t0,\n").unwrap();

        let pooled = dir.join("pooled.bed");
        let rows = vec![
            SampleRow {
                sample: "S1".to_owned(),
                group: Some("control".to_owned()),
                reads: s1.clone(),
            },
            SampleRow {
                sample: "S2".to_owned(),
                group: Some("treated".to_owned()),
                reads: s2.clone(),
            },
        ];
        let written = write_pooled_reads(&rows, &pooled).unwrap();
        assert_eq!(written, 2);

        let pooled_content = fs::read_to_string(&pooled).unwrap();
        assert!(pooled_content.contains("\tS1::r1\t"));
        assert!(pooled_content.contains("\tS2::r2\t"));
    }
}
