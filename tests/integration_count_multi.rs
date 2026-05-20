use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn repo_path(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn fresh_temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time went backwards")
        .as_nanos();
    let mut dir = std::env::temp_dir();
    dir.push(format!(
        "trackcluster_rs_{}_{}_{}",
        prefix,
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

#[test]
fn count_multi_outputs_expected_long_and_matrix_tables() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("count_multi");
    let prefix = "pooled";

    let flow_status = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run flow to generate pooled isoforms");
    assert!(flow_status.success());

    let isoform = out_dir.join(format!("{prefix}_isoform.bed"));
    let out_prefix = out_dir.join("recount");
    let count_status = Command::new(exe)
        .args([
            "count-multi",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("run count-multi");
    assert!(count_status.success());

    let long_tsv = out_dir.join("recount.isoform_usage.long.tsv");
    let matrix_tsv = out_dir.join("recount.isoform_counts.matrix.tsv");
    let group_tsv = out_dir.join("recount.isoform_usage.group.tsv");
    assert!(long_tsv.exists());
    assert!(matrix_tsv.exists());
    assert!(group_tsv.exists());

    let matrix = fs::read_to_string(matrix_tsv).expect("read matrix tsv");
    let mut lines = matrix.lines();
    assert_eq!(lines.next().unwrap(), "gene\tisoform_id\tS1\tS2");
    let mut counts = std::collections::HashMap::<String, (f64, f64)>::new();
    let mut seen_rows = 0usize;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        assert_eq!(fields[0], "GENEA");
        let c1: f64 = fields[2].parse().expect("parse S1 count");
        let c2: f64 = fields[3].parse().expect("parse S2 count");
        counts.insert(fields[1].to_owned(), (c1, c2));
        seen_rows += 1;
    }
    assert_eq!(seen_rows, 2);
    assert_eq!(counts.get("ref_a"), Some(&(1.0, 1.0)));
    assert_eq!(counts.get("ref_b"), Some(&(0.0, 0.0)));

    let long = fs::read_to_string(long_tsv).expect("read long tsv");
    let mut sums = std::collections::HashMap::<String, f64>::new();
    for (idx, line) in long.lines().enumerate() {
        if idx == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let key = format!("{}\t{}", fields[0], fields[2]); // gene + sample
        let proportion: f64 = fields[5].parse().expect("parse proportion");
        *sums.entry(key).or_insert(0.0) += proportion;
    }
    for total in sums.values() {
        assert!((*total - 1.0).abs() < 1e-9);
    }

    let fractional_out_prefix = out_dir.join("recount_fractional");
    let fractional_status = Command::new(exe)
        .args([
            "count-multi",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "--assignment-mode",
            "fractional",
            "-o",
            fractional_out_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("run fractional count-multi");
    assert!(fractional_status.success());

    let fractional_matrix =
        fs::read_to_string(out_dir.join("recount_fractional.isoform_counts.matrix.tsv"))
            .expect("read fractional matrix tsv");
    let mut lines = fractional_matrix.lines();
    assert_eq!(lines.next().unwrap(), "gene\tisoform_id\tS1\tS2");
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        let c1: f64 = fields[2].parse().expect("parse S1 fractional count");
        let c2: f64 = fields[3].parse().expect("parse S2 fractional count");
        assert!((c1 - 0.5).abs() < 1e-9);
        assert!((c2 - 0.5).abs() < 1e-9);
    }
}
