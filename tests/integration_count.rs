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
fn count_matches_golden_output() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let isoform = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_csv = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("count");
    let out_csv = out_dir.join("isoform_count.csv");

    let status = Command::new(exe)
        .args([
            "count",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "-o",
            out_csv.to_str().unwrap(),
        ])
        .status()
        .expect("run count");
    assert!(status.success());

    let produced = fs::read_to_string(out_csv).expect("read produced csv");
    let golden = fs::read_to_string(golden_csv).expect("read golden csv");
    assert_eq!(produced, golden);
}

#[test]
fn count_fractional_mode_preserves_split_counts() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let isoform = repo_path("tests/golden/clusterj/isoform.bed");

    let out_dir = fresh_temp_dir("count_fractional");
    let out_csv = out_dir.join("isoform_count.csv");

    let status = Command::new(exe)
        .args([
            "count",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "--assignment-mode",
            "fractional",
            "-o",
            out_csv.to_str().unwrap(),
        ])
        .status()
        .expect("run count");
    assert!(status.success());

    let produced = fs::read_to_string(out_csv).expect("read produced csv");
    assert_eq!(produced, "isoform_id,count\nref_a,0.5\nref_b,0.5\n");
}
