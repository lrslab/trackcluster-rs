mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{assert_success, TestDir};

fn repo_path(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn fresh_temp_dir(prefix: &str) -> TestDir {
    TestDir::new(prefix)
}

fn normalized_lines(path: &Path) -> Vec<String> {
    let lines: Vec<String> = fs::read_to_string(path)
        .expect("read text file")
        .lines()
        .map(|line| line.to_owned())
        .collect();
    lines
}

#[test]
fn count_output_root_mode_matches_golden_output() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_csv = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("count_output_root");
    let prefix = "sample";

    let output = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
        ])
        .output()
        .expect("run flow");
    assert_success(&output, "flow before output-root count");

    let output = Command::new(exe)
        .args([
            "count",
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
        ])
        .output()
        .expect("run output-root count");
    assert_success(&output, "output-root count");

    let count_csv = out_dir.join(format!("{prefix}_isoform_count.csv"));
    assert_eq!(normalized_lines(&count_csv), normalized_lines(&golden_csv));

    assert!(out_dir
        .join(format!("{prefix}_read_to_isoform.unique.tsv"))
        .exists());
}

#[test]
fn legacy_count_matches_golden_output() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let isoform = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_csv = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("count");
    let out_csv = out_dir.join("isoform_count.csv");

    let output = Command::new(exe)
        .args([
            "count",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "--unique-assignment-junction-offset",
            "8",
            "--out",
            out_csv.to_str().unwrap(),
        ])
        .output()
        .expect("run count");
    assert_success(&output, "legacy count");

    let produced = fs::read_to_string(&out_csv).expect("read produced csv");
    let golden = fs::read_to_string(golden_csv).expect("read golden csv");
    assert_eq!(produced, golden);
    let provenance = fs::read_to_string(out_csv.with_extension("provenance.tsv"))
        .expect("read unique-assignment provenance");
    assert!(provenance.contains("unique_assignment_junction_offset\t8\n"));
}

#[test]
fn legacy_count_rejects_an_output_that_aliases_an_input() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let root = fresh_temp_dir("count_output_alias");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    let isoform = root.join("isoform.bed");
    fs::copy(repo_path("tests/fixtures/reads.bed"), &reads).unwrap();
    fs::copy(repo_path("tests/fixtures/ref.bed"), &reference).unwrap();
    fs::copy(repo_path("tests/golden/clusterj/isoform.bed"), &isoform).unwrap();
    let original = fs::read(&isoform).unwrap();

    let output = Command::new(exe)
        .args(["count", "--reads"])
        .arg(&reads)
        .arg("--reference")
        .arg(&reference)
        .arg("--isoform")
        .arg(&isoform)
        .arg("--out")
        .arg(&isoform)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
    assert_eq!(fs::read(&isoform).unwrap(), original);
}

#[test]
fn count_fractional_mode_preserves_split_counts() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let isoform = repo_path("tests/golden/clusterj/isoform.bed");

    let out_dir = fresh_temp_dir("count_fractional");
    let out_csv = out_dir.join("isoform_count.csv");
    fs::write(out_csv.with_extension("provenance.tsv"), "stale\n")
        .expect("write stale unique-assignment provenance");

    let output = Command::new(exe)
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
            "--out",
            out_csv.to_str().unwrap(),
        ])
        .output()
        .expect("run count");
    assert_success(&output, "fractional count");

    let produced = fs::read_to_string(out_csv).expect("read produced csv");
    assert_eq!(
        produced,
        "gene,isoform_id,count\nGENEA,ref_a,0.5\nGENEA,ref_b,0.5\n"
    );
    assert!(!out_dir.join("isoform_count.provenance.tsv").exists());
}
