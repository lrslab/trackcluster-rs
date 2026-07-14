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
    let content = fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<String> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect();
    lines
}

#[test]
fn cluster_overlap_plain_bed12_reference_and_unmatched_reads_match_goldens() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/cluster_overlap/plain_reads.bed");
    let reference = repo_path("tests/fixtures/cluster_overlap/plain_ref.bed");
    let golden_isoforms = repo_path("tests/golden/cluster/plain_isoform.bed");
    let golden_mapping = repo_path("tests/golden/cluster/plain_isoform.read_to_isoform.tsv");
    let golden_unused = repo_path("tests/golden/cluster/plain_isoform.unused.bed");

    let out_dir = fresh_temp_dir("cluster_overlap_plain");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "cluster",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run overlap cluster");
    assert_success(&output, "overlap-cluster golden run");

    let produced_mapping = out_bed.with_extension("read_to_isoform.tsv");
    let produced_unused = out_bed.with_extension("unused.bed");

    assert_eq!(
        normalized_lines(&out_bed),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&produced_mapping),
        normalized_lines(&golden_mapping)
    );
    assert_eq!(
        normalized_lines(&produced_unused),
        normalized_lines(&golden_unused)
    );
}
