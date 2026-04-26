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

fn normalized_lines(path: &Path) -> Vec<String> {
    let content = fs::read_to_string(path).unwrap_or_default();
    let mut lines: Vec<String> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect();
    lines.sort();
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

    let status = Command::new(exe)
        .args([
            "cluster",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .status()
        .expect("run overlap cluster");
    assert!(status.success());

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
