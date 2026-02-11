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
fn clusterj_matches_golden_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_mapping = repo_path("tests/golden/clusterj/isoform.read_to_isoform.tsv");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");

    let out_dir = fresh_temp_dir("clusterj");
    let out_bed = out_dir.join("isoform.bed");

    let status = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .status()
        .expect("run clusterj");
    assert!(status.success());

    let produced_isoforms = out_bed.clone();
    let produced_mapping = out_bed.with_extension("read_to_isoform.tsv");
    let produced_unused = out_bed.with_extension("unused.bed");

    assert_eq!(
        normalized_lines(&produced_isoforms),
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
