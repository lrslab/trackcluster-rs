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
fn flow_runs_end_to_end_and_matches_goldens() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");
    let golden_count = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("flow");
    let prefix = "sample";

    let status = Command::new(exe)
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
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run flow");
    assert!(status.success());

    let isoform_out = out_dir.join(format!("{prefix}_isoform.bed"));
    let unused_out = out_dir.join(format!("{prefix}_unused.bed"));
    let count_out = out_dir.join(format!("{prefix}_isoform_count.csv"));

    assert_eq!(
        normalized_lines(&isoform_out),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&unused_out),
        normalized_lines(&golden_unused)
    );
    assert_eq!(
        normalized_lines(&count_out),
        normalized_lines(&golden_count)
    );

    assert!(out_dir.join("GENEA/GENEA_simple_coveragej.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_unused.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_read_to_isoform.tsv").exists());

    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class4.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_fusion.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class12.txt")).exists());
}
