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
fn clusterj_batch_keeps_compat_behavior() {
    let exe = env!("CARGO_BIN_EXE_clusterj_batch");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");

    let out_dir = fresh_temp_dir("clusterj_batch");
    let prefix = "sample";

    let status = Command::new(exe)
        .args([
            "--prepare-reads",
            reads.to_str().unwrap(),
            "--prepare-reference",
            reference.to_str().unwrap(),
            "--prepare-prefix",
            prefix,
            "--input-root",
            out_dir.to_str().unwrap(),
            "--output-root",
            out_dir.to_str().unwrap(),
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run clusterj_batch");
    assert!(status.success());

    let gene_isoform = out_dir.join("GENEA/GENEA_simple_coveragej.bed");
    let gene_unused = out_dir.join("GENEA/GENEA_unused.bed");
    assert!(gene_isoform.exists());
    assert!(gene_unused.exists());
    assert_eq!(
        normalized_lines(&gene_isoform),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&gene_unused),
        normalized_lines(&golden_unused)
    );

    let summary = out_dir.join("clusterj_batch_summary.txt");
    assert!(summary.exists());
    let summary_text = fs::read_to_string(summary).unwrap();
    assert!(summary_text.contains("platform_preset\tgeneric"));
    assert!(summary_text.contains("junction_correction_offset\t10"));
    assert!(summary_text.contains("junction_correction_min_support\t2"));
    assert!(summary_text.contains("total_genes\t1"));
    assert!(summary_text.contains("errors\t0"));
}
