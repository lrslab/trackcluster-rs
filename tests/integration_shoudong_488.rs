mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{assert_success, TestDir};

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
#[ignore]
fn clusterj_matches_python_flow_outputs_on_shoudong_488_subset() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let Some(root) = std::env::var_os("TRACKCLUSTER_SHOUDONG_488_ROOT") else {
        eprintln!("skipping: set TRACKCLUSTER_SHOUDONG_488_ROOT to the dataset root");
        return;
    };
    let root = PathBuf::from(root);
    let tracktest = root.join("tracktest");
    if !tracktest.exists() {
        eprintln!("skipping: missing dataset at {tracktest:?}");
        return;
    }

    // These folders are created by the Python flow in `test/flow_test.py` and contain:
    // - `<gene>_gff.bed` reference
    // - `<gene>_nano.bed` reads
    // - `<gene>_simple_coveragej.bed` expected isoform output
    // - `<gene>_unused.bed` expected rare-junction reads
    let genes = ["ATC24-1G10010", "ATC24-1G10020"];

    for gene in genes {
        let gene_dir = tracktest.join(gene);
        let reads = gene_dir.join(format!("{gene}_nano.bed"));
        let reference = gene_dir.join(format!("{gene}_gff.bed"));
        let expected_isoforms = gene_dir.join(format!("{gene}_simple_coveragej.bed"));
        let expected_unused = gene_dir.join(format!("{gene}_unused.bed"));

        for required in [&reads, &reference, &expected_isoforms, &expected_unused] {
            assert!(
                required.exists(),
                "missing required input {required:?} for gene {gene}"
            );
        }

        let out_dir = fresh_temp_dir(&format!("clusterj_shoudong_488_{gene}"));
        let out_bed = out_dir.join("isoform.bed");

        // Treat reads as having no SW-supported 5' signal while keeping ordinary merge behavior.
        let output = Command::new(exe)
            .args([
                "clusterj",
                "-s",
                reads.to_str().unwrap(),
                "-r",
                reference.to_str().unwrap(),
                "-o",
                out_bed.to_str().unwrap(),
                "--sw-score",
                "-1",
            ])
            .output()
            .expect("run clusterj");
        assert_success(&output, &format!("clusterj for {gene}"));

        let produced_unused = out_bed.with_extension("unused.bed");

        assert_eq!(
            normalized_lines(&out_bed),
            normalized_lines(&expected_isoforms),
            "isoforms mismatch for {gene}"
        );
        assert_eq!(
            normalized_lines(&produced_unused),
            normalized_lines(&expected_unused),
            "unused mismatch for {gene}"
        );
    }
}
