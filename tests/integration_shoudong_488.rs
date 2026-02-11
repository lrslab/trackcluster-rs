use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

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
#[ignore]
fn clusterj_matches_python_flow_outputs_on_shoudong_488_subset() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let root = std::env::var("TRACKCLUSTER_SHOUDONG_488_ROOT")
        .unwrap_or_else(|_| "/t1/shoudong_488/test".to_owned());
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

        // Match Python `flow_clusterj_all_gene_novel(..., sw_score=-1)` behavior: disable collapsing.
        let status = Command::new(exe)
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
            .status()
            .expect("run clusterj");
        assert!(status.success(), "clusterj failed for {gene}");

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
