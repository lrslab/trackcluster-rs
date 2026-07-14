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

fn bed_names(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .expect("read BED")
        .lines()
        .filter_map(|line| line.split('\t').nth(3).map(ToOwned::to_owned))
        .collect()
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

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run clusterj");
    assert_success(&output, "clusterj golden run");

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

#[test]
fn single_gene_cluster_commands_reject_an_output_that_aliases_an_input() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let root = fresh_temp_dir("single_gene_output_alias");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::copy(repo_path("tests/fixtures/reads.bed"), &reads).unwrap();
    fs::copy(repo_path("tests/fixtures/ref.bed"), &reference).unwrap();
    let original = fs::read(&reads).unwrap();

    for command in ["clusterj", "cluster"] {
        let output = Command::new(exe)
            .args([command, "--reads"])
            .arg(&reads)
            .arg("--reference")
            .arg(&reference)
            .arg("--out")
            .arg(&reads)
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "{command} accepted an input alias"
        );
        assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
        assert_eq!(fs::read(&reads).unwrap(), original);
    }
}

#[test]
fn clusterj_plain_bed_reference_is_protected_and_unmatched_reads_are_auditable() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/cluster_overlap/plain_reads.bed");
    let reference = repo_path("tests/fixtures/cluster_overlap/plain_ref.bed");
    let out_dir = fresh_temp_dir("clusterj_plain_bed");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run clusterj with plain BED12 inputs");

    assert!(
        output.status.success(),
        "clusterj failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(bed_names(&out_bed), vec!["ref_plain"]);
    assert_eq!(
        normalized_lines(&out_bed.with_extension("read_to_isoform.tsv")),
        vec!["read_match\tref_plain"]
    );
    assert_eq!(
        bed_names(&out_bed.with_extension("unused.bed")),
        vec!["read_wrong_strand", "read_disjoint", "read_wrong_chrom"]
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("input_reads=4"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("represented_reads=1"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("mapping_rows=1"),
        "missing summary: {stderr}"
    );
    assert!(stderr.contains("rare_reads=0"), "missing summary: {stderr}");
    assert!(
        stderr.contains("unmatched_reads=3"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("unused_reads=3"),
        "missing summary: {stderr}"
    );
}

#[test]
fn clusterj_rejects_correction_that_would_create_an_empty_exon() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("clusterj_invalid_snap_input");
    let reads = input_dir.join("reads.bed");
    let reference = input_dir.join("reference.bed");
    fs::write(
        &reference,
        "chr1\t80\t220\tref\t0\t+\t80\t220\t0\t2\t10,20,\t0,120,\n",
    )
    .expect("write correction reference");
    fs::write(
        &reads,
        "chr1\t100\t210\tread\t100\t+\t100\t210\t0\t2\t1,9,\t0,101,\n",
    )
    .expect("write correction read");
    let out_dir = fresh_temp_dir("clusterj_invalid_snap_output");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
            "--junction-correction-offset",
            "15",
            "--junction-correction-min-support",
            "5",
            "--sw-score",
            "11",
            "--sl-5prime-min-support",
            "1",
            "--sl-same-junction-5prime-offset",
            "0",
        ])
        .output()
        .expect("run clusterj invalid-snap regression");
    assert_success(&output, "clusterj invalid-snap regression");

    let unused = out_bed.with_extension("unused.bed");
    assert_eq!(bed_names(&out_bed), vec!["ref"]);
    assert_eq!(bed_names(&unused), vec!["read"]);
    assert!(String::from_utf8_lossy(&output.stderr).contains("rare_reads=1"));

    for path in [&out_bed, &unused] {
        let validation = Command::new(exe)
            .args(["validate-bed", "--input", path.to_str().unwrap()])
            .output()
            .expect("strictly validate clusterj regression output");
        assert_success(
            &validation,
            "strict validation of clusterj regression output",
        );
    }
}

#[test]
fn single_gene_cluster_commands_skip_only_bad_read_tracks() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("single_gene_bad_read_input");
    let reads = input_dir.join("dirty.bed");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    fs::write(&reads, format!("not-a-bed-record\n{good}")).unwrap();
    let reference = repo_path("tests/fixtures/ref.bed");

    for command in ["clusterj", "cluster"] {
        let out_dir = fresh_temp_dir(&format!("single_gene_{command}_bad_read"));
        let out_bed = out_dir.join("isoform.bed");
        let output = Command::new(exe)
            .args([
                command,
                "-s",
                reads.to_str().unwrap(),
                "-r",
                reference.to_str().unwrap(),
                "-o",
                out_bed.to_str().unwrap(),
            ])
            .output()
            .expect("run single-gene clustering");
        assert!(
            output.status.success(),
            "{command} stopped on one bad read: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let mapping = fs::read_to_string(out_bed.with_extension("read_to_isoform.tsv")).unwrap();
        assert!(mapping.contains("read_trunc\t"), "{mapping}");
        let rejected = fs::read_to_string(out_bed.with_extension("rejected_reads.tsv")).unwrap();
        assert_eq!(rejected.lines().count(), 2, "{rejected}");
        assert!(String::from_utf8_lossy(&output.stderr).contains("excluded 1 malformed read"));
    }

    let strict_out = input_dir.join("strict.bed");
    let strict = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--invalid-read-policy",
            "fail",
        ])
        .output()
        .expect("run strict single-gene clustering");
    assert!(!strict.status.success());
    assert!(!strict_out.exists());
}
