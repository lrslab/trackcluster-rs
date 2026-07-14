use std::process::Command;

fn assert_trackcluster_rejects(args: &[&str], expected: &str) {
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(args)
        .output()
        .expect("run trackcluster");
    assert!(
        !output.status.success(),
        "unexpected success for arguments {args:?}"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(expected),
        "stderr did not contain {expected:?} for arguments {args:?}:\n{stderr}"
    );
}

#[test]
fn cluster_rejects_nonfinite_out_of_range_and_zero_values() {
    let base = [
        "cluster",
        "--reads",
        "reads.bed",
        "--reference",
        "reference.bed",
    ];
    for (option, value, expected) in [
        ("--threads", "0", "worker threads must be at least 1"),
        ("--cutoff1", "NaN", "must be finite"),
        ("--cutoff1", "inf", "must be finite"),
        ("--cutoff2", "1.01", "must be within [0, 1]"),
        ("--intron-weight", "-0.1", "must be nonnegative"),
        ("--intron-weight", "inf", "must be finite"),
    ] {
        let mut args = base.to_vec();
        args.extend([option, value]);
        assert_trackcluster_rejects(&args, expected);
    }
}

#[test]
fn flow_and_preparedir_reject_invalid_fractions() {
    assert_trackcluster_rejects(
        &[
            "flow",
            "--reads",
            "reads.bed",
            "--reference",
            "reference.bed",
            "--output-root",
            "out",
            "--prefix",
            "sample",
            "--prepare-fraction-read",
            "-0.01",
        ],
        "must be within [0, 1]",
    );
    assert_trackcluster_rejects(
        &[
            "preparedir",
            "--reads",
            "reads.bed",
            "--reference",
            "reference.bed",
            "--output-root",
            "out",
            "--prefix",
            "sample",
            "--fraction-ref",
            "NaN",
        ],
        "must be finite",
    );
}

#[test]
fn clusterj_rejects_zero_support_and_negative_offsets() {
    let base = [
        "clusterj",
        "--reads",
        "reads.bed",
        "--reference",
        "reference.bed",
    ];
    for (option, value, expected) in [
        (
            "--junction-correction-min-support",
            "0",
            "minimum support must be at least 1",
        ),
        (
            "--sl-5prime-min-support",
            "0",
            "minimum support must be at least 1",
        ),
        (
            "--3prime-min-support",
            "0",
            "minimum support must be at least 1",
        ),
        (
            "--junction-correction-offset",
            "-1",
            "base-pair offset must be a nonnegative integer",
        ),
    ] {
        let mut args = base.to_vec();
        args.extend([option, value]);
        assert_trackcluster_rejects(&args, expected);
    }
}

#[test]
fn description_rejects_invalid_fusion_fraction() {
    assert_trackcluster_rejects(
        &[
            "desc",
            "--isoform",
            "isoform.bed",
            "--reference",
            "reference.bed",
            "--fusion-fraction-read",
            "2",
        ],
        "must be within [0, 1]",
    );
}

#[test]
fn batch_binary_rejects_zero_threads_before_io() {
    let output = Command::new(env!("CARGO_BIN_EXE_clusterj_batch"))
        .args([
            "--input-root",
            "input",
            "--output-root",
            "output",
            "--threads",
            "0",
        ])
        .output()
        .expect("run clusterj-batch");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("worker threads must be at least 1"));
}

#[test]
fn count_only_rejects_irrelevant_strict_gene_error_flag() {
    assert_trackcluster_rejects(
        &[
            "flow",
            "--count-only",
            "--strict-gene-errors",
            "--reference",
            "reference.bed",
            "--output-root",
            "out",
            "--prefix",
            "sample",
        ],
        "cannot be used with --count-only",
    );
}

#[test]
fn count_only_rejects_force_instead_of_ignoring_it() {
    assert_trackcluster_rejects(
        &[
            "flow",
            "--count-only",
            "--force",
            "--reference",
            "reference.bed",
            "--output-root",
            "out",
            "--prefix",
            "sample",
        ],
        "cannot be used with --count-only",
    );
}
