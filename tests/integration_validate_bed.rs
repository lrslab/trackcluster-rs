mod common;

use std::fs;
use std::path::Path;
use std::process::Command;

use common::TestDir;

fn temp_dir(label: &str) -> TestDir {
    TestDir::new(&format!("validate-bed-{label}"))
}

fn run_validate(path: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["validate-bed", "--input"])
        .arg(path)
        .output()
        .unwrap()
}

#[test]
fn strict_validate_bed_rejects_each_malformed_structure() {
    let cases = [
        (
            "invalid-score",
            "chr1\t100\t200\ttx\tbad\t+\t100\t200\t0\t1\t100,\t0,\n",
        ),
        (
            "range-score",
            "chr1\t100\t200\ttx\t1001\t+\t100\t200\t0\t1\t100,\t0,\n",
        ),
        (
            "empty-token",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t50,,30,\t0,70,\n",
        ),
        (
            "outside",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t1\t101,\t0,\n",
        ),
        (
            "zero-exon",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t0,100,\t0,0,\n",
        ),
        (
            "overlap",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t70,60,\t0,40,\n",
        ),
        (
            "ordering",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t40,40,\t60,0,\n",
        ),
        (
            "count",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t100,\t0,\n",
        ),
        (
            "span",
            "chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t1\t80,\t10,\n",
        ),
        (
            "thick",
            "chr1\t100\t200\ttx\t0\t+\t90\t200\t0\t1\t100,\t0,\n",
        ),
    ];

    let root = temp_dir("strict");
    for (label, contents) in cases {
        let input = root.join(format!("{label}.bed"));
        fs::write(&input, contents).unwrap();
        let output = run_validate(&input);
        assert!(
            !output.status.success(),
            "{label} unexpectedly passed: {}",
            String::from_utf8_lossy(&output.stdout)
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("validation_error"),
            "{label} did not emit a structured error: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn lenient_validate_bed_reports_repairs_and_writes_report() {
    let root = temp_dir("lenient");
    let input = root.join("legacy.bed");
    let report = root.join("report.tsv");
    fs::write(
        &input,
        "chr1\t100\t200\ttx\tbad\t+\t100\t200\t0\t2\t50,,50,\t0,70,\n",
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["validate-bed", "--input"])
        .arg(&input)
        .arg("--lenient")
        .arg("--report")
        .arg(&report)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(stderr.matches("repair\t").count(), 3, "{stderr}");
    let report = fs::read_to_string(report).unwrap();
    assert!(report.contains("schema\ttrackcluster-bed-validation-v1"));
    assert!(report.contains("mode\tlenient"));
    assert!(report.contains("repairs\t3"));
    assert!(report.contains("normalized_records\t1"));
    assert!(report.contains("errors\t0"));
}
