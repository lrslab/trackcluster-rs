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

#[test]
fn count_multi_outputs_expected_long_and_matrix_tables() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("count_multi");
    let prefix = "pooled";

    let flow_output = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
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
        .output()
        .expect("run flow to generate pooled isoforms");
    assert_success(&flow_output, "flow before multi-sample count");

    let isoform = out_dir.join(format!("{prefix}_isoform.bed"));
    let out_prefix = out_dir.join("recount.v1");
    let count_output = Command::new(exe)
        .args([
            "count-multi",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
            "--unique-assignment-junction-offset",
            "11",
        ])
        .output()
        .expect("run count-multi");
    assert_success(&count_output, "unique multi-sample count");

    let long_tsv = out_dir.join("recount.v1.isoform_usage.long.tsv");
    let matrix_tsv = out_dir.join("recount.v1.isoform_counts.matrix.tsv");
    let count_csv = out_dir.join("recount.v1.isoform_count.csv");
    let group_tsv = out_dir.join("recount.v1.isoform_usage.group.tsv");
    assert!(long_tsv.exists());
    assert!(matrix_tsv.exists());
    assert!(count_csv.exists());
    assert!(group_tsv.exists());
    let provenance =
        fs::read_to_string(out_dir.join("recount.v1.unique_assignment.provenance.tsv"))
            .expect("read unique-assignment provenance");
    assert!(provenance.contains("unique_assignment_junction_offset\t11\n"));
    assert!(!out_dir
        .join("recount.unique_assignment.provenance.tsv")
        .exists());

    let matrix = fs::read_to_string(matrix_tsv).expect("read matrix tsv");
    let mut lines = matrix.lines();
    assert_eq!(lines.next().unwrap(), "gene\tisoform_id\tS1\tS2");
    let mut counts = std::collections::HashMap::<String, (f64, f64)>::new();
    let mut seen_rows = 0usize;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        assert_eq!(fields[0], "GENEA");
        let c1: f64 = fields[2].parse().expect("parse S1 count");
        let c2: f64 = fields[3].parse().expect("parse S2 count");
        counts.insert(fields[1].to_owned(), (c1, c2));
        seen_rows += 1;
    }
    assert_eq!(seen_rows, 2);
    assert_eq!(counts.get("ref_a"), Some(&(1.0, 1.0)));
    assert_eq!(counts.get("ref_b"), Some(&(0.0, 0.0)));

    let mut count_reader = csv::Reader::from_path(count_csv).expect("read aggregate count csv");
    assert_eq!(
        count_reader
            .headers()
            .expect("count CSV header")
            .iter()
            .collect::<Vec<_>>(),
        ["gene", "isoform_id", "count"]
    );
    let mut aggregate_counts = std::collections::HashMap::<String, f64>::new();
    for record in count_reader.records() {
        let fields = record.expect("count CSV row");
        assert_eq!(fields.len(), 3);
        assert_eq!(&fields[0], "GENEA");
        aggregate_counts.insert(
            fields[1].to_owned(),
            fields[2].parse().expect("parse aggregate count"),
        );
    }
    for (isoform_id, (s1, s2)) in &counts {
        let aggregate = aggregate_counts
            .get(isoform_id)
            .copied()
            .expect("aggregate row for isoform");
        assert!((aggregate - (s1 + s2)).abs() < 1e-9);
    }

    let long = fs::read_to_string(long_tsv).expect("read long tsv");
    let mut sums = std::collections::HashMap::<String, f64>::new();
    for (idx, line) in long.lines().enumerate() {
        if idx == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let key = format!("{}\t{}", fields[0], fields[2]); // gene + sample
        let proportion: f64 = fields[5].parse().expect("parse proportion");
        *sums.entry(key).or_insert(0.0) += proportion;
    }
    for total in sums.values() {
        assert!((*total - 1.0).abs() < 1e-9);
    }

    let fractional_out_prefix = out_dir.join("recount_fractional");
    let fractional_output = Command::new(exe)
        .args([
            "count-multi",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-i",
            isoform.to_str().unwrap(),
            "--assignment-mode",
            "fractional",
            "-o",
            fractional_out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("run fractional count-multi");
    assert_success(&fractional_output, "fractional multi-sample count");

    let fractional_matrix =
        fs::read_to_string(out_dir.join("recount_fractional.isoform_counts.matrix.tsv"))
            .expect("read fractional matrix tsv");
    let mut lines = fractional_matrix.lines();
    assert_eq!(lines.next().unwrap(), "gene\tisoform_id\tS1\tS2");
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        let c1: f64 = fields[2].parse().expect("parse S1 fractional count");
        let c2: f64 = fields[3].parse().expect("parse S2 fractional count");
        assert!((c1 - 0.5).abs() < 1e-9);
        assert!((c2 - 0.5).abs() < 1e-9);
    }
}

#[test]
fn count_multi_unique_round_trips_a_raw_read_id_with_the_sample_delimiter() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let dir = fresh_temp_dir("count_multi_delimiter_read_id");
    let reads = dir.join("reads.bed");
    let manifest = dir.join("samples.tsv");
    let isoforms = dir.join("isoforms.bed");
    let mapping = dir.join("mapping.tsv");
    let out = dir.join("result");

    fs::write(
        &reads,
        "chr1\t100\t200\tS1::r1\t0\t+\t100\t200\t0\t1\t100,\t0,\n",
    )
    .expect("write delimiter-bearing read");
    fs::write(&manifest, "sample\tgroup\treads\nS1\tcontrol\treads.bed\n")
        .expect("write sample manifest");
    fs::write(
        &isoforms,
        concat!(
            "chr1\t100\t200\tiso1\t0\t+\t0\t0\t0\t1\t100,\t0,\tnone\tnone\tnone\t-1,\tisoform_anno\tGENEA\tnone\tnone\n",
            "chr1\t100\t210\tiso2\t0\t+\t0\t0\t0\t1\t110,\t0,\tnone\tnone\tnone\t-1,\tisoform_anno\tGENEA\tnone\tnone\n"
        ),
    )
    .expect("write candidate isoforms");
    fs::write(&mapping, "S1::S1::r1\tiso1\nS1::S1::r1\tiso2\n").expect("write pooled mapping");

    let output = Command::new(exe)
        .args([
            "count-multi",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference",
            isoforms.to_str().unwrap(),
            "--isoform",
            isoforms.to_str().unwrap(),
            "--read-to-isoform",
            mapping.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
        ])
        .output()
        .expect("run unique count-multi delimiter regression");
    assert_success(&output, "unique count-multi delimiter regression");

    let matrix = fs::read_to_string(dir.join("result.isoform_counts.matrix.tsv"))
        .expect("read delimiter regression matrix");
    assert!(matrix.lines().any(|line| line == "GENEA\tiso1\t1"));
    assert!(matrix.lines().any(|line| line == "GENEA\tiso2\t0"));
}

#[test]
fn count_multi_rejects_a_generated_output_that_aliases_an_input() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let dir = fresh_temp_dir("count_multi_output_alias");
    let out = dir.join("result");
    let isoforms = dir.join("result.isoform_count.csv");
    fs::copy(repo_path("tests/golden/clusterj/isoform.bed"), &isoforms).unwrap();
    let original = fs::read(&isoforms).unwrap();

    let output = Command::new(exe)
        .args(["count-multi", "--manifest"])
        .arg(repo_path("tests/fixtures/samples.tsv"))
        .arg("--reference")
        .arg(repo_path("tests/fixtures/ref.bed"))
        .arg("--isoform")
        .arg(&isoforms)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
    assert_eq!(fs::read(&isoforms).unwrap(), original);
}
