mod common;

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{assert_success, TestDir};

const CATALOG: &str = "chr1\t100\t400\tshort\t0\t+\t100\t400\t0\t2\t100,100\t0,200\nchr1\t100\t500\tlong\t0\t+\t100\t500\t0\t2\t100,200\t0,200\nchr2\t100\t200\tzero\t0\t+\t100\t200\t0\t1\t100\t0\n";
const SHORT_READ: &str = "chr1\t100\t400\tr_short\t0\t+\t100\t400\t0\t2\t100,100\t0,200\n";
const LONG_READ: &str = "chr1\t100\t500\tr_long\t0\t+\t100\t500\t0\t2\t100,200\t0,200\n";
const FUZZY_READ: &str = "chr1\t100\t390\tr_fuzzy\t0\t+\t100\t390\t0\t2\t105,80\t0,210\n";
const UNASSIGNED_READ: &str = "chr9\t100\t200\tr_unknown\t0\t+\t100\t200\t0\t1\t100\t0\n";

fn write_inputs(root: &Path) -> (PathBuf, PathBuf) {
    let reads = root.join("reads.bed");
    let catalog = root.join("catalog.bed");
    fs::write(&catalog, CATALOG).unwrap();
    fs::write(
        &reads,
        format!("{SHORT_READ}{LONG_READ}{FUZZY_READ}{UNASSIGNED_READ}{SHORT_READ}"),
    )
    .unwrap();
    // This malformed neighboring discovery mapping must never be consulted.
    fs::write(
        catalog.with_extension("read_to_isoform.tsv"),
        "stale\tmissing_isoform\n",
    )
    .unwrap();
    (reads, catalog)
}

fn count_command(reads: &Path, catalog: &Path, out: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    command
        .args(["count", "--assign-against-catalog", "--reads"])
        .arg(reads)
        .arg("--isoform")
        .arg(catalog)
        .arg("--out")
        .arg(out);
    command
}

fn counts(path: &Path) -> BTreeMap<String, f64> {
    csv::Reader::from_path(path)
        .unwrap()
        .records()
        .map(|row| {
            let row = row.unwrap();
            (row[1].to_owned(), row[2].parse().unwrap())
        })
        .collect()
}

#[test]
fn count_directly_from_external_catalog_without_discovery_or_reference() {
    let root = TestDir::new("catalog_count");
    let (reads, catalog) = write_inputs(&root);
    let out = root.join("nested/quant.v1.csv");
    let output = count_command(&reads, &catalog, &out)
        .args(["--unique-assignment-junction-offset", "8"])
        .output()
        .unwrap();
    assert_success(&output, "direct catalog count");
    assert_eq!(
        counts(&out),
        BTreeMap::from([
            ("short".into(), 2.0),
            ("long".into(), 1.0),
            ("zero".into(), 0.0)
        ])
    );
    let pairs =
        trackcluster_rs::count::read_read_to_isoform_tsv(out.with_extension("read_to_isoform.tsv"))
            .unwrap();
    assert_eq!(
        pairs,
        vec![
            ("r_fuzzy".into(), "short".into()),
            ("r_long".into(), "long".into()),
            ("r_short".into(), "short".into())
        ]
    );
    assert_eq!(
        fs::read_to_string(out.with_extension("unassigned_reads.tsv")).unwrap(),
        "read_id\tgene\treason\nr_unknown\tnone\tno_candidate_locus\n"
    );
    assert_eq!(
        fs::read_to_string(out.with_extension("assignment_stats.tsv")).unwrap(),
        "metric\tcount\ninput_records\t5\ntotal_reads\t4\nassigned_reads\t3\nunassigned_reads\t1\n"
    );
    let provenance = fs::read_to_string(out.with_extension("provenance.tsv")).unwrap();
    assert!(provenance.contains("unique_assignment_junction_offset\t8\n"));
    assert!(provenance.contains("assignment_source\tfixed_catalog\n"));
    assert_eq!(fs::read_to_string(catalog).unwrap(), CATALOG);
    assert_eq!(fs::read_dir(out.parent().unwrap()).unwrap().count(), 5);
}

#[test]
fn count_multi_keeps_sample_identities_zero_columns_and_group_totals() {
    let root = TestDir::new("catalog_count_multi");
    let (_, catalog) = write_inputs(&root);
    // Same original read name as S1, but a different alignment in S2.
    fs::write(root.join("s2.bed"), LONG_READ.replace("r_long", "r_short")).unwrap();
    fs::write(root.join("empty.bed"), "").unwrap();
    let manifest = root.join("samples.tsv");
    fs::write(
        &manifest,
        "sample\treads\tgroup\nS1\treads.bed\tcontrol\nS2\ts2.bed\tcase\nS3\tempty.bed\tcase\n",
    )
    .unwrap();
    let prefix = root.join("nested/quant.v1");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["count-multi", "--assign-against-catalog", "--manifest"])
        .arg(&manifest)
        .arg("--isoform")
        .arg(&catalog)
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert_success(&output, "direct multi-sample catalog count");
    assert_eq!(fs::read_to_string(prefix.with_file_name("quant.v1.isoform_counts.matrix.tsv")).unwrap(),
        "gene\tisoform_id\tS1\tS2\tS3\nnone\tlong\t1\t1\t0\nnone\tshort\t2\t0\t0\nnone\tzero\t0\t0\t0\n");
    assert_eq!(
        counts(&prefix.with_file_name("quant.v1.isoform_count.csv")),
        BTreeMap::from([
            ("short".into(), 2.0),
            ("long".into(), 2.0),
            ("zero".into(), 0.0)
        ])
    );
    let pairs = trackcluster_rs::count::read_read_to_isoform_tsv(
        prefix.with_file_name("quant.v1.read_to_isoform.tsv"),
    )
    .unwrap();
    assert_eq!(pairs.len(), 4);
    assert!(pairs.contains(&("S1::r_short".into(), "short".into())));
    assert!(pairs.contains(&("S2::r_short".into(), "long".into())));
    let group_path = prefix.with_file_name("quant.v1.isoform_usage.group.tsv");
    let mut group_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(group_path)
        .unwrap();
    let group_rows: Vec<_> = group_reader.records().map(Result::unwrap).collect();
    let long_case = group_rows
        .iter()
        .find(|row| &row[1] == "long" && &row[2] == "case")
        .unwrap();
    assert_eq!(&long_case[3], "1");
    assert_eq!(&long_case[4], "1");
    assert!(
        fs::read_to_string(prefix.with_file_name("quant.v1.unassigned_reads.tsv"))
            .unwrap()
            .contains("S1::r_unknown\tnone\tno_candidate_locus")
    );
    assert!(
        fs::read_to_string(prefix.with_file_name("quant.v1.unique_assignment.provenance.tsv"))
            .unwrap()
            .contains("assignment_source\tfixed_catalog\n")
    );
    assert!(!root.join("nested/quant.read_to_isoform.tsv").exists());
}

#[test]
fn catalog_count_preserves_zero_counts_when_all_reads_are_unassigned_or_empty() {
    let root = TestDir::new("catalog_count_empty");
    let (reads, catalog) = write_inputs(&root);
    let out = root.join("quant.csv");
    for contents in [UNASSIGNED_READ, ""] {
        fs::write(&reads, contents).unwrap();
        let output = count_command(&reads, &catalog, &out)
            .arg("--reference")
            .arg(&catalog)
            .output()
            .unwrap();
        assert_success(&output, "zero-count catalog");
        let records = counts(&out);
        assert_eq!(records.len(), 3);
        assert!(records.values().all(|&value| value == 0.0));
        assert_eq!(
            fs::read_to_string(out.with_extension("read_to_isoform.tsv")).unwrap(),
            ""
        );
    }
    assert_eq!(
        fs::read_to_string(out.with_extension("unassigned_reads.tsv")).unwrap(),
        "read_id\tgene\treason\n"
    );
}

#[test]
fn invalid_catalog_and_conflicting_read_alignments_fail_before_publishing() {
    let root = TestDir::new("catalog_count_invalid");
    let (reads, catalog) = write_inputs(&root);
    let out = root.join("quant.csv");
    fs::write(&out, "previous result\n").unwrap();
    for (read_input, catalog_input, error) in [
        (SHORT_READ.to_owned(), "".to_owned(), "catalog is empty"),
        (
            SHORT_READ.to_owned(),
            format!("{CATALOG}{CATALOG}"),
            "duplicate isoform id",
        ),
        (
            format!("{SHORT_READ}{}", LONG_READ.replace("r_long", "r_short")),
            CATALOG.to_owned(),
            "conflicting alignments",
        ),
    ] {
        fs::write(&reads, read_input).unwrap();
        fs::write(&catalog, catalog_input).unwrap();
        let output = count_command(&reads, &catalog, &out).output().unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains(error));
        assert_eq!(fs::read_to_string(&out).unwrap(), "previous result\n");
        assert!(!out.with_extension("assignment_stats.tsv").exists());
    }
}

#[test]
fn catalog_mode_rejects_fractional_counting_and_existing_mapping_inputs() {
    let root = TestDir::new("catalog_count_flags");
    let (reads, catalog) = write_inputs(&root);
    let out = root.join("quant.csv");
    let manifest = root.join("samples.tsv");
    fs::write(&manifest, "sample\treads\nS1\treads.bed\n").unwrap();
    for extra in [
        vec!["--assignment-mode", "fractional"],
        vec!["--read-to-isoform", "old.tsv"],
    ] {
        let single = count_command(&reads, &catalog, &out)
            .args(&extra)
            .output()
            .unwrap();
        let multi = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
            .args(["count-multi", "--assign-against-catalog", "--manifest"])
            .arg(&manifest)
            .arg("--isoform")
            .arg(&catalog)
            .arg("--out")
            .arg(root.join("multi"))
            .args(&extra)
            .output()
            .unwrap();
        assert!(!single.status.success());
        assert!(!multi.status.success());
    }
    let conflicting = count_command(&reads, &catalog, &out)
        .arg("--output-root")
        .arg(&*root)
        .output()
        .unwrap();
    assert!(!conflicting.status.success());
    assert!(!out.exists());
    assert!(!root.join("multi.isoform_count.csv").exists());
}

#[test]
fn catalog_reports_cannot_overwrite_inputs() {
    let root = TestDir::new("catalog_report_alias");
    let (_, catalog) = write_inputs(&root);
    let out = root.join("quant.csv");
    for extension in [
        "csv",
        "read_to_isoform.tsv",
        "unassigned_reads.tsv",
        "assignment_stats.tsv",
        "provenance.tsv",
    ] {
        let reads = out.with_extension(extension);
        fs::write(&reads, SHORT_READ).unwrap();
        let output = count_command(&reads, &catalog, &out).output().unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
        assert_eq!(fs::read_to_string(reads).unwrap(), SHORT_READ);
    }
    let manifest = root.join("samples.tsv");
    fs::write(&manifest, "sample\treads\nS1\tquant.read_to_isoform.tsv\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["count-multi", "--assign-against-catalog", "--manifest"])
        .arg(&manifest)
        .arg("--isoform")
        .arg(&catalog)
        .arg("--out")
        .arg(root.join("quant"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
    assert_eq!(
        fs::read_to_string(root.join("quant.read_to_isoform.tsv")).unwrap(),
        SHORT_READ
    );
}

#[test]
fn count_multi_rejects_empty_raw_read_ids_even_without_catalog_candidates() {
    let root = TestDir::new("catalog_empty_read_id");
    let (reads, catalog) = write_inputs(&root);
    fs::write(&reads, UNASSIGNED_READ.replace("r_unknown", "")).unwrap();
    let manifest = root.join("samples.tsv");
    fs::write(&manifest, "sample\treads\nS1\treads.bed\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["count-multi", "--assign-against-catalog", "--manifest"])
        .arg(manifest)
        .arg("--isoform")
        .arg(catalog)
        .arg("--out")
        .arg(root.join("quant"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("read id must not be empty"));
    assert!(!root.join("quant.isoform_count.csv").exists());
}
