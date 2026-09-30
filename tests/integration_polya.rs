mod common;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use common::{assert_success, TestDir};
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;
use sam::alignment::record::{data::field::Tag, Flags};
use sam::alignment::record_buf::{data::field::Value, Data, Sequence};

fn write_bam(path: &Path, reads: &[(&str, Option<Value>, Flags)]) {
    use sam::header::record::value::{
        map::{program::tag, Program},
        Map,
    };
    let program = Map::<Program>::builder()
        .insert(tag::NAME, "dorado")
        .insert(tag::VERSION, "0.9.1")
        .insert(
            tag::COMMAND_LINE,
            "dorado basecaller sup reads --estimate-poly-a",
        )
        .build()
        .unwrap();
    let header = sam::Header::builder()
        .add_program("dorado", program)
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for (name, value, flags) in reads {
        let data: Data = value
            .clone()
            .map(|value| (Tag::new(b'p', b't'), value))
            .into_iter()
            .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_name(*name)
            .set_flags(*flags)
            .set_sequence(Sequence::from(b"ACGTACGTAA".to_vec()))
            .set_data(data)
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn primary(name: &str, length: i32) -> (&str, Option<Value>, Flags) {
    (name, Some(Value::Int32(length)), Flags::UNMAPPED)
}

fn catalog(root: &Path) -> PathBuf {
    let path = root.join("isoforms.bed");
    let rows = ["long", "short", "zero"].map(|name| {
        format!("chr1\t100\t110\t{name}\t0\t+\t100\t110\t0\t1\t10,\t0,\tnone\tnone\tnone\t-1,\tisoform_anno\tGENE1\tnone\tnone\n")
    }).concat();
    fs::write(&path, rows).unwrap();
    path
}

fn tsv(path: &Path) -> Vec<HashMap<String, String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let header = reader.headers().unwrap().clone();
    reader
        .records()
        .map(|record| {
            header
                .iter()
                .zip(record.unwrap().iter())
                .map(|(k, v)| (k.to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

fn aggregate(root: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["polya-aggregate", "--bam"])
        .arg(root.join("reads.bam"))
        .args(["--sample", "S1", "--isoforms"])
        .arg(root.join("isoforms.bed"))
        .arg("--read-to-isoform")
        .arg(root.join("mapping.tsv"))
        .arg("--out")
        .arg(root.join("result.v1"))
        .output()
        .unwrap()
}

#[test]
fn dorado_tail_statistics_have_independent_expected_values_and_auditable_denominators() {
    let root = TestDir::new("polya-statistics");
    catalog(&root);
    let mut reads = vec![
        primary("a", 10),
        primary("b", 20),
        primary("c", 40),
        primary("d", 70),
        primary("anchor", -1),
        primary("failed", 0),
        ("no_tag", None, Flags::UNMAPPED),
        primary("short", 5),
        primary("a", 10),
        primary("unassigned", 999),
        (
            "a",
            Some(Value::Int32(999)),
            Flags::UNMAPPED | Flags::SECONDARY,
        ),
        (
            "a",
            Some(Value::Int32(999)),
            Flags::UNMAPPED | Flags::SUPPLEMENTARY,
        ),
    ];
    // Tail estimates come from signal; a valid estimate can exceed sequence length.
    reads[3].2 |= Flags::REVERSE_COMPLEMENTED;
    write_bam(&root.join("reads.bam"), &reads);
    fs::write(root.join("mapping.tsv"), "a\tlong\nb\tlong\nc\tlong\nd\tlong\nanchor\tlong\nfailed\tlong\nno_tag\tlong\nabsent\tlong\nshort\tshort\na\tlong\n").unwrap();
    assert_success(&aggregate(&root), "Dorado poly(A) aggregate");
    let rows = tsv(&root.join("result.v1.isoform_polya.tsv"));
    assert_eq!(rows.len(), 3);
    let long = rows.iter().find(|row| row["isoform_id"] == "long").unwrap();
    for (column, expected) in [
        ("gene", "GENE1"),
        ("assigned_reads", "8"),
        ("observed_reads", "7"),
        ("polya_reads", "4"),
        ("missing_bam_reads", "1"),
        ("missing_pt_reads", "1"),
        ("anchor_not_found_reads", "1"),
        ("estimation_failed_reads", "1"),
        ("polya_fraction", "0.500000"),
        ("polya_mean_nt", "35.000000"),
        ("polya_median_nt", "30.000000"),
        ("polya_q25_nt", "17.500000"),
        ("polya_q75_nt", "47.500000"),
        ("polya_min_nt", "10"),
        ("polya_max_nt", "70"),
        ("polya_stddev_nt", "26.457513"),
    ] {
        assert_eq!(long[column], expected, "{column}");
    }
    let short = rows
        .iter()
        .find(|row| row["isoform_id"] == "short")
        .unwrap();
    assert_eq!(short["polya_median_nt"], "5.000000");
    assert_eq!(short["polya_stddev_nt"], "NA");
    let zero = rows.iter().find(|row| row["isoform_id"] == "zero").unwrap();
    assert_eq!(zero["assigned_reads"], "0");
    assert_eq!(zero["polya_mean_nt"], "NA");
    assert_eq!(zero["polya_fraction"], "NA");
    let calls = tsv(&root.join("result.v1.read_polya.tsv"));
    assert_eq!(calls.len(), 9);
    assert_eq!(
        calls.iter().find(|row| row["read_id"] == "anchor").unwrap()["dorado_pt"],
        "-1"
    );
    assert_eq!(
        calls.iter().find(|row| row["read_id"] == "failed").unwrap()["polya_length_nt"],
        "NA"
    );
    assert_eq!(
        calls.iter().find(|row| row["read_id"] == "absent").unwrap()["status"],
        "missing_bam_read"
    );
    assert_eq!(
        calls.iter().find(|row| row["read_id"] == "no_tag").unwrap()["status"],
        "missing_pt_tag"
    );
    let qc = &tsv(&root.join("result.v1.polya_qc.tsv"))[0];
    assert_eq!(qc["dorado_versions"], "0.9.1");
    assert_eq!(qc["bam_records"], "12");
    assert_eq!(qc["skipped_secondary"], "1");
    assert_eq!(qc["skipped_supplementary"], "1");
    assert_eq!(qc["unassigned_primary_records"], "1");
    assert_eq!(qc["duplicate_assigned_primary_records"], "1");
    assert_eq!(qc["assigned_reads"], "9");
    assert_eq!(qc["read_join_rate"], "0.888889");
}

#[test]
fn manifest_keeps_identical_read_names_in_different_samples_separate() {
    let root = TestDir::new("polya-multi");
    let isoforms = catalog(&root);
    write_bam(&root.join("s1.bam"), &[primary("shared", 15)]);
    write_bam(&root.join("s2.bam"), &[primary("shared", 150)]);
    fs::write(
        root.join("polya.tsv"),
        "sample\tgroup\tbam\nS1\tcontrol\ts1.bam\nS2\ttreated\ts2.bam\n",
    )
    .unwrap();
    fs::write(
        root.join("mapping.tsv"),
        "S1::shared\tlong\nS2::shared\tlong\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["polya-aggregate", "--polya-manifest"])
        .arg(root.join("polya.tsv"))
        .arg("--isoforms")
        .arg(isoforms)
        .arg("--read-to-isoform")
        .arg(root.join("mapping.tsv"))
        .arg("--out")
        .arg(root.join("result"))
        .current_dir(root.path().parent().unwrap())
        .output()
        .unwrap();
    assert_success(&output, "poly(A) manifest join");
    let rows = tsv(&root.join("result.isoform_polya.tsv"));
    assert_eq!(rows.len(), 6);
    for (sample, group, length) in [
        ("S1", "control", "15.000000"),
        ("S2", "treated", "150.000000"),
    ] {
        let row = rows
            .iter()
            .find(|row| row["sample"] == sample && row["isoform_id"] == "long")
            .unwrap();
        assert_eq!(row["group"], group);
        assert_eq!(row["polya_median_nt"], length);
        assert_eq!(row["assigned_reads"], "1");
    }
}

#[test]
fn invalid_tags_and_conflicting_primary_estimates_fail_before_publication() {
    for (reads, expected_error) in [
        (vec![primary("a", -2)], "invalid Dorado pt value"),
        (
            vec![("a", Some(Value::Float(10.0)), Flags::UNMAPPED)],
            "pt tag must have integer type",
        ),
        (vec![("a", None, Flags::UNMAPPED)], "no Dorado pt:i tags"),
        (
            vec![primary("a", 10), primary("a", 20)],
            "conflicting Dorado pt values",
        ),
        (
            vec![primary("a", 10), ("a", None, Flags::UNMAPPED)],
            "conflicting Dorado pt values",
        ),
    ] {
        let root = TestDir::new("polya-invalid-tag");
        catalog(&root);
        fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
        write_bam(&root.join("reads.bam"), &reads);
        let output = aggregate(&root);
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(expected_error),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        for suffix in ["isoform_polya.tsv", "read_polya.tsv", "polya_qc.tsv"] {
            assert!(!root.join(format!("result.v1.{suffix}")).exists());
        }
    }
}

#[test]
fn ambiguous_or_unknown_isoform_assignments_are_rejected() {
    for (mapping, expected_error) in [
        ("a\tlong\na\tshort\n", "unique read-to-isoform mapping"),
        ("a\tunknown\n", "unknown isoform"),
    ] {
        let root = TestDir::new("polya-invalid-assignment");
        catalog(&root);
        write_bam(&root.join("reads.bam"), &[primary("a", 10)]);
        fs::write(root.join("mapping.tsv"), mapping).unwrap();
        let output = aggregate(&root);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains(expected_error));
        assert!(!root.join("result.v1.isoform_polya.tsv").exists());
    }
}

#[test]
fn failed_estimates_and_zero_read_id_join_are_explicitly_missing_lengths() {
    let root = TestDir::new("polya-no-valid-tails");
    catalog(&root);
    write_bam(
        &root.join("reads.bam"),
        &[primary("a", -1), primary("b", 0)],
    );
    fs::write(root.join("mapping.tsv"), "a\tlong\nb\tlong\n").unwrap();
    assert_success(&aggregate(&root), "failed-estimate summary");
    let long = tsv(&root.join("result.v1.isoform_polya.tsv"))
        .into_iter()
        .find(|row| row["isoform_id"] == "long")
        .unwrap();
    assert_eq!(long["polya_reads"], "0");
    assert_eq!(long["polya_median_nt"], "NA");
    assert_eq!(long["polya_fraction"], "0.000000");
    fs::write(root.join("mapping.tsv"), "renamed_read\tlong\n").unwrap();
    assert_success(&aggregate(&root), "zero-join audit");
    let qc = &tsv(&root.join("result.v1.polya_qc.tsv"))[0];
    assert_eq!(qc["missing_bam_reads"], "1");
    assert_eq!(qc["read_join_rate"], "0.000000");
}

fn fixture(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(path)
}

#[test]
fn flow_adds_polya_after_unique_assignment_and_removes_stale_summaries() {
    let root = TestDir::new("flow-polya-single");
    write_bam(&root.join("reads.bam"), &[primary("read_trunc", 80)]);
    let output_root = root.join("flow");
    let run_flow = |polya: bool| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
        command
            .args(["flow", "--reads"])
            .arg(fixture("reads.bed"))
            .arg("--reference")
            .arg(fixture("ref.bed"))
            .arg("--output-root")
            .arg(&output_root)
            .args(["--prefix", "sample", "--heartbeat-seconds", "0"]);
        if polya {
            command
                .arg("--polya-bam")
                .arg(root.join("reads.bam"))
                .args(["--polya-sample", "S1"]);
        }
        command.output().unwrap()
    };
    assert_success(&run_flow(true), "single-sample flow poly(A)");
    let calls = tsv(&output_root.join("sample.read_polya.tsv"));
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0]["read_id"], "read_trunc");
    assert_eq!(calls[0]["polya_length_nt"], "80");
    let mapping =
        fs::read_to_string(output_root.join("sample_read_to_isoform.unique.tsv")).unwrap();
    assert!(mapping.contains(&format!("read_trunc\t{}\n", calls[0]["isoform_id"])));
    let original = fs::read(output_root.join("sample.isoform_polya.tsv")).unwrap();
    let recount = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["flow", "--count-only", "--reference"])
        .arg(fixture("ref.bed"))
        .arg("--output-root")
        .arg(&output_root)
        .args(["--prefix", "sample", "--heartbeat-seconds", "0"])
        .arg("--polya-bam")
        .arg(root.join("reads.bam"))
        .args(["--polya-sample", "S1"])
        .output()
        .unwrap();
    assert_success(&recount, "poly(A) count-only flow");
    assert_eq!(
        fs::read(output_root.join("sample.isoform_polya.tsv")).unwrap(),
        original
    );
    assert_success(&run_flow(false), "flow without poly(A)");
    for suffix in ["isoform_polya.tsv", "read_polya.tsv", "polya_qc.tsv"] {
        assert!(!output_root.join(format!("sample.{suffix}")).exists());
    }
}

#[test]
fn multisample_flow_inherits_groups_from_reads_manifest() {
    let root = TestDir::new("flow-polya-multi");
    write_bam(&root.join("s1.bam"), &[primary("read_s1", 25)]);
    write_bam(&root.join("s2.bam"), &[primary("read_s2", 125)]);
    fs::write(
        root.join("polya.tsv"),
        "sample\tbam\nS1\ts1.bam\nS2\ts2.bam\n",
    )
    .unwrap();
    let output_root = root.join("flow");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["flow", "--manifest"])
        .arg(fixture("samples.tsv"))
        .arg("--reference")
        .arg(fixture("ref.bed"))
        .arg("--output-root")
        .arg(&output_root)
        .args(["--prefix", "pooled", "--heartbeat-seconds", "0"])
        .arg("--polya-manifest")
        .arg(root.join("polya.tsv"))
        .output()
        .unwrap();
    assert_success(&output, "multi-sample flow poly(A)");
    let calls = tsv(&output_root.join("pooled.read_polya.tsv"));
    assert_eq!(calls.len(), 2);
    assert_eq!(calls[0]["read_id"], "S1::read_s1");
    assert_eq!(calls[0]["group"], "control");
    assert_eq!(calls[1]["read_id"], "S2::read_s2");
    assert_eq!(calls[1]["group"], "treated");
    assert_eq!(calls[1]["polya_length_nt"], "125");
}

#[test]
fn polya_outputs_cannot_overwrite_input_bam() {
    let root = TestDir::new("polya-output-alias");
    catalog(&root);
    fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
    let bam = root.join("result.v1.read_polya.tsv");
    write_bam(&bam, &[primary("a", 10)]);
    let original = fs::read(&bam).unwrap();
    fs::hard_link(&bam, root.join("reads.bam")).unwrap();
    let output = aggregate(&root);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
    assert_eq!(fs::read(&bam).unwrap(), original);
}

#[test]
fn flow_rejects_fractional_polya_before_creating_output() {
    let root = TestDir::new("flow-polya-fractional");
    write_bam(&root.join("reads.bam"), &[primary("read_trunc", 80)]);
    let output_root = root.join("flow");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["flow", "--reads"])
        .arg(fixture("reads.bed"))
        .arg("--reference")
        .arg(fixture("ref.bed"))
        .arg("--output-root")
        .arg(&output_root)
        .args(["--prefix", "sample", "--assignment-mode", "fractional"])
        .arg("--polya-bam")
        .arg(root.join("reads.bam"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("poly(A) aggregation requires --assignment-mode unique"));
    assert!(!output_root.exists());
}

#[test]
fn flow_polya_inputs_follow_the_existing_output_root_ownership_contract() {
    for hard_link in [false, true] {
        let root = TestDir::new("flow-polya-input-ownership");
        let output_root = root.join("flow");
        fs::create_dir(&output_root).unwrap();
        let bam = if hard_link {
            root.join("reads.bam")
        } else {
            output_root.join("reads.bam")
        };
        write_bam(&bam, &[primary("read_trunc", 80)]);
        let original = fs::read(&bam).unwrap();
        if hard_link {
            fs::hard_link(&bam, output_root.join("existing_artifact")).unwrap();
        }
        let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
            .args(["flow", "--reads"])
            .arg(fixture("reads.bed"))
            .arg("--reference")
            .arg(fixture("ref.bed"))
            .arg("--output-root")
            .arg(&output_root)
            .args(["--prefix", "sample"])
            .arg("--polya-bam")
            .arg(&bam)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("pipeline-owned output_root"));
        assert_eq!(fs::read(&bam).unwrap(), original);
        assert_eq!(fs::read_dir(&output_root).unwrap().count(), 1);
    }
}
