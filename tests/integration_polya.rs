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

const NANOPOLISH_HEADER: &str = "readname\tcontig\tposition\tleader_start\tadapter_start\tpolya_start\ttranscript_start\tread_rate\tpolya_length\tqc_tag\n";

fn write_nanopolish(path: &Path, rows: &[(&str, &str, &str)], header: bool) {
    let mut text = if header {
        NANOPOLISH_HEADER.to_owned()
    } else {
        String::new()
    };
    for (read, length, qc) in rows {
        // Source alignment coordinates do not affect the molecule-level join.
        text.push_str(&format!(
            "{read}\tother_contig\t999\t0\t1\t2\t3\t100.0\t{length}\t{qc}\n"
        ));
    }
    fs::write(path, text).unwrap();
}

fn aggregate_nanopolish(root: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["polya-aggregate", "--nanopolish"])
        .arg(root.join("tails.tsv"))
        .args(["--sample", "S1", "--isoforms"])
        .arg(root.join("isoforms.bed"))
        .arg("--read-to-isoform")
        .arg(root.join("mapping.tsv"))
        .arg("--out")
        .arg(root.join("result.v1"))
        .output()
        .unwrap()
}

fn aggregate_manifest(root: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["polya-aggregate", "--polya-manifest"])
        .arg(root.join("polya.tsv"))
        .arg("--isoforms")
        .arg(root.join("isoforms.bed"))
        .arg("--read-to-isoform")
        .arg(root.join("mapping.tsv"))
        .arg("--out")
        .arg(root.join("result.v1"))
        .current_dir(root.parent().unwrap())
        .output()
        .unwrap()
}

#[test]
fn nanopolish_float_statistics_use_pass_molecules_and_keep_raw_qc_audits() {
    let root = TestDir::new("polya-nanopolish-statistics");
    catalog(&root);
    // Five PASS lengths from the Nanopolish documentation; expectations below
    // were computed independently, including the sample standard deviation.
    let mut rows = vec![
        ("a", "38.22", "PASS"),
        ("b", "27.48", "PASS"),
        ("c", "29.16", "PASS"),
        ("d", "32.43", "PASS"),
        ("e", "39.65", "PASS"),
        ("clip", "23.00", "SUFFCLIP"),
        ("load", "-1.00", "READ_FAILED_LOAD"),
        ("unassigned", "999.00", "PASS"),
        ("a", "38.220", "PASS"),
        ("a", "38.22", "SUFFCLIP"),
        ("clip", "-1.00", "ADAPTER"),
    ];
    write_nanopolish(&root.join("tails.tsv"), &rows, true);
    fs::write(root.join("mapping.tsv"), "a\tlong\nb\tlong\nc\tlong\nd\tlong\ne\tlong\nclip\tlong\nload\tlong\nabsent\tlong\na\tlong\n").unwrap();
    assert_success(&aggregate_nanopolish(&root), "Nanopolish statistics");
    let long = tsv(&root.join("result.v1.isoform_polya.tsv"))
        .into_iter()
        .find(|row| row["isoform_id"] == "long")
        .unwrap();
    for (column, expected) in [
        ("caller", "nanopolish"),
        ("assigned_reads", "8"),
        ("observed_reads", "7"),
        ("polya_reads", "5"),
        ("missing_bam_reads", "0"),
        ("missing_pt_reads", "0"),
        ("anchor_not_found_reads", "0"),
        ("estimation_failed_reads", "0"),
        ("missing_nanopolish_reads", "1"),
        ("qc_failed_reads", "2"),
        ("polya_fraction", "0.625000"),
        ("polya_mean_nt", "33.388000"),
        ("polya_median_nt", "32.430000"),
        ("polya_q25_nt", "29.160000"),
        ("polya_q75_nt", "38.220000"),
        ("polya_min_nt", "27.48"),
        ("polya_max_nt", "39.65"),
        ("polya_stddev_nt", "5.391175"),
    ] {
        assert_eq!(long[column], expected, "{column}");
    }
    let calls = tsv(&root.join("result.v1.read_polya.tsv"));
    assert_eq!(calls.len(), 8);
    let a = calls.iter().find(|row| row["read_id"] == "a").unwrap();
    assert_eq!(a["dorado_pt"], "NA");
    assert_eq!(a["polya_length_nt"], "38.22");
    assert_eq!(a["nanopolish_polya_length_nt"], "38.22");
    assert_eq!(a["nanopolish_qc_tag"], "PASS");
    let failed = calls.iter().find(|row| row["read_id"] == "load").unwrap();
    assert_eq!(failed["status"], "qc_failed");
    assert_eq!(failed["polya_length_nt"], "NA");
    assert_eq!(failed["nanopolish_polya_length_nt"], "-1.00");
    assert_eq!(failed["nanopolish_qc_tag"], "READ_FAILED_LOAD");
    let clip = calls.iter().find(|row| row["read_id"] == "clip").unwrap();
    assert_eq!(clip["nanopolish_qc_tag"], "ADAPTER");
    let absent = calls.iter().find(|row| row["read_id"] == "absent").unwrap();
    assert_eq!(absent["status"], "missing_nanopolish_read");
    assert_eq!(absent["nanopolish_qc_tag"], "NA");
    let qc = &tsv(&root.join("result.v1.polya_qc.tsv"))[0];
    for (column, expected) in [
        ("caller", "nanopolish"),
        ("bam", "NA"),
        ("bam_records", "NA"),
        ("nanopolish_rows", "11"),
        ("nanopolish_pass_rows", "7"),
        ("nanopolish_qc_failed_rows", "4"),
        ("unassigned_nanopolish_rows", "1"),
        ("duplicate_assigned_nanopolish_rows", "3"),
        ("assigned_reads", "8"),
        ("read_join_rate", "0.875000"),
    ] {
        assert_eq!(qc[column], expected, "{column}");
    }
    let tag_counts: serde_json::Value =
        serde_json::from_str(&qc["nanopolish_qc_tag_counts"]).unwrap();
    assert_eq!(
        tag_counts,
        serde_json::json!({"ADAPTER": 1, "PASS": 7, "READ_FAILED_LOAD": 1, "SUFFCLIP": 2})
    );
    let snapshots = ["isoform_polya.tsv", "read_polya.tsv", "polya_qc.tsv"].map(|suffix| {
        (
            suffix,
            fs::read(root.join(format!("result.v1.{suffix}"))).unwrap(),
        )
    });
    rows.reverse();
    write_nanopolish(&root.join("tails.tsv"), &rows, true);
    assert_success(&aggregate_nanopolish(&root), "reordered Nanopolish rows");
    for (suffix, original) in snapshots {
        assert_eq!(
            fs::read(root.join(format!("result.v1.{suffix}"))).unwrap(),
            original
        );
    }
}

#[test]
fn nanopolish_pass_zero_and_fractional_lengths_are_valid_with_reordered_columns() {
    let root = TestDir::new("polya-nanopolish-zero");
    catalog(&root);
    fs::write(root.join("mapping.tsv"), "zero\tlong\nfraction\tlong\n").unwrap();
    fs::write(
        root.join("tails.tsv"),
        "qc_tag\tpolya_length\treadname\nPASS\t-0.00\tzero\nPASS\t12.5\tfraction\n",
    )
    .unwrap();
    assert_success(&aggregate_nanopolish(&root), "Nanopolish PASS zero");
    let long = tsv(&root.join("result.v1.isoform_polya.tsv"))
        .into_iter()
        .find(|row| row["isoform_id"] == "long")
        .unwrap();
    assert_eq!(long["polya_reads"], "2");
    assert_eq!(long["polya_mean_nt"], "6.250000");
    assert_eq!(long["polya_min_nt"], "0");
    assert_eq!(long["polya_max_nt"], "12.5");
    assert_eq!(long["estimation_failed_reads"], "0");
    let zero = tsv(&root.join("result.v1.read_polya.tsv"))
        .into_iter()
        .find(|row| row["read_id"] == "zero")
        .unwrap();
    assert_eq!(zero["status"], "estimated");
    assert_eq!(zero["polya_length_nt"], "0");
    assert_eq!(zero["nanopolish_polya_length_nt"], "-0.00");
}

#[test]
fn nanopolish_accepts_headerless_pass_exports_and_empty_results() {
    let root = TestDir::new("polya-nanopolish-headerless");
    catalog(&root);
    fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
    write_nanopolish(&root.join("tails.tsv"), &[("a", "10.25", "PASS")], false);
    assert_success(
        &aggregate_nanopolish(&root),
        "headerless Nanopolish PASS export",
    );
    assert_eq!(
        tsv(&root.join("result.v1.read_polya.tsv"))[0]["polya_length_nt"],
        "10.25"
    );
    for contents in ["", NANOPOLISH_HEADER] {
        fs::write(root.join("tails.tsv"), contents).unwrap();
        assert_success(&aggregate_nanopolish(&root), "empty Nanopolish result");
        let qc = &tsv(&root.join("result.v1.polya_qc.tsv"))[0];
        assert_eq!(qc["nanopolish_rows"], "0");
        assert_eq!(qc["missing_nanopolish_reads"], "1");
        assert_eq!(qc["read_join_rate"], "0.000000");
    }
}

#[test]
fn invalid_nanopolish_pass_lengths_and_conflicting_estimates_fail_before_publication() {
    for length in ["NaN", "inf", "-1", "broken"] {
        let root = TestDir::new("polya-nanopolish-invalid-length");
        catalog(&root);
        fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
        write_nanopolish(&root.join("tails.tsv"), &[("a", length, "PASS")], true);
        let output = aggregate_nanopolish(&root);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("Nanopolish PASS polya_length"));
        for suffix in ["isoform_polya.tsv", "read_polya.tsv", "polya_qc.tsv"] {
            assert!(!root.join(format!("result.v1.{suffix}")).exists());
        }
    }
    let root = TestDir::new("polya-nanopolish-conflict");
    catalog(&root);
    fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
    write_nanopolish(
        &root.join("tails.tsv"),
        &[("a", "10.25", "PASS"), ("a", "20.5", "PASS")],
        true,
    );
    let output = aggregate_nanopolish(&root);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("conflicting Nanopolish PASS"));
    assert!(!root.join("result.v1.isoform_polya.tsv").exists());
    fs::write(root.join("tails.tsv"), "readname\tpolya_length\na\t10.0\n").unwrap();
    let output = aggregate_nanopolish(&root);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("header must include"));
    // Failed-QC rows can contain the caller's missing/invalid numeric sentinels.
    write_nanopolish(
        &root.join("tails.tsv"),
        &[("a", "nan", "NOREGION"), ("a", "NA", "READ_FAILED_LOAD")],
        true,
    );
    assert_success(
        &aggregate_nanopolish(&root),
        "failed-QC Nanopolish sentinels",
    );
    assert_eq!(
        tsv(&root.join("result.v1.polya_qc.tsv"))[0]["qc_failed_reads"],
        "1"
    );
}

#[test]
fn mixed_caller_manifest_keeps_samples_separate_and_rejects_multiple_sources() {
    let root = TestDir::new("polya-mixed-callers");
    catalog(&root);
    write_bam(&root.join("reads.bam"), &[primary("shared", 10)]);
    write_nanopolish(&root.join("tails.tsv"), &[("shared", "10.5", "PASS")], true);
    fs::write(
        root.join("mapping.tsv"),
        "S1::shared\tlong\nS2::shared\tlong\n",
    )
    .unwrap();
    fs::write(
        root.join("polya.tsv"),
        "sample\tbam\tnanopolish\tgroup\nS1\treads.bam\tNA\tcontrol\nS2\t\ttails.tsv\ttreated\n",
    )
    .unwrap();
    assert_success(&aggregate_manifest(&root), "mixed-caller manifest");
    let rows = tsv(&root.join("result.v1.isoform_polya.tsv"));
    for (sample, caller, length) in [
        ("S1", "dorado", "10.000000"),
        ("S2", "nanopolish", "10.500000"),
    ] {
        let row = rows
            .iter()
            .find(|row| row["sample"] == sample && row["isoform_id"] == "long")
            .unwrap();
        assert_eq!(row["caller"], caller);
        assert_eq!(row["polya_mean_nt"], length);
        assert_eq!(row["assigned_reads"], "1");
    }
    let valid_summary = fs::read(root.join("result.v1.isoform_polya.tsv")).unwrap();
    for contents in [
        "sample\tbam\tnanopolish\nS1\treads.bam\ttails.tsv\n",
        "sample\tbam\tnanopolish\nS1\tNA\t\n",
        "sample\tnanopolish\nS1\ttails.tsv\nS1\ttails.tsv\n",
    ] {
        fs::write(root.join("polya.tsv"), contents).unwrap();
        let output = aggregate_manifest(&root);
        assert!(!output.status.success());
        assert_eq!(
            fs::read(root.join("result.v1.isoform_polya.tsv")).unwrap(),
            valid_summary
        );
    }
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
fn failed_flow_reruns_cannot_keep_polya_tables_for_replaced_assignments() {
    for with_polya in [true, false] {
        let root = TestDir::new("flow-polya-failed-rerun");
        let reads = root.join("reads.bed");
        fs::copy(fixture("reads.bed"), &reads).unwrap();
        write_nanopolish(
            &root.join("tails.tsv"),
            &[("read_trunc", "80.25", "PASS")],
            true,
        );
        let output_root = root.join("flow");
        let run_flow = |polya: bool| {
            let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
            command
                .args(["flow", "--reads"])
                .arg(&reads)
                .arg("--reference")
                .arg(fixture("ref.bed"))
                .arg("--output-root")
                .arg(&output_root)
                .args([
                    "--prefix",
                    "sample",
                    "--threads",
                    "1",
                    "--heartbeat-seconds",
                    "0",
                ]);
            if polya {
                command
                    .arg("--polya-nanopolish")
                    .arg(root.join("tails.tsv"))
                    .args(["--polya-sample", "S1"]);
            }
            command.output().unwrap()
        };
        assert_success(&run_flow(true), "initial flow with poly(A)");
        assert!(output_root.join("sample.isoform_polya.tsv").is_file());

        // Force a late publication failure after the selected mapping changes.
        let description = output_root.join("sample_desc.txt");
        fs::remove_file(&description).unwrap();
        fs::create_dir(&description).unwrap();
        let new_reads = fs::read_to_string(&reads)
            .unwrap()
            .replace("read_trunc", "read_changed");
        fs::write(&reads, new_reads).unwrap();
        let output = run_flow(with_polya);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("publish derived output"));
        let mapping =
            fs::read_to_string(output_root.join("sample_read_to_isoform.unique.tsv")).unwrap();
        assert!(mapping.starts_with("read_changed\t"));
        for suffix in ["isoform_polya.tsv", "read_polya.tsv", "polya_qc.tsv"] {
            assert!(
                !output_root.join(format!("sample.{suffix}")).exists(),
                "stale {suffix} survived failed rerun (poly(A) requested: {with_polya})"
            );
        }
    }
}

#[test]
fn nanopolish_flow_joins_final_assignments_and_regenerates_in_count_only_mode() {
    let root = TestDir::new("flow-polya-nanopolish");
    write_nanopolish(
        &root.join("tails.tsv"),
        &[("read_trunc", "80.25", "PASS")],
        true,
    );
    let output_root = root.join("flow");
    let run_flow = |count_only: bool| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
        command.arg("flow");
        if count_only {
            command.arg("--count-only");
        } else {
            command.arg("--reads").arg(fixture("reads.bed"));
        }
        command
            .arg("--reference")
            .arg(fixture("ref.bed"))
            .arg("--output-root")
            .arg(&output_root)
            .args(["--prefix", "sample", "--heartbeat-seconds", "0"])
            .arg("--polya-nanopolish")
            .arg(root.join("tails.tsv"))
            .args(["--polya-sample", "S1"])
            .output()
            .unwrap()
    };
    assert_success(&run_flow(false), "single-sample Nanopolish flow");
    let calls = tsv(&output_root.join("sample.read_polya.tsv"));
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0]["sample"], "S1");
    assert_eq!(calls[0]["caller"], "nanopolish");
    assert_eq!(calls[0]["polya_length_nt"], "80.25");
    let mapping =
        fs::read_to_string(output_root.join("sample_read_to_isoform.unique.tsv")).unwrap();
    assert!(mapping.contains(&format!("read_trunc\t{}\n", calls[0]["isoform_id"])));
    let original = fs::read(output_root.join("sample.isoform_polya.tsv")).unwrap();
    assert_success(&run_flow(true), "count-only Nanopolish flow");
    assert_eq!(
        fs::read(output_root.join("sample.isoform_polya.tsv")).unwrap(),
        original
    );
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
fn mixed_caller_flow_inherits_groups_and_rejects_conflicting_manifest_groups() {
    let root = TestDir::new("flow-polya-mixed-callers");
    write_bam(&root.join("s1.bam"), &[primary("read_s1", 25)]);
    write_nanopolish(&root.join("s2.tsv"), &[("read_s2", "125.5", "PASS")], true);
    fs::write(
        root.join("polya.tsv"),
        "sample\tbam\tnanopolish\nS1\ts1.bam\t\nS2\t\ts2.tsv\n",
    )
    .unwrap();
    let output_root = root.join("flow");
    let run_flow = || {
        Command::new(env!("CARGO_BIN_EXE_trackcluster"))
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
            .unwrap()
    };
    assert_success(&run_flow(), "mixed-caller flow");
    let calls = tsv(&output_root.join("pooled.read_polya.tsv"));
    assert_eq!(calls.len(), 2);
    assert_eq!(calls[0]["read_id"], "S1::read_s1");
    assert_eq!(calls[0]["caller"], "dorado");
    assert_eq!(calls[0]["group"], "control");
    assert_eq!(calls[1]["read_id"], "S2::read_s2");
    assert_eq!(calls[1]["caller"], "nanopolish");
    assert_eq!(calls[1]["group"], "treated");
    assert_eq!(calls[1]["polya_length_nt"], "125.5");
    let original = fs::read(output_root.join("pooled.isoform_polya.tsv")).unwrap();
    fs::write(
        root.join("polya.tsv"),
        "sample\tbam\tnanopolish\tgroup\nS1\ts1.bam\t\tcontrol\nS2\t\ts2.tsv\tcontrol\n",
    )
    .unwrap();
    let output = run_flow();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("conflicts with --manifest"));
    assert_eq!(
        fs::read(output_root.join("pooled.isoform_polya.tsv")).unwrap(),
        original
    );
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
fn nanopolish_output_alias_and_multiple_caller_flags_fail_without_overwriting_inputs() {
    let root = TestDir::new("polya-nanopolish-output-alias");
    catalog(&root);
    fs::write(root.join("mapping.tsv"), "a\tlong\n").unwrap();
    let input = root.join("result.v1.read_polya.tsv");
    write_nanopolish(&input, &[("a", "10.5", "PASS")], true);
    fs::hard_link(&input, root.join("tails.tsv")).unwrap();
    let original = fs::read(&input).unwrap();
    let output = aggregate_nanopolish(&root);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
    assert_eq!(fs::read(&input).unwrap(), original);
    for args in [
        vec![
            "polya-aggregate",
            "--bam",
            "reads.bam",
            "--nanopolish",
            "tails.tsv",
            "--sample",
            "S1",
            "--isoforms",
            "isoforms.bed",
            "--read-to-isoform",
            "mapping.tsv",
            "--out",
            "rejected",
        ],
        vec![
            "flow",
            "--reads",
            "reads.bed",
            "--reference",
            "ref.bed",
            "--output-root",
            "rejected",
            "--polya-bam",
            "reads.bam",
            "--polya-nanopolish",
            "tails.tsv",
        ],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
            .args(args)
            .current_dir(&root)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(String::from_utf8_lossy(&output.stderr).contains("cannot be used with"));
        assert!(!root.join("rejected").exists());
    }
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
