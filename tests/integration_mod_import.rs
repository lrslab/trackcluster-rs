mod common;

use std::fs;
use std::num::NonZero;
use std::process::Command;

use common::TestDir;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;
use sam::alignment::record::cigar::{op::Kind, Op};

#[test]
fn m6anet_import_cli_projects_and_writes_auditable_outputs() {
    let root = TestDir::new("m6anet-import-e2e");
    let gtf = root.join("reference.gtf");
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\ttranscript\t101\t220\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1.1\";\n",
            "chr1\ttest\texon\t101\t120\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1.1\";\n",
            "chr1\ttest\texon\t201\t220\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1.1\";\n",
        ),
    )
    .unwrap();
    let indiv = root.join("data.indiv_proba.csv");
    fs::write(
        &indiv,
        concat!(
            "transcript_id,transcript_position,read_index,probability_modified\n",
            "tx1.1,5,7,0.9\n",
            "tx1.1,25,8,0.2\n",
        ),
    )
    .unwrap();
    let read_map = root.join("read_map.tsv");
    fs::write(&read_map, "read_index\tread_id\n7\tread-a\n8\tread-b\n").unwrap();
    let data_info = root.join("data.info");
    fs::write(
        &data_info,
        concat!(
            "transcript_id,transcript_position,start,end,n_reads\n",
            "tx1.1,5,0,1,1\n",
            "tx1.1,25,1,2,1\n",
        ),
    )
    .unwrap();
    let site_proba = root.join("data.site_proba.csv");
    fs::write(
        &site_proba,
        concat!(
            "transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio\n",
            "tx1.1,5,1,0.95,GGACT,1\n",
            "tx1.1,25,1,0.2,GGACT,1\n",
        ),
    )
    .unwrap();

    let prefix = root.join("sample.mod");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "mod-import-m6anet",
            "--sample",
            "S1",
            "--assay-id",
            "m6anet_hct116_rna002",
            "--model-id",
            "HCT116_RNA002",
            "--input-format",
            "gtf",
            "--min-reads",
            "1",
            "--indiv",
        ])
        .arg(&indiv)
        .arg("--data-info")
        .arg(&data_info)
        .arg("--site-proba")
        .arg(&site_proba)
        .arg("--read-map")
        .arg(&read_map)
        .arg("--reference")
        .arg(&gtf)
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "mod-import-m6anet: sample=S1 assay=m6anet_hct116_rna002 input_rows=2 observations=2 sites=2\n"
    );

    let observations = fs::read_to_string(root.join("sample.mod.observations.tsv")).unwrap();
    assert!(observations.contains(
        "m6anet_hct116_rna002\tS1\tS1::read-a\tchr1\t105\t+\tA+a\t0.9\texplicit_probability\tNA\ttx1.1\t5\n"
    ));
    assert!(observations.contains(
        "m6anet_hct116_rna002\tS1\tS1::read-b\tchr1\t205\t+\tA+a\t0.2\texplicit_probability\tNA\ttx1.1\t25\n"
    ));

    let assay = fs::read_to_string(root.join("sample.mod.assay.json")).unwrap();
    assert!(assay.contains("\"candidate_observations_complete\": true"));
    assert!(assay.contains("\"chemistry\": \"RNA002\""));
    let qc = fs::read_to_string(root.join("sample.mod.import_qc.tsv")).unwrap();
    assert!(qc.contains("data_info_retained_sites\t2\n"));
    assert!(qc.contains("read_map_entries_used\t2\n"));
    assert!(qc.contains("site_probability_sites\t2\n"));
    assert!(qc.contains("site_probability_sites_at_or_above_threshold\t1\n"));
    assert!(qc.contains("read_probability_threshold\t0.033379376\n"));
}

#[test]
fn dorado_import_cli_emits_explicit_and_implicit_candidates() {
    use sam::alignment::record_buf::{
        data::field::{value::Array as BufArray, Value as BufValue},
        Cigar, Data, Sequence,
    };
    use sam::header::record::value::{map::ReferenceSequence, Map};

    let root = TestDir::new("dorado-import-e2e");
    let bam_path = root.join("calls.bam");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
    let data: Data = [
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATIONS,
            BufValue::from("A+a?,0;"),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_PROBABILITIES,
            BufValue::Array(BufArray::UInt8(vec![255])),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            BufValue::UInt32(4),
        ),
    ]
    .into_iter()
    .collect();
    let record = sam::alignment::RecordBuf::builder()
        .set_name("read1")
        .set_flags(sam::alignment::record::Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start("101".parse().unwrap())
        .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(b"AAAA".to_vec()))
        .set_data(data)
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();

    let prefix = root.join("sample.mod");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "mod-import-dorado",
            "--sample",
            "S1",
            "--assay-id",
            "dorado_rna004_m6a",
            "--mod-code",
            "A+a",
            "--model-id",
            "rna004-test-m6a",
            "--source-emission-threshold",
            "0.05",
            "--question-mark-policy",
            "below-emission-threshold",
            "--bam",
        ])
        .arg(&bam_path)
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "mod-import-dorado: sample=S1 assay=dorado_rna004_m6a records=1 retained=1 observations=4 complete=true\n"
    );

    let observations = fs::read_to_string(root.join("sample.mod.observations.tsv")).unwrap();
    assert_eq!(observations.lines().count(), 5);
    assert_eq!(observations.matches("\texplicit_probability\t").count(), 1);
    assert_eq!(
        observations
            .matches("\timplicit_below_emission_threshold\t")
            .count(),
        3
    );
    let metadata = fs::read_to_string(root.join("sample.mod.assay.json")).unwrap();
    assert!(metadata.contains("\"candidate_observations_complete\": true"));
    assert!(metadata.contains("\"implicit_skip_policy\": \"low_probability\""));
    assert!(metadata.contains("mm_question_mark_policy=below_emission_threshold"));
    let qc = fs::read_to_string(root.join("sample.mod.import_qc.tsv")).unwrap();
    assert!(qc.contains("target_groups_unknown\t1\n"));
    assert!(qc.contains("unknown_candidates\t0\n"));
    assert!(qc.contains("mm_question_mark_policy\tbelow_emission_threshold\n"));
    assert!(qc.contains("emitted_explicit_observations\t1\n"));
    assert!(qc.contains("emitted_implicit_observations\t3\n"));
}

#[test]
fn dorado_import_cli_applies_centered_iupac_candidate_rule() {
    use sam::alignment::record_buf::{
        data::field::{value::Array as BufArray, Value as BufValue},
        Cigar, Data, Sequence,
    };
    use sam::header::record::value::{map::ReferenceSequence, Map};

    let root = TestDir::new("dorado-motif-e2e");
    let bam_path = root.join("calls.bam");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let cigar: Cigar = [Op::new(Kind::Match, 5)].into_iter().collect();
    let data: Data = [
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATIONS,
            BufValue::from("A+a.,2;"),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_PROBABILITIES,
            BufValue::Array(BufArray::UInt8(vec![255])),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            BufValue::UInt32(5),
        ),
    ]
    .into_iter()
    .collect();
    let record = sam::alignment::RecordBuf::builder()
        .set_name("motif-read")
        .set_flags(sam::alignment::record::Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start("101".parse().unwrap())
        .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(b"AAACA".to_vec()))
        .set_data(data)
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();

    let prefix = root.join("sample.mod");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "mod-import-dorado",
            "--sample",
            "S1",
            "--assay-id",
            "dorado_rna004_drach",
            "--mod-code",
            "A+a",
            "--model-id",
            "rna004-test-drach",
            "--candidate-rule",
            "drach",
            "--bam",
        ])
        .arg(&bam_path)
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let observations = fs::read_to_string(root.join("sample.mod.observations.tsv")).unwrap();
    assert_eq!(observations.lines().count(), 2);
    assert!(observations.contains(
        "dorado_rna004_drach\tS1\tS1::motif-read\tchr1\t102\t+\tA+a\t0.998046875\texplicit_probability\tDRACH\tNA\tNA\n"
    ));
    let metadata = fs::read_to_string(root.join("sample.mod.assay.json")).unwrap();
    assert!(metadata.contains("\"candidate_rule\": \"DRACH\""));
    let qc = fs::read_to_string(root.join("sample.mod.import_qc.tsv")).unwrap();
    assert!(qc.contains("target_canonical_bases\t4\n"));
    assert!(qc.contains("target_candidate_bases\t1\n"));
    assert!(qc.contains("candidate_observations_complete\ttrue\n"));
}

#[test]
fn dorado_import_cli_treats_sam_n_groups_as_any_base_and_preserves_ml_order() {
    use sam::alignment::record_buf::{
        data::field::{value::Array as BufArray, Value as BufValue},
        Cigar, Data, Sequence,
    };
    use sam::header::record::value::{map::ReferenceSequence, Map};

    let root = TestDir::new("dorado-n-canonical-e2e");
    let bam_path = root.join("calls.bam");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
    let data: Data = [
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATIONS,
            BufValue::from("A+a.,0;N+n.,0;"),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_PROBABILITIES,
            BufValue::Array(BufArray::UInt8(vec![255, 200])),
        ),
        (
            sam::alignment::record::data::field::Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            BufValue::UInt32(4),
        ),
    ]
    .into_iter()
    .collect();
    let record = sam::alignment::RecordBuf::builder()
        .set_name("mixed-groups")
        .set_flags(sam::alignment::record::Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start("101".parse().unwrap())
        .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_data(data)
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();

    let prefix = root.join("sample.mod");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "mod-import-dorado",
            "--sample",
            "S1",
            "--assay-id",
            "dorado_any_base",
            "--mod-code",
            "N+n",
            "--model-id",
            "any-base-test",
            "--source-emission-threshold",
            "0",
            "--bam",
        ])
        .arg(&bam_path)
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let observations = fs::read_to_string(root.join("sample.mod.observations.tsv")).unwrap();
    assert_eq!(observations.lines().count(), 5);
    assert_eq!(observations.matches("\tN+n\t").count(), 4);
    assert!(observations.contains(
        "dorado_any_base\tS1\tS1::mixed-groups\tchr1\t100\t+\tN+n\t0.783203125\texplicit_probability\tNA\tNA\tNA\n"
    ));
    assert_eq!(
        observations
            .matches("\timplicit_below_emission_threshold\t")
            .count(),
        3
    );
    let qc = fs::read_to_string(root.join("sample.mod.import_qc.tsv")).unwrap();
    assert!(qc.contains("target_canonical_bases\t4\n"));
    assert!(qc.contains("target_candidate_bases\t4\n"));
    assert!(qc.contains("ml_values_consumed\t2\n"));
    assert!(qc.contains("explicit_target_calls\t1\n"));
    assert!(qc.contains("emitted_implicit_observations\t3\n"));
}
