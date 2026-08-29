mod common;

use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_sam::alignment::io::Write as _;
use noodles_sam::{
    alignment::{
        record::{
            cigar::{op::Kind as CigarKind, Op as CigarOp},
            Flags, MappingQuality,
        },
        record_buf::{Cigar, QualityScores, Sequence},
        RecordBuf,
    },
    header::record::value::{map::ReferenceSequence, Map},
};

use common::TestDir;

fn temp_dir(label: &str) -> TestDir {
    TestDir::new(&format!("aux-cli-{label}"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn bam_record(name: &str, start1: usize, cigar: Cigar, mapq: u8) -> RecordBuf {
    let read_len = cigar.read_length();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(start1.to_string().parse().unwrap())
        .set_mapping_quality(MappingQuality::try_from(mapq).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; read_len]))
        .set_quality_scores(QualityScores::from(vec![30; read_len]))
        .build()
}

fn match_cigar(len: usize) -> Cigar {
    vec![CigarOp::new(CigarKind::Match, len)]
        .into_iter()
        .collect()
}

fn spliced_cigar(left: usize, intron: usize, right: usize) -> Cigar {
    vec![
        CigarOp::new(CigarKind::Match, left),
        CigarOp::new(CigarKind::Skip, intron),
        CigarOp::new(CigarKind::Match, right),
    ]
    .into_iter()
    .collect()
}

fn deletion_only_middle_block_cigar() -> Cigar {
    vec![
        CigarOp::new(CigarKind::Match, 10),
        CigarOp::new(CigarKind::Skip, 20),
        CigarOp::new(CigarKind::Deletion, 5),
        CigarOp::new(CigarKind::Skip, 20),
        CigarOp::new(CigarKind::Match, 10),
    ]
    .into_iter()
    .collect()
}

fn zero_length_intron_cigar() -> Cigar {
    vec![
        CigarOp::new(CigarKind::Match, 10),
        CigarOp::new(CigarKind::Skip, 0),
        CigarOp::new(CigarKind::Match, 10),
    ]
    .into_iter()
    .collect()
}

#[test]
fn addgene_and_desc_cli_publish_valid_outputs() {
    let root = temp_dir("annotate");
    let annotated = root.join("reads_gene.bed");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let output = Command::new(executable)
        .args(["addgene", "--reads"])
        .arg(fixture("reads.bed"))
        .arg("--reference")
        .arg(fixture("ref.bed"))
        .arg("--out")
        .arg(&annotated)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let records = trackcluster_rs::io::bed::read_bed12(&annotated)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].extra_fields[5], "GENEA");

    let prefix = root.join("catalog");
    let output = Command::new(executable)
        .args(["desc", "--isoform"])
        .arg(fixture("ref.bed"))
        .arg("--reference")
        .arg(fixture("ref.bed"))
        .arg("--out")
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let desc = fs::read_to_string(root.join("catalog_desc.txt")).unwrap();
    assert!(desc.starts_with("#schema\ttrackcluster-description-v2\tdesc\n"));
    assert!(desc.contains("isoform_id\treference_id\tgene_id"));
    assert!(!root.join("catalog_sqanti_structural_category.tsv").exists());
}

#[test]
fn export_cli_writes_standard_formats() {
    let root = temp_dir("interchange");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let gtf = root.join("catalog.gtf");
    let gff3 = root.join("catalog.gff3");
    let sqanti = root.join("catalog.sqanti.tsv");
    let output = Command::new(executable)
        .args(["export", "--input"])
        .arg(fixture("ref.bed"))
        .arg("--gtf")
        .arg(&gtf)
        .arg("--gff3")
        .arg(&gff3)
        .arg("--sqanti-input")
        .arg(&sqanti)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(fs::read_to_string(&gtf).unwrap().contains("\ttranscript\t"));
    assert!(fs::read_to_string(&gff3)
        .unwrap()
        .starts_with("##gff-version 3"));
    assert!(fs::read_to_string(sqanti)
        .unwrap()
        .contains("#schema\ttrackcluster-sqanti-input-v1"));

    let original = trackcluster_rs::io::bed::read_bed12(fixture("ref.bed"))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    for (annotation, format) in [(&gtf, "gtf"), (&gff3, "gff3")] {
        let round_trip = root.join(format!("round-trip.{format}.bed"));
        let output = Command::new(executable)
            .args(["gff2bigg", "--gff"])
            .arg(annotation)
            .arg("--out")
            .arg(&round_trip)
            .args(["--input-format", format])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let converted = trackcluster_rs::io::bed::read_bed12(&round_trip)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(converted.len(), original.len());
        for expected in &original {
            let actual = converted
                .iter()
                .find(|record| record.name == expected.name)
                .unwrap();
            assert_eq!(actual.chrom, expected.chrom);
            assert_eq!(actual.strand, expected.strand);
            assert_eq!(actual.exons, expected.exons);
            assert_eq!(actual.metadata().gene_id(), expected.metadata().gene_id());
        }
    }
}

#[test]
fn standalone_aux_commands_reject_output_aliases_without_clobbering_inputs() {
    let root = temp_dir("output-aliases");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let reference = root.join("reference.bed");
    fs::copy(fixture("ref.bed"), &reference).unwrap();
    let reference_before = fs::read(&reference).unwrap();

    let export = Command::new(executable)
        .args(["export", "--input"])
        .arg(&reference)
        .arg("--gtf")
        .arg(&reference)
        .output()
        .unwrap();
    assert!(!export.status.success());
    assert!(String::from_utf8_lossy(&export.stderr).contains("refer to the same file"));
    assert_eq!(fs::read(&reference).unwrap(), reference_before);

    let shared = root.join("shared-output.txt");
    fs::write(&shared, "previous output\n").unwrap();
    let export = Command::new(executable)
        .args(["export", "--input"])
        .arg(&reference)
        .arg("--gtf")
        .arg(&shared)
        .arg("--gff3")
        .arg(&shared)
        .output()
        .unwrap();
    assert!(!export.status.success());
    assert!(String::from_utf8_lossy(&export.stderr).contains("refer to the same file"));
    assert_eq!(fs::read_to_string(&shared).unwrap(), "previous output\n");

    let reads = root.join("reads.bed");
    fs::copy(fixture("reads.bed"), &reads).unwrap();
    let reads_before = fs::read(&reads).unwrap();
    let addgene = Command::new(executable)
        .args(["addgene", "--reads"])
        .arg(&reads)
        .arg("--reference")
        .arg(&reference)
        .arg("--out")
        .arg(&reads)
        .output()
        .unwrap();
    assert!(!addgene.status.success());
    assert_eq!(fs::read(&reads).unwrap(), reads_before);

    let description_input = root.join("catalog_desc.txt");
    fs::copy(&reference, &description_input).unwrap();
    let description_before = fs::read(&description_input).unwrap();
    let desc = Command::new(executable)
        .args(["desc", "--isoform"])
        .arg(&description_input)
        .arg("--reference")
        .arg(&reference)
        .arg("--out")
        .arg(root.join("catalog"))
        .output()
        .unwrap();
    assert!(!desc.status.success());
    assert_eq!(fs::read(&description_input).unwrap(), description_before);

    let validate = Command::new(executable)
        .args(["validate-bed", "--input"])
        .arg(&reference)
        .arg("--report")
        .arg(&reference)
        .output()
        .unwrap();
    assert!(!validate.status.success());
    assert_eq!(fs::read(&reference).unwrap(), reference_before);
}

#[test]
fn gff2bigg_cli_converts_minimal_gff3_and_gtf() {
    let root = temp_dir("gff2bigg");
    let executable = env!("CARGO_BIN_EXE_trackcluster");

    let gff3 = root.join("annotation.gff3");
    let gff3_bed = root.join("annotation.gff3.bed");
    fs::write(
        &gff3,
        concat!(
            "##gff-version 3\n",
            "chr1\ttest\tgene\t101\t250\t.\t+\t.\tID=gene1;Name=GENEA\n",
            "chr1\ttest\tmRNA\t101\t250\t.\t+\t.\tID=tx1;Parent=gene1\n",
            "chr1\ttest\texon\t101\t120\t.\t+\t.\tParent=tx1\n",
            "chr1\ttest\texon\t201\t250\t.\t+\t.\tParent=tx1\n",
        ),
    )
    .unwrap();
    let output = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&gff3)
        .arg("--out")
        .arg(&gff3_bed)
        .args(["--key", "Name", "--input-format", "gff3"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "gff2bigg: transcripts=1 rejected_records=0 rejected_transcripts=0\n"
    );
    assert_eq!(
        fs::read_to_string(format!("{}.rejected.tsv", gff3_bed.display())).unwrap(),
        "source_path\tanchor_line\ttranscript_ids_json\tkind\treason\n"
    );
    let records = trackcluster_rs::io::bed::read_bed12(&gff3_bed)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name, "tx1");
    assert_eq!(records[0].strand.as_char(), '+');
    assert_eq!(records[0].extra_fields[4], "isoform_anno");
    assert_eq!(records[0].extra_fields[5], "GENEA");
    assert_eq!(
        records[0]
            .exons
            .iter()
            .map(|exon| (exon.start.get(), exon.end.get()))
            .collect::<Vec<_>>(),
        vec![(100, 120), (200, 250)]
    );

    let gtf = root.join("annotation.gtf");
    let gtf_bed = root.join("annotation.gtf.bed");
    fs::write(
        &gtf,
        concat!(
            "chr2\ttest\ttranscript\t301\t450\t.\t-\t.\tgene_id \"GENEB\"; transcript_id \"tx2\";\n",
            "chr2\ttest\texon\t401\t450\t.\t-\t.\tgene_id \"GENEB\"; transcript_id \"tx2\";\n",
            "chr2\ttest\texon\t301\t350\t.\t-\t.\tgene_id \"GENEB\"; transcript_id \"tx2\";\n",
        ),
    )
    .unwrap();
    let output = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&gtf)
        .arg("--out")
        .arg(&gtf_bed)
        .args(["--input-format", "gtf"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let records = trackcluster_rs::io::bed::read_bed12(&gtf_bed)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name, "tx2");
    assert_eq!(records[0].strand.as_char(), '-');
    assert_eq!(records[0].extra_fields[5], "GENEB");
    assert_eq!(
        records[0]
            .exons
            .iter()
            .map(|exon| (exon.start.get(), exon.end.get()))
            .collect::<Vec<_>>(),
        vec![(300, 350), (400, 450)]
    );
}

#[test]
fn gff2bigg_cli_quarantines_identifiable_bad_models_and_audits_them() {
    let root = temp_dir("gff2bigg-recovery");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let gtf = root.join("annotation.gtf");
    let bed = root.join("annotation.bed");
    let rejected = root.join("rejected.tsv");
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"bad,tx\";\n",
            "chr1\ttest\texon\tnot-a-coordinate\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"bad,tx\";\n",
            "chr1\ttest\ttranscript\t201\t250\t.\t+\t.\tgene_id \"G2\"; transcript_id \"good_tx\";\n",
            "chr1\ttest\texon\t201\t250\t.\t+\t.\tgene_id \"G2\"; transcript_id \"good_tx\";\n",
        ),
    )
    .unwrap();

    let output = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&gtf)
        .arg("--out")
        .arg(&bed)
        .args(["--input-format", "gtf", "--rejected-records"])
        .arg(&rejected)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("excluded 1 invalid record(s) affecting 1 transcript model(s)"));
    assert!(stderr.contains("transcripts=1 rejected_records=1 rejected_transcripts=1"));
    let records = trackcluster_rs::io::bed::read_bed12(&bed)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name, "good_tx");
    let audit = fs::read_to_string(&rejected).unwrap();
    assert!(audit.contains("\t2\t[\"bad,tx\"]\tparse\t"), "{audit}");
    assert!(audit.contains("not-a-coordinate"), "{audit}");

    fs::write(&bed, "previous output\n").unwrap();
    let strict = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&gtf)
        .arg("--out")
        .arg(&bed)
        .args(["--input-format", "gtf", "--invalid-record-policy", "fail"])
        .output()
        .unwrap();
    assert!(!strict.status.success());
    assert!(String::from_utf8_lossy(&strict.stderr).contains("not-a-coordinate"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");

    let all_bad = root.join("all-bad.gtf");
    fs::write(
        &all_bad,
        "chr1\ttest\texon\tnot-a-coordinate\t100\t.\t+\t.\tgene_id \"G3\"; transcript_id \"only_bad\";\n",
    )
    .unwrap();
    let recovered = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&all_bad)
        .arg("--out")
        .arg(&bed)
        .args(["--input-format", "gtf"])
        .output()
        .unwrap();
    assert!(!recovered.status.success());
    assert!(String::from_utf8_lossy(&recovered.stderr)
        .contains("annotation contains no valid transcript models"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");
}

#[test]
fn bam2bigg_cli_converts_spliced_bam_and_filters_low_mapq() {
    let root = temp_dir("bam2bigg");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let bam = root.join("reads.bam");
    let bed = root.join("reads.bed");

    let header = noodles_sam::Header::builder()
        .add_reference_sequence(
            "chr1".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build();
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &bam_record("retained", 101, spliced_cigar(10, 20, 10), 60),
        )
        .unwrap();
    writer
        .write_alignment_record(&header, &bam_record("low_mapq", 201, match_cigar(10), 20))
        .unwrap();
    writer.try_finish().unwrap();

    let output = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .args(["--group", "sample-A"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains(
        "bam2bigg: total=2 records=1 skipped_unmapped=0 skipped_secondary=0 skipped_supplementary=0 skipped_below_mapq=1 skipped_invalid=0"
    ));
    let imported = trackcluster_rs::io::bed::read_bed12(&bed)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(imported.len(), 1);
    assert_eq!(imported[0].name, "retained");
    assert_eq!(imported[0].score, 60);
    assert_eq!(imported[0].item_rgb, "250,128,114");
    assert_eq!(imported[0].extra_fields[6], "sample-A");
    assert_eq!(
        imported[0]
            .exons
            .iter()
            .map(|exon| (exon.start.get(), exon.end.get()))
            .collect::<Vec<_>>(),
        vec![(100, 110), (130, 140)]
    );
}

#[test]
fn bam2bigg_skips_only_invalid_decoded_records_and_reports_bounded_diagnostics() {
    let root = temp_dir("bam2bigg-invalid-record");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let bam = root.join("mixed.bam");
    let bed = root.join("mixed.bed");

    let header = noodles_sam::Header::builder()
        .add_reference_sequence(
            "chr1".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build();
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(&header, &bam_record("before", 101, match_cigar(10), 60))
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &bam_record("zero-length", 121, zero_length_intron_cigar(), 60),
        )
        .unwrap();
    for name in ["invalid-1", "invalid-2"] {
        writer
            .write_alignment_record(
                &header,
                &bam_record(name, 151, deletion_only_middle_block_cigar(), 60),
            )
            .unwrap();
    }
    writer
        .write_alignment_record(
            &header,
            &bam_record("invalid-name", 251, match_cigar(10), 60),
        )
        .unwrap();
    writer
        .write_alignment_record(&header, &bam_record("after", 301, match_cigar(10), 60))
        .unwrap();
    writer.try_finish().unwrap();
    drop(writer);

    // The typed BAM writer correctly rejects non-ASCII names, so patch the
    // uncompressed BAM payload to exercise a legacy/malformed producer without
    // changing any record boundary.
    let mut raw_bam = Vec::new();
    let mut decoder = noodles_bgzf::io::Reader::new(fs::File::open(&bam).unwrap());
    std::io::Read::read_to_end(&mut decoder, &mut raw_bam).unwrap();
    let marker = b"invalid-name\0";
    let marker_offsets = raw_bam
        .windows(marker.len())
        .enumerate()
        .filter_map(|(offset, window)| (window == marker).then_some(offset))
        .collect::<Vec<_>>();
    assert_eq!(marker_offsets.len(), 1);
    raw_bam[marker_offsets[0]] = 0xff;
    let mut encoder = noodles_bgzf::io::Writer::new(Vec::new());
    std::io::Write::write_all(&mut encoder, &raw_bam).unwrap();
    fs::write(&bam, encoder.finish().unwrap()).unwrap();

    let output = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(
        "bam2bigg: total=6 records=2 skipped_unmapped=0 skipped_secondary=0 skipped_supplementary=0 skipped_below_mapq=0 skipped_invalid=4"
    ));
    assert!(
        stderr.contains("invalid_record_reason=invalid_cigar_structure records=3 first_record=2")
    );
    assert!(stderr.contains("zero-length CIGAR operation"));
    assert!(stderr.contains("zero-length"));
    assert!(!stderr.contains("invalid-1"));
    assert!(!stderr.contains("invalid-2"));
    assert!(stderr.contains("invalid_record_reason=invalid_query_name records=1 first_record=5"));
    assert!(stderr.contains("non-UTF-8 query name"));
    assert!(!stderr.contains('\u{fffd}'));

    let records = trackcluster_rs::io::bed::read_bed12(&bed)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(
        records
            .iter()
            .map(|record| record.name.as_str())
            .collect::<Vec<_>>(),
        vec!["before", "after"]
    );

    fs::write(&bed, "previous output\n").unwrap();
    let strict = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .args(["--invalid-record-policy", "fail"])
        .output()
        .unwrap();
    assert!(!strict.status.success());
    assert!(String::from_utf8_lossy(&strict.stderr).contains("convert BAM record 2"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");
}

#[test]
fn bam2bigg_does_not_recover_from_a_truncated_record_stream() {
    let root = temp_dir("bam2bigg-truncated-record");
    let executable = env!("CARGO_BIN_EXE_trackcluster");
    let bam = root.join("truncated.bam");
    let bed = root.join("existing.bed");

    let header = noodles_sam::Header::builder()
        .add_reference_sequence(
            "chr1".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build();
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    std::io::Write::flush(writer.get_mut()).unwrap();
    writer
        .write_alignment_record(&header, &bam_record("valid", 101, match_cigar(10), 60))
        .unwrap();
    std::io::Write::flush(writer.get_mut()).unwrap();
    let complete_prefix_len = writer.get_ref().position();
    writer
        .write_alignment_record(&header, &bam_record("truncated", 201, match_cigar(10), 60))
        .unwrap();
    std::io::Write::flush(writer.get_mut()).unwrap();
    let second_record_end = writer.get_ref().position();
    writer.try_finish().unwrap();
    drop(writer);

    assert!(second_record_end > complete_prefix_len + 18);
    let truncated_len = second_record_end - 1;
    assert!(fs::metadata(&bam).unwrap().len() > truncated_len);
    fs::OpenOptions::new()
        .write(true)
        .open(&bam)
        .unwrap()
        .set_len(truncated_len)
        .unwrap();
    fs::write(&bed, "previous output\n").unwrap();

    let output = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("read BAM record 2"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");
}

#[test]
fn converter_failures_preserve_inputs_and_previous_outputs() {
    let root = temp_dir("converter-failures");
    let executable = env!("CARGO_BIN_EXE_trackcluster");

    let gff3 = root.join("same.gff3");
    let gff3_text = concat!(
        "##gff-version 3\n",
        "chr1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=tx1\n",
        "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=tx1\n",
    );
    fs::write(&gff3, gff3_text).unwrap();
    let alias = root.join(".").join("same.gff3");
    let output = Command::new(executable)
        .args(["gff2bigg", "--gff"])
        .arg(&gff3)
        .arg("--out")
        .arg(&alias)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
    assert_eq!(fs::read_to_string(&gff3).unwrap(), gff3_text);

    let bam = root.join("invalid.bam");
    let bed = root.join("existing.bed");
    let header = noodles_sam::Header::builder()
        .add_reference_sequence(
            "chr1".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build();
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &bam_record("invalid", 101, deletion_only_middle_block_cigar(), 60),
        )
        .unwrap();
    writer.try_finish().unwrap();
    fs::write(&bed, "previous output\n").unwrap();

    let output = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("produced no valid records after skipping 1 invalid decoded record"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");

    let output = Command::new(executable)
        .args(["bam2bigg", "--bamfile"])
        .arg(&bam)
        .arg("--out")
        .arg(&bed)
        .args(["--invalid-record-policy", "fail"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("only deletions"));
    assert_eq!(fs::read_to_string(&bed).unwrap(), "previous output\n");
}
