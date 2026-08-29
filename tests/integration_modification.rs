mod common;

use std::fs;
use std::num::NonZero;
use std::path::Path;
use std::process::Command;

use common::TestDir;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write;
use sam::alignment::record::cigar::{op::Kind, Op};
use sam::alignment::record_buf::{Cigar, Sequence};
use sam::header::record::value::{map::ReferenceSequence, Map};
use trackcluster_rs::io::mod_calls::{
    write_assay_metadata_to_writer, write_observations_tsv_to_writer,
};
use trackcluster_rs::model::Strand;
use trackcluster_rs::modification::{
    AssayMetadata, ImplicitSkipPolicy, ModObservation, ModObservationKey, ModSiteKey,
    ObservationState,
};

fn observation(sample: &str, read: &str, probability: f64) -> ModObservation {
    ModObservation {
        key: ModObservationKey {
            assay_id: "a1".to_owned(),
            sample: sample.to_owned(),
            read_id: format!("{sample}::{read}"),
            site: ModSiteKey {
                chrom: "chr1".to_owned(),
                pos0: 110,
                strand: Strand::Plus,
                mod_code: "A+a".to_owned(),
            },
        },
        probability: Some(probability),
        observation_state: ObservationState::ExplicitProbability,
        context: Some("DRACH".to_owned()),
        source_transcript_id: None,
        source_pos0: None,
    }
}

fn write_coverage_bam(path: &Path, read_names: &[&str]) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for name in read_names {
        let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_name(*name)
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start("101".parse().unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 20]))
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_reference_fasta(path: &Path) {
    fs::write(path, format!(">chr1\n{}\n", "A".repeat(1000))).unwrap();
    fs::write(
        Path::new(&format!("{}.fai", path.display())),
        "chr1\t1000\t6\t1000\t1001\n",
    )
    .unwrap();
}

#[test]
fn aggregate_and_effect_only_contrast_run_end_to_end() {
    let root = TestDir::new("modification-e2e");
    let isoforms = root.join("isoforms.bed");
    fs::write(
        &isoforms,
        concat!(
            "chr1\t100\t220\tiso1\t0\t+\t100\t220\t0\t2\t20,20,\t0,100,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
            "chr1\t100\t320\tiso2\t0\t+\t100\t320\t0\t2\t20,20,\t0,200,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
        ),
    )
    .unwrap();

    let s1_reads = root.join("S1.reads.bed");
    let s2_reads = root.join("S2.reads.bed");
    fs::write(
        &s1_reads,
        "chr1\t100\t120\tr1\t0\t+\t100\t120\t0\t1\t20,\t0,\n",
    )
    .unwrap();
    fs::write(
        &s2_reads,
        "chr1\t100\t120\tr3\t0\t+\t100\t120\t0\t1\t20,\t0,\n",
    )
    .unwrap();
    let samples = root.join("samples.tsv");
    fs::write(
        &samples,
        format!(
            "sample\tgroup\treads\nS1\tcontrol\t{}\nS2\ttreated\t{}\n",
            s1_reads.display(),
            s2_reads.display()
        ),
    )
    .unwrap();

    let mapping = root.join("read_to_isoform.unique.tsv");
    fs::write(
        &mapping,
        concat!(
            "S1::r1\tiso1\n",
            "S1::r2\tiso2\n",
            "S2::r3\tiso1\n",
            "S2::r4\tiso2\n",
        ),
    )
    .unwrap();

    let metadata = AssayMetadata {
        schema_version: 1,
        assay_id: "a1".to_owned(),
        caller: "test".to_owned(),
        caller_version: "1".to_owned(),
        model_id: "test-model".to_owned(),
        chemistry: "RNA004".to_owned(),
        candidate_rule: "DRACH".to_owned(),
        source_emission_threshold: None,
        source_site_filter: "none".to_owned(),
        candidate_observations_complete: true,
        provenance_status: trackcluster_rs::modification::ProvenanceStatus::UserDeclared,
        implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
        coordinate_source: "synthetic_genomic".to_owned(),
        read_id_mapping: "explicit".to_owned(),
        source_files: Vec::new(),
    };
    let assay = root.join("a1.assay.json");
    write_assay_metadata_to_writer(fs::File::create(&assay).unwrap(), &metadata).unwrap();
    let s1_observations = root.join("S1.observations.tsv");
    write_observations_tsv_to_writer(
        fs::File::create(&s1_observations).unwrap(),
        &[observation("S1", "r1", 0.9), observation("S1", "r2", 0.1)],
    )
    .unwrap();
    let s2_observations = root.join("S2.observations.tsv");
    write_observations_tsv_to_writer(
        fs::File::create(&s2_observations).unwrap(),
        &[observation("S2", "r3", 0.4), observation("S2", "r4", 0.3)],
    )
    .unwrap();
    let s1_coverage = root.join("S1.coverage.bam");
    let s2_coverage = root.join("S2.coverage.bam");
    write_coverage_bam(&s1_coverage, &["r1", "r2"]);
    write_coverage_bam(&s2_coverage, &["r3", "r4"]);
    let mod_manifest = root.join("mod_samples.tsv");
    fs::write(
        &mod_manifest,
        format!(
            concat!(
                "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n",
                "S1\ta1\t{}\t{}\t{}\n",
                "S2\ta1\t{}\t{}\t{}\n",
            ),
            s1_observations.display(),
            assay.display(),
            s1_coverage.display(),
            s2_observations.display(),
            assay.display(),
            s2_coverage.display(),
        ),
    )
    .unwrap();
    let reference_fasta = root.join("reference.fa");
    write_reference_fasta(&reference_fasta);

    let prefix = root.join("result");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-aggregate", "--manifest"])
        .arg(&samples)
        .arg("--isoforms")
        .arg(&isoforms)
        .arg("--read-to-isoform")
        .arg(&mapping)
        .arg("--mod-manifest")
        .arg(&mod_manifest)
        .arg("--reference-fasta")
        .arg(&reference_fasta)
        .args([
            "--analysis-threshold",
            "a1=0.5",
            "--min-read-join-rate",
            "1",
        ])
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
        "mod-aggregate: assays=1 join_qc_rows=2 site_join_qc_rows=2 site_rows=4 design_rows=4\n"
    );

    let join_qc = fs::read_to_string(root.join("result.mod_join_qc.tsv")).unwrap();
    assert_eq!(join_qc.lines().count(), 3);
    for line in join_qc.lines().skip(1) {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(
            &fields[3..],
            &["2", "2", "2", "2", "2", "2", "1", "1", "0", "0", "0", "0", "0", "0", "0", "true"],
            "{line}"
        );
    }
    let site_join_qc = fs::read_to_string(root.join("result.mod_site_join_qc.tsv")).unwrap();
    assert_eq!(site_join_qc.lines().count(), 3);
    assert!(site_join_qc
        .lines()
        .skip(1)
        .all(|line| line.ends_with("\ttrue")));
    let sites = fs::read_to_string(root.join("result.isoform_mod_sites.tsv")).unwrap();
    assert_eq!(sites.lines().count(), 5);
    assert!(sites.contains("\tbam_exact\texploratory\t1\t1\t0\t1\t1\t1\t1\t1\t0\t0\t1\t0.9"));
    let summary_output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-site-summary", "--sites"])
        .arg(root.join("result.isoform_mod_sites.tsv"))
        .arg("--out")
        .arg(root.join("summary"))
        .output()
        .unwrap();
    assert!(
        summary_output.status.success(),
        "{}",
        String::from_utf8_lossy(&summary_output.stderr)
    );
    let design = root.join("result.isoform_mod_design.tsv");
    let design_text = fs::read_to_string(&design).unwrap();
    assert!(design_text.starts_with(
        "assay_id\tanalysis_threshold\tsample\tgroup\tgene\tsite_id\tmod_code\tisoform_id"
    ));

    let contrasts = root.join("contrasts.tsv");
    fs::write(
        &contrasts,
        concat!(
            "contrast_type\tassay_id\tgene\tsite_id\tmod_code\tisoform_a\tisoform_b\tgroup_a\tgroup_b\n",
            "isoform_effect\ta1\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\tNA\n",
            "condition_effect\ta1\tG1\tchr1:110:+\tA+a\tiso1\tNA\tcontrol\ttreated\n",
            "isoform_condition_interaction\ta1\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\ttreated\n",
        ),
    )
    .unwrap();
    let contrast_prefix = root.join("effects");
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-contrast", "--design"])
        .arg(&design)
        .arg("--contrasts")
        .arg(&contrasts)
        .arg("--out")
        .arg(&contrast_prefix)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let effects = fs::read_to_string(root.join("effects.isoform_mod_contrasts.tsv")).unwrap();
    assert_eq!(effects.lines().count(), 4);
    assert!(effects
        .lines()
        .skip(1)
        .all(|line| line.contains("\tNA\tNA\teffect_only\tok")));
    assert!(effects
        .contains("isoform_effect\ta1\t0.5\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\tNA\t1\t1"));
    assert!(effects.contains(
        "condition_effect\ta1\t0.5\tG1\tchr1:110:+\tA+a\tiso1\tNA\tcontrol\ttreated\t2\t-1"
    ));
    assert!(effects.contains("isoform_condition_interaction\ta1\t0.5\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\ttreated\t2\tNA\tNA\t-1"));
}
