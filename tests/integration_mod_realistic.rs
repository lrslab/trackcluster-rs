mod common;

use std::collections::BTreeMap;
use std::fs;
use std::num::NonZero;
use std::path::Path;
use std::process::Command;

use common::TestDir;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;
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

const READS_PER_ISOFORM: usize = 40;

#[derive(Clone, Copy)]
struct SamplePlan {
    sample: &'static str,
    group: &'static str,
    iso1_modified: usize,
    iso2_modified: usize,
}

const SAMPLE_PLANS: [SamplePlan; 6] = [
    SamplePlan {
        sample: "C1",
        group: "control",
        iso1_modified: 10,
        iso2_modified: 8,
    },
    SamplePlan {
        sample: "C2",
        group: "control",
        iso1_modified: 12,
        iso2_modified: 8,
    },
    SamplePlan {
        sample: "C3",
        group: "control",
        iso1_modified: 8,
        iso2_modified: 8,
    },
    SamplePlan {
        sample: "T1",
        group: "treated",
        iso1_modified: 30,
        iso2_modified: 16,
    },
    SamplePlan {
        sample: "T2",
        group: "treated",
        iso1_modified: 28,
        iso2_modified: 18,
    },
    SamplePlan {
        sample: "T3",
        group: "treated",
        iso1_modified: 32,
        iso2_modified: 14,
    },
];

fn observation(sample: &str, raw_read_id: &str, probability: f64) -> ModObservation {
    ModObservation {
        key: ModObservationKey {
            assay_id: "realistic_a1".to_owned(),
            sample: sample.to_owned(),
            read_id: format!("{sample}::{raw_read_id}"),
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

fn write_spliced_coverage_bam(path: &Path, sample: &str, reads: &[(&str, &str)]) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for &(raw_read_id, isoform) in reads {
        let skipped = match isoform {
            "iso1" => 80,
            "iso2" => 180,
            _ => panic!("unexpected isoform {isoform}"),
        };
        let cigar: Cigar = [
            Op::new(Kind::Match, 20),
            Op::new(Kind::Skip, skipped),
            Op::new(Kind::Match, 20),
        ]
        .into_iter()
        .collect();
        let record = sam::alignment::RecordBuf::builder()
            // Already-tagged names exercise synchronized TrackCluster/coverage inputs.
            .set_name(format!("{sample}::{raw_read_id}"))
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start("101".parse().unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 40]))
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn read_tsv(path: &Path) -> Vec<BTreeMap<String, String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let headers = reader.headers().unwrap().clone();
    reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            headers
                .iter()
                .zip(record.iter())
                .map(|(key, value)| (key.to_owned(), value.to_owned()))
                .collect()
        })
        .collect()
}

fn value_as_f64(row: &BTreeMap<String, String>, field: &str) -> f64 {
    row[field].parse().unwrap()
}

#[test]
fn multi_sample_isoform_condition_effects_recover_known_integer_truth() {
    let root = TestDir::new("modification-realistic");
    let isoforms = root.join("isoforms.bed");
    fs::write(
        &isoforms,
        concat!(
            "chr1\t100\t220\tiso1\t0\t+\t100\t220\t0\t2\t20,20,\t0,100,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
            "chr1\t100\t320\tiso2\t0\t+\t100\t320\t0\t2\t20,20,\t0,200,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
        ),
    )
    .unwrap();

    let reference_fasta = root.join("reference.fa");
    fs::write(&reference_fasta, format!(">chr1\n{}\n", "A".repeat(1000))).unwrap();
    fs::write(
        Path::new(&format!("{}.fai", reference_fasta.display())),
        "chr1\t1000\t6\t1000\t1001\n",
    )
    .unwrap();

    let metadata = AssayMetadata {
        schema_version: 1,
        assay_id: "realistic_a1".to_owned(),
        caller: "deterministic_simulation".to_owned(),
        caller_version: "1".to_owned(),
        model_id: "known_integer_truth".to_owned(),
        chemistry: "RNA004".to_owned(),
        candidate_rule: "DRACH".to_owned(),
        source_emission_threshold: None,
        source_site_filter: "shared_exon_site".to_owned(),
        candidate_observations_complete: true,
        implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
        coordinate_source: "synthetic_genomic_truth".to_owned(),
        read_id_mapping: "sample_prefixed_simulated_read".to_owned(),
        source_files: Vec::new(),
    };
    let assay_path = root.join("assay.json");
    write_assay_metadata_to_writer(fs::File::create(&assay_path).unwrap(), &metadata).unwrap();

    let mut sample_manifest = "sample\tgroup\treads\n".to_owned();
    let mut mod_manifest =
        "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n".to_owned();
    let mut mapping = String::new();
    let mut expected_modified = BTreeMap::new();

    for plan in SAMPLE_PLANS {
        let reads_path = root.join(format!("{}.reads.bed", plan.sample));
        let observations_path = root.join(format!("{}.observations.tsv", plan.sample));
        let coverage_path = root.join(format!("{}.coverage.bam", plan.sample));
        let mut reads_bed = String::new();
        let mut observations = Vec::new();
        let mut coverage_reads = Vec::new();

        for (isoform, modified) in [("iso1", plan.iso1_modified), ("iso2", plan.iso2_modified)] {
            expected_modified.insert(
                (plan.sample.to_owned(), isoform.to_owned()),
                u64::try_from(modified).unwrap(),
            );
            for read_index in 0..READS_PER_ISOFORM {
                let raw_read_id = format!("{}_{}_{read_index:03}", plan.sample, isoform);
                let (end, block_start) = if isoform == "iso1" {
                    (220, 100)
                } else {
                    (320, 200)
                };
                reads_bed.push_str(&format!(
                    "chr1\t100\t{end}\t{raw_read_id}\t0\t+\t100\t{end}\t0\t2\t20,20,\t0,{block_start},\n"
                ));
                mapping.push_str(&format!("{}::{}\t{isoform}\n", plan.sample, raw_read_id));
                observations.push(observation(
                    plan.sample,
                    &raw_read_id,
                    if read_index < modified { 0.9 } else { 0.1 },
                ));
                coverage_reads.push((raw_read_id, isoform));
            }
        }

        // Real callers contain reads that cannot be assigned to a final isoform.
        observations.push(observation(plan.sample, "unassigned_high", 0.9));
        observations.push(observation(plan.sample, "unassigned_low", 0.1));
        fs::write(&reads_path, reads_bed).unwrap();
        write_observations_tsv_to_writer(
            fs::File::create(&observations_path).unwrap(),
            &observations,
        )
        .unwrap();
        let coverage_refs = coverage_reads
            .iter()
            .map(|(read, isoform)| (read.as_str(), *isoform))
            .collect::<Vec<_>>();
        write_spliced_coverage_bam(&coverage_path, plan.sample, &coverage_refs);

        sample_manifest.push_str(&format!(
            "{}\t{}\t{}\n",
            plan.sample,
            plan.group,
            reads_path.display()
        ));
        mod_manifest.push_str(&format!(
            "{}\trealistic_a1\t{}\t{}\t{}\n",
            plan.sample,
            observations_path.display(),
            assay_path.display(),
            coverage_path.display()
        ));
    }

    let sample_manifest_path = root.join("samples.tsv");
    let mod_manifest_path = root.join("mod_samples.tsv");
    let mapping_path = root.join("read_to_isoform.unique.tsv");
    fs::write(&sample_manifest_path, sample_manifest).unwrap();
    fs::write(&mod_manifest_path, mod_manifest).unwrap();
    fs::write(&mapping_path, mapping).unwrap();

    let aggregate_prefix = root.join("result");
    let aggregate = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-aggregate", "--manifest"])
        .arg(&sample_manifest_path)
        .arg("--isoforms")
        .arg(&isoforms)
        .arg("--read-to-isoform")
        .arg(&mapping_path)
        .arg("--mod-manifest")
        .arg(&mod_manifest_path)
        .arg("--reference-fasta")
        .arg(&reference_fasta)
        .args([
            "--analysis-threshold",
            "realistic_a1=0.5",
            "--min-callable",
            "20",
            "--min-read-join-rate",
            "0.95",
        ])
        .arg("--out")
        .arg(&aggregate_prefix)
        .output()
        .unwrap();
    assert!(
        aggregate.status.success(),
        "{}",
        String::from_utf8_lossy(&aggregate.stderr)
    );

    let join_qc = read_tsv(&root.join("result.mod_join_qc.tsv"));
    assert_eq!(join_qc.len(), 6);
    for row in &join_qc {
        assert_eq!(row["input_rows"], "82");
        assert_eq!(row["joined_rows"], "80");
        assert_eq!(row["joined_reads"], "80");
        assert_eq!(row["unknown_read"], "2");
        assert!(value_as_f64(row, "read_join_rate") > 0.97);
        assert_eq!(row["candidate_observations_complete"], "true");
    }

    let sites = read_tsv(&root.join("result.isoform_mod_sites.tsv"));
    assert_eq!(sites.len(), 12);
    for row in &sites {
        let key = (row["sample"].clone(), row["isoform_id"].clone());
        let modified = expected_modified[&key];
        assert_eq!(row["coverage_basis"], "bam_exact");
        assert_eq!(row["n_assigned"], "40");
        assert_eq!(row["n_covering"], "40");
        assert_eq!(row["n_candidate"], "40");
        assert_eq!(row["n_callable"], "40");
        assert_eq!(row["n_modified"], modified.to_string());
        assert_eq!(row["n_unmodified"], (40 - modified).to_string());
        assert_eq!(row["n_unknown"], "0");
        assert_eq!(row["eligibility_reason"], "ok");
        assert!((value_as_f64(row, "mod_fraction") - modified as f64 / 40.0).abs() < 1e-12);
    }

    let contrast_spec = root.join("contrasts.tsv");
    fs::write(
        &contrast_spec,
        concat!(
            "contrast_type\tassay_id\tgene\tsite_id\tmod_code\tisoform_a\tisoform_b\tgroup_a\tgroup_b\n",
            "isoform_effect\trealistic_a1\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\tNA\n",
            "condition_effect\trealistic_a1\tG1\tchr1:110:+\tA+a\tiso1\tNA\tcontrol\ttreated\n",
            "isoform_condition_interaction\trealistic_a1\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tcontrol\ttreated\n",
        ),
    )
    .unwrap();
    let contrast_prefix = root.join("effects");
    let contrast = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-contrast", "--design"])
        .arg(root.join("result.isoform_mod_design.tsv"))
        .arg("--contrasts")
        .arg(&contrast_spec)
        .arg("--out")
        .arg(&contrast_prefix)
        .output()
        .unwrap();
    assert!(
        contrast.status.success(),
        "{}",
        String::from_utf8_lossy(&contrast.stderr)
    );

    let contrasts = read_tsv(&root.join("effects.isoform_mod_contrasts.tsv"));
    assert_eq!(contrasts.len(), 3);
    for row in &contrasts {
        assert_eq!(row["p_value"], "NA");
        assert_eq!(row["q_value"], "NA");
        assert_eq!(row["method"], "effect_only");
        assert_eq!(row["eligibility_reason"], "ok");
        match row["contrast_type"].as_str() {
            "isoform_effect" => {
                assert_eq!(row["n_eligible_samples"], "3");
                assert!((value_as_f64(row, "delta_fraction") - 0.05).abs() < 1e-12);
            }
            "condition_effect" => {
                assert_eq!(row["n_eligible_samples"], "6");
                assert!((value_as_f64(row, "delta_fraction") - 0.5).abs() < 1e-12);
            }
            "isoform_condition_interaction" => {
                assert_eq!(row["n_eligible_samples"], "6");
                assert!((value_as_f64(row, "interaction_delta") - 0.3).abs() < 1e-12);
            }
            other => panic!("unexpected contrast type {other}"),
        }
    }
}
