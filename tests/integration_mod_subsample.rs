mod common;

use std::fs;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{assert_success, TestDir};
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;
use sam::alignment::record::cigar::{op::Kind, Op};
use sam::alignment::record::data::field::Tag;
use sam::alignment::record_buf::{
    data::field::{value::Array as BufArray, Value as BufValue},
    Cigar, Data, Sequence,
};
use sam::header::record::value::{map::ReferenceSequence, Map};
use trackcluster_rs::io::mod_calls::{
    write_assay_metadata_to_writer, write_observations_tsv_to_writer,
};
use trackcluster_rs::model::Strand;
use trackcluster_rs::modification::{
    AssayMetadata, ImplicitSkipPolicy, ModObservation, ModObservationKey, ModSiteKey,
    ObservationState,
};

struct Fixture {
    samples: PathBuf,
    isoforms: PathBuf,
    mapping: PathBuf,
    mod_manifest: PathBuf,
}

fn observation(read_name: &str, modified: bool) -> ModObservation {
    ModObservation {
        key: ModObservationKey {
            assay_id: "a1".to_owned(),
            sample: "parent".to_owned(),
            read_id: format!("parent::{read_name}"),
            site: ModSiteKey {
                chrom: "chr1".to_owned(),
                pos0: 110,
                strand: Strand::Plus,
                mod_code: "A+a".to_owned(),
            },
        },
        probability: Some(if modified { 0.9 } else { 0.1 }),
        observation_state: ObservationState::ExplicitProbability,
        context: Some("DRACH".to_owned()),
        source_transcript_id: None,
        source_pos0: None,
    }
}

fn write_modbam(path: &Path, read_names: &[String]) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
        )
        .build();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for read_name in read_names {
        let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
        let data: Data = [
            (Tag::BASE_MODIFICATIONS, BufValue::from("A+a.,0;")),
            (
                Tag::BASE_MODIFICATION_PROBABILITIES,
                BufValue::Array(BufArray::UInt8(vec![255])),
            ),
            (Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, BufValue::UInt32(20)),
        ]
        .into_iter()
        .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_name(read_name.as_str())
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start("101".parse().unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 20]))
            .set_data(data)
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_fixture(root: &Path, read_count: usize, assigned_count: usize) -> Fixture {
    assert!(assigned_count <= read_count);
    let isoforms = root.join("isoforms.bed");
    fs::write(
        &isoforms,
        concat!(
            "chr1\t100\t220\tiso1\t0\t+\t100\t220\t0\t2\t20,20,\t0,100,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
            "chr1\t100\t320\tiso2\t0\t+\t100\t320\t0\t2\t20,20,\t0,200,\tnone\tnone\tnone\t-1,\tisoform\tG1\tnone\tnone\n",
        ),
    )
    .unwrap();

    let read_names = (0..read_count)
        .map(|index| format!("read-{index:03}"))
        .collect::<Vec<_>>();
    let reads = root.join("parent.reads.bed");
    let reads_text = read_names
        .iter()
        .map(|read_name| format!("chr1\t100\t120\t{read_name}\t0\t+\t100\t120\t0\t1\t20,\t0,\n"))
        .collect::<String>();
    fs::write(&reads, reads_text).unwrap();
    let samples = root.join("samples.tsv");
    fs::write(
        &samples,
        format!("sample\tgroup\treads\nparent\tUHRR\t{}\n", reads.display()),
    )
    .unwrap();

    let mapping = root.join("read_to_isoform.unique.tsv");
    let mapping_text = read_names
        .iter()
        .take(assigned_count)
        .enumerate()
        .map(|(index, read_name)| {
            let isoform = if index % 2 == 0 { "iso1" } else { "iso2" };
            format!("parent::{read_name}\t{isoform}\n")
        })
        .collect::<String>();
    fs::write(&mapping, mapping_text).unwrap();

    let observations = root.join("parent.observations.tsv");
    let observation_rows = read_names
        .iter()
        .enumerate()
        .map(|(index, read_name)| observation(read_name, index % 2 == 0))
        .collect::<Vec<_>>();
    write_observations_tsv_to_writer(fs::File::create(&observations).unwrap(), &observation_rows)
        .unwrap();
    let metadata = AssayMetadata {
        schema_version: 1,
        assay_id: "a1".to_owned(),
        caller: "dorado".to_owned(),
        caller_version: "test".to_owned(),
        model_id: "test-model".to_owned(),
        chemistry: "RNA004".to_owned(),
        candidate_rule: "DRACH".to_owned(),
        source_emission_threshold: None,
        source_site_filter: "none".to_owned(),
        candidate_observations_complete: true,
        provenance_status: trackcluster_rs::modification::ProvenanceStatus::UserDeclared,
        implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
        coordinate_source: "synthetic_genomic".to_owned(),
        read_id_mapping: "sample_prefix".to_owned(),
        source_files: Vec::new(),
    };
    let assay_metadata = root.join("parent.assay.json");
    write_assay_metadata_to_writer(fs::File::create(&assay_metadata).unwrap(), &metadata).unwrap();
    let coverage_bam = root.join("parent.modbam.bam");
    write_modbam(&coverage_bam, &read_names);
    let mod_manifest = root.join("mod_samples.tsv");
    fs::write(
        &mod_manifest,
        format!(
            concat!(
                "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n",
                "parent\ta1\t{}\t{}\t{}\n"
            ),
            observations.display(),
            assay_metadata.display(),
            coverage_bam.display()
        ),
    )
    .unwrap();

    Fixture {
        samples,
        isoforms,
        mapping,
        mod_manifest,
    }
}

fn run_subsample(fixture: &Fixture, out_dir: &Path, reads_per_sample: usize) {
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-subsample", "--manifest"])
        .arg(&fixture.samples)
        .arg("--read-to-isoform")
        .arg(&fixture.mapping)
        .arg("--mod-manifest")
        .arg(&fixture.mod_manifest)
        .args([
            "--source-sample",
            "parent",
            "--sample-prefix",
            "pseudo",
            "--replicates",
            "2",
            "--reads-per-sample",
            &reads_per_sample.to_string(),
            "--mode",
            "disjoint",
            "--seed",
            "17",
            "--out-dir",
        ])
        .arg(out_dir)
        .output()
        .unwrap();
    assert_success(&output, "mod-subsample");
}

fn run_aggregate(fixture: &Fixture, bundle: &Path, prefix: &Path, allow_low_join: bool) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    command
        .args(["mod-aggregate", "--manifest"])
        .arg(bundle.join("samples.tsv"))
        .arg("--isoforms")
        .arg(&fixture.isoforms)
        .arg("--read-to-isoform")
        .arg(bundle.join("read_to_isoform.unique.tsv"))
        .arg("--mod-manifest")
        .arg(bundle.join("mod_samples.tsv"))
        .args([
            "--analysis-threshold",
            "a1=0.5",
            "--min-callable",
            "1",
            "--min-read-join-rate",
            "1",
        ]);
    if allow_low_join {
        command.arg("--allow-low-join");
    }
    let output = command.arg("--out").arg(prefix).output().unwrap();
    assert_success(&output, "aggregate pseudo-samples");
}

fn tsv_rows(path: &Path) -> Vec<csv::StringRecord> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap()
        .records()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
}

#[test]
fn disjoint_pseudo_samples_run_through_isoform_quantification_and_contrast() {
    let root = TestDir::new("mod-subsample-quant");
    let fixture = write_fixture(root.path(), 20, 20);
    let bundle = root.join("bundle");
    run_subsample(&fixture, &bundle, 10);

    let selected = fs::read_to_string(bundle.join("subsample_read_ids.tsv")).unwrap();
    assert_eq!(selected.lines().count(), 21);
    assert!(selected
        .lines()
        .skip(1)
        .all(|line| !line.contains("\tNA\t")));
    let samples = fs::read_to_string(bundle.join("samples.tsv")).unwrap();
    assert!(samples.contains("pseudo_001\t\tsamples/pseudo_001.reads.bed"));
    let provenance = fs::read_to_string(bundle.join("subsample_provenance.json")).unwrap();
    assert!(provenance.contains("\"biological_replicates\": false"));
    assert!(provenance.contains("\"parent_group\": \"UHRR\""));
    let overlap = fs::read_to_string(bundle.join("overlap_qc.tsv")).unwrap();
    assert!(overlap.contains("pseudo_001\tpseudo_002\t0\t20\t0"));
    let checksums = fs::read_to_string(bundle.join("SHA256SUMS")).unwrap();
    assert!(checksums.contains("  subsample_provenance.json\n"));
    assert!(!checksums.contains("  SHA256SUMS\n"));

    let reimport = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "mod-import-dorado",
            "--sample",
            "pseudo_001",
            "--assay-id",
            "reimport",
            "--model-id",
            "test-model",
            "--source-emission-threshold",
            "0.05",
            "--bam",
        ])
        .arg(bundle.join("coverage/pseudo_001.assay_001.bam"))
        .arg("--out")
        .arg(root.join("reimport"))
        .output()
        .unwrap();
    assert_success(&reimport, "reimport filtered modBAM");
    assert_eq!(
        fs::read_to_string(root.join("reimport.observations.tsv"))
            .unwrap()
            .lines()
            .count(),
        201
    );

    let repeat_bundle = root.join("bundle-repeat");
    run_subsample(&fixture, &repeat_bundle, 10);
    assert_eq!(
        fs::read(bundle.join("subsample_read_ids.tsv")).unwrap(),
        fs::read(repeat_bundle.join("subsample_read_ids.tsv")).unwrap()
    );

    let aggregate_prefix = root.join("quant");
    run_aggregate(&fixture, &bundle, &aggregate_prefix, false);
    let join_rows = tsv_rows(&root.join("quant.mod_join_qc.tsv"));
    assert_eq!(join_rows.len(), 2);
    assert!(join_rows.iter().all(|row| row.get(10) == Some("1")));

    let site_rows = tsv_rows(&root.join("quant.isoform_mod_sites.tsv"));
    assert_eq!(site_rows.len(), 4);
    assert!(site_rows
        .iter()
        .all(|row| row.get(30) == Some("ok") && row.get(29) == Some("eligible")));

    let contrasts = root.join("contrasts.tsv");
    fs::write(
        &contrasts,
        concat!(
            "contrast_type\tassay_id\tgene\tsite_id\tmod_code\tisoform_a\tisoform_b\tgroup_a\tgroup_b\n",
            "isoform_effect\ta1\tG1\tchr1:110:+\tA+a\tiso1\tiso2\tNA\tNA\n",
        ),
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-contrast", "--design"])
        .arg(root.join("quant.isoform_mod_design.tsv"))
        .arg("--contrasts")
        .arg(&contrasts)
        .arg("--out")
        .arg(root.join("effect"))
        .output()
        .unwrap();
    assert_success(&output, "contrast pseudo-samples");
    let effects = tsv_rows(&root.join("effect.isoform_mod_contrasts.tsv"));
    assert_eq!(effects.len(), 1);
    assert_eq!(effects[0].get(10), Some("2"));
    assert_eq!(effects[0].get(11), Some("1"));
    assert_eq!(effects[0].get(17), Some("ok"));
}

#[test]
fn selected_unassigned_observations_remain_visible_to_the_join_gate() {
    let root = TestDir::new("mod-subsample-unassigned");
    let fixture = write_fixture(root.path(), 4, 3);
    let bundle = root.join("bundle");
    run_subsample(&fixture, &bundle, 2);

    let qc_rows = tsv_rows(&bundle.join("subsample_qc.tsv"));
    assert_eq!(qc_rows.len(), 2);
    assert_eq!(
        qc_rows
            .iter()
            .map(|row| row.get(9).unwrap().parse::<usize>().unwrap())
            .sum::<usize>(),
        3
    );
    assert_eq!(
        qc_rows
            .iter()
            .map(|row| row.get(10).unwrap().parse::<usize>().unwrap())
            .sum::<usize>(),
        4
    );

    let prefix = root.join("quant");
    run_aggregate(&fixture, &bundle, &prefix, false);
    let join_rows = tsv_rows(&root.join("quant.mod_join_qc.tsv"));
    assert_eq!(
        join_rows
            .iter()
            .map(|row| row.get(11).unwrap().parse::<usize>().unwrap())
            .sum::<usize>(),
        1
    );
    for row in &join_rows {
        let assigned_reads = row.get(8).unwrap().parse::<usize>().unwrap();
        let joined_reads = row.get(7).unwrap().parse::<usize>().unwrap();
        let read_join_rate = row.get(9).unwrap().parse::<f64>().unwrap();
        if assigned_reads > 0 {
            assert_eq!(joined_reads, assigned_reads);
            assert_eq!(read_join_rate, 1.0);
        }
    }
}
