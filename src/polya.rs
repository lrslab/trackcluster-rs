//! Join Dorado `pt:i` estimates to final unique read assignments.
//!
//! Tail length is a molecule attribute, independent of alignment coordinates.
//! Unmapped primary records can supply it; secondary/supplementary records cannot
//! create additional molecules. Dorado's -1 and 0 sentinels are failed estimates.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::data::field::Tag;

use crate::model::Transcript;

pub(crate) const OUTPUT_SUFFIXES: [&str; 3] =
    [".read_polya.tsv", ".polya_qc.tsv", ".isoform_polya.tsv"];

#[derive(Clone, Debug)]
pub(crate) struct SampleInput {
    pub sample: String,
    pub group: Option<String>,
    pub bam: PathBuf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadIdMode {
    Raw,
    SamplePrefixed,
}

fn validate_identifier(kind: &str, value: &str) -> anyhow::Result<()> {
    if value.trim().is_empty() || value.chars().any(char::is_control) {
        anyhow::bail!("{kind} must not be empty or contain control characters: {value:?}");
    }
    Ok(())
}

pub(crate) fn validate_sample(sample: &str) -> anyhow::Result<()> {
    validate_identifier("poly(A) sample", sample)?;
    if sample.contains(crate::sample::SAMPLE_DELIM) {
        anyhow::bail!("poly(A) sample must not contain '::': {sample:?}");
    }
    Ok(())
}

/// Manifest paths are relative to the manifest, just like the reads manifest.
pub(crate) fn read_manifest(path: &Path) -> anyhow::Result<Vec<SampleInput>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .from_path(path)
        .with_context(|| format!("open poly(A) manifest {path:?}"))?;
    let header = reader.headers()?.clone();
    let mut columns = HashSet::new();
    for field in &header {
        if !columns.insert(field.trim().to_ascii_lowercase()) {
            anyhow::bail!("duplicate poly(A) manifest column {field:?}");
        }
    }
    let column = |name: &str| {
        header
            .iter()
            .position(|value| value.trim().eq_ignore_ascii_case(name))
    };
    let sample_column = column("sample")
        .context("poly(A) manifest header must include 'sample' and 'bam' columns")?;
    let bam_column =
        column("bam").context("poly(A) manifest header must include 'sample' and 'bam' columns")?;
    let group_column = column("group");
    let mut samples = Vec::new();
    let mut seen = HashSet::new();
    for (index, row) in reader.records().enumerate() {
        let row = row.with_context(|| format!("parse poly(A) manifest row {}", index + 2))?;
        let sample = row[sample_column].trim();
        validate_sample(sample)?;
        if !seen.insert(sample.to_owned()) {
            anyhow::bail!("duplicate poly(A) sample {sample:?}");
        }
        let bam_field = row[bam_column].trim();
        validate_identifier("poly(A) BAM path", bam_field)?;
        let bam = PathBuf::from(bam_field);
        let bam = if bam.is_absolute() {
            bam
        } else {
            path.parent().unwrap_or_else(|| Path::new(".")).join(bam)
        };
        let group = group_column
            .map(|column| row[column].trim())
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        if let Some(group) = &group {
            validate_identifier("poly(A) group", group)?;
        }
        samples.push(SampleInput {
            sample: sample.to_owned(),
            group,
            bam,
        });
    }
    if samples.is_empty() {
        anyhow::bail!("poly(A) manifest contains no samples");
    }
    Ok(samples)
}

/// Use the flow's sample labels/groups and reject a mismatched manifest early.
pub(crate) fn match_flow_samples(
    inputs: &mut [SampleInput],
    samples: &[crate::io::manifest::SampleRow],
) -> anyhow::Result<()> {
    if inputs.len() != samples.len() {
        anyhow::bail!("flow: poly(A) manifest must contain exactly the samples in --manifest");
    }
    for input in inputs {
        let sample = samples
            .iter()
            .find(|sample| sample.sample == input.sample)
            .with_context(|| format!("unknown poly(A) sample {:?} in flow", input.sample))?;
        if input.group.is_some() && input.group != sample.group {
            anyhow::bail!(
                "poly(A) group for sample {:?} conflicts with --manifest",
                input.sample
            );
        }
        input.group = sample.group.clone();
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum TailCall {
    #[default]
    MissingBam,
    MissingTag,
    AnchorNotFound,
    EstimationFailed,
    Estimated(u32),
}

impl TailCall {
    fn status(self) -> &'static str {
        match self {
            Self::MissingBam => "missing_bam_read",
            Self::MissingTag => "missing_pt_tag",
            Self::AnchorNotFound => "anchor_not_found",
            Self::EstimationFailed => "estimation_failed",
            Self::Estimated(_) => "estimated",
        }
    }

    fn raw_pt(self) -> String {
        match self {
            Self::MissingBam | Self::MissingTag => "NA".to_owned(),
            Self::AnchorNotFound => "-1".to_owned(),
            Self::EstimationFailed => "0".to_owned(),
            Self::Estimated(length) => length.to_string(),
        }
    }

    fn length(self) -> Option<u32> {
        match self {
            Self::Estimated(length) => Some(length),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct AssignedRead {
    sample: usize,
    isoform: usize,
    call: TailCall,
}

#[derive(Debug, Default)]
struct BamQc {
    records: u64,
    primary_records: u64,
    skipped_secondary: u64,
    skipped_supplementary: u64,
    unassigned_primary_records: u64,
    duplicate_assigned_primary_records: u64,
    primary_records_with_pt: u64,
    dorado_versions: String,
}

#[derive(Debug, Default)]
struct TailSummary {
    assigned: u64,
    missing_bam: u64,
    missing_tag: u64,
    anchor_not_found: u64,
    estimation_failed: u64,
    lengths: Vec<u32>,
}

impl TailSummary {
    fn add(&mut self, call: TailCall) {
        self.assigned += 1;
        match call {
            TailCall::MissingBam => self.missing_bam += 1,
            TailCall::MissingTag => self.missing_tag += 1,
            TailCall::AnchorNotFound => self.anchor_not_found += 1,
            TailCall::EstimationFailed => self.estimation_failed += 1,
            TailCall::Estimated(length) => self.lengths.push(length),
        }
    }

    fn count_fields(&self) -> Vec<String> {
        vec![
            self.assigned.to_string(),
            (self.assigned - self.missing_bam).to_string(),
            self.lengths.len().to_string(),
            self.missing_bam.to_string(),
            self.missing_tag.to_string(),
            self.anchor_not_found.to_string(),
            self.estimation_failed.to_string(),
            number(rate(self.lengths.len() as u64, self.assigned)),
        ]
    }

    /// Exact quantiles use linear interpolation at (n - 1) * p (type 7).
    fn length_fields(&mut self) -> Vec<String> {
        if self.lengths.is_empty() {
            return vec!["NA".to_owned(); 7];
        }
        self.lengths.sort_unstable();
        let mut mean = 0.0;
        let mut m2 = 0.0;
        for (index, &length) in self.lengths.iter().enumerate() {
            let delta = f64::from(length) - mean;
            mean += delta / (index + 1) as f64;
            m2 += delta * (f64::from(length) - mean);
        }
        let quantile = |p: f64| {
            let h = (self.lengths.len() - 1) as f64 * p;
            let left = h.floor() as usize;
            let right = h.ceil() as usize;
            f64::from(self.lengths[left])
                + (h - left as f64)
                    * (f64::from(self.lengths[right]) - f64::from(self.lengths[left]))
        };
        vec![
            number(Some(mean)),
            number(Some(quantile(0.5))),
            number(Some(quantile(0.25))),
            number(Some(quantile(0.75))),
            self.lengths[0].to_string(),
            self.lengths[self.lengths.len() - 1].to_string(),
            number(
                (self.lengths.len() > 1)
                    .then(|| (m2.max(0.0) / (self.lengths.len() - 1) as f64).sqrt()),
            ),
        ]
    }
}

fn rate(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| numerator as f64 / denominator as f64)
}

fn number(value: Option<f64>) -> String {
    value.map_or_else(|| "NA".to_owned(), |value| format!("{value:.6}"))
}

const COUNT_COLUMNS: [&str; 8] = [
    "assigned_reads",
    "observed_reads",
    "polya_reads",
    "missing_bam_reads",
    "missing_pt_reads",
    "anchor_not_found_reads",
    "estimation_failed_reads",
    "polya_fraction",
];
const LENGTH_COLUMNS: [&str; 7] = [
    "polya_mean_nt",
    "polya_median_nt",
    "polya_q25_nt",
    "polya_q75_nt",
    "polya_min_nt",
    "polya_max_nt",
    "polya_stddev_nt",
];

fn read_bam(
    input: &SampleInput,
    sample_index: usize,
    mode: ReadIdMode,
    assignments: &mut BTreeMap<String, AssignedRead>,
) -> anyhow::Result<BamQc> {
    let mut reader = bam::io::Reader::new(
        File::open(&input.bam).with_context(|| format!("open poly(A) BAM {:?}", input.bam))?,
    );
    let header = reader
        .read_header()
        .with_context(|| format!("read poly(A) BAM header {:?}", input.bam))?;
    use sam::header::record::value::map::program::tag;
    let versions = header
        .programs()
        .as_ref()
        .values()
        .filter(|program| {
            program
                .other_fields()
                .get(&tag::NAME)
                .is_some_and(|v| v == "dorado")
        })
        .filter_map(|program| program.other_fields().get(&tag::VERSION))
        .map(|value| String::from_utf8_lossy(value.as_ref()).into_owned())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join(",");
    let mut qc = BamQc {
        dorado_versions: versions,
        ..BamQc::default()
    };
    for result in reader.records() {
        let record = result.with_context(|| format!("decode poly(A) BAM {:?}", input.bam))?;
        qc.records += 1;
        if record.flags().is_secondary() {
            qc.skipped_secondary += 1;
            continue;
        }
        if record.flags().is_supplementary() {
            qc.skipped_supplementary += 1;
            continue;
        }
        qc.primary_records += 1;
        let name = record
            .name()
            .context("poly(A) BAM primary record has no query name")?;
        let name =
            std::str::from_utf8(name.as_ref()).context("poly(A) BAM query name is not UTF-8")?;
        validate_identifier("poly(A) BAM query name", name)?;
        if name == "*" {
            anyhow::bail!("poly(A) BAM primary record has no query name");
        }
        // Check all primary pt tags, including unassigned reads, so a BAM with
        // no estimates is distinguished from a read-ID join mismatch.
        let call = match record.data().get(&Tag::new(b'p', b't')) {
            None => TailCall::MissingTag,
            Some(value) => {
                let value = value.context("decode Dorado pt tag")?;
                let length = value
                    .as_int()
                    .context("Dorado pt tag must have integer type (pt:i)")?;
                qc.primary_records_with_pt += 1;
                match length {
                    -1 => TailCall::AnchorNotFound,
                    0 => TailCall::EstimationFailed,
                    1.. => TailCall::Estimated(u32::try_from(length).context("Dorado pt exceeds u32")?),
                    _ => anyhow::bail!("invalid Dorado pt value {length} for read {name:?}; expected -1, 0, or positive tail length"),
                }
            }
        };
        let read_id = match mode {
            ReadIdMode::Raw => name.to_owned(),
            ReadIdMode::SamplePrefixed => crate::sample::tagged_read_name(&input.sample, name),
        };
        let Some(assigned) = assignments.get_mut(&read_id) else {
            qc.unassigned_primary_records += 1;
            continue;
        };
        if assigned.sample != sample_index {
            anyhow::bail!("poly(A) read {read_id:?} joined to the wrong sample");
        }
        if assigned.call != TailCall::MissingBam {
            if assigned.call != call {
                anyhow::bail!(
                    "conflicting Dorado pt values in primary records for read {read_id:?}"
                );
            }
            qc.duplicate_assigned_primary_records += 1;
        } else {
            assigned.call = call;
        }
    }
    if qc.primary_records > 0 && qc.primary_records_with_pt == 0 {
        anyhow::bail!("poly(A) BAM {:?} has no Dorado pt:i tags on primary reads; run Dorado with --estimate-poly-a and preserve tags during alignment", input.bam);
    }
    Ok(qc)
}

pub(crate) struct AggregateResult {
    samples: Vec<SampleInput>,
    isoforms: Vec<Transcript>,
    assignments: BTreeMap<String, AssignedRead>,
    qc: Vec<BamQc>,
}

pub(crate) fn aggregate(
    inputs: &[SampleInput],
    isoforms: Vec<Transcript>,
    pairs: Vec<(String, String)>,
    mode: ReadIdMode,
) -> anyhow::Result<AggregateResult> {
    if inputs.is_empty() || (mode == ReadIdMode::Raw && inputs.len() != 1) {
        anyhow::bail!("poly(A) requires samples; raw read IDs require exactly one sample");
    }
    if isoforms.is_empty() {
        anyhow::bail!("poly(A) isoform catalog is empty");
    }
    crate::identity::validate_isoform_ids(&isoforms)?;
    let mut sample_indices = HashMap::new();
    for (index, input) in inputs.iter().enumerate() {
        validate_sample(&input.sample)?;
        if let Some(group) = &input.group {
            validate_identifier("poly(A) group", group)?;
        }
        if sample_indices
            .insert(input.sample.as_str(), index)
            .is_some()
        {
            anyhow::bail!("duplicate poly(A) sample {:?}", input.sample);
        }
    }
    let isoform_indices = isoforms
        .iter()
        .enumerate()
        .map(|(index, isoform)| (isoform.name.as_str(), index))
        .collect::<HashMap<_, _>>();
    let mut assignments: BTreeMap<String, AssignedRead> = BTreeMap::new();
    for (read_id, isoform_id) in pairs {
        validate_identifier("assigned read ID", &read_id)?;
        let isoform = *isoform_indices
            .get(isoform_id.as_str())
            .with_context(|| format!("poly(A) mapping contains unknown isoform {isoform_id:?}"))?;
        let sample = match mode {
            ReadIdMode::Raw => 0,
            ReadIdMode::SamplePrefixed => {
                let (sample, _) =
                    crate::sample::split_tagged_read_name(&read_id).with_context(|| {
                        format!("poly(A) manifest mode requires sample::read IDs; got {read_id:?}")
                    })?;
                *sample_indices.get(sample).with_context(|| {
                    format!("poly(A) mapping contains unknown sample {sample:?}")
                })?
            }
        };
        if let Some(existing) = assignments.get(&read_id) {
            if existing.isoform != isoform {
                anyhow::bail!("poly(A) requires a unique read-to-isoform mapping; read {read_id:?} is assigned to multiple isoforms");
            }
            continue;
        }
        assignments.insert(
            read_id,
            AssignedRead {
                sample,
                isoform,
                call: TailCall::MissingBam,
            },
        );
    }
    let mut qc = Vec::new();
    for (index, input) in inputs.iter().enumerate() {
        qc.push(
            read_bam(input, index, mode, &mut assignments)
                .with_context(|| format!("import Dorado poly(A) for sample {:?}", input.sample))?,
        );
    }
    Ok(AggregateResult {
        samples: inputs.to_vec(),
        isoforms,
        assignments,
        qc,
    })
}

pub(crate) fn output_paths(prefix: &Path) -> [PathBuf; 3] {
    OUTPUT_SUFFIXES.map(|suffix| {
        let mut path = prefix.as_os_str().to_os_string();
        path.push(suffix);
        PathBuf::from(path)
    })
}

impl AggregateResult {
    fn write_reads(&self, writer: impl Write) -> anyhow::Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        writer.write_record([
            "sample",
            "group",
            "gene",
            "isoform_id",
            "read_id",
            "dorado_pt",
            "polya_length_nt",
            "status",
        ])?;
        for (read_id, assigned) in &self.assignments {
            let input = &self.samples[assigned.sample];
            let isoform = &self.isoforms[assigned.isoform];
            writer.write_record([
                input.sample.as_str(),
                input.group.as_deref().unwrap_or("NA"),
                crate::identity::gene_id(isoform),
                isoform.name.as_str(),
                read_id.as_str(),
                &assigned.call.raw_pt(),
                &assigned
                    .call
                    .length()
                    .map_or_else(|| "NA".to_owned(), |length| length.to_string()),
                assigned.call.status(),
            ])?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_qc(&self, writer: impl Write) -> anyhow::Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        let mut header = vec![
            "sample",
            "group",
            "bam",
            "dorado_versions",
            "bam_records",
            "primary_records",
            "skipped_secondary",
            "skipped_supplementary",
            "unassigned_primary_records",
            "duplicate_assigned_primary_records",
            "primary_records_with_pt",
        ];
        header.extend(COUNT_COLUMNS);
        header.push("read_join_rate");
        writer.write_record(header)?;
        let mut summaries = (0..self.samples.len())
            .map(|_| TailSummary::default())
            .collect::<Vec<_>>();
        for assigned in self.assignments.values() {
            summaries[assigned.sample].add(assigned.call);
        }
        for (index, input) in self.samples.iter().enumerate() {
            let qc = &self.qc[index];
            let summary = &summaries[index];
            let mut row = vec![
                input.sample.clone(),
                input.group.clone().unwrap_or_else(|| "NA".to_owned()),
                input.bam.to_string_lossy().into_owned(),
                if qc.dorado_versions.is_empty() {
                    "NA".to_owned()
                } else {
                    qc.dorado_versions.clone()
                },
                qc.records.to_string(),
                qc.primary_records.to_string(),
                qc.skipped_secondary.to_string(),
                qc.skipped_supplementary.to_string(),
                qc.unassigned_primary_records.to_string(),
                qc.duplicate_assigned_primary_records.to_string(),
                qc.primary_records_with_pt.to_string(),
            ];
            row.extend(summary.count_fields());
            row.push(number(rate(
                summary.assigned - summary.missing_bam,
                summary.assigned,
            )));
            writer.write_record(row)?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_isoforms(&self, writer: impl Write) -> anyhow::Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        let mut header = vec!["sample", "group", "gene", "isoform_id"];
        header.extend(COUNT_COLUMNS);
        header.extend(LENGTH_COLUMNS);
        writer.write_record(header)?;
        let mut summaries: HashMap<(usize, usize), TailSummary> = HashMap::new();
        for assigned in self.assignments.values() {
            summaries
                .entry((assigned.sample, assigned.isoform))
                .or_default()
                .add(assigned.call);
        }
        let mut indices = (0..self.isoforms.len()).collect::<Vec<_>>();
        indices.sort_by(|&left, &right| {
            (
                crate::identity::gene_id(&self.isoforms[left]),
                &self.isoforms[left].name,
            )
                .cmp(&(
                    crate::identity::gene_id(&self.isoforms[right]),
                    &self.isoforms[right].name,
                ))
        });
        for (sample, input) in self.samples.iter().enumerate() {
            for &index in &indices {
                let isoform = &self.isoforms[index];
                let mut summary = summaries.remove(&(sample, index)).unwrap_or_default();
                let mut row = vec![
                    input.sample.clone(),
                    input.group.clone().unwrap_or_else(|| "NA".to_owned()),
                    crate::identity::gene_id(isoform).to_owned(),
                    isoform.name.clone(),
                ];
                row.extend(summary.count_fields());
                row.extend(summary.length_fields());
                writer.write_record(row)?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    pub(crate) fn write(&self, prefix: &Path) -> anyhow::Result<()> {
        let [reads, qc, isoforms] = output_paths(prefix);
        // Parse and validate the complete input before publishing any table.
        // The isoform summary is published last.
        crate::flow::artifact_manifest::atomic_write_with(&reads, |writer| {
            self.write_reads(writer)
        })?;
        crate::flow::artifact_manifest::atomic_write_with(&qc, |writer| self.write_qc(writer))?;
        crate::flow::artifact_manifest::atomic_write_with(&isoforms, |writer| {
            self.write_isoforms(writer)
        })?;
        eprintln!(
            "polya-aggregate: samples={} isoforms={} assigned_reads={} summary={isoforms:?}",
            self.samples.len(),
            self.isoforms.len(),
            self.assignments.len()
        );
        Ok(())
    }
}
