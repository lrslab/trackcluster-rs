//! Explicit, effect-only isoform and condition modification contrasts.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;

use anyhow::Context;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};

/// Supported descriptive contrast families.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ContrastType {
    /// Within-sample difference between two isoforms at a shared site.
    IsoformEffect,
    /// Difference between sample-group means for one isoform/site.
    ConditionEffect,
    /// Difference of within-sample isoform differences between two groups.
    IsoformConditionInteraction,
}

impl ContrastType {
    /// Return the stable TSV token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::IsoformEffect => "isoform_effect",
            Self::ConditionEffect => "condition_effect",
            Self::IsoformConditionInteraction => "isoform_condition_interaction",
        }
    }
}

impl std::str::FromStr for ContrastType {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "isoform_effect" => Ok(Self::IsoformEffect),
            "condition_effect" => Ok(Self::ConditionEffect),
            "isoform_condition_interaction" => Ok(Self::IsoformConditionInteraction),
            _ => Err(format!(
                "invalid contrast_type {value:?}; expected isoform_effect, condition_effect, or isoform_condition_interaction"
            )),
        }
    }
}

/// One explicit contrast request.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContrastSpec {
    /// Contrast family.
    pub contrast_type: ContrastType,
    /// Assay identifier.
    pub assay_id: String,
    /// Gene identifier.
    pub gene: String,
    /// Genomic site identifier.
    pub site_id: String,
    /// SAM-style modification code.
    pub mod_code: String,
    /// Primary isoform.
    pub isoform_a: String,
    /// Comparator isoform for isoform/interaction contrasts.
    pub isoform_b: Option<String>,
    /// Baseline or optional filter group.
    pub group_a: Option<String>,
    /// Comparator group for condition/interaction contrasts.
    pub group_b: Option<String>,
}

/// One effect-only contrast result.
#[derive(Clone, Debug, PartialEq)]
pub struct ContrastRow {
    /// Contrast family.
    pub contrast_type: ContrastType,
    /// Assay identifier.
    pub assay_id: String,
    /// Applied hard-call threshold.
    pub analysis_threshold: f64,
    /// Gene identifier.
    pub gene: String,
    /// Genomic site identifier.
    pub site_id: String,
    /// SAM-style modification code.
    pub mod_code: String,
    /// Primary isoform.
    pub isoform_a: String,
    /// Comparator isoform.
    pub isoform_b: Option<String>,
    /// Baseline or filter group.
    pub group_a: Option<String>,
    /// Comparator group.
    pub group_b: Option<String>,
    /// Samples or paired samples contributing to the effect.
    pub n_eligible_samples: usize,
    /// Mean fraction difference for isoform or condition effects.
    pub delta_fraction: Option<f64>,
    /// Descriptive pooled odds ratio with a 0.5 continuity correction.
    pub odds_ratio: Option<f64>,
    /// Difference-in-differences for interaction contrasts.
    pub interaction_delta: Option<f64>,
    /// Inferential p-value; V1 always emits NA.
    pub p_value: Option<f64>,
    /// Multiple-testing adjusted p-value; V1 always emits NA.
    pub q_value: Option<f64>,
    /// Stable method identifier.
    pub method: String,
    /// Stable eligibility reason.
    pub eligibility_reason: String,
}

#[derive(Clone, Debug)]
struct DesignRow {
    assay_id: String,
    threshold: f64,
    sample: String,
    group: Option<String>,
    gene: String,
    site_id: String,
    mod_code: String,
    isoform_id: String,
    modified: u64,
    unmodified: u64,
    fraction: Option<f64>,
    eligible: bool,
}

const DESIGN_COLUMNS: [&str; 13] = [
    "assay_id",
    "analysis_threshold",
    "sample",
    "group",
    "gene",
    "site_id",
    "mod_code",
    "isoform_id",
    "n_modified",
    "n_unmodified",
    "mod_fraction",
    "eligibility",
    "eligibility_reason",
];

const SPEC_COLUMNS: [&str; 9] = [
    "contrast_type",
    "assay_id",
    "gene",
    "site_id",
    "mod_code",
    "isoform_a",
    "isoform_b",
    "group_a",
    "group_b",
];

fn parse_optional(value: &str, field: &str) -> anyhow::Result<Option<String>> {
    match value {
        "NA" => Ok(None),
        "" => anyhow::bail!("{field} must use NA for a missing value"),
        value if value.chars().any(char::is_control) => {
            anyhow::bail!("{field} contains a control character")
        }
        value => Ok(Some(value.to_owned())),
    }
}

fn parse_required(record: &StringRecord, index: usize, field: &str) -> anyhow::Result<String> {
    let value = record.get(index).unwrap_or_default();
    if value.is_empty() || value == "NA" || value.chars().any(char::is_control) {
        anyhow::bail!("{field} must not be empty, NA, or contain control characters");
    }
    Ok(value.to_owned())
}

fn validate_contrast_spec(spec: &ContrastSpec) -> Result<(), String> {
    match spec.contrast_type {
        ContrastType::IsoformEffect if spec.isoform_b.is_none() || spec.group_b.is_some() => {
            return Err("isoform_effect requires isoform_b and group_b=NA".to_owned());
        }
        ContrastType::ConditionEffect
            if spec.isoform_b.is_some() || spec.group_a.is_none() || spec.group_b.is_none() =>
        {
            return Err("condition_effect requires isoform_b=NA and both groups".to_owned());
        }
        ContrastType::IsoformConditionInteraction
            if spec.isoform_b.is_none() || spec.group_a.is_none() || spec.group_b.is_none() =>
        {
            return Err(
                "isoform_condition_interaction requires isoform_b and both groups".to_owned(),
            );
        }
        _ => {}
    }

    if matches!(
        spec.contrast_type,
        ContrastType::IsoformEffect | ContrastType::IsoformConditionInteraction
    ) && spec.isoform_b.as_deref() == Some(spec.isoform_a.as_str())
    {
        return Err(format!(
            "{} requires isoform_a and isoform_b to differ",
            spec.contrast_type.as_str()
        ));
    }
    if matches!(
        spec.contrast_type,
        ContrastType::ConditionEffect | ContrastType::IsoformConditionInteraction
    ) && spec.group_a == spec.group_b
    {
        return Err(format!(
            "{} requires group_a and group_b to differ",
            spec.contrast_type.as_str()
        ));
    }
    Ok(())
}

fn read_design(path: &Path) -> anyhow::Result<Vec<DesignRow>> {
    let file = std::fs::File::open(path).with_context(|| format!("open design {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(file);
    let expected = StringRecord::from(DESIGN_COLUMNS.to_vec());
    if reader.headers()? != &expected {
        anyhow::bail!("design {path:?} header mismatch; expected {DESIGN_COLUMNS:?}");
    }
    let mut rows = Vec::new();
    let mut seen = BTreeMap::new();
    for (index, result) in reader.records().enumerate() {
        let line = index + 2;
        let record = result.with_context(|| format!("parse design {path:?}:{line}"))?;
        let threshold = parse_required(&record, 1, "analysis_threshold")?
            .parse::<f64>()
            .with_context(|| format!("invalid analysis_threshold at {path:?}:{line}"))?;
        if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
            anyhow::bail!("analysis_threshold at {path:?}:{line} must be in [0, 1]");
        }
        let modified = parse_required(&record, 8, "n_modified")?
            .parse::<u64>()
            .with_context(|| format!("invalid n_modified at {path:?}:{line}"))?;
        let unmodified = parse_required(&record, 9, "n_unmodified")?
            .parse::<u64>()
            .with_context(|| format!("invalid n_unmodified at {path:?}:{line}"))?;
        let fraction = match record.get(10).unwrap_or_default() {
            "NA" => None,
            value => {
                let value = value
                    .parse::<f64>()
                    .with_context(|| format!("invalid mod_fraction at {path:?}:{line}"))?;
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    anyhow::bail!("mod_fraction at {path:?}:{line} must be in [0, 1] or NA");
                }
                Some(value)
            }
        };
        if let Some(fraction) = fraction {
            let callable = modified + unmodified;
            if callable == 0 || (fraction - modified as f64 / callable as f64).abs() > 1e-12 {
                anyhow::bail!(
                    "mod_fraction at {path:?}:{line} is inconsistent with integer counts"
                );
            }
        }
        let eligibility = parse_required(&record, 11, "eligibility")?;
        let reason = parse_required(&record, 12, "eligibility_reason")?;
        let eligible = match eligibility.as_str() {
            "eligible" if reason == "ok" && fraction.is_some() => true,
            "ineligible" if reason != "ok" => false,
            _ => {
                anyhow::bail!("inconsistent eligibility/reason/fraction at design {path:?}:{line}")
            }
        };
        let row = DesignRow {
            assay_id: parse_required(&record, 0, "assay_id")?,
            threshold,
            sample: parse_required(&record, 2, "sample")?,
            group: parse_optional(record.get(3).unwrap_or_default(), "group")?,
            gene: parse_required(&record, 4, "gene")?,
            site_id: parse_required(&record, 5, "site_id")?,
            mod_code: parse_required(&record, 6, "mod_code")?,
            isoform_id: parse_required(&record, 7, "isoform_id")?,
            modified,
            unmodified,
            fraction,
            eligible,
        };
        let key = (
            row.assay_id.clone(),
            row.sample.clone(),
            row.gene.clone(),
            row.site_id.clone(),
            row.mod_code.clone(),
            row.isoform_id.clone(),
        );
        if seen.insert(key, line).is_some() {
            anyhow::bail!("duplicate design key at {path:?}:{line}");
        }
        rows.push(row);
    }
    if rows.is_empty() {
        anyhow::bail!("design {path:?} has no rows");
    }
    Ok(rows)
}

/// Read and validate an explicit contrast specification TSV.
pub fn read_contrast_specs(path: &Path) -> anyhow::Result<Vec<ContrastSpec>> {
    let file = std::fs::File::open(path).with_context(|| format!("open contrasts {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(file);
    let expected = StringRecord::from(SPEC_COLUMNS.to_vec());
    if reader.headers()? != &expected {
        anyhow::bail!("contrasts {path:?} header mismatch; expected {SPEC_COLUMNS:?}");
    }
    let mut specs = Vec::new();
    for (index, result) in reader.records().enumerate() {
        let line = index + 2;
        let record = result.with_context(|| format!("parse contrasts {path:?}:{line}"))?;
        let contrast_type = parse_required(&record, 0, "contrast_type")?
            .parse::<ContrastType>()
            .map_err(anyhow::Error::msg)?;
        let spec = ContrastSpec {
            contrast_type,
            assay_id: parse_required(&record, 1, "assay_id")?,
            gene: parse_required(&record, 2, "gene")?,
            site_id: parse_required(&record, 3, "site_id")?,
            mod_code: parse_required(&record, 4, "mod_code")?,
            isoform_a: parse_required(&record, 5, "isoform_a")?,
            isoform_b: parse_optional(record.get(6).unwrap_or_default(), "isoform_b")?,
            group_a: parse_optional(record.get(7).unwrap_or_default(), "group_a")?,
            group_b: parse_optional(record.get(8).unwrap_or_default(), "group_b")?,
        };
        validate_contrast_spec(&spec)
            .map_err(|error| anyhow::anyhow!("{error} at {path:?}:{line}"))?;
        specs.push(spec);
    }
    if specs.is_empty() {
        anyhow::bail!("contrasts {path:?} has no rows");
    }
    Ok(specs)
}

fn pooled_odds_ratio(left: &[&DesignRow], right: &[&DesignRow]) -> Option<f64> {
    if left.is_empty() || right.is_empty() {
        return None;
    }
    let counts = |rows: &[&DesignRow]| {
        rows.iter()
            .fold((0u64, 0u64), |(modified, unmodified), row| {
                (modified + row.modified, unmodified + row.unmodified)
            })
    };
    let (left_modified, left_unmodified) = counts(left);
    let (right_modified, right_unmodified) = counts(right);
    Some(
        ((left_modified as f64 + 0.5) * (right_unmodified as f64 + 0.5))
            / ((left_unmodified as f64 + 0.5) * (right_modified as f64 + 0.5)),
    )
}

fn mean(values: &[f64]) -> Option<f64> {
    (!values.is_empty()).then(|| values.iter().sum::<f64>() / values.len() as f64)
}

/// Calculate explicit descriptive contrasts without inferential p-values.
pub fn calculate_contrasts(
    design_path: &Path,
    specs: &[ContrastSpec],
) -> anyhow::Result<Vec<ContrastRow>> {
    if specs.is_empty() {
        anyhow::bail!("at least one contrast specification is required");
    }
    for (index, spec) in specs.iter().enumerate() {
        validate_contrast_spec(spec).map_err(|error| {
            anyhow::anyhow!("invalid contrast specification {}: {error}", index + 1)
        })?;
    }
    let design = read_design(design_path)?;
    let mut thresholds = HashMap::new();
    let mut eligible_by_site: HashMap<(&str, &str, &str, &str), Vec<&DesignRow>> = HashMap::new();
    for row in &design {
        match thresholds.insert(row.assay_id.as_str(), row.threshold) {
            Some(existing) if existing != row.threshold => anyhow::bail!(
                "assay {:?} has inconsistent analysis thresholds in design",
                row.assay_id
            ),
            _ => {}
        }
        if row.eligible {
            eligible_by_site
                .entry((
                    row.assay_id.as_str(),
                    row.gene.as_str(),
                    row.site_id.as_str(),
                    row.mod_code.as_str(),
                ))
                .or_default()
                .push(row);
        }
    }

    let mut results = Vec::with_capacity(specs.len());
    for spec in specs {
        let threshold = *thresholds
            .get(spec.assay_id.as_str())
            .with_context(|| format!("contrast references unknown assay {:?}", spec.assay_id))?;
        let matching = eligible_by_site
            .get(&(
                spec.assay_id.as_str(),
                spec.gene.as_str(),
                spec.site_id.as_str(),
                spec.mod_code.as_str(),
            ))
            .map(Vec::as_slice)
            .unwrap_or_default();

        let mut delta_fraction = None;
        let mut odds_ratio = None;
        let mut interaction_delta = None;
        let n_eligible_samples;
        let eligibility_reason;

        match spec.contrast_type {
            ContrastType::IsoformEffect => {
                let isoform_b = spec.isoform_b.as_deref().expect("validated spec");
                let by_key = matching
                    .iter()
                    .map(|row| ((row.sample.as_str(), row.isoform_id.as_str()), *row))
                    .collect::<HashMap<_, _>>();
                let mut deltas = Vec::new();
                let mut rows_a = Vec::new();
                let mut rows_b = Vec::new();
                for row_a in matching.iter().copied().filter(|row| {
                    row.isoform_id == spec.isoform_a
                        && spec
                            .group_a
                            .as_deref()
                            .is_none_or(|group| row.group.as_deref() == Some(group))
                }) {
                    if let Some(row_b) = by_key.get(&(row_a.sample.as_str(), isoform_b)).copied() {
                        if row_b.group == row_a.group {
                            deltas.push(row_a.fraction.unwrap() - row_b.fraction.unwrap());
                            rows_a.push(row_a);
                            rows_b.push(row_b);
                        }
                    }
                }
                n_eligible_samples = deltas.len();
                delta_fraction = mean(&deltas);
                odds_ratio = pooled_odds_ratio(&rows_a, &rows_b);
                eligibility_reason = if deltas.is_empty() {
                    "no_shared_eligible_samples"
                } else {
                    "ok"
                };
            }
            ContrastType::ConditionEffect => {
                let group_a = spec.group_a.as_deref().expect("validated spec");
                let group_b = spec.group_b.as_deref().expect("validated spec");
                let rows_a = matching
                    .iter()
                    .copied()
                    .filter(|row| {
                        row.isoform_id == spec.isoform_a && row.group.as_deref() == Some(group_a)
                    })
                    .collect::<Vec<_>>();
                let rows_b = matching
                    .iter()
                    .copied()
                    .filter(|row| {
                        row.isoform_id == spec.isoform_a && row.group.as_deref() == Some(group_b)
                    })
                    .collect::<Vec<_>>();
                n_eligible_samples = rows_a.len() + rows_b.len();
                let mean_a = mean(
                    &rows_a
                        .iter()
                        .map(|row| row.fraction.unwrap())
                        .collect::<Vec<_>>(),
                );
                let mean_b = mean(
                    &rows_b
                        .iter()
                        .map(|row| row.fraction.unwrap())
                        .collect::<Vec<_>>(),
                );
                if let (Some(mean_a), Some(mean_b)) = (mean_a, mean_b) {
                    delta_fraction = Some(mean_b - mean_a);
                    odds_ratio = pooled_odds_ratio(&rows_b, &rows_a);
                    eligibility_reason = "ok";
                } else {
                    eligibility_reason = "missing_eligible_group";
                }
            }
            ContrastType::IsoformConditionInteraction => {
                let isoform_b = spec.isoform_b.as_deref().expect("validated spec");
                let group_a = spec.group_a.as_deref().expect("validated spec");
                let group_b = spec.group_b.as_deref().expect("validated spec");
                let by_key = matching
                    .iter()
                    .map(|row| ((row.sample.as_str(), row.isoform_id.as_str()), *row))
                    .collect::<HashMap<_, _>>();
                let mut deltas_a = Vec::new();
                let mut deltas_b = Vec::new();
                for row_a in matching
                    .iter()
                    .copied()
                    .filter(|row| row.isoform_id == spec.isoform_a)
                {
                    let Some(row_b) = by_key.get(&(row_a.sample.as_str(), isoform_b)).copied()
                    else {
                        continue;
                    };
                    let delta = row_a.fraction.unwrap() - row_b.fraction.unwrap();
                    match row_a.group.as_deref() {
                        Some(group) if group == group_a && row_b.group == row_a.group => {
                            deltas_a.push(delta)
                        }
                        Some(group) if group == group_b && row_b.group == row_a.group => {
                            deltas_b.push(delta)
                        }
                        _ => {}
                    }
                }
                n_eligible_samples = deltas_a.len() + deltas_b.len();
                if let (Some(mean_a), Some(mean_b)) = (mean(&deltas_a), mean(&deltas_b)) {
                    interaction_delta = Some(mean_b - mean_a);
                    eligibility_reason = "ok";
                } else {
                    eligibility_reason = "missing_paired_group";
                }
            }
        }

        results.push(ContrastRow {
            contrast_type: spec.contrast_type,
            assay_id: spec.assay_id.clone(),
            analysis_threshold: threshold,
            gene: spec.gene.clone(),
            site_id: spec.site_id.clone(),
            mod_code: spec.mod_code.clone(),
            isoform_a: spec.isoform_a.clone(),
            isoform_b: spec.isoform_b.clone(),
            group_a: spec.group_a.clone(),
            group_b: spec.group_b.clone(),
            n_eligible_samples,
            delta_fraction,
            odds_ratio,
            interaction_delta,
            p_value: None,
            q_value: None,
            method: "effect_only".to_owned(),
            eligibility_reason: eligibility_reason.to_owned(),
        });
    }
    results.sort_by(|left, right| {
        left.contrast_type
            .cmp(&right.contrast_type)
            .then_with(|| left.assay_id.cmp(&right.assay_id))
            .then_with(|| left.gene.cmp(&right.gene))
            .then_with(|| left.site_id.cmp(&right.site_id))
            .then_with(|| left.mod_code.cmp(&right.mod_code))
            .then_with(|| left.isoform_a.cmp(&right.isoform_a))
            .then_with(|| left.isoform_b.cmp(&right.isoform_b))
            .then_with(|| left.group_a.cmp(&right.group_a))
            .then_with(|| left.group_b.cmp(&right.group_b))
    });
    Ok(results)
}

fn optional<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

/// Write deterministic effect-only contrast rows.
pub fn write_contrasts_tsv<W: Write>(writer: W, rows: &[ContrastRow]) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record([
        "contrast_type",
        "assay_id",
        "analysis_threshold",
        "gene",
        "site_id",
        "mod_code",
        "isoform_a",
        "isoform_b",
        "group_a",
        "group_b",
        "n_eligible_samples",
        "delta_fraction",
        "odds_ratio",
        "interaction_delta",
        "p_value",
        "q_value",
        "method",
        "eligibility_reason",
    ])?;
    for row in rows {
        output.write_record([
            row.contrast_type.as_str().to_owned(),
            row.assay_id.clone(),
            row.analysis_threshold.to_string(),
            row.gene.clone(),
            row.site_id.clone(),
            row.mod_code.clone(),
            row.isoform_a.clone(),
            optional(row.isoform_b.clone()),
            optional(row.group_a.clone()),
            optional(row.group_b.clone()),
            row.n_eligible_samples.to_string(),
            optional(row.delta_fraction),
            optional(row.odds_ratio),
            optional(row.interaction_delta),
            optional(row.p_value),
            optional(row.q_value),
            row.method.clone(),
            row.eligibility_reason.clone(),
        ])?;
    }
    output.flush().context("flush modification contrasts")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn temp_file(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "trackcluster-contrast-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    fn design_row(
        sample: &str,
        group: &str,
        isoform: &str,
        modified: u64,
        unmodified: u64,
    ) -> String {
        let fraction = modified as f64 / (modified + unmodified) as f64;
        format!("a1\t0.5\t{sample}\t{group}\tG1\tchr1:10:+\tA+a\t{isoform}\t{modified}\t{unmodified}\t{fraction}\teligible\tok\n")
    }

    #[test]
    fn calculates_all_three_effect_only_contrasts_without_p_values() {
        let design = temp_file("design.tsv");
        let mut text = format!("{}\n", DESIGN_COLUMNS.join("\t"));
        text.push_str(&design_row("S1", "control", "i1", 8, 2));
        text.push_str(&design_row("S1", "control", "i2", 2, 8));
        text.push_str(&design_row("S2", "treated", "i1", 4, 6));
        text.push_str(&design_row("S2", "treated", "i2", 3, 7));
        fs::write(&design, text).unwrap();
        let common = |contrast_type| ContrastSpec {
            contrast_type,
            assay_id: "a1".to_owned(),
            gene: "G1".to_owned(),
            site_id: "chr1:10:+".to_owned(),
            mod_code: "A+a".to_owned(),
            isoform_a: "i1".to_owned(),
            isoform_b: Some("i2".to_owned()),
            group_a: Some("control".to_owned()),
            group_b: Some("treated".to_owned()),
        };
        let specs = vec![
            ContrastSpec {
                group_b: None,
                ..common(ContrastType::IsoformEffect)
            },
            ContrastSpec {
                isoform_b: None,
                ..common(ContrastType::ConditionEffect)
            },
            common(ContrastType::IsoformConditionInteraction),
        ];
        let rows = calculate_contrasts(&design, &specs).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].delta_fraction, Some(0.6000000000000001));
        assert_eq!(rows[1].delta_fraction, Some(-0.4));
        assert!((rows[2].interaction_delta.unwrap() + 0.5).abs() < 1e-12);
        assert!(rows.iter().all(|row| row.p_value.is_none()
            && row.q_value.is_none()
            && row.method == "effect_only"));
        let _ = fs::remove_file(design);
    }

    #[test]
    fn rejects_isoform_and_group_self_comparisons_from_tsv() {
        for (label, row, expected) in [
            (
                "same-isoform",
                "isoform_effect\ta1\tG1\tchr1:10:+\tA+a\ti1\ti1\tNA\tNA\n",
                "isoform_a and isoform_b to differ",
            ),
            (
                "same-interaction-isoform",
                "isoform_condition_interaction\ta1\tG1\tchr1:10:+\tA+a\ti1\ti1\tcontrol\ttreated\n",
                "isoform_a and isoform_b to differ",
            ),
            (
                "same-condition-group",
                "condition_effect\ta1\tG1\tchr1:10:+\tA+a\ti1\tNA\tcontrol\tcontrol\n",
                "group_a and group_b to differ",
            ),
            (
                "same-interaction-group",
                "isoform_condition_interaction\ta1\tG1\tchr1:10:+\tA+a\ti1\ti2\tcontrol\tcontrol\n",
                "group_a and group_b to differ",
            ),
        ] {
            let path = temp_file(label);
            fs::write(&path, format!("{}\n{row}", SPEC_COLUMNS.join("\t"))).unwrap();
            let error = read_contrast_specs(&path).unwrap_err();
            assert!(format!("{error:#}").contains(expected), "{error:#}");
            let _ = fs::remove_file(path);
        }
    }

    #[test]
    fn calculate_contrasts_revalidates_programmatic_specs_before_reading_design() {
        let missing_design = temp_file("missing-design.tsv");
        let spec = ContrastSpec {
            contrast_type: ContrastType::IsoformEffect,
            assay_id: "a1".to_owned(),
            gene: "G1".to_owned(),
            site_id: "chr1:10:+".to_owned(),
            mod_code: "A+a".to_owned(),
            isoform_a: "i1".to_owned(),
            isoform_b: Some("i1".to_owned()),
            group_a: None,
            group_b: None,
        };
        let error = calculate_contrasts(&missing_design, &[spec]).unwrap_err();
        assert!(
            format!("{error:#}").contains("isoform_a and isoform_b to differ"),
            "{error:#}"
        );
    }
}
