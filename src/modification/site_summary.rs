//! Streaming reduction of isoform/site audit TSVs into per-site QC summaries.

use std::collections::{btree_map::Entry, BTreeMap};
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;
use csv::{Reader, ReaderBuilder, StringRecord, WriterBuilder};

use crate::model::Strand;

const SITE_COLUMNS: [&str; 27] = [
    "assay_id",
    "analysis_threshold",
    "sample",
    "group",
    "gene",
    "isoform_id",
    "site_id",
    "chrom",
    "pos0",
    "strand",
    "mod_code",
    "context",
    "site_state",
    "coverage_basis",
    "n_assigned",
    "n_covering",
    "n_candidate",
    "n_callable",
    "n_modified",
    "n_unmodified",
    "n_unknown",
    "mod_fraction",
    "mean_probability",
    "ci_low",
    "ci_high",
    "eligibility",
    "eligibility_reason",
];

pub(crate) const SITE_SUMMARY_COLUMNS: [&str; 27] = [
    "assay_id",
    "analysis_threshold",
    "sample",
    "group",
    "gene",
    "site_id",
    "chrom",
    "pos0",
    "strand",
    "mod_code",
    "context",
    "coverage_basis",
    "n_isoforms_total",
    "n_isoforms_assigned",
    "n_isoforms_present",
    "n_isoforms_eligible",
    "n_isoforms_site_absent",
    "n_isoforms_context_dependent",
    "n_isoforms_reference_base_mismatch",
    "n_isoforms_unprojectable",
    "n_isoforms_incomplete_candidate_universe",
    "n_isoforms_join_rate_low",
    "n_isoforms_low_callable",
    "n_isoforms_other_ineligible",
    "min_eligible_n_covering",
    "min_eligible_n_callable",
    "summary_state",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ParsedSiteState {
    Present,
    StructurallyAbsent,
    ContextDependent,
    ReferenceBaseMismatch,
    Unprojectable,
}

impl ParsedSiteState {
    fn parse(value: &str) -> anyhow::Result<Self> {
        match value {
            "present" => Ok(Self::Present),
            "structurally_absent" => Ok(Self::StructurallyAbsent),
            "context_dependent" => Ok(Self::ContextDependent),
            "reference_base_mismatch" => Ok(Self::ReferenceBaseMismatch),
            "unprojectable" => Ok(Self::Unprojectable),
            _ => anyhow::bail!("invalid site_state {value:?}"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ParsedCoverageBasis {
    Unavailable,
    BedApproximate,
    BamExact,
}

impl ParsedCoverageBasis {
    fn parse(value: &str) -> anyhow::Result<Self> {
        match value {
            "unavailable" => Ok(Self::Unavailable),
            "bed_approximate" => Ok(Self::BedApproximate),
            "bam_exact" => Ok(Self::BamExact),
            _ => anyhow::bail!("invalid coverage_basis {value:?}"),
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::Unavailable => "unavailable",
            Self::BedApproximate => "bed_approximate",
            Self::BamExact => "bam_exact",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ParsedEligibilityReason(String);

impl ParsedEligibilityReason {
    fn parse(value: &str) -> anyhow::Result<Self> {
        if !value
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_')
        {
            anyhow::bail!(
                "eligibility_reason must be a lowercase ASCII token containing letters, digits, or underscores"
            );
        }
        Ok(Self(value.to_owned()))
    }

    fn as_str(&self) -> &str {
        &self.0
    }

    fn is_eligible(&self) -> bool {
        self.as_str() == "ok"
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SummaryState {
    SharedEligible,
    SingleEligible,
    NoEligibleIsoform,
}

impl SummaryState {
    const fn from_count(count: u64) -> Self {
        match count {
            0 => Self::NoEligibleIsoform,
            1 => Self::SingleEligible,
            _ => Self::SharedEligible,
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::SharedEligible => "shared_eligible",
            Self::SingleEligible => "single_eligible",
            Self::NoEligibleIsoform => "no_eligible_isoform",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct FullRowKey {
    assay_id: String,
    sample: String,
    gene: String,
    isoform_id: String,
    chrom: String,
    pos0: u32,
    strand: Strand,
    mod_code: String,
}

impl FullRowKey {
    fn summary_key(&self) -> SummaryKey {
        SummaryKey {
            assay_id: self.assay_id.clone(),
            sample: self.sample.clone(),
            gene: self.gene.clone(),
            chrom: self.chrom.clone(),
            pos0: self.pos0,
            strand: self.strand,
            mod_code: self.mod_code.clone(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SummaryKey {
    assay_id: String,
    sample: String,
    gene: String,
    chrom: String,
    pos0: u32,
    strand: Strand,
    mod_code: String,
}

impl SummaryKey {
    fn site_id(&self) -> String {
        format!("{}:{}:{}", self.chrom, self.pos0, self.strand.as_char())
    }
}

#[derive(Clone, Debug)]
struct ParsedSiteRow {
    key: FullRowKey,
    analysis_threshold: f64,
    group: Option<String>,
    context: Option<String>,
    site_state: ParsedSiteState,
    coverage_basis: ParsedCoverageBasis,
    n_assigned: u64,
    n_covering: Option<u64>,
    n_callable: u64,
    eligibility_reason: ParsedEligibilityReason,
    line: u64,
}

fn required<'a>(record: &'a StringRecord, index: usize, field: &str) -> anyhow::Result<&'a str> {
    let value = record.get(index).unwrap_or_default();
    if value.trim().is_empty() || value == "NA" || value.chars().any(char::is_control) {
        anyhow::bail!(
            "{field} must not be empty, whitespace-only, NA, or contain control characters"
        );
    }
    Ok(value)
}

fn optional_token(
    record: &StringRecord,
    index: usize,
    field: &str,
) -> anyhow::Result<Option<String>> {
    let value = record.get(index).unwrap_or_default();
    match value {
        "NA" => Ok(None),
        "" => anyhow::bail!("{field} must use NA for a missing value"),
        value if value.trim().is_empty() || value.chars().any(char::is_control) => {
            anyhow::bail!("{field} must not be whitespace-only or contain control characters")
        }
        value => Ok(Some(value.to_owned())),
    }
}

fn parse_u64(record: &StringRecord, index: usize, field: &str) -> anyhow::Result<u64> {
    required(record, index, field)?
        .parse::<u64>()
        .with_context(|| format!("invalid {field}"))
}

fn parse_optional_u64(
    record: &StringRecord,
    index: usize,
    field: &str,
) -> anyhow::Result<Option<u64>> {
    match record.get(index).unwrap_or_default() {
        "NA" => Ok(None),
        "" => anyhow::bail!("{field} must use NA for a missing value"),
        value => value
            .parse::<u64>()
            .map(Some)
            .with_context(|| format!("invalid {field}")),
    }
}

fn parse_unit_interval(value: &str, field: &str) -> anyhow::Result<f64> {
    let value = value
        .parse::<f64>()
        .with_context(|| format!("invalid {field}"))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("{field} must be finite and in [0, 1]");
    }
    Ok(if value == 0.0 { 0.0 } else { value })
}

fn parse_optional_unit_interval(
    record: &StringRecord,
    index: usize,
    field: &str,
) -> anyhow::Result<Option<f64>> {
    match record.get(index).unwrap_or_default() {
        "NA" => Ok(None),
        "" => anyhow::bail!("{field} must use NA for a missing value"),
        value => parse_unit_interval(value, field).map(Some),
    }
}

fn validate_mod_code(value: &str) -> anyhow::Result<()> {
    let bytes = value.as_bytes();
    if bytes.len() < 3
        || !matches!(bytes[0], b'A' | b'C' | b'G' | b'T' | b'U' | b'N')
        || !matches!(bytes[1], b'+' | b'-')
        || !bytes[2..].iter().all(u8::is_ascii_alphanumeric)
    {
        anyhow::bail!(
            "mod_code must contain canonical base, +/- strand, and an ASCII alphanumeric code"
        );
    }
    Ok(())
}

fn parse_site_row(record: &StringRecord, line: u64) -> anyhow::Result<ParsedSiteRow> {
    if record.len() != SITE_COLUMNS.len() {
        anyhow::bail!(
            "site row has {} columns; expected {}",
            record.len(),
            SITE_COLUMNS.len()
        );
    }

    let assay_id = required(record, 0, "assay_id")?.to_owned();
    let analysis_threshold = parse_unit_interval(
        required(record, 1, "analysis_threshold")?,
        "analysis_threshold",
    )?;
    let sample = required(record, 2, "sample")?.to_owned();
    let group = optional_token(record, 3, "group")?;
    let gene = required(record, 4, "gene")?.to_owned();
    let isoform_id = required(record, 5, "isoform_id")?.to_owned();
    let site_id = required(record, 6, "site_id")?;
    let chrom = required(record, 7, "chrom")?.to_owned();
    let pos0 = required(record, 8, "pos0")?
        .parse::<u32>()
        .context("invalid pos0")?;
    let strand = Strand::try_from(required(record, 9, "strand")?)?;
    if strand == Strand::Unknown {
        anyhow::bail!("site strand must be + or -");
    }
    let mod_code = required(record, 10, "mod_code")?.to_owned();
    validate_mod_code(&mod_code)?;
    let expected_site_id = format!("{chrom}:{pos0}:{}", strand.as_char());
    if site_id != expected_site_id {
        anyhow::bail!("site_id {site_id:?} does not match {expected_site_id:?}");
    }
    let context = optional_token(record, 11, "context")?;
    let site_state = ParsedSiteState::parse(required(record, 12, "site_state")?)?;
    let coverage_basis = ParsedCoverageBasis::parse(required(record, 13, "coverage_basis")?)?;
    let n_assigned = parse_u64(record, 14, "n_assigned")?;
    let n_covering = parse_optional_u64(record, 15, "n_covering")?;
    let n_candidate = parse_u64(record, 16, "n_candidate")?;
    let n_callable = parse_u64(record, 17, "n_callable")?;
    let n_modified = parse_u64(record, 18, "n_modified")?;
    let n_unmodified = parse_u64(record, 19, "n_unmodified")?;
    let n_unknown = parse_u64(record, 20, "n_unknown")?;
    let mod_fraction = parse_optional_unit_interval(record, 21, "mod_fraction")?;
    let _mean_probability = parse_optional_unit_interval(record, 22, "mean_probability")?;
    let ci_low = parse_optional_unit_interval(record, 23, "ci_low")?;
    let ci_high = parse_optional_unit_interval(record, 24, "ci_high")?;
    let eligibility = required(record, 25, "eligibility")?;
    let eligibility_reason =
        ParsedEligibilityReason::parse(required(record, 26, "eligibility_reason")?)?;

    if n_modified
        .checked_add(n_unmodified)
        .context("n_modified + n_unmodified overflows u64")?
        != n_callable
    {
        anyhow::bail!("n_callable must equal n_modified + n_unmodified");
    }
    if n_callable
        .checked_add(n_unknown)
        .context("n_callable + n_unknown overflows u64")?
        != n_candidate
    {
        anyhow::bail!("n_candidate must equal n_callable + n_unknown");
    }
    if n_candidate > n_assigned {
        anyhow::bail!("n_candidate must not exceed n_assigned");
    }
    match (coverage_basis, n_covering) {
        (ParsedCoverageBasis::Unavailable, None) => {}
        (ParsedCoverageBasis::Unavailable, Some(_)) => {
            anyhow::bail!("coverage_basis=unavailable requires n_covering=NA")
        }
        (_, None) => anyhow::bail!("numeric coverage basis requires numeric n_covering"),
        (_, Some(covering)) if covering > n_assigned => {
            anyhow::bail!("n_covering must not exceed n_assigned")
        }
        (_, Some(covering)) if n_candidate > covering => {
            anyhow::bail!("n_candidate must not exceed n_covering")
        }
        _ => {}
    }

    if let Some(fraction) = mod_fraction {
        if n_callable == 0 {
            anyhow::bail!("numeric mod_fraction requires n_callable > 0");
        }
        let expected = n_modified as f64 / n_callable as f64;
        if (fraction - expected).abs() > 1e-12 {
            anyhow::bail!("mod_fraction is inconsistent with integer counts");
        }
    }
    match (mod_fraction, ci_low, ci_high) {
        (Some(fraction), Some(low), Some(high))
            if low <= fraction && fraction <= high && low <= high => {}
        (Some(_), _, _) => {
            anyhow::bail!("numeric mod_fraction requires an ordered CI containing the fraction")
        }
        (None, None, None) => {}
        (None, _, _) => anyhow::bail!("NA mod_fraction requires ci_low=NA and ci_high=NA"),
    }

    let expected_reason = match site_state {
        ParsedSiteState::StructurallyAbsent => Some("site_absent"),
        ParsedSiteState::ContextDependent => Some("context_dependent"),
        ParsedSiteState::ReferenceBaseMismatch => Some("reference_base_mismatch"),
        ParsedSiteState::Unprojectable => Some("unprojectable"),
        ParsedSiteState::Present => None,
    };
    if expected_reason.is_some_and(|expected| eligibility_reason.as_str() != expected) {
        anyhow::bail!("site_state and eligibility_reason are inconsistent");
    }
    if site_state == ParsedSiteState::Present
        && matches!(
            eligibility_reason.as_str(),
            "site_absent" | "context_dependent" | "reference_base_mismatch" | "unprojectable"
        )
    {
        anyhow::bail!("present site has incompatible eligibility_reason");
    }
    match eligibility {
        "eligible" if eligibility_reason.is_eligible() => {}
        "ineligible" if !eligibility_reason.is_eligible() => {}
        "eligible" | "ineligible" => {
            anyhow::bail!("eligibility and eligibility_reason are inconsistent")
        }
        _ => anyhow::bail!("eligibility must be eligible or ineligible"),
    }
    if eligibility_reason.as_str() == "ok" && (n_callable == 0 || mod_fraction.is_none()) {
        anyhow::bail!("eligible row requires a callable denominator and mod_fraction");
    }
    if eligibility_reason.as_str() == "incomplete_candidate_universe" && mod_fraction.is_some() {
        anyhow::bail!("incomplete candidate universe requires mod_fraction=NA");
    }

    Ok(ParsedSiteRow {
        key: FullRowKey {
            assay_id,
            sample,
            gene,
            isoform_id,
            chrom,
            pos0,
            strand,
            mod_code,
        },
        analysis_threshold,
        group,
        context,
        site_state,
        coverage_basis,
        n_assigned,
        n_covering,
        n_callable,
        eligibility_reason,
        line,
    })
}

type SiteCsvReader = Reader<BufReader<File>>;

struct SiteStream {
    path: PathBuf,
    reader: SiteCsvReader,
    current: Option<ParsedSiteRow>,
    previous_key: Option<FullRowKey>,
}

impl SiteStream {
    fn open(path: &Path) -> anyhow::Result<Self> {
        let file = File::open(path).with_context(|| format!("open site table {path:?}"))?;
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(false)
            .from_reader(BufReader::new(file));
        let expected = StringRecord::from(SITE_COLUMNS.to_vec());
        let actual = reader
            .headers()
            .with_context(|| format!("read site table header {path:?}"))?;
        if actual != &expected {
            anyhow::bail!(
                "site table {path:?} header mismatch; expected {:?}",
                SITE_COLUMNS
            );
        }
        let mut stream = Self {
            path: path.to_path_buf(),
            reader,
            current: None,
            previous_key: None,
        };
        stream.advance()?;
        Ok(stream)
    }

    fn advance(&mut self) -> anyhow::Result<()> {
        let mut record = StringRecord::new();
        if !self
            .reader
            .read_record(&mut record)
            .with_context(|| format!("read site table {:?}", self.path))?
        {
            self.current = None;
            return Ok(());
        }
        let line = record.position().map_or(0, |position| position.line());
        let row = parse_site_row(&record, line)
            .with_context(|| format!("parse site table {:?}:{line}", self.path))?;
        if let Some(previous) = &self.previous_key {
            match row.key.cmp(previous) {
                std::cmp::Ordering::Less => anyhow::bail!(
                    "site table {:?}:{line} is not sorted by assay, sample, gene, isoform, and genomic site",
                    self.path
                ),
                std::cmp::Ordering::Equal => anyhow::bail!(
                    "duplicate isoform-site key in site table {:?}:{line}",
                    self.path
                ),
                std::cmp::Ordering::Greater => {}
            }
        }
        self.previous_key = Some(row.key.clone());
        self.current = Some(row);
        Ok(())
    }
}

#[derive(Default)]
struct GlobalAudit {
    assay_thresholds: BTreeMap<String, f64>,
    sample_groups: BTreeMap<(String, String), Option<String>>,
    gene_coverage_basis: BTreeMap<(String, String, String), ParsedCoverageBasis>,
    assigned_counts: BTreeMap<(String, String, String, String), u64>,
}

impl GlobalAudit {
    fn validate(&mut self, row: &ParsedSiteRow) -> anyhow::Result<()> {
        match self
            .assay_thresholds
            .insert(row.key.assay_id.clone(), row.analysis_threshold)
        {
            Some(existing) if existing != row.analysis_threshold => anyhow::bail!(
                "assay {:?} has inconsistent analysis thresholds",
                row.key.assay_id
            ),
            _ => {}
        }
        let sample_key = (row.key.assay_id.clone(), row.key.sample.clone());
        match self.sample_groups.insert(sample_key, row.group.clone()) {
            Some(existing) if existing != row.group => anyhow::bail!(
                "assay {:?}, sample {:?} has inconsistent groups",
                row.key.assay_id,
                row.key.sample
            ),
            _ => {}
        }
        let coverage_key = (
            row.key.assay_id.clone(),
            row.key.sample.clone(),
            row.key.gene.clone(),
        );
        match self
            .gene_coverage_basis
            .insert(coverage_key, row.coverage_basis)
        {
            Some(existing) if existing != row.coverage_basis => anyhow::bail!(
                "assay {:?}, sample {:?}, gene {:?} has inconsistent coverage_basis",
                row.key.assay_id,
                row.key.sample,
                row.key.gene
            ),
            _ => {}
        }
        let assigned_key = (
            row.key.assay_id.clone(),
            row.key.sample.clone(),
            row.key.gene.clone(),
            row.key.isoform_id.clone(),
        );
        match self.assigned_counts.insert(assigned_key, row.n_assigned) {
            Some(existing) if existing != row.n_assigned => anyhow::bail!(
                "isoform {:?} has inconsistent n_assigned values",
                row.key.isoform_id
            ),
            _ => {}
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct SummaryAccumulator {
    analysis_threshold: f64,
    group: Option<String>,
    context: Option<String>,
    coverage_basis: ParsedCoverageBasis,
    n_isoforms_total: u64,
    n_isoforms_assigned: u64,
    n_isoforms_present: u64,
    n_isoforms_eligible: u64,
    n_isoforms_site_absent: u64,
    n_isoforms_context_dependent: u64,
    n_isoforms_reference_base_mismatch: u64,
    n_isoforms_unprojectable: u64,
    n_isoforms_incomplete_candidate_universe: u64,
    n_isoforms_join_rate_low: u64,
    n_isoforms_low_callable: u64,
    n_isoforms_other_ineligible: u64,
    min_eligible_n_covering: Option<u64>,
    min_eligible_n_callable: Option<u64>,
}

fn increment(value: &mut u64, field: &str) -> anyhow::Result<()> {
    *value = value
        .checked_add(1)
        .with_context(|| format!("{field} overflows u64"))?;
    Ok(())
}

fn update_min(target: &mut Option<u64>, value: u64) {
    *target = Some(target.map_or(value, |current| current.min(value)));
}

impl SummaryAccumulator {
    fn new(row: &ParsedSiteRow) -> Self {
        Self {
            analysis_threshold: row.analysis_threshold,
            group: row.group.clone(),
            context: row.context.clone(),
            coverage_basis: row.coverage_basis,
            n_isoforms_total: 0,
            n_isoforms_assigned: 0,
            n_isoforms_present: 0,
            n_isoforms_eligible: 0,
            n_isoforms_site_absent: 0,
            n_isoforms_context_dependent: 0,
            n_isoforms_reference_base_mismatch: 0,
            n_isoforms_unprojectable: 0,
            n_isoforms_incomplete_candidate_universe: 0,
            n_isoforms_join_rate_low: 0,
            n_isoforms_low_callable: 0,
            n_isoforms_other_ineligible: 0,
            min_eligible_n_covering: None,
            min_eligible_n_callable: None,
        }
    }

    fn add(&mut self, row: &ParsedSiteRow) -> anyhow::Result<()> {
        if self.analysis_threshold != row.analysis_threshold {
            anyhow::bail!("summary site has inconsistent analysis_threshold");
        }
        if self.group != row.group {
            anyhow::bail!("summary site has inconsistent group");
        }
        if self.context != row.context {
            anyhow::bail!("summary site has conflicting context");
        }
        if self.coverage_basis != row.coverage_basis {
            anyhow::bail!("summary site has inconsistent coverage_basis");
        }

        increment(&mut self.n_isoforms_total, "n_isoforms_total")?;
        if row.n_assigned > 0 {
            increment(&mut self.n_isoforms_assigned, "n_isoforms_assigned")?;
        }
        match row.site_state {
            ParsedSiteState::Present => {
                increment(&mut self.n_isoforms_present, "n_isoforms_present")?
            }
            ParsedSiteState::StructurallyAbsent => {
                increment(&mut self.n_isoforms_site_absent, "n_isoforms_site_absent")?
            }
            ParsedSiteState::ContextDependent => increment(
                &mut self.n_isoforms_context_dependent,
                "n_isoforms_context_dependent",
            )?,
            ParsedSiteState::ReferenceBaseMismatch => increment(
                &mut self.n_isoforms_reference_base_mismatch,
                "n_isoforms_reference_base_mismatch",
            )?,
            ParsedSiteState::Unprojectable => increment(
                &mut self.n_isoforms_unprojectable,
                "n_isoforms_unprojectable",
            )?,
        }
        match row.eligibility_reason.as_str() {
            "ok" => {
                increment(&mut self.n_isoforms_eligible, "n_isoforms_eligible")?;
                update_min(&mut self.min_eligible_n_callable, row.n_callable);
                if let Some(n_covering) = row.n_covering {
                    update_min(&mut self.min_eligible_n_covering, n_covering);
                }
            }
            "incomplete_candidate_universe" => increment(
                &mut self.n_isoforms_incomplete_candidate_universe,
                "n_isoforms_incomplete_candidate_universe",
            )?,
            "join_rate_low" => increment(
                &mut self.n_isoforms_join_rate_low,
                "n_isoforms_join_rate_low",
            )?,
            "low_callable" => {
                increment(&mut self.n_isoforms_low_callable, "n_isoforms_low_callable")?
            }
            "site_absent" | "context_dependent" | "reference_base_mismatch" | "unprojectable" => {}
            _ => increment(
                &mut self.n_isoforms_other_ineligible,
                "n_isoforms_other_ineligible",
            )?,
        }
        Ok(())
    }

    fn validate(&self) -> anyhow::Result<()> {
        let state_total = self
            .n_isoforms_present
            .checked_add(self.n_isoforms_site_absent)
            .and_then(|value| value.checked_add(self.n_isoforms_context_dependent))
            .and_then(|value| value.checked_add(self.n_isoforms_reference_base_mismatch))
            .and_then(|value| value.checked_add(self.n_isoforms_unprojectable))
            .context("site-state summary counts overflow u64")?;
        if state_total != self.n_isoforms_total {
            anyhow::bail!("site-state summary counts do not equal n_isoforms_total");
        }
        let present_total = self
            .n_isoforms_eligible
            .checked_add(self.n_isoforms_incomplete_candidate_universe)
            .and_then(|value| value.checked_add(self.n_isoforms_join_rate_low))
            .and_then(|value| value.checked_add(self.n_isoforms_low_callable))
            .and_then(|value| value.checked_add(self.n_isoforms_other_ineligible))
            .context("present-site eligibility counts overflow u64")?;
        if present_total != self.n_isoforms_present {
            anyhow::bail!("present-site eligibility counts do not equal n_isoforms_present");
        }
        if self.n_isoforms_eligible == 0 {
            if self.min_eligible_n_callable.is_some() || self.min_eligible_n_covering.is_some() {
                anyhow::bail!("site without eligible isoforms has eligible minima");
            }
        } else if self.min_eligible_n_callable.is_none() {
            anyhow::bail!("eligible site is missing min_eligible_n_callable");
        }
        if self.coverage_basis == ParsedCoverageBasis::Unavailable
            && self.min_eligible_n_covering.is_some()
        {
            anyhow::bail!("unavailable coverage cannot have min_eligible_n_covering");
        }
        Ok(())
    }
}

pub(crate) struct SiteSummaryResult {
    input_rows: u64,
    summaries: BTreeMap<SummaryKey, SummaryAccumulator>,
}

impl SiteSummaryResult {
    pub(crate) fn input_rows(&self) -> u64 {
        self.input_rows
    }

    pub(crate) fn site_count(&self) -> usize {
        self.summaries.len()
    }
}

fn queue_insert(
    queue: &mut BTreeMap<FullRowKey, Vec<usize>>,
    streams: &[SiteStream],
    stream_index: usize,
) {
    if let Some(row) = streams[stream_index].current.as_ref() {
        queue.entry(row.key.clone()).or_default().push(stream_index);
    }
}

fn queue_pop(queue: &mut BTreeMap<FullRowKey, Vec<usize>>) -> Option<(FullRowKey, usize)> {
    let key = queue.keys().next()?.clone();
    let (stream_index, remove) = {
        let indices = queue.get_mut(&key).expect("key came from queue");
        let stream_index = indices.pop().expect("queue entries are non-empty");
        (stream_index, indices.is_empty())
    };
    if remove {
        queue.remove(&key);
    }
    Some((key, stream_index))
}

pub(crate) fn summarize_site_files(paths: &[PathBuf]) -> anyhow::Result<SiteSummaryResult> {
    if paths.is_empty() {
        anyhow::bail!("at least one site table is required");
    }
    let mut streams = paths
        .iter()
        .map(|path| SiteStream::open(path))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut queue = BTreeMap::new();
    for stream_index in 0..streams.len() {
        queue_insert(&mut queue, &streams, stream_index);
    }

    let mut audit = GlobalAudit::default();
    let mut summaries = BTreeMap::new();
    let mut input_rows = 0u64;
    let mut previous_merged_key: Option<FullRowKey> = None;
    let mut previous_source: Option<(PathBuf, u64)> = None;

    while let Some((queued_key, stream_index)) = queue_pop(&mut queue) {
        let row = streams[stream_index]
            .current
            .take()
            .expect("queued stream has a current row");
        debug_assert_eq!(queued_key, row.key);
        if previous_merged_key.as_ref() == Some(&row.key) {
            let (previous_path, previous_line) =
                previous_source.as_ref().expect("previous key has source");
            anyhow::bail!(
                "duplicate isoform-site key across inputs at {:?}:{} and {:?}:{}",
                previous_path,
                previous_line,
                streams[stream_index].path,
                row.line
            );
        }

        audit.validate(&row).with_context(|| {
            format!(
                "validate site table {:?}:{}",
                streams[stream_index].path, row.line
            )
        })?;
        let summary_key = row.key.summary_key();
        match summaries.entry(summary_key) {
            Entry::Vacant(entry) => {
                let mut accumulator = SummaryAccumulator::new(&row);
                accumulator.add(&row)?;
                entry.insert(accumulator);
            }
            Entry::Occupied(mut entry) => entry.get_mut().add(&row)?,
        }
        input_rows = input_rows
            .checked_add(1)
            .context("input row count overflows u64")?;
        previous_merged_key = Some(row.key);
        previous_source = Some((streams[stream_index].path.clone(), row.line));

        streams[stream_index].advance()?;
        queue_insert(&mut queue, &streams, stream_index);
    }

    if input_rows == 0 {
        anyhow::bail!("site tables contain no data rows");
    }
    for accumulator in summaries.values() {
        accumulator.validate()?;
    }
    Ok(SiteSummaryResult {
        input_rows,
        summaries,
    })
}

fn optional<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

pub(crate) fn write_site_summary_tsv<W: Write>(
    writer: W,
    result: &SiteSummaryResult,
) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record(SITE_SUMMARY_COLUMNS)?;
    for (key, summary) in &result.summaries {
        summary.validate()?;
        output.write_record([
            key.assay_id.clone(),
            summary.analysis_threshold.to_string(),
            key.sample.clone(),
            optional(summary.group.clone()),
            key.gene.clone(),
            key.site_id(),
            key.chrom.clone(),
            key.pos0.to_string(),
            key.strand.as_char().to_string(),
            key.mod_code.clone(),
            optional(summary.context.clone()),
            summary.coverage_basis.as_str().to_owned(),
            summary.n_isoforms_total.to_string(),
            summary.n_isoforms_assigned.to_string(),
            summary.n_isoforms_present.to_string(),
            summary.n_isoforms_eligible.to_string(),
            summary.n_isoforms_site_absent.to_string(),
            summary.n_isoforms_context_dependent.to_string(),
            summary.n_isoforms_reference_base_mismatch.to_string(),
            summary.n_isoforms_unprojectable.to_string(),
            summary.n_isoforms_incomplete_candidate_universe.to_string(),
            summary.n_isoforms_join_rate_low.to_string(),
            summary.n_isoforms_low_callable.to_string(),
            summary.n_isoforms_other_ineligible.to_string(),
            optional(summary.min_eligible_n_covering),
            optional(summary.min_eligible_n_callable),
            SummaryState::from_count(summary.n_isoforms_eligible)
                .as_str()
                .to_owned(),
        ])?;
    }
    output.flush().context("flush modification site summary")?;
    Ok(())
}
