mod common;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use common::{assert_success, TestDir};

const SITE_HEADER: &str = concat!(
    "assay_id\tanalysis_threshold\tsample\tgroup\tgene\tisoform_id\tsite_id\t",
    "chrom\tpos0\tstrand\tmod_code\tcontext\tsite_state\tcoverage_basis\t",
    "eligibility_profile\tn_assigned\tn_covering\tn_not_candidate\tn_candidate\t",
    "candidate_rate\tn_callable\tcallable_rate\tn_modified\tn_unmodified\t",
    "n_unknown\tmod_fraction\tmean_probability\tci_low\tci_high\teligibility\t",
    "eligibility_reason\n"
);

const V0_3_0_SITE_HEADER: &str = concat!(
    "assay_id\tanalysis_threshold\tsample\tgroup\tgene\tisoform_id\tsite_id\t",
    "chrom\tpos0\tstrand\tmod_code\tcontext\tsite_state\tcoverage_basis\t",
    "n_assigned\tn_covering\tn_candidate\tn_callable\tn_modified\t",
    "n_unmodified\tn_unknown\tmod_fraction\tmean_probability\tci_low\tci_high\t",
    "eligibility\teligibility_reason\n"
);

const SUMMARY_HEADER: &str = concat!(
    "assay_id\tanalysis_threshold\tsample\tgroup\tgene\tsite_id\tchrom\tpos0\t",
    "strand\tmod_code\tcontext\tcoverage_basis\teligibility_profile\t",
    "n_isoforms_total\t",
    "n_isoforms_assigned\tn_isoforms_present\tn_isoforms_eligible\t",
    "n_isoforms_site_absent\tn_isoforms_context_dependent\t",
    "n_isoforms_reference_base_mismatch\tn_isoforms_unprojectable\t",
    "n_isoforms_incomplete_candidate_universe\tn_isoforms_join_rate_low\t",
    "n_isoforms_site_join_rate_low\tn_isoforms_unknown_denominator\t",
    "n_isoforms_coverage_unavailable\tn_isoforms_reference_unvalidated\t",
    "n_isoforms_low_covering\tn_isoforms_low_callable\t",
    "n_isoforms_low_candidate_rate\tn_isoforms_low_callable_rate\t",
    "n_isoforms_provenance_unverified\tn_isoforms_other_ineligible\t",
    "min_eligible_n_covering\tmin_eligible_n_callable\tsummary_state\n"
);

#[derive(Clone, Debug)]
struct SiteRow {
    assay_id: &'static str,
    analysis_threshold: &'static str,
    sample: &'static str,
    group: &'static str,
    gene: &'static str,
    isoform_id: &'static str,
    chrom: &'static str,
    pos0: u32,
    strand: &'static str,
    mod_code: &'static str,
    context: &'static str,
    site_state: &'static str,
    coverage_basis: &'static str,
    eligibility_profile: &'static str,
    n_assigned: u64,
    n_covering: Option<u64>,
    n_candidate: u64,
    n_callable: u64,
    n_modified: u64,
    n_unmodified: u64,
    n_unknown: u64,
    mod_fraction: Option<&'static str>,
    mean_probability: Option<&'static str>,
    ci_low: Option<&'static str>,
    ci_high: Option<&'static str>,
    eligibility: &'static str,
    eligibility_reason: &'static str,
}

impl SiteRow {
    fn new(gene: &'static str, isoform_id: &'static str, pos0: u32) -> Self {
        Self {
            assay_id: "a1",
            analysis_threshold: "0.5",
            sample: "S1",
            group: "control",
            gene,
            isoform_id,
            chrom: "chr1",
            pos0,
            strand: "+",
            mod_code: "A+a",
            context: "DRACH",
            site_state: "present",
            coverage_basis: "bam_exact",
            eligibility_profile: "exploratory",
            n_assigned: 1,
            n_covering: Some(0),
            n_candidate: 0,
            n_callable: 0,
            n_modified: 0,
            n_unmodified: 0,
            n_unknown: 0,
            mod_fraction: None,
            mean_probability: None,
            ci_low: None,
            ci_high: None,
            eligibility: "ineligible",
            eligibility_reason: "low_callable",
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn exact_eligible(
        mut self,
        n_assigned: u64,
        n_covering: u64,
        n_modified: u64,
        n_unmodified: u64,
        fraction: &'static str,
        ci_low: &'static str,
        ci_high: &'static str,
    ) -> Self {
        self.n_assigned = n_assigned;
        self.n_covering = Some(n_covering);
        self.n_candidate = n_modified + n_unmodified;
        self.n_callable = self.n_candidate;
        self.n_modified = n_modified;
        self.n_unmodified = n_unmodified;
        self.mod_fraction = Some(fraction);
        self.mean_probability = Some(fraction);
        self.ci_low = Some(ci_low);
        self.ci_high = Some(ci_high);
        self.eligibility = "eligible";
        self.eligibility_reason = "ok";
        self
    }

    fn unavailable_eligible(
        mut self,
        n_assigned: u64,
        n_modified: u64,
        n_unmodified: u64,
        fraction: &'static str,
    ) -> Self {
        self = self.exact_eligible(
            n_assigned,
            n_assigned,
            n_modified,
            n_unmodified,
            fraction,
            "0",
            "1",
        );
        self.coverage_basis = "unavailable";
        self.n_covering = None;
        self
    }

    fn zero_count_reason(
        mut self,
        n_assigned: u64,
        n_covering: u64,
        site_state: &'static str,
        reason: &'static str,
    ) -> Self {
        self.n_assigned = n_assigned;
        self.n_covering = Some(n_covering);
        self.site_state = site_state;
        self.eligibility_reason = reason;
        self
    }

    fn render(&self) -> String {
        let site_id = format!("{}:{}:{}", self.chrom, self.pos0, self.strand);
        let (n_not_candidate, candidate_rate, callable_rate) = match self.n_covering {
            None => (None, None, None),
            Some(0) => (Some(0), None, None),
            Some(covering) => (
                Some(covering - self.n_candidate),
                Some(self.n_candidate as f64 / covering as f64),
                Some(self.n_callable as f64 / covering as f64),
            ),
        };
        [
            self.assay_id.to_owned(),
            self.analysis_threshold.to_owned(),
            self.sample.to_owned(),
            self.group.to_owned(),
            self.gene.to_owned(),
            self.isoform_id.to_owned(),
            site_id,
            self.chrom.to_owned(),
            self.pos0.to_string(),
            self.strand.to_owned(),
            self.mod_code.to_owned(),
            self.context.to_owned(),
            self.site_state.to_owned(),
            self.coverage_basis.to_owned(),
            self.eligibility_profile.to_owned(),
            self.n_assigned.to_string(),
            optional(self.n_covering),
            optional(n_not_candidate),
            self.n_candidate.to_string(),
            optional(candidate_rate),
            self.n_callable.to_string(),
            optional(callable_rate),
            self.n_modified.to_string(),
            self.n_unmodified.to_string(),
            self.n_unknown.to_string(),
            optional(self.mod_fraction),
            optional(self.mean_probability),
            optional(self.ci_low),
            optional(self.ci_high),
            self.eligibility.to_owned(),
            self.eligibility_reason.to_owned(),
        ]
        .join("\t")
    }

    fn render_v0_3_0(&self) -> String {
        let site_id = format!("{}:{}:{}", self.chrom, self.pos0, self.strand);
        [
            self.assay_id.to_owned(),
            self.analysis_threshold.to_owned(),
            self.sample.to_owned(),
            self.group.to_owned(),
            self.gene.to_owned(),
            self.isoform_id.to_owned(),
            site_id,
            self.chrom.to_owned(),
            self.pos0.to_string(),
            self.strand.to_owned(),
            self.mod_code.to_owned(),
            self.context.to_owned(),
            self.site_state.to_owned(),
            self.coverage_basis.to_owned(),
            self.n_assigned.to_string(),
            optional(self.n_covering),
            self.n_candidate.to_string(),
            self.n_callable.to_string(),
            self.n_modified.to_string(),
            self.n_unmodified.to_string(),
            self.n_unknown.to_string(),
            optional(self.mod_fraction),
            optional(self.mean_probability),
            optional(self.ci_low),
            optional(self.ci_high),
            self.eligibility.to_owned(),
            self.eligibility_reason.to_owned(),
        ]
        .join("\t")
    }
}

fn optional<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

fn write_sites(path: &Path, rows: &[SiteRow]) {
    let mut text = SITE_HEADER.to_owned();
    for row in rows {
        text.push_str(&row.render());
        text.push('\n');
    }
    fs::write(path, text).unwrap();
}

fn summary_path(prefix: &Path) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(".mod_site_summary.tsv");
    PathBuf::from(value)
}

fn run_summary(inputs: &[&Path], prefix: &Path) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    command.arg("mod-site-summary");
    for input in inputs {
        command.arg("--sites").arg(input);
    }
    command.arg("--out").arg(prefix).output().unwrap()
}

#[test]
fn accepts_v0_3_0_site_tables_as_exploratory_inputs() {
    let root = TestDir::new("mod-site-summary-v0-3-0");
    let sites = root.join("legacy.tsv");
    let row = SiteRow::new("G1", "iso1", 110).exact_eligible(2, 2, 1, 1, "0.5", "0", "1");
    fs::write(
        &sites,
        format!("{V0_3_0_SITE_HEADER}{}\n", row.render_v0_3_0()),
    )
    .unwrap();

    let prefix = root.join("summary");
    let output = run_summary(&[sites.as_path()], &prefix);
    assert_success(&output, "v0.3.0 site summary");

    let text = fs::read_to_string(summary_path(&prefix)).unwrap();
    let mut lines = text.lines();
    let columns = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let values = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let value = |column: &str| {
        let index = columns
            .iter()
            .position(|candidate| *candidate == column)
            .unwrap();
        values[index]
    };
    assert_eq!(value("eligibility_profile"), "exploratory");
    assert_eq!(value("n_isoforms_eligible"), "1");
    assert_eq!(value("min_eligible_n_covering"), "2");
    assert_eq!(value("min_eligible_n_callable"), "2");
}

#[test]
fn summarizes_sites_deterministically_with_all_qc_states() {
    let root = TestDir::new("mod-site-summary-success");
    let left = root.join("left.tsv");
    let right = root.join("right.tsv");

    let mut incomplete = SiteRow::new("G1", "iso4", 130);
    incomplete.n_assigned = 3;
    incomplete.n_covering = Some(1);
    incomplete.n_candidate = 1;
    incomplete.n_callable = 1;
    incomplete.n_modified = 1;
    incomplete.mean_probability = Some("0.9");
    incomplete.eligibility_reason = "incomplete_candidate_universe";

    write_sites(
        &left,
        &[
            SiteRow::new("G1", "iso1", 110).exact_eligible(10, 8, 3, 2, "0.6", "0.2", "0.9"),
            SiteRow::new("G1", "iso1", 120).zero_count_reason(10, 5, "present", "low_callable"),
            SiteRow::new("G1", "iso1", 130).zero_count_reason(
                10,
                0,
                "context_dependent",
                "context_dependent",
            ),
            SiteRow::new("G1", "iso1", 140).exact_eligible(10, 4, 1, 1, "0.5", "0.1", "0.9"),
            SiteRow::new("G2", "isoA", 10).unavailable_eligible(2, 1, 0, "1"),
        ],
    );
    write_sites(
        &right,
        &[
            SiteRow::new("G1", "iso2", 110).exact_eligible(7, 6, 1, 3, "0.25", "0.05", "0.7"),
            SiteRow::new("G1", "iso2", 120).zero_count_reason(
                7,
                0,
                "structurally_absent",
                "site_absent",
            ),
            SiteRow::new("G1", "iso2", 130).zero_count_reason(
                7,
                0,
                "reference_base_mismatch",
                "reference_base_mismatch",
            ),
            SiteRow::new("G1", "iso2", 140).zero_count_reason(7, 3, "present", "low_callable"),
            SiteRow::new("G1", "iso3", 130).zero_count_reason(
                3,
                0,
                "unprojectable",
                "unprojectable",
            ),
            incomplete,
            SiteRow::new("G1", "iso5", 130).zero_count_reason(0, 0, "present", "join_rate_low"),
        ],
    );

    let first_prefix = root.join("first");
    let second_prefix = root.join("second");
    let first = run_summary(&[left.as_path(), right.as_path()], &first_prefix);
    assert_success(&first, "site summary");
    assert_eq!(
        String::from_utf8(first.stderr).unwrap(),
        "mod-site-summary: inputs=2 input_rows=12 sites=5\n"
    );
    let second = run_summary(&[right.as_path(), left.as_path()], &second_prefix);
    assert_success(&second, "reordered site summary");

    let first_text = fs::read_to_string(summary_path(&first_prefix)).unwrap();
    let second_text = fs::read_to_string(summary_path(&second_prefix)).unwrap();
    assert_eq!(first_text, second_text);
    assert!(first_text.starts_with(SUMMARY_HEADER));
    assert_eq!(first_text.lines().count(), 6);
    assert!(first_text.contains("\tshared_eligible\n"));
    assert!(first_text.contains("\tsingle_eligible\n"));
    assert!(first_text.contains("\tno_eligible_isoform\n"));
}

#[test]
fn counts_new_and_future_eligibility_reasons_in_catch_all() {
    let root = TestDir::new("mod-site-summary-open-reasons");
    let sites = root.join("sites.tsv");

    let site_join =
        SiteRow::new("G1", "iso1", 110).zero_count_reason(1, 0, "present", "site_join_rate_low");
    let mut unknown_denominator =
        SiteRow::new("G1", "iso2", 110).zero_count_reason(1, 1, "present", "unknown_denominator");
    unknown_denominator.n_candidate = 1;
    unknown_denominator.n_unknown = 1;
    let coverage_unavailable =
        SiteRow::new("G1", "iso3", 110).zero_count_reason(1, 0, "present", "coverage_unavailable");
    let reference_unvalidated =
        SiteRow::new("G1", "iso4", 110).zero_count_reason(1, 1, "present", "reference_unvalidated");
    let low_covering =
        SiteRow::new("G1", "iso5", 110).zero_count_reason(1, 1, "present", "low_covering");
    let low_candidate_rate =
        SiteRow::new("G1", "iso6", 110).zero_count_reason(1, 1, "present", "low_candidate_rate");
    let low_callable_rate =
        SiteRow::new("G1", "iso7", 110).zero_count_reason(1, 1, "present", "low_callable_rate");
    let provenance_unverified =
        SiteRow::new("G1", "iso8", 110).zero_count_reason(1, 1, "present", "provenance_unverified");
    let future =
        SiteRow::new("G1", "iso9", 110).zero_count_reason(1, 0, "present", "future_qc_gate");
    let mut rows = vec![
        site_join,
        unknown_denominator,
        coverage_unavailable,
        reference_unvalidated,
        low_covering,
        low_candidate_rate,
        low_callable_rate,
        provenance_unverified,
        future,
    ];
    for row in &mut rows {
        row.coverage_basis = "unavailable";
        row.n_covering = None;
    }
    write_sites(&sites, &rows);

    let prefix = root.join("summary");
    let output = run_summary(&[sites.as_path()], &prefix);
    assert_success(&output, "open eligibility reasons");
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        "mod-site-summary: inputs=1 input_rows=9 sites=1\n"
    );

    let text = fs::read_to_string(summary_path(&prefix)).unwrap();
    let mut lines = text.lines();
    let columns = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let values = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    assert!(lines.next().is_none());
    let value = |column: &str| {
        let index = columns
            .iter()
            .position(|candidate| *candidate == column)
            .unwrap();
        values[index]
    };
    assert_eq!(value("n_isoforms_present"), "9");
    assert_eq!(value("n_isoforms_eligible"), "0");
    assert_eq!(value("n_isoforms_reference_base_mismatch"), "0");
    assert_eq!(value("n_isoforms_site_join_rate_low"), "1");
    assert_eq!(value("n_isoforms_unknown_denominator"), "1");
    assert_eq!(value("n_isoforms_coverage_unavailable"), "1");
    assert_eq!(value("n_isoforms_reference_unvalidated"), "1");
    assert_eq!(value("n_isoforms_low_covering"), "1");
    assert_eq!(value("n_isoforms_low_candidate_rate"), "1");
    assert_eq!(value("n_isoforms_low_callable_rate"), "1");
    assert_eq!(value("n_isoforms_provenance_unverified"), "1");
    assert_eq!(value("n_isoforms_other_ineligible"), "1");
    assert_eq!(value("summary_state"), "no_eligible_isoform");
}

#[test]
fn rejects_duplicate_and_unsorted_isoform_site_keys() {
    let root = TestDir::new("mod-site-summary-order");
    let left = root.join("left.tsv");
    let duplicate = root.join("duplicate.tsv");
    let row = SiteRow::new("G1", "iso1", 110).exact_eligible(2, 2, 1, 0, "1", "0", "1");
    write_sites(&left, std::slice::from_ref(&row));
    write_sites(&duplicate, std::slice::from_ref(&row));
    let output = run_summary(
        &[left.as_path(), duplicate.as_path()],
        &root.join("duplicate-out"),
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("duplicate isoform-site key across inputs"));

    let unsorted = root.join("unsorted.tsv");
    write_sites(
        &unsorted,
        &[
            SiteRow::new("G1", "iso1", 120),
            SiteRow::new("G1", "iso1", 110),
        ],
    );
    let output = run_summary(&[unsorted.as_path()], &root.join("unsorted-out"));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("is not sorted"));
}

#[test]
fn rejects_schema_count_and_cross_row_inconsistencies() {
    let root = TestDir::new("mod-site-summary-invalid");

    let bad_header = root.join("bad-header.tsv");
    fs::write(&bad_header, "assay_id\twrong\n").unwrap();
    let output = run_summary(&[bad_header.as_path()], &root.join("bad-header-out"));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("header mismatch"));

    let bad_counts = root.join("bad-counts.tsv");
    let mut count_row = SiteRow::new("G1", "iso1", 110);
    count_row.n_assigned = 1;
    count_row.n_covering = Some(1);
    count_row.n_candidate = 1;
    count_row.n_callable = 1;
    write_sites(&bad_counts, &[count_row]);
    let output = run_summary(&[bad_counts.as_path()], &root.join("bad-counts-out"));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("n_callable must equal n_modified + n_unmodified"));

    let thresholds = root.join("thresholds.tsv");
    let first = SiteRow::new("G1", "iso1", 110);
    let mut second = SiteRow::new("G1", "iso1", 120);
    second.analysis_threshold = "0.6";
    write_sites(&thresholds, &[first, second]);
    let output = run_summary(&[thresholds.as_path()], &root.join("thresholds-out"));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("inconsistent analysis thresholds"));

    let coverage = root.join("coverage.tsv");
    let exact = SiteRow::new("G1", "iso1", 110);
    let mut unavailable = SiteRow::new("G1", "iso1", 120);
    unavailable.coverage_basis = "unavailable";
    unavailable.n_covering = None;
    write_sites(&coverage, &[exact, unavailable]);
    let output = run_summary(&[coverage.as_path()], &root.join("coverage-out"));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("inconsistent coverage_basis"));
}

#[test]
fn conflicting_context_does_not_replace_existing_output() {
    let root = TestDir::new("mod-site-summary-atomic");
    let left = root.join("left.tsv");
    let right = root.join("right.tsv");
    write_sites(&left, &[SiteRow::new("G1", "iso1", 110)]);
    let mut conflicting = SiteRow::new("G1", "iso2", 110);
    conflicting.context = "RRACH";
    write_sites(&right, &[conflicting]);

    let prefix = root.join("result");
    let output_path = summary_path(&prefix);
    fs::write(&output_path, "sentinel\n").unwrap();
    let output = run_summary(&[left.as_path(), right.as_path()], &prefix);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("conflicting context"));
    assert_eq!(fs::read_to_string(output_path).unwrap(), "sentinel\n");
}

#[test]
fn help_exposes_repeatable_sites_and_output_prefix() {
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["mod-site-summary", "--help"])
        .output()
        .unwrap();
    assert_success(&output, "mod-site-summary help");
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Usage: trackcluster mod-site-summary"));
    assert!(stdout.contains("--sites <SITES>"));
    assert!(stdout.contains("--out <OUT>"));
}
