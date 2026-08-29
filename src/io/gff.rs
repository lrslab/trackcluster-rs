//! GFF3/GTF transcript annotation parsing and bigGenePred conversion.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;
use flate2::read::MultiGzDecoder;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Open a plain-text or gzip-compressed input as a buffered reader.
///
/// Compression is detected from the gzip magic bytes so compressed inputs do
/// not need a `.gz` suffix. A `.gz` suffix without gzip content is rejected.
pub(crate) fn open_maybe_gzip(path: &Path) -> anyhow::Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("open text input {path:?}"))?;
    let mut reader = BufReader::new(file);
    let is_gzip = reader
        .fill_buf()
        .with_context(|| format!("inspect text input {path:?}"))?
        .starts_with(&GZIP_MAGIC);
    let has_gzip_suffix = path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"));
    if has_gzip_suffix && !is_gzip {
        anyhow::bail!("input {path:?} has a .gz suffix but is not gzip-compressed");
    }
    if is_gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

/// Annotation attribute syntax accepted by [`read_annotation_transcripts`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum AnnotationFormat {
    /// Detect GFF3 (`key=value`) or GTF (`key "value"`) attributes per row.
    #[default]
    Auto,
    /// Require GFF3 attributes and `ID`/`Parent` relationships.
    Gff3,
    /// Require GTF attributes such as `gene_id` and `transcript_id`.
    Gtf,
}

/// Options controlling annotation-to-bigGenePred conversion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GffToBiggOptions {
    /// Input attribute syntax.
    pub format: AnnotationFormat,
    /// GFF3 gene-feature attribute used as the output gene label.
    ///
    /// Graph relationships always use canonical `ID`/`Parent`; this option
    /// changes only the gene label written to BED.
    pub gene_key: String,
}

impl Default for GffToBiggOptions {
    fn default() -> Self {
        Self {
            format: AnnotationFormat::Auto,
            gene_key: "ID".to_owned(),
        }
    }
}

/// Policy for annotation-record failures in general-purpose GFF/GTF conversion.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum InvalidAnnotationPolicy {
    /// Quarantine an identifiable invalid transcript and continue with other models.
    #[default]
    Skip,
    /// Stop at the first invalid annotation record or transcript model.
    Fail,
}

/// Auditable description of an annotation record or transcript model excluded
/// by [`InvalidAnnotationPolicy::Skip`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RejectedAnnotationRecord {
    /// Source annotation path.
    pub path: PathBuf,
    /// One-based triggering source line. For a multi-row `model` rejection,
    /// this is the earliest row in that model and therefore an anchor line;
    /// the reason retains exact participating line numbers when available.
    pub line: usize,
    /// Transcript IDs quarantined by this rejection, when identifiable.
    pub transcript_ids: Vec<String>,
    /// Stable rejection class: `ignored_feature`, `parse`, `gene`, or `model`.
    pub kind: &'static str,
    /// Human-readable reason for exclusion.
    pub reason: String,
}

/// Result of policy-aware annotation conversion.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotationReadResult {
    /// Valid transcript models assembled from the input.
    pub transcripts: Vec<Transcript>,
    /// Source records and transcript models excluded in recovering mode.
    pub rejected_records: Vec<RejectedAnnotationRecord>,
    /// Number of exon-bearing transcript models quarantined in recovering mode.
    pub rejected_transcripts: usize,
}

#[derive(Clone, Debug)]
struct FeatureRow {
    line: usize,
    seqid: String,
    feature_type: String,
    interval: Interval,
    strand: Strand,
    format: AnnotationFormat,
    attributes: BTreeMap<String, String>,
    parents: Vec<String>,
}

#[derive(Clone, Debug)]
struct TranscriptAssembly {
    seqid: String,
    strand: Strand,
    exons: Vec<Interval>,
    gene_hints: BTreeSet<String>,
}

#[derive(Clone, Debug)]
struct TranscriptFeature {
    line: usize,
    seqid: String,
    interval: Interval,
    strand: Strand,
    gene_hints: BTreeSet<String>,
}

#[derive(Clone, Debug)]
struct ParsedAttributes {
    format: AnnotationFormat,
    values: BTreeMap<String, String>,
    parents: Vec<String>,
}

fn parse_u32(field: &str, value: &str) -> anyhow::Result<u32> {
    value
        .parse::<u32>()
        .with_context(|| format!("invalid integer for {field}: {value:?}"))
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

fn percent_decode(value: &str) -> anyhow::Result<String> {
    let input = value.as_bytes();
    let mut decoded = Vec::with_capacity(input.len());
    let mut index = 0usize;
    while index < input.len() {
        if input[index] != b'%' {
            decoded.push(input[index]);
            index += 1;
            continue;
        }
        let high = input
            .get(index + 1)
            .copied()
            .and_then(hex_value)
            .with_context(|| {
                format!("invalid percent encoding in attribute {value:?} at byte {index}")
            })?;
        let low = input
            .get(index + 2)
            .copied()
            .and_then(hex_value)
            .with_context(|| {
                format!("invalid percent encoding in attribute {value:?} at byte {index}")
            })?;
        decoded.push((high << 4) | low);
        index += 3;
    }
    String::from_utf8(decoded)
        .with_context(|| format!("percent-decoded attribute is not UTF-8: {value:?}"))
}

fn insert_attribute(
    attributes: &mut BTreeMap<String, String>,
    key: String,
    value: String,
) -> anyhow::Result<()> {
    if attributes.insert(key.clone(), value).is_some() {
        anyhow::bail!("duplicate annotation attribute {key:?}");
    }
    Ok(())
}

fn insert_gtf_attribute(
    attributes: &mut BTreeMap<String, String>,
    key: String,
    value: String,
) -> anyhow::Result<()> {
    use std::collections::btree_map::Entry;

    match attributes.entry(key) {
        Entry::Vacant(entry) => {
            entry.insert(value);
        }
        Entry::Occupied(mut entry)
            if !matches!(entry.key().as_str(), "gene_id" | "transcript_id") =>
        {
            entry.get_mut().push(',');
            entry.get_mut().push_str(&value);
        }
        Entry::Occupied(entry) => {
            anyhow::bail!("duplicate annotation identity attribute {:?}", entry.key());
        }
    }
    Ok(())
}

fn parse_gff3_attributes(raw: &str) -> anyhow::Result<ParsedAttributes> {
    let mut attributes = BTreeMap::new();
    let mut parents = Vec::new();
    if raw == "." {
        return Ok(ParsedAttributes {
            format: AnnotationFormat::Gff3,
            values: attributes,
            parents,
        });
    }
    for field in raw
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
    {
        let (key, value) = field
            .split_once('=')
            .with_context(|| format!("expected GFF3 key=value attribute, got {field:?}"))?;
        let key = percent_decode(key.trim())?;
        let value = value.trim();
        if key == "Parent" {
            for parent in value.split(',') {
                let parent = parent.trim();
                if parent.is_empty() {
                    anyhow::bail!("GFF3 Parent contains an empty identity in {field:?}");
                }
                parents.push(percent_decode(parent)?);
            }
        }
        insert_attribute(&mut attributes, key, percent_decode(value)?)?;
    }
    Ok(ParsedAttributes {
        format: AnnotationFormat::Gff3,
        values: attributes,
        parents,
    })
}

fn parse_gtf_attributes(raw: &str) -> anyhow::Result<ParsedAttributes> {
    let mut attributes = BTreeMap::new();
    if raw == "." {
        return Ok(ParsedAttributes {
            format: AnnotationFormat::Gtf,
            values: attributes,
            parents: Vec::new(),
        });
    }

    let bytes = raw.as_bytes();
    let mut index = 0usize;
    while index < bytes.len() {
        while index < bytes.len() && (bytes[index].is_ascii_whitespace() || bytes[index] == b';') {
            index += 1;
        }
        if index == bytes.len() {
            break;
        }

        let key_start = index;
        while index < bytes.len() && !bytes[index].is_ascii_whitespace() && bytes[index] != b';' {
            index += 1;
        }
        if key_start == index {
            anyhow::bail!("invalid GTF attribute near byte {index}: {raw:?}");
        }
        let key = &raw[key_start..index];
        if index == bytes.len() || !bytes[index].is_ascii_whitespace() {
            anyhow::bail!("expected whitespace after GTF attribute key {key:?}");
        }
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        if index == bytes.len() {
            anyhow::bail!("missing value for GTF attribute {key:?}");
        }

        let value = if bytes[index] == b'"' {
            index += 1;
            let mut value = Vec::new();
            let mut terminated = false;
            while index < bytes.len() {
                match bytes[index] {
                    b'"' => {
                        index += 1;
                        terminated = true;
                        break;
                    }
                    b'\\' => {
                        index += 1;
                        let escaped = bytes.get(index).copied().with_context(|| {
                            format!("truncated escape in GTF attribute {key:?}")
                        })?;
                        if matches!(escaped, b'"' | b'\\') {
                            value.push(escaped);
                        } else {
                            value.push(b'\\');
                            value.push(escaped);
                        }
                        index += 1;
                    }
                    byte => {
                        value.push(byte);
                        index += 1;
                    }
                }
            }
            if !terminated {
                anyhow::bail!("unterminated quoted GTF attribute {key:?}");
            }
            String::from_utf8(value)
                .with_context(|| format!("GTF attribute {key:?} is not valid UTF-8"))?
        } else {
            let value_start = index;
            while index < bytes.len() && bytes[index] != b';' {
                index += 1;
            }
            raw[value_start..index].trim().to_owned()
        };

        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        if index < bytes.len() && bytes[index] != b';' {
            anyhow::bail!("unexpected text after GTF attribute {key:?}");
        }
        if index < bytes.len() {
            index += 1;
        }
        insert_gtf_attribute(&mut attributes, key.to_owned(), value)?;
    }

    Ok(ParsedAttributes {
        format: AnnotationFormat::Gtf,
        values: attributes,
        parents: Vec::new(),
    })
}

fn parse_attributes(
    raw: &str,
    requested_format: AnnotationFormat,
) -> anyhow::Result<ParsedAttributes> {
    let format = match requested_format {
        AnnotationFormat::Auto
            if raw
                .split(';')
                .next()
                .and_then(|field| field.find('=').map(|index| (&field[..index], index)))
                .is_some_and(|(prefix, _)| !prefix.chars().any(char::is_whitespace)) =>
        {
            AnnotationFormat::Gff3
        }
        AnnotationFormat::Auto => AnnotationFormat::Gtf,
        format => format,
    };
    match format {
        AnnotationFormat::Auto => unreachable!("auto format is resolved above"),
        AnnotationFormat::Gff3 => parse_gff3_attributes(raw),
        AnnotationFormat::Gtf => parse_gtf_attributes(raw),
    }
}

fn parse_strand(value: &str) -> anyhow::Result<Strand> {
    match value {
        "?" => Ok(Strand::Unknown),
        value => Strand::try_from(value).map_err(Into::into),
    }
}

fn parse_feature_row(
    line: &str,
    line_number: usize,
    format: AnnotationFormat,
) -> anyhow::Result<FeatureRow> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 9 {
        anyhow::bail!(
            "expected exactly 9 tab-separated columns, got {}",
            fields.len()
        );
    }
    if fields[0].is_empty() || fields[0] == "." {
        anyhow::bail!("annotation seqid must not be empty or '.'");
    }
    if fields[0].chars().any(char::is_control) {
        anyhow::bail!("annotation seqid contains a control character");
    }
    if fields[2].is_empty() || fields[2] == "." {
        anyhow::bail!("annotation feature type must not be empty or '.'");
    }
    let start = parse_u32("start", fields[3])?;
    let end = parse_u32("end", fields[4])?;
    if start == 0 {
        anyhow::bail!("annotation start must be at least 1");
    }
    if end < start {
        anyhow::bail!("annotation end {end} is before start {start}");
    }
    let parsed_attributes = parse_attributes(fields[8], format)?;
    Ok(FeatureRow {
        line: line_number,
        seqid: fields[0].to_owned(),
        feature_type: fields[2].to_owned(),
        interval: Interval::new(Coord::new(start - 1), Coord::new(end))?,
        strand: parse_strand(fields[6])?,
        format: parsed_attributes.format,
        attributes: parsed_attributes.values,
        parents: parsed_attributes.parents,
    })
}

fn is_transcript_feature(feature_type: &str) -> bool {
    let feature_type = feature_type.to_ascii_lowercase();
    matches!(
        feature_type.as_str(),
        "mrna"
            | "transcript"
            | "rna"
            | "lnc_rna"
            | "ncrna"
            | "mirna"
            | "rrna"
            | "trna"
            | "snrna"
            | "snorna"
            | "circrna"
            | "pseudogenic_transcript"
    ) || feature_type.ends_with("_rna")
        || feature_type.ends_with("rna")
}

fn is_gene_feature(feature_type: &str) -> bool {
    feature_type.to_ascii_lowercase().ends_with("gene")
}

fn transcript_id(row: &FeatureRow) -> Option<&str> {
    match row.format {
        AnnotationFormat::Gff3 => row.attributes.get("ID"),
        AnnotationFormat::Gtf => row.attributes.get("transcript_id"),
        AnnotationFormat::Auto => unreachable!("feature rows always have a resolved format"),
    }
    .map(String::as_str)
    .filter(|value| !value.trim().is_empty())
}

fn validate_bed_identity(kind: &str, value: &str, line: usize) -> anyhow::Result<()> {
    if value.trim().is_empty() || value == "." {
        anyhow::bail!("annotation line {line} has an empty {kind}");
    }
    if value.chars().any(char::is_control) {
        anyhow::bail!("annotation line {line} {kind} {value:?} contains a control character");
    }
    Ok(())
}

fn add_gene_hints(target: &mut BTreeSet<String>, values: impl IntoIterator<Item = String>) {
    target.extend(
        values
            .into_iter()
            .map(|value| value.trim().to_owned())
            .filter(|value| !value.is_empty() && value != "none"),
    );
}

fn build_transcripts(
    rows: &[FeatureRow],
    options: &GffToBiggOptions,
) -> anyhow::Result<Vec<Transcript>> {
    let gene_key = options.gene_key.trim();
    if gene_key.is_empty() {
        anyhow::bail!("gene attribute key must not be empty");
    }

    let gff_transcript_ids: BTreeSet<String> = rows
        .iter()
        .filter(|row| {
            row.format == AnnotationFormat::Gff3 && row.feature_type.eq_ignore_ascii_case("exon")
        })
        .flat_map(|row| row.parents.iter().cloned())
        .collect();
    let mut gene_labels: HashMap<String, String> = HashMap::new();
    let mut transcript_features: HashMap<String, TranscriptFeature> = HashMap::new();

    for row in rows {
        if is_gene_feature(&row.feature_type) {
            let id = match row.format {
                AnnotationFormat::Gff3 => row.attributes.get("ID"),
                AnnotationFormat::Gtf => row.attributes.get("gene_id"),
                AnnotationFormat::Auto => {
                    unreachable!("feature rows always have a resolved format")
                }
            };
            if let Some(id) = id {
                validate_bed_identity("gene identity", id, row.line)?;
                let label = row
                    .attributes
                    .get(gene_key)
                    .or_else(|| row.attributes.get("ID"))
                    .or_else(|| row.attributes.get("gene_id"))
                    .cloned()
                    .unwrap_or_else(|| id.clone());
                validate_bed_identity("gene label", &label, row.line)?;
                if gene_labels.insert(id.clone(), label).is_some() {
                    anyhow::bail!(
                        "duplicate gene identity {id:?} at annotation line {}",
                        row.line
                    );
                }
            }
            continue;
        }

        let id = match row.format {
            AnnotationFormat::Gff3 if !row.feature_type.eq_ignore_ascii_case("exon") => row
                .attributes
                .get("ID")
                .filter(|id| gff_transcript_ids.contains(*id))
                .map(String::as_str),
            AnnotationFormat::Gtf if is_transcript_feature(&row.feature_type) => transcript_id(row),
            AnnotationFormat::Auto => unreachable!("feature rows always have a resolved format"),
            _ => None,
        };
        if let Some(id) = id {
            validate_bed_identity("transcript identity", id, row.line)?;
            let mut gene_hints = BTreeSet::new();
            add_gene_hints(&mut gene_hints, row.parents.iter().cloned());
            if let Some(gene_id) = row.attributes.get("gene_id") {
                add_gene_hints(&mut gene_hints, std::iter::once(gene_id.clone()));
            }
            let feature = TranscriptFeature {
                line: row.line,
                seqid: row.seqid.clone(),
                interval: row.interval,
                strand: row.strand,
                gene_hints,
            };
            if let Some(previous) = transcript_features.insert(id.to_owned(), feature) {
                anyhow::bail!(
                    "duplicate transcript identity {id:?} at annotation lines {} and {}",
                    previous.line,
                    row.line
                );
            }
        }
    }

    let mut assemblies: BTreeMap<String, TranscriptAssembly> = BTreeMap::new();
    for row in rows
        .iter()
        .filter(|row| row.feature_type.eq_ignore_ascii_case("exon"))
    {
        let parents = match row.format {
            AnnotationFormat::Gff3 => row.parents.clone(),
            AnnotationFormat::Gtf => row
                .attributes
                .get("transcript_id")
                .cloned()
                .into_iter()
                .collect(),
            AnnotationFormat::Auto => unreachable!("feature rows always have a resolved format"),
        };
        if parents.is_empty() {
            anyhow::bail!(
                "annotation line {} exon has no Parent or transcript_id",
                row.line
            );
        }

        for parent in parents {
            validate_bed_identity("exon parent identity", &parent, row.line)?;
            if row.format == AnnotationFormat::Gff3 && !transcript_features.contains_key(&parent) {
                if gene_labels.contains_key(&parent) {
                    anyhow::bail!(
                        "annotation line {} exon Parent {parent:?} refers to a gene, not a transcript",
                        row.line
                    );
                }
                anyhow::bail!(
                    "annotation line {} exon Parent {parent:?} has no declared transcript feature",
                    row.line
                );
            }
            let entry = assemblies
                .entry(parent.clone())
                .or_insert_with(|| TranscriptAssembly {
                    seqid: row.seqid.clone(),
                    strand: row.strand,
                    exons: Vec::new(),
                    gene_hints: BTreeSet::new(),
                });
            if entry.seqid != row.seqid {
                anyhow::bail!(
                    "transcript {parent:?} has exons on multiple sequences: {:?} and {:?}",
                    entry.seqid,
                    row.seqid
                );
            }
            if entry.strand == Strand::Unknown {
                entry.strand = row.strand;
            } else if row.strand != Strand::Unknown && entry.strand != row.strand {
                anyhow::bail!(
                    "transcript {parent:?} has conflicting exon strands {} and {}",
                    entry.strand.as_char(),
                    row.strand.as_char()
                );
            }
            entry.exons.push(row.interval);
            if let Some(gene_id) = row.attributes.get("gene_id") {
                add_gene_hints(&mut entry.gene_hints, std::iter::once(gene_id.clone()));
            }
        }
    }

    if assemblies.is_empty() {
        anyhow::bail!("annotation contains no exon features with transcript identities");
    }

    let mut transcripts = Vec::with_capacity(assemblies.len());
    for (id, mut assembly) in assemblies {
        assembly.exons.sort_unstable();
        assembly.exons.dedup();
        let tx_start = assembly.exons[0].start;
        let tx_end = assembly.exons[assembly.exons.len() - 1].end;

        let mut raw_genes = BTreeSet::new();
        if let Some(feature) = transcript_features.remove(&id) {
            if feature.seqid != assembly.seqid {
                anyhow::bail!(
                    "transcript {id:?} feature at line {} is on {:?}, but its exons are on {:?}",
                    feature.line,
                    feature.seqid,
                    assembly.seqid
                );
            }
            if feature.interval.start > tx_start || feature.interval.end < tx_end {
                anyhow::bail!(
                    "transcript {id:?} feature span [{}, {}) does not contain exon span [{}, {})",
                    feature.interval.start.get(),
                    feature.interval.end.get(),
                    tx_start.get(),
                    tx_end.get()
                );
            }
            if assembly.strand == Strand::Unknown {
                assembly.strand = feature.strand;
            } else if feature.strand != Strand::Unknown && assembly.strand != feature.strand {
                anyhow::bail!(
                    "transcript {id:?} feature strand {} conflicts with exon strand {}",
                    feature.strand.as_char(),
                    assembly.strand.as_char()
                );
            }
            raw_genes.extend(feature.gene_hints);
        }
        raw_genes.extend(assembly.gene_hints);
        let mut genes = BTreeSet::new();
        for raw_gene in raw_genes {
            genes.insert(gene_labels.get(&raw_gene).cloned().unwrap_or(raw_gene));
        }
        let gene = if genes.is_empty() {
            "none".to_owned()
        } else {
            genes.into_iter().collect::<Vec<_>>().join("||")
        };
        if gene != "none" {
            for value in gene.split("||") {
                crate::flow::path_key::GeneId::parse(value).with_context(|| {
                    format!("transcript {id:?} has an unsafe output gene ID {value:?}")
                })?;
            }
        }

        let exon_frames = format!(
            "{},",
            std::iter::repeat_n("-1", assembly.exons.len())
                .collect::<Vec<_>>()
                .join(",")
        );
        let transcript = Transcript::new(
            assembly.seqid,
            assembly.strand,
            tx_start,
            tx_end,
            id.clone(),
            assembly.exons,
            Bed12Attrs {
                score: 100,
                thick_start: Coord::new(0),
                thick_end: Coord::new(0),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    exon_frames,
                    "isoform_anno".to_owned(),
                    gene,
                    "none".to_owned(),
                    "none".to_owned(),
                ],
            },
        )
        .with_context(|| format!("build transcript {id:?} from annotation exons"))?;
        transcripts.push(transcript);
    }

    crate::identity::validate_reference_ids(&transcripts)
        .context("validate GFF/GTF transcript identities")?;
    transcripts.sort_by(crate::identity::transcript_order);
    Ok(transcripts)
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum ParseRecoveryScope {
    IgnoredFeature(String),
    Gene(String),
    Transcripts { ids: Vec<String>, is_exon: bool },
}

fn row_transcript_ids(row: &FeatureRow) -> Vec<String> {
    if row.feature_type.eq_ignore_ascii_case("exon") {
        return match row.format {
            AnnotationFormat::Gff3 => row.parents.clone(),
            AnnotationFormat::Gtf => row
                .attributes
                .get("transcript_id")
                .cloned()
                .into_iter()
                .collect(),
            AnnotationFormat::Auto => {
                unreachable!("feature rows always have a resolved format")
            }
        };
    }
    if is_transcript_feature(&row.feature_type) {
        return transcript_id(row).map(str::to_owned).into_iter().collect();
    }
    Vec::new()
}

fn recover_parse_scope(
    line: &str,
    requested_format: AnnotationFormat,
) -> anyhow::Result<ParseRecoveryScope> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 9 {
        anyhow::bail!(
            "record has {} columns, so its feature type and model ownership cannot be trusted",
            fields.len()
        );
    }
    let feature_type = fields[2];
    if feature_type.is_empty() || feature_type == "." || feature_type.chars().any(char::is_control)
    {
        anyhow::bail!("record has an invalid feature type {feature_type:?}");
    }
    if !feature_type.eq_ignore_ascii_case("exon")
        && !is_transcript_feature(feature_type)
        && !is_gene_feature(feature_type)
    {
        return Ok(ParseRecoveryScope::IgnoredFeature(feature_type.to_owned()));
    }

    let parsed = parse_attributes(fields[8], requested_format)
        .context("cannot recover identity from malformed record attributes")?;
    let identity = |key: &str| {
        parsed
            .values
            .get(key)
            .map(String::as_str)
            .map(str::trim)
            .filter(|value| !value.is_empty() && *value != ".")
            .map(str::to_owned)
    };
    if is_gene_feature(feature_type) {
        let id = match parsed.format {
            AnnotationFormat::Gff3 => identity("ID"),
            AnnotationFormat::Gtf => identity("gene_id"),
            AnnotationFormat::Auto => unreachable!("attribute format is resolved"),
        }
        .context("malformed gene record has no recoverable gene identity")?;
        return Ok(ParseRecoveryScope::Gene(id));
    }

    let ids = if feature_type.eq_ignore_ascii_case("exon") {
        match parsed.format {
            AnnotationFormat::Gff3 => parsed.parents,
            AnnotationFormat::Gtf => identity("transcript_id").into_iter().collect(),
            AnnotationFormat::Auto => unreachable!("attribute format is resolved"),
        }
    } else {
        match parsed.format {
            AnnotationFormat::Gff3 => identity("ID").into_iter().collect(),
            AnnotationFormat::Gtf => identity("transcript_id").into_iter().collect(),
            AnnotationFormat::Auto => unreachable!("attribute format is resolved"),
        }
    };
    if ids.is_empty() {
        anyhow::bail!("malformed {feature_type:?} record has no recoverable transcript identity");
    }
    Ok(ParseRecoveryScope::Transcripts {
        ids,
        is_exon: feature_type.eq_ignore_ascii_case("exon"),
    })
}

fn gene_identity(row: &FeatureRow) -> Option<&str> {
    match row.format {
        AnnotationFormat::Gff3 => row.attributes.get("ID"),
        AnnotationFormat::Gtf => row.attributes.get("gene_id"),
        AnnotationFormat::Auto => unreachable!("feature rows always have a resolved format"),
    }
    .map(String::as_str)
}

fn transcript_gene_hints(rows: &[FeatureRow]) -> BTreeSet<String> {
    let mut hints = BTreeSet::new();
    for row in rows {
        if is_transcript_feature(&row.feature_type) {
            add_gene_hints(&mut hints, row.parents.iter().cloned());
        }
        if let Some(gene_id) = row.attributes.get("gene_id") {
            add_gene_hints(&mut hints, std::iter::once(gene_id.clone()));
        }
    }
    hints
}

fn build_transcripts_recovering(
    rows: &[FeatureRow],
    options: &GffToBiggOptions,
    path: &Path,
    poisoned_transcripts: &BTreeSet<String>,
    poisoned_exon_transcripts: &BTreeSet<String>,
    poisoned_genes: &BTreeMap<String, usize>,
    rejected_records: &mut Vec<RejectedAnnotationRecord>,
) -> (Vec<Transcript>, usize) {
    let candidate_ids: BTreeSet<String> = rows
        .iter()
        .filter(|row| row.feature_type.eq_ignore_ascii_case("exon"))
        .flat_map(row_transcript_ids)
        .collect();
    let mut groups: BTreeMap<String, Vec<FeatureRow>> = candidate_ids
        .iter()
        .cloned()
        .map(|id| (id, Vec::new()))
        .collect();
    let mut gene_rows: BTreeMap<String, Vec<FeatureRow>> = BTreeMap::new();

    for row in rows {
        if is_gene_feature(&row.feature_type) {
            if let Some(id) = gene_identity(row) {
                gene_rows
                    .entry(id.to_owned())
                    .or_default()
                    .push(row.clone());
            }
            continue;
        }
        if row.feature_type.eq_ignore_ascii_case("exon") {
            for id in row_transcript_ids(row) {
                let Some(group) = groups.get_mut(&id) else {
                    continue;
                };
                let mut owned = row.clone();
                if owned.format == AnnotationFormat::Gff3 {
                    owned.parents = vec![id];
                }
                group.push(owned);
            }
            continue;
        }
        if is_transcript_feature(&row.feature_type) {
            if let Some(id) = transcript_id(row) {
                if let Some(group) = groups.get_mut(id) {
                    group.push(row.clone());
                }
            }
        }
    }

    let mut transcripts = Vec::new();
    let mut rejected_transcript_ids = poisoned_exon_transcripts.clone();
    for (id, mut group) in groups {
        if poisoned_transcripts.contains(&id) {
            rejected_transcript_ids.insert(id);
            continue;
        }
        let line = group.iter().map(|row| row.line).min().unwrap_or(1);
        let gene_hints = transcript_gene_hints(&group);
        let invalid_gene = gene_hints
            .iter()
            .find_map(|gene| poisoned_genes.get(gene).map(|bad_line| (gene, *bad_line)));
        if let Some((gene, bad_line)) = invalid_gene {
            rejected_records.push(RejectedAnnotationRecord {
                path: path.to_path_buf(),
                line: bad_line,
                transcript_ids: vec![id.clone()],
                kind: "model",
                reason: format!(
                    "transcript depends on malformed gene {gene:?} from annotation line {bad_line}"
                ),
            });
            rejected_transcript_ids.insert(id);
            continue;
        }
        for gene in gene_hints {
            if let Some(features) = gene_rows.get(&gene) {
                group.extend(features.iter().cloned());
            }
        }
        match build_transcripts(&group, options) {
            Ok(mut built) if built.len() == 1 => transcripts.append(&mut built),
            Ok(built) => {
                rejected_records.push(RejectedAnnotationRecord {
                    path: path.to_path_buf(),
                    line,
                    transcript_ids: vec![id.clone()],
                    kind: "model",
                    reason: format!(
                        "isolated transcript assembly produced {} models instead of one",
                        built.len()
                    ),
                });
                rejected_transcript_ids.insert(id);
            }
            Err(error) => {
                rejected_records.push(RejectedAnnotationRecord {
                    path: path.to_path_buf(),
                    line,
                    transcript_ids: vec![id.clone()],
                    kind: "model",
                    reason: format!("{error:#}"),
                });
                rejected_transcript_ids.insert(id);
            }
        }
    }
    transcripts.sort_by(crate::identity::transcript_order);
    (transcripts, rejected_transcript_ids.len())
}

fn read_annotation_transcripts_recovering(
    path: &Path,
    options: &GffToBiggOptions,
) -> anyhow::Result<AnnotationReadResult> {
    let reader =
        open_maybe_gzip(path).with_context(|| format!("open annotation input {path:?}"))?;
    let mut rows = Vec::new();
    let mut rejected_records = Vec::new();
    let mut poisoned_transcripts = BTreeSet::new();
    let mut poisoned_exon_transcripts = BTreeSet::new();
    let mut poisoned_genes = BTreeMap::new();
    for (line_index, result) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = result.with_context(|| format!("read annotation {path:?}:{line_number}"))?;
        let line = line.trim();
        if line == "##FASTA" {
            break;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let row = match parse_feature_row(line, line_number, options.format) {
            Ok(row) => row,
            Err(parse_error) => {
                let scope = recover_parse_scope(line, options.format).map_err(|scope_error| {
                    anyhow::anyhow!(
                        "parse annotation {path:?}:{line_number}: {parse_error:#}; refusing unsafe recovery: {scope_error:#}"
                    )
                })?;
                match scope {
                    ParseRecoveryScope::IgnoredFeature(feature_type) => {
                        rejected_records.push(RejectedAnnotationRecord {
                            path: path.to_path_buf(),
                            line: line_number,
                            transcript_ids: Vec::new(),
                            kind: "ignored_feature",
                            reason: format!(
                                "ignored malformed non-model feature {feature_type:?}: {parse_error:#}"
                            ),
                        });
                    }
                    ParseRecoveryScope::Gene(id) => {
                        poisoned_genes.insert(id.clone(), line_number);
                        rejected_records.push(RejectedAnnotationRecord {
                            path: path.to_path_buf(),
                            line: line_number,
                            transcript_ids: Vec::new(),
                            kind: "gene",
                            reason: format!("malformed gene {id:?}: {parse_error:#}"),
                        });
                    }
                    ParseRecoveryScope::Transcripts { mut ids, is_exon } => {
                        ids.sort();
                        ids.dedup();
                        poisoned_transcripts.extend(ids.iter().cloned());
                        if is_exon {
                            poisoned_exon_transcripts.extend(ids.iter().cloned());
                        }
                        rejected_records.push(RejectedAnnotationRecord {
                            path: path.to_path_buf(),
                            line: line_number,
                            transcript_ids: ids,
                            kind: "parse",
                            reason: format!("{parse_error:#}"),
                        });
                    }
                }
                continue;
            }
        };
        if row.feature_type.eq_ignore_ascii_case("exon") && row_transcript_ids(&row).is_empty() {
            anyhow::bail!(
                "parse annotation {path:?}:{line_number}: exon has no Parent or transcript_id; refusing unsafe recovery because the affected transcript is unknown"
            );
        }
        if is_gene_feature(&row.feature_type)
            || is_transcript_feature(&row.feature_type)
            || row.feature_type.eq_ignore_ascii_case("exon")
        {
            rows.push(row);
        }
    }

    let (transcripts, rejected_transcripts) = build_transcripts_recovering(
        &rows,
        options,
        path,
        &poisoned_transcripts,
        &poisoned_exon_transcripts,
        &poisoned_genes,
        &mut rejected_records,
    );
    rejected_records.sort_by(|left, right| {
        left.line
            .cmp(&right.line)
            .then_with(|| left.kind.cmp(right.kind))
            .then_with(|| left.transcript_ids.cmp(&right.transcript_ids))
    });
    Ok(AnnotationReadResult {
        transcripts,
        rejected_records,
        rejected_transcripts,
    })
}

/// Read GFF3 or GTF exon annotations and build validated BED12+8 transcripts.
///
/// GFF3 relationships use `ID`/`Parent`; GTF relationships use
/// `transcript_id`/`gene_id`. Input order does not affect output order. This
/// compatibility API is strict; use [`read_annotation_transcripts_with_policy`]
/// for explicitly auditable record recovery.
pub fn read_annotation_transcripts(
    path: &Path,
    options: &GffToBiggOptions,
) -> anyhow::Result<Vec<Transcript>> {
    let transcripts = read_annotation_transcripts_where(path, options, true, |_| true)?;
    if transcripts.is_empty() {
        anyhow::bail!("annotation contains no exon features with transcript identities");
    }
    Ok(transcripts)
}

/// Read an annotation under an explicit invalid-record policy.
///
/// Recovering mode only skips malformed non-model features or quarantines an
/// entire transcript whose identity is known. A malformed exon/transcript row
/// with unknown ownership, decompression error, or line-I/O error remains
/// fatal so recovery cannot silently publish a truncated transcript model.
pub fn read_annotation_transcripts_with_policy(
    path: &Path,
    options: &GffToBiggOptions,
    policy: InvalidAnnotationPolicy,
) -> anyhow::Result<AnnotationReadResult> {
    let result = match policy {
        InvalidAnnotationPolicy::Skip => read_annotation_transcripts_recovering(path, options)?,
        InvalidAnnotationPolicy::Fail => AnnotationReadResult {
            transcripts: read_annotation_transcripts(path, options)?,
            rejected_records: Vec::new(),
            rejected_transcripts: 0,
        },
    };
    if result.transcripts.is_empty() {
        if result.rejected_transcripts == 0 {
            anyhow::bail!("annotation contains no exon features with transcript identities");
        }
        anyhow::bail!(
            "annotation contains no valid transcript models; rejected_records={} rejected_transcripts={}",
            result.rejected_records.len(),
            result.rejected_transcripts
        );
    }
    Ok(result)
}

/// Header used by rejected-annotation audit files.
pub const REJECTED_ANNOTATIONS_TSV_HEADER: &str =
    "source_path\tanchor_line\ttranscript_ids_json\tkind\treason";

fn escape_rejected_annotation_field(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

/// Write rejected annotation diagnostics as stable, escaped TSV.
pub fn write_rejected_annotations_tsv_to_writer<W: Write>(
    writer: &mut W,
    records: &[RejectedAnnotationRecord],
) -> std::io::Result<()> {
    writeln!(writer, "{REJECTED_ANNOTATIONS_TSV_HEADER}")?;
    for record in records {
        let transcript_ids = serde_json::to_string(&record.transcript_ids)
            .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidData, error))?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            escape_rejected_annotation_field(&record.path.to_string_lossy()),
            record.line,
            transcript_ids,
            record.kind,
            escape_rejected_annotation_field(&record.reason)
        )?;
    }
    writer.flush()
}

/// Read only annotation rows needed to assemble transcripts accepted by
/// `keep_transcript`.
///
/// `retain_gene_features` preserves configured GFF3 gene labels for general
/// conversion. Coordinate-only callers can disable it: transcript and exon
/// gene IDs remain available, while unrelated genome-wide gene rows are not
/// materialized. Irrelevant feature kinds are validated but released
/// immediately instead of being retained for the whole annotation.
pub(crate) fn read_annotation_transcripts_where<F>(
    path: &Path,
    options: &GffToBiggOptions,
    retain_gene_features: bool,
    keep_transcript: F,
) -> anyhow::Result<Vec<Transcript>>
where
    F: Fn(&str) -> bool,
{
    let reader =
        open_maybe_gzip(path).with_context(|| format!("open annotation input {path:?}"))?;
    let mut rows = Vec::new();
    for (line_index, result) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = result.with_context(|| format!("read annotation {path:?}:{line_number}"))?;
        let line = line.trim();
        if line == "##FASTA" {
            break;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut row = parse_feature_row(line, line_number, options.format)
            .with_context(|| format!("parse annotation {path:?}:{line_number}"))?;
        let retain = if is_gene_feature(&row.feature_type) {
            retain_gene_features
        } else if is_transcript_feature(&row.feature_type) {
            transcript_id(&row).is_some_and(&keep_transcript)
        } else if row.feature_type.eq_ignore_ascii_case("exon") {
            match row.format {
                AnnotationFormat::Gff3 => {
                    row.parents.retain(|parent| keep_transcript(parent));
                    !row.parents.is_empty()
                }
                AnnotationFormat::Gtf => row
                    .attributes
                    .get("transcript_id")
                    .is_some_and(|id| keep_transcript(id)),
                AnnotationFormat::Auto => {
                    unreachable!("feature rows always have a resolved format")
                }
            }
        } else {
            false
        };
        if retain {
            rows.push(row);
        }
    }
    if !rows
        .iter()
        .any(|row| row.feature_type.eq_ignore_ascii_case("exon"))
    {
        return Ok(Vec::new());
    }
    build_transcripts(&rows, options).with_context(|| format!("convert annotation {path:?}"))
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::*;

    fn temp_dir(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-gff-{label}-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&path).unwrap();
        path
    }

    fn parse_rows(text: &str, format: AnnotationFormat) -> Vec<FeatureRow> {
        text.lines()
            .enumerate()
            .map(|(index, line)| parse_feature_row(line, index + 1, format).unwrap())
            .collect()
    }

    #[test]
    fn parses_gff3_percent_encoding_and_gtf_attributes() {
        let gff = parse_gff3_attributes("ID=tx%2D1;Parent=tx%2C1,gene%3A1").unwrap();
        assert_eq!(gff.values.get("ID").map(String::as_str), Some("tx-1"));
        assert_eq!(gff.parents, ["tx,1", "gene:1"]);

        let gtf = parse_gtf_attributes(
            "gene_id \"G1\"; transcript_id \"T1\"; description \"left;right=A\"; tag \"basic\"; tag \"CCDS\";",
        )
        .unwrap();
        assert_eq!(gtf.values.get("gene_id").map(String::as_str), Some("G1"));
        assert_eq!(
            gtf.values.get("transcript_id").map(String::as_str),
            Some("T1")
        );
        assert_eq!(
            gtf.values.get("description").map(String::as_str),
            Some("left;right=A")
        );
        assert_eq!(
            gtf.values.get("tag").map(String::as_str),
            Some("basic,CCDS")
        );

        let auto = parse_feature_row(
            "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; description \"A=B;C\";",
            1,
            AnnotationFormat::Auto,
        )
        .unwrap();
        assert_eq!(auto.format, AnnotationFormat::Gtf);
        assert_eq!(
            auto.attributes.get("description").map(String::as_str),
            Some("A=B;C")
        );
    }

    #[test]
    fn gff3_graph_preserves_encoded_commas_and_rejects_invalid_nodes() {
        let valid = parse_rows(
            concat!(
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1\n",
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx%2C1;Parent=gene1\n",
                "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=tx%2C1\n",
            ),
            AnnotationFormat::Gff3,
        );
        let transcripts = build_transcripts(&valid, &GffToBiggOptions::default()).unwrap();
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].name, "tx,1");
        assert_eq!(transcripts[0].metadata().gene_id(), Some("gene1"));

        let orphan = parse_rows(
            "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=missing\n",
            AnnotationFormat::Gff3,
        );
        assert!(build_transcripts(&orphan, &GffToBiggOptions::default())
            .unwrap_err()
            .to_string()
            .contains("no declared transcript feature"));

        let duplicate = parse_rows(
            concat!(
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx1\n",
                "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=tx1\n",
                "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=tx1\n",
            ),
            AnnotationFormat::Gff3,
        );
        assert!(build_transcripts(&duplicate, &GffToBiggOptions::default())
            .unwrap_err()
            .to_string()
            .contains("duplicate transcript identity"));

        let control = parse_rows(
            concat!(
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx%09bad\n",
                "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=tx%09bad\n",
            ),
            AnnotationFormat::Gff3,
        );
        assert!(build_transcripts(&control, &GffToBiggOptions::default())
            .unwrap_err()
            .to_string()
            .contains("control character"));
    }

    #[test]
    fn annotation_reader_accepts_plain_and_gzip_without_changing_api() {
        let root = temp_dir("plain-gzip");
        let plain = root.join("annotation.gtf");
        let gzip = root.join("annotation.gtf.gz");
        let annotation = concat!(
            "chr1\ttest\texon\t101\t103\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1.1\";\n",
            "chr1\ttest\texon\t201\t203\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1.1\";\n",
        );
        fs::write(&plain, annotation).unwrap();
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(annotation.as_bytes()).unwrap();
        fs::write(&gzip, encoder.finish().unwrap()).unwrap();

        let plain_transcripts =
            read_annotation_transcripts(&plain, &GffToBiggOptions::default()).unwrap();
        let gzip_transcripts =
            read_annotation_transcripts(&gzip, &GffToBiggOptions::default()).unwrap();
        assert_eq!(plain_transcripts, gzip_transcripts);
        assert_eq!(plain_transcripts[0].name, "TX1.1");
        assert_eq!(plain_transcripts[0].exons.len(), 2);

        let fake_gzip = root.join("not-gzip.gtf.gz");
        fs::write(&fake_gzip, annotation).unwrap();
        let error =
            read_annotation_transcripts(&fake_gzip, &GffToBiggOptions::default()).unwrap_err();
        assert!(format!("{error:#}").contains("not gzip-compressed"));
        let recovering_error = read_annotation_transcripts_with_policy(
            &fake_gzip,
            &GffToBiggOptions::default(),
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap_err();
        assert!(format!("{recovering_error:#}").contains("not gzip-compressed"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn filtered_reader_keeps_only_selected_gff3_transcript_parents() {
        let root = temp_dir("filtered");
        let path = root.join("annotation.gff3");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1\n",
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n",
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx2;Parent=gene1\n",
                "chr1\ttest\texon\t1\t10\t.\t+\t.\tParent=tx1,tx2\n",
                "chr1\ttest\tCDS\t2\t9\t.\t+\t0\tParent=tx1,tx2\n",
            ),
        )
        .unwrap();

        let transcripts = read_annotation_transcripts_where(
            &path,
            &GffToBiggOptions {
                format: AnnotationFormat::Gff3,
                gene_key: "ID".to_owned(),
            },
            true,
            |id| id == "tx2",
        )
        .unwrap();
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].name, "tx2");
        assert_eq!(transcripts[0].metadata().gene_id(), Some("gene1"));

        let empty =
            read_annotation_transcripts_where(&path, &GffToBiggOptions::default(), false, |_| {
                false
            })
            .unwrap();
        assert!(empty.is_empty());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_skips_only_irrelevant_bad_rows() {
        let root = temp_dir("recover-irrelevant");
        let path = root.join("annotation.gff3");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx1\n",
                "chr1\ttest\tCDS\tnot-a-coordinate\t90\t.\t+\t0\tParent=tx1\n",
                "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=tx1\n",
            ),
        )
        .unwrap();

        let recovered = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions::default(),
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap();
        assert_eq!(recovered.transcripts.len(), 1);
        assert_eq!(recovered.transcripts[0].name, "tx1");
        assert_eq!(recovered.rejected_transcripts, 0);
        assert_eq!(recovered.rejected_records.len(), 1);
        assert_eq!(recovered.rejected_records[0].line, 2);
        assert_eq!(recovered.rejected_records[0].kind, "ignored_feature");

        let strict = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions::default(),
            InvalidAnnotationPolicy::Fail,
        )
        .unwrap_err();
        assert!(format!("{strict:#}").contains("not-a-coordinate"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_quarantines_whole_identifiable_bad_transcript() {
        let root = temp_dir("recover-model");
        let path = root.join("annotation.gtf");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"bad_tx\";\n",
                "chr1\ttest\texon\tnot-a-coordinate\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"bad_tx\";\n",
                "chr1\ttest\ttranscript\t201\t250\t.\t+\t.\tgene_id \"G2\"; transcript_id \"good_tx\";\n",
                "chr1\ttest\texon\t201\t250\t.\t+\t.\tgene_id \"G2\"; transcript_id \"good_tx\";\n",
            ),
        )
        .unwrap();

        let recovered = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions {
                format: AnnotationFormat::Gtf,
                gene_key: "ID".to_owned(),
            },
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap();
        assert_eq!(recovered.transcripts.len(), 1);
        assert_eq!(recovered.transcripts[0].name, "good_tx");
        assert_eq!(recovered.rejected_transcripts, 1);
        assert_eq!(recovered.rejected_records.len(), 1);
        assert_eq!(recovered.rejected_records[0].line, 2);
        assert_eq!(recovered.rejected_records[0].kind, "parse");
        assert_eq!(recovered.rejected_records[0].transcript_ids, ["bad_tx"]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_refuses_unowned_bad_exon() {
        let root = temp_dir("recover-unsafe");
        let path = root.join("annotation.gtf");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\texon\tnot-a-coordinate\t20\t.\t+\t.\tgene_id \"G1\";\n",
                "chr1\ttest\texon\t101\t120\t.\t+\t.\tgene_id \"G2\"; transcript_id \"good_tx\";\n",
            ),
        )
        .unwrap();
        let error = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions {
                format: AnnotationFormat::Gtf,
                gene_key: "ID".to_owned(),
            },
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("refusing unsafe recovery"), "{message}");
        assert!(
            message.contains("no recoverable transcript identity"),
            "{message}"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_refuses_invalid_feature_type() {
        let root = temp_dir("recover-feature-type");
        let path = root.join("annotation.gff3");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\t.\tnot-a-coordinate\t20\t.\t+\t.\tParent=good_tx\n",
                "chr1\ttest\tmRNA\t101\t120\t.\t+\t.\tID=good_tx\n",
                "chr1\ttest\texon\t101\t120\t.\t+\t.\tParent=good_tx\n",
            ),
        )
        .unwrap();
        let error = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions::default(),
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("refusing unsafe recovery"), "{message}");
        assert!(message.contains("invalid feature type"), "{message}");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_isolates_graph_invalid_model() {
        let root = temp_dir("recover-graph");
        let path = root.join("annotation.gff3");
        fs::write(
            &path,
            concat!(
                "chr1\ttest\texon\t1\t20\t.\t+\t.\tParent=orphan\n",
                "chr1\ttest\tmRNA\t101\t120\t.\t+\t.\tID=good_tx\n",
                "chr1\ttest\texon\t101\t120\t.\t+\t.\tParent=good_tx\n",
            ),
        )
        .unwrap();
        let recovered = read_annotation_transcripts_with_policy(
            &path,
            &GffToBiggOptions {
                format: AnnotationFormat::Gff3,
                gene_key: "ID".to_owned(),
            },
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap();
        assert_eq!(recovered.transcripts.len(), 1);
        assert_eq!(recovered.transcripts[0].name, "good_tx");
        assert_eq!(recovered.rejected_transcripts, 1);
        assert_eq!(recovered.rejected_records.len(), 1);
        assert_eq!(recovered.rejected_records[0].kind, "model");
        assert_eq!(recovered.rejected_records[0].transcript_ids, ["orphan"]);
        assert!(recovered.rejected_records[0]
            .reason
            .contains("no declared transcript feature"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn recovering_reader_rejects_empty_or_all_quarantined_catalogs() {
        let root = temp_dir("recover-empty");
        let no_exons = root.join("no-exons.gff3");
        fs::write(&no_exons, "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1\n").unwrap();
        let error = read_annotation_transcripts_with_policy(
            &no_exons,
            &GffToBiggOptions::default(),
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains("no exon features"));

        let all_bad = root.join("all-bad.gtf");
        fs::write(
            &all_bad,
            "chr1\ttest\texon\tnot-a-coordinate\t20\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\n",
        )
        .unwrap();
        let error = read_annotation_transcripts_with_policy(
            &all_bad,
            &GffToBiggOptions {
                format: AnnotationFormat::Gtf,
                gene_key: "ID".to_owned(),
            },
            InvalidAnnotationPolicy::Skip,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("no valid transcript models"), "{message}");
        assert!(message.contains("rejected_transcripts=1"), "{message}");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rejected_annotation_tsv_escapes_free_text() {
        let records = [RejectedAnnotationRecord {
            path: PathBuf::from("annotation.gtf"),
            line: 7,
            transcript_ids: vec!["tx,1".to_owned()],
            kind: "parse",
            reason: "bad\tfield\ncontinued".to_owned(),
        }];
        let mut output = Vec::new();
        write_rejected_annotations_tsv_to_writer(&mut output, &records).unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            concat!(
                "source_path\tanchor_line\ttranscript_ids_json\tkind\treason\n",
                "annotation.gtf\t7\t[\"tx,1\"]\tparse\tbad\\tfield\\ncontinued\n",
            )
        );
    }
}
