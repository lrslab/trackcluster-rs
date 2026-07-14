//! GFF3/GTF transcript annotation parsing and bigGenePred conversion.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Context;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

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

/// Read GFF3 or GTF exon annotations and build validated BED12+8 transcripts.
///
/// GFF3 relationships use `ID`/`Parent`; GTF relationships use
/// `transcript_id`/`gene_id`. Input order does not affect output order.
pub fn read_annotation_transcripts(
    path: &Path,
    options: &GffToBiggOptions,
) -> anyhow::Result<Vec<Transcript>> {
    let file = File::open(path).with_context(|| format!("open annotation input {path:?}"))?;
    let reader = BufReader::new(file);
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
        rows.push(
            parse_feature_row(line, line_number, options.format)
                .with_context(|| format!("parse annotation {path:?}:{line_number}"))?,
        );
    }
    build_transcripts(&rows, options).with_context(|| format!("convert annotation {path:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
