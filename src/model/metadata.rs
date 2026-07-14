//! Typed access to TrackCluster's bigGenePred extension columns.

use super::{Coord, Interval, Strand, Transcript};

const NAME2: usize = 0;
const CDS_START_STAT: usize = 1;
const CDS_END_STAT: usize = 2;
const EXON_FRAMES: usize = 3;
const TRANSCRIPT_TYPE: usize = 4;
const GENE_ID: usize = 5;
const SAMPLE_GROUP: usize = 6;
const RESERVED: usize = 7;
const STANDARD_FIELD_COUNT: usize = 8;

fn optional_field(fields: &[String], index: usize) -> Option<&str> {
    fields
        .get(index)
        .map(String::as_str)
        .filter(|value| !value.is_empty() && *value != "none")
}

fn raw_field(fields: &[String], index: usize) -> Option<&str> {
    fields.get(index).map(String::as_str)
}

/// Immutable view of transcript extension metadata.
#[derive(Clone, Copy, Debug)]
pub struct TrackMetadataRef<'a> {
    fields: &'a [String],
}

impl<'a> TrackMetadataRef<'a> {
    pub(crate) fn new(fields: &'a [String]) -> Self {
        Self { fields }
    }

    /// Serialized subread/coverage payload.
    pub fn name2(self) -> Option<&'a str> {
        optional_field(self.fields, NAME2)
    }

    /// Serialized subread/coverage column when it is physically present,
    /// including legacy sentinel values such as `none`.
    pub fn name2_field(self) -> Option<&'a str> {
        raw_field(self.fields, NAME2)
    }

    /// CDS start-status token.
    pub fn cds_start_status(self) -> Option<&'a str> {
        optional_field(self.fields, CDS_START_STAT)
    }

    /// CDS end-status token.
    pub fn cds_end_status(self) -> Option<&'a str> {
        optional_field(self.fields, CDS_END_STAT)
    }

    /// Serialized exon-frame list.
    pub fn exon_frames(self) -> Option<&'a str> {
        optional_field(self.fields, EXON_FRAMES)
    }

    /// Biological transcript-type annotation, such as `isoform_anno`.
    pub fn transcript_type(self) -> Option<&'a str> {
        optional_field(self.fields, TRANSCRIPT_TYPE)
    }

    /// Transcript-type column when it is physically present.
    pub fn transcript_type_field(self) -> Option<&'a str> {
        raw_field(self.fields, TRANSCRIPT_TYPE)
    }

    /// Biological gene identifier field.
    pub fn gene_id(self) -> Option<&'a str> {
        optional_field(self.fields, GENE_ID)
    }

    /// Gene-ID column when it is physically present, including `none`.
    pub fn gene_id_field(self) -> Option<&'a str> {
        raw_field(self.fields, GENE_ID)
    }

    /// Sample/group annotation field.
    pub fn sample_group(self) -> Option<&'a str> {
        optional_field(self.fields, SAMPLE_GROUP)
    }

    /// Reserved TrackCluster annotation field.
    pub fn reserved(self) -> Option<&'a str> {
        optional_field(self.fields, RESERVED)
    }
}

/// Mutable view that changes annotations without exposing numeric column indices.
#[derive(Debug)]
pub struct TrackMetadataMut<'a> {
    fields: &'a mut Vec<String>,
}

impl<'a> TrackMetadataMut<'a> {
    pub(crate) fn new(fields: &'a mut Vec<String>) -> Self {
        Self { fields }
    }

    fn set(&mut self, index: usize, value: impl Into<String>) {
        if self.fields.len() <= index {
            self.fields.resize(index + 1, "none".to_owned());
        }
        self.fields[index] = value.into();
    }

    /// Set the subread/coverage payload.
    pub fn set_name2(&mut self, value: impl Into<String>) {
        self.set(NAME2, value);
    }

    /// Set the biological transcript-type annotation.
    pub fn set_transcript_type(&mut self, value: impl Into<String>) {
        self.set(TRANSCRIPT_TYPE, value);
    }

    /// Set the biological gene identifier field.
    pub fn set_gene_id(&mut self, value: impl Into<String>) {
        self.set(GENE_ID, value);
    }

    /// Set the sample/group annotation field.
    pub fn set_sample_group(&mut self, value: impl Into<String>) {
        self.set(SAMPLE_GROUP, value);
    }
}

/// Owned, round-trippable codec for all standard bigGenePred extension fields.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BigGenePredAttrs {
    /// Subread/coverage payload.
    pub name2: Option<String>,
    /// CDS start-status token.
    pub cds_start_status: Option<String>,
    /// CDS end-status token.
    pub cds_end_status: Option<String>,
    /// Serialized exon-frame list.
    pub exon_frames: Option<String>,
    /// Biological transcript-type annotation.
    pub transcript_type: Option<String>,
    /// Biological gene identifier.
    pub gene_id: Option<String>,
    /// Sample/group annotation.
    pub sample_group: Option<String>,
    /// Reserved TrackCluster field.
    pub reserved: Option<String>,
    /// Extension columns after the eight standard TrackCluster fields.
    pub trailing: Vec<String>,
    original_len: usize,
    decoded_missing_fields: Vec<Option<String>>,
}

impl BigGenePredAttrs {
    /// Decode preserved BED extension columns into named fields.
    pub fn decode(fields: &[String]) -> Self {
        let get = |index| optional_field(fields, index).map(str::to_owned);
        let decoded_missing_fields = fields
            .iter()
            .take(STANDARD_FIELD_COUNT)
            .map(|value| (value.is_empty() || value == "none").then(|| value.clone()))
            .collect();
        Self {
            name2: get(NAME2),
            cds_start_status: get(CDS_START_STAT),
            cds_end_status: get(CDS_END_STAT),
            exon_frames: get(EXON_FRAMES),
            transcript_type: get(TRANSCRIPT_TYPE),
            gene_id: get(GENE_ID),
            sample_group: get(SAMPLE_GROUP),
            reserved: get(RESERVED),
            trailing: fields
                .get(STANDARD_FIELD_COUNT..)
                .unwrap_or_default()
                .to_vec(),
            original_len: fields.len(),
            decoded_missing_fields,
        }
    }

    /// Encode named fields back to extension columns, preserving absent columns
    /// and any trailing fields from [`Self::decode`].
    pub fn encode(&self) -> Vec<String> {
        let values = [
            &self.name2,
            &self.cds_start_status,
            &self.cds_end_status,
            &self.exon_frames,
            &self.transcript_type,
            &self.gene_id,
            &self.sample_group,
            &self.reserved,
        ];
        let last_set = values
            .iter()
            .rposition(|value| value.is_some())
            .map_or(0, |index| index + 1);
        let decoded_standard_len = self.original_len.min(STANDARD_FIELD_COUNT);
        let standard_len = if self.trailing.is_empty() {
            decoded_standard_len.max(last_set)
        } else {
            STANDARD_FIELD_COUNT
        };
        let mut fields: Vec<String> = values[..standard_len]
            .iter()
            .enumerate()
            .map(|(index, value)| {
                (**value).clone().unwrap_or_else(|| {
                    self.decoded_missing_fields
                        .get(index)
                        .and_then(Clone::clone)
                        .unwrap_or_else(|| "none".to_owned())
                })
            })
            .collect();
        fields.extend(self.trailing.iter().cloned());
        fields
    }
}

/// Immutable geometry view, separate from mutable pipeline annotations.
#[derive(Clone, Copy, Debug)]
pub struct TranscriptGeometry<'a> {
    /// Reference sequence name.
    pub chrom: &'a str,
    /// Genomic strand.
    pub strand: Strand,
    /// Transcript start.
    pub tx_start: Coord,
    /// Transcript end.
    pub tx_end: Coord,
    /// Ordered exon intervals.
    pub exons: &'a [Interval],
}

impl<'a> From<&'a Transcript> for TranscriptGeometry<'a> {
    fn from(transcript: &'a Transcript) -> Self {
        Self {
            chrom: &transcript.chrom,
            strand: transcript.strand,
            tx_start: transcript.tx_start,
            tx_end: transcript.tx_end,
            exons: &transcript.exons,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codec_round_trips_missing_and_trailing_fields() {
        let empty = BigGenePredAttrs::decode(&[]);
        assert!(empty.encode().is_empty());

        let fields = vec![
            "reads,|2".to_owned(),
            "none".to_owned(),
            "none".to_owned(),
            "-1,".to_owned(),
            "isoform_anno".to_owned(),
            "GENE1".to_owned(),
            "group".to_owned(),
            "none".to_owned(),
            "future".to_owned(),
        ];
        assert_eq!(BigGenePredAttrs::decode(&fields).encode(), fields);
    }

    #[test]
    fn codec_preserves_explicitly_empty_standard_fields() {
        let fields = vec![
            "".to_owned(),
            "none".to_owned(),
            "".to_owned(),
            "-1,".to_owned(),
            "".to_owned(),
            "GENE1".to_owned(),
            "".to_owned(),
            "".to_owned(),
        ];

        assert_eq!(BigGenePredAttrs::decode(&fields).encode(), fields);
    }

    #[test]
    fn codec_places_new_trailing_fields_after_all_standard_fields() {
        let attrs = BigGenePredAttrs {
            trailing: vec!["future".to_owned()],
            ..BigGenePredAttrs::default()
        };

        let encoded = attrs.encode();
        assert_eq!(encoded.len(), STANDARD_FIELD_COUNT + 1);
        assert!(encoded[..STANDARD_FIELD_COUNT]
            .iter()
            .all(|field| field == "none"));
        assert_eq!(encoded[STANDARD_FIELD_COUNT], "future");
    }

    #[test]
    fn mutable_view_materializes_only_through_requested_field() {
        let mut fields = Vec::new();
        TrackMetadataMut::new(&mut fields).set_gene_id("GENE1");
        assert_eq!(fields.len(), GENE_ID + 1);
        assert_eq!(TrackMetadataRef::new(&fields).gene_id(), Some("GENE1"));
    }
}
