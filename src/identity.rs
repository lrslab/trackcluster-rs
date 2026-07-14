//! Stable catalog identities and the versioned embedded-read payload codec.

use std::collections::HashMap;

use thiserror::Error;

use crate::model::{Strand, Transcript};

pub(crate) const NOVEL_ISOFORM_PREFIX: &str = "tc_novel_v1:";
const NAME2_PREFIX: &str = "tc_name2_v1:";

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub(crate) enum IdentityError {
    #[error("{kind} id must not be empty at index {index}")]
    EmptyId { kind: &'static str, index: usize },
    #[error(
        "duplicate {kind} id {id:?} at indices {first_index} and {second_index}; ids must be globally unique"
    )]
    DuplicateId {
        kind: &'static str,
        id: String,
        first_index: usize,
        second_index: usize,
    },
    #[error(
        "reference id {id:?} at index {index} uses reserved novel-isoform namespace {NOVEL_ISOFORM_PREFIX:?}"
    )]
    ReservedReferenceId { id: String, index: usize },
    #[error("embedded read id at position {index} is empty")]
    EmptyEmbeddedReadId { index: usize },
    #[error("malformed percent escape in name2 token {token:?} at byte {offset}")]
    MalformedPercentEscape { token: String, offset: usize },
    #[error("name2 token is not valid UTF-8 after percent decoding: {token:?}")]
    InvalidUtf8 { token: String },
}

pub(crate) fn gene_id(tx: &Transcript) -> &str {
    tx.metadata()
        .gene_id()
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
        .unwrap_or("none")
}

fn hex_bytes(value: &str) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(value.len() * 2);
    for byte in value.as_bytes() {
        out.push(HEX[(byte >> 4) as usize] as char);
        out.push(HEX[(byte & 0x0f) as usize] as char);
    }
    out
}

/// Build a collision-free, human-auditable structural ID for a novel isoform.
///
/// The gene and chromosome are UTF-8 hex encoded, while strand and every exon
/// boundary are serialized literally. Unlike a truncated hash, this encoding
/// cannot silently collide for two distinct `(gene, chromosome, strand, exons)`
/// tuples.
pub(crate) fn novel_isoform_id(tx: &Transcript) -> String {
    let strand = match tx.strand {
        Strand::Plus => 'p',
        Strand::Minus => 'm',
        Strand::Unknown => 'u',
    };
    let exons = tx
        .exons
        .iter()
        .map(|exon| format!("{}-{}", exon.start.get(), exon.end.get()))
        .collect::<Vec<_>>()
        .join(";");
    format!(
        "{NOVEL_ISOFORM_PREFIX}{}:{}:{strand}:{exons}",
        hex_bytes(gene_id(tx)),
        hex_bytes(&tx.chrom)
    )
}

fn validate_unique_ids<'a>(
    kind: &'static str,
    ids: impl IntoIterator<Item = &'a str>,
) -> Result<(), IdentityError> {
    let mut seen: HashMap<&str, usize> = HashMap::new();
    for (index, id) in ids.into_iter().enumerate() {
        if id.is_empty() {
            return Err(IdentityError::EmptyId { kind, index });
        }
        if let Some(first_index) = seen.insert(id, index) {
            return Err(IdentityError::DuplicateId {
                kind,
                id: id.to_owned(),
                first_index,
                second_index: index,
            });
        }
    }
    Ok(())
}

pub(crate) fn validate_reference_ids(references: &[Transcript]) -> Result<(), IdentityError> {
    validate_unique_ids(
        "reference isoform",
        references.iter().map(|tx| tx.name.as_str()),
    )?;
    for (index, reference) in references.iter().enumerate() {
        if reference.name.starts_with(NOVEL_ISOFORM_PREFIX) {
            return Err(IdentityError::ReservedReferenceId {
                id: reference.name.clone(),
                index,
            });
        }
    }
    Ok(())
}

pub(crate) fn validate_isoform_ids(isoforms: &[Transcript]) -> Result<(), IdentityError> {
    validate_unique_ids("isoform", isoforms.iter().map(|tx| tx.name.as_str()))
}

pub(crate) fn validate_read_ids(reads: &[Transcript]) -> Result<(), IdentityError> {
    for (index, read) in reads.iter().enumerate() {
        if read.name.is_empty() {
            return Err(IdentityError::EmptyId {
                kind: "read",
                index,
            });
        }
    }
    Ok(())
}

pub(crate) fn transcript_order(left: &Transcript, right: &Transcript) -> std::cmp::Ordering {
    left.chrom
        .cmp(&right.chrom)
        .then_with(|| left.tx_start.cmp(&right.tx_start))
        .then_with(|| left.tx_end.cmp(&right.tx_end))
        .then_with(|| left.strand.cmp(&right.strand))
        .then_with(|| left.name.cmp(&right.name))
        .then_with(|| left.exons.cmp(&right.exons))
}

fn is_unreserved(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'.' | b'_' | b'~' | b':')
}

fn percent_encode(value: &str) -> String {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut out = String::with_capacity(value.len());
    for byte in value.as_bytes() {
        if is_unreserved(*byte) {
            out.push(*byte as char);
        } else {
            out.push('%');
            out.push(HEX[(byte >> 4) as usize] as char);
            out.push(HEX[(byte & 0x0f) as usize] as char);
        }
    }
    out
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

fn percent_decode(token: &str) -> Result<String, IdentityError> {
    let bytes = token.as_bytes();
    let mut decoded = Vec::with_capacity(bytes.len());
    let mut index = 0usize;
    while index < bytes.len() {
        if bytes[index] != b'%' {
            decoded.push(bytes[index]);
            index += 1;
            continue;
        }
        let Some(high) = bytes.get(index + 1).copied().and_then(hex_value) else {
            return Err(IdentityError::MalformedPercentEscape {
                token: token.to_owned(),
                offset: index,
            });
        };
        let Some(low) = bytes.get(index + 2).copied().and_then(hex_value) else {
            return Err(IdentityError::MalformedPercentEscape {
                token: token.to_owned(),
                offset: index,
            });
        };
        decoded.push((high << 4) | low);
        index += 3;
    }
    String::from_utf8(decoded).map_err(|_| IdentityError::InvalidUtf8 {
        token: token.to_owned(),
    })
}

pub(crate) fn encode_name2<'a>(
    read_ids: impl IntoIterator<Item = &'a str>,
    coverage: f64,
) -> Result<String, IdentityError> {
    let mut encoded = Vec::new();
    for (index, read_id) in read_ids.into_iter().enumerate() {
        if read_id.is_empty() {
            return Err(IdentityError::EmptyEmbeddedReadId { index });
        }
        encoded.push(percent_encode(read_id));
    }
    Ok(format!("{NAME2_PREFIX}{},|{coverage}", encoded.join(",")))
}

/// Decode both the safe versioned payload and the legacy comma-separated form.
pub(crate) fn decode_name2(payload: &str) -> Result<Vec<String>, IdentityError> {
    if payload == "none" || payload.is_empty() || payload.starts_with('|') {
        return Ok(Vec::new());
    }

    if let Some(versioned) = payload.strip_prefix(NAME2_PREFIX) {
        let read_part = versioned.split_once(",|").map_or(versioned, |(ids, _)| ids);
        if read_part.is_empty() {
            return Ok(Vec::new());
        }
        return read_part
            .split(',')
            .enumerate()
            .map(|(index, token)| {
                if token.is_empty() {
                    Err(IdentityError::EmptyEmbeddedReadId { index })
                } else {
                    percent_decode(token)
                }
            })
            .collect();
    }

    // Compatibility reader for pre-v1 payloads. Legacy payloads cannot
    // represent a comma inside an ID; new writers always emit the safe form.
    if !payload.contains(',') {
        return Ok(Vec::new());
    }
    let read_part = payload.split_once(",|").map_or(payload, |(ids, _)| ids);
    Ok(read_part
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(str::to_owned)
        .collect())
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval};

    use super::*;

    fn tx(name: &str, gene: &str) -> Transcript {
        Transcript::new(
            "chr/一".to_owned(),
            Strand::Minus,
            Coord::new(10),
            Coord::new(40),
            name.to_owned(),
            vec![
                Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
                Interval::new(Coord::new(30), Coord::new(40)).unwrap(),
            ],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(10),
                thick_end: Coord::new(40),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    gene.to_owned(),
                ],
            },
        )
        .unwrap()
    }

    #[test]
    fn novel_id_is_structural_namespaced_and_unicode_safe() {
        let left = novel_isoform_id(&tx("representative-a", "基因/A"));
        let right = novel_isoform_id(&tx("representative-b", "基因/A"));
        assert_eq!(left, right);
        assert!(left.starts_with(NOVEL_ISOFORM_PREFIX));
        assert!(!left.contains("representative"));
        assert!(!left.contains('基'));
    }

    #[test]
    fn name2_round_trips_reserved_delimiters_and_unicode() {
        let ids = ["plain", "comma,id", "pipe|id", "percent%id", "讀取::一"];
        let payload = encode_name2(ids, 2.5).unwrap();
        assert!(payload.starts_with(NAME2_PREFIX));
        assert!(!payload.contains("comma,id"));
        assert_eq!(decode_name2(&payload).unwrap(), ids);
    }

    #[test]
    fn rejects_duplicate_and_reserved_reference_ids() {
        let duplicate = vec![tx("same", "G"), tx("same", "G")];
        assert!(matches!(
            validate_reference_ids(&duplicate),
            Err(IdentityError::DuplicateId { .. })
        ));

        let reserved = vec![tx(&format!("{NOVEL_ISOFORM_PREFIX}claimed"), "G")];
        assert!(matches!(
            validate_reference_ids(&reserved),
            Err(IdentityError::ReservedReferenceId { .. })
        ));
    }
}
