use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::model::Transcript;

pub mod multi;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct CountRecord {
    pub isoform_id: String,
    pub count: f64,
}

pub(crate) fn parse_subreads(tx: &Transcript) -> Vec<&str> {
    let Some(name2) = tx.extra_fields.first() else {
        return Vec::new();
    };
    if !name2.contains(',') {
        return Vec::new();
    }
    let sub_part = name2.split(",|").next().unwrap_or(name2.as_str());
    sub_part
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .collect()
}

pub fn count_by_subreads(isoforms: &[Transcript], references: &[Transcript]) -> Vec<CountRecord> {
    let ref_names: HashSet<&str> = references.iter().map(|tx| tx.name.as_str()).collect();

    let mut occ: HashMap<&str, u32> = HashMap::new();
    for isoform in isoforms {
        for name in parse_subreads(isoform) {
            if ref_names.contains(name) {
                continue;
            }
            *occ.entry(name).or_insert(0) += 1;
        }
    }

    isoforms
        .iter()
        .map(|isoform| {
            let mut coverage = 0.0f64;
            let mut subreads = parse_subreads(isoform);
            subreads.sort_unstable();
            for name in subreads {
                if ref_names.contains(name) {
                    continue;
                }
                let denom = occ.get(name).copied().unwrap_or(0);
                if denom > 0 {
                    coverage += 1.0f64 / denom as f64;
                }
            }

            CountRecord {
                isoform_id: isoform.name.clone(),
                count: coverage,
            }
        })
        .collect()
}

pub fn read_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<(String, String)>, std::io::Error> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);

    let mut pairs: Vec<(String, String)> = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let Some((read_id, isoform_id)) = line.split_once('\t') else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid read_to_isoform line {}: {:?}", line_no + 1, line),
            ));
        };

        let read_id = read_id.trim();
        let isoform_id = isoform_id.trim();
        if read_id.is_empty() || isoform_id.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "empty read/isoform value at read_to_isoform line {}",
                    line_no + 1
                ),
            ));
        }

        pairs.push((read_id.to_owned(), isoform_id.to_owned()));
    }

    Ok(pairs)
}

pub fn count_by_read_to_isoform(
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
) -> Vec<CountRecord> {
    let mut read_occurrence: HashMap<&str, u32> = HashMap::new();
    for (read_id, _) in read_to_isoform {
        *read_occurrence.entry(read_id.as_str()).or_insert(0) += 1;
    }

    let mut counts: HashMap<&str, f64> = HashMap::new();
    for (read_id, isoform_id) in read_to_isoform {
        let denom = read_occurrence.get(read_id.as_str()).copied().unwrap_or(0);
        if denom == 0 {
            continue;
        }
        *counts.entry(isoform_id.as_str()).or_insert(0.0) += 1.0f64 / denom as f64;
    }

    isoforms
        .iter()
        .map(|isoform| CountRecord {
            isoform_id: isoform.name.clone(),
            count: counts.get(isoform.name.as_str()).copied().unwrap_or(0.0),
        })
        .collect()
}

pub fn write_counts_csv<P: AsRef<Path>>(
    path: P,
    records: &[CountRecord],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(&mut writer, "isoform_id,count")?;
    for record in records {
        writeln!(&mut writer, "{},{}", record.isoform_id, record.count)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(name: &str, exons: &[(u32, u32)], name2: &str) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: vec![name2.to_owned()],
            },
        )
        .unwrap()
    }

    #[test]
    fn counts_split_duplicates_across_isoforms() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,|0"),
            make_tx("iso3", &[(0, 10)], ",|0"), // empty
        ];

        let records = count_by_subreads(&isoforms, &references);
        let iso1 = records.iter().find(|r| r.isoform_id == "iso1").unwrap();
        let iso2 = records.iter().find(|r| r.isoform_id == "iso2").unwrap();
        let iso3 = records.iter().find(|r| r.isoform_id == "iso3").unwrap();

        assert!((iso1.count - 1.5).abs() < 1e-9);
        assert!((iso2.count - 0.5).abs() < 1e-9);
        assert!((iso3.count - 0.0).abs() < 1e-9);
    }

    #[test]
    fn mapping_counts_match_subread_counts() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,|0"),
        ];
        let pairs = vec![
            ("r1".to_owned(), "iso1".to_owned()),
            ("r2".to_owned(), "iso1".to_owned()),
            ("r2".to_owned(), "iso2".to_owned()),
        ];

        let by_subreads = count_by_subreads(&isoforms, &references);
        let by_mapping = count_by_read_to_isoform(&isoforms, &pairs);
        assert_eq!(by_subreads.len(), by_mapping.len());

        for (left, right) in by_subreads.iter().zip(by_mapping.iter()) {
            assert_eq!(left.isoform_id, right.isoform_id);
            assert!((left.count - right.count).abs() < 1e-9);
        }
    }
}
