use std::cmp::Ordering;

use thiserror::Error;

use super::metadata::{TrackMetadataMut, TrackMetadataRef, TranscriptGeometry};
use super::{Coord, Interval, Strand};

#[derive(Clone, Debug, PartialEq, Eq)]
/// A BED12-compatible transcript with ordered, non-overlapping exons.
///
/// [`Transcript::new`] validates the BED score, thick span, and transcript
/// geometry. Fields remain public for compatibility, so callers that mutate
/// them directly are responsible for preserving those invariants.
pub struct Transcript {
    /// Reference sequence name.
    pub chrom: String,
    /// Transcript strand.
    pub strand: Strand,
    /// Inclusive transcript start, equal to the first exon start.
    pub tx_start: Coord,
    /// Exclusive transcript end, equal to the last exon end.
    pub tx_end: Coord,
    /// Transcript identifier.
    pub name: String,
    /// BED score; validated constructors require the inclusive range 0–1000.
    pub score: u32,
    /// Coding/thick-region start.
    pub thick_start: Coord,
    /// Coding/thick-region end.
    pub thick_end: Coord,
    /// BED item RGB field.
    pub item_rgb: String,
    /// Ordered, non-empty, non-overlapping exon intervals.
    pub exons: Vec<Interval>,
    /// Optional columns following the twelve standard BED fields.
    pub extra_fields: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// Non-geometric BED12 attributes supplied to [`Transcript::new`].
pub struct Bed12Attrs {
    /// BED score in the inclusive range 0–1000.
    pub score: u32,
    /// Coding/thick-region start, or zero as part of the legacy `0, 0` sentinel.
    pub thick_start: Coord,
    /// Coding/thick-region end, or zero as part of the legacy `0, 0` sentinel.
    pub thick_end: Coord,
    /// BED item RGB value.
    pub item_rgb: String,
    /// Optional columns following the standard BED fields.
    pub extra_fields: Vec<String>,
}

#[derive(Error, Debug)]
/// Error returned when a transcript violates BED12 or geometry invariants.
pub enum TranscriptError {
    /// Transcript start is greater than transcript end.
    #[error("invalid transcript span: start {start} > end {end}")]
    InvalidSpan {
        /// Supplied transcript start.
        start: Coord,
        /// Supplied transcript end.
        end: Coord,
    },

    /// BED score is outside the standard range.
    #[error("BED score must be between 0 and 1000, got {value}")]
    ScoreOutOfRange {
        /// Supplied BED score.
        value: u32,
    },

    /// Coding/thick coordinates are ordered incorrectly or outside the transcript.
    #[error(
        "invalid thick span [{thick_start}, {thick_end}) for transcript [{tx_start}, {tx_end})"
    )]
    InvalidThickSpan {
        /// Supplied coding/thick start.
        thick_start: Coord,
        /// Supplied coding/thick end.
        thick_end: Coord,
        /// Declared transcript start.
        tx_start: Coord,
        /// Declared transcript end.
        tx_end: Coord,
    },

    /// No exons were supplied.
    #[error("expected at least 1 exon")]
    EmptyExons,

    /// An exon lies outside the declared transcript span.
    #[error("exon is outside transcript span: exon {exon:?}, transcript [{tx_start}, {tx_end})")]
    ExonOutsideSpan {
        /// Exon outside the transcript.
        exon: Interval,
        /// Declared transcript start.
        tx_start: Coord,
        /// Declared transcript end.
        tx_end: Coord,
    },

    /// An exon contains zero bases.
    #[error("transcript contains an empty exon at {exon:?}")]
    EmptyExon {
        /// Empty exon interval.
        exon: Interval,
    },

    /// Two exons overlap.
    #[error("transcript contains overlapping exons: {left:?} and {right:?}")]
    OverlappingExons {
        /// Earlier exon.
        left: Interval,
        /// Later exon that overlaps `left`.
        right: Interval,
    },

    /// The outer exon bounds do not equal the transcript span.
    #[error(
        "transcript span [{tx_start}, {tx_end}) does not match exon bounds [{exon_start}, {exon_end})"
    )]
    ExonBoundsDoNotMatchSpan {
        /// Declared transcript start.
        tx_start: Coord,
        /// Declared transcript end.
        tx_end: Coord,
        /// First exon start.
        exon_start: Coord,
        /// Last exon end.
        exon_end: Coord,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// Chromosome, strand, and ordered intron chain used for exact comparisons.
pub struct JunctionSignature {
    /// Reference sequence name.
    pub chrom: String,
    /// Transcript strand.
    pub strand: Strand,
    /// Ordered half-open intron intervals.
    pub introns: Vec<Interval>,
}

impl Transcript {
    /// Borrow typed biological/pipeline annotations.
    pub fn metadata(&self) -> TrackMetadataRef<'_> {
        TrackMetadataRef::new(&self.extra_fields)
    }

    /// Mutably borrow typed biological/pipeline annotations.
    pub fn metadata_mut(&mut self) -> TrackMetadataMut<'_> {
        TrackMetadataMut::new(&mut self.extra_fields)
    }

    /// Borrow immutable transcript geometry independently of annotations.
    pub fn geometry(&self) -> TranscriptGeometry<'_> {
        TranscriptGeometry::from(self)
    }

    /// Construct a transcript and validate its BED12 and geometry invariants.
    ///
    /// The thick span must be ordered and contained within the transcript. The
    /// legacy `0, 0` non-coding sentinel is also accepted, including for
    /// transcripts whose genomic span does not contain zero.
    pub fn new(
        chrom: String,
        strand: Strand,
        tx_start: Coord,
        tx_end: Coord,
        name: String,
        mut exons: Vec<Interval>,
        bed: Bed12Attrs,
    ) -> Result<Self, TranscriptError> {
        if tx_start > tx_end {
            return Err(TranscriptError::InvalidSpan {
                start: tx_start,
                end: tx_end,
            });
        }
        if bed.score > 1000 {
            return Err(TranscriptError::ScoreOutOfRange { value: bed.score });
        }
        let thick_is_legacy_sentinel = bed.thick_start.get() == 0 && bed.thick_end.get() == 0;
        if !thick_is_legacy_sentinel
            && (bed.thick_start > bed.thick_end
                || bed.thick_start < tx_start
                || bed.thick_end > tx_end)
        {
            return Err(TranscriptError::InvalidThickSpan {
                thick_start: bed.thick_start,
                thick_end: bed.thick_end,
                tx_start,
                tx_end,
            });
        }
        if exons.is_empty() {
            return Err(TranscriptError::EmptyExons);
        }

        exons.sort_by(|left, right| match left.start.cmp(&right.start) {
            Ordering::Equal => left.end.cmp(&right.end),
            ordering => ordering,
        });

        for exon in &exons {
            if exon.is_empty() {
                return Err(TranscriptError::EmptyExon { exon: *exon });
            }
            if exon.start < tx_start || exon.end > tx_end {
                return Err(TranscriptError::ExonOutsideSpan {
                    exon: *exon,
                    tx_start,
                    tx_end,
                });
            }
        }

        for pair in exons.windows(2) {
            if pair[1].start < pair[0].end {
                return Err(TranscriptError::OverlappingExons {
                    left: pair[0],
                    right: pair[1],
                });
            }
        }

        let exon_start = exons[0].start;
        let exon_end = exons[exons.len() - 1].end;
        if exon_start != tx_start || exon_end != tx_end {
            return Err(TranscriptError::ExonBoundsDoNotMatchSpan {
                tx_start,
                tx_end,
                exon_start,
                exon_end,
            });
        }

        Ok(Self {
            chrom,
            strand,
            tx_start,
            tx_end,
            name,
            score: bed.score,
            thick_start: bed.thick_start,
            thick_end: bed.thick_end,
            item_rgb: bed.item_rgb,
            exons,
            extra_fields: bed.extra_fields,
        })
    }

    /// Return the ordered non-empty intervals between adjacent exons.
    pub fn introns(&self) -> Vec<Interval> {
        let mut introns = Vec::new();
        for window in self.exons.windows(2) {
            let left = window[0];
            let right = window[1];
            if left.end < right.start {
                introns.push(Interval {
                    start: left.end,
                    end: right.start,
                });
            }
        }
        introns
    }

    /// Build an exact intron-chain signature.
    pub fn junction_signature(&self) -> JunctionSignature {
        JunctionSignature {
            chrom: self.chrom.clone(),
            strand: self.strand,
            introns: self.introns(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn introns_from_exons() {
        let exons = vec![
            Interval::new(Coord::new(100), Coord::new(150)).unwrap(),
            Interval::new(Coord::new(170), Coord::new(200)).unwrap(),
        ];
        let transcript = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "tx1".to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(100),
                thick_end: Coord::new(200),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap();

        assert_eq!(
            transcript.introns(),
            vec![Interval::new(Coord::new(150), Coord::new(170)).unwrap()]
        );
    }

    #[test]
    fn junction_signature_is_stable() {
        let exons_a = vec![
            Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            Interval::new(Coord::new(30), Coord::new(40)).unwrap(),
        ];
        let exons_b = vec![
            Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            Interval::new(Coord::new(30), Coord::new(40)).unwrap(),
        ];

        let a = Transcript::new(
            "chr1".to_owned(),
            Strand::Minus,
            Coord::new(10),
            Coord::new(40),
            "a".to_owned(),
            exons_a,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(10),
                thick_end: Coord::new(40),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap();

        let b = Transcript::new(
            "chr1".to_owned(),
            Strand::Minus,
            Coord::new(10),
            Coord::new(40),
            "b".to_owned(),
            exons_b,
            Bed12Attrs {
                score: 999,
                thick_start: Coord::new(10),
                thick_end: Coord::new(40),
                item_rgb: "0".to_owned(),
                extra_fields: vec!["extra".to_owned()],
            },
        )
        .unwrap();

        assert_eq!(a.junction_signature(), b.junction_signature());
    }

    fn attrs(start: u32, end: u32) -> Bed12Attrs {
        Bed12Attrs {
            score: 0,
            thick_start: Coord::new(start),
            thick_end: Coord::new(end),
            item_rgb: "0".to_owned(),
            extra_fields: Vec::new(),
        }
    }

    #[test]
    fn rejects_empty_and_overlapping_exons() {
        let empty = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "empty".to_owned(),
            vec![
                Interval::new(Coord::new(100), Coord::new(100)).unwrap(),
                Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
            ],
            attrs(100, 200),
        );
        assert!(matches!(empty, Err(TranscriptError::EmptyExon { .. })));

        let overlapping = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "overlap".to_owned(),
            vec![
                Interval::new(Coord::new(100), Coord::new(170)).unwrap(),
                Interval::new(Coord::new(140), Coord::new(200)).unwrap(),
            ],
            attrs(100, 200),
        );
        assert!(matches!(
            overlapping,
            Err(TranscriptError::OverlappingExons { .. })
        ));
    }

    #[test]
    fn rejects_exon_bounds_that_do_not_match_transcript_span() {
        let result = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "bad-bounds".to_owned(),
            vec![Interval::new(Coord::new(110), Coord::new(190)).unwrap()],
            attrs(100, 200),
        );
        assert!(matches!(
            result,
            Err(TranscriptError::ExonBoundsDoNotMatchSpan { .. })
        ));
    }

    #[test]
    fn rejects_score_outside_bed_range() {
        let mut bed = attrs(100, 200);
        bed.score = 1001;
        let result = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "bad-score".to_owned(),
            vec![Interval::new(Coord::new(100), Coord::new(200)).unwrap()],
            bed,
        );

        assert!(matches!(
            result,
            Err(TranscriptError::ScoreOutOfRange { value: 1001 })
        ));
    }

    #[test]
    fn rejects_invalid_thick_spans_and_accepts_legacy_sentinel() {
        for (thick_start, thick_end) in [(99, 200), (100, 201), (180, 170)] {
            let result = Transcript::new(
                "chr1".to_owned(),
                Strand::Plus,
                Coord::new(100),
                Coord::new(200),
                "bad-thick".to_owned(),
                vec![Interval::new(Coord::new(100), Coord::new(200)).unwrap()],
                attrs(thick_start, thick_end),
            );

            match result {
                Err(TranscriptError::InvalidThickSpan {
                    thick_start: actual_start,
                    thick_end: actual_end,
                    tx_start,
                    tx_end,
                }) => {
                    assert_eq!(actual_start, Coord::new(thick_start));
                    assert_eq!(actual_end, Coord::new(thick_end));
                    assert_eq!(tx_start, Coord::new(100));
                    assert_eq!(tx_end, Coord::new(200));
                }
                other => panic!("expected InvalidThickSpan, got {other:?}"),
            }
        }

        let sentinel = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "non-coding".to_owned(),
            vec![Interval::new(Coord::new(100), Coord::new(200)).unwrap()],
            attrs(0, 0),
        );
        assert!(sentinel.is_ok());
    }
}
