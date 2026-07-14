#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct OrderedAlignmentScore {
    matches: usize,
    total_delta: u64,
}

impl OrderedAlignmentScore {
    fn with_match(self, delta: u64) -> Self {
        Self {
            matches: self.matches + 1,
            total_delta: self.total_delta.saturating_add(delta),
        }
    }

    fn is_better_than(self, other: Self) -> bool {
        self.matches > other.matches
            || (self.matches == other.matches && self.total_delta < other.total_delta)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum OrderedAlignmentStep {
    Match,
    SkipRight,
    #[default]
    SkipLeft,
}

/// Return a maximum-cardinality, minimum-total-delta ordered alignment.
///
/// `match_delta(left_idx, right_idx)` returns the distance for a permitted
/// match and `None` for an incompatible pair. Each item can occur in at most
/// one result pair. Exact score ties keep the previously established prefix
/// alignment, making ambiguous repeated offsets deterministic.
pub(crate) fn ordered_one_to_one_matches_by<F>(
    left_len: usize,
    right_len: usize,
    mut match_delta: F,
) -> Vec<(usize, usize)>
where
    F: FnMut(usize, usize) -> Option<u64>,
{
    let columns = right_len + 1;
    let mut scores = vec![OrderedAlignmentScore::default(); (left_len + 1) * columns];
    let mut steps = vec![OrderedAlignmentStep::default(); scores.len()];

    for step in steps.iter_mut().take(columns).skip(1) {
        *step = OrderedAlignmentStep::SkipRight;
    }

    for left_idx in 1..=left_len {
        for right_idx in 1..=right_len {
            let cell = left_idx * columns + right_idx;
            let mut best_score = scores[(left_idx - 1) * columns + right_idx];
            let mut best_step = OrderedAlignmentStep::SkipLeft;

            let skip_right_score = scores[left_idx * columns + (right_idx - 1)];
            if skip_right_score.is_better_than(best_score) {
                best_score = skip_right_score;
                best_step = OrderedAlignmentStep::SkipRight;
            }

            if let Some(delta) = match_delta(left_idx - 1, right_idx - 1) {
                let match_score =
                    scores[(left_idx - 1) * columns + (right_idx - 1)].with_match(delta);
                if match_score.is_better_than(best_score) {
                    best_score = match_score;
                    best_step = OrderedAlignmentStep::Match;
                }
            }

            scores[cell] = best_score;
            steps[cell] = best_step;
        }
    }

    let mut matches = Vec::with_capacity(scores.last().map_or(0, |score| score.matches));
    let mut left_idx = left_len;
    let mut right_idx = right_len;
    while left_idx > 0 && right_idx > 0 {
        match steps[left_idx * columns + right_idx] {
            OrderedAlignmentStep::Match => {
                matches.push((left_idx - 1, right_idx - 1));
                left_idx -= 1;
                right_idx -= 1;
            }
            OrderedAlignmentStep::SkipRight => right_idx -= 1,
            OrderedAlignmentStep::SkipLeft => left_idx -= 1,
        }
    }
    matches.reverse();
    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maximizes_cardinality_before_minimizing_distance() {
        let distances = [[Some(0), Some(1)], [None, Some(100)]];
        let matches = ordered_one_to_one_matches_by(2, 2, |left, right| distances[left][right]);
        assert_eq!(matches, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn exact_ties_keep_the_earliest_prefix_match() {
        let matches = ordered_one_to_one_matches_by(2, 1, |_, _| Some(5));
        assert_eq!(matches, vec![(0, 0)]);
    }
}
