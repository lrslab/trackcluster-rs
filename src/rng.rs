//! Deterministic, unbiased sampling helpers shared by flow and clusterj.

#[derive(Clone, Copy, Debug)]
pub(crate) struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    pub(crate) fn gen_below(&mut self, upper: usize) -> usize {
        debug_assert!(upper > 0);
        let upper = upper as u64;
        // Reject the short leading interval before taking the remainder. This
        // removes modulo bias while preserving the sampler's deterministic
        // seed and bounded-memory behavior.
        let threshold = upper.wrapping_neg() % upper;
        loop {
            let value = self.next_u64();
            if value >= threshold {
                return (value % upper) as usize;
            }
        }
    }
}

pub(crate) fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 14695981039346656037;
    update_fnv1a64(&mut hash, bytes);
    hash
}

pub(crate) fn update_fnv1a64(hash: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *hash ^= b as u64;
        *hash = hash.wrapping_mul(1099511628211);
    }
}

/// Reservoir-sample `target` indices from `0..len` in encounter order.
pub(crate) fn reservoir_sample_indices(len: usize, target: usize, seed: u64) -> Vec<usize> {
    if target == 0 || len == 0 {
        return Vec::new();
    }
    if len <= target {
        return (0..len).collect();
    }

    let mut rng = Lcg64::new(seed);
    let mut sampled: Vec<usize> = (0..target).collect();
    for next in target..len {
        let idx = rng.gen_below(next + 1);
        if idx < target {
            sampled[idx] = next;
        }
    }
    sampled.sort_unstable();
    sampled
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reservoir_keeps_everything_when_under_the_cap() {
        assert_eq!(reservoir_sample_indices(3, 10, 1), vec![0, 1, 2]);
    }

    #[test]
    fn reservoir_is_deterministic_and_sized_to_the_cap() {
        let first = reservoir_sample_indices(20, 5, 7);
        let second = reservoir_sample_indices(20, 5, 7);
        assert_eq!(first, second);
        assert_eq!(first.len(), 5);
        assert!(first.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(first.iter().all(|&idx| idx < 20));
    }
}
