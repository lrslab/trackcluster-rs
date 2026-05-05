use std::time::Instant;

use trackcluster_rs::cluster::{cluster_overlap, clusterj::Name2Mode};
use trackcluster_rs::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

#[derive(Clone)]
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    fn gen_range_u32(&mut self, low: u32, high: u32) -> u32 {
        let span = (high - low) as u64;
        low + (self.next_u64() % span) as u32
    }
}

fn bed12_attrs(score: u32, tx_start: u32, tx_end: u32, ttype: &str) -> Bed12Attrs {
    let mut extra_fields = vec!["none".to_owned(); 8];
    extra_fields[4] = ttype.to_owned();
    Bed12Attrs {
        score,
        thick_start: Coord::new(tx_start),
        thick_end: Coord::new(tx_end),
        item_rgb: "0".to_owned(),
        extra_fields,
    }
}

fn make_tx(
    chrom: &str,
    strand: Strand,
    name: String,
    exons: Vec<(u32, u32)>,
    ttype: &str,
    score: u32,
) -> Transcript {
    let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
    let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
    let exon_intervals = exons
        .into_iter()
        .map(|(s, e)| Interval::new(Coord::new(s), Coord::new(e)).expect("valid exon"))
        .collect::<Vec<_>>();
    Transcript::new(
        chrom.to_owned(),
        strand,
        Coord::new(tx_start),
        Coord::new(tx_end),
        name,
        exon_intervals,
        bed12_attrs(score, tx_start, tx_end, ttype),
    )
    .expect("valid transcript")
}

fn exon_chain(tx_start: u32, exon_len: u32, gap_len: u32, exon_count: usize) -> Vec<(u32, u32)> {
    let mut exons = Vec::with_capacity(exon_count);
    let mut cursor = tx_start;
    for _ in 0..exon_count {
        exons.push((cursor, cursor + exon_len));
        cursor = cursor + exon_len + gap_len;
    }
    exons
}

fn make_cluster_overlap_inputs(
    seed: u64,
    refs_len: usize,
    reads_len: usize,
    locus_start: u32,
    locus_span: u32,
) -> (Vec<Transcript>, Vec<Transcript>) {
    let mut rng = Lcg64::new(seed);
    let mut refs = Vec::with_capacity(refs_len);
    let mut reads = Vec::with_capacity(reads_len);

    for i in 0..refs_len {
        let tx_start = locus_start + rng.gen_range_u32(0, locus_span / 4);
        let exons = exon_chain(tx_start, 60, 90, 3);
        refs.push(make_tx(
            "chr1",
            Strand::Plus,
            format!("ref{i}"),
            exons,
            "isoform_anno",
            100,
        ));
    }

    for i in 0..reads_len {
        let tx_start = locus_start + rng.gen_range_u32(0, locus_span / 4);
        let exon_len = 40 + rng.gen_range_u32(0, 30);
        let gap_len = 70 + rng.gen_range_u32(0, 50);
        let exons = exon_chain(tx_start, exon_len, gap_len, 3);
        reads.push(make_tx(
            "chr1",
            Strand::Plus,
            format!("read{i}"),
            exons,
            "nanopore_read",
            0,
        ));
    }

    (refs, reads)
}

fn env_parse<T: std::str::FromStr>(key: &str, default: T) -> T {
    std::env::var(key)
        .ok()
        .and_then(|value| value.parse::<T>().ok())
        .unwrap_or(default)
}

fn env_batch_sizes() -> Vec<usize> {
    std::env::var("TRACKCLUSTER_RS_PERF_BATCHES")
        .ok()
        .map(|value| {
            value
                .split(',')
                .filter_map(|part| part.trim().parse::<usize>().ok())
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
        .unwrap_or_else(|| vec![0, 500, 100_000])
}

#[test]
#[ignore]
fn overlap_cluster_120k_speed_probe() {
    let refs_len = env_parse("TRACKCLUSTER_RS_PERF_REFS", 500usize);
    let reads_len = env_parse("TRACKCLUSTER_RS_PERF_READS", 120_000usize);
    let locus_span = env_parse("TRACKCLUSTER_RS_PERF_LOCUS_SPAN", 500_000u32);
    let seed = env_parse("TRACKCLUSTER_RS_PERF_SEED", 17u64);
    let batch_sizes = env_batch_sizes();

    assert!(
        reads_len > 100_000,
        "speed probe should exercise more than 100k simulated reads"
    );

    let (refs, reads) = make_cluster_overlap_inputs(seed, refs_len, reads_len, 10_000, locus_span);

    for batch_size in batch_sizes {
        let started = Instant::now();
        let result = cluster_overlap::cluster_with_options(
            &reads,
            Some(&refs),
            1,
            cluster_overlap::ClusterOptions {
                batch_size,
                name2_mode: Name2Mode::Coverage,
                ..cluster_overlap::ClusterOptions::default()
            },
        );
        let elapsed = started.elapsed();

        eprintln!(
            "overlap_speed_probe\trefs={refs_len}\treads={reads_len}\tlocus_span={locus_span}\tbatch_size={batch_size}\telapsed_s={:.6}\tisoforms={}\tmappings={}\tunused={}",
            elapsed.as_secs_f64(),
            result.isoforms.len(),
            result.read_to_isoform.len(),
            result.unused.len()
        );

        assert!(!result.isoforms.is_empty());
    }
}
