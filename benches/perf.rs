use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use trackcluster_rs::cluster::{cluster_overlap, clusterj};
use trackcluster_rs::interval::{sort_by_coord, sweep_intersect_pairs, IntersectOpts, StrandMode};
use trackcluster_rs::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

#[cfg(feature = "index-binned")]
use trackcluster_rs::interval::index_binned::{BinnedIntersectIndex, BinnedIntersectScratch};

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
        debug_assert!(low < high);
        let span = (high - low) as u64;
        low + (self.next_u64() % span) as u32
    }

    fn gen_bool(&mut self) -> bool {
        (self.next_u64() & 1) == 1
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

fn make_span_tx(chrom: &str, strand: Strand, name: String, start: u32, end: u32) -> Transcript {
    make_tx(chrom, strand, name, vec![(start, end)], "none", 0)
}

fn make_intersect_inputs(
    seed: u64,
    a_len: usize,
    b_len: usize,
) -> (Vec<Transcript>, Vec<Transcript>) {
    let mut rng = Lcg64::new(seed);
    let chroms = ["chr1", "chr2", "chr3", "chrX"];
    let mut a = Vec::with_capacity(a_len);
    let mut b = Vec::with_capacity(b_len);

    for i in 0..a_len {
        let chrom = chroms[(rng.next_u64() as usize) % chroms.len()];
        let strand = if rng.gen_bool() {
            Strand::Plus
        } else {
            Strand::Minus
        };
        let start = rng.gen_range_u32(0, 1_000_000);
        let len = rng.gen_range_u32(20, 500);
        let end = start.saturating_add(len);
        a.push(make_span_tx(chrom, strand, format!("a{i}"), start, end));
    }

    for i in 0..b_len {
        let chrom = chroms[(rng.next_u64() as usize) % chroms.len()];
        let strand = if rng.gen_bool() {
            Strand::Plus
        } else {
            Strand::Minus
        };
        let start = rng.gen_range_u32(0, 1_000_000);
        let len = rng.gen_range_u32(20, 500);
        let end = start.saturating_add(len);
        b.push(make_span_tx(chrom, strand, format!("b{i}"), start, end));
    }

    sort_by_coord(&mut a);
    sort_by_coord(&mut b);
    (a, b)
}

fn make_clusterj_inputs(
    seed: u64,
    refs_len: usize,
    reads_len: usize,
) -> (Vec<Transcript>, Vec<Transcript>) {
    let mut rng = Lcg64::new(seed);
    let mut refs = Vec::with_capacity(refs_len);
    let mut reads = Vec::with_capacity(reads_len);

    for i in 0..refs_len {
        let tx_start = rng.gen_range_u32(0, 200_000);
        let exons = exon_chain(tx_start, 50, 100, 4);
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
        let tx_start = rng.gen_range_u32(0, 200_000);
        let shift = rng.gen_range_u32(0, 5);
        let exons = exon_chain(tx_start + shift, 50, 100, 4);
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

fn make_clusterj_single_locus_inputs(
    seed: u64,
    refs_len: usize,
    reads_len: usize,
    locus_start: u32,
) -> (Vec<Transcript>, Vec<Transcript>) {
    let mut rng = Lcg64::new(seed);
    let mut refs = Vec::with_capacity(refs_len);
    let mut reads = Vec::with_capacity(reads_len);

    for i in 0..refs_len {
        let exon_len = 45 + rng.gen_range_u32(0, 10);
        let gap_len = 90 + rng.gen_range_u32(0, 30);
        let exons = exon_chain(locus_start, exon_len, gap_len, 4);
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
        let exon_len = 30 + rng.gen_range_u32(0, 40);
        let gap_len = 60 + rng.gen_range_u32(0, 80);
        let exon_count = if rng.gen_bool() { 3 } else { 4 };
        let exons = exon_chain(locus_start, exon_len, gap_len, exon_count);
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

fn make_clusterj_high_diversity_inputs(reads_len: usize, locus_start: u32) -> Vec<Transcript> {
    (0..reads_len)
        .map(|index| {
            let middle_start = locus_start + 100 + (index as u32 * 7);
            make_tx(
                "chr1",
                Strand::Plus,
                format!("diverse-read-{index}"),
                vec![
                    (locus_start, locus_start + 40),
                    (middle_start, middle_start + 30),
                    (locus_start + 100_000, locus_start + 100_040),
                ],
                "nanopore_read",
                0,
            )
        })
        .collect()
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

fn bench_interval_sweep_intersect(c: &mut Criterion) {
    let mut group = c.benchmark_group("interval");
    group.sample_size(10);

    let opts = IntersectOpts {
        strand_mode: StrandMode::Ignore,
        min_overlap_bp: None,
    };

    for (a_len, b_len) in [(5_000usize, 5_000usize), (20_000, 20_000)] {
        let (a, b) = make_intersect_inputs(1, a_len, b_len);
        group.bench_with_input(
            BenchmarkId::new("sweep_intersect_pairs", format!("{a_len}x{b_len}")),
            &(a, b),
            |bench, (a, b)| {
                bench.iter(|| {
                    let pairs = sweep_intersect_pairs(black_box(a), black_box(b), black_box(&opts));
                    black_box(pairs.len());
                });
            },
        );
    }

    group.finish();
}

fn bench_clusterj_grouping(c: &mut Criterion) {
    let mut group = c.benchmark_group("clusterj");
    group.sample_size(10);

    for (refs_len, reads_len) in [(200usize, 2_000usize), (500, 5_000)] {
        let (refs, reads) = make_clusterj_inputs(2, refs_len, reads_len);
        group.bench_with_input(
            BenchmarkId::new(
                "clusterj_default",
                format!("{refs_len}_refs_{reads_len}_reads"),
            ),
            &(refs, reads),
            |bench, (refs, reads)| {
                bench.iter(|| {
                    let result = clusterj::clusterj_with_name2_mode(
                        black_box(reads),
                        Some(black_box(refs)),
                        1,
                        clusterj::DEFAULT_SW_SCORE,
                        500,
                        100,
                        clusterj::Name2Mode::Coverage,
                    );
                    black_box(result.isoforms.len());
                });
            },
        );
    }

    let (refs, reads) = make_clusterj_inputs(22, 100, 1_000);
    group.bench_with_input(
        BenchmarkId::new("clusterj_sl_full_payload", "100_refs_1000_reads"),
        &(refs, reads),
        |bench, (refs, reads)| {
            bench.iter(|| {
                let result = clusterj::clusterj_with_name2_mode(
                    black_box(reads),
                    Some(black_box(refs)),
                    1,
                    11,
                    500,
                    100,
                    clusterj::Name2Mode::Full,
                );
                black_box(result.isoforms.len());
            });
        },
    );

    group.finish();
}

fn bench_clusterj_large_single_locus(c: &mut Criterion) {
    let mut group = c.benchmark_group("clusterj_large");
    group.sample_size(10);

    let refs_len = 200usize;
    let reads_len = 20_000usize;
    let (refs, reads) = make_clusterj_single_locus_inputs(4, refs_len, reads_len, 100_000);
    group.bench_with_input(
        BenchmarkId::new(
            "clusterj_simple_merge",
            format!("{refs_len}_refs_{reads_len}_reads"),
        ),
        &(refs, reads),
        |bench, (refs, reads)| {
            bench.iter(|| {
                let result = clusterj::clusterj_with_name2_mode(
                    black_box(reads),
                    Some(black_box(refs)),
                    1,
                    clusterj::DEFAULT_SW_SCORE,
                    0,
                    1,
                    clusterj::Name2Mode::Coverage,
                );
                black_box(result.isoforms.len());
            });
        },
    );

    group.finish();
}

fn bench_clusterj_high_diversity_single_locus(c: &mut Criterion) {
    let mut group = c.benchmark_group("clusterj_high_diversity");
    group.sample_size(10);

    let reads = make_clusterj_high_diversity_inputs(2_000, 100_000);
    let refs = vec![make_tx(
        "chr1",
        Strand::Plus,
        "diverse-reference".to_owned(),
        vec![(100_000, 100_040), (100_100, 100_130), (200_000, 200_040)],
        "isoform_anno",
        100,
    )];
    group.bench_with_input(
        BenchmarkId::new("mostly_nonmergeable_default", reads.len()),
        &(refs, reads),
        |bench, (refs, reads)| {
            bench.iter(|| {
                let result = clusterj::clusterj_with_name2_mode(
                    black_box(reads),
                    Some(black_box(refs)),
                    1,
                    clusterj::DEFAULT_SW_SCORE,
                    500,
                    100,
                    clusterj::Name2Mode::Coverage,
                );
                black_box(result.isoforms.len() + result.unused.len());
            });
        },
    );

    group.finish();
}

fn bench_cluster_overlap_synthetic_locus(c: &mut Criterion) {
    let mut group = c.benchmark_group("cluster_overlap");
    group.sample_size(10);

    for (refs_len, reads_len) in [(50usize, 400usize), (200, 2_000), (500, 5_000)] {
        let (refs, reads) = make_cluster_overlap_inputs(3, refs_len, reads_len, 10_000, 50_000);
        group.bench_with_input(
            BenchmarkId::new(
                "cluster_overlap",
                format!("{refs_len}_refs_{reads_len}_reads"),
            ),
            &(refs, reads),
            |bench, (refs, reads)| {
                bench.iter(|| {
                    let result =
                        cluster_overlap::cluster(black_box(reads), Some(black_box(refs)), 1);
                    black_box(result.isoforms.len());
                });
            },
        );
    }

    group.finish();
}

fn bench_cluster_overlap_batch_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("cluster_overlap_batch_size");
    group.sample_size(10);

    let refs_len = 500usize;
    let reads_len = 5_000usize;
    let (refs, reads) = make_cluster_overlap_inputs(5, refs_len, reads_len, 10_000, 50_000);

    for batch_size in [0usize, 250, 500, 1_000, 100_000] {
        group.bench_with_input(
            BenchmarkId::new("batch_size", batch_size),
            &batch_size,
            |bench, &batch_size| {
                bench.iter(|| {
                    let result = cluster_overlap::cluster_with_options(
                        black_box(&reads),
                        Some(black_box(&refs)),
                        1,
                        cluster_overlap::ClusterOptions {
                            batch_size,
                            name2_mode: clusterj::Name2Mode::Coverage,
                            ..cluster_overlap::ClusterOptions::default()
                        },
                    );
                    black_box(result.isoforms.len());
                });
            },
        );
    }

    group.finish();
}

#[cfg(feature = "index-binned")]
fn bench_interval_binned_reuse(c: &mut Criterion) {
    let mut group = c.benchmark_group("interval_reuse");
    group.sample_size(10);

    let opts = IntersectOpts {
        strand_mode: StrandMode::Ignore,
        min_overlap_bp: None,
    };

    let b_len = 50_000usize;
    let set_count = 50usize;
    let set_len = 200usize;

    let (_, b) = make_intersect_inputs(42, 0, b_len);
    let mut a_sets: Vec<Vec<Transcript>> = Vec::with_capacity(set_count);
    for i in 0..set_count {
        let (a, _) = make_intersect_inputs(1_000 + i as u64, set_len, 0);
        a_sets.push(a);
    }

    let index = BinnedIntersectIndex::build(&b, opts.strand_mode);
    let mut scratch = BinnedIntersectScratch::new(b.len());
    let mut out: Vec<(usize, usize)> = Vec::new();

    group.bench_function("binned_index_reuse", |bench| {
        bench.iter(|| {
            let mut total = 0usize;
            for a in &a_sets {
                index.intersect_pairs_into(
                    black_box(a),
                    black_box(&b),
                    black_box(&opts),
                    &mut scratch,
                    &mut out,
                );
                total = total.wrapping_add(out.len());
            }
            black_box(total);
        });
    });

    group.bench_function("sweep_repeated_queries", |bench| {
        bench.iter(|| {
            let mut total = 0usize;
            for a in &a_sets {
                let pairs = sweep_intersect_pairs(black_box(a), black_box(&b), black_box(&opts));
                total = total.wrapping_add(pairs.len());
            }
            black_box(total);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_interval_sweep_intersect,
    bench_clusterj_grouping,
    bench_clusterj_large_single_locus,
    bench_clusterj_high_diversity_single_locus,
    bench_cluster_overlap_synthetic_locus,
    bench_cluster_overlap_batch_sizes
);

#[cfg(feature = "index-binned")]
criterion_group!(benches_index_binned, bench_interval_binned_reuse);

#[cfg(feature = "index-binned")]
criterion_main!(benches, benches_index_binned);

#[cfg(not(feature = "index-binned"))]
criterion_main!(benches);
