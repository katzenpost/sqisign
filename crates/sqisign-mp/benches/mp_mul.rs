//! Benchmark co-located with the `mp_mul` vector test, as the plan
//! requires. Schoolbook low-half multiply is O(n^2); the fixed
//! `INNER`-fold repeat lifts the smaller widths above the host noise
//! floor (see mp_add.rs for the rationale).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_mul;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_mul(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_mul");
    for &n in &[4usize, 8, 16, 64] {
        let a = vec![0x0123_4567_89ab_cdefu64; n];
        let b = vec![0xfedc_ba98_7654_3210u64; n];
        let mut out = vec![0u64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_mul(black_box(&mut out), black_box(&a), black_box(&b));
                }
                black_box(&out);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_mul);
criterion_main!(benches);
