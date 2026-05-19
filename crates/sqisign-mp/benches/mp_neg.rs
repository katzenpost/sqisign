//! Benchmark co-located with the `mp_neg` vector test, as the plan
//! requires. Fixed `INNER`-fold repeat lifts this short complement loop
//! above the host noise floor (see mp_add.rs for the rationale). Cost is
//! value-independent, so no reset is needed (negating twice returns the
//! input for a[0] != 0, which keeps the operand non-degenerate).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_neg;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_neg(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_neg");
    for &n in &[4usize, 8, 16, 64] {
        let mut a = vec![0x0123_4567_89ab_cdefu64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_neg(black_box(&mut a));
                }
                black_box(&a);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_neg);
criterion_main!(benches);
