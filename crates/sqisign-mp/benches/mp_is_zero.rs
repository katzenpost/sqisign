//! Benchmark co-located with the `mp_is_zero` vector test, as the plan
//! requires. All-zero input is the worst case for the result-faithful
//! `iter().all` (it scans every limb), and value-independent for the
//! reference's constant-time reduction; the fixed `INNER`-fold repeat
//! lifts it above the host noise floor (see mp_add.rs).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_is_zero;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_is_zero(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_is_zero");
    for &n in &[4usize, 8, 16, 64] {
        let a = vec![0u64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                let mut acc = 0u64;
                for _ in 0..reps {
                    acc += mp_is_zero(black_box(&a)) as u64;
                }
                black_box(acc);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_is_zero);
criterion_main!(benches);
