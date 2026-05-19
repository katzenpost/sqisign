//! Benchmark co-located with the `mp_compare` vector test, as the plan
//! requires. Worst case is full-length equality (the scan visits every
//! limb); the fixed `INNER`-fold repeat lifts it above the host noise
//! floor (see mp_add.rs for the rationale).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_compare;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_compare(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_compare");
    for &n in &[4usize, 8, 16, 64] {
        // Equal operands: the scan does not early-return, so this is the
        // worst case and value-independent in cost.
        let a = vec![0x0123_4567_89ab_cdefu64; n];
        let b = a.clone();
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                let mut acc = 0i64;
                for _ in 0..reps {
                    acc += mp_compare(black_box(&a), black_box(&b)) as i64;
                }
                black_box(acc);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_compare);
criterion_main!(benches);
