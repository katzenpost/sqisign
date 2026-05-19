//! Benchmark co-located with the `mp_inv_2e` vector test, as the plan
//! requires. The Hensel iteration count grows with e, so this is the
//! heaviest mp op; the fixed `INNER`-fold repeat keeps it comfortably
//! above the host noise floor (see mp_add.rs for the rationale).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_inv_2e;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_inv_2e(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_inv_2e");
    for &n in &[4usize, 8, 16, 64] {
        let mut a = vec![0x0123_4567_89ab_cdefu64; n];
        a[0] |= 1; // the reference requires odd a
        let e = (64 * n) as i32; // full-width inverse, the heaviest case
        let mut b = vec![0u64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_inv_2e(black_box(&mut b), black_box(&a), black_box(e));
                }
                black_box(&b);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_inv_2e);
criterion_main!(benches);
