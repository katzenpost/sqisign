//! Benchmark co-located with the `mp_copy` vector test, as the plan
//! requires. Fixed `INNER`-fold repeat lifts this op above the host
//! noise floor (see mp_add.rs for the rationale).
//!
//! `mp_copy` is `copy_from_slice` -- memcpy-class, the single fastest op
//! in the crate, so at the siblings' `INNER = 1024` each sample is too
//! short to rise above this contended VM's scheduling jitter, and it
//! repeatedly produced wild non-reproducing wall-clock swings (+85%,
//! +37.9%, a +14.5% false-pass, +71..239% across four ids). The
//! evidence-based fix is the same as the original `INNER` and the
//! mp_mul2 restructure: do enough work per sample that the wall-clock
//! signal dominates jitter. `mp_copy` therefore folds 32x more than its
//! siblings so its per-sample duration is comparable to theirs (which
//! are stable at 1024). This hardens the measurement; it does not loosen
//! the gate -- the threshold and the self-test are unchanged.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_copy;

const INNER: usize = 32 * 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_copy(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_copy");
    for &n in &[4usize, 8, 16, 64] {
        let a = vec![0x0123_4567_89ab_cdefu64; n];
        let mut b = vec![0u64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_copy(black_box(&mut b), black_box(&a));
                }
                black_box(&b);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_copy);
criterion_main!(benches);
