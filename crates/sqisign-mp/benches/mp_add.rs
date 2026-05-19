//! Benchmark co-located with the `mp_add` vector test, as the plan
//! requires. Times addition across limb counts spanning the field/order
//! sizes SQIsign actually uses.
//!
//! A single word-array add of a few limbs is a few nanoseconds, below
//! this host's wall-clock noise floor: gating one op at a time
//! false-failed on pure host jitter even under criterion's
//! baseline-differential. Each sample therefore performs a fixed
//! `INNER`-fold repeat so the measured duration is microseconds, well
//! above the floor. This does not weaken the gate: a real >10% regression
//! still shows as >10%, and the slowdown self-test multiplies the same
//! repeat so it still trips.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_add;

/// Fixed, code-independent repeat per timed sample. Chosen so the smallest
/// benched width is still microsecond-scale.
const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_add(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_add");
    for &n in &[4usize, 8, 16, 64] {
        let a = vec![0x0123_4567_89ab_cdefu64; n];
        let b = vec![0xfedc_ba98_7654_3210u64; n];
        let mut out = vec![0u64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_add(black_box(&mut out), black_box(&a), black_box(&b));
                }
                black_box(&out);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_add);
criterion_main!(benches);
