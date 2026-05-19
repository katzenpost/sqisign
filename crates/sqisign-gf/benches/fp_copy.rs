//! Benchmark co-located with the `fp_copy` vector test, as the plan
//! requires. Times one plain five-limb copy on the level-1 field.
//!
//! `fp_copy` is the cheapest gf op (a single 40-byte memcpy with no
//! reduction); a single call is sub-nanosecond. Each sample therefore
//! performs a fixed `INNER`-fold repeat so the measured duration is
//! microseconds, well above the host noise floor. This does not weaken
//! the gate: a real >25% regression still shows as >25%, and the
//! slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_neg.rs`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_copy;

/// Fixed, code-independent repeat per timed sample.
const INNER: usize = 8192;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_copy(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_copy");
    let a = [
        0x0123_4567_89ab_cdefu64,
        0xfedc_ba98_7654_3210,
        0xdead_beef_cafe_babe,
        0x1357_9bdf_2468_ace0,
        0x0000_5000_0000_0000,
    ];
    let mut out = [0u64; 5];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_copy(black_box(&mut out), black_box(&a));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_copy);
criterion_main!(benches);
