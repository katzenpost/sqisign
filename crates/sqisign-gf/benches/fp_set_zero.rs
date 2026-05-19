//! Benchmark co-located with the `fp_set_zero` vector test, as the plan
//! requires. Times one plain five-limb zero-fill on the level-1 field.
//!
//! `fp_set_zero` is the cheapest possible field op (a single 40-byte
//! zero store with no reduction); a single call is sub-nanosecond. Each
//! sample therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the host noise floor. This does
//! not weaken the gate: a real >25% regression still shows as >25%, and
//! the slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_copy.rs`.
//!
//! The destination is given a non-trivial striped pre-fill so a port
//! that quietly skipped some limbs would be timed honestly rather than
//! benefiting from a buffer already at zero.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_set_zero;

/// Fixed, code-independent repeat per timed sample.
const INNER: usize = 8192;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_set_zero(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_set_zero");
    let mut out = [
        0x0123_4567_89ab_cdefu64,
        0xfedc_ba98_7654_3210,
        0xdead_beef_cafe_babe,
        0x1357_9bdf_2468_ace0,
        0x0000_5000_0000_0000,
    ];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_set_zero(black_box(&mut out));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_set_zero);
criterion_main!(benches);
