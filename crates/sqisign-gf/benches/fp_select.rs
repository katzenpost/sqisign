//! Benchmark co-located with the `fp_select` vector test, as the plan
//! requires. Times one branchless conditional select on the level-1
//! field at the `ctl == 0xFFFFFFFF` endpoint, the path that actually
//! exercises the per-limb XOR blend.
//!
//! `fp_select` is among the cheapest fp ops (five XOR-AND-XOR limb
//! steps, no reduction); a single call is sub-nanosecond. Each sample
//! therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the host noise floor. This does
//! not weaken the gate: a real >25% regression still shows as >25%,
//! and the slowdown self-test multiplies the same repeat so it still
//! trips. Mirrors `sqisign-gf/benches/fp_copy.rs`.
//!
//! `ctl` is fixed at `0xFFFFFFFF` so the blend path is timed honestly
//! (`ctl == 0` would also blend, but the all-zero `cw` means the inner
//! `cw & (a0 ^ a1)` term is identically zero; with `cw == u64::MAX`
//! the term is `a0 ^ a1`, the genuine blend). The destination is given
//! a non-trivial striped pre-fill so a port that quietly skipped any
//! limb would be timed honestly rather than benefiting from a buffer
//! already close to the target value.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_select;

/// Fixed, code-independent repeat per timed sample.
const INNER: usize = 8192;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_select(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_select");
    let a0 = [
        0x0123_4567_89ab_cdefu64,
        0xfedc_ba98_7654_3210,
        0xdead_beef_cafe_babe,
        0x1357_9bdf_2468_ace0,
        0x0000_5000_0000_0000,
    ];
    let a1 = [
        0xfedc_ba98_7654_3210u64,
        0x0123_4567_89ab_cdef,
        0xcafe_babe_dead_beef,
        0x2468_ace0_1357_9bdf,
        0x0000_3000_0000_0000,
    ];
    let mut d = [
        0xa5a5_a5a5_a5a5_a5a5u64,
        0x5a5a_5a5a_5a5a_5a5a,
        0xa5a5_a5a5_a5a5_a5a5,
        0x5a5a_5a5a_5a5a_5a5a,
        0xa5a5_a5a5_a5a5_a5a5,
    ];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_select(
                    black_box(&mut d),
                    black_box(&a0),
                    black_box(&a1),
                    black_box(0xFFFFFFFFu32),
                );
            }
            black_box(&d);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_select);
criterion_main!(benches);
