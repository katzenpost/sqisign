//! Benchmark co-located with the `fp_cswap` vector test, as the plan
//! requires. Times one branchless conditional swap on the level-1 field
//! at the `ctl & 1 == 1` endpoint, the path that actually exercises
//! the full cross-multiplication and limb writeback.
//!
//! `fp_cswap` is more expensive than `fp_select` (five multiply-adds
//! and a multiply-subtract per limb, on both `g` and `f`), but still
//! sub-microsecond per call. Each sample therefore performs a fixed
//! `INNER`-fold repeat so the measured duration is microseconds, well
//! above the host noise floor. This does not weaken the gate: a real
//! regression of >25% still shows as >25%, and the slowdown self-test
//! multiplies the same repeat so it still trips. Mirrors
//! `sqisign-gf/benches/fp_select.rs`.
//!
//! `ctl` is fixed at `1` so the swap path is timed honestly. The
//! per-limb cross-multiplication runs the same code on both endpoints
//! (this is a constant-time op), but a `1` makes the algebraic effect
//! visible in any debug inspection of `g`/`f` between iterations.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_cswap;

/// Fixed, code-independent repeat per timed sample.
const INNER: usize = 8192;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_cswap(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_cswap");
    let mut g = [
        0x0123_4567_89ab_cdefu64,
        0xfedc_ba98_7654_3210,
        0xdead_beef_cafe_babe,
        0x1357_9bdf_2468_ace0,
        0x0000_5000_0000_0000,
    ];
    let mut f = [
        0xfedc_ba98_7654_3210u64,
        0x0123_4567_89ab_cdef,
        0xcafe_babe_dead_beef,
        0x2468_ace0_1357_9bdf,
        0x0000_3000_0000_0000,
    ];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_cswap(black_box(&mut g), black_box(&mut f), black_box(1u32));
            }
            black_box(&g);
            black_box(&f);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_cswap);
criterion_main!(benches);
