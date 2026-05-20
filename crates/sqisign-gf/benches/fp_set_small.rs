//! Benchmark co-located with the `fp_set_small` vector test, as the
//! plan requires. Times one int-to-Montgomery setter on the level-1
//! field.
//!
//! `fp_set_small` threads through `modint` and `nres`, the latter being
//! a single `modmul` call with the precomputed [`NRES_C`] constant. A
//! single call is therefore on the order of a single `modmul` on
//! redundant five-limb operands (low tens of nanoseconds), already well
//! above the host noise floor on its own, but the per-sample `INNER`
//! fold is kept identical to the other gf setters
//! (`fp_set_zero`/`fp_set_one`) so the slowdown self-test trips with the
//! same SLOWDOWN multiplier as the rest of the gf battery and the bench
//! reads at microsecond scale. This does not weaken the gate: a real
//! 25%+ regression still shows as 25%+, and the slowdown self-test
//! multiplies the same repeat so it still trips. Mirrors
//! `sqisign-gf/benches/fp_set_zero.rs`.
//!
//! The destination is given a non-trivial striped pre-fill so a port
//! that quietly skipped some limbs would be timed honestly rather than
//! benefiting from a buffer already at the target value, and the `val`
//! is rotated through three representative widths per iteration to
//! exercise both the small-positive and the high-bits-ignored
//! narrowing paths.
//!
//! [`NRES_C`]: ../sqisign_gf/index.html

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_set_small;

/// Fixed, code-independent repeat per timed sample.
const INNER: usize = 8192;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_set_small(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_set_small");
    let mut out = [
        0x0123_4567_89ab_cdefu64,
        0xfedc_ba98_7654_3210,
        0xdead_beef_cafe_babe,
        0x1357_9bdf_2468_ace0,
        0x0000_5000_0000_0000,
    ];
    // Three representative vals: a small positive int (the intended
    // domain), a value above 2^31 - 1 (forces the sign-extending
    // narrowing branch), and a value with non-zero high bits (forces
    // the high-bits-ignored narrowing branch). Rotating through them
    // per iteration keeps the branch predictor honest without
    // perturbing the timing scale, which is dominated by the modmul
    // inside nres.
    let vals: [u64; 3] = [7, 0x80000000, 0xffffffff_00000000];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for i in 0..reps {
                fp_set_small(black_box(&mut out), black_box(vals[i % vals.len()]));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_set_small);
criterion_main!(benches);
