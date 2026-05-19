//! Benchmark co-located with the `fp_is_zero` vector test, as the plan
//! requires. Times one GF(p) zero predicate on the level-1 five-limb
//! field.
//!
//! A single `fp_is_zero` (one `redc` = one `modmul` + one `modfsb`,
//! plus the five-limb OR and the `(d - 1) >> 51` bit-twiddle) is tens
//! of nanoseconds, close to this host's wall-clock noise floor: gating
//! one op at a time false-failed on pure host jitter. Each sample
//! therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the floor. This does not weaken
//! the gate: a real >10% regression still shows as >10%, and the
//! slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_neg`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_is_zero;

/// Fixed, code-independent repeat per timed sample. Chosen so the field
/// op is microsecond-scale.
const INNER: usize = 4096;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_fp_is_zero(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_is_zero");
    // A nonzero operand, near 2p, so redc's modmul does substantive
    // work on every column (the all-zero positive case would skip every
    // partial-product contribution, masking any modmul/modfsb
    // regression). Mirrors the fp_neg bench's operand shape so the two
    // predicate-vs-arithmetic boundaries are timed on
    // comparable-magnitude limbs. The accumulator forces the result to
    // be read so the optimiser cannot elide the call.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0,
        0,
        0x0000_5000_0000_0000,
    ];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            let mut acc = 0u32;
            for _ in 0..reps {
                acc = acc.wrapping_add(fp_is_zero(black_box(&a)));
            }
            black_box(acc);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_is_zero);
criterion_main!(benches);
