//! Benchmark co-located with the `fp_sqr` vector test, as the plan
//! requires. Times one GF(p) Montgomery squaring on the level-1 five-limb
//! field.
//!
//! A single `modsqr` (15 distinct partial products plus the five inline
//! Montgomery `v_k * p4` folds, with five `tot *= 2;` doublings) is a few
//! tens of nanoseconds, close to this host's wall-clock noise floor:
//! gating one op at a time would false-fail on pure host jitter. Each
//! sample therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the floor. This does not weaken
//! the gate: a real >10% regression still shows as >10%, and the
//! slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_mul`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_sqr;

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

fn bench_fp_sqr(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_sqr");
    // A near-2p operand so every column accumulator is wide and every
    // v_k * p4 Montgomery fold contributes a near-maximum term,
    // exercising the heaviest path through modsqr. Mirrors the fp_mul
    // bench's operand shape so the two boundaries are timed on
    // comparable-magnitude limbs.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0000_4fff_ffff_ffff,
    ];
    let mut out = [0u64; 5];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_sqr(black_box(&mut out), black_box(&a));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_sqr);
criterion_main!(benches);
