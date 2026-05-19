//! Benchmark co-located with the `fp_neg` vector test, as the plan
//! requires. Times one GF(p) negation on the level-1 five-limb field.
//!
//! A single `modneg` (five limbwise negations plus two `prop` passes) is
//! a few nanoseconds, below this host's wall-clock noise floor: gating
//! one op at a time false-failed on pure host jitter. Each sample
//! therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the floor. This does not weaken
//! the gate: a real >10% regression still shows as >10%, and the
//! slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_sub`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_neg;

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

fn bench_fp_neg(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_neg");
    // A nonzero operand so the limbwise 0 - a[i] borrows and prop's
    // negative-carry path fires, taking the modneg 2p-correction branch.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0,
        0,
        0x0000_5000_0000_0000,
    ];
    let mut out = [0u64; 5];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_neg(black_box(&mut out), black_box(&a));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_neg);
criterion_main!(benches);
