//! Benchmark co-located with the `fp_mul` vector test, as the plan
//! requires. Times one GF(p) Montgomery multiplication on the level-1
//! five-limb field.
//!
//! A single `modmul` (the 25-product schoolbook with five inline
//! Montgomery `v_k * p4` folds) is a few tens of nanoseconds, still close
//! to this host's wall-clock noise floor: gating one op at a time would
//! false-fail on pure host jitter. Each sample therefore performs a fixed
//! `INNER`-fold repeat so the measured duration is microseconds, well
//! above the floor. This does not weaken the gate: a real >10% regression
//! still shows as >10%, and the slowdown self-test multiplies the same
//! repeat so it still trips. Mirrors `sqisign-gf/benches/fp_add`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_mul;

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

fn bench_fp_mul(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_mul");
    // Two near-2p operands so every column accumulator is wide and every
    // v_k * p4 Montgomery fold contributes a near-maximum term, exercising
    // the heaviest path through modmul.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0000_4fff_ffff_ffff,
    ];
    let b = [
        0x0007_ffff_ffff_fffdu64,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0000_5000_0000_0000,
    ];
    let mut out = [0u64; 5];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                fp_mul(black_box(&mut out), black_box(&a), black_box(&b));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_mul);
criterion_main!(benches);
