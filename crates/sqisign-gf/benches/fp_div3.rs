//! Benchmark co-located with the `fp_div3` vector test, as the plan
//! requires. Times one GF(p) Montgomery division-by-three on the
//! level-1 five-limb field.
//!
//! `fp_div3` is exactly one `modmul` call with the precomputed constant
//! `THREE_INV` as the first operand, so its cost is the cost of a single
//! Montgomery multiplication. A single `modmul` is a few tens of
//! nanoseconds, close to this host's wall-clock noise floor: gating one
//! op at a time would false-fail on pure host jitter. Each sample
//! therefore performs a fixed `INNER`-fold repeat so the measured
//! duration is microseconds, well above the floor. This does not weaken
//! the gate: a real >10% regression still shows as >10%, and the
//! slowdown self-test multiplies the same repeat so it still trips.
//! Mirrors `sqisign-gf/benches/fp_half` and `sqisign-gf/benches/fp_mul`.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_div3;

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

fn bench_fp_div3(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_div3");
    // A near-2p operand so every column accumulator in the underlying
    // modmul is wide and every v_k * p4 Montgomery fold contributes a
    // near-maximum term, exercising the heaviest path. Mirrors the
    // fp_half / fp_sqr / fp_mul bench operand shape so the four
    // boundaries are timed on comparable-magnitude limbs.
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
                fp_div3(black_box(&mut out), black_box(&a));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_div3);
criterion_main!(benches);
