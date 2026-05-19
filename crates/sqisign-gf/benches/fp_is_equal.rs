//! Benchmark co-located with the `fp_is_equal` vector test, as the plan
//! requires. Times one GF(p) equality predicate on the level-1 five-limb
//! field.
//!
//! A single `fp_is_equal` (two `redc`s = two `modmul` + two `modfsb`,
//! plus the five-limb per-limb XOR / zero-detect / AND-fold) is tens of
//! nanoseconds, close to this host's wall-clock noise floor: gating one
//! op at a time false-failed on pure host jitter. Each sample therefore
//! performs a fixed `INNER`-fold repeat so the measured duration is
//! microseconds, well above the floor. This does not weaken the gate: a
//! real >10% regression still shows as >10%, and the slowdown self-test
//! multiplies the same repeat so it still trips. Mirrors
//! `sqisign-gf/benches/fp_is_zero`.
//!
//! Two distinct near-2p operands so each `redc`'s `modmul` does
//! substantive work on every column (an all-zero operand or an
//! identical-operand pair would skip partial-product contributions or
//! short-circuit the per-limb zero detect at the first limb, masking
//! any modmul/modcmp regression). The two operands are also picked to
//! be unequal under canonical reduction, so the AND-fold runs over a
//! mix of zero and nonzero per-limb XORs.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_is_equal;

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

fn bench_fp_is_equal(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_is_equal");
    // Two near-2p operands. `a` mirrors the fp_is_zero bench's operand
    // shape so the two predicate boundaries are timed on
    // comparable-magnitude limbs; `b` perturbs the low limb so the two
    // are unequal after canonical reduction. The accumulator forces
    // the result to be read so the optimiser cannot elide the call.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0,
        0,
        0x0000_5000_0000_0000,
    ];
    let b = [
        0x0007_ffff_ffff_fffdu64,
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
                acc = acc.wrapping_add(fp_is_equal(black_box(&a), black_box(&b)));
            }
            black_box(acc);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_is_equal);
criterion_main!(benches);
