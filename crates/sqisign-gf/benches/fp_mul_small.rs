//! Benchmark co-located with the `fp_mul_small` vector test, as the
//! plan requires. Times one GF(p) Montgomery multiplication by a small
//! integer on the level-1 five-limb field.
//!
//! `fp_mul_small` is the thin wrapper around `modmli`, which builds the
//! Montgomery representative of the integer scalar through `modint`
//! (one `modmul` against `NRES_C`) and then folds that into a second
//! `modmul` against the field input. Two `modmul` calls per invocation,
//! so the per-op cost is roughly twice a plain `fp_mul`; the per-sample
//! `INNER` fold is kept identical to the other gf binops so the
//! slowdown self-test trips with the same SLOWDOWN multiplier as the
//! rest of the gf battery and the bench reads at microsecond scale.
//! Mirrors `sqisign-gf/benches/fp_mul.rs`.
//!
//! The `a` operand is a near-2p pattern so every column accumulator is
//! wide and every `v_k * p4` Montgomery fold contributes a near-maximum
//! term, exercising the heaviest path through the inner `modmul`. The
//! `val` is rotated through three representative widths per iteration
//! (a small positive int, a value above `2^31 - 1` forcing the
//! sign-extending narrowing branch, and a value with non-zero high
//! bits forcing the high-bits-ignored narrowing branch) to keep the
//! branch predictor honest without perturbing the timing scale.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_gf::fp_mul_small;

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

fn bench_fp_mul_small(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("fp_mul_small");
    // Near-2p operand: every column accumulator is wide and every
    // v_k * p4 Montgomery fold contributes a near-maximum term,
    // exercising the heaviest path through the inner modmul.
    let a = [
        0x0007_ffff_ffff_fffeu64,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0007_ffff_ffff_ffff,
        0x0000_4fff_ffff_ffff,
    ];
    // Three representative vals: a small positive int (the intended
    // domain), a value above 2^31 - 1 (forces the sign-extending
    // narrowing branch), and a value with bit 31 set as a u32 (forces
    // the same sign-extending narrowing). Rotating through them per
    // iteration keeps the branch predictor honest without perturbing
    // the timing scale, which is dominated by the two modmul calls
    // inside modmli.
    let vals: [u32; 3] = [7, 0x8000_0000, 0xffff_ffff];
    let mut out = [0u64; 5];
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("level1", |bch| {
        bch.iter(|| {
            for i in 0..reps {
                fp_mul_small(
                    black_box(&mut out),
                    black_box(&a),
                    black_box(vals[i % vals.len()]),
                );
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fp_mul_small);
criterion_main!(benches);
