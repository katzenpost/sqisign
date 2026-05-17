//! Phase 0 placeholder benchmark, co-located with the `shake256` vector
//! test as the plan requires (every test ships a benchmark in the same
//! commit).
//!
//! It currently times the `sha3` oracle, not ported code, because nothing
//! is ported yet. When `sqisign-common`'s SHAKE lands in Phase 1 this body
//! swaps to the ported implementation with the same inputs, so the recorded
//! baseline transfers and the regression gate stays meaningful across the
//! handover.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sha3::digest::{ExtendableOutput, Update, XofReader};

fn shake256(input: &[u8], out_len: usize) -> Vec<u8> {
    let mut hasher = sha3::Shake256::default();
    hasher.update(input);
    let mut out = vec![0u8; out_len];
    hasher.finalize_xof().read(&mut out);
    out
}

fn bench_shake256(c: &mut Criterion) {
    let mut group = c.benchmark_group("shake256");
    for &len in &[32usize, 256, 4096] {
        let input = vec![0xa5u8; len];
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("absorb_{len}_squeeze_64"), |b| {
            b.iter(|| shake256(black_box(&input), black_box(64)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_shake256);
criterion_main!(benches);
