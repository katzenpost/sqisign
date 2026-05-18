//! Phase 0 placeholder benchmark, co-located with the `shake256` vector
//! test as the plan requires (every test ships a benchmark in the same
//! commit).
//!
//! As of Phase 1 this times the ported `sqisign_common::shake256` with the
//! same inputs the Phase 0 placeholder used, so the recorded baseline
//! transfers across the handover and the regression gate stays meaningful.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_common::hash::shake256_vec as shake256;

/// Regression-gate self-test hook. `SQISIGN_BENCH_SLOWDOWN=N` makes every
/// iteration do N extra hash passes, deliberately regressing the benchmark
/// so CI can assert that `tools/bench-gate.sh` actually fails on a
/// regression rather than only being green. Unset (the normal case) it is a
/// no-op. This is gate plumbing, never a real measurement.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_shake256(c: &mut Criterion) {
    let slowdown = extra_passes();
    let mut group = c.benchmark_group("shake256");
    for &len in &[32usize, 256, 4096] {
        let input = vec![0xa5u8; len];
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("absorb_{len}_squeeze_64"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(shake256(black_box(&input), black_box(64)));
                }
                shake256(black_box(&input), black_box(64))
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_shake256);
criterion_main!(benches);
