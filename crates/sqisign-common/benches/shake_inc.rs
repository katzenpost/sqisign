//! Benchmark co-located with the incremental SHAKE vector test, as the plan
//! requires. It times the realistic incremental usage pattern: absorb the
//! input in fixed-size chunks, finalize, then squeeze 64 bytes. The chunked
//! shape is what SQIsign's callers actually do; a one-shot bench would not
//! exercise the per-absorb-call overhead the gate is meant to track.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_common::{Shake128Absorb, Shake256Absorb};

/// Regression-gate self-test hook; see `benches/shake256.rs`.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn s256(input: &[u8], chunk: usize) -> [u8; 64] {
    let mut a = Shake256Absorb::new();
    for c in input.chunks(chunk) {
        a.absorb(c);
    }
    let mut sq = a.finalize();
    let mut out = [0u8; 64];
    sq.squeeze(&mut out);
    out
}

fn s128(input: &[u8], chunk: usize) -> [u8; 64] {
    let mut a = Shake128Absorb::new();
    for c in input.chunks(chunk) {
        a.absorb(c);
    }
    let mut sq = a.finalize();
    let mut out = [0u8; 64];
    sq.squeeze(&mut out);
    out
}

fn bench_shake_inc(c: &mut Criterion) {
    let slowdown = extra_passes();
    let mut group = c.benchmark_group("shake_inc");
    for &len in &[32usize, 256, 4096] {
        let input = vec![0xa5u8; len];
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("shake256_absorb32_{len}_squeeze_64"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(s256(black_box(&input), 32));
                }
                s256(black_box(&input), 32)
            })
        });
        group.bench_function(format!("shake128_absorb32_{len}_squeeze_64"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(s128(black_box(&input), 32));
                }
                s128(black_box(&input), 32)
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_shake_inc);
criterion_main!(benches);
