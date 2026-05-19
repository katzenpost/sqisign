//! Benchmark co-located with the CTR-DRBG vector test, as the plan
//! requires. Times a seed plus a single draw of varying size; the per-call
//! update and AES key schedule dominate the small sizes, the counter stream
//! the large ones.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_common::CtrDrbg;

/// Regression-gate self-test hook; see `benches/shake256.rs`.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn seed_and_draw(entropy: &[u8; 48], n: usize) -> Vec<u8> {
    let mut drbg = CtrDrbg::new(entropy, None);
    let mut out = vec![0u8; n];
    drbg.fill(&mut out);
    out
}

fn bench_ctr_drbg(c: &mut Criterion) {
    let slowdown = extra_passes();
    let entropy = [0xa5u8; 48];
    let mut group = c.benchmark_group("ctr_drbg");
    for &len in &[32usize, 256, 4096] {
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("seed_then_draw_{len}"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(seed_and_draw(black_box(&entropy), len));
                }
                seed_and_draw(black_box(&entropy), len)
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_ctr_drbg);
criterion_main!(benches);
