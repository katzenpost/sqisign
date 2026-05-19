//! Benchmark co-located with the `secure_clear` vector test, as the plan
//! requires. Times the optimiser-resistant wipe across buffer sizes; the
//! `black_box` around the buffer is essential here, since the entire point
//! of the routine is that the write cannot be elided.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_common::secure_clear;

/// Regression-gate self-test hook; see `benches/shake256.rs`.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_secure_clear(c: &mut Criterion) {
    let slowdown = extra_passes();
    let mut group = c.benchmark_group("secure_clear");
    for &len in &[32usize, 256, 4096] {
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("clear_{len}"), |b| {
            let mut buf = vec![0xa5u8; len];
            b.iter(|| {
                for _ in 0..=slowdown {
                    secure_clear(black_box(&mut buf));
                }
                black_box(&buf);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_secure_clear);
criterion_main!(benches);
