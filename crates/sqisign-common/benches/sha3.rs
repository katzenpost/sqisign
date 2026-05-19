//! Benchmark co-located with the SHA3 vector test, as the plan requires.
//! Times the three fixed-output digests over inputs straddling their
//! distinct Keccak rates; the gate tracks each independently.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_common::{sha3_256, sha3_384, sha3_512};

/// Regression-gate self-test hook; see `benches/shake256.rs`.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_sha3(c: &mut Criterion) {
    let slowdown = extra_passes();
    let mut group = c.benchmark_group("sha3");
    for &len in &[32usize, 256, 4096] {
        let input = vec![0xa5u8; len];
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_function(format!("sha3_256_{len}"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(sha3_256(black_box(&input)));
                }
                sha3_256(black_box(&input))
            })
        });
        group.bench_function(format!("sha3_384_{len}"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(sha3_384(black_box(&input)));
                }
                sha3_384(black_box(&input))
            })
        });
        group.bench_function(format!("sha3_512_{len}"), |b| {
            b.iter(|| {
                for _ in 0..slowdown {
                    black_box(sha3_512(black_box(&input)));
                }
                sha3_512(black_box(&input))
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_sha3);
criterion_main!(benches);
