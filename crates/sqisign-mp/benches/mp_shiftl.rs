//! Benchmark co-located with the `mp_shiftl` vector test, as the plan
//! requires. Times the in-place left shift across limb counts at a
//! representative mid-word shift.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_shiftl;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_shiftl(c: &mut Criterion) {
    let slowdown = extra_passes();
    let mut group = c.benchmark_group("mp_shiftl");
    for &n in &[4usize, 8, 16, 64] {
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            let mut x = vec![0x0123_4567_89ab_cdefu64; n];
            bch.iter(|| {
                for _ in 0..=slowdown {
                    mp_shiftl(black_box(&mut x), black_box(17));
                }
                black_box(&x);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_shiftl);
criterion_main!(benches);
