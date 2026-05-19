//! Benchmark co-located with the `mp_mul2` vector test, as the plan
//! requires. Fixed two-digit operands, so a single width; the fixed
//! `INNER`-fold repeat lifts this few-instruction op above the host noise
//! floor (see mp_add.rs for the rationale).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_mul2;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_mul2(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let a = [0x0123_4567_89ab_cdefu64, 0xfedc_ba98_7654_3210u64];
    let b = [0x1111_2222_3333_4444u64, 0x5555_6666_7777_8888u64];
    let mut out = [0u64; 4];
    let mut group = c.benchmark_group("mp_mul2");
    group.throughput(Throughput::Elements(INNER as u64));
    group.bench_function("two_digit", |bch| {
        bch.iter(|| {
            for _ in 0..reps {
                mp_mul2(black_box(&mut out), black_box(&a), black_box(&b));
            }
            black_box(&out);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_mp_mul2);
criterion_main!(benches);
