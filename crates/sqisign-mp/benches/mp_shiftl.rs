//! Benchmark co-located with the `mp_shiftl` vector test, as the plan
//! requires. Fixed `INNER`-fold repeat lifts the op above the host's
//! wall-clock noise floor (see mp_add.rs for the rationale). `mp_shiftl`
//! does a fixed per-limb loop independent of the operand value, so
//! repeating the shift in place (the value saturating toward zero) costs
//! the same per call and needs no reset.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_shiftl;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_shiftl(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_shiftl");
    for &n in &[4usize, 8, 16, 64] {
        let mut x = vec![0x0123_4567_89ab_cdefu64; n];
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
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
