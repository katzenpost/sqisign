//! Benchmark co-located with the `mp_mod_2exp` vector test, as the plan
//! requires. Fixed `INNER`-fold repeat lifts this short masking op above
//! the host noise floor (see mp_add.rs for the rationale); the op is
//! value-independent in cost so no reset is needed.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_mod_2exp;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_mod_2exp(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_mod_2exp");
    for &n in &[4usize, 8, 16, 64] {
        let mut a = vec![0xffff_ffff_ffff_ffffu64; n];
        // A mid-array reduction point: keeps the mask + zero-fill loop.
        let e = (n as u32 * 64) / 2 + 7;
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_mod_2exp(black_box(&mut a), black_box(e));
                }
                black_box(&a);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_mod_2exp);
criterion_main!(benches);
