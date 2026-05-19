//! Benchmark co-located with the `mp_mul2` vector test, as the plan
//! requires. `mp_mul2` is fixed two-digit, so there is no `nwords` to
//! vary; instead it is benched over four operand patterns as four group
//! ids. This matches the structure of the other `mp` benches (four
//! `nwords_*` ids): a lone single-id group at the host noise floor was
//! the most volatile of all and tripped the wall-clock gate on one run
//! of a pair while its sibling and the self-test stayed clean. Averaging
//! over several ids restores parity; the fixed `INNER`-fold repeat lifts
//! each above the floor (see mp_add.rs for the rationale).

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
    let m = u64::MAX;
    let cases: [(&str, [u64; 2], [u64; 2]); 4] = [
        (
            "random",
            [0x0123_4567_89ab_cdef, 0xfedc_ba98_7654_3210],
            [0x1111_2222_3333_4444, 0x5555_6666_7777_8888],
        ),
        ("max", [m, m], [m, m]),
        (
            "low_only",
            [0x9e37_79b9_7f4a_7c15, 0],
            [0xc2b2_ae3d_27d4_eb4f, 0],
        ),
        (
            "cross",
            [0, 0xa5a5_a5a5_a5a5_a5a5],
            [0x5a5a_5a5a_5a5a_5a5a, 0],
        ),
    ];
    let mut out = [0u64; 4];
    let mut group = c.benchmark_group("mp_mul2");
    group.throughput(Throughput::Elements(INNER as u64));
    for (name, a, b) in cases {
        group.bench_function(name, |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    mp_mul2(black_box(&mut out), black_box(&a), black_box(&b));
                }
                black_box(&out);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_mul2);
criterion_main!(benches);
