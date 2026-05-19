//! Benchmark co-located with the `mp_invert_matrix` vector test, as the
//! plan requires. It is the heaviest mp composition: a determinant, one
//! full-width `mp_inv_2e`, four `mp_mul`s and the reductions. The fixed
//! `INNER`-fold repeat keeps it above the host noise floor (see
//! mp_add.rs for the rationale).

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqisign_mp::mp_invert_matrix;

const INNER: usize = 1024;

/// Regression-gate self-test hook; see sqisign-common/benches/shake256.rs.
fn extra_passes() -> usize {
    std::env::var("SQISIGN_BENCH_SLOWDOWN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn bench_mp_invert_matrix(c: &mut Criterion) {
    let reps = INNER * (1 + extra_passes());
    let mut group = c.benchmark_group("mp_invert_matrix");
    for &n in &[4usize, 8, 16, 64] {
        // A fixed invertible matrix: r1,s2 odd; r2,s1 even => det odd.
        // Re-inverting in place can drive the determinant even, so each
        // rep restores from these immutable sources (cheap copies, no
        // per-rep allocation).
        let mut sr1 = vec![0x0123_4567_89ab_cdefu64; n];
        sr1[0] |= 1;
        let mut sr2 = vec![0x0f1e_2d3c_4b5a_6978u64; n];
        sr2[0] &= !1;
        let mut ss1 = vec![0x1122_3344_5566_7788u64; n];
        ss1[0] &= !1;
        let mut ss2 = vec![0xfedc_ba98_7654_3210u64; n];
        ss2[0] |= 1;
        let e = (64 * n) as i32; // full-width, the heaviest case
        let (mut r1, mut r2) = (sr1.clone(), sr2.clone());
        let (mut s1, mut s2) = (ss1.clone(), ss2.clone());
        group.throughput(Throughput::Elements((n * INNER) as u64));
        group.bench_function(format!("nwords_{n}"), |bch| {
            bch.iter(|| {
                for _ in 0..reps {
                    r1.copy_from_slice(&sr1);
                    r2.copy_from_slice(&sr2);
                    s1.copy_from_slice(&ss1);
                    s2.copy_from_slice(&ss2);
                    mp_invert_matrix(
                        black_box(&mut r1),
                        black_box(&mut r2),
                        black_box(&mut s1),
                        black_box(&mut s2),
                        black_box(e),
                    );
                }
                black_box((&r1, &r2, &s1, &s2));
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mp_invert_matrix);
criterion_main!(benches);
