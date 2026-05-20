//! Top-level benchmark for [`sqisign_sign::protocols_keygen`].
//!
//! Seeds a deterministic [`CtrDrbg`] with a fixed entropy block and
//! measures wall-clock time for one keypair generation. Held high-level
//! per sir's divide-and-conquer policy: bench the public API first;
//! drill into primitives only if this gate trips.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sqisign_common::CtrDrbg;
use sqisign_sign::{protocols_keygen, SecretKey};
use sqisign_verify::PublicKey;

fn bench_keygen(c: &mut Criterion) {
    let entropy = [0x5au8; 48];
    let mut group = c.benchmark_group("sqisign_sign");
    // Each iteration constructs fresh DRBG, pk, sk; the keypair routine
    // dominates the iteration cost (~30 ms in release on a modest host).
    group.sample_size(20);
    group.bench_function("protocols_keygen_lvl1", |b| {
        b.iter(|| {
            let mut drbg = CtrDrbg::new(&entropy, None);
            let mut pk = PublicKey::zero();
            let mut sk = SecretKey::new();
            let ok = protocols_keygen(&mut drbg, &mut pk, &mut sk);
            black_box(ok);
            black_box(&pk);
            black_box(&sk);
        });
    });
    group.finish();
}

criterion_group!(benches, bench_keygen);
criterion_main!(benches);
