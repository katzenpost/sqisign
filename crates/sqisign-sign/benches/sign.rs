//! Top-level benchmark for [`sqisign_sign::protocols_sign`].
//!
//! Generates a keypair once (outside the timed loop) and then measures
//! one sign call per iteration. The per-iteration DRBG is reseeded from
//! the same entropy block as the keygen so the byte stream the signer
//! sees is deterministic across runs.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sqisign_common::CtrDrbg;
use sqisign_sign::{protocols_keygen, protocols_sign, SecretKey};
use sqisign_verify::{PublicKey, Signature};

fn bench_sign(c: &mut Criterion) {
    let entropy = [0xa5u8; 48];
    let mut drbg = CtrDrbg::new(&entropy, None);
    let mut pk = PublicKey::zero();
    let mut sk_template = SecretKey::new();
    assert_eq!(protocols_keygen(&mut drbg, &mut pk, &mut sk_template), 1);

    // The DRBG now sits in the state protocols_sign would see immediately
    // after a real keypair generation. Snapshot it; each iteration clones
    // this template so the signer always consumes the byte stream it
    // expects. Re-seeding from scratch would feed protocols_sign the
    // keygen's byte budget, which trips the divisibility debug_assert
    // deep in the signing loop.
    let drbg_template = drbg.clone();

    let message = b"sqisign bench: deterministic 32-byte payload     ";

    let mut group = c.benchmark_group("sqisign_sign");
    group.sample_size(20);
    group.bench_function("protocols_sign_lvl1", |b| {
        b.iter(|| {
            let mut drbg = drbg_template.clone();
            let mut sk = sk_template.clone();
            let mut sig = Signature::zero();
            let ok = protocols_sign(&mut drbg, &mut sig, &pk, &mut sk, message);
            black_box(ok);
            black_box(&sig);
        });
    });
    group.finish();
}

criterion_group!(benches, bench_sign);
criterion_main!(benches);
