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

    let message = b"sqisign bench: deterministic 32-byte payload     ";

    let mut group = c.benchmark_group("sqisign_sign");
    group.sample_size(20);
    group.bench_function("protocols_sign_lvl1", |b| {
        b.iter(|| {
            // Fresh DRBG per iteration so the signer's byte stream is
            // identical each pass. Cloning the secret key keeps the
            // template untouched.
            let mut drbg = CtrDrbg::new(&entropy, None);
            // Burn the keygen's worth of DRBG draws so the signer is
            // seeded the way protocols_sign would see it after a real
            // keypair generation. (We skip the keypair work, but mirror
            // its byte budget by drawing the same prefix.) For the bench
            // we just reseed and accept that the byte stream differs
            // from a full keypair-then-sign run; what we measure is the
            // sign-call hot loop, not the DRBG advance.
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
