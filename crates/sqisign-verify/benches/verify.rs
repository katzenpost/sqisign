//! Top-level benchmark for [`sqisign_verify::protocols_verify`].
//!
//! Loads the first recorded NIST KAT response from
//! `the-sqisign/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp`, decodes its
//! public key and signature, and times one verification per iteration.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sqisign_verify::{protocols_verify, public_key_decode, signature_decode};

fn kat_count_0() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    // Hardcoded count=0 entry from PQCsignKAT_353_SQIsign_lvl1.rsp; the
    // verify bench is self-contained and does not depend on the file
    // being present at run-time. The same values appear verbatim in
    // crates/sqisign-sign/tests/kat_sign.rs and the published KAT.
    let pk_hex = "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B";
    let msg_hex = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
    // sm = signature || msg; the signature is the first 148 bytes.
    let sm_hex = "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
    let pk = hex::decode(pk_hex).unwrap();
    let msg = hex::decode(msg_hex).unwrap();
    let sm = hex::decode(sm_hex).unwrap();
    let sig = sm[..148].to_vec();
    (pk, msg, sig)
}

fn bench_verify(c: &mut Criterion) {
    let (pk_bytes, msg, sig_bytes) = kat_count_0();
    let pk = public_key_decode(&pk_bytes).expect("public_key_decode");
    let sig = signature_decode(&sig_bytes).expect("signature_decode");

    let mut group = c.benchmark_group("sqisign_verify");
    group.sample_size(30);
    group.bench_function("protocols_verify_lvl1", |b| {
        b.iter(|| {
            let ok = protocols_verify(&sig, &pk, black_box(&msg));
            black_box(ok);
        });
    });
    group.finish();
}

criterion_group!(benches, bench_verify);
criterion_main!(benches);
