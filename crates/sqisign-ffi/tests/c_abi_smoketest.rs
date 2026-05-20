//! C ABI smoke test: drive `sqisign_lvl1_verify` from Rust using the
//! same KAT data that `sqisign-verify::tests::kat_verify` consumes, and
//! confirm the FFI returns identically to a direct in-Rust verify call.
//!
//! This is the FFI-side parallel of `crates/sqisign-verify/tests/kat_verify.rs`.
//! Verifying 100 entries through the full protocol is slow, so this test
//! samples 10 entries: enough to cover the common decode and verification
//! paths without inflating CI runtime. The verify crate's own test still
//! covers all 100.

use std::fs;
use std::path::PathBuf;

use sqisign_ffi::{
    sqisign_lvl1_verify, SQISIGN_LVL1_PUBLIC_KEY_BYTES, SQISIGN_LVL1_SIGNATURE_BYTES,
};

const SAMPLE_SIZE: usize = 10;

struct KatEntry {
    count: u64,
    msg: Vec<u8>,
    pk: Vec<u8>,
    sig: Vec<u8>,
}

fn kat_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("vendor")
        .join("the-sqisign")
        .join("KAT")
        .join("PQCsignKAT_353_SQIsign_lvl1.rsp")
}

fn parse_kat(path: &PathBuf) -> Vec<KatEntry> {
    let raw = fs::read_to_string(path).expect("KAT file must be readable");
    let mut entries = Vec::new();
    let mut count: Option<u64> = None;
    let mut mlen: Option<usize> = None;
    let mut msg: Option<Vec<u8>> = None;
    let mut pk: Option<Vec<u8>> = None;
    let mut smlen: Option<usize> = None;

    for line in raw.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((k, v)) = line.split_once('=') else {
            continue;
        };
        let k = k.trim();
        let v = v.trim();
        match k {
            "count" => count = Some(v.parse().expect("count int")),
            "mlen" => mlen = Some(v.parse().expect("mlen int")),
            "msg" => msg = Some(hex::decode(v).expect("msg hex")),
            "pk" => pk = Some(hex::decode(v).expect("pk hex")),
            "smlen" => smlen = Some(v.parse().expect("smlen int")),
            "sm" => {
                let sm = hex::decode(v).expect("sm hex");
                let count = count.expect("count present");
                let mlen = mlen.expect("mlen present");
                let smlen = smlen.expect("smlen present");
                let msg = msg.take().expect("msg present");
                let pk = pk.take().expect("pk present");
                assert_eq!(msg.len(), mlen, "msg/mlen mismatch in count {count}");
                assert_eq!(sm.len(), smlen, "sm/smlen mismatch in count {count}");
                assert_eq!(pk.len(), SQISIGN_LVL1_PUBLIC_KEY_BYTES, "pk length");
                assert!(
                    sm.len() >= SQISIGN_LVL1_SIGNATURE_BYTES,
                    "sm shorter than SIGNATURE_BYTES",
                );
                let sig_bytes = sm[..SQISIGN_LVL1_SIGNATURE_BYTES].to_vec();
                entries.push(KatEntry {
                    count,
                    msg,
                    pk,
                    sig: sig_bytes,
                });
            }
            _ => { /* sk, seed, etc. */ }
        }
    }
    entries
}

fn ffi_verify(sig: &[u8], pk: &[u8], msg: &[u8]) -> bool {
    let r = unsafe {
        sqisign_lvl1_verify(
            sig.as_ptr(),
            sig.len(),
            pk.as_ptr(),
            pk.len(),
            if msg.is_empty() {
                core::ptr::null()
            } else {
                msg.as_ptr()
            },
            msg.len(),
        )
    };
    match r {
        0 => false,
        1 => true,
        other => panic!("sqisign_lvl1_verify returned {other}, expected 0 or 1"),
    }
}

#[test]
fn ffi_verifies_sample_of_kat_entries() {
    let path = kat_path();
    let entries = parse_kat(&path);
    assert_eq!(entries.len(), 100, "lvl1 KAT must have 100 entries");

    let sample = &entries[..SAMPLE_SIZE];
    let mut passed = 0;
    let mut failed = Vec::new();
    for e in sample {
        if ffi_verify(&e.sig, &e.pk, &e.msg) {
            passed += 1;
        } else {
            failed.push(e.count);
        }
    }
    assert!(
        failed.is_empty(),
        "{} of {} FFI verifications failed: {:?}",
        failed.len(),
        SAMPLE_SIZE,
        failed
    );
    assert_eq!(passed, SAMPLE_SIZE);
}

#[test]
fn ffi_rejects_tampered_signature() {
    let path = kat_path();
    let entries = parse_kat(&path);

    // Tamper with one entry: flip a bit in the auxiliary curve's encoded
    // A coordinate. The verifier should reject.
    let e = &entries[0];
    let mut bad = e.sig.clone();
    bad[0] ^= 0x80;
    assert!(
        !ffi_verify(&bad, &e.pk, &e.msg),
        "tampered signature unexpectedly verified",
    );
}

#[test]
fn ffi_rejects_length_mismatch() {
    let path = kat_path();
    let entries = parse_kat(&path);
    let e = &entries[0];

    // One byte short.
    let short = &e.sig[..SQISIGN_LVL1_SIGNATURE_BYTES - 1];
    assert!(!ffi_verify(short, &e.pk, &e.msg));

    // One byte short on the public key side.
    let short_pk = &e.pk[..SQISIGN_LVL1_PUBLIC_KEY_BYTES - 1];
    assert!(!ffi_verify(&e.sig, short_pk, &e.msg));
}
