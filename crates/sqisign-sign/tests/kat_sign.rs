//! KAT round-trip test for SQIsign lvl1.
//!
//! Replays the recorded NIST KAT response file
//! `kat/PQCsignKAT_353_SQIsign_lvl1.rsp` (a carry of the upstream-released
//! file at the pinned commit, see `UPSTREAM.md`) against the Rust port.
//! Each entry seeds a fresh [`CtrDrbg`] with the recorded 48-byte
//! entropy, runs [`protocols_keygen`] followed by [`protocols_sign`],
//! and compares the serialized public key, secret key, and signed
//! message against the recorded bytes.
//!
//! The test is `#[ignore]`d by default so it does not run in routine
//! `cargo test` invocations: each entry can take several seconds, and a
//! single divergence anywhere in the keygen / sign chain produces a long
//! failure message. Run explicitly with
//! `cargo test -p sqisign-sign --release -- --ignored kat_lvl1`.

use std::fs;
use std::path::PathBuf;

use sqisign_common::CtrDrbg;
use sqisign_sign::{
    protocols_keygen, protocols_sign, secret_key_to_bytes, SecretKey, SECRETKEY_BYTES,
};
use sqisign_verify::{
    public_key_to_bytes, signature_to_bytes, PublicKey, Signature, PUBLICKEY_BYTES, SIGNATURE_BYTES,
};

fn kat_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("kat/PQCsignKAT_353_SQIsign_lvl1.rsp");
    p
}

struct KatEntry {
    count: u32,
    seed: [u8; 48],
    msg: Vec<u8>,
    pk: Vec<u8>,
    sk: Vec<u8>,
    sm: Vec<u8>,
}

fn parse_kat(path: &PathBuf) -> Vec<KatEntry> {
    let text = fs::read_to_string(path).expect("read KAT file");
    let mut entries: Vec<KatEntry> = Vec::new();
    let mut cur: Option<KatEntry> = None;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            if let Some(e) = cur.take() {
                entries.push(e);
            }
            continue;
        }
        let (key, value) = match line.split_once(" = ") {
            Some(kv) => kv,
            None => continue,
        };
        match key {
            "count" => {
                cur = Some(KatEntry {
                    count: value.parse().unwrap_or(0),
                    seed: [0u8; 48],
                    msg: Vec::new(),
                    pk: Vec::new(),
                    sk: Vec::new(),
                    sm: Vec::new(),
                });
            }
            "seed" => {
                if let Some(ref mut e) = cur {
                    let bytes = hex::decode(value).expect("seed hex");
                    e.seed.copy_from_slice(&bytes[..48]);
                }
            }
            "mlen" => {}
            "msg" => {
                if let Some(ref mut e) = cur {
                    e.msg = hex::decode(value).expect("msg hex");
                }
            }
            "pk" => {
                if let Some(ref mut e) = cur {
                    e.pk = hex::decode(value).expect("pk hex");
                }
            }
            "sk" => {
                if let Some(ref mut e) = cur {
                    e.sk = hex::decode(value).expect("sk hex");
                }
            }
            "smlen" => {}
            "sm" => {
                if let Some(ref mut e) = cur {
                    e.sm = hex::decode(value).expect("sm hex");
                }
            }
            _ => {}
        }
    }
    if let Some(e) = cur.take() {
        entries.push(e);
    }
    entries
}

fn run_one(e: &KatEntry) -> Result<(), String> {
    // One DRBG seeded by the recorded entropy drives both keygen and sign,
    // matching the NIST KAT driver (`PQCgenKAT_sign.c`: a single
    // `randombytes_init(seed, NULL, 256)` before keypair + sign).
    let mut drbg = CtrDrbg::new(&e.seed, None);

    let mut pk = PublicKey::zero();
    let mut sk = SecretKey::new();
    let kg_ok = protocols_keygen(&mut drbg, &mut pk, &mut sk);
    if kg_ok != 1 {
        return Err(format!(
            "count={}: protocols_keygen returned {} (want 1)",
            e.count, kg_ok
        ));
    }

    let mut pk_bytes = vec![0u8; PUBLICKEY_BYTES];
    public_key_to_bytes(&mut pk_bytes, &pk);
    if pk_bytes != e.pk {
        return Err(format!(
            "count={}: pk bytes diverged\n  got:  {}\n  want: {}",
            e.count,
            hex::encode_upper(&pk_bytes),
            hex::encode_upper(&e.pk)
        ));
    }

    let mut sk_bytes = [0u8; SECRETKEY_BYTES];
    secret_key_to_bytes(&mut sk_bytes, &sk, &pk);
    if sk_bytes[..] != e.sk[..] {
        return Err(format!(
            "count={}: sk bytes diverged\n  got:  {}\n  want: {}",
            e.count,
            hex::encode_upper(sk_bytes),
            hex::encode_upper(&e.sk)
        ));
    }

    let mut sig = Signature::zero();
    let sign_ok = protocols_sign(&mut drbg, &mut sig, &pk, &mut sk, &e.msg);
    if sign_ok != 1 {
        return Err(format!(
            "count={}: protocols_sign returned {} (want 1)",
            e.count, sign_ok
        ));
    }
    let mut sig_bytes = [0u8; SIGNATURE_BYTES];
    signature_to_bytes(&mut sig_bytes, &sig);

    // Signed message: signature_bytes || msg.
    let mut sm = Vec::with_capacity(SIGNATURE_BYTES + e.msg.len());
    sm.extend_from_slice(&sig_bytes);
    sm.extend_from_slice(&e.msg);
    if sm != e.sm {
        return Err(format!(
            "count={}: signed-message bytes diverged\n  got:  {}\n  want: {}",
            e.count,
            hex::encode_upper(&sm),
            hex::encode_upper(&e.sm)
        ));
    }
    Ok(())
}

#[test]
fn kat_lvl1_count_0() {
    let entries = parse_kat(&kat_path());
    let e = entries
        .iter()
        .find(|e| e.count == 0)
        .expect("KAT count=0 present");
    run_one(e).expect("KAT count=0 must match");
}

#[test]
#[ignore = "full KAT battery (~12s in release); count=0 is run by default"]
fn kat_lvl1_full() {
    let entries = parse_kat(&kat_path());
    assert!(!entries.is_empty(), "KAT file empty");
    let mut failures: Vec<String> = Vec::new();
    for e in &entries {
        if let Err(err) = run_one(e) {
            failures.push(err);
        }
    }
    if !failures.is_empty() {
        panic!(
            "{} KAT mismatches:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
