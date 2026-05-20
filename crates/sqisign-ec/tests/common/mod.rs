//! Shared decode helpers for the ec differential vector tests.
//!
//! Each boundary records its inputs and outputs decomposed into fp2
//! `_re` / `_im` limb pairs (see the `put_*_named` helpers in
//! `tools/cdump/src/dump_main.c`); this module reassembles those into
//! the corresponding Rust structs.

#![allow(dead_code)]

use std::collections::BTreeMap;

use sqisign_ec::{AddComponents, EcBasis, EcCurve, EcKps2, EcKps4, EcPoint, JacPoint};
use sqisign_gf::{Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::decode;

pub fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be exactly {NWORDS_FIELD} u64 limbs"
    );
    let mut limbs = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        limbs[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    limbs
}

pub fn fp2_from(prefix: &str, fields: &BTreeMap<String, String>) -> Fp2 {
    let re_key = format!("{prefix}_re");
    let im_key = format!("{prefix}_im");
    let re_bytes = decode(&re_key, &fields[&re_key]).expect("re hex");
    let im_bytes = decode(&im_key, &fields[&im_key]).expect("im hex");
    Fp2 {
        re: fp_from(&re_key, &re_bytes),
        im: fp_from(&im_key, &im_bytes),
    }
}

pub fn ec_point_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcPoint {
    EcPoint {
        x: fp2_from(&format!("{prefix}_x"), fields),
        z: fp2_from(&format!("{prefix}_z"), fields),
    }
}

pub fn jac_point_from(prefix: &str, fields: &BTreeMap<String, String>) -> JacPoint {
    JacPoint {
        x: fp2_from(&format!("{prefix}_x"), fields),
        y: fp2_from(&format!("{prefix}_y"), fields),
        z: fp2_from(&format!("{prefix}_z"), fields),
    }
}

pub fn ec_curve_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcCurve {
    let is_a24_key = format!("{prefix}_is_A24");
    let is_a24_bytes = decode(&is_a24_key, &fields[&is_a24_key]).expect("is_A24 hex");
    assert_eq!(is_a24_bytes.len(), 4);
    let is_a24 = u32::from_le_bytes(is_a24_bytes.as_slice().try_into().unwrap()) != 0;
    EcCurve {
        A: fp2_from(&format!("{prefix}_A"), fields),
        C: fp2_from(&format!("{prefix}_C"), fields),
        A24: ec_point_from(&format!("{prefix}_A24"), fields),
        is_A24_computed_and_normalized: is_a24,
    }
}

pub fn ec_basis_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcBasis {
    EcBasis {
        P: ec_point_from(&format!("{prefix}_P"), fields),
        Q: ec_point_from(&format!("{prefix}_Q"), fields),
        PmQ: ec_point_from(&format!("{prefix}_PmQ"), fields),
    }
}

pub fn add_components_from(prefix: &str, fields: &BTreeMap<String, String>) -> AddComponents {
    AddComponents {
        u: fp2_from(&format!("{prefix}_u"), fields),
        v: fp2_from(&format!("{prefix}_v"), fields),
        w: fp2_from(&format!("{prefix}_w"), fields),
    }
}

pub fn ec_kps2_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcKps2 {
    EcKps2 {
        K: ec_point_from(&format!("{prefix}_K"), fields),
    }
}

pub fn ec_kps4_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcKps4 {
    EcKps4 {
        K: [
            ec_point_from(&format!("{prefix}_K0"), fields),
            ec_point_from(&format!("{prefix}_K1"), fields),
            ec_point_from(&format!("{prefix}_K2"), fields),
        ],
    }
}

pub fn u32_field(key: &str, fields: &BTreeMap<String, String>) -> u32 {
    let bytes = decode(key, &fields[key]).expect("u32 hex");
    assert_eq!(bytes.len(), 4);
    u32::from_le_bytes(bytes.as_slice().try_into().unwrap())
}

pub fn i32_field(key: &str, fields: &BTreeMap<String, String>) -> i32 {
    u32_field(key, fields) as i32
}

pub fn u64_field(key: &str, fields: &BTreeMap<String, String>) -> u64 {
    let bytes = decode(key, &fields[key]).expect("u64 hex");
    assert_eq!(bytes.len(), 8);
    u64::from_le_bytes(bytes.as_slice().try_into().unwrap())
}

pub fn digits_field(key: &str, n: usize, fields: &BTreeMap<String, String>) -> Vec<u64> {
    let bytes = decode(key, &fields[key]).expect("digits hex");
    assert_eq!(bytes.len(), n * 8);
    let mut out = vec![0u64; n];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        out[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    out
}
