//! Shared decode helpers for the id2iso differential vector tests.
//!
//! Mirrors the layout of `crates/sqisign-ec/tests/common/mod.rs` and
//! `crates/sqisign-quaternion/tests/common.rs`, lifted into one module
//! so the per-boundary test files stay minimal.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use std::collections::BTreeMap;

use sqisign_ec::{EcBasis, EcCurve, EcPoint};
use sqisign_gf::{Fp, Fp2, NWORDS_FIELD};
use sqisign_quaternion::dim2::ibz_mat_2x2_new;
use sqisign_quaternion::{
    ibz_vec_2_new, Ibz, IbzMat2x2, IbzVec2, QuatAlgElem, QuatLattice, QuatLeftIdeal,
};
use sqisign_vectors::decode;

pub fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be {NWORDS_FIELD} u64 limbs"
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

pub fn read_ibz(name: &str, fields: &BTreeMap<String, String>) -> Ibz {
    Ibz::from_canonical_bytes(&decode(name, &fields[name]).unwrap()).unwrap()
}

pub fn read_vec2(prefix: &str, fields: &BTreeMap<String, String>) -> IbzVec2 {
    let mut v = ibz_vec_2_new();
    v[0] = read_ibz(&format!("{prefix}_0"), fields);
    v[1] = read_ibz(&format!("{prefix}_1"), fields);
    v
}

pub fn read_mat2x2(prefix: &str, fields: &BTreeMap<String, String>) -> IbzMat2x2 {
    let mut m = ibz_mat_2x2_new();
    m[0][0] = read_ibz(&format!("{prefix}_00"), fields);
    m[0][1] = read_ibz(&format!("{prefix}_01"), fields);
    m[1][0] = read_ibz(&format!("{prefix}_10"), fields);
    m[1][1] = read_ibz(&format!("{prefix}_11"), fields);
    m
}

pub fn read_mat4x4(
    prefix: &str,
    fields: &BTreeMap<String, String>,
) -> sqisign_quaternion::IbzMat4x4 {
    let mut m = sqisign_quaternion::ibz_mat_4x4_new();
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = read_ibz(&format!("{prefix}_{i}_{j}"), fields);
        }
    }
    m
}

pub fn read_lattice(prefix: &str, fields: &BTreeMap<String, String>) -> QuatLattice {
    let denom = read_ibz(&format!("{prefix}_denom"), fields);
    let basis = read_mat4x4(&format!("{prefix}_basis"), fields);
    QuatLattice { denom, basis }
}

pub fn read_lideal(prefix: &str, fields: &BTreeMap<String, String>) -> QuatLeftIdeal {
    QuatLeftIdeal {
        norm: read_ibz(&format!("{prefix}_norm"), fields),
        lattice: read_lattice(&format!("{prefix}_lat"), fields),
    }
}

pub fn read_quat_elem(prefix: &str, fields: &BTreeMap<String, String>) -> QuatAlgElem {
    let mut e = QuatAlgElem::new();
    e.denom = read_ibz(&format!("{prefix}_denom"), fields);
    e.coord[0] = read_ibz(&format!("{prefix}_c0"), fields);
    e.coord[1] = read_ibz(&format!("{prefix}_c1"), fields);
    e.coord[2] = read_ibz(&format!("{prefix}_c2"), fields);
    e.coord[3] = read_ibz(&format!("{prefix}_c3"), fields);
    e
}

pub fn read_i32(key: &str, fields: &BTreeMap<String, String>) -> i32 {
    let b = decode(key, &fields[key]).expect("i32 hex");
    assert_eq!(b.len(), 4);
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

pub fn ibz_eq(a: &Ibz, b: &Ibz) -> bool {
    a.0 == b.0
}

pub fn mat2_eq(a: &IbzMat2x2, b: &IbzMat2x2) -> bool {
    (0..2).all(|i| (0..2).all(|j| ibz_eq(&a[i][j], &b[i][j])))
}

pub fn vec2_eq(a: &IbzVec2, b: &IbzVec2) -> bool {
    ibz_eq(&a[0], &b[0]) && ibz_eq(&a[1], &b[1])
}

pub fn mat4x4_eq(a: &sqisign_quaternion::IbzMat4x4, b: &sqisign_quaternion::IbzMat4x4) -> bool {
    (0..4).all(|i| (0..4).all(|j| ibz_eq(&a[i][j], &b[i][j])))
}

pub fn lattice_eq(a: &QuatLattice, b: &QuatLattice) -> bool {
    ibz_eq(&a.denom, &b.denom) && mat4x4_eq(&a.basis, &b.basis)
}

pub fn lideal_eq(a: &QuatLeftIdeal, b: &QuatLeftIdeal) -> bool {
    ibz_eq(&a.norm, &b.norm) && lattice_eq(&a.lattice, &b.lattice)
}
