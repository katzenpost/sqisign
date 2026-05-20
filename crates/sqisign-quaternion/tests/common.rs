//! Shared helpers for the closeout differential tests.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use std::collections::BTreeMap;

use sqisign_quaternion::{ibz_mat_4x4_new, ibz_vec_4_new, Ibz, IbzMat4x4, IbzVec4, QuatLattice};
use sqisign_vectors::decode;

pub fn read_ibz(name: &str, inputs: &BTreeMap<String, String>) -> Ibz {
    Ibz::from_canonical_bytes(&decode(name, &inputs[name]).unwrap()).unwrap()
}

pub fn read_vec4(prefix: &str, inputs: &BTreeMap<String, String>) -> IbzVec4 {
    let mut v = ibz_vec_4_new();
    for i in 0..4 {
        let k = format!("{prefix}_{i}");
        v[i] = read_ibz(&k, inputs);
    }
    v
}

pub fn read_mat4x4(prefix: &str, inputs: &BTreeMap<String, String>) -> IbzMat4x4 {
    let mut m = ibz_mat_4x4_new();
    for i in 0..4 {
        for j in 0..4 {
            let k = format!("{prefix}_{i}_{j}");
            m[i][j] = read_ibz(&k, inputs);
        }
    }
    m
}

pub fn read_lattice(prefix: &str, inputs: &BTreeMap<String, String>) -> QuatLattice {
    let denom_key = format!("{prefix}_denom");
    let denom = read_ibz(&denom_key, inputs);
    let basis = read_mat4x4(&format!("{prefix}_basis"), inputs);
    QuatLattice { denom, basis }
}

pub fn ibz_eq(a: &Ibz, b: &Ibz) -> bool {
    a.0 == b.0
}

pub fn vec4_eq(a: &IbzVec4, b: &IbzVec4) -> bool {
    (0..4).all(|i| ibz_eq(&a[i], &b[i]))
}

pub fn mat4x4_eq(a: &IbzMat4x4, b: &IbzMat4x4) -> bool {
    (0..4).all(|i| (0..4).all(|j| ibz_eq(&a[i][j], &b[i][j])))
}

pub fn read_i32(name: &str, inputs: &BTreeMap<String, String>) -> i32 {
    let b = decode(name, &inputs[name]).unwrap();
    assert!(b.len() == 4, "expected 4-byte i32 for {name}");
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

/// Re-export `ibz_mat_4x4_new` under a local name so tests can build
/// scratch matrices without re-importing it everywhere.
pub fn ibz_mat_4x4_new_local() -> IbzMat4x4 {
    ibz_mat_4x4_new()
}
