//! Differential verification of the ported lvl1 precomputed constants
//! against the committed canonical-bytes vectors under
//! `vectors/precomp/`. For every constant the test loads the vector,
//! re-encodes the `LazyLock`-materialised Rust value to canonical bytes
//! (the very encoding `tools/cdump` emitted on the C side via
//! `mpz_export`), and asserts bit-equality with the recorded hex.
//!
//! The integrity gate the plan demands ("add a generator-script
//! verification test that re-derives them and bit-compares") is therefore
//! anchored at the same value level as every other quaternion vector in
//! this workspace: the C side and the Rust side meet at canonical bytes.
//! Drift on either side is caught by this test.

use sqisign_ec::{EcBasis, EcCurve, EcPoint};
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_precomp::{
    COM_DEGREE, CONJUGATING_ELEMENTS, CONNECTING_IDEALS, CURVES_WITH_ENDOMORPHISMS,
    EXTREMAL_ORDERS, QUATALG_PINFTY, QUAT_PRIME_COFACTOR, SEC_DEGREE, TORSION_PLUS_2POWER,
    TWO_TO_SECURITY_BITS,
};
use sqisign_quaternion::dim2::IbzMat2x2;
use sqisign_quaternion::{Ibz, QuatAlgElem, QuatLattice, QuatLeftIdeal};
use sqisign_vectors::{encode, load, Vector};

const VEC_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/precomp");

fn vec_path(name: &str) -> String {
    format!("{VEC_DIR}/{name}.json")
}

fn assert_ibz(rec: &Vector, name: &str, val: &Ibz) {
    let expected = rec
        .outputs
        .get(name)
        .unwrap_or_else(|| panic!("vector missing output {name}"));
    let got = encode(&val.to_canonical_bytes());
    assert_eq!(
        &got, expected,
        "ibz field {name} (vector id {}) diverged",
        rec.id
    );
}

fn assert_u32(rec: &Vector, name: &str, val: u32) {
    let expected = rec
        .outputs
        .get(name)
        .unwrap_or_else(|| panic!("vector missing output {name}"));
    let got = encode(&val.to_le_bytes());
    assert_eq!(&got, expected, "u32 field {name} diverged");
}

fn fp_to_bytes(fp: &[u64; NWORDS_FIELD]) -> Vec<u8> {
    let mut out = Vec::with_capacity(NWORDS_FIELD * 8);
    for limb in fp {
        out.extend_from_slice(&limb.to_le_bytes());
    }
    out
}

fn assert_fp(rec: &Vector, name: &str, val: &[u64; NWORDS_FIELD]) {
    let expected = rec
        .outputs
        .get(name)
        .unwrap_or_else(|| panic!("vector missing output {name}"));
    let got = encode(&fp_to_bytes(val));
    assert_eq!(&got, expected, "fp field {name} diverged");
}

fn assert_fp2(rec: &Vector, prefix: &str, val: &Fp2) {
    assert_fp(rec, &format!("{prefix}_re"), &val.re);
    assert_fp(rec, &format!("{prefix}_im"), &val.im);
}

fn assert_ec_point(rec: &Vector, prefix: &str, p: &EcPoint) {
    assert_fp2(rec, &format!("{prefix}_x"), &p.x);
    assert_fp2(rec, &format!("{prefix}_z"), &p.z);
}

fn assert_ec_curve(rec: &Vector, prefix: &str, c: &EcCurve) {
    assert_fp2(rec, &format!("{prefix}_A"), &c.A);
    assert_fp2(rec, &format!("{prefix}_C"), &c.C);
    assert_ec_point(rec, &format!("{prefix}_A24"), &c.A24);
    assert_u32(
        rec,
        &format!("{prefix}_is_A24"),
        if c.is_A24_computed_and_normalized {
            1
        } else {
            0
        },
    );
}

fn assert_ec_basis(rec: &Vector, prefix: &str, b: &EcBasis) {
    assert_ec_point(rec, &format!("{prefix}_P"), &b.P);
    assert_ec_point(rec, &format!("{prefix}_Q"), &b.Q);
    assert_ec_point(rec, &format!("{prefix}_PmQ"), &b.PmQ);
}

fn assert_ibz_mat_2x2(rec: &Vector, prefix: &str, m: &IbzMat2x2) {
    assert_ibz(rec, &format!("{prefix}_00"), &m[0][0]);
    assert_ibz(rec, &format!("{prefix}_01"), &m[0][1]);
    assert_ibz(rec, &format!("{prefix}_10"), &m[1][0]);
    assert_ibz(rec, &format!("{prefix}_11"), &m[1][1]);
}

fn assert_lattice(rec: &Vector, prefix: &str, l: &QuatLattice) {
    assert_ibz(rec, &format!("{prefix}_denom"), &l.denom);
    for i in 0..4 {
        for j in 0..4 {
            assert_ibz(rec, &format!("{prefix}_basis_{i}_{j}"), &l.basis[i][j]);
        }
    }
}

fn assert_alg_elem(rec: &Vector, denom_field: &str, coord_prefix: &str, e: &QuatAlgElem) {
    assert_ibz(rec, denom_field, &e.denom);
    for (i, c) in e.coord.iter().enumerate() {
        assert_ibz(rec, &format!("{coord_prefix}_{i}"), c);
    }
}

fn assert_lideal(rec: &Vector, prefix: &str, id: &QuatLeftIdeal) {
    assert_ibz(rec, &format!("{prefix}_norm"), &id.norm);
    assert_lattice(rec, &format!("{prefix}_lat"), &id.lattice);
}

fn single_record(boundary_short: &str, expected_boundary: &str) -> Vector {
    let f = load(vec_path(boundary_short)).expect("vector file must load");
    assert_eq!(f.boundary, expected_boundary);
    assert_eq!(f.vectors.len(), 1);
    f.vectors.into_iter().next().unwrap()
}

#[test]
fn torsion_constants_match_reference_vectors() {
    let r = single_record(
        "TWO_TO_SECURITY_BITS",
        "sqisign_precomp::TWO_TO_SECURITY_BITS",
    );
    assert_ibz(&r, "value", &TWO_TO_SECURITY_BITS);

    let r = single_record(
        "TORSION_PLUS_2POWER",
        "sqisign_precomp::TORSION_PLUS_2POWER",
    );
    assert_ibz(&r, "value", &TORSION_PLUS_2POWER);

    let r = single_record("SEC_DEGREE", "sqisign_precomp::SEC_DEGREE");
    assert_ibz(&r, "value", &SEC_DEGREE);

    let r = single_record("COM_DEGREE", "sqisign_precomp::COM_DEGREE");
    assert_ibz(&r, "value", &COM_DEGREE);
}

#[test]
fn quaternion_data_simple_match_reference_vectors() {
    let r = single_record(
        "QUAT_prime_cofactor",
        "sqisign_precomp::QUAT_prime_cofactor",
    );
    assert_ibz(&r, "value", &QUAT_PRIME_COFACTOR);

    let r = single_record("QUATALG_PINFTY", "sqisign_precomp::QUATALG_PINFTY");
    assert_ibz(&r, "p", &QUATALG_PINFTY.p);
}

#[test]
fn extremal_orders_match_reference_vectors() {
    let f = load(vec_path("EXTREMAL_ORDERS")).expect("vector file must load");
    assert_eq!(f.boundary, "sqisign_precomp::EXTREMAL_ORDERS");
    assert_eq!(f.vectors.len(), 7);
    for (i, rec) in f.vectors.iter().enumerate() {
        let eo = &EXTREMAL_ORDERS[i];
        assert_lattice(rec, "order", &eo.order);
        assert_alg_elem(rec, "z_denom", "z", &eo.z);
        assert_alg_elem(rec, "t_denom", "t", &eo.t);
        assert_u32(rec, "q", eo.q);
    }
}

#[test]
fn connecting_ideals_match_reference_vectors() {
    let f = load(vec_path("CONNECTING_IDEALS")).expect("vector file must load");
    assert_eq!(f.boundary, "sqisign_precomp::CONNECTING_IDEALS");
    assert_eq!(f.vectors.len(), 7);
    for (i, rec) in f.vectors.iter().enumerate() {
        assert_lideal(rec, "ideal", &CONNECTING_IDEALS[i]);
    }
}

#[test]
fn conjugating_elements_match_reference_vectors() {
    let f = load(vec_path("CONJUGATING_ELEMENTS")).expect("vector file must load");
    assert_eq!(f.boundary, "sqisign_precomp::CONJUGATING_ELEMENTS");
    assert_eq!(f.vectors.len(), 7);
    for (i, rec) in f.vectors.iter().enumerate() {
        assert_alg_elem(rec, "denom", "c", &CONJUGATING_ELEMENTS[i]);
    }
}

#[test]
fn curves_with_endomorphisms_match_reference_vectors() {
    let f = load(vec_path("CURVES_WITH_ENDOMORPHISMS")).expect("vector file must load");
    assert_eq!(f.boundary, "sqisign_precomp::CURVES_WITH_ENDOMORPHISMS");
    assert_eq!(f.vectors.len(), 7);
    for (i, rec) in f.vectors.iter().enumerate() {
        let c = &CURVES_WITH_ENDOMORPHISMS[i];
        assert_ec_curve(rec, "curve", &c.curve);
        assert_ec_basis(rec, "basis", &c.basis_even);
        assert_ibz_mat_2x2(rec, "action_i", &c.action_i);
        assert_ibz_mat_2x2(rec, "action_j", &c.action_j);
        assert_ibz_mat_2x2(rec, "action_k", &c.action_k);
        assert_ibz_mat_2x2(rec, "action_gen2", &c.action_gen2);
        assert_ibz_mat_2x2(rec, "action_gen3", &c.action_gen3);
        assert_ibz_mat_2x2(rec, "action_gen4", &c.action_gen4);
    }
}
