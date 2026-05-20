//! Norm-equation helpers (deterministic subset).
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/normeq.c`. The
//! `quat_represent_integer` and `quat_sampling_random_ideal_O0_given_norm`
//! entry points are **RNG-driven** and out of scope for this batch (no
//! differential boundary). We port only the deterministic helpers:
//! `quat_lattice_O0_set`, `quat_lattice_O0_set_extremal`,
//! `quat_order_elem_create`, and `quat_change_to_O0_basis`.

#![allow(non_snake_case)]

use crate::algebra::{quat_alg_add, quat_alg_mul, quat_alg_scalar, QuatAlg, QuatAlgElem};
use crate::dim4::IbzVec4;
use crate::ibz::{ibz_add, ibz_const_one, ibz_div, ibz_divides, ibz_set, ibz_sub, Ibz};
use crate::lattice::QuatLattice;

/// `quat_p_extremal_maximal_order_t`.
#[derive(Clone, Debug, Default)]
pub struct QuatPExtremalMaximalOrder {
    pub order: QuatLattice,
    pub z: QuatAlgElem,
    pub t: QuatAlgElem,
    pub q: u32,
}

impl QuatPExtremalMaximalOrder {
    pub fn new() -> Self {
        Self {
            order: QuatLattice::new(),
            z: QuatAlgElem::new(),
            t: QuatAlgElem::new(),
            q: 0,
        }
    }
}

/// `quat_lattice_O0_set(O0)`: the order `(1, i, (i+j)/2, (1+ij)/2)`.
pub fn quat_lattice_O0_set(o0: &mut QuatLattice) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut o0.basis[i][j], 0);
        }
    }
    ibz_set(&mut o0.denom, 2);
    ibz_set(&mut o0.basis[0][0], 2);
    ibz_set(&mut o0.basis[1][1], 2);
    ibz_set(&mut o0.basis[2][2], 1);
    ibz_set(&mut o0.basis[1][2], 1);
    ibz_set(&mut o0.basis[3][3], 1);
    ibz_set(&mut o0.basis[0][3], 1);
}

/// `quat_lattice_O0_set_extremal(O0)`.
pub fn quat_lattice_O0_set_extremal(o0: &mut QuatPExtremalMaximalOrder) {
    ibz_set(&mut o0.z.coord[1], 1);
    ibz_set(&mut o0.t.coord[2], 1);
    ibz_set(&mut o0.z.denom, 1);
    ibz_set(&mut o0.t.denom, 1);
    o0.q = 1;
    quat_lattice_O0_set(&mut o0.order);
}

/// `quat_order_elem_create(elem, order, coeffs, alg)`: build an algebra
/// element from its `(1, z, t, t*z)`-basis coefficients.
pub fn quat_order_elem_create(
    elem: &mut QuatAlgElem,
    order: &QuatPExtremalMaximalOrder,
    coeffs: &IbzVec4,
    bpoo: &QuatAlg,
) {
    let mut quat_temp = QuatAlgElem::new();

    quat_alg_scalar(elem, &coeffs[0], &ibz_const_one());

    quat_alg_scalar(&mut quat_temp, &coeffs[1], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.z, &tmp_clone, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);

    quat_alg_scalar(&mut quat_temp, &coeffs[2], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &tmp_clone, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);

    quat_alg_scalar(&mut quat_temp, &coeffs[3], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &tmp_clone, bpoo);
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &tmp_clone, &order.z, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);
}

/// `quat_change_to_O0_basis(vec, el)`.
pub fn quat_change_to_O0_basis(vec: &mut IbzVec4, el: &QuatAlgElem) {
    vec[2] = el.coord[2].clone();
    let v2_clone = vec[2].clone();
    ibz_add(&mut vec[2], &v2_clone, &v2_clone);
    vec[3] = el.coord[3].clone();
    let v3_clone = vec[3].clone();
    ibz_add(&mut vec[3], &v3_clone, &v3_clone);
    ibz_sub(&mut vec[0], &el.coord[0], &el.coord[3]);
    ibz_sub(&mut vec[1], &el.coord[1], &el.coord[2]);

    assert!(ibz_divides(&vec[0], &el.denom) != 0);
    assert!(ibz_divides(&vec[1], &el.denom) != 0);
    assert!(ibz_divides(&vec[2], &el.denom) != 0);
    assert!(ibz_divides(&vec[3], &el.denom) != 0);

    let mut tmp = Ibz::zero();
    for i in 0..4 {
        let saved = vec[i].clone();
        ibz_div(&mut vec[i], &mut tmp, &saved, &el.denom);
    }
}
