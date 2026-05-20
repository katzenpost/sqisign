//! Quaternion lattices.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/lattice.c`.
//! A lattice is a 4x4 integer basis divided by a positive denominator.
//! All bases are assumed to be in HNF unless explicitly noted otherwise.

use crate::algebra::{quat_alg_coord_mul, QuatAlg, QuatAlgElem};
use crate::dim4::{
    ibz_mat_4x4_copy, ibz_mat_4x4_equal, ibz_mat_4x4_eval, ibz_mat_4x4_gcd,
    ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_new, ibz_mat_4x4_scalar_div,
    ibz_mat_4x4_scalar_mul, ibz_mat_4x4_transpose, ibz_vec_4_copy_ibz, ibz_vec_4_new,
    ibz_vec_4_scalar_div, ibz_vec_4_scalar_mul, IbzMat4x4, IbzVec4,
};
use crate::hnf::ibz_mat_4xn_hnf_mod_core;
use crate::ibz::{ibz_abs, ibz_cmp, ibz_div, ibz_gcd, ibz_is_zero, ibz_mul, ibz_neg, ibz_set, Ibz};

/// `quat_lattice_t`: a 4x4 integer basis with a non-zero denominator.
#[derive(Clone, Debug)]
pub struct QuatLattice {
    pub denom: Ibz,
    pub basis: IbzMat4x4,
}

impl QuatLattice {
    pub fn new() -> Self {
        let mut basis = ibz_mat_4x4_new();
        // C reference's quat_lattice_init relies on ibz_init zeroing.
        // Match that. denom defaults to zero too (the C code does so).
        let _ = &mut basis;
        Self {
            denom: Ibz::zero(),
            basis,
        }
    }
}

impl Default for QuatLattice {
    fn default() -> Self {
        Self::new()
    }
}

/// `quat_lattice_reduce_denom(reduced, lat)`: divide basis and denom by
/// their joint gcd, force `denom > 0`.
pub fn quat_lattice_reduce_denom(reduced: &mut QuatLattice, lat: &QuatLattice) {
    let mut gcd = Ibz::zero();
    ibz_mat_4x4_gcd(&mut gcd, &lat.basis);
    let mut t = Ibz::zero();
    ibz_gcd(&mut t, &gcd, &lat.denom);
    gcd = t;
    let _ = ibz_mat_4x4_scalar_div(&mut reduced.basis, &gcd, &lat.basis);
    let mut q = Ibz::zero();
    let mut rem = Ibz::zero();
    ibz_div(&mut q, &mut rem, &lat.denom, &gcd);
    reduced.denom = q;
    let mut a = Ibz::zero();
    ibz_abs(&mut a, &reduced.denom);
    reduced.denom = a;
}

/// `quat_lattice_conjugate_without_hnf(conj, lat)`: negate the `i`, `j`,
/// `ij` rows of the basis. The result is NOT in HNF.
pub fn quat_lattice_conjugate_without_hnf(conj: &mut QuatLattice, lat: &QuatLattice) {
    ibz_mat_4x4_copy(&mut conj.basis, &lat.basis);
    conj.denom = lat.denom.clone();
    for row in 1..4 {
        for col in 0..4 {
            let mut t = Ibz::zero();
            ibz_neg(&mut t, &conj.basis[row][col]);
            conj.basis[row][col] = t;
        }
    }
}

/// `quat_lattice_dual_without_hnf(dual, lat)`: dual lattice; result NOT
/// in HNF. Caller must apply `quat_lattice_hnf` before use.
pub fn quat_lattice_dual_without_hnf(dual: &mut QuatLattice, lat: &QuatLattice) {
    let mut inv = ibz_mat_4x4_new();
    let mut det = Ibz::zero();
    let _ = ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), Some(&mut det), &lat.basis);
    let inv_copy = clone_mat(&inv);
    ibz_mat_4x4_transpose(&mut inv, &inv_copy);
    ibz_mat_4x4_scalar_mul(&mut dual.basis, &lat.denom, &inv);
    dual.denom = det;
}

/// `quat_lattice_add(res, lat1, lat2)`.
pub fn quat_lattice_add(res: &mut QuatLattice, lat1: &QuatLattice, lat2: &QuatLattice) {
    let mut generators: Vec<IbzVec4> = Vec::with_capacity(8);
    for _ in 0..8 {
        generators.push(ibz_vec_4_new());
    }
    let mut tmp = ibz_mat_4x4_new();
    let mut det1 = Ibz::zero();
    let mut det2 = Ibz::zero();
    let mut detprod = Ibz::zero();

    ibz_mat_4x4_scalar_mul(&mut tmp, &lat1.denom, &lat2.basis);
    for i in 0..4 {
        for j in 0..4 {
            generators[j][i] = tmp[i][j].clone();
        }
    }
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det1), &tmp);

    let mut tmp2 = ibz_mat_4x4_new();
    ibz_mat_4x4_scalar_mul(&mut tmp2, &lat2.denom, &lat1.basis);
    for i in 0..4 {
        for j in 0..4 {
            generators[4 + j][i] = tmp2[i][j].clone();
        }
    }
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det2), &tmp2);

    assert!(ibz_is_zero(&det1) == 0);
    assert!(ibz_is_zero(&det2) == 0);
    ibz_gcd(&mut detprod, &det1, &det2);
    ibz_mat_4xn_hnf_mod_core(&mut res.basis, 8, &generators, &detprod);
    ibz_mul(&mut res.denom, &lat1.denom, &lat2.denom);
    let cloned = clone_lat(res);
    quat_lattice_reduce_denom(res, &cloned);
}

/// `quat_lattice_intersect(res, lat1, lat2)`.
pub fn quat_lattice_intersect(res: &mut QuatLattice, lat1: &QuatLattice, lat2: &QuatLattice) {
    let mut dual1 = QuatLattice::new();
    let mut dual2 = QuatLattice::new();
    let mut dual_res = QuatLattice::new();
    quat_lattice_dual_without_hnf(&mut dual1, lat1);
    quat_lattice_dual_without_hnf(&mut dual2, lat2);
    quat_lattice_add(&mut dual_res, &dual1, &dual2);
    quat_lattice_dual_without_hnf(res, &dual_res);
    quat_lattice_hnf(res);
}

/// `quat_lattice_hnf(lat)`: in-place HNF normalisation.
pub fn quat_lattice_hnf(lat: &mut QuatLattice) {
    let mut mod_ = Ibz::zero();
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut mod_), &lat.basis);
    let mut a = Ibz::zero();
    ibz_abs(&mut a, &mod_);
    mod_ = a;

    let mut generators: Vec<IbzVec4> = (0..4).map(|_| ibz_vec_4_new()).collect();
    for i in 0..4 {
        for j in 0..4 {
            generators[j][i] = lat.basis[i][j].clone();
        }
    }
    ibz_mat_4xn_hnf_mod_core(&mut lat.basis, 4, &generators, &mod_);
    let cloned = clone_lat(lat);
    quat_lattice_reduce_denom(lat, &cloned);
}

/// `quat_lattice_equal(lat1, lat2)`.
pub fn quat_lattice_equal(lat1: &QuatLattice, lat2: &QuatLattice) -> i32 {
    let mut a = clone_lat(lat1);
    let mut b = clone_lat(lat2);
    let a_clone = clone_lat(&a);
    quat_lattice_reduce_denom(&mut a, &a_clone);
    let b_clone = clone_lat(&b);
    quat_lattice_reduce_denom(&mut b, &b_clone);
    let mut ad = Ibz::zero();
    ibz_abs(&mut ad, &a.denom);
    a.denom = ad;
    let mut bd = Ibz::zero();
    ibz_abs(&mut bd, &b.denom);
    b.denom = bd;
    quat_lattice_hnf(&mut a);
    quat_lattice_hnf(&mut b);
    let denom_eq = ibz_cmp(&a.denom, &b.denom) == 0;
    let basis_eq = ibz_mat_4x4_equal(&a.basis, &b.basis) != 0;
    if denom_eq && basis_eq {
        1
    } else {
        0
    }
}

/// `quat_lattice_inclusion(sublat, overlat)`: 1 if `sublat <= overlat`.
pub fn quat_lattice_inclusion(sublat: &QuatLattice, overlat: &QuatLattice) -> i32 {
    let mut sum = QuatLattice::new();
    quat_lattice_add(&mut sum, overlat, sublat);
    quat_lattice_equal(&sum, overlat)
}

/// `quat_lattice_mat_alg_coord_mul_without_hnf(prod, lat, coord, alg)`.
pub fn quat_lattice_mat_alg_coord_mul_without_hnf(
    prod: &mut IbzMat4x4,
    lat: &IbzMat4x4,
    coord: &IbzVec4,
    alg: &QuatAlg,
) {
    let mut p = ibz_vec_4_new();
    let mut a = ibz_vec_4_new();
    for i in 0..4 {
        ibz_vec_4_copy_ibz(&mut a, &lat[0][i], &lat[1][i], &lat[2][i], &lat[3][i]);
        quat_alg_coord_mul(&mut p, &a, coord, alg);
        prod[0][i] = p[0].clone();
        prod[1][i] = p[1].clone();
        prod[2][i] = p[2].clone();
        prod[3][i] = p[3].clone();
    }
}

/// `quat_lattice_alg_elem_mul(prod, lat, elem, alg)`.
pub fn quat_lattice_alg_elem_mul(
    prod: &mut QuatLattice,
    lat: &QuatLattice,
    elem: &QuatAlgElem,
    alg: &QuatAlg,
) {
    quat_lattice_mat_alg_coord_mul_without_hnf(&mut prod.basis, &lat.basis, &elem.coord, alg);
    ibz_mul(&mut prod.denom, &lat.denom, &elem.denom);
    quat_lattice_hnf(prod);
}

/// `quat_lattice_mul(res, lat1, lat2, alg)`.
pub fn quat_lattice_mul(
    res: &mut QuatLattice,
    lat1: &QuatLattice,
    lat2: &QuatLattice,
    alg: &QuatAlg,
) {
    let mut elem1 = ibz_vec_4_new();
    let mut elem2 = ibz_vec_4_new();
    let mut elem_res = ibz_vec_4_new();
    let mut generators: Vec<IbzVec4> = (0..16).map(|_| ibz_vec_4_new()).collect();
    let mut detmat = ibz_mat_4x4_new();
    let mut det = Ibz::zero();

    for k in 0..4 {
        ibz_vec_4_copy_ibz(
            &mut elem1,
            &lat1.basis[0][k],
            &lat1.basis[1][k],
            &lat1.basis[2][k],
            &lat1.basis[3][k],
        );
        for i in 0..4 {
            ibz_vec_4_copy_ibz(
                &mut elem2,
                &lat2.basis[0][i],
                &lat2.basis[1][i],
                &lat2.basis[2][i],
                &lat2.basis[3][i],
            );
            quat_alg_coord_mul(&mut elem_res, &elem1, &elem2, alg);
            for j in 0..4 {
                if k == 0 {
                    detmat[i][j] = elem_res[j].clone();
                }
                generators[4 * k + i][j] = elem_res[j].clone();
            }
        }
    }
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det), &detmat);
    let mut adet = Ibz::zero();
    ibz_abs(&mut adet, &det);
    det = adet;
    ibz_mat_4xn_hnf_mod_core(&mut res.basis, 16, &generators, &det);
    ibz_mul(&mut res.denom, &lat1.denom, &lat2.denom);
    let cloned = clone_lat(res);
    quat_lattice_reduce_denom(res, &cloned);
}

/// `quat_lattice_contains(coord, lat, x)`: returns 1 if `x` is in the
/// lattice and (if not NULL) writes its basis-coordinate vector.
pub fn quat_lattice_contains(
    coord: Option<&mut IbzVec4>,
    lat: &QuatLattice,
    x: &QuatAlgElem,
) -> i32 {
    let mut inv = ibz_mat_4x4_new();
    let mut det = Ibz::zero();
    let _ = ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), Some(&mut det), &lat.basis);
    assert!(ibz_is_zero(&det) == 0);
    let mut work_coord = ibz_vec_4_new();
    ibz_mat_4x4_eval(&mut work_coord, &inv, &x.coord);
    let cloned = clone_vec(&work_coord);
    ibz_vec_4_scalar_mul(&mut work_coord, &lat.denom, &cloned);
    let mut prod = Ibz::zero();
    ibz_mul(&mut prod, &x.denom, &det);
    let cloned = clone_vec(&work_coord);
    let divisible = ibz_vec_4_scalar_div(&mut work_coord, &prod, &cloned);
    if divisible != 0 {
        if let Some(out) = coord {
            for i in 0..4 {
                out[i] = work_coord[i].clone();
            }
        }
    }
    divisible
}

/// `quat_lattice_index(index, sublat, overlat)`.
pub fn quat_lattice_index(index: &mut Ibz, sublat: &QuatLattice, overlat: &QuatLattice) {
    let mut det = Ibz::zero();
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det), &sublat.basis);
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &overlat.denom, &overlat.denom);
    let cloned = tmp.clone();
    ibz_mul(&mut tmp, &cloned, &cloned);
    ibz_mul(index, &det, &tmp);
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &sublat.denom, &sublat.denom);
    let cloned = tmp.clone();
    ibz_mul(&mut tmp, &cloned, &cloned);
    let mut det = Ibz::zero();
    let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det), &overlat.basis);
    let cloned_tmp = tmp.clone();
    ibz_mul(&mut tmp, &cloned_tmp, &det);
    let saved_index = index.clone();
    let mut r = Ibz::zero();
    ibz_div(index, &mut r, &saved_index, &tmp);
    assert!(ibz_is_zero(&r) != 0);
    let mut a = Ibz::zero();
    ibz_abs(&mut a, index);
    *index = a;
}

/// `quat_lattice_gram(G, lattice, alg)`.
pub fn quat_lattice_gram(g: &mut IbzMat4x4, lattice: &QuatLattice, alg: &QuatAlg) {
    let mut tmp = Ibz::zero();
    for i in 0..4 {
        for j in 0..=i {
            ibz_set(&mut g[i][j], 0);
            for k in 0..4 {
                ibz_mul(&mut tmp, &lattice.basis[k][i], &lattice.basis[k][j]);
                if k >= 2 {
                    let cloned = tmp.clone();
                    ibz_mul(&mut tmp, &cloned, &alg.p);
                }
                let mut s = Ibz::zero();
                let g_ij = g[i][j].clone();
                use crate::ibz::ibz_add;
                ibz_add(&mut s, &g_ij, &tmp);
                g[i][j] = s;
            }
            let g_ij = g[i][j].clone();
            ibz_mul(&mut g[i][j], &g_ij, &crate::ibz::ibz_const_two());
        }
    }
    for i in 0..4 {
        for j in (i + 1)..4 {
            g[i][j] = g[j][i].clone();
        }
    }
}

fn clone_mat(m: &IbzMat4x4) -> IbzMat4x4 {
    let mut out = ibz_mat_4x4_new();
    ibz_mat_4x4_copy(&mut out, m);
    out
}
fn clone_vec(v: &IbzVec4) -> IbzVec4 {
    [v[0].clone(), v[1].clone(), v[2].clone(), v[3].clone()]
}
fn clone_lat(l: &QuatLattice) -> QuatLattice {
    QuatLattice {
        denom: l.denom.clone(),
        basis: clone_mat(&l.basis),
    }
}
