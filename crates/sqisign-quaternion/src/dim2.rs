//! Dimension-2 ibz vectors and 2x2 ibz matrices.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/dim2.c`.

use crate::ibz::{ibz_add, ibz_invmod, ibz_mod, ibz_mul, ibz_neg, ibz_set, ibz_sub, Ibz};

/// `ibz_vec_2_t`.
pub type IbzVec2 = [Ibz; 2];
/// `ibz_mat_2x2_t`.
pub type IbzMat2x2 = [[Ibz; 2]; 2];

pub fn ibz_vec_2_new() -> IbzVec2 {
    [Ibz::zero(), Ibz::zero()]
}

pub fn ibz_mat_2x2_new() -> IbzMat2x2 {
    [[Ibz::zero(), Ibz::zero()], [Ibz::zero(), Ibz::zero()]]
}

/// `ibz_vec_2_set(vec, a0, a1)`.
pub fn ibz_vec_2_set(vec: &mut IbzVec2, a0: i32, a1: i32) {
    ibz_set(&mut vec[0], a0);
    ibz_set(&mut vec[1], a1);
}

/// `ibz_mat_2x2_set(mat, a00, a01, a10, a11)`.
pub fn ibz_mat_2x2_set(mat: &mut IbzMat2x2, a00: i32, a01: i32, a10: i32, a11: i32) {
    ibz_set(&mut mat[0][0], a00);
    ibz_set(&mut mat[0][1], a01);
    ibz_set(&mut mat[1][0], a10);
    ibz_set(&mut mat[1][1], a11);
}

/// `ibz_mat_2x2_copy(copy, copied)`.
pub fn ibz_mat_2x2_copy(copy: &mut IbzMat2x2, copied: &IbzMat2x2) {
    for i in 0..2 {
        for j in 0..2 {
            copy[i][j] = copied[i][j].clone();
        }
    }
}

/// `ibz_mat_2x2_add(sum, a, b)`.
pub fn ibz_mat_2x2_add(sum: &mut IbzMat2x2, a: &IbzMat2x2, b: &IbzMat2x2) {
    for i in 0..2 {
        for j in 0..2 {
            let mut t = Ibz::zero();
            ibz_add(&mut t, &a[i][j], &b[i][j]);
            sum[i][j] = t;
        }
    }
}

/// `ibz_mat_2x2_det_from_ibz(det, a11, a12, a21, a22)`: `det = a11*a22 - a12*a21`.
pub fn ibz_mat_2x2_det_from_ibz(det: &mut Ibz, a11: &Ibz, a12: &Ibz, a21: &Ibz, a22: &Ibz) {
    let mut prod = Ibz::zero();
    ibz_mul(&mut prod, a12, a21);
    let mut t = Ibz::zero();
    ibz_mul(&mut t, a11, a22);
    let mut d = Ibz::zero();
    ibz_sub(&mut d, &t, &prod);
    *det = d;
}

/// `ibz_mat_2x2_eval(res, mat, vec)`.
pub fn ibz_mat_2x2_eval(res: &mut IbzVec2, mat: &IbzMat2x2, vec: &IbzVec2) {
    let mut prod = Ibz::zero();
    let mut matvec = ibz_vec_2_new();
    ibz_mul(&mut prod, &mat[0][0], &vec[0]);
    matvec[0] = prod.clone();
    ibz_mul(&mut prod, &mat[0][1], &vec[1]);
    let mut t = Ibz::zero();
    ibz_add(&mut t, &matvec[0], &prod);
    matvec[0] = t;

    ibz_mul(&mut prod, &mat[1][0], &vec[0]);
    matvec[1] = prod.clone();
    ibz_mul(&mut prod, &mat[1][1], &vec[1]);
    let mut t = Ibz::zero();
    ibz_add(&mut t, &matvec[1], &prod);
    matvec[1] = t;

    res[0] = matvec[0].clone();
    res[1] = matvec[1].clone();
}

/// `ibz_2x2_mul_mod(prod, a, b, m)`: 2x2 matmul reduced mod `m`.
pub fn ibz_2x2_mul_mod(prod: &mut IbzMat2x2, mat_a: &IbzMat2x2, mat_b: &IbzMat2x2, m: &Ibz) {
    let mut mul_ = Ibz::zero();
    let mut sums = ibz_mat_2x2_new();
    for i in 0..2 {
        for j in 0..2 {
            ibz_set(&mut sums[i][j], 0);
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                ibz_mul(&mut mul_, &mat_a[i][k], &mat_b[k][j]);
                let mut t = Ibz::zero();
                ibz_add(&mut t, &sums[i][j], &mul_);
                sums[i][j] = t;
                let mut r = Ibz::zero();
                ibz_mod(&mut r, &sums[i][j], m);
                sums[i][j] = r;
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            prod[i][j] = sums[i][j].clone();
        }
    }
}

/// `ibz_mat_2x2_inv_mod(inv, mat, m)`: 2x2 modular inverse via the
/// classical adjugate formula. Returns 1 on success, 0 on failure.
///
/// Mirrors the C reference exactly, including the quirk that on failure
/// the C code writes a 0 matrix into `inv` via the `ibz_set(&prod, res)`
/// dance.
pub fn ibz_mat_2x2_inv_mod(inv: &mut IbzMat2x2, mat: &IbzMat2x2, m: &Ibz) -> i32 {
    let mut det = Ibz::zero();
    let mut prod = Ibz::zero();
    ibz_mul(&mut det, &mat[0][0], &mat[1][1]);
    let mut tmp = Ibz::zero();
    ibz_mod(&mut tmp, &det, m);
    det = tmp;
    ibz_mul(&mut prod, &mat[0][1], &mat[1][0]);
    let mut tmp = Ibz::zero();
    ibz_sub(&mut tmp, &det, &prod);
    det = tmp;
    let mut tmp = Ibz::zero();
    ibz_mod(&mut tmp, &det, m);
    det = tmp;

    let mut det_inv = Ibz::zero();
    let res = ibz_invmod(&mut det_inv, &det, m);

    // The C reference does `ibz_set(&prod, res)` and multiplies det_inv
    // by that, zeroing it on failure.
    ibz_set(&mut prod, res);
    let mut t = Ibz::zero();
    ibz_mul(&mut t, &det_inv, &prod);
    det = t;

    // Compute adjugate then scale by det^{-1} mod m.
    let saved_00 = mat[0][0].clone();
    let new00 = mat[1][1].clone();
    inv[0][0] = new00;
    inv[1][1] = saved_00;
    let mut t = Ibz::zero();
    ibz_neg(&mut t, &mat[1][0]);
    inv[1][0] = t;
    let mut t = Ibz::zero();
    ibz_neg(&mut t, &mat[0][1]);
    inv[0][1] = t;

    for i in 0..2 {
        for j in 0..2 {
            let mut t = Ibz::zero();
            ibz_mul(&mut t, &inv[i][j], &det);
            inv[i][j] = t;
            let mut r = Ibz::zero();
            ibz_mod(&mut r, &inv[i][j], m);
            inv[i][j] = r;
        }
    }
    res
}
