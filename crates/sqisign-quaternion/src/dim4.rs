//! Dimension-4 ibz vectors and 4x4 ibz matrices.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/dim4.c`. All
//! operations are deterministic, pure-`Ibz`, and free of randomness.
//!
//! The C reference uses fixed-size C arrays `ibz_vec_4_t` and
//! `ibz_mat_4x4_t`. We mirror them as `[Ibz; 4]` and `[[Ibz; 4]; 4]`,
//! preserving row-major indexing `mat[row][col]`. Function names mirror
//! the C symbols (snake_case) so the differential boundary is easy to
//! audit.

use crate::ibz::{
    ibz_abs, ibz_add, ibz_cmp, ibz_const_zero, ibz_div, ibz_gcd, ibz_is_one, ibz_is_zero, ibz_mul,
    ibz_neg, ibz_set, ibz_sub, Ibz,
};

/// `ibz_vec_4_t`.
pub type IbzVec4 = [Ibz; 4];
/// `ibz_mat_4x4_t`.
pub type IbzMat4x4 = [[Ibz; 4]; 4];

/// Allocate a zero `IbzVec4`.
pub fn ibz_vec_4_new() -> IbzVec4 {
    [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()]
}

/// Allocate a zero `IbzMat4x4`.
pub fn ibz_mat_4x4_new() -> IbzMat4x4 {
    [
        ibz_vec_4_new(),
        ibz_vec_4_new(),
        ibz_vec_4_new(),
        ibz_vec_4_new(),
    ]
}

/// `ibz_vec_4_set(vec, c0, c1, c2, c3)`.
pub fn ibz_vec_4_set(vec: &mut IbzVec4, c0: i32, c1: i32, c2: i32, c3: i32) {
    ibz_set(&mut vec[0], c0);
    ibz_set(&mut vec[1], c1);
    ibz_set(&mut vec[2], c2);
    ibz_set(&mut vec[3], c3);
}

/// `ibz_vec_4_copy(new, vec)`.
pub fn ibz_vec_4_copy(new_: &mut IbzVec4, vec: &IbzVec4) {
    for i in 0..4 {
        new_[i] = vec[i].clone();
    }
}

/// `ibz_vec_4_copy_ibz(res, c0, c1, c2, c3)`.
pub fn ibz_vec_4_copy_ibz(res: &mut IbzVec4, c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) {
    res[0] = c0.clone();
    res[1] = c1.clone();
    res[2] = c2.clone();
    res[3] = c3.clone();
}

/// `ibz_vec_4_content(content, v)`: gcd of all four entries.
pub fn ibz_vec_4_content(content: &mut Ibz, v: &IbzVec4) {
    let mut t = Ibz::zero();
    ibz_gcd(&mut t, &v[0], &v[1]);
    let mut t2 = Ibz::zero();
    ibz_gcd(&mut t2, &v[2], &t);
    let mut t3 = Ibz::zero();
    ibz_gcd(&mut t3, &v[3], &t2);
    *content = t3;
}

/// `ibz_vec_4_negate(neg, vec)`.
pub fn ibz_vec_4_negate(neg: &mut IbzVec4, vec: &IbzVec4) {
    for i in 0..4 {
        let mut t = Ibz::zero();
        ibz_neg(&mut t, &vec[i]);
        neg[i] = t;
    }
}

/// `ibz_vec_4_add(res, a, b)`.
pub fn ibz_vec_4_add(res: &mut IbzVec4, a: &IbzVec4, b: &IbzVec4) {
    for i in 0..4 {
        let mut t = Ibz::zero();
        ibz_add(&mut t, &a[i], &b[i]);
        res[i] = t;
    }
}

/// `ibz_vec_4_sub(res, a, b)`.
pub fn ibz_vec_4_sub(res: &mut IbzVec4, a: &IbzVec4, b: &IbzVec4) {
    for i in 0..4 {
        let mut t = Ibz::zero();
        ibz_sub(&mut t, &a[i], &b[i]);
        res[i] = t;
    }
}

/// `ibz_vec_4_is_zero(x)`: 1 iff all four entries are zero.
pub fn ibz_vec_4_is_zero(x: &IbzVec4) -> i32 {
    let mut res = 1;
    for i in 0..4 {
        res &= ibz_is_zero(&x[i]);
    }
    res
}

/// `ibz_vec_4_linear_combination(lc, ca, va, cb, vb)`: `lc = ca*va + cb*vb`.
pub fn ibz_vec_4_linear_combination(
    lc: &mut IbzVec4,
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
) {
    let mut sums = ibz_vec_4_new();
    let mut prod = Ibz::zero();
    for i in 0..4 {
        ibz_mul(&mut sums[i], coeff_a, &vec_a[i]);
        ibz_mul(&mut prod, coeff_b, &vec_b[i]);
        let mut t = Ibz::zero();
        ibz_add(&mut t, &sums[i], &prod);
        sums[i] = t;
    }
    for i in 0..4 {
        lc[i] = sums[i].clone();
    }
}

/// `ibz_vec_4_scalar_mul(prod, scalar, vec)`.
pub fn ibz_vec_4_scalar_mul(prod: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4) {
    for i in 0..4 {
        let mut t = Ibz::zero();
        ibz_mul(&mut t, &vec[i], scalar);
        prod[i] = t;
    }
}

/// `ibz_vec_4_scalar_div(quot, scalar, vec)`: returns 1 if exact, 0 if any
/// remainder is non-zero. Division is always performed.
pub fn ibz_vec_4_scalar_div(quot: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4) -> i32 {
    let mut ok = 1;
    let mut r = Ibz::zero();
    for i in 0..4 {
        let mut q = Ibz::zero();
        ibz_div(&mut q, &mut r, &vec[i], scalar);
        quot[i] = q;
        if ibz_is_zero(&r) == 0 {
            ok = 0;
        }
    }
    ok
}

/// `ibz_mat_4x4_copy(new, mat)`.
pub fn ibz_mat_4x4_copy(new_: &mut IbzMat4x4, mat: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            new_[i][j] = mat[i][j].clone();
        }
    }
}

/// `ibz_mat_4x4_negate(neg, mat)`.
pub fn ibz_mat_4x4_negate(neg: &mut IbzMat4x4, mat: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            let mut t = Ibz::zero();
            ibz_neg(&mut t, &mat[i][j]);
            neg[i][j] = t;
        }
    }
}

/// `ibz_mat_4x4_transpose(transposed, mat)`.
pub fn ibz_mat_4x4_transpose(transposed: &mut IbzMat4x4, mat: &IbzMat4x4) {
    let mut work = ibz_mat_4x4_new();
    for i in 0..4 {
        for j in 0..4 {
            work[i][j] = mat[j][i].clone();
        }
    }
    ibz_mat_4x4_copy(transposed, &work);
}

/// `ibz_mat_4x4_zero(zero)`.
pub fn ibz_mat_4x4_zero(zero_: &mut IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut zero_[i][j], 0);
        }
    }
}

/// `ibz_mat_4x4_identity(id)`.
pub fn ibz_mat_4x4_identity(id: &mut IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut id[i][j], 0);
        }
        ibz_set(&mut id[i][i], 1);
    }
}

/// `ibz_mat_4x4_is_identity(mat)`.
///
/// Note: this mirrors the C reference's quirky implementation exactly,
/// which uses `ibz_is_one(&mat[i][j]) == (i == j)`. That comparison is
/// 1 iff `mat[i][j] == 1` and we're on the diagonal, OR `mat[i][j] != 1`
/// and we're off-diagonal. So strictly speaking this returns 1 iff every
/// diagonal entry is exactly 1 and every off-diagonal entry is anything
/// other than 1. It is **not** a strict identity check (an off-diagonal
/// of 5 still satisfies the predicate). We preserve the C semantics
/// verbatim for the differential boundary.
pub fn ibz_mat_4x4_is_identity(mat: &IbzMat4x4) -> i32 {
    let mut res = 1;
    for i in 0..4 {
        for j in 0..4 {
            let is_one = ibz_is_one(&mat[i][j]);
            let expected = if i == j { 1 } else { 0 };
            res = (res != 0 && is_one == expected) as i32;
        }
    }
    res
}

/// `ibz_mat_4x4_equal(mat1, mat2)`.
pub fn ibz_mat_4x4_equal(mat1: &IbzMat4x4, mat2: &IbzMat4x4) -> i32 {
    let mut diff = 0;
    for i in 0..4 {
        for j in 0..4 {
            diff |= ibz_cmp(&mat1[i][j], &mat2[i][j]);
        }
    }
    if diff == 0 {
        1
    } else {
        0
    }
}

/// `ibz_mat_4x4_scalar_mul(prod, scalar, mat)`.
pub fn ibz_mat_4x4_scalar_mul(prod: &mut IbzMat4x4, scalar: &Ibz, mat: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            let mut t = Ibz::zero();
            ibz_mul(&mut t, &mat[i][j], scalar);
            prod[i][j] = t;
        }
    }
}

/// `ibz_mat_4x4_gcd(gcd, mat)`: gcd of all sixteen entries.
///
/// Matches the C reference: it starts with `d = mat[0][0]` and folds
/// `gcd(d, mat[i][j])` over the whole matrix. The initial `mat[0][0]`
/// fold-in is redundant on iteration `(0,0)` but harmless.
pub fn ibz_mat_4x4_gcd(gcd: &mut Ibz, mat: &IbzMat4x4) {
    let mut d = mat[0][0].clone();
    for i in 0..4 {
        for j in 0..4 {
            let mut t = Ibz::zero();
            ibz_gcd(&mut t, &d, &mat[i][j]);
            d = t;
        }
    }
    *gcd = d;
}

/// `ibz_mat_4x4_scalar_div(quot, scalar, mat)`: returns 1 if exact.
pub fn ibz_mat_4x4_scalar_div(quot: &mut IbzMat4x4, scalar: &Ibz, mat: &IbzMat4x4) -> i32 {
    let mut ok = 1;
    let mut r = Ibz::zero();
    for i in 0..4 {
        for j in 0..4 {
            let mut q = Ibz::zero();
            ibz_div(&mut q, &mut r, &mat[i][j], scalar);
            quot[i][j] = q;
            if ibz_is_zero(&r) == 0 {
                ok = 0;
            }
        }
    }
    ok
}

/// `ibz_inv_dim4_make_coeff_pmp(coeff, a1, a2, b1, b2, c1, c2)`:
/// `coeff = a1*a2 - b1*b2 + c1*c2`.
pub fn ibz_inv_dim4_make_coeff_pmp(
    coeff: &mut Ibz,
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) {
    let mut sum = Ibz::zero();
    let mut prod = Ibz::zero();
    ibz_mul(&mut sum, a1, a2);
    ibz_mul(&mut prod, b1, b2);
    let mut t = Ibz::zero();
    ibz_sub(&mut t, &sum, &prod);
    ibz_mul(&mut prod, c1, c2);
    ibz_add(coeff, &t, &prod);
}

/// `ibz_inv_dim4_make_coeff_mpm(coeff, a1, a2, b1, b2, c1, c2)`:
/// `coeff = b1*b2 - a1*a2 - c1*c2`.
pub fn ibz_inv_dim4_make_coeff_mpm(
    coeff: &mut Ibz,
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) {
    let mut sum = Ibz::zero();
    let mut prod = Ibz::zero();
    ibz_mul(&mut sum, b1, b2);
    ibz_mul(&mut prod, a1, a2);
    let mut t = Ibz::zero();
    ibz_sub(&mut t, &sum, &prod);
    ibz_mul(&mut prod, c1, c2);
    ibz_sub(coeff, &t, &prod);
}

/// `ibz_mat_2x2_det_from_ibz(det, a11, a12, a21, a22)`:
/// `det = a11*a22 - a12*a21`. (Used as a helper of the 4x4 inverse.)
pub fn ibz_mat_2x2_det_from_ibz(det: &mut Ibz, a11: &Ibz, a12: &Ibz, a21: &Ibz, a22: &Ibz) {
    let mut prod = Ibz::zero();
    ibz_mul(&mut prod, a12, a21);
    let mut t = Ibz::zero();
    ibz_mul(&mut t, a11, a22);
    let mut d = Ibz::zero();
    ibz_sub(&mut d, &t, &prod);
    *det = d;
}

/// `ibz_mat_4x4_inv_with_det_as_denom(inv, det, mat)`: 4x4 determinant
/// and adjugate, structured as `mat * inv = det * I`. Returns 1 if
/// `det != 0` (the inverse exists), 0 otherwise.
///
/// Strictly transcribes the C reference's 2x2-minor Laplace expansion.
pub fn ibz_mat_4x4_inv_with_det_as_denom(
    inv: Option<&mut IbzMat4x4>,
    det: Option<&mut Ibz>,
    mat: &IbzMat4x4,
) -> i32 {
    let mut work = ibz_mat_4x4_new();
    let mut work_det = Ibz::zero();
    let mut s: [Ibz; 6] = [
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
    ];
    let mut c: [Ibz; 6] = [
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
    ];

    for i in 0..3 {
        ibz_mat_2x2_det_from_ibz(
            &mut s[i],
            &mat[0][0],
            &mat[0][i + 1],
            &mat[1][0],
            &mat[1][i + 1],
        );
        ibz_mat_2x2_det_from_ibz(
            &mut c[i],
            &mat[2][0],
            &mat[2][i + 1],
            &mat[3][0],
            &mat[3][i + 1],
        );
    }
    for i in 0..2 {
        ibz_mat_2x2_det_from_ibz(
            &mut s[3 + i],
            &mat[0][1],
            &mat[0][2 + i],
            &mat[1][1],
            &mat[1][2 + i],
        );
        ibz_mat_2x2_det_from_ibz(
            &mut c[3 + i],
            &mat[2][1],
            &mat[2][2 + i],
            &mat[3][1],
            &mat[3][2 + i],
        );
    }
    ibz_mat_2x2_det_from_ibz(&mut s[5], &mat[0][2], &mat[0][3], &mat[1][2], &mat[1][3]);
    ibz_mat_2x2_det_from_ibz(&mut c[5], &mat[2][2], &mat[2][3], &mat[3][2], &mat[3][3]);

    let mut prod = Ibz::zero();
    ibz_set(&mut work_det, 0);
    for i in 0..6 {
        ibz_mul(&mut prod, &s[i], &c[5 - i]);
        if i != 1 && i != 4 {
            let mut t = Ibz::zero();
            ibz_add(&mut t, &work_det, &prod);
            work_det = t;
        } else {
            let mut t = Ibz::zero();
            ibz_sub(&mut t, &work_det, &prod);
            work_det = t;
        }
    }

    // Transposed adjugate. Mirror the C indexing exactly, including the
    // `(j == 0)`, `(j > 1)`, etc. ternary patterns. The C uses bool->int
    // conversion (true => 1, false => 0) inside index expressions.
    for j in 0..4i32 {
        for k in 0..2i32 {
            let row = (1 - k) as usize;
            let i1 = (j == 0) as i32;
            let i2 = 2 - if j > 1 { 1 } else { 0 };
            let i3 = 3 - if j == 3 { 1 } else { 0 };
            let sj = 6 - j - (if j == 0 { 1 } else { 0 });
            let sj_mid = 4 - j - (if j == 1 { 1 } else { 0 });
            let sj_lo = 3 - j - (if j == 1 { 1 } else { 0 }) - (if j == 2 { 1 } else { 0 });
            let sign_pos = (k + j + 1) % 2 == 1;
            if sign_pos {
                ibz_inv_dim4_make_coeff_pmp(
                    &mut work[j as usize][k as usize],
                    &mat[row][i1 as usize],
                    &c[sj as usize],
                    &mat[row][i2 as usize],
                    &c[sj_mid as usize],
                    &mat[row][i3 as usize],
                    &c[sj_lo as usize],
                );
            } else {
                ibz_inv_dim4_make_coeff_mpm(
                    &mut work[j as usize][k as usize],
                    &mat[row][i1 as usize],
                    &c[sj as usize],
                    &mat[row][i2 as usize],
                    &c[sj_mid as usize],
                    &mat[row][i3 as usize],
                    &c[sj_lo as usize],
                );
            }
        }
        for k in 2..4i32 {
            let row = (3 - if k == 3 { 1 } else { 0 }) as usize;
            let i1 = (j == 0) as i32;
            let i2 = 2 - if j > 1 { 1 } else { 0 };
            let i3 = 3 - if j == 3 { 1 } else { 0 };
            let sj = 6 - j - (if j == 0 { 1 } else { 0 });
            let sj_mid = 4 - j - (if j == 1 { 1 } else { 0 });
            let sj_lo = 3 - j - (if j == 1 { 1 } else { 0 }) - (if j == 2 { 1 } else { 0 });
            let sign_pos = (k + j + 1) % 2 == 1;
            if sign_pos {
                ibz_inv_dim4_make_coeff_pmp(
                    &mut work[j as usize][k as usize],
                    &mat[row][i1 as usize],
                    &s[sj as usize],
                    &mat[row][i2 as usize],
                    &s[sj_mid as usize],
                    &mat[row][i3 as usize],
                    &s[sj_lo as usize],
                );
            } else {
                ibz_inv_dim4_make_coeff_mpm(
                    &mut work[j as usize][k as usize],
                    &mat[row][i1 as usize],
                    &s[sj as usize],
                    &mat[row][i2 as usize],
                    &s[sj_mid as usize],
                    &mat[row][i3 as usize],
                    &s[sj_lo as usize],
                );
            }
        }
    }

    let det_is_zero = ibz_is_zero(&work_det);
    if let Some(inv) = inv {
        let mut nz = Ibz::zero();
        ibz_set(&mut nz, if det_is_zero == 0 { 1 } else { 0 });
        ibz_mat_4x4_scalar_mul(inv, &nz, &work);
    }
    if let Some(det) = det {
        *det = work_det.clone();
    }
    // The C reference's return value is `!ibz_is_zero(det)` where `det`
    // is the *output pointer*; if the caller passed NULL, the reference
    // reads through the NULL pointer (undefined behaviour). Match the
    // observable case where det != NULL by returning the work_det test.
    if det_is_zero == 0 {
        1
    } else {
        0
    }
}

/// `ibz_mat_4x4_mul(res, a, b)`.
pub fn ibz_mat_4x4_mul(res: &mut IbzMat4x4, a: &IbzMat4x4, b: &IbzMat4x4) {
    let mut work = ibz_mat_4x4_new();
    let mut prod = Ibz::zero();
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut work[i][j], 0);
            for k in 0..4 {
                ibz_mul(&mut prod, &a[i][k], &b[k][j]);
                let mut t = Ibz::zero();
                ibz_add(&mut t, &work[i][j], &prod);
                work[i][j] = t;
            }
        }
    }
    ibz_mat_4x4_copy(res, &work);
}

/// `ibz_mat_4x4_eval(res, mat, vec)`: matrix * column-vector.
pub fn ibz_mat_4x4_eval(res: &mut IbzVec4, mat: &IbzMat4x4, vec: &IbzVec4) {
    let mut sum = ibz_vec_4_new();
    let mut prod = Ibz::zero();
    for i in 0..4 {
        for j in 0..4 {
            ibz_mul(&mut prod, &mat[i][j], &vec[j]);
            let mut t = Ibz::zero();
            ibz_add(&mut t, &sum[i], &prod);
            sum[i] = t;
        }
    }
    ibz_vec_4_copy(res, &sum);
}

/// `ibz_mat_4x4_eval_t(res, vec, mat)`: row-vector * matrix.
pub fn ibz_mat_4x4_eval_t(res: &mut IbzVec4, vec: &IbzVec4, mat: &IbzMat4x4) {
    let mut sum = ibz_vec_4_new();
    let mut prod = Ibz::zero();
    for i in 0..4 {
        for j in 0..4 {
            ibz_mul(&mut prod, &mat[j][i], &vec[j]);
            let mut t = Ibz::zero();
            ibz_add(&mut t, &sum[i], &prod);
            sum[i] = t;
        }
    }
    ibz_vec_4_copy(res, &sum);
}

/// `quat_qf_eval(res, qf, coord)`: scalar = coord^T * qf * coord.
///
/// The C reference reuses `sum[0]` as the accumulator on iteration `i==0`
/// (overwriting it instead of adding), so we mirror that.
pub fn quat_qf_eval(res: &mut Ibz, qf: &IbzMat4x4, coord: &IbzVec4) {
    let mut sum = ibz_vec_4_new();
    let mut prod = Ibz::zero();
    ibz_mat_4x4_eval(&mut sum, qf, coord);
    for i in 0..4 {
        ibz_mul(&mut prod, &sum[i], &coord[i]);
        if i > 0 {
            let mut t = Ibz::zero();
            ibz_add(&mut t, &sum[0], &prod);
            sum[0] = t;
        } else {
            sum[0] = prod.clone();
        }
    }
    *res = sum[0].clone();
    let _ = ibz_abs;
    let _ = ibz_const_zero;
}
