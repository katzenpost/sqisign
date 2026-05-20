//! LLL verification (oracle for non-unique L2 output).
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/lll/lll_verification.c`.
//! Given a candidate reduced basis, asserts that it satisfies the
//! size-reduction and Lovász conditions for the L2 parameters `(delta, eta)`
//! defined in `lll_internals.h`.
//!
//! This is the primary differential boundary for [`crate::lll::quat_lll_core`]:
//! the LLL output is not unique, so we cross-check both the Rust and the C
//! outputs against this exact-rational oracle. Bit-exact agreement is a
//! bonus; structural validity is the binding contract.

use crate::algebra::QuatAlg;
use crate::dim4::IbzMat4x4;
use crate::ibz::{ibz_const_one, ibz_const_two, ibz_set, Ibz};
use crate::lll::{DELTA_DENOM, DELTA_NUM, EPSILON_DENOM, EPSILON_NUM};
use crate::rationals::{
    ibq_abs, ibq_add, ibq_cmp, ibq_copy, ibq_inv, ibq_mul, ibq_set, ibq_sub, ibq_vec_4_new, Ibq,
    IbqMat4x4, IbqVec4,
};

/// `quat_lll_set_ibq_parameters(delta, eta)`. Returns L2 parameters as
/// rationals: `delta = 99/100` (the rational lower bound for the L2
/// constant `delta-bar`), and `eta = 1/2 + 1/100 = 51/100`.
pub fn quat_lll_set_ibq_parameters(delta: &mut Ibq, eta: &mut Ibq) {
    // delta initially set to 1/2 (the C reference uses it as a temporary).
    ibq_set(delta, &ibz_const_one(), &ibz_const_two());
    let mut num = Ibz::zero();
    let mut denom = Ibz::zero();
    ibz_set(&mut num, EPSILON_NUM);
    ibz_set(&mut denom, EPSILON_DENOM);
    ibq_set(eta, &num, &denom);
    // eta = eta + delta = (1/100) + (1/2). The C does ibq_add(eta, eta, delta).
    let eta_copy = eta.clone();
    ibq_add(eta, &eta_copy, delta);
    // Now set delta to the final value DELTA_NUM / DELTA_DENOM.
    ibz_set(&mut num, DELTA_NUM);
    ibz_set(&mut denom, DELTA_DENOM);
    ibq_set(delta, &num, &denom);
}

/// `ibq_vec_4_copy_ibz(vec, c0, c1, c2, c3)`: build an Ibq 4-vector from
/// four `Ibz` coefficients, all with denominator 1.
pub fn ibq_vec_4_copy_ibz(vec: &mut IbqVec4, c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) {
    let one = ibz_const_one();
    ibq_set(&mut vec[0], c0, &one);
    ibq_set(&mut vec[1], c1, &one);
    ibq_set(&mut vec[2], c2, &one);
    ibq_set(&mut vec[3], c3, &one);
}

/// `quat_lll_bilinear(b, v0, v1, q)`:
/// `b = v00*v10 + v01*v11 + q*(v02*v12 + v03*v13)`.
pub fn quat_lll_bilinear(b: &mut Ibq, vec0: &IbqVec4, vec1: &IbqVec4, q: &Ibz) {
    let one = ibz_const_one();
    let mut norm_q = Ibq::new();
    ibq_set(&mut norm_q, q, &one);

    let mut sum = Ibq::new();
    let mut prod = Ibq::new();
    ibq_mul(&mut sum, &vec0[0], &vec1[0]);
    ibq_mul(&mut prod, &vec0[1], &vec1[1]);
    let sum_copy = sum.clone();
    ibq_add(&mut sum, &sum_copy, &prod);
    ibq_mul(&mut prod, &vec0[2], &vec1[2]);
    let prod_copy = prod.clone();
    ibq_mul(&mut prod, &prod_copy, &norm_q);
    let sum_copy = sum.clone();
    ibq_add(&mut sum, &sum_copy, &prod);
    ibq_mul(&mut prod, &vec0[3], &vec1[3]);
    let prod_copy = prod.clone();
    ibq_mul(&mut prod, &prod_copy, &norm_q);
    ibq_add(b, &sum, &prod);
}

/// `quat_lll_gram_schmidt_transposed_with_ibq`: Gram-Schmidt orthogonalisation
/// over the rationals, written as four row vectors. The C reference uses an
/// auxiliary `ibq_mat_4x4_t work` indexed by row.
pub fn quat_lll_gram_schmidt_transposed_with_ibq(
    orthogonalised_transposed: &mut IbqMat4x4,
    mat: &IbzMat4x4,
    q: &Ibz,
) {
    let mut work: IbqMat4x4 = [
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
    ];

    // Initialise work[i] := column i of mat (i.e. transpose).
    for i in 0..4 {
        ibq_vec_4_copy_ibz(&mut work[i], &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
    }

    let mut norm = Ibq::new();
    let mut b = Ibq::new();
    let mut coeff = Ibq::new();
    let mut prod = Ibq::new();

    for i in 0..4 {
        let work_i = work[i].clone();
        quat_lll_bilinear(&mut norm, &work_i, &work_i, q);
        let norm_copy = norm.clone();
        ibq_inv(&mut norm, &norm_copy);
        for j in (i + 1)..4 {
            let mut vec = ibq_vec_4_new();
            ibq_vec_4_copy_ibz(&mut vec, &mat[0][j], &mat[1][j], &mat[2][j], &mat[3][j]);
            let work_i = work[i].clone();
            quat_lll_bilinear(&mut b, &work_i, &vec, q);
            ibq_mul(&mut coeff, &norm, &b);
            for k in 0..4 {
                ibq_mul(&mut prod, &coeff, &work[i][k]);
                let old = work[j][k].clone();
                ibq_sub(&mut work[j][k], &old, &prod);
            }
        }
    }

    for i in 0..4 {
        for j in 0..4 {
            ibq_copy(&mut orthogonalised_transposed[i][j], &work[i][j]);
        }
    }
}

/// `quat_lll_verify(mat, delta, eta, alg)`: returns `1` iff `mat` is a
/// `(delta, eta)`-LLL-reduced basis for the bilinear form determined by
/// `alg.p`.
///
/// The two conditions checked, exactly as in the C reference:
///  * size-reduction: `|<b_i*, b_j>| / <b_j*, b_j*> <= eta` for all `j < i`;
///  * Lovász: `<b_i*, b_i*> >= (delta - mu^2) * <b_{i-1}*, b_{i-1}*>`
///    where `mu = <b_{i-1}*, b_i>/<b_{i-1}*, b_{i-1}*>`.
pub fn quat_lll_verify(mat: &IbzMat4x4, delta: &Ibq, eta: &Ibq, alg: &QuatAlg) -> i32 {
    let mut res: i32 = 1;
    let mut orth = [
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
    ];
    quat_lll_gram_schmidt_transposed_with_ibq(&mut orth, mat, &alg.p);

    let mut tmp_vec = ibq_vec_4_new();
    let mut div = Ibq::new();
    let mut tmp = Ibq::new();
    let mut mu = Ibq::new();
    let mut norm = Ibq::new();
    let mut b = Ibq::new();

    // Size-reduction check
    for i in 0..4 {
        for j in 0..i {
            ibq_vec_4_copy_ibz(&mut tmp_vec, &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
            quat_lll_bilinear(&mut b, &orth[j], &tmp_vec, &alg.p);
            quat_lll_bilinear(&mut norm, &orth[j], &orth[j], &alg.p);
            ibq_inv(&mut tmp, &norm);
            ibq_mul(&mut mu, &b, &tmp);
            let mu_copy = mu.clone();
            ibq_abs(&mut mu, &mu_copy);
            res &= (ibq_cmp(&mu, eta) <= 0) as i32;
        }
    }
    // Lovász check
    for i in 1..4 {
        ibq_vec_4_copy_ibz(&mut tmp_vec, &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
        quat_lll_bilinear(&mut b, &orth[i - 1], &tmp_vec, &alg.p);
        quat_lll_bilinear(&mut norm, &orth[i - 1], &orth[i - 1], &alg.p);
        ibq_inv(&mut tmp, &norm);
        ibq_mul(&mut mu, &b, &tmp);
        // tmp = mu^2
        let mu_copy = mu.clone();
        ibq_mul(&mut tmp, &mu_copy, &mu_copy);
        // mu = delta - mu^2
        ibq_sub(&mut mu, delta, &tmp);
        quat_lll_bilinear(&mut tmp, &orth[i], &orth[i], &alg.p);
        // div = norm * (delta - mu^2)
        ibq_mul(&mut div, &norm, &mu);
        res &= (ibq_cmp(&tmp, &div) >= 0) as i32;
    }
    res
}
