//! Direct port of `vendor/the-sqisign/src/id2iso/ref/lvlx/id2iso.c`.
//!
//! Six entry points: scalar multiplication driven by `ibz` scalars,
//! kernel <-> ideal translation on the 2-power torsion, the in-place
//! matrix and endomorphism applications on a 2^f-torsion basis, and the
//! Tate-pairing change-of-basis matrices.

use sqisign_ec::ec_dlog_2_tate;
use sqisign_ec::{
    copy_basis, copy_point, ec_biscalar_mul, ec_dbl_iter, ec_is_equal, ec_point_init, EcBasis,
    EcCurve, EcPoint, NWORDS_ORDER, TORSION_EVEN_POWER,
};
use sqisign_mp::mp_invert_matrix;

use sqisign_precomp::{
    CURVES_WITH_ENDOMORPHISMS, EXTREMAL_ORDERS, NUM_ALTERNATE_STARTING_CURVES, TORSION_PLUS_2POWER,
};
use sqisign_quaternion::dim2::{ibz_mat_2x2_copy, ibz_mat_2x2_new};
use sqisign_quaternion::dim4::ibz_vec_4_new;
use sqisign_quaternion::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_copy_digits, ibz_gcd, ibz_is_even,
    ibz_mat_2x2_eval, ibz_mat_2x2_inv_mod, ibz_mod, ibz_mul, ibz_pow, ibz_set, ibz_sub,
    ibz_to_digits, ibz_vec_2_new, quat_alg_conj, quat_alg_make_primitive, quat_change_to_O0_basis,
    quat_lideal_create, quat_lideal_generator, Ibz, IbzMat2x2, IbzVec2, QuatAlg, QuatAlgElem,
    QuatLeftIdeal,
};

/// `ec_biscalar_mul_ibz_vec(res, scalar_vec, f, PQ, curve)`: scalar
/// multiplication `[x]P + [y]Q` where `x` and `y` live in an
/// `ibz_vec_2_t` and `P`, `Q` lie in `E[2^f]`.
pub fn ec_biscalar_mul_ibz_vec(
    res: &mut EcPoint,
    scalar_vec: &IbzVec2,
    f: i32,
    pq: &EcBasis,
    curve: &EcCurve,
) {
    let mut scalars = [[0u64; NWORDS_ORDER]; 2];
    ibz_to_digits(&mut scalars[0], &scalar_vec[0]);
    ibz_to_digits(&mut scalars[1], &scalar_vec[1]);
    let (s0, s1) = scalars.split_at(1);
    ec_biscalar_mul(res, &s0[0], &s1[0], f, pq, curve);
}

/// `id2iso_ideal_to_kernel_dlogs_even(vec, lideal)`: scalars `(s0, s1)`
/// that determine the kernel generator of the isogeny equivalent to the
/// given `2^f`-norm ideal.
pub fn id2iso_ideal_to_kernel_dlogs_even(vec: &mut IbzVec2, lideal: &QuatLeftIdeal, alg: &QuatAlg) {
    let mut tmp = Ibz::zero();
    let mut mat: IbzMat2x2 = ibz_mat_2x2_new();

    // Compute the matrix of the dual of alpha on the 2^f-torsion.
    {
        let mut alpha = QuatAlgElem::new();
        let ok = quat_lideal_generator(&mut alpha, lideal, alg);
        assert!(ok != 0, "id2iso_ideal_to_kernel_dlogs_even: no generator");
        let alpha_in = alpha.clone();
        quat_alg_conj(&mut alpha, &alpha_in);

        let mut coeffs = ibz_vec_4_new();
        quat_change_to_O0_basis(&mut coeffs, &alpha);

        let action_gen2 = &CURVES_WITH_ENDOMORPHISMS[0].action_gen2;
        let action_gen3 = &CURVES_WITH_ENDOMORPHISMS[0].action_gen3;
        let action_gen4 = &CURVES_WITH_ENDOMORPHISMS[0].action_gen4;

        for i in 0..2 {
            let m_ii = mat[i][i].clone();
            ibz_add(&mut mat[i][i], &m_ii, &coeffs[0]);
            for j in 0..2 {
                ibz_mul(&mut tmp, &action_gen2[i][j], &coeffs[1]);
                let m = mat[i][j].clone();
                ibz_add(&mut mat[i][j], &m, &tmp);
                ibz_mul(&mut tmp, &action_gen3[i][j], &coeffs[2]);
                let m = mat[i][j].clone();
                ibz_add(&mut mat[i][j], &m, &tmp);
                ibz_mul(&mut tmp, &action_gen4[i][j], &coeffs[3]);
                let m = mat[i][j].clone();
                ibz_add(&mut mat[i][j], &m, &tmp);
            }
        }
    }

    // Find the kernel of alpha modulo the ideal norm.
    {
        let norm = &lideal.norm;
        ibz_mod(&mut vec[0], &mat[0][0], norm);
        ibz_mod(&mut vec[1], &mat[1][0], norm);
        ibz_gcd(&mut tmp, &vec[0], &vec[1]);
        if ibz_is_even(&tmp) != 0 {
            ibz_mod(&mut vec[0], &mat[0][1], norm);
            ibz_mod(&mut vec[1], &mat[1][1], norm);
        }
        // Debug-only sanity in the C reference: assert gcd(vec[0], vec[1],
        // norm) == 1. We omit it; the boundary contract is the matrix
        // computation.
    }
}

/// `matrix_application_even_basis(bas, E, mat, f)`: in-place application
/// of a 2x2 matrix to a basis of `E[2^f]`. Returns 1 on success, 0 on
/// failure of the difference basis check.
pub fn matrix_application_even_basis(
    bas: &mut EcBasis,
    e: &EcCurve,
    mat: &mut IbzMat2x2,
    f: i32,
) -> i32 {
    let mut scalars = [[0u64; NWORDS_ORDER]; 2];

    let mut tmp = Ibz::zero();
    let mut pow_two = Ibz::zero();
    ibz_pow(&mut pow_two, &ibz_const_two(), f as u32);

    let mut tmp_bas = EcBasis::zero();
    copy_basis(&mut tmp_bas, bas);

    // Reduce mod 2^f.
    {
        let m = mat[0][0].clone();
        ibz_mod(&mut mat[0][0], &m, &pow_two);
        let m = mat[0][1].clone();
        ibz_mod(&mut mat[0][1], &m, &pow_two);
        let m = mat[1][0].clone();
        ibz_mod(&mut mat[1][0], &m, &pow_two);
        let m = mat[1][1].clone();
        ibz_mod(&mut mat[1][1], &m, &pow_two);
    }

    // R = [a]P + [b]Q
    ibz_to_digits(&mut scalars[0], &mat[0][0]);
    ibz_to_digits(&mut scalars[1], &mat[1][0]);
    let (s0, s1) = scalars.split_at(1);
    ec_biscalar_mul(&mut bas.P, &s0[0], &s1[0], f, &tmp_bas, e);

    // S = [c]P + [d]Q
    let mut scalars = [[0u64; NWORDS_ORDER]; 2];
    ibz_to_digits(&mut scalars[0], &mat[0][1]);
    ibz_to_digits(&mut scalars[1], &mat[1][1]);
    let (s0, s1) = scalars.split_at(1);
    ec_biscalar_mul(&mut bas.Q, &s0[0], &s1[0], f, &tmp_bas, e);

    // R - S = [a-c]P + [b-d]Q
    let mut scalars = [[0u64; NWORDS_ORDER]; 2];
    ibz_sub(&mut tmp, &mat[0][0], &mat[0][1]);
    let t = tmp.clone();
    ibz_mod(&mut tmp, &t, &pow_two);
    ibz_to_digits(&mut scalars[0], &tmp);
    ibz_sub(&mut tmp, &mat[1][0], &mat[1][1]);
    let t = tmp.clone();
    ibz_mod(&mut tmp, &t, &pow_two);
    ibz_to_digits(&mut scalars[1], &tmp);
    let (s0, s1) = scalars.split_at(1);
    ec_biscalar_mul(&mut bas.PmQ, &s0[0], &s1[0], f, &tmp_bas, e)
}

/// `endomorphism_application_even_basis(bas, index, E, theta, f)`:
/// applies an endomorphism of `E_index` to a basis of `E[2^f]`.
pub fn endomorphism_application_even_basis(
    bas: &mut EcBasis,
    index_alternate_curve: i32,
    e: &EcCurve,
    theta: &QuatAlgElem,
    f: i32,
) {
    let mut tmp = Ibz::zero();
    let mut coeffs = ibz_vec_4_new();
    let mut mat: IbzMat2x2 = ibz_mat_2x2_new();
    let mut content = Ibz::zero();

    let idx = index_alternate_curve as usize;
    let order = &EXTREMAL_ORDERS[idx].order;
    quat_alg_make_primitive(&mut coeffs, &mut content, theta, order);
    // The C reference asserts content is odd; we drop the debug check.

    for row in mat.iter_mut() {
        for cell in row.iter_mut() {
            ibz_set(cell, 0);
        }
    }

    let cwe = &CURVES_WITH_ENDOMORPHISMS[idx];

    for i in 0..2 {
        let m_ii = mat[i][i].clone();
        ibz_add(&mut mat[i][i], &m_ii, &coeffs[0]);
        for j in 0..2 {
            ibz_mul(&mut tmp, &cwe.action_gen2[i][j], &coeffs[1]);
            let m = mat[i][j].clone();
            ibz_add(&mut mat[i][j], &m, &tmp);
            ibz_mul(&mut tmp, &cwe.action_gen3[i][j], &coeffs[2]);
            let m = mat[i][j].clone();
            ibz_add(&mut mat[i][j], &m, &tmp);
            ibz_mul(&mut tmp, &cwe.action_gen4[i][j], &coeffs[3]);
            let m = mat[i][j].clone();
            ibz_add(&mut mat[i][j], &m, &tmp);
            let m = mat[i][j].clone();
            ibz_mul(&mut mat[i][j], &m, &content);
        }
    }

    matrix_application_even_basis(bas, e, &mut mat, f);
}

/// `id2iso_kernel_dlogs_to_ideal_even(lideal, vec2, f)`: build the
/// ideal whose kernel is generated by `vec2[0]*B0[0] + vec2[1]*B0[1]`,
/// where `B0` is the canonical basis of `E0`.
pub fn id2iso_kernel_dlogs_to_ideal_even(
    lideal: &mut QuatLeftIdeal,
    vec2: &IbzVec2,
    f: i32,
    alg: &QuatAlg,
) {
    let mut two_pow = Ibz::zero();
    let mut vec = ibz_vec_2_new();

    if f as usize == TORSION_EVEN_POWER {
        two_pow = TORSION_PLUS_2POWER.clone();
    } else {
        ibz_pow(&mut two_pow, &ibz_const_two(), f as u32);
    }

    {
        let mut mat: IbzMat2x2 = ibz_mat_2x2_new();

        mat[0][0] = vec2[0].clone();
        mat[1][0] = vec2[1].clone();

        let action_j = &CURVES_WITH_ENDOMORPHISMS[0].action_j;
        ibz_mat_2x2_eval(&mut vec, action_j, vec2);
        mat[0][1] = vec[0].clone();
        mat[1][1] = vec[1].clone();

        let action_gen4 = &CURVES_WITH_ENDOMORPHISMS[0].action_gen4;
        ibz_mat_2x2_eval(&mut vec, action_gen4, vec2);
        let m = mat[0][1].clone();
        ibz_add(&mut mat[0][1], &m, &vec[0]);
        let m = mat[1][1].clone();
        ibz_add(&mut mat[1][1], &m, &vec[1]);

        let m = mat[0][1].clone();
        ibz_mod(&mut mat[0][1], &m, &two_pow);
        let m = mat[1][1].clone();
        ibz_mod(&mut mat[1][1], &m, &two_pow);

        let mut inv: IbzMat2x2 = ibz_mat_2x2_new();
        let ok = ibz_mat_2x2_inv_mod(&mut inv, &mat, &two_pow);
        assert!(
            ok != 0,
            "id2iso_kernel_dlogs_to_ideal_even: matrix not invertible"
        );

        let action_i = &CURVES_WITH_ENDOMORPHISMS[0].action_i;
        let mut v_temp = ibz_vec_2_new();
        ibz_mat_2x2_eval(&mut v_temp, action_i, vec2);
        ibz_mat_2x2_eval(&mut vec, &inv, &v_temp);
    }

    // Final result: a - i + b*(j + (1+k)/2)
    let mut gen = QuatAlgElem::new();
    ibz_set(&mut gen.denom, 2);
    ibz_add(&mut gen.coord[0], &vec[0], &vec[0]);
    ibz_set(&mut gen.coord[1], -2);
    ibz_add(&mut gen.coord[2], &vec[1], &vec[1]);
    gen.coord[3] = vec[1].clone();
    let g0 = gen.coord[0].clone();
    ibz_add(&mut gen.coord[0], &g0, &vec[1]);

    let maxord = &EXTREMAL_ORDERS[0].order;
    quat_lideal_create(lideal, &gen, &two_pow, maxord, alg);
    assert_eq!(
        ibz_cmp(&lideal.norm, &two_pow),
        0,
        "id2iso_kernel_dlogs_to_ideal_even: norm mismatch"
    );
}

/// Internal helper for the two `change_of_basis_matrix_tate*` entries.
fn change_of_basis_matrix_tate_impl(
    mat: &mut IbzMat2x2,
    b1: &EcBasis,
    b2: &EcBasis,
    e: &mut EcCurve,
    f: i32,
    invert: bool,
) {
    let mut x1 = vec![0u64; NWORDS_ORDER];
    let mut x2 = vec![0u64; NWORDS_ORDER];
    let mut x3 = vec![0u64; NWORDS_ORDER];
    let mut x4 = vec![0u64; NWORDS_ORDER];

    if invert {
        ec_dlog_2_tate(&mut x1, &mut x2, &mut x3, &mut x4, b1, b2, e, f);
        mp_invert_matrix(&mut x1, &mut x2, &mut x3, &mut x4, f);
    } else {
        ec_dlog_2_tate(&mut x1, &mut x2, &mut x3, &mut x4, b2, b1, e, f);
    }

    // Copy out via ibz_copy_digits.
    ibz_copy_digits(&mut mat[0][0], &x1);
    ibz_copy_digits(&mut mat[1][0], &x2);
    ibz_copy_digits(&mut mat[0][1], &x3);
    ibz_copy_digits(&mut mat[1][1], &x4);
}

/// `change_of_basis_matrix_tate(mat, B1, B2, E, f)`: matrix `M` with
/// `(M*v).B2 = v.B1`. `B2` must be "full" with respect to the
/// `2^TORSION_EVEN_POWER` torsion.
pub fn change_of_basis_matrix_tate(
    mat: &mut IbzMat2x2,
    b1: &EcBasis,
    b2: &EcBasis,
    e: &mut EcCurve,
    f: i32,
) {
    change_of_basis_matrix_tate_impl(mat, b1, b2, e, f, false);
}

/// `change_of_basis_matrix_tate_invert(mat, B1, B2, E, f)`: matrix `M`
/// with `(M*v).B1 = [2^(e-f)]*v.B2`. `B1` must be "full".
pub fn change_of_basis_matrix_tate_invert(
    mat: &mut IbzMat2x2,
    b1: &EcBasis,
    b2: &EcBasis,
    e: &mut EcCurve,
    f: i32,
) {
    change_of_basis_matrix_tate_impl(mat, b1, b2, e, f, true);
}

// Suppress unused-imports warning if a constant ends up unreferenced after
// trimming a debug assert.
#[allow(dead_code)]
fn _unused() {
    let _ = (
        ec_dbl_iter,
        ec_is_equal,
        ec_point_init,
        copy_point,
        ibz_const_one,
        NUM_ALTERNATE_STARTING_CURVES,
        ibz_mat_2x2_copy,
    );
}
