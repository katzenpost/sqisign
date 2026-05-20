//! `lat_ball.c`: bounding parallelogram and lattice ball sampling.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/lat_ball.c`.
//! Both the deterministic [`quat_lattice_bound_parallelogram`] and the
//! RNG-driven [`quat_lattice_sample_from_ball`] are in scope. The
//! sampler takes `&mut impl RngSource` as its first formal parameter,
//! replacing the reference's thread-local DRBG global; bytes drawn from
//! the trait match what the reference would have drawn under an
//! identically seeded DRBG.

use sqisign_common::RngSource;

use crate::algebra::{quat_alg_normalize, QuatAlg, QuatAlgElem};
use crate::dim4::{
    ibz_mat_4x4_eval, ibz_mat_4x4_eval_t, ibz_mat_4x4_identity, ibz_mat_4x4_inv_with_det_as_denom,
    ibz_mat_4x4_scalar_mul, ibz_vec_4_new, quat_qf_eval, IbzMat4x4, IbzVec4,
};
use crate::ibz::{
    ibz_abs, ibz_add, ibz_cmp, ibz_const_two, ibz_const_zero, ibz_div, ibz_is_one, ibz_is_zero,
    ibz_mul, ibz_sqrt_floor, ibz_sub, Ibz,
};
use crate::ibz_rand::ibz_rand_interval;
use crate::lattice::{quat_lattice_gram, QuatLattice};
use crate::lll::quat_lll_core;

/// `quat_lattice_bound_parallelogram(box, U, G, radius)`.
///
/// Given a Gram matrix `G` and a non-negative radius, computes a bounding
/// parallelogram for the ball of that radius in the lattice and the
/// integer change-of-basis matrix `U`. Returns `1` if the bounding box is
/// non-trivial (some `box[i]` non-zero), `0` if the box collapses to the
/// origin (in which case the only lattice point in the ball is `0`).
///
/// The C reference relies on three intermediate facts:
///  1. the dual Gram matrix has the same lattice geometry as the primal;
///  2. running [`quat_lll_core`] on the dual Gram matrix yields a
///     reduced basis with diagonal entries that are reciprocals of the
///     squared GS norms (up to the shared denominator);
///  3. `U` is unitary, so its inverse multiplied by the denominator
///     reconstructs the integer change-of-basis without losing precision.
pub fn quat_lattice_bound_parallelogram(
    box_: &mut IbzVec4,
    u_mat: &mut IbzMat4x4,
    g: &IbzMat4x4,
    radius: &Ibz,
) -> i32 {
    let mut denom = Ibz::zero();
    let mut rem = Ibz::zero();
    let mut dual_g = crate::dim4::ibz_mat_4x4_new();

    // Compute the Gram matrix of the dual lattice.
    let _inv_check = ibz_mat_4x4_inv_with_det_as_denom(Some(&mut dual_g), Some(&mut denom), g);
    debug_assert!(_inv_check != 0, "dual Gram inverse failed");

    // Initialize the dual lattice basis to the identity matrix.
    ibz_mat_4x4_identity(u_mat);

    // Reduce the dual lattice.
    quat_lll_core(&mut dual_g, u_mat);

    // Compute the parallelogram's bounds.
    let mut trivial: i32 = 1;
    for i in 0..4 {
        let mut prod = Ibz::zero();
        ibz_mul(&mut prod, &dual_g[i][i], radius);
        let mut q = Ibz::zero();
        ibz_div(&mut q, &mut rem, &prod, &denom);
        let q_clone = q.clone();
        ibz_sqrt_floor(&mut box_[i], &q_clone);
        trivial &= ibz_is_zero(&box_[i]);
    }

    // Compute the transpose transformation matrix.
    let u_clone = u_mat.clone();
    let _inv = ibz_mat_4x4_inv_with_det_as_denom(Some(u_mat), Some(&mut denom), &u_clone);
    // U is unitary, det(U) = +/-1.
    let u_clone = u_mat.clone();
    ibz_mat_4x4_scalar_mul(u_mat, &denom, &u_clone);
    debug_assert!(_inv != 0, "U inversion failed");
    let mut denom_abs = Ibz::zero();
    ibz_abs(&mut denom_abs, &denom);
    debug_assert!(
        ibz_is_one(&denom_abs) != 0,
        "U should be unitary (det = +/-1)"
    );

    if trivial != 0 {
        0
    } else {
        1
    }
}

/// `quat_lattice_sample_from_ball(rng, res, lattice, alg, radius)`.
///
/// Samples a non-zero lattice element of (twice) the squared norm at most
/// `radius * lattice.denom^2 * 2` (the C reference's "Gram matrix is
/// twice the norm" correction is preserved verbatim). Returns 1 on
/// success, 0 if the bounding box collapses to the origin (only the zero
/// element lies inside the ball).
///
/// The rejection loop mirrors the reference exactly: per coordinate, if
/// `box[i] == 0` then `x[i] = 0`, else draw `x[i]` from
/// `ibz_rand_interval(rng, 0, 2 * box[i])` and shift by `-box[i]`. The
/// `rng` argument replaces the reference's thread-local DRBG global; the
/// byte stream consumed is identical for an identically seeded DRBG.
pub fn quat_lattice_sample_from_ball<R: RngSource>(
    rng: &mut R,
    res: &mut QuatAlgElem,
    lattice: &QuatLattice,
    alg: &QuatAlg,
    radius: &Ibz,
) -> i32 {
    assert!(
        ibz_cmp(radius, &ibz_const_zero()) > 0,
        "quat_lattice_sample_from_ball: radius must be positive"
    );

    let mut box_ = ibz_vec_4_new();
    let mut u_mat = crate::dim4::ibz_mat_4x4_new();
    let mut g = crate::dim4::ibz_mat_4x4_new();
    let mut x = ibz_vec_4_new();
    let mut rad = Ibz::zero();
    let mut tmp = Ibz::zero();

    // Compute the Gram matrix of the lattice.
    quat_lattice_gram(&mut g, lattice, alg);

    // Correct ball radius by the denominator.
    ibz_mul(&mut rad, radius, &lattice.denom);
    let rad_clone = rad.clone();
    ibz_mul(&mut rad, &rad_clone, &lattice.denom);
    // Correct by 2 (Gram matrix corresponds to twice the norm).
    let rad_clone = rad.clone();
    ibz_mul(&mut rad, &rad_clone, &ibz_const_two());

    // Compute a bounding parallelogram for the ball, stop if it only
    // contains the origin.
    let mut ok = quat_lattice_bound_parallelogram(&mut box_, &mut u_mat, &g, &rad);
    if ok == 0 {
        return 0;
    }

    // Rejection sampling from the parallelogram.
    loop {
        // Sample vector.
        for i in 0..4 {
            if ibz_is_zero(&box_[i]) != 0 {
                x[i] = ibz_const_zero();
            } else {
                let box_i = box_[i].clone();
                ibz_add(&mut tmp, &box_i, &box_i);
                ok &= ibz_rand_interval(rng, &mut x[i], &ibz_const_zero(), &tmp);
                let cur = x[i].clone();
                ibz_sub(&mut x[i], &cur, &box_[i]);
                if ok == 0 {
                    return ok;
                }
            }
        }
        // Map to parallelogram.
        let x_clone = x.clone();
        ibz_mat_4x4_eval_t(&mut x, &x_clone, &u_mat);
        // Evaluate quadratic form.
        quat_qf_eval(&mut tmp, &g, &x);
        if ibz_is_zero(&tmp) == 0 && ibz_cmp(&tmp, &rad) <= 0 {
            break;
        }
    }

    // Evaluate linear combination.
    ibz_mat_4x4_eval(&mut res.coord, &lattice.basis, &x);
    res.denom = lattice.denom.clone();
    quat_alg_normalize(res);

    ok
}
