//! `lat_ball.c`: bounding parallelogram and lattice ball sampling.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/lat_ball.c`. Only
//! the deterministic boundary [`quat_lattice_bound_parallelogram`] is in
//! scope this batch: the surrounding `quat_lattice_sample_from_ball`
//! depends on `ibz_rand_interval` (RNG-driven) and joins the deferred
//! RNG-driven boundary set alongside `ibz_rand_*` and
//! `ibz_generate_random_prime`.

use crate::dim4::{
    ibz_mat_4x4_identity, ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_scalar_mul, IbzMat4x4,
    IbzVec4,
};
use crate::ibz::{ibz_abs, ibz_div, ibz_is_one, ibz_is_zero, ibz_mul, ibz_sqrt_floor, Ibz};
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

// `quat_lattice_sample_from_ball` is deferred: it depends on
// `ibz_rand_interval`, which is part of the RNG-driven boundary set the
// quaternion port has not yet brought in.
