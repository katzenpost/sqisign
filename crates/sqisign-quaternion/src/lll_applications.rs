//! LLL-derived helpers on lattices and ideals.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/lll/lll_applications.c`
//! for the deterministic boundaries. The RNG-driven boundary
//! `quat_lideal_prime_norm_reduced_equivalent` is deferred to the broader
//! RNG-driven boundary handling (the rest of the quaternion module follows
//! the same convention: `ibz_rand_*` paths are not in scope for this batch).

use crate::algebra::QuatAlg;
use crate::dim4::{ibz_mat_4x4_copy, ibz_mat_4x4_scalar_mul, IbzMat4x4};
use crate::ibz::{ibz_div_2exp, ibz_mul, ibz_set, Ibz};
use crate::ideal::{quat_lideal_class_gram, quat_lideal_norm, QuatLeftIdeal};
use crate::lattice::quat_lattice_mul;
use crate::lll::quat_lll_core;

/// `quat_lideal_reduce_basis(reduced, gram, lideal, alg)`.
///
/// Runs the L2 core on the lideal's class Gram matrix and stores the
/// reduced basis. The Gram matrix returned is the **reduced** Gram matrix
/// scaled by `lideal.lattice.denom^2`, with its diagonal halved and the
/// strict upper triangle zeroed, matching the C reference exactly.
pub fn quat_lideal_reduce_basis(
    reduced: &mut IbzMat4x4,
    gram: &mut IbzMat4x4,
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    let mut gram_corrector = Ibz::zero();
    ibz_mul(
        &mut gram_corrector,
        &lideal.lattice.denom,
        &lideal.lattice.denom,
    );
    quat_lideal_class_gram(gram, lideal, alg);
    ibz_mat_4x4_copy(reduced, &lideal.lattice.basis);
    quat_lll_core(gram, reduced);
    let gram_clone = gram.clone();
    ibz_mat_4x4_scalar_mul(gram, &gram_corrector, &gram_clone);
    for i in 0..4 {
        let cell = gram[i][i].clone();
        ibz_div_2exp(&mut gram[i][i], &cell, 1);
        for j in (i + 1)..4 {
            ibz_set(&mut gram[i][j], 0);
        }
    }
}

/// `quat_lideal_lideal_mul_reduced(prod, gram, lideal1, lideal2, alg)`.
///
/// Forms the ideal product `lideal1 * lideal2` (lattice multiplication
/// followed by norm recomputation) and reduces its basis in place. The
/// reduced Gram matrix is returned in `gram` with the same convention as
/// [`quat_lideal_reduce_basis`] (scaled, halved diagonal, upper triangle
/// zeroed).
pub fn quat_lideal_lideal_mul_reduced(
    prod: &mut QuatLeftIdeal,
    gram: &mut IbzMat4x4,
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    parent_order: &crate::lattice::QuatLattice,
    alg: &QuatAlg,
) {
    let mut red = crate::dim4::ibz_mat_4x4_new();
    quat_lattice_mul(&mut prod.lattice, &lideal1.lattice, &lideal2.lattice, alg);
    quat_lideal_norm(prod, parent_order);
    quat_lideal_reduce_basis(&mut red, gram, prod, alg);
    ibz_mat_4x4_copy(&mut prod.lattice.basis, &red);
}

// The C reference also defines `quat_lideal_prime_norm_reduced_equivalent`,
// which depends on `ibz_rand_interval_minm_m` and the surrounding RNG
// machinery (not yet ported). It belongs to the deferred RNG-driven boundary
// batch alongside `ibz_generate_random_prime`.
