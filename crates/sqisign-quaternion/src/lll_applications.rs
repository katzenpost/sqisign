//! LLL-derived helpers on lattices and ideals.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/lll/lll_applications.c`.
//! The RNG-driven boundary [`quat_lideal_prime_norm_reduced_equivalent`]
//! takes `&mut impl RngSource` as its first formal parameter, replacing
//! the reference's thread-local DRBG global.

use sqisign_common::RngSource;

use crate::algebra::{quat_alg_conj, QuatAlg, QuatAlgElem};
use crate::dim4::{
    ibz_mat_4x4_copy, ibz_mat_4x4_eval, ibz_mat_4x4_new, ibz_mat_4x4_scalar_mul, quat_qf_eval,
    IbzMat4x4,
};
use crate::ibz::{ibz_div, ibz_div_2exp, ibz_is_zero, ibz_mul, ibz_probab_prime, ibz_set, Ibz};
use crate::ibz_rand::ibz_rand_interval_minm_m;
use crate::ideal::{quat_lideal_class_gram, quat_lideal_mul, quat_lideal_norm, QuatLeftIdeal};
use crate::lattice::{quat_lattice_contains, quat_lattice_mul, QuatLattice};
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

/// `quat_lideal_prime_norm_reduced_equivalent(rng, lideal, alg,
/// primality_num_iter, equiv_bound_coeff, parent_order)`.
///
/// Searches a small box of `[-equiv_bound_coeff, equiv_bound_coeff]^4`
/// integer combinations of the LLL-reduced basis of `lideal` for an
/// element whose (lideal-adjusted) norm is prime. On success, replaces
/// `lideal` with the equivalent left ideal of prime norm and returns 1.
/// On exhaustion of the bounded search returns 0; the reference instead
/// `assert(found)`s, which we surface as a return-value-only signal so
/// the caller may decide.
///
/// `parent_order` is the (maximal) order the ideal lives in; it is
/// needed by [`quat_lideal_norm`], which the C reference invokes via the
/// pointer stored in `lideal.parent_order`. The Rust port carries that
/// reference at the call site rather than on the struct, matching the
/// rest of the ideal/lattice ports.
///
/// The `rng` argument replaces the reference's thread-local DRBG global;
/// bytes drawn match what the reference would have drawn under an
/// identically seeded DRBG.
pub fn quat_lideal_prime_norm_reduced_equivalent<R: RngSource>(
    rng: &mut R,
    lideal: &mut QuatLeftIdeal,
    alg: &QuatAlg,
    primality_num_iter: i32,
    equiv_bound_coeff: i32,
    parent_order: &QuatLattice,
) -> i32 {
    let mut gram = ibz_mat_4x4_new();
    let mut red = ibz_mat_4x4_new();

    let mut found = 0i32;

    // Compute the reduced basis.
    quat_lideal_reduce_basis(&mut red, &mut gram, lideal, alg);

    let mut new_alpha = QuatAlgElem::new();
    let mut tmp = Ibz::zero();
    let mut remainder = Ibz::zero();
    let mut adjusted_norm = Ibz::zero();

    ibz_mul(&mut adjusted_norm, &lideal.lattice.denom, &lideal.lattice.denom);

    let mut ctr = 0i32;

    // equiv_num_iter = (2 * equiv_bound_coeff + 1)^4.
    assert!(
        equiv_bound_coeff < (1 << 20),
        "quat_lideal_prime_norm_reduced_equivalent: equiv_bound_coeff too large"
    );
    let mut equiv_num_iter = 2 * equiv_bound_coeff + 1;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;

    while found == 0 && ctr < equiv_num_iter {
        ctr += 1;
        // Select linear combination uniformly at random.
        ibz_rand_interval_minm_m(rng, &mut new_alpha.coord[0], equiv_bound_coeff);
        ibz_rand_interval_minm_m(rng, &mut new_alpha.coord[1], equiv_bound_coeff);
        ibz_rand_interval_minm_m(rng, &mut new_alpha.coord[2], equiv_bound_coeff);
        ibz_rand_interval_minm_m(rng, &mut new_alpha.coord[3], equiv_bound_coeff);

        // Compute the (Gram-form) norm of the sampled vector.
        quat_qf_eval(&mut tmp, &gram, &new_alpha.coord);

        // Compute the norm of the equivalent ideal.
        let tmp_clone = tmp.clone();
        ibz_div(&mut tmp, &mut remainder, &tmp_clone, &adjusted_norm);

        // Debug: remainder must be zero.
        debug_assert!(
            ibz_is_zero(&remainder) != 0,
            "quat_lideal_prime_norm_reduced_equivalent: adjusted-norm divisibility broken"
        );

        // Pseudo-primality test.
        if ibz_probab_prime(&tmp, primality_num_iter) != 0 {
            // Compute the generator via matrix-vector product.
            let coord_clone = new_alpha.coord.clone();
            ibz_mat_4x4_eval(&mut new_alpha.coord, &red, &coord_clone);
            new_alpha.denom = lideal.lattice.denom.clone();
            debug_assert!(
                quat_lattice_contains(None, &lideal.lattice, &new_alpha) != 0,
                "quat_lideal_prime_norm_reduced_equivalent: generator not in lattice"
            );

            let alpha_clone = new_alpha.clone();
            quat_alg_conj(&mut new_alpha, &alpha_clone);
            let denom_clone = new_alpha.denom.clone();
            ibz_mul(&mut new_alpha.denom, &denom_clone, &lideal.norm);
            let lideal_clone = lideal.clone();
            quat_lideal_mul(lideal, &lideal_clone, &new_alpha, alg);
            // Recompute the norm against the parent order; the reference
            // does this implicitly via the borrowed `parent_order`
            // pointer carried on the struct.
            quat_lideal_norm(lideal, parent_order);
            debug_assert!(
                ibz_probab_prime(&lideal.norm, primality_num_iter) != 0,
                "quat_lideal_prime_norm_reduced_equivalent: post-mul norm not prime"
            );

            found = 1;
            break;
        }
    }

    found
}
