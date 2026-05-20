//! Direct port of `the-sqisign/src/id2iso/ref/lvlx/dim2id2iso.c`.
//!
//! Three families of entry points:
//!
//! 1. [`find_uv`]: the Clapotis enumerator that picks the small-norm
//!    quaternion pair `(beta1, beta2)` and the integer coefficients
//!    `(u, v, d1, d2)` driving the dimension-two Kani diagram. This is
//!    deterministic given the input ideal; the helpers
//!    `post_LLL_basis_treatment`, `enumerate_hypercube`, and
//!    `find_uv_from_lists` live inside this module.
//!
//! 2. [`fixed_degree_isogeny_and_eval`]: builds an isogeny of fixed degree
//!    `u` and evaluates it at a list of points. Requires
//!    `quat_represent_integer` (RNG-driven), which is **not yet ported**.
//!    The Rust signature compiles but the body panics with a deferral
//!    notice when invoked. To remove the panic, port the RNG-driven
//!    norm-equation solvers in `sqisign-quaternion::normeq` first.
//!
//! 3. [`dim2id2iso_ideal_to_isogeny_clapotis`] and
//!    [`dim2id2iso_arbitrary_isogeny_evaluation`]: the public driver that
//!    composes `find_uv` with two `fixed_degree_isogeny_and_eval` calls
//!    and a dimension-two `theta_chain_compute_and_eval_randomized` step.
//!    Inherits the deferral above.

use sqisign_common::RngSource;
use sqisign_ec::{
    copy_basis, copy_curve, copy_point, ec_curve_normalize_a24, ec_dbl_iter_basis, ec_point_init,
    EcBasis, EcCurve, NWORDS_ORDER, TORSION_EVEN_POWER,
};
use sqisign_hd::{
    copy_bases_to_kernel, double_couple_point_iter, theta_chain_compute_and_eval,
    theta_chain_compute_and_eval_randomized, ThetaCoupleCurve, ThetaCouplePoint,
    ThetaKernelCouplePoints, HD_EXTRA_TORSION,
};
use sqisign_precomp::{
    CONNECTING_IDEALS, CURVES_WITH_ENDOMORPHISMS, EXTREMAL_ORDERS, NUM_ALTERNATE_EXTREMAL_ORDERS,
    QUATALG_PINFTY, TORSION_PLUS_2POWER,
};
use sqisign_quaternion::dim4::{
    ibz_mat_4x4_copy, ibz_mat_4x4_eval, ibz_mat_4x4_new, ibz_vec_4_new, quat_qf_eval, IbzMat4x4,
    IbzVec4,
};
use sqisign_quaternion::{
    ibz_add, ibz_bitsize, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero,
    ibz_cornacchia_prime, ibz_div, ibz_get, ibz_invmod, ibz_is_even, ibz_is_odd, ibz_is_zero,
    ibz_mod, ibz_mod_ui, ibz_mul, ibz_neg, ibz_pow, ibz_set, ibz_sub, ibz_two_adic, quat_alg_conj,
    quat_alg_mul, quat_alg_normalize, quat_lattice_alg_elem_mul, quat_lattice_contains,
    quat_lideal_conjugate_without_hnf, quat_lideal_copy, quat_lideal_create,
    quat_lideal_lideal_mul_reduced, quat_lideal_reduce_basis, quat_represent_integer, Ibz, QuatAlg,
    QuatAlgElem, QuatLattice, QuatLeftIdeal, QuatRepresentIntegerParams,
};

use crate::id2iso::endomorphism_application_even_basis;

/// Local zero constructor for [`Fp2`], mirroring the inline literal the
/// rest of `sqisign-hd` uses (`Fp2 { re: [0; NWORDS_FIELD], im: [0; NWORDS_FIELD] }`).
#[inline]
fn fp2_zero() -> sqisign_gf::Fp2 {
    sqisign_gf::Fp2 {
        re: [0u64; sqisign_gf::NWORDS_FIELD],
        im: [0u64; sqisign_gf::NWORDS_FIELD],
    }
}

/// `QUAT_repres_bound_input` (lvl1). Used by [`fixed_degree_isogeny_and_eval`]
/// to size the dimension-two isogeny step.
pub const QUAT_REPRES_BOUND_INPUT: i32 = 20;

/// `QUAT_primality_num_iter` (lvl1). Iteration count for the
/// Miller-Rabin probable-prime tests inside `quat_represent_integer`.
pub const QUAT_REPRESENT_INTEGER_PRIMALITY_ITER: i32 = 32;

/// `QUAT_equiv_bound_coeff` (lvl1). Bound passed to the
/// `quat_lideal_prime_norm_reduced_equivalent` search in keygen/sign.
pub const QUAT_EQUIV_BOUND_COEFF: i32 = 64;

/// `FINDUV_box_size` (lvl1).
pub const FINDUV_BOX_SIZE: i32 = 2;

/// `FINDUV_cube_size` (lvl1): pre-sized upper bound on the number of
/// short vectors used by [`enumerate_hypercube`].
pub const FINDUV_CUBE_SIZE: usize = 624;

// ----------------------------------------------------------------------
// Internal helpers (file-private in the C reference).
// ----------------------------------------------------------------------

/// Mirrors the static `post_LLL_basis_treatment` in the C reference. The
/// `is_special_order` flag triggers the special-order specific reordering
/// and sign-flipping pass.
fn post_lll_basis_treatment(gram: &mut IbzMat4x4, reduced: &mut IbzMat4x4, is_special_order: bool) {
    if !is_special_order {
        return;
    }

    if ibz_cmp(&gram[0][0], &gram[2][2]) == 0 {
        for i in 0..4 {
            reduced[i].swap(1, 2);
        }
        gram_swap(gram, (0, 2), (0, 1));
        gram_swap(gram, (2, 0), (1, 0));
        gram_swap(gram, (3, 2), (3, 1));
        gram_swap(gram, (2, 3), (1, 3));
        gram_swap(gram, (2, 2), (1, 1));
    } else if ibz_cmp(&gram[0][0], &gram[3][3]) == 0 {
        for i in 0..4 {
            reduced[i].swap(1, 3);
        }
        gram_swap(gram, (0, 3), (0, 1));
        gram_swap(gram, (3, 0), (1, 0));
        gram_swap(gram, (2, 3), (2, 1));
        gram_swap(gram, (3, 2), (1, 2));
        gram_swap(gram, (3, 3), (1, 1));
    } else if ibz_cmp(&gram[1][1], &gram[3][3]) == 0 {
        for i in 0..4 {
            reduced[i].swap(1, 2);
        }
        gram_swap(gram, (0, 2), (0, 1));
        gram_swap(gram, (2, 0), (1, 0));
        gram_swap(gram, (3, 2), (3, 1));
        gram_swap(gram, (2, 3), (1, 3));
        gram_swap(gram, (2, 2), (1, 1));
    }

    if ibz_cmp(&reduced[0][0], &reduced[1][1]) != 0 {
        for i in 0..4 {
            let r = reduced[i][1].clone();
            ibz_neg(&mut reduced[i][1], &r);
            let g = gram[i][1].clone();
            ibz_neg(&mut gram[i][1], &g);
            let g = gram[1][i].clone();
            ibz_neg(&mut gram[1][i], &g);
        }
    }
    if ibz_cmp(&reduced[0][2], &reduced[1][3]) != 0 {
        for i in 0..4 {
            let r = reduced[i][3].clone();
            ibz_neg(&mut reduced[i][3], &r);
            let g = gram[i][3].clone();
            ibz_neg(&mut gram[i][3], &g);
            let g = gram[3][i].clone();
            ibz_neg(&mut gram[3][i], &g);
        }
    }
}

fn gram_swap(m: &mut IbzMat4x4, a: (usize, usize), b: (usize, usize)) {
    let tmp = m[a.0][a.1].clone();
    m[a.0][a.1] = m[b.0][b.1].clone();
    m[b.0][b.1] = tmp;
}

/// Mirrors the static `enumerate_hypercube` helper. Enumerates non-trivial
/// integer hypercube points of `[-m, m]^4` whose quadratic form divides
/// the adjusted norm to an odd quotient. Returns `count - 1` to match the
/// C reference (the caller treats `0` as "no vectors found").
fn enumerate_hypercube(
    vecs: &mut [IbzVec4],
    norms: &mut [Ibz],
    m: i32,
    gram: &IbzMat4x4,
    adjusted_norm: &Ibz,
) -> i32 {
    let mut remain = Ibz::zero();
    let mut norm = Ibz::zero();
    let mut point: IbzVec4 = ibz_vec_4_new();

    assert!(m > 0);

    let mut count: usize = 0;
    let dim = 2 * m + 1;
    let dim2 = dim * dim;
    let dim3 = dim2 * dim;

    let need_remove_symmetry =
        ibz_cmp(&gram[0][0], &gram[1][1]) == 0 && ibz_cmp(&gram[3][3], &gram[2][2]) == 0;

    for x in -m..=0 {
        for y in -m..=m {
            if x == 0 && y > 0 {
                break;
            }
            for z in -m..=m {
                if x == 0 && y == 0 && z > 0 {
                    break;
                }
                for w in -m..=m {
                    if x == 0 && y == 0 && z == 0 && w >= 0 {
                        break;
                    }
                    if ((x | y | z | w) & 1) == 0 {
                        continue;
                    }
                    if x % 3 == 0 && y % 3 == 0 && z % 3 == 0 && w % 3 == 0 {
                        continue;
                    }
                    let check1 = (m + w) + dim * (m + z) + dim2 * (m + y) + dim3 * (m + x);
                    let check2 = (m - z) + dim * (m + w) + dim2 * (m - x) + dim3 * (m + y);
                    let check3 = (m + z) + dim * (m - w) + dim2 * (m + x) + dim3 * (m - y);

                    if !need_remove_symmetry || (check1 <= check2 && check1 <= check3) {
                        ibz_set(&mut point[0], x);
                        ibz_set(&mut point[1], y);
                        ibz_set(&mut point[2], z);
                        ibz_set(&mut point[3], w);

                        quat_qf_eval(&mut norm, gram, &point);
                        let saved = norm.clone();
                        ibz_div(&mut norm, &mut remain, &saved, adjusted_norm);
                        debug_assert!(ibz_is_zero(&remain) != 0);

                        if ibz_mod_ui(&norm, 2) == 1 {
                            ibz_set(&mut vecs[count][0], x);
                            ibz_set(&mut vecs[count][1], y);
                            ibz_set(&mut vecs[count][2], z);
                            ibz_set(&mut vecs[count][3], w);
                            norms[count] = norm.clone();
                            count += 1;
                        }
                    }
                }
            }
        }
    }

    (count as i32) - 1
}

/// Mirrors the static `find_uv_from_lists` helper.
#[allow(clippy::too_many_arguments)]
fn find_uv_from_lists(
    au: &mut Ibz,
    bu: &mut Ibz,
    av: &mut Ibz,
    bv: &mut Ibz,
    u: &mut Ibz,
    v: &mut Ibz,
    index_sol1: &mut i32,
    index_sol2: &mut i32,
    target: &Ibz,
    small_norms1: &[Ibz],
    small_norms2: &[Ibz],
    quotients: &[Ibz],
    index1: i32,
    index2: i32,
    is_diagonal: bool,
    number_sum_square: i32,
) -> i32 {
    let mut n = target.clone();
    let mut remain = Ibz::zero();
    let mut adjusted_norm = Ibz::zero();
    let mut found = 0i32;
    let mut cmp;

    for i1 in 0..(index1 as usize) {
        ibz_mod(&mut adjusted_norm, &n, &small_norms1[i1]);
        let starting_index2 = if is_diagonal { i1 as i32 } else { 0 };
        for i2 in (starting_index2 as usize)..(index2 as usize) {
            if ibz_invmod(&mut remain, &small_norms2[i2], &small_norms1[i1]) == 0 {
                continue;
            }
            ibz_mul(v, &remain, &adjusted_norm);
            let v_saved = v.clone();
            ibz_mod(v, &v_saved, &small_norms1[i1]);
            cmp = ibz_cmp(v, &quotients[i2]);
            while found == 0 && cmp < 0 {
                if number_sum_square > 0 {
                    found = ibz_cornacchia_prime(av, bv, &ibz_const_one(), v);
                } else if number_sum_square == 0 {
                    found = 1;
                }
                if found != 0 {
                    ibz_mul(&mut remain, v, &small_norms2[i2]);
                    *au = n.clone();
                    ibz_sub(u, au, &remain);
                    assert!(ibz_cmp(u, &ibz_const_zero()) > 0);
                    let u_saved = u.clone();
                    ibz_div(u, &mut remain, &u_saved, &small_norms1[i1]);
                    assert!(ibz_is_zero(&remain) != 0);
                    let nonzero_u = ibz_get(u) != 0;
                    let nonzero_v = ibz_get(v) != 0;
                    if !(nonzero_u && nonzero_v) {
                        found = 0;
                    }
                    if number_sum_square == 2 {
                        found = ibz_cornacchia_prime(au, bu, &ibz_const_one(), u);
                    }
                }
                if found == 0 {
                    let v_saved = v.clone();
                    ibz_add(v, &v_saved, &small_norms1[i1]);
                    cmp = ibz_cmp(v, &quotients[i2]);
                }
            }
            if found != 0 {
                *index_sol1 = i1 as i32;
                *index_sol2 = i2 as i32;
                break;
            }
        }
        if found != 0 {
            // Update n to satisfy the assertion below; n was unchanged in C.
            let _ = &mut n;
            break;
        }
    }

    found
}

// ----------------------------------------------------------------------
// Public entry points.
// ----------------------------------------------------------------------

/// `find_uv(u, v, beta1, beta2, d1, d2, idx1, idx2, target, lideal,
/// Bpoo, num_alternate_order)`. Mirrors the C reference verbatim.
///
/// Returns 1 on success, 0 on failure (no suitable `(u, v, d1, d2)`
/// found across the alternate orders).
#[allow(clippy::too_many_arguments)]
pub fn find_uv(
    u: &mut Ibz,
    v: &mut Ibz,
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    d1: &mut Ibz,
    d2: &mut Ibz,
    index_alternate_order_1: &mut i32,
    index_alternate_order_2: &mut i32,
    target: &Ibz,
    lideal: &QuatLeftIdeal,
    bpoo: &QuatAlg,
    num_alternate_order: i32,
) -> i32 {
    let mut au = Ibz::zero();
    let mut bu = Ibz::zero();
    let mut av = Ibz::zero();
    let mut bv = Ibz::zero();
    let mut norm_d = Ibz::zero();
    let mut remain = Ibz::zero();
    let mut n = target.clone();

    let no = (num_alternate_order + 1) as usize;
    let mut adjusted_norm: Vec<Ibz> = (0..no).map(|_| Ibz::zero()).collect();
    let mut gram: Vec<IbzMat4x4> = (0..no).map(|_| ibz_mat_4x4_new()).collect();
    let mut reduced: Vec<IbzMat4x4> = (0..no).map(|_| ibz_mat_4x4_new()).collect();
    let mut ideal: Vec<QuatLeftIdeal> = (0..no).map(|_| QuatLeftIdeal::new()).collect();

    // Reduce the input ideal in place at index 0.
    quat_lideal_copy(&mut ideal[0], lideal);
    quat_lideal_reduce_basis(&mut reduced[0], &mut gram[0], &ideal[0], bpoo);
    ibz_mat_4x4_copy(&mut ideal[0].lattice.basis, &reduced[0]);
    ibz_set(&mut adjusted_norm[0], 1);
    let an = adjusted_norm[0].clone();
    ibz_mul(&mut adjusted_norm[0], &an, &ideal[0].lattice.denom);
    let an = adjusted_norm[0].clone();
    ibz_mul(&mut adjusted_norm[0], &an, &ideal[0].lattice.denom);
    {
        let (g_split, r_split) = (&mut gram[0], &mut reduced[0]);
        post_lll_basis_treatment(g_split, r_split, true);
    }

    // reduced_id = ideal[0] * \overline{delta} / n(ideal[0])
    let mut reduced_id = QuatLeftIdeal::new();
    quat_lideal_copy(&mut reduced_id, &ideal[0]);
    let mut delta = QuatAlgElem::new();
    ibz_set(&mut delta.coord[0], 1);
    ibz_set(&mut delta.coord[1], 0);
    ibz_set(&mut delta.coord[2], 0);
    ibz_set(&mut delta.coord[3], 0);
    delta.denom = reduced_id.lattice.denom.clone();
    let saved = delta.coord.clone();
    ibz_mat_4x4_eval(&mut delta.coord, &reduced[0], &saved);
    debug_assert!(quat_lattice_contains(None, &reduced_id.lattice, &delta) != 0);

    let delta_in = delta.clone();
    quat_alg_conj(&mut delta, &delta_in);
    let saved = delta.denom.clone();
    ibz_mul(&mut delta.denom, &saved, &ideal[0].norm);
    let lat = reduced_id.lattice.clone();
    quat_lattice_alg_elem_mul(&mut reduced_id.lattice, &lat, &delta, bpoo);
    reduced_id.norm = gram[0][0][0].clone();
    let norm_saved = reduced_id.norm.clone();
    ibz_div(
        &mut reduced_id.norm,
        &mut remain,
        &norm_saved,
        &adjusted_norm[0],
    );
    debug_assert!(ibz_cmp(&remain, &ibz_const_zero()) == 0);

    // conj_ideal is the conjugate of reduced_id; we also produce the new
    // parent (right) order.
    let mut right_order = QuatLattice::new();
    let mut conj_ideal = QuatLeftIdeal::new();
    quat_lideal_conjugate_without_hnf(&mut conj_ideal, &mut right_order, &reduced_id, bpoo);

    // Build the alternate ideals.
    let maxord_o0 = &EXTREMAL_ORDERS[0].order;
    for i in 1..no {
        // ALTERNATE_CONNECTING_IDEALS[i - 1] = CONNECTING_IDEALS[i]
        let alt_id = &CONNECTING_IDEALS[i];
        // The Rust port of `quat_lideal_lideal_mul_reduced` expects an
        // explicit `parent_order`; lideal1 is `conj_ideal`, with parent
        // order `right_order` in the C reference.
        quat_lideal_lideal_mul_reduced(
            &mut ideal[i],
            &mut gram[i],
            &conj_ideal,
            alt_id,
            maxord_o0,
            bpoo,
        );
        let r = ideal[i].lattice.basis.clone();
        ibz_mat_4x4_copy(&mut reduced[i], &r);
        ibz_set(&mut adjusted_norm[i], 1);
        let an = adjusted_norm[i].clone();
        ibz_mul(&mut adjusted_norm[i], &an, &ideal[i].lattice.denom);
        let an = adjusted_norm[i].clone();
        ibz_mul(&mut adjusted_norm[i], &an, &ideal[i].lattice.denom);
        let (g_split, r_split) = (&mut gram[i], &mut reduced[i]);
        post_lll_basis_treatment(g_split, r_split, false);
    }

    // Short-vector enumeration parameters (lvl-specific).
    let m = FINDUV_BOX_SIZE;
    let m4 = FINDUV_CUBE_SIZE;

    let mut small_vecs: Vec<Vec<IbzVec4>> = (0..no)
        .map(|_| (0..m4).map(|_| ibz_vec_4_new()).collect())
        .collect();
    let mut small_norms: Vec<Vec<Ibz>> = (0..no)
        .map(|_| (0..m4).map(|_| Ibz::zero()).collect())
        .collect();
    let mut quotients: Vec<Vec<Ibz>> = (0..no)
        .map(|_| (0..m4).map(|_| Ibz::zero()).collect())
        .collect();
    let mut indices: Vec<i32> = vec![0; no];

    for j in 0..no {
        indices[j] = enumerate_hypercube(
            &mut small_vecs[j],
            &mut small_norms[j],
            m,
            &gram[j],
            &adjusted_norm[j],
        );

        // Sort the (vec, norm) pairs by norm, breaking ties by original index.
        let n_used = indices[j] as usize;
        if n_used > 0 {
            let mut tagged: Vec<(IbzVec4, Ibz, usize)> = (0..n_used)
                .map(|i| (small_vecs[j][i].clone(), small_norms[j][i].clone(), i))
                .collect();
            tagged.sort_by(|a, b| {
                let c = ibz_cmp(&a.1, &b.1);
                if c != 0 {
                    c.cmp(&0)
                } else {
                    a.2.cmp(&b.2)
                }
            });
            for (i, (v, n_, _)) in tagged.into_iter().enumerate() {
                small_vecs[j][i] = v;
                small_norms[j][i] = n_;
            }
        }

        for i in 0..(indices[j] as usize) {
            let mut r = Ibz::zero();
            ibz_div(&mut quotients[j][i], &mut r, &n, &small_norms[j][i]);
        }
    }

    let mut found = 0i32;
    let mut i1: i32 = 0;
    let mut i2: i32 = 0;
    let mut chosen_j1 = 0usize;
    let mut chosen_j2 = 0usize;

    'outer: for j1 in 0..no {
        for j2 in j1..no {
            let is_diago = j1 == j2;
            found = find_uv_from_lists(
                &mut au,
                &mut bu,
                &mut av,
                &mut bv,
                u,
                v,
                &mut i1,
                &mut i2,
                target,
                &small_norms[j1],
                &small_norms[j2],
                &quotients[j2],
                indices[j1],
                indices[j2],
                is_diago,
                0,
            );

            if found != 0 {
                beta1.denom = ideal[j1].lattice.denom.clone();
                beta2.denom = ideal[j2].lattice.denom.clone();
                *d1 = small_norms[j1][i1 as usize].clone();
                *d2 = small_norms[j2][i2 as usize].clone();
                ibz_mat_4x4_eval(&mut beta1.coord, &reduced[j1], &small_vecs[j1][i1 as usize]);
                ibz_mat_4x4_eval(&mut beta2.coord, &reduced[j2], &small_vecs[j2][i2 as usize]);

                if j1 != 0 || j2 != 0 {
                    let saved = delta.denom.clone();
                    ibz_div(&mut delta.denom, &mut remain, &saved, &lideal.norm);
                    debug_assert!(ibz_cmp(&remain, &ibz_const_zero()) == 0);
                    let saved = delta.denom.clone();
                    ibz_mul(&mut delta.denom, &saved, &conj_ideal.norm);
                }
                if j1 != 0 {
                    let saved = beta1.clone();
                    quat_alg_mul(beta1, &delta, &saved, bpoo);
                    quat_alg_normalize(beta1);
                }
                if j2 != 0 {
                    let saved = beta2.clone();
                    quat_alg_mul(beta2, &delta, &saved, bpoo);
                    quat_alg_normalize(beta2);
                }
                if j1 != 0 {
                    let saved = beta1.clone();
                    quat_alg_conj(beta1, &saved);
                }
                if j2 != 0 {
                    let saved = beta2.clone();
                    quat_alg_conj(beta2, &saved);
                }

                chosen_j1 = j1;
                chosen_j2 = j2;
                break 'outer;
            }
        }
    }

    if found != 0 {
        *index_alternate_order_1 = chosen_j1 as i32;
        *index_alternate_order_2 = chosen_j2 as i32;
    }

    // Touch unused variables to silence unused-mut warnings on the rare
    // path where the function fails (matches the C reference's cleanup).
    let _ = (&mut n, &mut norm_d, &mut au, &mut bu, &mut av, &mut bv);

    found
}

/// `fixed_degree_isogeny_and_eval(rng, lideal, u, small, E34, P12, numP,
/// index_alternate_order)`. Mirrors the C entry point.
///
/// Builds an isogeny of fixed degree `u` from an alternate maximal
/// order and evaluates it at the points of `p12`. Returns the chain
/// length on success and 0 if the underlying `quat_represent_integer`
/// search exhausts its bounded counter.
///
/// The `rng` argument replaces the reference's thread-local DRBG global;
/// the byte stream consumed is identical for an identically seeded DRBG.
#[allow(clippy::too_many_arguments)]
pub fn fixed_degree_isogeny_and_eval<R: RngSource>(
    rng: &mut R,
    lideal: &mut QuatLeftIdeal,
    u: &Ibz,
    small: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    num_p: usize,
    index_alternate_order: i32,
) -> i32 {
    fixed_degree_isogeny_impl(
        rng,
        lideal,
        u,
        small,
        e34,
        p12,
        num_p,
        index_alternate_order,
    )
}

/// Internal worker for [`fixed_degree_isogeny_and_eval`]. Mirrors the
/// static `_fixed_degree_isogeny_impl` in
/// `the-sqisign/src/id2iso/ref/lvlx/dim2id2iso.c`.
#[allow(clippy::too_many_arguments)]
fn fixed_degree_isogeny_impl<R: RngSource>(
    rng: &mut R,
    lideal: &mut QuatLeftIdeal,
    u: &Ibz,
    small: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    num_p: usize,
    index_alternate_order: i32,
) -> i32 {
    let mut e0 = EcCurve::zero();
    copy_curve(
        &mut e0,
        &CURVES_WITH_ENDOMORPHISMS[index_alternate_order as usize].curve,
    );
    ec_curve_normalize_a24(&mut e0);

    let u_bitsize = ibz_bitsize(u);

    // Decide the dimension-two step's power-of-two length. Smaller is
    // faster but risks `quat_represent_integer` failing the search.
    let length: u32 = if !small {
        (TORSION_EVEN_POWER as u32).saturating_sub(HD_EXTRA_TORSION)
    } else {
        let l = ibz_bitsize(&QUATALG_PINFTY.p) + QUAT_REPRES_BOUND_INPUT - u_bitsize;
        debug_assert!(
            u_bitsize < l,
            "fixed_degree_isogeny_and_eval: bitsize bound"
        );
        debug_assert!(
            (l as u32) < (TORSION_EVEN_POWER as u32).saturating_sub(HD_EXTRA_TORSION),
            "fixed_degree_isogeny_and_eval: length under torsion bound"
        );
        l as u32
    };
    debug_assert!(length > 0, "fixed_degree_isogeny_and_eval: length is zero");

    // theta = quat_represent_integer with target norm u * (2^length - u).
    let mut two_pow = Ibz::zero();
    ibz_pow(&mut two_pow, &ibz_const_two(), length);
    let mut tmp = u.clone();
    debug_assert!(ibz_cmp(&two_pow, &tmp) > 0);
    debug_assert!(ibz_is_even(&tmp) == 0, "u must be odd");

    let cur = tmp.clone();
    ibz_sub(&mut tmp, &two_pow, &cur);
    let cur = tmp.clone();
    ibz_mul(&mut tmp, &cur, u);
    debug_assert!(ibz_is_even(&tmp) == 0, "target norm must be odd");

    let extremal = &EXTREMAL_ORDERS[index_alternate_order as usize];
    let ri_params = QuatRepresentIntegerParams {
        primality_test_iterations: QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
        order: extremal,
        algebra: &QUATALG_PINFTY,
    };
    let mut theta = QuatAlgElem::new();
    let ret_ri = quat_represent_integer(rng, &mut theta, &tmp, 1, &ri_params);
    if ret_ri == 0 {
        return 0;
    }

    quat_lideal_create(lideal, &theta, u, &extremal.order, &QUATALG_PINFTY);

    // Pull down the precomputed even-torsion basis to order length+HD_extra.
    let mut b0_two = EcBasis::zero();
    copy_basis(
        &mut b0_two,
        &CURVES_WITH_ENDOMORPHISMS[index_alternate_order as usize].basis_even,
    );
    let drop_iters = (TORSION_EVEN_POWER as i32) - (length as i32) - (HD_EXTRA_TORSION as i32);
    let b0_two_clone = b0_two;
    ec_dbl_iter_basis(&mut b0_two, drop_iters, &b0_two_clone, &mut e0);

    // theta *= u^{-1} mod 2^(length + 2).
    let cur = two_pow.clone();
    ibz_mul(&mut two_pow, &cur, &ibz_const_two());
    let cur = two_pow.clone();
    ibz_mul(&mut two_pow, &cur, &ibz_const_two());
    let mut inv_u = u.clone();
    let _ok_inv = ibz_invmod(&mut inv_u, &u.clone(), &two_pow);
    debug_assert!(ibz_is_even(&inv_u) == 0, "u^{{-1}} mod 2^(L+2) must be odd");
    let inv_u_clone = inv_u.clone();
    for k in 0..4 {
        let cur = theta.coord[k].clone();
        ibz_mul(&mut theta.coord[k], &cur, &inv_u_clone);
    }

    // Apply theta to the basis; record the image as the dim-two kernel
    // second basis.
    let mut b0_two_theta = EcBasis::zero();
    copy_basis(&mut b0_two_theta, &b0_two);
    endomorphism_application_even_basis(
        &mut b0_two_theta,
        index_alternate_order,
        &e0,
        &theta,
        (length + HD_EXTRA_TORSION) as i32,
    );

    // Domain: E0 x E0.
    let mut e00 = ThetaCoupleCurve::zero();
    e00.e1 = e0;
    e00.e2 = e0;

    let mut dim_two_ker = ThetaKernelCouplePoints::zero();
    copy_bases_to_kernel(&mut dim_two_ker, &b0_two, &b0_two_theta);

    let ret_chain =
        theta_chain_compute_and_eval(length, &mut e00, &dim_two_ker, true, e34, &mut p12[..num_p]);
    if ret_chain == 0 {
        return 0;
    }
    length as i32
}

/// `dim2id2iso_ideal_to_isogeny_clapotis(rng, beta1, beta2, u, v, d1, d2,
/// codomain, basis, lideal, Bpoo)`. Mirrors the C entry point.
///
/// Composes [`find_uv`] (deterministic) with two
/// [`fixed_degree_isogeny_and_eval`] calls and one
/// [`theta_chain_compute_and_eval_randomized`] step. All RNG-derived
/// bytes flow through `rng`; the byte stream consumed matches what the
/// C reference would have drawn under an identically seeded DRBG.
#[allow(clippy::too_many_arguments)]
pub fn dim2id2iso_ideal_to_isogeny_clapotis<R: RngSource>(
    rng: &mut R,
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    u: &mut Ibz,
    v: &mut Ibz,
    d1: &mut Ibz,
    d2: &mut Ibz,
    codomain: &mut EcCurve,
    basis: &mut EcBasis,
    lideal: &QuatLeftIdeal,
    bpoo: &QuatAlg,
) -> i32 {
    let mut target = Ibz::zero();
    let mut tmp = Ibz::zero();
    let mut two_pow = Ibz::zero();
    let mut theta = QuatAlgElem::new();
    let _ = (&mut target, &mut two_pow); // (mirror C scaffolding; some scratch is local)

    // 1. Pick u, v, d1, d2, beta1, beta2 and the alternate-order indices.
    let mut index_order1: i32 = 0;
    let mut index_order2: i32 = 0;
    let ret_uv = find_uv(
        u,
        v,
        beta1,
        beta2,
        d1,
        d2,
        &mut index_order1,
        &mut index_order2,
        &TORSION_PLUS_2POWER,
        lideal,
        bpoo,
        NUM_ALTERNATE_EXTREMAL_ORDERS as i32,
    );
    if ret_uv == 0 {
        return 0;
    }
    debug_assert!(
        ibz_is_odd(d1) != 0 && ibz_is_odd(d2) != 0,
        "dim2id2iso_ideal_to_isogeny_clapotis: d1, d2 must be odd"
    );

    // 2. Strip the power-of-two GCD of (u, v).
    let mut gcd_uv = Ibz::zero();
    sqisign_quaternion::ibz_gcd(&mut gcd_uv, u, v);
    debug_assert!(ibz_cmp(&gcd_uv, &ibz_const_zero()) != 0);
    let exp_gcd = ibz_two_adic(&gcd_uv);
    let exp: i32 = (TORSION_EVEN_POWER as i32) - exp_gcd;
    let mut rem = Ibz::zero();
    let u_clone = u.clone();
    ibz_div(u, &mut rem, &u_clone, &gcd_uv);
    debug_assert!(ibz_is_zero(&rem) != 0);
    let v_clone = v.clone();
    ibz_div(v, &mut rem, &v_clone, &gcd_uv);
    debug_assert!(ibz_is_zero(&rem) != 0);

    // 3. Pull down the precomputed bases at the two alternate orders.
    let mut e1 = EcCurve::zero();
    let mut e2 = EcCurve::zero();
    copy_curve(
        &mut e1,
        &CURVES_WITH_ENDOMORPHISMS[index_order1 as usize].curve,
    );
    copy_curve(
        &mut e2,
        &CURVES_WITH_ENDOMORPHISMS[index_order2 as usize].curve,
    );
    let mut bas1 = EcBasis::zero();
    let mut bas2 = EcBasis::zero();
    copy_basis(
        &mut bas1,
        &CURVES_WITH_ENDOMORPHISMS[index_order1 as usize].basis_even,
    );
    copy_basis(
        &mut bas2,
        &CURVES_WITH_ENDOMORPHISMS[index_order2 as usize].basis_even,
    );

    // 4. theta = beta2 * conj(beta1) / n(lideal).
    ibz_set(&mut theta.denom, 1);
    let beta1_clone = beta1.clone();
    quat_alg_conj(&mut theta, &beta1_clone);
    let theta_clone = theta.clone();
    quat_alg_mul(&mut theta, beta2, &theta_clone, &QUATALG_PINFTY);
    let cur = theta.denom.clone();
    ibz_mul(&mut theta.denom, &cur, &lideal.norm);

    let mut idealu = QuatLeftIdeal::new();
    let mut idealv = QuatLeftIdeal::new();
    let mut fu_codomain = ThetaCoupleCurve::zero();
    let mut fv_codomain = ThetaCoupleCurve::zero();

    // pushed_points := { P, Q, P-Q } on the dim-2 product.
    let mut pushed_points: [ThetaCouplePoint; 3] = [
        ThetaCouplePoint::zero(),
        ThetaCouplePoint::zero(),
        ThetaCouplePoint::zero(),
    ];
    copy_point(&mut pushed_points[0].p1, &bas1.P);
    copy_point(&mut pushed_points[1].p1, &bas1.Q);
    copy_point(&mut pushed_points[2].p1, &bas1.PmQ);
    ec_point_init(&mut pushed_points[0].p2);
    ec_point_init(&mut pushed_points[1].p2);
    ec_point_init(&mut pushed_points[2].p2);

    // 5. phi_u: fixed-degree isogeny from index_order1 of degree u.
    let ret_u = fixed_degree_isogeny_and_eval(
        rng,
        &mut idealu,
        u,
        true,
        &mut fu_codomain,
        &mut pushed_points,
        3,
        index_order1,
    );
    if ret_u == 0 {
        return 0;
    }

    // Capture phi_u(bas1) as the new bas_u on Fu_codomain.E1.
    let mut bas_u = EcBasis::zero();
    copy_point(&mut bas_u.P, &pushed_points[0].p1);
    copy_point(&mut bas_u.Q, &pushed_points[1].p1);
    copy_point(&mut bas_u.PmQ, &pushed_points[2].p1);

    // Half of the dim-2 kernel: (phi_u(bas1), 0) on Fu_codomain.E1.
    let mut ker = ThetaKernelCouplePoints::zero();
    copy_point(&mut ker.t1.p1, &bas_u.P);
    copy_point(&mut ker.t2.p1, &bas_u.Q);
    copy_point(&mut ker.t1m2.p1, &bas_u.PmQ);
    let mut e01 = ThetaCoupleCurve::zero();
    copy_curve(&mut e01.e1, &fu_codomain.e1);

    // Reset pushed_points to (bas2, 0); call phi_v at degree v.
    copy_point(&mut pushed_points[0].p1, &bas2.P);
    copy_point(&mut pushed_points[1].p1, &bas2.Q);
    copy_point(&mut pushed_points[2].p1, &bas2.PmQ);
    ec_point_init(&mut pushed_points[0].p2);
    ec_point_init(&mut pushed_points[1].p2);
    ec_point_init(&mut pushed_points[2].p2);
    let ret_v = fixed_degree_isogeny_and_eval(
        rng,
        &mut idealv,
        v,
        true,
        &mut fv_codomain,
        &mut pushed_points,
        3,
        index_order2,
    );
    if ret_v == 0 {
        return 0;
    }

    // bas2 := phi_v(bas2) on Fv_codomain.E1.
    copy_point(&mut bas2.P, &pushed_points[0].p1);
    copy_point(&mut bas2.Q, &pushed_points[1].p1);
    copy_point(&mut bas2.PmQ, &pushed_points[2].p1);

    // 6. theta *= 1 / (d1 * n(connecting_ideal2)) mod 2^TORSION_EVEN_POWER.
    ibz_pow(&mut two_pow, &ibz_const_two(), TORSION_EVEN_POWER as u32);
    let mut inv_factor = d1.clone();
    if index_order2 > 0 {
        let cur = inv_factor.clone();
        ibz_mul(
            &mut inv_factor,
            &cur,
            &CONNECTING_IDEALS[index_order2 as usize].norm,
        );
    }
    let mut inv_tmp = Ibz::zero();
    let _ = ibz_invmod(&mut inv_tmp, &inv_factor, &two_pow);
    for k in 0..4 {
        let cur = theta.coord[k].clone();
        ibz_mul(&mut theta.coord[k], &cur, &inv_tmp);
    }

    // Apply theta to bas2 on Fv_codomain.E1.
    endomorphism_application_even_basis(
        &mut bas2,
        0,
        &fv_codomain.e1,
        &theta,
        TORSION_EVEN_POWER as i32,
    );

    // Second half of the dim-2 kernel: (0, theta(phi_v(bas2))) on
    // Fv_codomain.E1; set E01.E2 to that curve.
    copy_point(&mut ker.t1.p2, &bas2.P);
    copy_point(&mut ker.t2.p2, &bas2.Q);
    copy_point(&mut ker.t1m2.p2, &bas2.PmQ);
    copy_curve(&mut e01.e2, &fv_codomain.e1);

    // 7. Drop the kernel down to order 2^exp.
    let dim2_drop = (TORSION_EVEN_POWER as i32 - exp) as u32;
    let mut tmp_pt;
    tmp_pt = ker.t1;
    double_couple_point_iter(&mut ker.t1, dim2_drop, &tmp_pt, &e01);
    tmp_pt = ker.t2;
    double_couple_point_iter(&mut ker.t2, dim2_drop, &tmp_pt, &e01);
    tmp_pt = ker.t1m2;
    double_couple_point_iter(&mut ker.t1m2, dim2_drop, &tmp_pt, &e01);

    debug_assert!(ibz_is_odd(u) != 0, "u must remain odd after gcd strip");

    // 8. Evaluate (phi_u(bas1), 0) through the dim-2 isogeny of degree u*d1.
    copy_point(&mut pushed_points[0].p1, &bas_u.P);
    copy_point(&mut pushed_points[1].p1, &bas_u.Q);
    copy_point(&mut pushed_points[2].p1, &bas_u.PmQ);
    ec_point_init(&mut pushed_points[0].p2);
    ec_point_init(&mut pushed_points[1].p2);
    ec_point_init(&mut pushed_points[2].p2);

    let mut theta_codomain = ThetaCoupleCurve::zero();
    let ret_theta = theta_chain_compute_and_eval_randomized(
        rng,
        exp as u32,
        &mut e01,
        &ker,
        false,
        &mut theta_codomain,
        &mut pushed_points,
    );
    if ret_theta == 0 {
        return 0;
    }
    let t1_out = pushed_points[0];
    let t2_out = pushed_points[1];
    let t1m2_out = pushed_points[2];

    // 9. Select the correct codomain curve via a Weil-pairing check.
    copy_point(&mut basis.P, &t1_out.p1);
    copy_point(&mut basis.Q, &t2_out.p1);
    copy_point(&mut basis.PmQ, &t1m2_out.p1);
    copy_curve(codomain, &theta_codomain.e1);

    let mut w0 = fp2_zero();
    let mut w1 = fp2_zero();
    sqisign_ec::weil(
        &mut w0,
        TORSION_EVEN_POWER as u32,
        &bas1.P,
        &bas1.Q,
        &bas1.PmQ,
        &mut e1,
    );
    {
        let mut codomain_tmp = *codomain;
        sqisign_ec::weil(
            &mut w1,
            TORSION_EVEN_POWER as u32,
            &basis.P,
            &basis.Q,
            &basis.PmQ,
            &mut codomain_tmp,
        );
    }
    // (d1 * u * u) mod 2^TORSION_EVEN_POWER.
    let mut digit_target = d1.clone();
    let cur = digit_target.clone();
    ibz_mul(&mut digit_target, &cur, u);
    let cur = digit_target.clone();
    ibz_mul(&mut digit_target, &cur, u);
    let cur = digit_target.clone();
    ibz_mod(&mut digit_target, &cur, &TORSION_PLUS_2POWER);
    let mut digit_u = vec![0u64; NWORDS_ORDER];
    sqisign_quaternion::ibz_to_digits(&mut digit_u, &digit_target);
    let mut test_pow = fp2_zero();
    sqisign_gf::fp2_pow_vartime(&mut test_pow, &w0, &digit_u);
    if sqisign_gf::fp2_is_equal(&w1, &test_pow) == 0 {
        copy_point(&mut basis.P, &t1_out.p2);
        copy_point(&mut basis.Q, &t2_out.p2);
        copy_point(&mut basis.PmQ, &t1m2_out.p2);
        copy_curve(codomain, &theta_codomain.e2);
    }

    // 10. Apply beta1 scaled by 1 / (u * d1 * n(connecting_ideal1)) mod 2^TORSION_EVEN_POWER.
    ibz_mul(&mut tmp, u, d1);
    if index_order1 != 0 {
        let cur = tmp.clone();
        ibz_mul(
            &mut tmp,
            &cur,
            &CONNECTING_IDEALS[index_order1 as usize].norm,
        );
    }
    let mut inv_ud1 = Ibz::zero();
    let _ = ibz_invmod(&mut inv_ud1, &tmp, &TORSION_PLUS_2POWER);
    for k in 0..4 {
        let cur = beta1.coord[k].clone();
        ibz_mul(&mut beta1.coord[k], &cur, &inv_ud1);
    }
    endomorphism_application_even_basis(basis, 0, codomain, beta1, TORSION_EVEN_POWER as i32);

    1
}

/// `dim2id2iso_arbitrary_isogeny_evaluation(rng, basis, codomain, lideal)`.
/// Mirrors the C entry point: wraps
/// [`dim2id2iso_ideal_to_isogeny_clapotis`] with fresh scratch buffers
/// and discards the intermediate `(beta1, beta2, u, v, d1, d2)`.
pub fn dim2id2iso_arbitrary_isogeny_evaluation<R: RngSource>(
    rng: &mut R,
    basis: &mut EcBasis,
    codomain: &mut EcCurve,
    lideal: &QuatLeftIdeal,
) -> i32 {
    let mut beta1 = QuatAlgElem::new();
    let mut beta2 = QuatAlgElem::new();
    let mut u = Ibz::zero();
    let mut v = Ibz::zero();
    let mut d1 = Ibz::zero();
    let mut d2 = Ibz::zero();
    dim2id2iso_ideal_to_isogeny_clapotis(
        rng,
        &mut beta1,
        &mut beta2,
        &mut u,
        &mut v,
        &mut d1,
        &mut d2,
        codomain,
        basis,
        lideal,
        &QUATALG_PINFTY,
    )
}
