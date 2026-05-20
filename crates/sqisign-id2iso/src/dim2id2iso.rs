//! Direct port of `vendor/the-sqisign/src/id2iso/ref/lvlx/dim2id2iso.c`.
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

use sqisign_ec::{
    copy_basis, copy_curve, copy_point, ec_curve_normalize_a24, ec_dbl_iter_basis, ec_point_init,
    EcBasis, EcCurve, NWORDS_ORDER, TORSION_EVEN_POWER,
};
use sqisign_hd::{
    copy_bases_to_kernel, double_couple_point_iter, theta_chain_compute_and_eval_randomized,
    ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints, HD_EXTRA_TORSION,
};
use sqisign_precomp::{
    CONNECTING_IDEALS, CURVES_WITH_ENDOMORPHISMS, EXTREMAL_ORDERS, NUM_ALTERNATE_EXTREMAL_ORDERS,
    QUATALG_PINFTY, TORSION_PLUS_2POWER,
};
use sqisign_quaternion::dim2::ibz_mat_2x2_new;
use sqisign_quaternion::dim4::{
    ibz_mat_4x4_copy, ibz_mat_4x4_eval, ibz_mat_4x4_new, ibz_vec_4_new, quat_qf_eval, IbzMat4x4,
    IbzVec4,
};
use sqisign_quaternion::{
    ibz_add, ibz_bitsize, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero,
    ibz_cornacchia_prime, ibz_div, ibz_get, ibz_invmod, ibz_is_even, ibz_is_odd, ibz_is_one,
    ibz_is_zero, ibz_mod, ibz_mod_ui, ibz_mul, ibz_neg, ibz_pow, ibz_set, ibz_sub, ibz_two_adic,
    quat_alg_conj, quat_alg_elem_copy, quat_alg_mul, quat_alg_normalize, quat_lattice_alg_elem_mul,
    quat_lattice_contains, quat_lideal_conjugate_without_hnf, quat_lideal_copy,
    quat_lideal_lideal_mul_reduced, quat_lideal_reduce_basis, Ibz, QuatAlg, QuatAlgElem,
    QuatLattice, QuatLeftIdeal,
};

use crate::id2iso::endomorphism_application_even_basis;

/// `QUAT_repres_bound_input` (lvl1). Used by [`fixed_degree_isogeny_and_eval`]
/// to size the dimension-two isogeny step.
pub const QUAT_REPRES_BOUND_INPUT: i32 = 20;

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

/// `fixed_degree_isogeny_and_eval(lideal, u, small, E34, P12, numP,
/// index_alternate_order)`. Mirrors the C entry point.
///
/// **Body deferred**: requires `quat_represent_integer`, which is
/// RNG-driven and not yet ported in `sqisign-quaternion::normeq`.
#[allow(clippy::too_many_arguments)]
pub fn fixed_degree_isogeny_and_eval(
    _lideal: &mut QuatLeftIdeal,
    _u: &Ibz,
    _small: bool,
    _e34: &mut ThetaCoupleCurve,
    _p12: &mut [ThetaCouplePoint],
    _num_p: usize,
    _index_alternate_order: i32,
) -> i32 {
    panic!(
        "fixed_degree_isogeny_and_eval is RNG-driven (quat_represent_integer); \
         port pending sqisign-quaternion::normeq RNG path."
    );
}

/// `dim2id2iso_ideal_to_isogeny_clapotis(beta1, beta2, u, v, d1, d2,
/// codomain, basis, lideal, Bpoo)`. Mirrors the C entry point.
///
/// **Body deferred**: composes [`find_uv`] (deterministic) with two
/// `fixed_degree_isogeny_and_eval` calls and one
/// `theta_chain_compute_and_eval_randomized` step, all RNG-driven. See
/// crate-level docs.
#[allow(clippy::too_many_arguments)]
pub fn dim2id2iso_ideal_to_isogeny_clapotis(
    _beta1: &mut QuatAlgElem,
    _beta2: &mut QuatAlgElem,
    _u: &mut Ibz,
    _v: &mut Ibz,
    _d1: &mut Ibz,
    _d2: &mut Ibz,
    _codomain: &mut EcCurve,
    _basis: &mut EcBasis,
    _lideal: &QuatLeftIdeal,
    _bpoo: &QuatAlg,
) -> i32 {
    panic!(
        "dim2id2iso_ideal_to_isogeny_clapotis is RNG-driven (composes \
         fixed_degree_isogeny_and_eval); port pending sqisign-quaternion::normeq."
    );
}

/// `dim2id2iso_arbitrary_isogeny_evaluation(basis, codomain, lideal)`.
/// Wrapper around [`dim2id2iso_ideal_to_isogeny_clapotis`]. **Deferred**
/// for the same reason.
pub fn dim2id2iso_arbitrary_isogeny_evaluation(
    _basis: &mut EcBasis,
    _codomain: &mut EcCurve,
    _lideal: &QuatLeftIdeal,
) -> i32 {
    panic!(
        "dim2id2iso_arbitrary_isogeny_evaluation: RNG-driven dependency chain \
         (quat_represent_integer); port pending."
    );
}

#[allow(dead_code)]
fn _unused() {
    // Symbols imported solely to allow the deferred bodies to compile
    // when they are filled in. The arms are exercised by the RNG-driven
    // ports in a follow-up unit.
    let _ = (
        copy_basis,
        copy_curve,
        copy_point,
        ec_curve_normalize_a24,
        ec_dbl_iter_basis,
        ec_point_init,
        TORSION_EVEN_POWER,
        NWORDS_ORDER,
        copy_bases_to_kernel,
        double_couple_point_iter,
        theta_chain_compute_and_eval_randomized,
        ThetaCoupleCurve::zero,
        ThetaCouplePoint::zero,
        ThetaKernelCouplePoints::zero,
        HD_EXTRA_TORSION,
        CURVES_WITH_ENDOMORPHISMS.len(),
        NUM_ALTERNATE_EXTREMAL_ORDERS,
        QUATALG_PINFTY.p.clone(),
        TORSION_PLUS_2POWER.clone(),
        ibz_mat_2x2_new(),
        ibz_bitsize(&Ibz::zero()),
        ibz_const_one(),
        ibz_const_two(),
        ibz_const_zero(),
        ibz_get(&Ibz::zero()),
        ibz_is_one(&Ibz::zero()),
        ibz_is_odd(&Ibz::zero()),
        ibz_is_even(&Ibz::zero()),
        ibz_pow,
        ibz_two_adic(&Ibz::zero()),
        ibz_set,
        ibz_neg,
        endomorphism_application_even_basis,
        QUAT_REPRES_BOUND_INPUT,
        quat_alg_elem_copy,
    );
}
