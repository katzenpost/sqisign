//! Norm-equation helpers.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/normeq.c`. The
//! deterministic helpers (`quat_lattice_O0_set`, `quat_lattice_O0_set_extremal`,
//! `quat_order_elem_create`, `quat_change_to_O0_basis`) are joined here by the
//! two RNG-driven entry points the signing path needs: [`quat_represent_integer`]
//! and [`quat_sampling_random_ideal_O0_given_norm`]. Both pull bytes through
//! `&mut impl RngSource` (replacing the reference's thread-local DRBG), so the
//! same fixed-seed entropy that drives the C harness drives the Rust port; the
//! differential vectors prove the two agree value-for-value.
//!
//! ## Upstream pinned quirks
//!
//! * `quat_represent_integer` reuses `temp` as both the upper-bound and the
//!   in-place output of `ibz_sqrt_floor`, then immediately divides
//!   `adjusted_n_gamma` by that floored square root with a final aliased
//!   division (`ibz_div(&counter, &temp, &adjusted_n_gamma, &temp)`). We mirror
//!   the aliased divisions verbatim via explicit clones, since `num-bigint`
//!   forbids the alias. See `normeq.c:138` and `normeq.c:154`.
//!
//! * `quat_sampling_random_ideal_O0_given_norm` performs a `found = found && !zero`
//!   short-circuit on the composite path then unconditionally clobbers `found`
//!   back to zero before entering the rerandomization loop. The first
//!   assignment is therefore a dead store. We keep it for side-by-side review
//!   parity. See `normeq.c:306-319`.
//!
//! * `quat_represent_integer` uses the C `int` arithmetic `(ibz_get(c0) -
//!   ibz_get(c3)) % 4 == 2` for the non-diagonal parity check. We mirror the
//!   wrapping-subtraction and truncated-remainder semantics with
//!   `i32::wrapping_sub` and Rust's `%`. The values in question are positive
//!   and well below `i32::MAX` for every input we test, so the wrap is
//!   inert in practice; the explicit `wrapping_sub` is for faithful porting
//!   rather than load-bearing correctness. See `normeq.c:187-188`.

#![allow(non_snake_case)]

use crate::algebra::{
    quat_alg_add, quat_alg_elem_is_zero, quat_alg_make_primitive, quat_alg_mul, quat_alg_norm,
    quat_alg_scalar, QuatAlg, QuatAlgElem,
};
use crate::dim4::{ibz_mat_4x4_eval, ibz_vec_4_new, IbzVec4};
use crate::ibz::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_div, ibz_divides, ibz_gcd,
    ibz_get, ibz_is_even, ibz_is_odd, ibz_is_one, ibz_is_zero, ibz_mod, ibz_mul, ibz_neg,
    ibz_probab_prime, ibz_set, ibz_sqrt_floor, ibz_sqrt_mod_p, ibz_sub, Ibz,
};
use crate::ibz_rand::ibz_rand_interval;
use crate::ideal::{quat_lideal_create, QuatLeftIdeal};
use crate::integers::ibz_cornacchia_prime;
use crate::lattice::QuatLattice;

use sqisign_common::RngSource;

/// `quat_p_extremal_maximal_order_t`.
#[derive(Clone, Debug, Default)]
pub struct QuatPExtremalMaximalOrder {
    pub order: QuatLattice,
    pub z: QuatAlgElem,
    pub t: QuatAlgElem,
    pub q: u32,
}

impl QuatPExtremalMaximalOrder {
    pub fn new() -> Self {
        Self {
            order: QuatLattice::new(),
            z: QuatAlgElem::new(),
            t: QuatAlgElem::new(),
            q: 0,
        }
    }
}

/// `quat_lattice_O0_set(O0)`: the order `(1, i, (i+j)/2, (1+ij)/2)`.
pub fn quat_lattice_O0_set(o0: &mut QuatLattice) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut o0.basis[i][j], 0);
        }
    }
    ibz_set(&mut o0.denom, 2);
    ibz_set(&mut o0.basis[0][0], 2);
    ibz_set(&mut o0.basis[1][1], 2);
    ibz_set(&mut o0.basis[2][2], 1);
    ibz_set(&mut o0.basis[1][2], 1);
    ibz_set(&mut o0.basis[3][3], 1);
    ibz_set(&mut o0.basis[0][3], 1);
}

/// `quat_lattice_O0_set_extremal(O0)`.
pub fn quat_lattice_O0_set_extremal(o0: &mut QuatPExtremalMaximalOrder) {
    ibz_set(&mut o0.z.coord[1], 1);
    ibz_set(&mut o0.t.coord[2], 1);
    ibz_set(&mut o0.z.denom, 1);
    ibz_set(&mut o0.t.denom, 1);
    o0.q = 1;
    quat_lattice_O0_set(&mut o0.order);
}

/// `quat_order_elem_create(elem, order, coeffs, alg)`: build an algebra
/// element from its `(1, z, t, t*z)`-basis coefficients.
pub fn quat_order_elem_create(
    elem: &mut QuatAlgElem,
    order: &QuatPExtremalMaximalOrder,
    coeffs: &IbzVec4,
    bpoo: &QuatAlg,
) {
    let mut quat_temp = QuatAlgElem::new();

    quat_alg_scalar(elem, &coeffs[0], &ibz_const_one());

    quat_alg_scalar(&mut quat_temp, &coeffs[1], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.z, &tmp_clone, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);

    quat_alg_scalar(&mut quat_temp, &coeffs[2], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &tmp_clone, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);

    quat_alg_scalar(&mut quat_temp, &coeffs[3], &ibz_const_one());
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &tmp_clone, bpoo);
    let tmp_clone = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &tmp_clone, &order.z, bpoo);

    let elem_clone = elem.clone();
    quat_alg_add(elem, &elem_clone, &quat_temp);
}

/// `quat_represent_integer_params_t`: parameters bundle for
/// [`quat_represent_integer`] and [`quat_sampling_random_ideal_O0_given_norm`].
///
/// Mirrors the C struct field-for-field; `order` and `algebra` are borrowed so
/// callers may pass a single immutable copy of each by reference.
#[derive(Clone, Debug)]
pub struct QuatRepresentIntegerParams<'a> {
    pub primality_test_iterations: i32,
    pub order: &'a QuatPExtremalMaximalOrder,
    pub algebra: &'a QuatAlg,
}

/// `quat_represent_integer(gamma, n_gamma, non_diag, params)`.
///
/// Finds a quaternion element `gamma` of (adjusted) norm `n_gamma` in
/// `params.order`. The non_diag flag enables the off-diagonal parity check used
/// when the resulting endomorphism must behave well for dim-2 computations.
/// Returns 1 on success; 0 if the bounded search exhausts its counter.
///
/// The reference draws two coordinates from the DRBG via `ibz_rand_interval`,
/// solves Cornacchia for the other two, then post-filters for primitivity and
/// (optionally) parity. We mirror the control flow verbatim so the byte stream
/// the underlying DRBG sees matches what C consumes.
pub fn quat_represent_integer<R: RngSource>(
    rng: &mut R,
    gamma: &mut QuatAlgElem,
    n_gamma: &Ibz,
    non_diag: i32,
    params: &QuatRepresentIntegerParams<'_>,
) -> i32 {
    if ibz_is_even(n_gamma) != 0 {
        return 0;
    }

    let mut found = 0i32;
    let mut bound = Ibz::zero();
    let mut temp = Ibz::zero();
    let mut q = Ibz::zero();
    let mut sq_bound = Ibz::zero();
    let mut coeffs: IbzVec4 = ibz_vec_4_new();
    let mut adjusted_n_gamma = Ibz::zero();
    let mut cornacchia_target = Ibz::zero();

    if non_diag != 0 {
        assert_eq!(params.order.q % 4, 1);
    }

    ibz_set(&mut q, params.order.q as i32);

    let standard_order = params.order.q == 1;

    // Adjusting the norm of gamma (multiplying by 4 to find a solution in
    // an order of odd level).
    if non_diag != 0 || standard_order {
        ibz_mul(&mut adjusted_n_gamma, n_gamma, &ibz_const_two());
        let cur = adjusted_n_gamma.clone();
        ibz_mul(&mut adjusted_n_gamma, &cur, &ibz_const_two());
    } else {
        adjusted_n_gamma = n_gamma.clone();
    }

    // First bound = sqrt(adjusted_n_gamma / p - q).
    let mut rmd = Ibz::zero();
    ibz_div(&mut sq_bound, &mut rmd, &adjusted_n_gamma, &params.algebra.p);
    ibz_set(&mut temp, params.order.q as i32);
    let cur = sq_bound.clone();
    ibz_sub(&mut sq_bound, &cur, &temp);
    ibz_sqrt_floor(&mut bound, &sq_bound);

    // counter = adjusted_n_gamma / sqrt(q * p * p).
    let mut counter = Ibz::zero();
    let cur = temp.clone();
    ibz_mul(&mut temp, &cur, &params.algebra.p);
    let cur = temp.clone();
    ibz_mul(&mut temp, &cur, &params.algebra.p);
    let cur = temp.clone();
    ibz_sqrt_floor(&mut temp, &cur);
    // C: ibz_div(&counter, &temp, &adjusted_n_gamma, &temp). The divisor and
    // remainder output alias the same `temp` slot; clone the divisor first.
    let divisor = temp.clone();
    ibz_div(&mut counter, &mut temp, &adjusted_n_gamma, &divisor);

    while found == 0 && ibz_cmp(&counter, &ibz_const_zero()) != 0 {
        let cur = counter.clone();
        ibz_sub(&mut counter, &cur, &ibz_const_one());

        // First coordinate: c[2] in [1, bound].
        ibz_rand_interval(rng, &mut coeffs[2], &ibz_const_one(), &bound);

        // Second-bound = sqrt((adjusted_n_gamma - p*c[2]^2) / (q*p)).
        ibz_mul(&mut cornacchia_target, &coeffs[2], &coeffs[2]);
        ibz_mul(&mut temp, &cornacchia_target, &params.algebra.p);
        let cur = temp.clone();
        ibz_sub(&mut temp, &adjusted_n_gamma, &cur);
        ibz_mul(&mut sq_bound, &q, &params.algebra.p);
        // C: ibz_div(&temp, &sq_bound, &temp, &sq_bound). Both temp and
        // sq_bound alias the divisor and output rem slots; clone both.
        let num = temp.clone();
        let den = sq_bound.clone();
        ibz_div(&mut temp, &mut sq_bound, &num, &den);
        let cur = temp.clone();
        ibz_sqrt_floor(&mut temp, &cur);

        if ibz_cmp(&temp, &ibz_const_zero()) == 0 {
            continue;
        }

        // Second coordinate: c[3] in [1, temp].
        ibz_rand_interval(rng, &mut coeffs[3], &ibz_const_one(), &temp);

        // cornacchia_target = adjusted_n_gamma - p * (c[2]^2 + q*c[3]^2).
        ibz_mul(&mut temp, &coeffs[3], &coeffs[3]);
        let cur = temp.clone();
        ibz_mul(&mut temp, &q, &cur);
        let cur = cornacchia_target.clone();
        ibz_add(&mut cornacchia_target, &cur, &temp);
        let cur = cornacchia_target.clone();
        ibz_mul(&mut cornacchia_target, &cur, &params.algebra.p);
        let cur = cornacchia_target.clone();
        ibz_sub(&mut cornacchia_target, &adjusted_n_gamma, &cur);
        debug_assert!(ibz_cmp(&cornacchia_target, &ibz_const_zero()) > 0);

        // Cornacchia: solve c[0]^2 + q * c[1]^2 = cornacchia_target.
        if ibz_probab_prime(&cornacchia_target, params.primality_test_iterations) != 0 {
            let mut c0 = Ibz::zero();
            let mut c1 = Ibz::zero();
            found = ibz_cornacchia_prime(&mut c0, &mut c1, &q, &cornacchia_target);
            coeffs[0] = c0;
            coeffs[1] = c1;
        } else {
            found = 0;
        }

        if found != 0 && non_diag != 0 && standard_order {
            // Parity adjustment to ensure x = t mod 2 and y = z mod 2: with
            // q=1 we may swap c[0] and c[1] freely.
            if ibz_is_odd(&coeffs[0]) != ibz_is_odd(&coeffs[3]) {
                let tmp = coeffs[0].clone();
                coeffs[0] = coeffs[1].clone();
                coeffs[1] = tmp;
            }
            // Further require (x-t)/2 odd and (y-z)/2 odd so the
            // resulting endomorphism is well-behaved for dim-2 work.
            let d03 = ibz_get(&coeffs[0]).wrapping_sub(ibz_get(&coeffs[3]));
            let d12 = ibz_get(&coeffs[1]).wrapping_sub(ibz_get(&coeffs[2]));
            let cond = (d03 % 4 == 2) && (d12 % 4 == 2);
            if !cond {
                found = 0;
            }
        }

        if found != 0 {
            // Translate (x, y, z, t) into the algebra element gamma.
            quat_order_elem_create(gamma, params.order, &coeffs, params.algebra);

            // Primitivize: coeffs is overwritten with the primitive
            // coordinates in the order's basis; temp receives the content.
            quat_alg_make_primitive(&mut coeffs, &mut temp, gamma, &params.order.order);

            if non_diag != 0 || standard_order {
                found = (ibz_cmp(&temp, &ibz_const_two()) == 0) as i32;
            } else {
                found = (ibz_cmp(&temp, &ibz_const_one()) == 0) as i32;
            }
        }
    }

    if found != 0 {
        // Recompose: gamma.coord = order.basis * coeffs, gamma.denom = order.denom.
        let saved = coeffs.clone();
        ibz_mat_4x4_eval(&mut coeffs, &params.order.order.basis, &saved);
        gamma.coord[0] = coeffs[0].clone();
        gamma.coord[1] = coeffs[1].clone();
        gamma.coord[2] = coeffs[2].clone();
        gamma.coord[3] = coeffs[3].clone();
        gamma.denom = params.order.order.denom.clone();
    }

    found
}

/// `quat_sampling_random_ideal_O0_given_norm(lideal, norm, is_prime, params, prime_cofactor)`.
///
/// Produces a uniformly random left `O0`-ideal of the requested norm. The fast
/// prime path samples a trace-zero generator and recovers the first coordinate
/// via a modular square root; the composite path delegates to
/// [`quat_represent_integer`] with `prime_cofactor * norm` as the norm target.
/// In both cases the ideal class is then rerandomized by left-multiplication
/// with a uniform element coprime to `norm`.
///
/// `prime_cofactor` must be `Some` on the composite path; the C reference
/// asserts `prime_cofactor != NULL` there.
pub fn quat_sampling_random_ideal_O0_given_norm<R: RngSource>(
    rng: &mut R,
    lideal: &mut QuatLeftIdeal,
    norm: &Ibz,
    is_prime: i32,
    params: &QuatRepresentIntegerParams<'_>,
    prime_cofactor: Option<&Ibz>,
) -> i32 {
    let mut n_temp = Ibz::zero();
    let mut norm_d = Ibz::zero();
    let mut disc = Ibz::zero();
    let mut gen = QuatAlgElem::new();
    let mut gen_rerand = QuatAlgElem::new();
    let mut found = 0i32;

    if is_prime != 0 {
        while found == 0 {
            // Trace-zero element with coordinates in [0, norm - 1].
            ibz_set(&mut gen.coord[0], 0);
            ibz_set(&mut gen.denom, 1);
            ibz_sub(&mut n_temp, norm, &ibz_const_one());
            for i in 1..4 {
                ibz_rand_interval(rng, &mut gen.coord[i], &ibz_const_zero(), &n_temp);
            }

            quat_alg_norm(&mut n_temp, &mut norm_d, &gen, params.algebra);
            debug_assert!(ibz_is_one(&norm_d) != 0);

            ibz_neg(&mut disc, &n_temp);
            let cur = disc.clone();
            ibz_mod(&mut disc, &cur, norm);

            found = ibz_sqrt_mod_p(&mut gen.coord[0], &disc, norm);
            if found != 0 && quat_alg_elem_is_zero(&gen) != 0 {
                found = 0;
            }
        }
    } else {
        let cof = prime_cofactor.expect("prime_cofactor must be set for composite norm");
        assert!(ibz_is_zero(norm) == 0);
        ibz_mul(&mut n_temp, cof, norm);
        found = quat_represent_integer(rng, &mut gen, &n_temp, 0, params);
        if found != 0 && quat_alg_elem_is_zero(&gen) != 0 {
            found = 0;
        }
        let _ = found; // dead store: the next line clobbers found to 0.
    }

    // Rerandomize the ideal class with a uniform element coprime to `norm`.
    found = 0;
    while found == 0 {
        for i in 0..4 {
            ibz_rand_interval(rng, &mut gen_rerand.coord[i], &ibz_const_one(), norm);
        }
        ibz_set(&mut gen_rerand.denom, 1);
        quat_alg_norm(&mut n_temp, &mut norm_d, &gen_rerand, params.algebra);
        debug_assert!(ibz_is_one(&norm_d) != 0);
        ibz_gcd(&mut disc, &n_temp, norm);
        found = ibz_is_one(&disc);
        if found != 0 && quat_alg_elem_is_zero(&gen_rerand) != 0 {
            found = 0;
        }
    }

    let gen_saved = gen.clone();
    quat_alg_mul(&mut gen, &gen_saved, &gen_rerand, params.algebra);
    quat_lideal_create(lideal, &gen, norm, &params.order.order, params.algebra);
    debug_assert!(ibz_cmp(norm, &lideal.norm) == 0);

    found
}

/// `quat_change_to_O0_basis(vec, el)`.
pub fn quat_change_to_O0_basis(vec: &mut IbzVec4, el: &QuatAlgElem) {
    vec[2] = el.coord[2].clone();
    let v2_clone = vec[2].clone();
    ibz_add(&mut vec[2], &v2_clone, &v2_clone);
    vec[3] = el.coord[3].clone();
    let v3_clone = vec[3].clone();
    ibz_add(&mut vec[3], &v3_clone, &v3_clone);
    ibz_sub(&mut vec[0], &el.coord[0], &el.coord[3]);
    ibz_sub(&mut vec[1], &el.coord[1], &el.coord[2]);

    assert!(ibz_divides(&vec[0], &el.denom) != 0);
    assert!(ibz_divides(&vec[1], &el.denom) != 0);
    assert!(ibz_divides(&vec[2], &el.denom) != 0);
    assert!(ibz_divides(&vec[3], &el.denom) != 0);

    let mut tmp = Ibz::zero();
    for i in 0..4 {
        let saved = vec[i].clone();
        ibz_div(&mut vec[i], &mut tmp, &saved, &el.denom);
    }
}
