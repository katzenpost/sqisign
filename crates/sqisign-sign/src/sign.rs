//! `sign.c` port: SQIsign Sign protocol.
//!
//! Mirrors `the-sqisign/src/signature/ref/lvlx/sign.c`. The
//! private helpers (`commit`, `compute_challenge_ideal_signature`, etc.)
//! are preserved as private functions in this module with the same
//! orchestration as the reference. RNG flows in through one
//! caller-supplied [`RngSource`] handle, threaded through every primitive
//! that consumes randomness.

use sqisign_common::RngSource;
use sqisign_ec::{
    copy_basis, copy_curve, copy_point, ec_dbl_iter, ec_dbl_iter_basis, ec_eval_even,
    ec_eval_small_chain, ec_iso_eval, ec_isomorphism, ec_ladder3pt, ec_mul, ec_normalize_curve,
    EcBasis, EcCurve, EcIsogEven, EcIsom, EcPoint, NWORDS_ORDER, TORSION_EVEN_POWER,
};
use sqisign_gf::fp2_copy;
use sqisign_hd::{
    copy_bases_to_kernel, double_couple_point_iter, theta_chain_compute_and_eval_randomized,
    ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints, HD_EXTRA_TORSION,
};
use sqisign_id2iso::{
    change_of_basis_matrix_tate, change_of_basis_matrix_tate_invert,
    dim2id2iso_arbitrary_isogeny_evaluation, ec_biscalar_mul_ibz_vec,
    id2iso_ideal_to_kernel_dlogs_even, id2iso_kernel_dlogs_to_ideal_even,
    matrix_application_even_basis, QUAT_EQUIV_BOUND_COEFF, QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
};
use sqisign_precomp::{
    COM_DEGREE, EXTREMAL_ORDERS, QUATALG_PINFTY, QUAT_PRIME_COFACTOR, TORSION_PLUS_2POWER,
};
use sqisign_quaternion::dim2::{ibz_vec_2_new, IbzVec2};
use sqisign_quaternion::{
    ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_copy_digits, ibz_div, ibz_invmod,
    ibz_is_one, ibz_mul, ibz_pow, ibz_set, ibz_sub, ibz_to_digits, ibz_two_adic, quat_alg_conj,
    quat_alg_make_primitive, quat_alg_norm, quat_lattice_conjugate_without_hnf,
    quat_lattice_intersect, quat_lattice_sample_from_ball, quat_lideal_create, quat_lideal_inter,
    quat_lideal_prime_norm_reduced_equivalent, quat_sampling_random_ideal_O0_given_norm, Ibz,
    QuatAlgElem, QuatLattice, QuatLeftIdeal, QuatRepresentIntegerParams,
};
use sqisign_verify::{
    hash_to_challenge, PublicKey, Signature, SECURITY_BITS, SQISIGN_RESPONSE_LENGTH,
};

use crate::keygen::SecretKey;

/// `theta_couple_curve_with_basis_t`. Mirrors the C struct of the same name.
#[derive(Clone, Debug)]
struct ThetaCoupleCurveWithBasis {
    pub e1: EcCurve,
    pub e2: EcCurve,
    pub b1: EcBasis,
    pub b2: EcBasis,
}

impl ThetaCoupleCurveWithBasis {
    fn zero() -> Self {
        Self {
            e1: EcCurve::zero(),
            e2: EcCurve::zero(),
            b1: EcBasis::zero(),
            b2: EcBasis::zero(),
        }
    }
}

/// `commit`. Sample a random ideal of norm COM_DEGREE, reduce to its
/// prime-norm equivalent, and run the Clapotis isogeny evaluator.
fn commit<R: RngSource>(
    rng: &mut R,
    e_com: &mut EcCurve,
    basis_even_com: &mut EcBasis,
    lideal_com: &mut QuatLeftIdeal,
) -> bool {
    let ri_params = QuatRepresentIntegerParams {
        primality_test_iterations: QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
        order: &EXTREMAL_ORDERS[0],
        algebra: &QUATALG_PINFTY,
    };
    let mut found =
        quat_sampling_random_ideal_O0_given_norm(rng, lideal_com, &COM_DEGREE, 1, &ri_params, None);
    if found != 0 {
        found = quat_lideal_prime_norm_reduced_equivalent(
            rng,
            lideal_com,
            &QUATALG_PINFTY,
            QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
            QUAT_EQUIV_BOUND_COEFF,
            &EXTREMAL_ORDERS[0].order,
        );
    }
    if found != 0 {
        found = dim2id2iso_arbitrary_isogeny_evaluation(rng, basis_even_com, e_com, lideal_com);
    }
    found != 0
}

/// `compute_challenge_ideal_signature`: build the lideal that the
/// challenge isogeny generates after pull-back through the secret-key
/// isogeny. Deterministic given the secret key and the recorded
/// challenge coefficient.
fn compute_challenge_ideal_signature(
    lideal_chall_two: &mut QuatLeftIdeal,
    sig: &Signature,
    sk: &SecretKey,
) {
    let mut vec = ibz_vec_2_new();
    ibz_set(&mut vec[0], 1);
    ibz_copy_digits(&mut vec[1], &sig.chall_coeff);

    let vec_clone = vec.clone();
    sqisign_quaternion::dim2::ibz_mat_2x2_eval(&mut vec, &sk.mat_BAcan_to_BA0_two, &vec_clone);

    id2iso_kernel_dlogs_to_ideal_even(
        lideal_chall_two,
        &vec,
        TORSION_EVEN_POWER as i32,
        &QUATALG_PINFTY,
    );
    debug_assert!(ibz_cmp(&lideal_chall_two.norm, &TORSION_PLUS_2POWER) == 0);
}

/// `sample_response`. Wrap [`quat_lattice_sample_from_ball`] with the
/// SQIsign-response bound.
fn sample_response<R: RngSource>(
    rng: &mut R,
    x: &mut QuatAlgElem,
    lattice: &QuatLattice,
    lattice_content: &Ibz,
) {
    let mut bound = Ibz::zero();
    ibz_pow(&mut bound, &ibz_const_two(), SQISIGN_RESPONSE_LENGTH as u32);
    let cur = bound.clone();
    ibz_sub(&mut bound, &cur, &ibz_const_one());
    let cur = bound.clone();
    ibz_mul(&mut bound, &cur, lattice_content);
    let ok = quat_lattice_sample_from_ball(rng, x, lattice, &QUATALG_PINFTY, &bound);
    debug_assert!(ok != 0, "sample_response: sampler returned 0");
    let _ = ok;
}

/// `compute_response_quat_element`. Construct the lattice
/// `dual(I_com) * (I_secret ∩ I_chall_two)` and sample a short element
/// from its ball; the resulting quaternion element is the response
/// generator.
fn compute_response_quat_element<R: RngSource>(
    rng: &mut R,
    resp_quat: &mut QuatAlgElem,
    lattice_content: &mut Ibz,
    sk: &SecretKey,
    lideal_chall_two: &QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
) {
    let mut lideal_chall_secret = QuatLeftIdeal::new();
    quat_lideal_inter(
        &mut lideal_chall_secret,
        lideal_chall_two,
        &sk.secret_ideal,
        &EXTREMAL_ORDERS[0].order,
    );

    let mut lat_commit = QuatLattice::new();
    quat_lattice_conjugate_without_hnf(&mut lat_commit, &lideal_commit.lattice);

    let mut lattice_hom_chall_to_com = QuatLattice::new();
    quat_lattice_intersect(
        &mut lattice_hom_chall_to_com,
        &lideal_chall_secret.lattice,
        &lat_commit,
    );

    ibz_mul(
        lattice_content,
        &lideal_chall_secret.norm,
        &lideal_commit.norm,
    );
    sample_response(rng, resp_quat, &lattice_hom_chall_to_com, lattice_content);
}

/// `compute_backtracking_signature`. Make `resp_quat` primitive and
/// record its 2-adic valuation as the signature's backtracking length.
fn compute_backtracking_signature(
    sig: &mut Signature,
    resp_quat: &mut QuatAlgElem,
    lattice_content: &mut Ibz,
    remain: &mut Ibz,
) {
    let mut tmp = Ibz::zero();
    let mut dummy_coord = sqisign_quaternion::dim4::ibz_vec_4_new();
    quat_alg_make_primitive(
        &mut dummy_coord,
        &mut tmp,
        resp_quat,
        &EXTREMAL_ORDERS[0].order,
    );
    let denom_clone = resp_quat.denom.clone();
    ibz_mul(&mut resp_quat.denom, &denom_clone, &tmp);

    let backtracking = ibz_two_adic(&tmp) as u8;
    sig.backtracking = backtracking;

    ibz_pow(&mut tmp, &ibz_const_two(), backtracking as u32);
    let cur = lattice_content.clone();
    ibz_div(lattice_content, remain, &cur, &tmp);
}

/// `compute_random_aux_norm_and_helpers`. Compute the response degree,
/// strip the power-of-two part, encode the response ideal, and produce
/// the random-aux-isogeny norm + degree-response inverse helpers.
#[allow(clippy::too_many_arguments)]
fn compute_random_aux_norm_and_helpers(
    sig: &mut Signature,
    random_aux_norm: &mut Ibz,
    degree_resp_inv: &mut Ibz,
    remain: &mut Ibz,
    lattice_content: &Ibz,
    resp_quat: &mut QuatAlgElem,
    lideal_com_resp: &mut QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
) -> u8 {
    let mut degree_full_resp = Ibz::zero();
    let mut degree_odd_resp = Ibz::zero();
    let mut norm_d = Ibz::zero();
    let mut tmp = Ibz::zero();

    quat_alg_norm(
        &mut degree_full_resp,
        &mut norm_d,
        resp_quat,
        &QUATALG_PINFTY,
    );
    debug_assert!(ibz_is_one(&norm_d) != 0);
    let cur = degree_full_resp.clone();
    ibz_div(&mut degree_full_resp, remain, &cur, lattice_content);
    debug_assert!(ibz_cmp(remain, &ibz_const_zero()) == 0);

    let exp_diadic = ibz_two_adic(&degree_full_resp) as u8;
    sig.two_resp_length = exp_diadic;

    ibz_pow(&mut tmp, &ibz_const_two(), exp_diadic as u32);
    ibz_div(&mut degree_odd_resp, remain, &degree_full_resp, &tmp);
    debug_assert!(ibz_cmp(remain, &ibz_const_zero()) == 0);

    // resp_quat = conj(resp_quat).
    let cur = resp_quat.clone();
    quat_alg_conj(resp_quat, &cur);

    ibz_mul(&mut tmp, &lideal_commit.norm, &degree_odd_resp);
    quat_lideal_create(
        lideal_com_resp,
        resp_quat,
        &tmp,
        &EXTREMAL_ORDERS[0].order,
        &QUATALG_PINFTY,
    );

    let pow_dim2_deg_resp =
        (SQISIGN_RESPONSE_LENGTH as i32) - (exp_diadic as i32) - (sig.backtracking as i32);
    debug_assert!(
        pow_dim2_deg_resp >= 0,
        "pow_dim2_deg_resp must be non-negative"
    );
    ibz_pow(remain, &ibz_const_two(), pow_dim2_deg_resp as u32);
    ibz_sub(random_aux_norm, remain, &degree_odd_resp);

    for _ in 0..HD_EXTRA_TORSION {
        let cur = remain.clone();
        ibz_mul(remain, &cur, &ibz_const_two());
    }

    let _ok = ibz_invmod(degree_resp_inv, &degree_odd_resp, remain);

    pow_dim2_deg_resp as u8
}

/// `evaluate_random_aux_isogeny_signature`. Sample a random ideal of a
/// given non-prime norm with the chosen cofactor, intersect with the
/// response ideal, and run the Clapotis evaluator.
fn evaluate_random_aux_isogeny_signature<R: RngSource>(
    rng: &mut R,
    e_aux: &mut EcCurve,
    b_aux: &mut EcBasis,
    norm: &Ibz,
    lideal_com_resp: &QuatLeftIdeal,
) -> bool {
    let mut lideal_aux = QuatLeftIdeal::new();
    let mut lideal_aux_resp_com = QuatLeftIdeal::new();

    let ri_params = QuatRepresentIntegerParams {
        primality_test_iterations: QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
        order: &EXTREMAL_ORDERS[0],
        algebra: &QUATALG_PINFTY,
    };
    let found = quat_sampling_random_ideal_O0_given_norm(
        rng,
        &mut lideal_aux,
        norm,
        0,
        &ri_params,
        Some(&QUAT_PRIME_COFACTOR),
    );
    if found == 0 {
        return false;
    }

    quat_lideal_inter(
        &mut lideal_aux_resp_com,
        lideal_com_resp,
        &lideal_aux,
        &EXTREMAL_ORDERS[0].order,
    );

    let evaluated =
        dim2id2iso_arbitrary_isogeny_evaluation(rng, b_aux, e_aux, &lideal_aux_resp_com);
    evaluated != 0
}

/// `compute_dim2_isogeny_challenge`. Build the dim-2 isogeny kernel from
/// the response-inverse-scaled commitment basis and run the randomised
/// theta chain.
fn compute_dim2_isogeny_challenge<R: RngSource>(
    rng: &mut R,
    codomain: &mut ThetaCoupleCurveWithBasis,
    domain: &ThetaCoupleCurveWithBasis,
    degree_resp_inv: &Ibz,
    pow_dim2_deg_resp: u8,
    exp_diadic_val_full_resp: u8,
    reduced_order: i32,
) -> bool {
    let mut ecom_x_eaux = ThetaCoupleCurve::zero();
    copy_curve(&mut ecom_x_eaux.e1, &domain.e1);
    copy_curve(&mut ecom_x_eaux.e2, &domain.e2);

    let mut dim_two_ker = ThetaKernelCouplePoints::zero();
    copy_bases_to_kernel(&mut dim_two_ker, &domain.b1, &domain.b2);

    let mut scalar = vec![0u64; NWORDS_ORDER];
    ibz_to_digits(&mut scalar, degree_resp_inv);

    let cur = dim_two_ker.t1.p2;
    ec_mul(
        &mut dim_two_ker.t1.p2,
        &scalar,
        reduced_order,
        &cur,
        &mut ecom_x_eaux.e2,
    );
    let cur = dim_two_ker.t2.p2;
    ec_mul(
        &mut dim_two_ker.t2.p2,
        &scalar,
        reduced_order,
        &cur,
        &mut ecom_x_eaux.e2,
    );
    let cur = dim_two_ker.t1m2.p2;
    ec_mul(
        &mut dim_two_ker.t1m2.p2,
        &scalar,
        reduced_order,
        &cur,
        &mut ecom_x_eaux.e2,
    );

    let cur = dim_two_ker.t1;
    double_couple_point_iter(
        &mut dim_two_ker.t1,
        exp_diadic_val_full_resp as u32,
        &cur,
        &ecom_x_eaux,
    );
    let cur = dim_two_ker.t2;
    double_couple_point_iter(
        &mut dim_two_ker.t2,
        exp_diadic_val_full_resp as u32,
        &cur,
        &ecom_x_eaux,
    );
    let cur = dim_two_ker.t1m2;
    double_couple_point_iter(
        &mut dim_two_ker.t1m2,
        exp_diadic_val_full_resp as u32,
        &cur,
        &ecom_x_eaux,
    );

    let mut pushed_points: [ThetaCouplePoint; 3] = [
        ThetaCouplePoint::zero(),
        ThetaCouplePoint::zero(),
        ThetaCouplePoint::zero(),
    ];
    copy_point(&mut pushed_points[0].p1, &domain.b1.P);
    copy_point(&mut pushed_points[1].p1, &domain.b1.Q);
    copy_point(&mut pushed_points[2].p1, &domain.b1.PmQ);
    sqisign_ec::ec_point_init(&mut pushed_points[0].p2);
    sqisign_ec::ec_point_init(&mut pushed_points[1].p2);
    sqisign_ec::ec_point_init(&mut pushed_points[2].p2);

    let mut codomain_product = ThetaCoupleCurve::zero();
    let ret = theta_chain_compute_and_eval_randomized(
        rng,
        pow_dim2_deg_resp as u32,
        &mut ecom_x_eaux,
        &dim_two_ker,
        true,
        &mut codomain_product,
        &mut pushed_points,
    );
    if ret == 0 {
        return false;
    }
    let _ = reduced_order;

    copy_curve(&mut codomain.e1, &codomain_product.e2);
    copy_curve(&mut codomain.e2, &codomain_product.e1);

    copy_point(&mut codomain.b1.P, &pushed_points[0].p2);
    copy_point(&mut codomain.b1.Q, &pushed_points[1].p2);
    copy_point(&mut codomain.b1.PmQ, &pushed_points[2].p2);

    copy_point(&mut codomain.b2.P, &pushed_points[0].p1);
    copy_point(&mut codomain.b2.Q, &pushed_points[1].p1);
    copy_point(&mut codomain.b2.PmQ, &pushed_points[2].p1);
    true
}

/// `compute_small_chain_isogeny_signature`. The 2-power-chain step that
/// finishes the response when `two_resp_length > 0`.
fn compute_small_chain_isogeny_signature(
    e_chall_2: &mut EcCurve,
    b_chall_2: &mut EcBasis,
    resp_quat: &QuatAlgElem,
    pow_dim2_deg_resp: u8,
    length: i32,
) -> bool {
    let mut two_pow = Ibz::zero();
    let mut vec_resp_two: IbzVec2 = ibz_vec_2_new();
    let mut lideal_resp_two = QuatLeftIdeal::new();

    ibz_pow(&mut two_pow, &ibz_const_two(), length as u32);
    quat_lideal_create(
        &mut lideal_resp_two,
        resp_quat,
        &two_pow,
        &EXTREMAL_ORDERS[0].order,
        &QUATALG_PINFTY,
    );

    id2iso_ideal_to_kernel_dlogs_even(&mut vec_resp_two, &lideal_resp_two, &QUATALG_PINFTY);

    let mut points: [EcPoint; 3] = [EcPoint::zero(); 3];
    copy_point(&mut points[0], &b_chall_2.P);
    copy_point(&mut points[1], &b_chall_2.Q);
    copy_point(&mut points[2], &b_chall_2.PmQ);

    let drop = (pow_dim2_deg_resp as i32) + (HD_EXTRA_TORSION as i32);
    let bas_clone = *b_chall_2;
    ec_dbl_iter_basis(b_chall_2, drop, &bas_clone, e_chall_2);

    let mut ker = EcPoint::zero();
    ec_biscalar_mul_ibz_vec(&mut ker, &vec_resp_two, length, b_chall_2, e_chall_2);

    if ec_eval_small_chain(e_chall_2, &ker, length, &mut points, true) != 0 {
        return false;
    }

    copy_point(&mut b_chall_2.P, &points[0]);
    copy_point(&mut b_chall_2.Q, &points[1]);
    copy_point(&mut b_chall_2.PmQ, &points[2]);
    true
}

/// `compute_challenge_codomain_signature`. Derive the challenge codomain
/// curve by walking the challenge isogeny from the secret-key curve;
/// then push the dim-2 basis through the isomorphism connecting the two
/// representations.
fn compute_challenge_codomain_signature(
    sig: &Signature,
    sk: &SecretKey,
    e_chall: &mut EcCurve,
    _e_chall_2: &EcCurve,
    b_chall_2: &mut EcBasis,
) -> bool {
    let mut phi_chall = EcIsogEven::zero();
    let mut bas_sk = EcBasis::zero();
    copy_basis(&mut bas_sk, &sk.canonical_basis);

    phi_chall.curve = sk.curve;
    phi_chall.length = (TORSION_EVEN_POWER as i32 - sig.backtracking as i32) as u32;

    let mut sk_curve = sk.curve;
    ec_ladder3pt(
        &mut phi_chall.kernel,
        &sig.chall_coeff,
        &bas_sk.P,
        &bas_sk.Q,
        &bas_sk.PmQ,
        &sk_curve,
    );

    let mut kernel = phi_chall.kernel;
    ec_dbl_iter(
        &mut phi_chall.kernel,
        sig.backtracking as i32,
        &kernel,
        &mut sk_curve,
    );
    let _ = kernel; // (suppress unused-var in release builds)
    kernel = phi_chall.kernel;
    let _ = kernel;

    if ec_eval_even(e_chall, &phi_chall, &mut []) != 0 {
        return false;
    }

    let mut isom = EcIsom::zero();
    if ec_isomorphism(&mut isom, _e_chall_2, e_chall) != 0 {
        return false;
    }
    ec_iso_eval(&mut b_chall_2.P, &isom);
    ec_iso_eval(&mut b_chall_2.Q, &isom);
    ec_iso_eval(&mut b_chall_2.PmQ, &isom);
    true
}

/// `set_aux_curve_signature`. Normalize the auxiliary curve and copy
/// its Montgomery A-coefficient into the signature.
fn set_aux_curve_signature(sig: &mut Signature, e_aux: &mut EcCurve) {
    ec_normalize_curve(e_aux);
    fp2_copy(&mut sig.E_aux_A, &e_aux.A);
}

/// `compute_and_set_basis_change_matrix`. Compute the two canonical
/// bases (with hints recorded into the signature), then the change-of-
/// basis matrix from `B_chall_can` to `B_chall_2`. Pack each cell as a
/// `NWORDS_ORDER`-long digit array.
fn compute_and_set_basis_change_matrix(
    sig: &mut Signature,
    b_aux_2: &EcBasis,
    b_chall_2: &mut EcBasis,
    e_aux_2: &mut EcCurve,
    e_chall: &mut EcCurve,
    f: i32,
) {
    let mut mat_baux2_to_can = sqisign_quaternion::dim2::ibz_mat_2x2_new();
    let mut mat_can_to_bchall = sqisign_quaternion::dim2::ibz_mat_2x2_new();

    let mut b_can_chall = EcBasis::zero();
    let mut b_aux_2_can = EcBasis::zero();
    sig.hint_chall = sqisign_ec::ec_curve_to_basis_2f_to_hint(
        &mut b_can_chall,
        e_chall,
        TORSION_EVEN_POWER as i32,
    );
    sig.hint_aux = sqisign_ec::ec_curve_to_basis_2f_to_hint(
        &mut b_aux_2_can,
        e_aux_2,
        TORSION_EVEN_POWER as i32,
    );

    change_of_basis_matrix_tate_invert(&mut mat_baux2_to_can, &b_aux_2_can, b_aux_2, e_aux_2, f);
    matrix_application_even_basis(b_chall_2, e_chall, &mut mat_baux2_to_can, f);

    change_of_basis_matrix_tate(&mut mat_can_to_bchall, b_chall_2, &b_can_chall, e_chall, f);

    for i in 0..2 {
        for j in 0..2 {
            let mut digits = [0u64; NWORDS_ORDER];
            ibz_to_digits(&mut digits, &mat_can_to_bchall[i][j]);
            sig.mat_Bchall_can_to_B_chall[i][j] = digits;
        }
    }
}

/// `protocols_sign(rng, sig, pk, sk, message)`. The top-level Sign
/// orchestration.
pub fn protocols_sign<R: RngSource>(
    rng: &mut R,
    sig: &mut Signature,
    pk: &PublicKey,
    sk: &mut SecretKey,
    message: &[u8],
) -> i32 {
    let mut ret = false;
    let mut reduced_order: i32 = 0;

    let mut remain = Ibz::zero();
    let mut lattice_content = Ibz::zero();
    let mut random_aux_norm = Ibz::zero();
    let mut degree_resp_inv = Ibz::zero();

    let mut resp_quat = QuatAlgElem::new();

    let mut lideal_commit = QuatLeftIdeal::new();
    let mut lideal_com_resp = QuatLeftIdeal::new();

    let mut ecom_eaux = ThetaCoupleCurveWithBasis::zero();
    let mut eaux2_echall2 = ThetaCoupleCurveWithBasis::zero();
    let mut e_chall = sk.curve;

    let mut pow_dim2_deg_resp: u8;

    while !ret {
        ret = commit(
            rng,
            &mut ecom_eaux.e1,
            &mut ecom_eaux.b1,
            &mut lideal_commit,
        );
        if !ret {
            continue;
        }

        // Hash the message into the challenge scalar.
        let mut chall = [0u64; NWORDS_ORDER];
        hash_to_challenge(&mut chall, pk, &ecom_eaux.e1, message);
        sig.chall_coeff = chall;

        {
            let mut lideal_chall_two = QuatLeftIdeal::new();
            compute_challenge_ideal_signature(&mut lideal_chall_two, sig, sk);
            compute_response_quat_element(
                rng,
                &mut resp_quat,
                &mut lattice_content,
                sk,
                &lideal_chall_two,
                &lideal_commit,
            );
        }

        compute_backtracking_signature(sig, &mut resp_quat, &mut lattice_content, &mut remain);

        pow_dim2_deg_resp = compute_random_aux_norm_and_helpers(
            sig,
            &mut random_aux_norm,
            &mut degree_resp_inv,
            &mut remain,
            &lattice_content,
            &mut resp_quat,
            &mut lideal_com_resp,
            &lideal_commit,
        );

        if pow_dim2_deg_resp > 0 {
            ret = evaluate_random_aux_isogeny_signature(
                rng,
                &mut ecom_eaux.e2,
                &mut ecom_eaux.b2,
                &random_aux_norm,
                &lideal_com_resp,
            );
            if !ret {
                continue;
            }

            reduced_order = (pow_dim2_deg_resp as i32)
                + (HD_EXTRA_TORSION as i32)
                + (sig.two_resp_length as i32);
            let drop = (TORSION_EVEN_POWER as i32) - reduced_order;
            let b1_clone = ecom_eaux.b1;
            ec_dbl_iter_basis(&mut ecom_eaux.b1, drop, &b1_clone, &mut ecom_eaux.e1);
            let b2_clone = ecom_eaux.b2;
            ec_dbl_iter_basis(&mut ecom_eaux.b2, drop, &b2_clone, &mut ecom_eaux.e2);

            ret = compute_dim2_isogeny_challenge(
                rng,
                &mut eaux2_echall2,
                &ecom_eaux,
                &degree_resp_inv,
                pow_dim2_deg_resp,
                sig.two_resp_length,
                reduced_order,
            );
            if !ret {
                continue;
            }
        } else {
            copy_curve(&mut eaux2_echall2.e1, &ecom_eaux.e1);
            copy_curve(&mut eaux2_echall2.e2, &ecom_eaux.e1);
            reduced_order = sig.two_resp_length as i32;
            let drop = (TORSION_EVEN_POWER as i32) - reduced_order;
            let b1_clone = ecom_eaux.b1;
            ec_dbl_iter_basis(&mut eaux2_echall2.b1, drop, &b1_clone, &mut ecom_eaux.e1);
            // (The C reference applies the same double-iter twice; we
            // mirror that verbatim to keep the exit basis identical.)
            let b1_clone = ecom_eaux.b1;
            ec_dbl_iter_basis(&mut eaux2_echall2.b1, drop, &b1_clone, &mut ecom_eaux.e1);
            copy_basis(&mut eaux2_echall2.b2, &eaux2_echall2.b1);
        }

        if sig.two_resp_length > 0 {
            let ok = compute_small_chain_isogeny_signature(
                &mut eaux2_echall2.e2,
                &mut eaux2_echall2.b2,
                &resp_quat,
                pow_dim2_deg_resp,
                sig.two_resp_length as i32,
            );
            debug_assert!(ok, "compute_small_chain_isogeny_signature failed");
            let _ = ok;
        }

        let _ok2 = compute_challenge_codomain_signature(
            sig,
            sk,
            &mut e_chall,
            &eaux2_echall2.e2,
            &mut eaux2_echall2.b2,
        );
        debug_assert!(_ok2, "compute_challenge_codomain_signature failed");
    }

    set_aux_curve_signature(sig, &mut eaux2_echall2.e1);
    compute_and_set_basis_change_matrix(
        sig,
        &eaux2_echall2.b1,
        &mut eaux2_echall2.b2,
        &mut eaux2_echall2.e1,
        &mut e_chall,
        reduced_order,
    );

    let _ = SECURITY_BITS;
    if ret {
        1
    } else {
        0
    }
}
