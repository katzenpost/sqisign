//! `keygen.c` port: sample a random secret ideal, reduce to prime-norm
//! equivalent, build the secret isogeny, derive the public key.

use sqisign_common::RngSource;
use sqisign_ec::{ec_curve_to_basis_2f_to_hint, EcBasis, EcCurve, TORSION_EVEN_POWER};
use sqisign_gf::fp2_is_one;
use sqisign_id2iso::{
    change_of_basis_matrix_tate, dim2id2iso_arbitrary_isogeny_evaluation, QUAT_EQUIV_BOUND_COEFF,
    QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
};
use sqisign_precomp::{EXTREMAL_ORDERS, QUATALG_PINFTY, SEC_DEGREE};
use sqisign_quaternion::dim2::{ibz_mat_2x2_new, IbzMat2x2};
use sqisign_quaternion::{
    quat_lideal_prime_norm_reduced_equivalent, quat_sampling_random_ideal_O0_given_norm,
    QuatLeftIdeal, QuatRepresentIntegerParams,
};
use sqisign_verify::PublicKey;

/// `secret_key_t`: the secret key bundle. Mirrors the C struct of the
/// same name in `the-sqisign/src/signature/ref/include/signature.h`.
#[derive(Clone, Debug)]
pub struct SecretKey {
    pub secret_ideal: QuatLeftIdeal,
    pub mat_BAcan_to_BA0_two: IbzMat2x2,
    pub curve: EcCurve,
    pub canonical_basis: EcBasis,
}

impl SecretKey {
    pub fn new() -> Self {
        Self {
            secret_ideal: QuatLeftIdeal::new(),
            mat_BAcan_to_BA0_two: ibz_mat_2x2_new(),
            curve: EcCurve::zero(),
            canonical_basis: EcBasis::zero(),
        }
    }
}

impl Default for SecretKey {
    fn default() -> Self {
        Self::new()
    }
}

/// `protocols_keygen(pk, sk)`. Mirrors the C entry point. Returns 1 on
/// success; the C reference's `int found` flag is exposed verbatim.
pub fn protocols_keygen<R: RngSource>(rng: &mut R, pk: &mut PublicKey, sk: &mut SecretKey) -> i32 {
    let mut found = 0i32;
    let mut b_0_two = EcBasis::zero();

    while found == 0 {
        let ri_params = QuatRepresentIntegerParams {
            primality_test_iterations: QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
            order: &EXTREMAL_ORDERS[0],
            algebra: &QUATALG_PINFTY,
        };

        found = quat_sampling_random_ideal_O0_given_norm(
            rng,
            &mut sk.secret_ideal,
            &SEC_DEGREE,
            1,
            &ri_params,
            None,
        );

        if found != 0 {
            found = quat_lideal_prime_norm_reduced_equivalent(
                rng,
                &mut sk.secret_ideal,
                &QUATALG_PINFTY,
                QUAT_REPRESENT_INTEGER_PRIMALITY_ITER,
                QUAT_EQUIV_BOUND_COEFF,
                &EXTREMAL_ORDERS[0].order,
            );
        }

        if found != 0 {
            found = dim2id2iso_arbitrary_isogeny_evaluation(
                rng,
                &mut b_0_two,
                &mut sk.curve,
                &sk.secret_ideal,
            );
        }
    }

    // Compute a deterministic canonical basis with a hint for verification.
    pk.hint_pk = ec_curve_to_basis_2f_to_hint(
        &mut sk.canonical_basis,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
    );

    // Basis-change matrix from the canonical basis to the freshly-evaluated B_0_two.
    change_of_basis_matrix_tate(
        &mut sk.mat_BAcan_to_BA0_two,
        &sk.canonical_basis,
        &b_0_two,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
    );

    // Public key is the codomain curve; clear precomputation flag.
    sqisign_ec::copy_curve(&mut pk.curve, &sk.curve);
    pk.curve.is_A24_computed_and_normalized = false;

    debug_assert_eq!(
        fp2_is_one(&pk.curve.C),
        0xFFFF_FFFF,
        "protocols_keygen: pk.curve.C must be one after copy"
    );

    found
}
