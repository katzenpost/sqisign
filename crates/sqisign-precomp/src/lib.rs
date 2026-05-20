//! SQIsign `precomp`: precomputed constants, copied verbatim from the
//! reference with a re-derivation verification test.
//!
//! Mirrors `the-sqisign/src/precomp`. Ported in **Phase 1, unit 5**.
//!
//! # Scope
//!
//! The reference's lvl1 precomputed-constant tables come from offline Sage
//! scripts and are checked in as enormous `.c` files with `#if
//! GMP_LIMB_BITS` cascades whose only deterministic content is a sequence
//! of `mp_limb_t` arrays. The Rust port records each constant once via the
//! `cdump` harness as a canonical-bytes vector under
//! `vectors/precomp/<NAME>.json`, and rebuilds the typed value at startup
//! from those bytes via [`std::sync::LazyLock`].
//!
//! Verbatim-ness is enforced by the differential test in
//! `tests/precomp_vectors.rs`: the test loads each vector, re-encodes the
//! Rust-constructed value to canonical bytes, and asserts bit-equality
//! with the recorded bytes. Regenerating the vectors from a fresh build
//! of `tools/cdump` therefore round-trips: any drift in the Rust-side
//! reconstruction would diverge from the recorded value, and any drift in
//! the C constants would diverge after regeneration.
//!
//! # Constants ported in this unit
//!
//! From `precomp/ref/lvl1/torsion_constants.c`:
//! [`TWO_TO_SECURITY_BITS`], [`TORSION_PLUS_2POWER`], [`SEC_DEGREE`],
//! [`COM_DEGREE`].
//!
//! From `precomp/ref/lvl1/quaternion_data.c`:
//! [`QUAT_PRIME_COFACTOR`], [`QUATALG_PINFTY`], [`EXTREMAL_ORDERS`],
//! [`CONNECTING_IDEALS`], [`CONJUGATING_ELEMENTS`].
//!
//! From `precomp/ref/lvl1/endomorphism_action.c`:
//! [`CURVES_WITH_ENDOMORPHISMS`].
//!
//! Constants previously absorbed into their consumers
//! (`p_cofactor_for_2f`, `TORSION_EVEN_POWER`, `BASIS_E0_PX`,
//! `BASIS_E0_QX`, `SPLITTING_TRANSFORMS`, etc.) live in the `sqisign-ec`
//! and `sqisign-hd` crates as `const` items.

#![forbid(unsafe_code)]

use std::sync::LazyLock;

use sqisign_ec::{EcBasis, EcCurve, EcPoint};
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_quaternion::dim2::IbzMat2x2;
use sqisign_quaternion::{Ibz, QuatAlg, QuatAlgElem, QuatLeftIdeal};

/// Re-export of [`sqisign_quaternion::QuatPExtremalMaximalOrder`] so the
/// precomputed [`EXTREMAL_ORDERS`] array stays directly usable by the
/// quaternion crate's RNG-driven `quat_represent_integer` and friends
/// without a wrapping conversion.
pub use sqisign_quaternion::QuatPExtremalMaximalOrder;

mod loader;

/// `curve_with_endomorphism_ring_t`: a starting curve together with a
/// precomputed even-torsion basis and the action of the endomorphism
/// ring's quaternion-algebra generators on that basis. Mirrors the C
/// struct in `precomp/ref/lvl1/include/endomorphism_action.h`.
#[derive(Clone, Debug)]
pub struct CurveWithEndomorphismRing {
    pub curve: EcCurve,
    pub basis_even: EcBasis,
    pub action_i: IbzMat2x2,
    pub action_j: IbzMat2x2,
    pub action_k: IbzMat2x2,
    pub action_gen2: IbzMat2x2,
    pub action_gen3: IbzMat2x2,
    pub action_gen4: IbzMat2x2,
}

// ---------------------------------------------------------------------------
// torsion_constants.c
// ---------------------------------------------------------------------------

/// `2^128`. Mirrors the C constant of the same name.
pub static TWO_TO_SECURITY_BITS: LazyLock<Ibz> =
    LazyLock::new(|| loader::load_single_ibz("TWO_TO_SECURITY_BITS"));

/// `2 ^ TORSION_EVEN_POWER`. The full 2-power torsion order.
pub static TORSION_PLUS_2POWER: LazyLock<Ibz> =
    LazyLock::new(|| loader::load_single_ibz("TORSION_PLUS_2POWER"));

/// Degree of the secret isogeny in signing.
pub static SEC_DEGREE: LazyLock<Ibz> = LazyLock::new(|| loader::load_single_ibz("SEC_DEGREE"));

/// Degree of the commitment isogeny in signing.
pub static COM_DEGREE: LazyLock<Ibz> = LazyLock::new(|| loader::load_single_ibz("COM_DEGREE"));

/// `TORSION_2POWER_BYTES`: byte length of a 2-power-torsion order
/// representative. Mirrors `#define TORSION_2POWER_BYTES 32` in
/// `precomp/ref/lvl1/include/torsion_constants.h`.
pub const TORSION_2POWER_BYTES: usize = 32;

// ---------------------------------------------------------------------------
// quaternion_data.c
// ---------------------------------------------------------------------------

/// `QUAT_prime_cofactor`. Mirrors the C constant.
#[allow(non_upper_case_globals)]
pub static QUAT_PRIME_COFACTOR: LazyLock<Ibz> =
    LazyLock::new(|| loader::load_single_ibz("QUAT_prime_cofactor"));

/// The quaternion algebra ramified at `p` and infinity, with `p` the
/// lvl1 prime. Mirrors `const quat_alg_t QUATALG_PINFTY`.
pub static QUATALG_PINFTY: LazyLock<QuatAlg> = LazyLock::new(|| {
    let p = loader::load_named_ibz("QUATALG_PINFTY", &[], 0, "p");
    QuatAlg::init_set(&p)
});

/// `EXTREMAL_ORDERS[0..7]`. Mirrors the C array; index 0 is the standard
/// extremal order `MAXORD_O0`, the remaining six are alternates.
pub static EXTREMAL_ORDERS: LazyLock<[QuatPExtremalMaximalOrder; 7]> =
    LazyLock::new(loader::load_extremal_orders);

/// `CONNECTING_IDEALS[0..7]`. Mirrors the C array; index 0 is the trivial
/// ideal at `MAXORD_O0`, the rest connect to the alternate orders.
pub static CONNECTING_IDEALS: LazyLock<[QuatLeftIdeal; 7]> =
    LazyLock::new(loader::load_connecting_ideals);

/// `CONJUGATING_ELEMENTS[0..7]`. Mirrors the C array.
pub static CONJUGATING_ELEMENTS: LazyLock<[QuatAlgElem; 7]> =
    LazyLock::new(loader::load_conjugating_elements);

/// Number of alternate extremal orders. Mirrors the macro
/// `#define NUM_ALTERNATE_EXTREMAL_ORDERS 6`.
pub const NUM_ALTERNATE_EXTREMAL_ORDERS: usize = 6;

// ---------------------------------------------------------------------------
// endomorphism_action.c
// ---------------------------------------------------------------------------

/// `CURVES_WITH_ENDOMORPHISMS[0..7]`. Mirrors the C array of starting
/// curves; index 0 is `E0` with its precomputed even-torsion basis and
/// the action of the quaternion-algebra generators.
pub static CURVES_WITH_ENDOMORPHISMS: LazyLock<[CurveWithEndomorphismRing; 7]> =
    LazyLock::new(loader::load_curves_with_endomorphisms);

/// Number of alternate starting curves. Mirrors
/// `#define NUM_ALTERNATE_STARTING_CURVES 6`.
pub const NUM_ALTERNATE_STARTING_CURVES: usize = 6;

// Re-export the underlying field-limb count so callers do not need a
// separate `sqisign-gf` import just to size raw buffers.
pub use sqisign_gf::NWORDS_FIELD as PRECOMP_NWORDS_FIELD;

// A short helper to construct a zero Fp2 with the lvl1 limb count; used
// by the loader and re-exported in case downstream consumers find it
// convenient.
pub fn fp2_zero() -> Fp2 {
    Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }
}

pub fn ec_point_zero() -> EcPoint {
    EcPoint::zero()
}
