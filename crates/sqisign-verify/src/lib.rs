//! SQIsign `verification`: the top-level Verify protocol.
//!
//! Mirrors `vendor/the-sqisign/src/verification/ref/lvlx`. Ported in
//! **Phase 2, unit 9**, which closes the verification half of SQIsign.
//!
//! This crate is the public entry point for verifying a SQIsign signature.
//! It is deliberately standalone: it does not depend on `sqisign-sign`,
//! and it does not pull in the LLL / dpe paths of the quaternion module,
//! so Katzenpost mix nodes and clients embed no signing code or
//! floating-point quaternion machinery.
//!
//! ## C-side correspondence
//!
//! The crate mirrors three reference sources:
//!
//! - `verification/ref/lvlx/common.c` (~88 lines): [`public_key_init`],
//!   [`public_key_finalize`], [`hash_to_challenge`].
//! - `verification/ref/lvlx/encode_verification.c` (~220 lines): the
//!   serialization helpers [`signature_to_bytes`], [`signature_from_bytes`],
//!   [`public_key_to_bytes`], [`public_key_from_bytes`], plus the safe
//!   [`signature_decode`] / [`public_key_decode`] front ends.
//! - `verification/ref/lvlx/verify.c` (~309 lines): [`protocols_verify`]
//!   and its private helpers (the canonical-basis check, challenge
//!   recomputation, basis recovery, two-response splitting, and the
//!   2D-isogeny commitment-curve recovery).
//!
//! ## Constants
//!
//! `SECURITY_BITS`, `SQISIGN_RESPONSE_LENGTH`, `HASH_ITERATIONS`,
//! `PUBLICKEY_BYTES`, and `SIGNATURE_BYTES` correspond to the macros in
//! `precomp/ref/lvl1/include/encoded_sizes.h`. They are exposed at the
//! crate root so callers can size buffers without re-deriving them.
#![forbid(unsafe_code)]
#![allow(non_snake_case)]

use sqisign_common::hash::Shake256Absorb;
use sqisign_ec::{
    copy_basis, copy_curve, copy_point, ec_biscalar_mul, ec_curve_init, ec_curve_init_from_a,
    ec_curve_to_basis_2f_from_hint, ec_curve_verify_a, ec_dbl_iter, ec_dbl_iter_basis,
    ec_eval_even, ec_eval_small_chain, ec_is_basis_four_torsion, ec_j_inv, ec_ladder3pt, EcBasis,
    EcCurve, EcIsogEven, EcPoint, NWORDS_ORDER, RADIX, TORSION_EVEN_POWER,
};
use sqisign_gf::{fp2_decode, fp2_encode, fp2_is_one, Fp2, FP2_ENCODED_BYTES};
use sqisign_hd::{
    copy_bases_to_kernel, theta_chain_compute_and_eval_verify, ThetaCoupleCurve, ThetaCouplePoint,
    ThetaKernelCouplePoints, HD_EXTRA_TORSION,
};
use sqisign_mp::{mp_compare, mp_mod_2exp, mp_sub, multiple_mp_shiftl};

// ---------------------------------------------------------------------------
// Level-1 constants from precomp/ref/lvl1/include/encoded_sizes.h.
// ---------------------------------------------------------------------------

/// `SECURITY_BITS`: the level-1 security parameter `lambda` in bits.
/// Mirrors the macro of the same name in the reference's
/// `encoded_sizes.h`.
pub const SECURITY_BITS: usize = 128;

/// `SQIsign_response_length`: the bit length of the SQIsign response,
/// used to size the dim-2 isogeny chain and the basis-change matrix.
pub const SQISIGN_RESPONSE_LENGTH: usize = 126;

/// `HASH_ITERATIONS`: the number of SHAKE256 rejection-sampling
/// iterations in [`hash_to_challenge`].
pub const HASH_ITERATIONS: usize = 64;

/// `PUBLICKEY_BYTES`: serialized public key size.
pub const PUBLICKEY_BYTES: usize = 65;

/// `SIGNATURE_BYTES`: serialized signature size.
pub const SIGNATURE_BYTES: usize = 148;

// ---------------------------------------------------------------------------
// Types: signature_t and public_key_t from verification.h.
// ---------------------------------------------------------------------------

/// Scalar type: `digit_t[NWORDS_ORDER]`. Mirrors the
/// `typedef digit_t scalar_t[NWORDS_ORDER];` in
/// `verification/ref/include/verification.h`.
pub type Scalar = [u64; NWORDS_ORDER];

/// `scalar_mtx_2x2_t`: a 2x2 matrix of [`Scalar`]s, laid out in C as
/// `scalar_t mat[2][2]` (row-major).
pub type ScalarMtx2x2 = [[Scalar; 2]; 2];

/// `signature_t` from `verification.h`.
///
/// The verification protocol consumes the signature as a passive
/// witness: the auxiliary curve `E_aux`, the backtracking and
/// two-response lengths, the canonical-to-arbitrary basis change matrix,
/// the challenge coefficient, and the basis hints for the two canonical
/// bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Signature {
    /// Montgomery A-coefficient for the auxiliary curve.
    pub E_aux_A: Fp2,
    pub backtracking: u8,
    pub two_resp_length: u8,
    /// `mat_Bchall_can_to_B_chall[i][j]` is the 2x2 change-of-basis
    /// matrix from the canonical challenge basis to the arbitrary one.
    pub mat_Bchall_can_to_B_chall: ScalarMtx2x2,
    pub chall_coeff: Scalar,
    pub hint_aux: u8,
    pub hint_chall: u8,
}

impl Signature {
    /// All-zero placeholder. The C reference does not zero-initialize
    /// signatures either; this constructor exists only so callers can
    /// allocate without `unsafe`.
    pub const fn zero() -> Self {
        Signature {
            E_aux_A: Fp2 {
                re: [0u64; sqisign_gf::NWORDS_FIELD],
                im: [0u64; sqisign_gf::NWORDS_FIELD],
            },
            backtracking: 0,
            two_resp_length: 0,
            mat_Bchall_can_to_B_chall: [[[0u64; NWORDS_ORDER]; 2]; 2],
            chall_coeff: [0u64; NWORDS_ORDER],
            hint_aux: 0,
            hint_chall: 0,
        }
    }
}

/// `public_key_t` from `verification.h`. Just the normalized
/// A-coefficient of the Montgomery curve and a basis hint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PublicKey {
    pub curve: EcCurve,
    pub hint_pk: u8,
}

impl PublicKey {
    pub const fn zero() -> Self {
        PublicKey {
            curve: EcCurve::zero(),
            hint_pk: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// common.c: init / finalize / hash_to_challenge.
// ---------------------------------------------------------------------------

/// Mirrors `public_key_init`. The C version just calls `ec_curve_init`
/// on the embedded curve.
pub fn public_key_init(pk: &mut PublicKey) {
    ec_curve_init(&mut pk.curve);
}

/// Mirrors `public_key_finalize`, which is a no-op in the reference but
/// exists as a destructor placeholder. Kept for symmetry.
pub fn public_key_finalize(_pk: &mut PublicKey) {}

/// `hash_to_challenge` from `common.c`: derive the challenge scalar
/// from `(public_key, commitment_curve, message)` via rejection-style
/// SHAKE256 iteration.
///
/// The reference squeezes `2*SECURITY_BITS` bits in each iteration,
/// re-absorbs them, and on the final round squeezes
/// `TORSION_EVEN_POWER - SQIsign_response_length` bits, then reduces
/// modulo `2^SECURITY_BITS`. We mirror this exactly.
pub fn hash_to_challenge(scalar: &mut Scalar, pk: &PublicKey, com_curve: &EcCurve, message: &[u8]) {
    let mut buf = [0u8; 2 * FP2_ENCODED_BYTES];
    let mut j1 = Fp2 {
        re: [0u64; sqisign_gf::NWORDS_FIELD],
        im: [0u64; sqisign_gf::NWORDS_FIELD],
    };
    let mut j2 = j1;
    ec_j_inv(&mut j1, &pk.curve);
    ec_j_inv(&mut j2, com_curve);

    let (first, second) = buf.split_at_mut(FP2_ENCODED_BYTES);
    fp2_encode(first.try_into().unwrap(), &j1);
    fp2_encode(second.try_into().unwrap(), &j2);

    let hash_bytes = (2 * SECURITY_BITS).div_ceil(8);
    let limbs = hash_bytes.div_ceil(8);
    let bits = (2 * SECURITY_BITS) % RADIX;
    let mask: u64 = if bits == 0 {
        !0u64
    } else {
        (!0u64) >> (RADIX - bits)
    };

    let mut tmp = [0u8; (2 * SECURITY_BITS).div_ceil(8)];

    // First iteration: absorb buf and message.
    let mut ctx = Shake256Absorb::new();
    ctx.absorb(&buf);
    ctx.absorb(message);
    let mut squeeze = ctx.finalize();
    squeeze.squeeze(&mut tmp[..hash_bytes]);
    scalar_from_le_bytes(scalar, &tmp[..hash_bytes], limbs, mask);

    // Iterations 2 .. HASH_ITERATIONS - 1.
    for _ in 2..HASH_ITERATIONS {
        scalar_to_le_bytes(&mut tmp[..hash_bytes], scalar, limbs);
        let mut ctx = Shake256Absorb::new();
        ctx.absorb(&tmp[..hash_bytes]);
        let mut squeeze = ctx.finalize();
        squeeze.squeeze(&mut tmp[..hash_bytes]);
        scalar_from_le_bytes(scalar, &tmp[..hash_bytes], limbs, mask);
    }

    // Final iteration: absorb the running scalar, then squeeze
    // (TORSION_EVEN_POWER - SQIsign_response_length) bits.
    scalar_to_le_bytes(&mut tmp[..hash_bytes], scalar, limbs);
    let mut ctx = Shake256Absorb::new();
    ctx.absorb(&tmp[..hash_bytes]);
    let mut squeeze = ctx.finalize();

    let hash_bytes2 = (TORSION_EVEN_POWER - SQISIGN_RESPONSE_LENGTH).div_ceil(8);
    let limbs2 = hash_bytes2.div_ceil(8);
    let bits2 = (TORSION_EVEN_POWER - SQISIGN_RESPONSE_LENGTH) % RADIX;
    let mask2: u64 = if bits2 == 0 {
        !0u64
    } else {
        (!0u64) >> (RADIX - bits2)
    };

    *scalar = [0u64; NWORDS_ORDER];
    let mut tmp2 = [0u8; 64];
    squeeze.squeeze(&mut tmp2[..hash_bytes2]);
    scalar_from_le_bytes(scalar, &tmp2[..hash_bytes2], limbs2, mask2);

    mp_mod_2exp(scalar, SECURITY_BITS as u32);
}

/// Read `bytes` into the low `limbs` words of `scalar` (little-endian)
/// and mask the high word. The remaining words of `scalar` are left
/// alone (the caller pre-zeroes when needed).
fn scalar_from_le_bytes(scalar: &mut Scalar, bytes: &[u8], limbs: usize, mask: u64) {
    debug_assert!(limbs <= NWORDS_ORDER);
    debug_assert!(bytes.len() <= limbs * 8);
    // Zero the limbs we will overwrite, plus the mask target.
    for w in scalar.iter_mut().take(limbs) {
        *w = 0;
    }
    for (i, chunk) in bytes.chunks(8).enumerate() {
        let mut buf = [0u8; 8];
        buf[..chunk.len()].copy_from_slice(chunk);
        scalar[i] = u64::from_le_bytes(buf);
    }
    if limbs > 0 {
        scalar[limbs - 1] &= mask;
    }
}

/// Write the low `limbs` words of `scalar` (little-endian) into `out`.
fn scalar_to_le_bytes(out: &mut [u8], scalar: &Scalar, limbs: usize) {
    debug_assert!(limbs <= NWORDS_ORDER);
    debug_assert!(out.len() <= limbs * 8);
    let total = out.len();
    for (i, limb) in scalar.iter().take(limbs).enumerate() {
        let bytes = limb.to_le_bytes();
        let off = i * 8;
        if off >= total {
            break;
        }
        let take = (total - off).min(8);
        out[off..off + take].copy_from_slice(&bytes[..take]);
    }
}

// ---------------------------------------------------------------------------
// encode_verification.c: signature / public-key serialization.
// ---------------------------------------------------------------------------

/// Encode a [`Scalar`] (or any digit-array) to `nbytes` little-endian
/// bytes. Mirrors the local `encode_digits` helper. The C version
/// branches on `TARGET_BIG_ENDIAN`; the Rust port is unconditionally
/// little-endian because `u64::to_le_bytes` is.
fn encode_digits(enc: &mut [u8], x: &[u64], nbytes: usize) {
    debug_assert!(enc.len() >= nbytes);
    debug_assert!(nbytes <= x.len() * 8);
    for i in 0..nbytes {
        let limb = x[i / 8];
        enc[i] = (limb >> (8 * (i % 8))) as u8;
    }
}

/// Decode `nbytes` little-endian bytes into `ndigits` `u64` limbs,
/// zero-padding the high limbs. Mirrors the local `decode_digits`
/// helper.
fn decode_digits(x: &mut [u64], enc: &[u8], nbytes: usize, ndigits: usize) {
    debug_assert!(nbytes <= ndigits * 8);
    debug_assert!(enc.len() >= nbytes);
    for w in x.iter_mut().take(ndigits) {
        *w = 0;
    }
    for i in 0..nbytes {
        let byte = enc[i] as u64;
        x[i / 8] |= byte << (8 * (i % 8));
    }
}

/// Encode an [`EcCurve`] (the `(A : C)` projective representation
/// normalized to `(A/C : 1)`) into 64 bytes of `Fp2` for the `A`
/// coordinate. Mirrors the local `proj_to_bytes` plus `ec_curve_to_bytes`
/// in `encode_verification.c`.
fn ec_curve_to_bytes(enc: &mut [u8; FP2_ENCODED_BYTES], curve: &EcCurve) {
    let mut tmp = curve.C;
    sqisign_gf::fp2_inv(&mut tmp);
    let tmp_z = tmp;
    sqisign_gf::fp2_mul(&mut tmp, &curve.A, &tmp_z);
    fp2_encode(enc, &tmp);
}

/// Decode an [`EcCurve`] from 64 bytes: read `A`, set `C = 1`, leave
/// `A24` and the cache flag at zero. Mirrors `ec_curve_from_bytes`.
fn ec_curve_from_bytes(curve: &mut EcCurve, enc: &[u8; FP2_ENCODED_BYTES]) {
    *curve = EcCurve::zero();
    fp2_decode(&mut curve.A, enc);
    sqisign_gf::fp2_set_one(&mut curve.C);
}

/// `public_key_to_bytes`: 64 bytes for the curve, 1 byte for the hint.
/// Returns the number of bytes written.
pub fn public_key_to_bytes(enc: &mut [u8], pk: &PublicKey) -> usize {
    assert!(enc.len() >= PUBLICKEY_BYTES);
    let (curve_slot, rest) = enc[..PUBLICKEY_BYTES].split_at_mut(FP2_ENCODED_BYTES);
    let curve_arr: &mut [u8; FP2_ENCODED_BYTES] = curve_slot.try_into().unwrap();
    ec_curve_to_bytes(curve_arr, &pk.curve);
    rest[0] = pk.hint_pk;
    PUBLICKEY_BYTES
}

/// `public_key_from_bytes`: inverse of [`public_key_to_bytes`]. Returns
/// the number of bytes consumed.
pub fn public_key_from_bytes(pk: &mut PublicKey, enc: &[u8]) -> usize {
    assert!(enc.len() >= PUBLICKEY_BYTES);
    let curve_arr: &[u8; FP2_ENCODED_BYTES] = enc[..FP2_ENCODED_BYTES].try_into().unwrap();
    ec_curve_from_bytes(&mut pk.curve, curve_arr);
    pk.hint_pk = enc[FP2_ENCODED_BYTES];
    PUBLICKEY_BYTES
}

/// `signature_to_bytes`: produce the canonical [`SIGNATURE_BYTES`]-long
/// serialization. Layout (matching `encode_verification.c`):
///
/// ```text
/// fp2(E_aux_A) | u8(backtracking) | u8(two_resp_length) |
/// 4 * digits((SQIsign_response_length + 9) / 8) for mat_Bchall_can_to_B_chall |
/// digits(SECURITY_BITS / 8) for chall_coeff |
/// u8(hint_aux) | u8(hint_chall)
/// ```
pub fn signature_to_bytes(enc: &mut [u8], sig: &Signature) {
    assert!(enc.len() >= SIGNATURE_BYTES);
    let mut pos = 0;

    {
        let dst: &mut [u8; FP2_ENCODED_BYTES] =
            (&mut enc[pos..pos + FP2_ENCODED_BYTES]).try_into().unwrap();
        fp2_encode(dst, &sig.E_aux_A);
    }
    pos += FP2_ENCODED_BYTES;

    enc[pos] = sig.backtracking;
    pos += 1;
    enc[pos] = sig.two_resp_length;
    pos += 1;

    let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
    for i in 0..2 {
        for j in 0..2 {
            encode_digits(
                &mut enc[pos..pos + nbytes],
                &sig.mat_Bchall_can_to_B_chall[i][j],
                nbytes,
            );
            pos += nbytes;
        }
    }

    let nbytes = SECURITY_BITS / 8;
    encode_digits(&mut enc[pos..pos + nbytes], &sig.chall_coeff, nbytes);
    pos += nbytes;

    enc[pos] = sig.hint_aux;
    pos += 1;
    enc[pos] = sig.hint_chall;
    pos += 1;

    debug_assert_eq!(pos, SIGNATURE_BYTES);
}

/// `signature_from_bytes`: inverse of [`signature_to_bytes`].
pub fn signature_from_bytes(sig: &mut Signature, enc: &[u8]) {
    assert!(enc.len() >= SIGNATURE_BYTES);
    let mut pos = 0;

    let src: &[u8; FP2_ENCODED_BYTES] = enc[pos..pos + FP2_ENCODED_BYTES].try_into().unwrap();
    fp2_decode(&mut sig.E_aux_A, src);
    pos += FP2_ENCODED_BYTES;

    sig.backtracking = enc[pos];
    pos += 1;
    sig.two_resp_length = enc[pos];
    pos += 1;

    let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
    for i in 0..2 {
        for j in 0..2 {
            decode_digits(
                &mut sig.mat_Bchall_can_to_B_chall[i][j],
                &enc[pos..pos + nbytes],
                nbytes,
                NWORDS_ORDER,
            );
            pos += nbytes;
        }
    }

    let nbytes = SECURITY_BITS / 8;
    decode_digits(
        &mut sig.chall_coeff,
        &enc[pos..pos + nbytes],
        nbytes,
        NWORDS_ORDER,
    );
    pos += nbytes;

    sig.hint_aux = enc[pos];
    pos += 1;
    sig.hint_chall = enc[pos];
    pos += 1;

    debug_assert_eq!(pos, SIGNATURE_BYTES);
}

/// Safe length-checked wrapper around [`signature_from_bytes`].
/// Returns `None` if the input is shorter than [`SIGNATURE_BYTES`].
pub fn signature_decode(bytes: &[u8]) -> Option<Signature> {
    if bytes.len() < SIGNATURE_BYTES {
        return None;
    }
    let mut sig = Signature::zero();
    signature_from_bytes(&mut sig, bytes);
    Some(sig)
}

/// Safe length-checked wrapper around [`public_key_from_bytes`].
/// Returns `None` if the input is shorter than [`PUBLICKEY_BYTES`].
pub fn public_key_decode(bytes: &[u8]) -> Option<PublicKey> {
    if bytes.len() < PUBLICKEY_BYTES {
        return None;
    }
    let mut pk = PublicKey::zero();
    public_key_init(&mut pk);
    public_key_from_bytes(&mut pk, bytes);
    Some(pk)
}

// ---------------------------------------------------------------------------
// verify.c: protocols_verify and its private helpers.
// ---------------------------------------------------------------------------

/// Mirrors the static `check_canonical_basis_change_matrix` in
/// `verify.c`. Verifies that every entry of the change-of-basis matrix
/// is strictly less than
/// `2^(SQIsign_response_length + HD_extra_torsion - backtracking)`.
/// Returns `1` if the matrix is canonical, `0` otherwise.
fn check_canonical_basis_change_matrix(sig: &Signature) -> i32 {
    let mut aux: Scalar = [0u64; NWORDS_ORDER];
    aux[0] = 1;
    let shift =
        (SQISIGN_RESPONSE_LENGTH as i32) + (HD_EXTRA_TORSION as i32) - (sig.backtracking as i32);
    if shift < 0 {
        // The C version assumes shift is non-negative; mirror the
        // failure mode (degenerate aux compares <= every matrix entry).
        return 0;
    }
    // The Rust port of mp_shiftl panics on a zero shift; the C call
    // would be a no-op. Skip the shift when shift == 0.
    if shift > 0 {
        multiple_mp_shiftl(&mut aux, shift as u32);
    }

    let mut ret = 1;
    for i in 0..2 {
        for j in 0..2 {
            if mp_compare(&aux, &sig.mat_Bchall_can_to_B_chall[i][j]) <= 0 {
                ret = 0;
            }
        }
    }
    ret
}

/// Mirrors `compute_challenge_verify`: construct the 2^n isogeny with
/// kernel `P + [chall_coeff]Q` and return its codomain in `e_chall`.
fn compute_challenge_verify(
    e_chall: &mut EcCurve,
    sig: &Signature,
    epk: &EcCurve,
    hint_pk: u8,
) -> i32 {
    let mut phi_chall = EcIsogEven::zero();
    copy_curve(&mut phi_chall.curve, epk);
    phi_chall.length = (TORSION_EVEN_POWER as u32) - (sig.backtracking as u32);

    let mut bas_ea = EcBasis::zero();
    if ec_curve_to_basis_2f_from_hint(
        &mut bas_ea,
        &mut phi_chall.curve,
        TORSION_EVEN_POWER as i32,
        hint_pk,
    ) == 0
    {
        return 0;
    }

    let bas_p = bas_ea.P;
    let bas_q = bas_ea.Q;
    let bas_pmq = bas_ea.PmQ;
    let chall = sig.chall_coeff;
    if ec_ladder3pt(
        &mut phi_chall.kernel,
        &chall,
        &bas_p,
        &bas_q,
        &bas_pmq,
        &phi_chall.curve,
    ) == 0
    {
        return 0;
    }

    let ker = phi_chall.kernel;
    ec_dbl_iter(
        &mut phi_chall.kernel,
        sig.backtracking as i32,
        &ker,
        &mut phi_chall.curve,
    );

    copy_curve(e_chall, &phi_chall.curve);
    // ec_eval_even returns 0 on success in the C reference (it returns
    // the early-exit boolean), so the verify path treats nonzero as
    // failure. Mirror that.
    if ec_eval_even(e_chall, &phi_chall, &mut []) != 0 {
        return 0;
    }
    1
}

/// Mirrors `matrix_scalar_application_even_basis`. The reference
/// dispenses with the per-entry `mp_mod_2exp` reduction because
/// [`check_canonical_basis_change_matrix`] has already verified the
/// entries are canonical.
fn matrix_scalar_application_even_basis(
    bas: &mut EcBasis,
    e: &EcCurve,
    mat: &ScalarMtx2x2,
    f: i32,
) -> i32 {
    let mut tmp_bas = EcBasis::zero();
    copy_basis(&mut tmp_bas, bas);

    if ec_biscalar_mul(&mut bas.P, &mat[0][0], &mat[1][0], f, &tmp_bas, e) == 0 {
        return 0;
    }
    if ec_biscalar_mul(&mut bas.Q, &mat[0][1], &mat[1][1], f, &tmp_bas, e) == 0 {
        return 0;
    }

    let mut scalar0: Scalar = [0u64; NWORDS_ORDER];
    let mut scalar1: Scalar = [0u64; NWORDS_ORDER];
    mp_sub(&mut scalar0, &mat[0][0], &mat[0][1]);
    mp_mod_2exp(&mut scalar0, f as u32);
    mp_sub(&mut scalar1, &mat[1][0], &mat[1][1]);
    mp_mod_2exp(&mut scalar1, f as u32);

    ec_biscalar_mul(&mut bas.PmQ, &scalar0, &scalar1, f, &tmp_bas, e)
}

/// Mirrors `challenge_and_aux_basis_verify`. Recovers the canonical
/// bases on both curves, multiplies each down to the working order, and
/// applies the change-of-basis matrix to the challenge basis.
fn challenge_and_aux_basis_verify(
    b_chall_can: &mut EcBasis,
    b_aux_can: &mut EcBasis,
    e_chall: &mut EcCurve,
    e_aux: &mut EcCurve,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> i32 {
    if ec_curve_to_basis_2f_from_hint(
        b_chall_can,
        e_chall,
        TORSION_EVEN_POWER as i32,
        sig.hint_chall,
    ) == 0
    {
        return 0;
    }

    let shift_chall = (TORSION_EVEN_POWER as i32)
        - pow_dim2_deg_resp
        - (HD_EXTRA_TORSION as i32)
        - (sig.two_resp_length as i32);
    let src = *b_chall_can;
    ec_dbl_iter_basis(b_chall_can, shift_chall, &src, e_chall);

    if ec_curve_to_basis_2f_from_hint(b_aux_can, e_aux, TORSION_EVEN_POWER as i32, sig.hint_aux)
        == 0
    {
        return 0;
    }

    let shift_aux = (TORSION_EVEN_POWER as i32) - pow_dim2_deg_resp - (HD_EXTRA_TORSION as i32);
    let src = *b_aux_can;
    ec_dbl_iter_basis(b_aux_can, shift_aux, &src, e_aux);

    matrix_scalar_application_even_basis(
        b_chall_can,
        e_chall,
        &sig.mat_Bchall_can_to_B_chall,
        pow_dim2_deg_resp + (HD_EXTRA_TORSION as i32) + (sig.two_resp_length as i32),
    )
}

/// Mirrors `two_response_isogeny_verify`. When `two_resp_length != 0`
/// the verifier must perform an extra small 2^r-isogeny step from the
/// challenge curve and push the canonical basis through it.
fn two_response_isogeny_verify(
    e_chall: &mut EcCurve,
    b_chall_can: &mut EcBasis,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> i32 {
    // Mirror mp_is_even: a digit-array is even iff its low limb is even.
    let m00_even = (sig.mat_Bchall_can_to_B_chall[0][0][0] & 1) == 0;
    let m10_even = (sig.mat_Bchall_can_to_B_chall[1][0][0] & 1) == 0;

    let mut ker = EcPoint::zero();
    if m00_even && m10_even {
        copy_point(&mut ker, &b_chall_can.Q);
    } else {
        copy_point(&mut ker, &b_chall_can.P);
    }

    let mut points = [b_chall_can.P, b_chall_can.Q, b_chall_can.PmQ];

    let ker_src = ker;
    ec_dbl_iter(
        &mut ker,
        pow_dim2_deg_resp + (HD_EXTRA_TORSION as i32),
        &ker_src,
        e_chall,
    );

    // ec_eval_small_chain returns 0 on success in the C reference, so
    // the verify path treats nonzero as failure.
    if ec_eval_small_chain(
        e_chall,
        &ker,
        sig.two_resp_length as i32,
        &mut points,
        false,
    ) != 0
    {
        return 0;
    }

    copy_point(&mut b_chall_can.P, &points[0]);
    copy_point(&mut b_chall_can.Q, &points[1]);
    copy_point(&mut b_chall_can.PmQ, &points[2]);
    1
}

/// Mirrors `compute_commitment_curve_verify`. Builds the 2D isogeny
/// from `E_chall x E_aux` and returns the commitment curve as the first
/// component of the codomain.
fn compute_commitment_curve_verify(
    e_com: &mut EcCurve,
    b_chall_can: &EcBasis,
    b_aux_can: &EcBasis,
    e_chall: &EcCurve,
    e_aux: &EcCurve,
    pow_dim2_deg_resp: i32,
) -> i32 {
    let mut e_chall_x_e_aux = ThetaCoupleCurve::zero();
    copy_curve(&mut e_chall_x_e_aux.e1, e_chall);
    copy_curve(&mut e_chall_x_e_aux.e2, e_aux);

    let mut dim_two_ker = ThetaKernelCouplePoints::zero();
    copy_bases_to_kernel(&mut dim_two_ker, b_chall_can, b_aux_can);

    let mut codomain = ThetaCoupleCurve::zero();
    ec_curve_init(&mut codomain.e1);
    ec_curve_init(&mut codomain.e2);

    let codomain_splits = if pow_dim2_deg_resp == 0 {
        // Special case: no dim-2 computation, but we still need to
        // confirm E_chall is supersingular by checking the basis is
        // four-torsion. The C version assumes HD_extra_torsion == 2.
        copy_curve(&mut codomain.e1, &e_chall_x_e_aux.e1);
        copy_curve(&mut codomain.e2, &e_chall_x_e_aux.e2);
        if ec_is_basis_four_torsion(b_chall_can, e_chall) == 0 {
            return 0;
        }
        1
    } else {
        let mut empty: [ThetaCouplePoint; 0] = [];
        theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp as u32,
            &mut e_chall_x_e_aux,
            &dim_two_ker,
            true,
            &mut codomain,
            &mut empty,
        )
    };

    copy_curve(e_com, &codomain.e1);
    codomain_splits
}

/// SQIsign verification. Mirrors `protocols_verify` in `verify.c`.
/// Returns `true` on a valid signature, `false` otherwise.
pub fn protocols_verify(sig: &Signature, pk: &PublicKey, m: &[u8]) -> bool {
    if check_canonical_basis_change_matrix(sig) == 0 {
        return false;
    }

    let pow_dim2_deg_resp =
        (SQISIGN_RESPONSE_LENGTH as i32) - (sig.two_resp_length as i32) - (sig.backtracking as i32);

    if pow_dim2_deg_resp < 0 {
        return false;
    }
    // The dim-2 isogeny embeds a dim-1 isogeny of odd degree, so it
    // can never be of length 2.
    if pow_dim2_deg_resp == 1 {
        return false;
    }

    if ec_curve_verify_a(&pk.curve.A) == 0 {
        return false;
    }

    let mut e_aux = EcCurve::zero();
    if ec_curve_init_from_a(&mut e_aux, &sig.E_aux_A) == 0 {
        return false;
    }

    // Reference asserts the public key is in (A : 1) form with no
    // A24 cache. The Rust counterpart relies on the decoder to
    // produce exactly that shape; assert in debug builds.
    debug_assert!(fp2_is_one(&pk.curve.C) == 0xFFFFFFFF);
    debug_assert!(!pk.curve.is_A24_computed_and_normalized);

    let mut e_chall = EcCurve::zero();
    if compute_challenge_verify(&mut e_chall, sig, &pk.curve, pk.hint_pk) == 0 {
        return false;
    }

    let mut b_chall_can = EcBasis::zero();
    let mut b_aux_can = EcBasis::zero();
    if challenge_and_aux_basis_verify(
        &mut b_chall_can,
        &mut b_aux_can,
        &mut e_chall,
        &mut e_aux,
        sig,
        pow_dim2_deg_resp,
    ) == 0
    {
        return false;
    }

    if sig.two_resp_length > 0
        && two_response_isogeny_verify(&mut e_chall, &mut b_chall_can, sig, pow_dim2_deg_resp) == 0
    {
        return false;
    }

    let mut e_com = EcCurve::zero();
    if compute_commitment_curve_verify(
        &mut e_com,
        &b_chall_can,
        &b_aux_can,
        &e_chall,
        &e_aux,
        pow_dim2_deg_resp,
    ) == 0
    {
        return false;
    }

    let mut chk_chall: Scalar = [0u64; NWORDS_ORDER];
    hash_to_challenge(&mut chk_chall, pk, &e_com, m);

    mp_compare(&sig.chall_coeff, &chk_chall) == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signature_round_trip_zero() {
        let mut sig = Signature::zero();
        sig.backtracking = 5;
        sig.two_resp_length = 3;
        sig.hint_aux = 0x7f;
        sig.hint_chall = 0x42;
        sig.chall_coeff[0] = 0xdeadbeef;
        sig.mat_Bchall_can_to_B_chall[0][0][0] = 1;
        sig.mat_Bchall_can_to_B_chall[1][1][0] = 0x1234_5678;

        let mut buf = [0u8; SIGNATURE_BYTES];
        signature_to_bytes(&mut buf, &sig);
        let mut decoded = Signature::zero();
        signature_from_bytes(&mut decoded, &buf);
        assert_eq!(sig, decoded);
    }

    #[test]
    fn signature_decode_rejects_short() {
        let buf = [0u8; SIGNATURE_BYTES - 1];
        assert!(signature_decode(&buf).is_none());
    }

    #[test]
    fn public_key_round_trip_zero() {
        let mut pk = PublicKey::zero();
        public_key_init(&mut pk);
        pk.hint_pk = 0xab;
        // Put something nonzero into the A coordinate.
        sqisign_gf::fp2_set_one(&mut pk.curve.A);
        let mut buf = [0u8; PUBLICKEY_BYTES];
        public_key_to_bytes(&mut buf, &pk);
        let decoded = public_key_decode(&buf).expect("decode");
        assert_eq!(decoded.hint_pk, pk.hint_pk);
        assert_eq!(decoded.curve.A, pk.curve.A);
        assert_eq!(fp2_is_one(&decoded.curve.C), 0xFFFFFFFF);
    }

    #[test]
    fn public_key_decode_rejects_short() {
        let buf = [0u8; PUBLICKEY_BYTES - 1];
        assert!(public_key_decode(&buf).is_none());
    }

    #[test]
    fn encode_decode_digits_are_inverses() {
        let x: [u64; 4] = [0x0102030405060708, 0x1112131415161718, 0xdeadbeef, 0];
        let mut buf = [0u8; 17];
        encode_digits(&mut buf, &x, 17);
        let mut y = [0u64; 4];
        decode_digits(&mut y, &buf, 17, 4);
        assert_eq!(x[0], y[0]);
        assert_eq!(x[1], y[1]);
        // High limb only partially restored.
        assert_eq!(y[2] & 0xff, x[2] & 0xff);
        assert_eq!(y[3], 0);
    }
}
