//! `encode_signature.c` port: secret-key (de)serialization. Despite the
//! file name in the C reference, this is only the secret-key path; the
//! public-key and signature encoders live in `sqisign-verify`.
//!
//! The on-wire layout (little-endian throughout) is:
//!
//! ```text
//! +--------------------+--------+--------+--------+--------+--------+
//! | public_key (65 B)  | norm   | g.c0   | g.c1   | g.c2   | g.c3   |
//! +--------------------+--------+--------+--------+--------+--------+
//! | m[0][0] | m[0][1] | m[1][0] | m[1][1]                            |
//! +-------- 4 x TORSION_2POWER_BYTES ----------------------+
//! ```
//!
//! `norm` is encoded unsigned, the four generator coordinates are signed
//! (two's complement, 32-byte width each), and the basis-change matrix
//! cells are unsigned. Sizes are taken verbatim from
//! `precomp/ref/lvl1/include/encoded_sizes.h`.

use sqisign_ec::{ec_curve_to_basis_2f_from_hint, TORSION_EVEN_POWER};
use sqisign_precomp::{EXTREMAL_ORDERS, QUATALG_PINFTY, TORSION_2POWER_BYTES};
use sqisign_quaternion::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_copy_digits, ibz_neg,
    ibz_pow, ibz_sub, ibz_to_digits, quat_lideal_create, quat_lideal_generator, Ibz, QuatAlgElem,
};
use sqisign_verify::{public_key_from_bytes, public_key_to_bytes, PublicKey, PUBLICKEY_BYTES};

use crate::keygen::SecretKey;

/// `FP_ENCODED_BYTES` mirrored from `encoded_sizes.h`.
pub const FP_ENCODED_BYTES: usize = 32;

/// `SECRETKEY_BYTES` mirrored from `encoded_sizes.h`.
pub const SECRETKEY_BYTES: usize =
    PUBLICKEY_BYTES + FP_ENCODED_BYTES * 5 + TORSION_2POWER_BYTES * 4;

const DIGIT_BYTES: usize = 8;

fn encode_digits(enc: &mut [u8], x: &[u64]) {
    // Little-endian host: bytes of each u64 in LE order. The C reference
    // copies the limbs verbatim on LE hosts; we do the same explicitly.
    let nbytes = enc.len();
    let mut off = 0;
    let mut i = 0;
    while off + DIGIT_BYTES <= nbytes {
        let w = x[i];
        enc[off..off + DIGIT_BYTES].copy_from_slice(&w.to_le_bytes());
        i += 1;
        off += DIGIT_BYTES;
    }
    if off < nbytes {
        let w = x[i];
        let bytes = w.to_le_bytes();
        let rem = nbytes - off;
        enc[off..off + rem].copy_from_slice(&bytes[..rem]);
    }
}

fn decode_digits(x: &mut [u64], enc: &[u8]) {
    let nbytes = enc.len();
    let mut off = 0;
    let mut i = 0;
    while off + DIGIT_BYTES <= nbytes {
        let mut b = [0u8; DIGIT_BYTES];
        b.copy_from_slice(&enc[off..off + DIGIT_BYTES]);
        x[i] = u64::from_le_bytes(b);
        i += 1;
        off += DIGIT_BYTES;
    }
    if off < nbytes {
        let mut b = [0u8; DIGIT_BYTES];
        let rem = nbytes - off;
        b[..rem].copy_from_slice(&enc[off..]);
        x[i] = u64::from_le_bytes(b);
    }
    // The C reference zeros the trailing words; our `x` slice is the
    // exact length of the encoded prefix so the caller is responsible
    // for sizing.
}

/// Encode a single `Ibz` to a fixed-width little-endian byte slice. The
/// `sgn` flag enables two's-complement for negative values. Mirrors the
/// reference's static `ibz_to_bytes`.
fn ibz_to_bytes(enc: &mut [u8], x: &Ibz, sgn: bool) {
    let nbytes = enc.len();
    let digits = (nbytes + DIGIT_BYTES - 1) / DIGIT_BYTES;
    let mut d = vec![0u64; digits];
    if ibz_cmp(x, &ibz_const_zero()) >= 0 {
        ibz_to_digits(&mut d, x);
    } else {
        assert!(sgn, "ibz_to_bytes: negative value without sgn flag");
        let mut tmp = Ibz::zero();
        ibz_neg(&mut tmp, x);
        let tmp_clone = tmp.clone();
        ibz_sub(&mut tmp, &tmp_clone, &ibz_const_one());
        ibz_to_digits(&mut d, &tmp);
        for v in d.iter_mut() {
            *v = !*v;
        }
    }
    encode_digits(enc, &d);
}

/// Decode a single `Ibz` from a fixed-width little-endian byte slice.
/// Mirrors the reference's static `ibz_from_bytes`.
fn ibz_from_bytes(x: &mut Ibz, enc: &[u8], sgn: bool) {
    let nbytes = enc.len();
    assert!(nbytes > 0);
    let digits = (nbytes + DIGIT_BYTES - 1) / DIGIT_BYTES;
    let mut d = vec![0u64; digits];
    decode_digits(&mut d, enc);
    if sgn && (enc[nbytes - 1] >> 7) != 0 {
        // Negative: two's-complement decode. The C reference sign-extends
        // the high word; we do the same.
        let s = DIGIT_BYTES - 1 - (digits * DIGIT_BYTES - nbytes);
        debug_assert!(s < DIGIT_BYTES);
        let shift = 8 * s;
        let mask: u64 = (!0u64).wrapping_shr(shift as u32).wrapping_shl(shift as u32);
        d[digits - 1] |= mask;
        for v in d.iter_mut() {
            *v = !*v;
        }
        let mut tmp = Ibz::zero();
        ibz_copy_digits(&mut tmp, &d);
        let cur = tmp.clone();
        ibz_add(&mut tmp, &cur, &ibz_const_one());
        ibz_neg(x, &tmp);
    } else {
        ibz_copy_digits(x, &d);
    }
}

/// `secret_key_to_bytes(enc, sk, pk)`. The C reference encodes the public
/// key first, then the secret-ideal generator and norm, then the basis
/// change matrix. The on-wire size is `SECRETKEY_BYTES`.
pub fn secret_key_to_bytes(enc: &mut [u8; SECRETKEY_BYTES], sk: &SecretKey, pk: &PublicKey) {
    let mut off = public_key_to_bytes(&mut enc[..PUBLICKEY_BYTES], pk);

    // Norm of the secret ideal.
    ibz_to_bytes(
        &mut enc[off..off + FP_ENCODED_BYTES],
        &sk.secret_ideal.norm,
        false,
    );
    off += FP_ENCODED_BYTES;

    // Generator of the secret ideal: encode the four coordinates only
    // (the C reference drops the denominator since it does not change
    // the generated ideal modulo the norm).
    let mut gen = QuatAlgElem::new();
    let _ok = quat_lideal_generator(&mut gen, &sk.secret_ideal, &QUATALG_PINFTY);
    for k in 0..4 {
        ibz_to_bytes(&mut enc[off..off + FP_ENCODED_BYTES], &gen.coord[k], true);
        off += FP_ENCODED_BYTES;
    }

    // Basis change matrix.
    for i in 0..2 {
        for j in 0..2 {
            ibz_to_bytes(
                &mut enc[off..off + TORSION_2POWER_BYTES],
                &sk.mat_BAcan_to_BA0_two[i][j],
                false,
            );
            off += TORSION_2POWER_BYTES;
        }
    }
    debug_assert_eq!(off, SECRETKEY_BYTES);
}

/// `secret_key_from_bytes(sk, pk, enc)`: the inverse of
/// [`secret_key_to_bytes`]. The C reference rebuilds the canonical basis
/// from the recorded hint.
pub fn secret_key_from_bytes(sk: &mut SecretKey, pk: &mut PublicKey, enc: &[u8; SECRETKEY_BYTES]) {
    public_key_from_bytes(pk, &enc[..PUBLICKEY_BYTES]);
    let mut off = PUBLICKEY_BYTES;

    let mut norm = Ibz::zero();
    ibz_from_bytes(&mut norm, &enc[off..off + FP_ENCODED_BYTES], false);
    off += FP_ENCODED_BYTES;

    let mut gen = QuatAlgElem::new();
    for k in 0..4 {
        ibz_from_bytes(&mut gen.coord[k], &enc[off..off + FP_ENCODED_BYTES], true);
        off += FP_ENCODED_BYTES;
    }
    quat_lideal_create(
        &mut sk.secret_ideal,
        &gen,
        &norm,
        &EXTREMAL_ORDERS[0].order,
        &QUATALG_PINFTY,
    );

    for i in 0..2 {
        for j in 0..2 {
            ibz_from_bytes(
                &mut sk.mat_BAcan_to_BA0_two[i][j],
                &enc[off..off + TORSION_2POWER_BYTES],
                false,
            );
            off += TORSION_2POWER_BYTES;
        }
    }
    debug_assert_eq!(off, SECRETKEY_BYTES);

    sk.curve = pk.curve.clone();
    ec_curve_to_basis_2f_from_hint(
        &mut sk.canonical_basis,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
        pk.hint_pk,
    );

    // Suppress unused-import warnings in release if any helpers are
    // consulted only via debug assertions.
    let _ = (ibz_pow, ibz_const_two);
}
