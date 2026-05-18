//! SHA-3 / SHAKE hashing.
//!
//! Boundary parity with the C reference: the reference exposes
//!
//! ```c
//! void shake256(uint8_t *output, size_t outlen,
//!               const uint8_t *input, size_t inlen);
//! ```
//!
//! a one-shot SHAKE256 extendable-output function. [`shake256`] is the
//! buffer-filling form with identical semantics (no allocation, caller owns
//! the output length); [`shake256_vec`] is an allocating convenience.

use sha3::digest::{ExtendableOutput, Update, XofReader};

/// SHAKE256 of `input`, squeezed into `output` (its full length is the XOF
/// output length, exactly as the reference's `outlen`).
///
/// SHAKE is an extendable-output function, so the first `n` bytes of any
/// longer output are a prefix of a shorter one; callers relying on that
/// (the reference does) get the same guarantee here.
pub fn shake256(input: &[u8], output: &mut [u8]) {
    let mut hasher = sha3::Shake256::default();
    hasher.update(input);
    hasher.finalize_xof().read(output);
}

/// Allocating convenience over [`shake256`]: returns `out_len` squeezed
/// bytes.
pub fn shake256_vec(input: &[u8], out_len: usize) -> Vec<u8> {
    let mut out = vec![0u8; out_len];
    shake256(input, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // FIPS-202 / NIST: SHAKE256("") first bytes. Also the value our C
    // dump recorded for vector id 0, so this anchors the unit test to the
    // same ground truth as the differential suite.
    #[test]
    fn empty_input_known_answer() {
        let got = shake256_vec(b"", 32);
        assert_eq!(
            hex::encode(got),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
    }

    #[test]
    fn vec_and_buffer_forms_agree() {
        let v = shake256_vec(b"katzenpost", 100);
        let mut b = [0u8; 100];
        shake256(b"katzenpost", &mut b);
        assert_eq!(v, b);
    }
}
