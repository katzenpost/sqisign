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
//!
//! The reference exposes the same shape for SHAKE128
//! (`void shake128(uint8_t *output, size_t outlen, const uint8_t *input,
//! size_t inlen)`); [`shake128`] / [`shake128_vec`] mirror it. The two
//! differ only in Keccak rate (168 vs 136 bytes); both are the same audited
//! `sha3` Keccak underneath, equivalence proven by the committed C-derived
//! vectors (`tests/shake128_vectors.rs`, `tests/shake256_vectors.rs`).

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

/// SHAKE128 of `input`, squeezed into `output`. Same one-shot XOF contract
/// as [`shake256`], differing only in Keccak rate (168 bytes); the prefix
/// guarantee the reference relies on holds identically.
pub fn shake128(input: &[u8], output: &mut [u8]) {
    let mut hasher = sha3::Shake128::default();
    hasher.update(input);
    hasher.finalize_xof().read(output);
}

/// Allocating convenience over [`shake128`]: returns `out_len` squeezed
/// bytes.
pub fn shake128_vec(input: &[u8], out_len: usize) -> Vec<u8> {
    let mut out = vec![0u8; out_len];
    shake128(input, &mut out);
    out
}

/// Incremental SHAKE256, mirroring the reference's
/// `shake256_inc_init` / `_absorb` / `_finalize` / `_squeeze` contract:
/// construct, [`absorb`](Shake256Absorb::absorb) any number of times,
/// [`finalize`](Shake256Absorb::finalize) once, then
/// [`squeeze`](Shake256Squeeze::squeeze) any number of times. The type
/// system enforces the one-way absorb -> squeeze transition the C API only
/// documents in a comment.
#[derive(Clone, Default)]
pub struct Shake256Absorb(sha3::Shake256);

/// The squeeze phase of an incremental SHAKE256, after finalize.
#[derive(Clone)]
pub struct Shake256Squeeze(sha3::Shake256Reader);

impl Shake256Absorb {
    /// A fresh sponge, equivalent to `shake256_inc_init`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Absorb more input, equivalent to one `shake256_inc_absorb` call.
    /// Calling it `k` times with chunks whose concatenation is `m` is
    /// indistinguishable from one call with `m`.
    pub fn absorb(&mut self, input: &[u8]) {
        self.0.update(input);
    }

    /// Finalize for squeezing, equivalent to `shake256_inc_finalize`.
    pub fn finalize(self) -> Shake256Squeeze {
        Shake256Squeeze(self.0.finalize_xof())
    }
}

impl Shake256Squeeze {
    /// Squeeze `out.len()` more bytes, equivalent to one
    /// `shake256_inc_squeeze`. The byte stream is continuous across calls,
    /// so chunking the squeeze does not change the output.
    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.0.read(out);
    }
}

/// Incremental SHAKE128. Identical contract to [`Shake256Absorb`], differing
/// only in Keccak rate (168 bytes).
#[derive(Clone, Default)]
pub struct Shake128Absorb(sha3::Shake128);

/// The squeeze phase of an incremental SHAKE128, after finalize.
#[derive(Clone)]
pub struct Shake128Squeeze(sha3::Shake128Reader);

impl Shake128Absorb {
    /// A fresh sponge, equivalent to `shake128_inc_init`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Absorb more input, equivalent to one `shake128_inc_absorb` call.
    pub fn absorb(&mut self, input: &[u8]) {
        self.0.update(input);
    }

    /// Finalize for squeezing, equivalent to `shake128_inc_finalize`.
    pub fn finalize(self) -> Shake128Squeeze {
        Shake128Squeeze(self.0.finalize_xof())
    }
}

impl Shake128Squeeze {
    /// Squeeze `out.len()` more bytes, equivalent to one
    /// `shake128_inc_squeeze`.
    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.0.read(out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // FIPS-202 / NIST: SHAKE256("") first bytes. Also the value our C
    // dump recorded for vector id 0, so this anchors the unit test to the
    // same ground truth as the differential suite.
    #[test]
    fn shake256_empty_input_known_answer() {
        let got = shake256_vec(b"", 32);
        assert_eq!(
            hex::encode(got),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
    }

    // FIPS-202 / NIST: SHAKE128("") first bytes. Anchors the unit test to
    // the same ground truth our C dump recorded for shake128 vector id 0.
    #[test]
    fn shake128_empty_input_known_answer() {
        let got = shake128_vec(b"", 32);
        assert_eq!(
            hex::encode(got),
            "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26"
        );
    }

    // Chunked absorb then chunked squeeze must reproduce the one-shot
    // result exactly; anchored to the same FIPS-202 empty-input answers.
    #[test]
    fn incremental_chunking_matches_one_shot() {
        let mut a = Shake256Absorb::new();
        a.absorb(b"");
        a.absorb(b"");
        let mut sq = a.finalize();
        let mut out = [0u8; 32];
        sq.squeeze(&mut out[..1]);
        sq.squeeze(&mut out[1..32]);
        assert_eq!(
            hex::encode(out),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );

        let mut a = Shake128Absorb::new();
        a.absorb(b"katzen");
        a.absorb(b"post");
        let mut sq = a.finalize();
        let mut inc = [0u8; 100];
        sq.squeeze(&mut inc[..7]);
        sq.squeeze(&mut inc[7..]);
        assert_eq!(inc.to_vec(), shake128_vec(b"katzenpost", 100));
    }

    #[test]
    fn vec_and_buffer_forms_agree() {
        let v = shake256_vec(b"katzenpost", 100);
        let mut b = [0u8; 100];
        shake256(b"katzenpost", &mut b);
        assert_eq!(v, b);

        let v = shake128_vec(b"katzenpost", 100);
        let mut b = [0u8; 100];
        shake128(b"katzenpost", &mut b);
        assert_eq!(v, b);
    }
}
