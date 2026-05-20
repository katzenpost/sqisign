//! SQIsign C ABI surface (Katzenpost-facing).
//!
//! This crate exposes a small C ABI over [`sqisign-verify`] and
//! [`sqisign-sign`]. It is the surface non-Rust consumers (the
//! `bindings/go` cgo wrapper, future Python bindings, the Katzenpost
//! PKI's hybrid Ed25519+SQIsign verifier in mix nodes and clients)
//! talk to.
//!
//! ## ABI
//!
//! Verification:
//!
//! - [`sqisign_lvl1_verify`]: pass-or-fail verification.
//!
//! Keypair generation and signing come in two flavours; both produce
//! the same wire output for the same byte stream, the difference is
//! which RNG drives them:
//!
//! - [`sqisign_lvl1_keygen`] / [`sqisign_lvl1_sign`] take a 48-byte
//!   entropy block, seed a NIST AES-256 CTR-DRBG with it, and drive
//!   the algorithm from that. These exist so KAT replay against the
//!   upstream NIST vectors works bit-for-bit; production callers
//!   should not use them.
//!
//! - [`sqisign_lvl1_keygen_with_rng`] /
//!   [`sqisign_lvl1_sign_with_rng`] take a caller-supplied
//!   `(callback, context)` pair: every byte of randomness the
//!   algorithm demands comes from the callback. No NIST DRBG, no
//!   hidden state. This is the surface non-Rust production callers
//!   are expected to use; it mirrors the highctidh callback shape so
//!   the same `mattn/go-pointer`-style trick (uintptr_t context
//!   smuggling a Go pointer) works here too.
//!
//! Every entry point returns `1` on success and `0` on any failure
//! (verification false, length mismatch, null pointer where data was
//! required, null callback, internal Rust panic caught at the
//! boundary, algorithmic non-success status). The C ABI never returns
//! any other value and never lets a Rust panic cross the boundary.
//!
//! ## Safety
//!
//! `#![forbid(unsafe_code)]` would be ideal, but every FFI entry point
//! in this crate must turn caller-supplied `*const u8` pointers into
//! Rust slices, which is only possible through `unsafe`. The unsafe
//! surface is limited to that pointer-to-slice conversion, a
//! `catch_unwind` around the body, and the indirect call through the
//! caller-supplied RNG callback function pointer in the `_with_rng`
//! entries. The verification crate itself forbids unsafe code.

use core::ffi::c_int;
use core::panic::AssertUnwindSafe;
use std::panic;

use sqisign_common::RngSource;
use sqisign_verify::{protocols_verify, public_key_decode, signature_decode};

/// Wire size of a serialized level-1 SQIsign public key, in bytes.
///
/// Mirrors `sqisign_verify::PUBLICKEY_BYTES`. Re-exported here so C
/// callers see a single, ABI-stable constant.
pub const SQISIGN_LVL1_PUBLIC_KEY_BYTES: usize = sqisign_verify::PUBLICKEY_BYTES;

/// Wire size of a serialized level-1 SQIsign signature, in bytes.
///
/// Mirrors `sqisign_verify::SIGNATURE_BYTES`.
pub const SQISIGN_LVL1_SIGNATURE_BYTES: usize = sqisign_verify::SIGNATURE_BYTES;

/// Verify a level-1 SQIsign signature.
///
/// # Parameters
///
/// - `sig` / `sig_len`: signature bytes; must be exactly
///   [`SQISIGN_LVL1_SIGNATURE_BYTES`] long.
/// - `pk` / `pk_len`: public key bytes; must be exactly
///   [`SQISIGN_LVL1_PUBLIC_KEY_BYTES`] long.
/// - `msg` / `msg_len`: message that was signed. Zero-length messages are
///   permitted; `msg` may be null iff `msg_len == 0`.
///
/// # Return value
///
/// - `1` on a valid signature.
/// - `0` on any failure: verification returned false, lengths do not
///   match, a non-empty buffer pointer was null, or the verifier panicked
///   (which would otherwise be undefined behaviour across an FFI
///   boundary). The function never returns any other value.
///
/// # Safety
///
/// The caller must ensure that `sig`, `pk`, and `msg` point to readable
/// regions of at least the corresponding `*_len` bytes (or are null and
/// the matching length is zero, for `msg` only; `sig` and `pk` must
/// always reference a populated buffer of the exact required length).
/// The buffers are not retained past the call.
///
/// # No-panic guarantee
///
/// The implementation wraps the verifier in `catch_unwind` and maps any
/// panic to a `0` return value. This makes the function safe to call
/// from C even if a future bug in the verify crate triggers a panic.
#[no_mangle]
pub unsafe extern "C" fn sqisign_lvl1_verify(
    sig: *const u8,
    sig_len: usize,
    pk: *const u8,
    pk_len: usize,
    msg: *const u8,
    msg_len: usize,
) -> c_int {
    // Length gate. The verify crate's decoders themselves accept any
    // buffer >= the exact size, so we tighten the contract here to make
    // length mismatches a hard failure rather than silently truncating.
    if sig_len != SQISIGN_LVL1_SIGNATURE_BYTES {
        return 0;
    }
    if pk_len != SQISIGN_LVL1_PUBLIC_KEY_BYTES {
        return 0;
    }

    // Null-pointer gate. A null pointer with a nonzero length would be
    // unsound to dereference; a null `msg` with `msg_len == 0` is fine
    // because we build an empty slice in that case.
    if sig.is_null() || pk.is_null() {
        return 0;
    }
    if msg.is_null() && msg_len != 0 {
        return 0;
    }

    // Materialize slices. Safe under the caller's preconditions above.
    let sig_slice = unsafe { core::slice::from_raw_parts(sig, sig_len) };
    let pk_slice = unsafe { core::slice::from_raw_parts(pk, pk_len) };
    let msg_slice: &[u8] = if msg_len == 0 {
        &[]
    } else {
        unsafe { core::slice::from_raw_parts(msg, msg_len) }
    };

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let Some(sig) = signature_decode(sig_slice) else {
            return false;
        };
        let Some(pk) = public_key_decode(pk_slice) else {
            return false;
        };
        protocols_verify(&sig, &pk, msg_slice)
    }));

    match result {
        Ok(true) => 1,
        Ok(false) => 0,
        Err(_) => 0,
    }
}

// Signing FFI: keypair and sign. Both consume a 48-byte caller-supplied
// entropy block to seed the KAT-compatible CTR-DRBG. Production callers
// who want their own RNG should use the Rust-level
// `sqisign_sign::protocols_keygen` / `protocols_sign` entry points
// directly with their own `RngSource` implementation.

/// Wire size of a serialized level-1 SQIsign secret key, in bytes.
///
/// Mirrors `sqisign_sign::SECRETKEY_BYTES`.
pub const SQISIGN_LVL1_SECRET_KEY_BYTES: usize = sqisign_sign::SECRETKEY_BYTES;

/// Length of the entropy block consumed by [`sqisign_lvl1_keygen`] and
/// [`sqisign_lvl1_sign`]. Matches the NIST KAT format (48 bytes).
pub const SQISIGN_LVL1_ENTROPY_BYTES: usize = 48;

/// Generate a level-1 SQIsign keypair from a 48-byte entropy seed.
///
/// # Parameters
///
/// - `pk` / `pk_len`: output public key buffer; must be exactly
///   [`SQISIGN_LVL1_PUBLIC_KEY_BYTES`] long.
/// - `sk` / `sk_len`: output secret key buffer; must be exactly
///   [`SQISIGN_LVL1_SECRET_KEY_BYTES`] long.
/// - `entropy` / `entropy_len`: 48 bytes used to seed the
///   KAT-compatible CTR-DRBG.
///
/// # Return value
///
/// - `1` on success.
/// - `0` on any failure: lengths do not match, any pointer is null, or
///   the keypair routine panicked.
///
/// # Safety
///
/// All pointers must reference readable / writable buffers of the
/// indicated length. The buffers are not retained past the call.
#[no_mangle]
pub unsafe extern "C" fn sqisign_lvl1_keygen(
    pk: *mut u8,
    pk_len: usize,
    sk: *mut u8,
    sk_len: usize,
    entropy: *const u8,
    entropy_len: usize,
) -> c_int {
    if pk_len != SQISIGN_LVL1_PUBLIC_KEY_BYTES
        || sk_len != SQISIGN_LVL1_SECRET_KEY_BYTES
        || entropy_len != SQISIGN_LVL1_ENTROPY_BYTES
    {
        return 0;
    }
    if pk.is_null() || sk.is_null() || entropy.is_null() {
        return 0;
    }

    let entropy_slice = unsafe { core::slice::from_raw_parts(entropy, entropy_len) };
    let mut entropy_arr = [0u8; SQISIGN_LVL1_ENTROPY_BYTES];
    entropy_arr.copy_from_slice(entropy_slice);

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut drbg = sqisign_common::CtrDrbg::new(&entropy_arr, None);
        let mut pkv = sqisign_verify::PublicKey::zero();
        let mut skv = sqisign_sign::SecretKey::new();
        let ok = sqisign_sign::protocols_keygen(&mut drbg, &mut pkv, &mut skv);
        if ok != 1 {
            return false;
        }
        let mut pk_buf = vec![0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        sqisign_verify::public_key_to_bytes(&mut pk_buf, &pkv);
        let mut sk_buf = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
        sqisign_sign::secret_key_to_bytes(&mut sk_buf, &skv, &pkv);
        unsafe {
            core::ptr::copy_nonoverlapping(pk_buf.as_ptr(), pk, SQISIGN_LVL1_PUBLIC_KEY_BYTES);
            core::ptr::copy_nonoverlapping(sk_buf.as_ptr(), sk, SQISIGN_LVL1_SECRET_KEY_BYTES);
        }
        true
    }));
    match result {
        Ok(true) => 1,
        Ok(false) => 0,
        Err(_) => 0,
    }
}

/// Sign a message with a level-1 SQIsign secret key, seeded by the given
/// 48-byte entropy block. The output `sig` is the signature only (length
/// [`SQISIGN_LVL1_SIGNATURE_BYTES`]); the NIST `sm = signature || msg`
/// concatenation is the caller's responsibility.
///
/// # Parameters
///
/// - `sig` / `sig_len`: output signature buffer; must be exactly
///   [`SQISIGN_LVL1_SIGNATURE_BYTES`] long.
/// - `msg` / `msg_len`: input message.
/// - `sk` / `sk_len`: secret key; must be exactly
///   [`SQISIGN_LVL1_SECRET_KEY_BYTES`] long.
/// - `entropy` / `entropy_len`: 48 bytes used to seed the
///   KAT-compatible CTR-DRBG.
///
/// # Return value
///
/// - `1` on success.
/// - `0` on any failure (length mismatch, null pointer with nonzero
///   length, signing returned failure, or panic).
///
/// # Safety
///
/// All pointers must reference readable / writable buffers of the
/// indicated length.
#[no_mangle]
pub unsafe extern "C" fn sqisign_lvl1_sign(
    sig: *mut u8,
    sig_len: usize,
    msg: *const u8,
    msg_len: usize,
    sk: *const u8,
    sk_len: usize,
    entropy: *const u8,
    entropy_len: usize,
) -> c_int {
    if sig_len != SQISIGN_LVL1_SIGNATURE_BYTES
        || sk_len != SQISIGN_LVL1_SECRET_KEY_BYTES
        || entropy_len != SQISIGN_LVL1_ENTROPY_BYTES
    {
        return 0;
    }
    if sig.is_null() || sk.is_null() || entropy.is_null() {
        return 0;
    }
    if msg.is_null() && msg_len != 0 {
        return 0;
    }
    let sk_slice = unsafe { core::slice::from_raw_parts(sk, sk_len) };
    let entropy_slice = unsafe { core::slice::from_raw_parts(entropy, entropy_len) };
    let msg_slice: &[u8] = if msg_len == 0 {
        &[]
    } else {
        unsafe { core::slice::from_raw_parts(msg, msg_len) }
    };
    let mut entropy_arr = [0u8; SQISIGN_LVL1_ENTROPY_BYTES];
    entropy_arr.copy_from_slice(entropy_slice);
    let mut sk_arr = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
    sk_arr.copy_from_slice(sk_slice);
    let msg_vec = msg_slice.to_vec();

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut drbg = sqisign_common::CtrDrbg::new(&entropy_arr, None);
        let mut pkv = sqisign_verify::PublicKey::zero();
        let mut skv = sqisign_sign::SecretKey::new();
        sqisign_sign::secret_key_from_bytes(&mut skv, &mut pkv, &sk_arr);
        let mut sigv = sqisign_verify::Signature::zero();
        let ok = sqisign_sign::protocols_sign(&mut drbg, &mut sigv, &pkv, &mut skv, &msg_vec);
        if ok != 1 {
            return false;
        }
        let mut sig_buf = [0u8; SQISIGN_LVL1_SIGNATURE_BYTES];
        sqisign_verify::signature_to_bytes(&mut sig_buf, &sigv);
        unsafe {
            core::ptr::copy_nonoverlapping(sig_buf.as_ptr(), sig, SQISIGN_LVL1_SIGNATURE_BYTES);
        }
        true
    }));
    match result {
        Ok(true) => 1,
        Ok(false) => 0,
        Err(_) => 0,
    }
}

// Caller-supplied-RNG signing FFI. These mirror the entropy-based
// entries above but take a `(callback, context)` pair instead of a
// 48-byte seed: every byte of randomness the algorithm demands comes
// from the callback. No NIST DRBG, no hidden state. The callback
// shape matches highctidh's `ctidh_fillrandom`, so the same
// `mattn/go-pointer`-style context-smuggling works on the Go side.

/// C-callable RNG callback. Implementations must, on each call, fill
/// exactly `len` bytes at `out` from the caller's randomness source.
/// `context` is an opaque value the FFI threads through unchanged;
/// callers typically use it to carry a pointer to RNG state through
/// the C boundary (e.g. a Go `gopointer.Save` handle).
pub type SqisignFillRandomFn =
    unsafe extern "C" fn(out: *mut u8, len: usize, context: usize);

/// Bridge from a C function pointer to the workspace's [`RngSource`]
/// trait. Holds no state of its own beyond the callback and context.
struct CallbackRng {
    callback: SqisignFillRandomFn,
    context: usize,
}

impl RngSource for CallbackRng {
    fn fill(&mut self, out: &mut [u8]) {
        // SAFETY: the C ABI contract is that the callback fills
        // exactly `len` bytes at `out` and does not retain the
        // pointer. Zero-length requests are still valid and must not
        // touch memory; we leave it to the callback to handle.
        unsafe {
            (self.callback)(out.as_mut_ptr(), out.len(), self.context);
        }
    }
}

/// Drive [`sqisign_sign::protocols_keygen`] and serialize the
/// resulting keypair into the caller's `pk` / `sk` buffers. Returns
/// `false` if the algorithm itself reports failure.
fn keygen_into_buffers<R: RngSource>(rng: &mut R, pk: *mut u8, sk: *mut u8) -> bool {
    let mut pkv = sqisign_verify::PublicKey::zero();
    let mut skv = sqisign_sign::SecretKey::new();
    if sqisign_sign::protocols_keygen(rng, &mut pkv, &mut skv) != 1 {
        return false;
    }
    let mut pk_buf = vec![0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
    sqisign_verify::public_key_to_bytes(&mut pk_buf, &pkv);
    let mut sk_buf = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
    sqisign_sign::secret_key_to_bytes(&mut sk_buf, &skv, &pkv);
    unsafe {
        core::ptr::copy_nonoverlapping(pk_buf.as_ptr(), pk, SQISIGN_LVL1_PUBLIC_KEY_BYTES);
        core::ptr::copy_nonoverlapping(sk_buf.as_ptr(), sk, SQISIGN_LVL1_SECRET_KEY_BYTES);
    }
    true
}

/// Drive [`sqisign_sign::protocols_sign`] and serialize the resulting
/// signature into the caller's `sig` buffer. Returns `false` if the
/// algorithm itself reports failure.
fn sign_into_buffer<R: RngSource>(
    rng: &mut R,
    sig: *mut u8,
    sk_bytes: &[u8; SQISIGN_LVL1_SECRET_KEY_BYTES],
    msg: &[u8],
) -> bool {
    let mut pkv = sqisign_verify::PublicKey::zero();
    let mut skv = sqisign_sign::SecretKey::new();
    sqisign_sign::secret_key_from_bytes(&mut skv, &mut pkv, sk_bytes);
    let mut sigv = sqisign_verify::Signature::zero();
    if sqisign_sign::protocols_sign(rng, &mut sigv, &pkv, &mut skv, msg) != 1 {
        return false;
    }
    let mut sig_buf = [0u8; SQISIGN_LVL1_SIGNATURE_BYTES];
    sqisign_verify::signature_to_bytes(&mut sig_buf, &sigv);
    unsafe {
        core::ptr::copy_nonoverlapping(sig_buf.as_ptr(), sig, SQISIGN_LVL1_SIGNATURE_BYTES);
    }
    true
}

/// Generate a level-1 SQIsign keypair, taking every byte of
/// randomness from a caller-supplied callback.
///
/// # Parameters
///
/// - `pk` / `pk_len`: output public key buffer; must be exactly
///   [`SQISIGN_LVL1_PUBLIC_KEY_BYTES`] long.
/// - `sk` / `sk_len`: output secret key buffer; must be exactly
///   [`SQISIGN_LVL1_SECRET_KEY_BYTES`] long.
/// - `fill_random`: function pointer the algorithm calls every time
///   it needs randomness. Must be non-null.
/// - `rng_context`: opaque value passed through to every call of
///   `fill_random`. Not interpreted by this crate.
///
/// # Return value
///
/// - `1` on success.
/// - `0` on any failure: length mismatch, null pointer (including a
///   null `fill_random`), keypair routine returned non-success, or an
///   internal panic caught at the FFI boundary.
///
/// # Safety
///
/// `pk` and `sk` must reference writable buffers of exactly the
/// indicated length; `fill_random` must satisfy the documented
/// contract on [`SqisignFillRandomFn`]. Buffers are not retained past
/// the call.
#[no_mangle]
pub unsafe extern "C" fn sqisign_lvl1_keygen_with_rng(
    pk: *mut u8,
    pk_len: usize,
    sk: *mut u8,
    sk_len: usize,
    fill_random: Option<SqisignFillRandomFn>,
    rng_context: usize,
) -> c_int {
    if pk_len != SQISIGN_LVL1_PUBLIC_KEY_BYTES || sk_len != SQISIGN_LVL1_SECRET_KEY_BYTES {
        return 0;
    }
    if pk.is_null() || sk.is_null() {
        return 0;
    }
    let Some(callback) = fill_random else { return 0 };

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut rng = CallbackRng {
            callback,
            context: rng_context,
        };
        keygen_into_buffers(&mut rng, pk, sk)
    }));
    match result {
        Ok(true) => 1,
        _ => 0,
    }
}

/// Sign a message with a level-1 SQIsign secret key, taking every
/// byte of randomness from a caller-supplied callback. The output
/// `sig` is the signature alone (length
/// [`SQISIGN_LVL1_SIGNATURE_BYTES`]); the NIST `sm = signature || msg`
/// concatenation is the caller's responsibility.
///
/// # Parameters
///
/// - `sig` / `sig_len`: output signature buffer; must be exactly
///   [`SQISIGN_LVL1_SIGNATURE_BYTES`] long.
/// - `msg` / `msg_len`: input message. `msg` may be null iff
///   `msg_len == 0`.
/// - `sk` / `sk_len`: secret key; must be exactly
///   [`SQISIGN_LVL1_SECRET_KEY_BYTES`] long.
/// - `fill_random`: function pointer the algorithm calls every time
///   it needs randomness. Must be non-null.
/// - `rng_context`: opaque value passed through to every call of
///   `fill_random`.
///
/// # Return value
///
/// - `1` on success.
/// - `0` on any failure.
///
/// # Safety
///
/// All non-null pointers must reference readable / writable buffers
/// of the indicated length; `fill_random` must satisfy the documented
/// contract on [`SqisignFillRandomFn`].
#[no_mangle]
pub unsafe extern "C" fn sqisign_lvl1_sign_with_rng(
    sig: *mut u8,
    sig_len: usize,
    msg: *const u8,
    msg_len: usize,
    sk: *const u8,
    sk_len: usize,
    fill_random: Option<SqisignFillRandomFn>,
    rng_context: usize,
) -> c_int {
    if sig_len != SQISIGN_LVL1_SIGNATURE_BYTES || sk_len != SQISIGN_LVL1_SECRET_KEY_BYTES {
        return 0;
    }
    if sig.is_null() || sk.is_null() {
        return 0;
    }
    if msg.is_null() && msg_len != 0 {
        return 0;
    }
    let Some(callback) = fill_random else { return 0 };

    let sk_slice = unsafe { core::slice::from_raw_parts(sk, sk_len) };
    let msg_slice: &[u8] = if msg_len == 0 {
        &[]
    } else {
        unsafe { core::slice::from_raw_parts(msg, msg_len) }
    };
    let mut sk_arr = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
    sk_arr.copy_from_slice(sk_slice);
    let msg_vec = msg_slice.to_vec();

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut rng = CallbackRng {
            callback,
            context: rng_context,
        };
        sign_into_buffer(&mut rng, sig, &sk_arr, &msg_vec)
    }));
    match result {
        Ok(true) => 1,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_track_verify_crate() {
        assert_eq!(SQISIGN_LVL1_PUBLIC_KEY_BYTES, 65);
        assert_eq!(SQISIGN_LVL1_SIGNATURE_BYTES, 148);
    }

    #[test]
    fn length_mismatch_returns_zero() {
        let sig = [0u8; SQISIGN_LVL1_SIGNATURE_BYTES - 1];
        let pk = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let msg = [0u8; 1];
        let r = unsafe {
            sqisign_lvl1_verify(
                sig.as_ptr(),
                sig.len(),
                pk.as_ptr(),
                pk.len(),
                msg.as_ptr(),
                msg.len(),
            )
        };
        assert_eq!(r, 0);
    }

    #[test]
    fn null_pointer_with_zero_msg_does_not_crash() {
        // A bogus zero buffer still fails verification, but must not UB.
        let sig = [0u8; SQISIGN_LVL1_SIGNATURE_BYTES];
        let pk = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_verify(
                sig.as_ptr(),
                sig.len(),
                pk.as_ptr(),
                pk.len(),
                core::ptr::null(),
                0,
            )
        };
        assert_eq!(r, 0);
    }

    #[test]
    fn null_sig_returns_zero() {
        let pk = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_verify(
                core::ptr::null(),
                SQISIGN_LVL1_SIGNATURE_BYTES,
                pk.as_ptr(),
                pk.len(),
                core::ptr::null(),
                0,
            )
        };
        assert_eq!(r, 0);
    }

    // The callback path uses a CTR-DRBG passed through the C ABI just
    // like the entropy path uses one passed by value. Reusing the
    // same DRBG here is the simplest way to assert the two paths
    // produce identical output for the same byte stream, which is
    // what the differential KAT replay relies on.

    use core::cell::RefCell;

    struct DrbgCell(RefCell<sqisign_common::CtrDrbg>);

    unsafe extern "C" fn drbg_fill(out: *mut u8, len: usize, context: usize) {
        let cell = unsafe { &*(context as *const DrbgCell) };
        let slice = unsafe { core::slice::from_raw_parts_mut(out, len) };
        cell.0.borrow_mut().fill(slice);
    }

    #[test]
    fn keygen_with_rng_null_callback_returns_zero() {
        let mut pk = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let mut sk = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_keygen_with_rng(pk.as_mut_ptr(), pk.len(), sk.as_mut_ptr(), sk.len(), None, 0)
        };
        assert_eq!(r, 0);
    }

    #[test]
    fn callback_keygen_and_sign_roundtrip() {
        let entropy = [0x42u8; SQISIGN_LVL1_ENTROPY_BYTES];
        let cell = DrbgCell(RefCell::new(sqisign_common::CtrDrbg::new(&entropy, None)));
        let context = &cell as *const _ as usize;

        let mut pk = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let mut sk = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_keygen_with_rng(
                pk.as_mut_ptr(),
                pk.len(),
                sk.as_mut_ptr(),
                sk.len(),
                Some(drbg_fill),
                context,
            )
        };
        assert_eq!(r, 1);

        let msg = b"callback rng smoke test";
        let mut sig = [0u8; SQISIGN_LVL1_SIGNATURE_BYTES];
        let r = unsafe {
            sqisign_lvl1_sign_with_rng(
                sig.as_mut_ptr(),
                sig.len(),
                msg.as_ptr(),
                msg.len(),
                sk.as_ptr(),
                sk.len(),
                Some(drbg_fill),
                context,
            )
        };
        assert_eq!(r, 1);

        let v = unsafe {
            sqisign_lvl1_verify(
                sig.as_ptr(),
                sig.len(),
                pk.as_ptr(),
                pk.len(),
                msg.as_ptr(),
                msg.len(),
            )
        };
        assert_eq!(v, 1);
    }

    #[test]
    fn callback_path_matches_entropy_path_bit_for_bit() {
        // Same entropy fed both ways must yield the same pk/sk and
        // (given identical message and DRBG state) the same sig.
        let entropy = [0xa5u8; SQISIGN_LVL1_ENTROPY_BYTES];

        let mut pk_a = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let mut sk_a = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_keygen(
                pk_a.as_mut_ptr(),
                pk_a.len(),
                sk_a.as_mut_ptr(),
                sk_a.len(),
                entropy.as_ptr(),
                entropy.len(),
            )
        };
        assert_eq!(r, 1);

        let cell = DrbgCell(RefCell::new(sqisign_common::CtrDrbg::new(&entropy, None)));
        let context = &cell as *const _ as usize;
        let mut pk_b = [0u8; SQISIGN_LVL1_PUBLIC_KEY_BYTES];
        let mut sk_b = [0u8; SQISIGN_LVL1_SECRET_KEY_BYTES];
        let r = unsafe {
            sqisign_lvl1_keygen_with_rng(
                pk_b.as_mut_ptr(),
                pk_b.len(),
                sk_b.as_mut_ptr(),
                sk_b.len(),
                Some(drbg_fill),
                context,
            )
        };
        assert_eq!(r, 1);
        assert_eq!(pk_a, pk_b, "callback keygen pk diverged from entropy keygen pk");
        assert_eq!(sk_a, sk_b, "callback keygen sk diverged from entropy keygen sk");
    }
}
