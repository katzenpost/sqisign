//! SQIsign C ABI surface (Katzenpost-facing).
//!
//! This crate exposes a minimal, verify-only C ABI over [`sqisign-verify`].
//! The immediate consumer is the Katzenpost PKI, where SQIsign is used in
//! a hybrid Ed25519+SQIsign signature on dirauth documents; the verifier
//! ships in mix nodes and clients, the signer does not. Signing is
//! deliberately deferred: it will not be wired here until the LLL / dpe
//! quaternion paths land.
//!
//! ## Safety
//!
//! `#![forbid(unsafe_code)]` would be ideal, but every FFI entry point in
//! this crate must turn caller-supplied `*const u8` pointers into Rust
//! slices, which is only possible through `unsafe`. The unsafe surface is
//! limited to that pointer-to-slice conversion (and a `catch_unwind`
//! around the body to keep panics from crossing the C boundary, which
//! would be undefined behaviour). The verification call itself is on a
//! crate that does forbid unsafe code.
//!
//! ## ABI
//!
//! The C-callable surface is intentionally tiny:
//!
//! - [`SQISIGN_LVL1_PUBLIC_KEY_BYTES`] / [`SQISIGN_LVL1_SIGNATURE_BYTES`]
//!   are the wire sizes of the level-1 public key and signature.
//! - [`sqisign_lvl1_verify`] returns `1` on a valid signature, `0` on any
//!   failure (verification false, length mismatch, null pointer with
//!   nonzero length, or panic).
//!
//! No keypair generation, no signing, no error codes beyond pass/fail.
//! That is sufficient for a hybrid verifier embedding and avoids
//! committing to a richer C API before the signing path is ready.

use core::ffi::c_int;
use core::panic::AssertUnwindSafe;
use std::panic;

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
}
