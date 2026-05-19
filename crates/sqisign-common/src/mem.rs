//! Optimiser-resistant memory clearing.
//!
//! Mirrors the reference's `sqisign_secure_clear(mem, size)` in
//! `src/common/generic/mem.c`, which zeroes `size` bytes through a
//! `volatile` function pointer to `memset` so a dead-store optimisation
//! cannot elide the wipe. The reference also exposes `sqisign_secure_free`
//! (the same clear followed by `free`); in Rust "free" is ownership and
//! `Drop`, so the boundary worth porting and proving is the clear itself.
//!
//! The optimiser-resistant write is delegated to the audited `zeroize`
//! crate rather than hand-rolled, for the same reason the SHAKE and AES
//! primitives are: this crate is `#![forbid(unsafe_code)]`, and a bespoke
//! `write_volatile` loop would both require `unsafe` and forgo audited
//! code. The observable contract (exactly `buf.len()` bytes become zero,
//! nothing beyond is touched) is proven bit-equal to the reference by the
//! committed C-derived vectors (`tests/secure_clear_vectors.rs`).

use zeroize::Zeroize;

/// Zero every byte of `buf` in a way the optimiser may not elide,
/// equivalent to `sqisign_secure_clear(buf.as_mut_ptr(), buf.len())`.
///
/// Only `buf` is affected; the caller controls the length by the slice it
/// passes, exactly as the reference's `size` argument selects how much of
/// an allocation to wipe.
pub fn secure_clear(buf: &mut [u8]) {
    buf.zeroize();
}
