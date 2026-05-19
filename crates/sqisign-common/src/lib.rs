//! SQIsign `common`: hashing (SHA-3 / SHAKE), seed expansion, PRNG, memory
//! utilities.
//!
//! Mirrors `vendor/the-sqisign/src/common`. **Phase 1, unit 1.**
//!
//! Ported so far:
//! - [`hash::shake256`] is the one-shot SHAKE256 XOF boundary
//!   (`sqisign_common::shake256`), corresponding to the reference's
//!   `shake256(out, outlen, in, inlen)` in `common/generic/fips202.c`.
//! - [`hash::shake128`] is the one-shot SHAKE128 XOF boundary
//!   (`sqisign_common::shake128`), the reference's
//!   `shake128(out, outlen, in, inlen)` in the same translation unit.
//! - [`hash::Shake256Absorb`] / [`hash::Shake128Absorb`] are the
//!   incremental absorb/finalize/squeeze API
//!   (`sqisign_common::shake256_inc`, `..::shake128_inc`), mirroring the
//!   reference's `shake*_inc_init/_absorb/_finalize/_squeeze`.
//!
//! Not yet ported (later units within Phase 1 `common`): SHA3-256/384/512,
//! seed expansion, the CTR-DRBG (`randombytes`), and the memory-zeroing
//! utilities.
//!
//! Correctness is established the way the whole port is: every committed
//! C-derived vector is replayed against this code (see the per-boundary
//! tests under `tests/`). Equivalence to the reference is proven, not
//! presumed.
#![forbid(unsafe_code)]

pub mod hash;

pub use hash::{
    shake128, shake256, Shake128Absorb, Shake128Squeeze, Shake256Absorb, Shake256Squeeze,
};
