//! SQIsign `common`: hashing (SHA-3 / SHAKE), seed expansion, PRNG, memory
//! utilities.
//!
//! Mirrors `vendor/the-sqisign/src/common`. **Phase 1, unit 1.**
//!
//! Ported so far:
//! - [`hash::shake256`] — the one-shot SHAKE256 XOF boundary
//!   (`sqisign_common::shake256`), corresponding to the reference's
//!   `shake256(out, outlen, in, inlen)` in `common/generic/fips202.c`.
//!
//! Not yet ported (later units within Phase 1 `common`): SHAKE128, the
//! incremental absorb/squeeze API, SHA3-256/384/512, seed expansion, the
//! CTR-DRBG (`randombytes`), and the memory-zeroing utilities.
//!
//! Correctness is established the way the whole port is: every committed
//! C-derived vector is replayed against this code (see
//! `tests/shake256_vectors.rs`). Equivalence to the reference is proven,
//! not presumed.
#![forbid(unsafe_code)]

pub mod hash;

pub use hash::shake256;
