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
//! - [`hash::sha3_256`] / [`hash::sha3_384`] / [`hash::sha3_512`] are the
//!   fixed-output SHA3 digests (`sqisign_common::sha3_256` etc.),
//!   mirroring the reference's `sha3_256(out, in, inlen)` family.
//! - [`rng::RngSource`] is the byte-source trait every RNG-driven
//!   primitive in this workspace takes; production callers wire any
//!   implementation they trust (a future Rust port of Katzenpost's
//!   hpqc/rand, an `OsRng` shim, a hardware RNG, ...).
//! - [`rng::CtrDrbg`] is the NIST AES-256 CTR-DRBG, present **only** so
//!   the differential tests can replay the upstream KAT seeds bit-for-bit.
//!   It implements [`rng::RngSource`] so a KAT test constructs one and
//!   hands `&mut drbg` to whatever RNG-driven boundary it exercises.
//!   Production builds construct a different [`rng::RngSource`] and never
//!   touch the NIST DRBG (its backdoor history is the reason).
//! - [`mem::secure_clear`] is the optimiser-resistant memory wipe
//!   (`sqisign_common::secure_clear`), mirroring the reference's
//!   `sqisign_secure_clear`.
//!
//! With this the `common` hashing/PRNG/memory surface SQIsign actually
//! uses is ported. Seed expansion (`nistseedexpander`) is not used by the
//! Round 2 reference paths we target and is intentionally out of scope
//! unless a later phase proves otherwise.
//!
//! Correctness is established the way the whole port is: every committed
//! C-derived vector is replayed against this code (see the per-boundary
//! tests under `tests/`). Equivalence to the reference is proven, not
//! presumed.
#![forbid(unsafe_code)]

pub mod hash;
pub mod mem;
pub mod rng;

pub use hash::{
    sha3_256, sha3_384, sha3_512, shake128, shake256, Shake128Absorb, Shake128Squeeze,
    Shake256Absorb, Shake256Squeeze,
};
pub use mem::secure_clear;
pub use rng::{CtrDrbg, RngSource};
