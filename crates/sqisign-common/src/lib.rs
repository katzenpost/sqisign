//! SQIsign `common`: hashing (SHA-3 / SHAKE), seed expansion, PRNG, memory
//! utilities.
//!
//! Mirrors `vendor/the-sqisign/src/common`. Ported in **Phase 1, unit 1**.
//!
//! Nothing is ported yet. Phase 0 only proves the differential harness end
//! to end; see `tests/shake256_vectors.rs`, which compares committed
//! reference vectors against a known-good SHAKE oracle (the `sha3` crate),
//! not against any code in this crate.
#![forbid(unsafe_code)]
