//! SQIsign C ABI surface.
//!
//! Phase 4 will expose the NIST sign/verify API
//! (`crypto_sign_keypair`, `crypto_sign`, `crypto_sign_open`) as
//! `#[no_mangle] extern "C"` for the Python (pyo3/maturin) and Go (cgo)
//! bindings. `unsafe` is permitted here, and only here, because an FFI
//! boundary is inherently unsafe; the safe core crates keep
//! `#![forbid(unsafe_code)]`. Not yet implemented.
