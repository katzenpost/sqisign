//! SQIsign `signature`: KeyGen and Sign protocols.
//!
//! Mirrors `the-sqisign/src/signature/ref/lvlx`. Pulls the full
//! quaternion / id2iso paths and the dim-two theta-chain randomised
//! variant. Production callers thread their own [`sqisign_common::RngSource`]
//! through [`protocols_keygen`] and [`protocols_sign`]; the KAT round-trip
//! tests use [`sqisign_common::CtrDrbg`] seeded by the recorded entropy.
#![forbid(unsafe_code)]
#![allow(non_snake_case)]
#![allow(clippy::needless_range_loop)]

pub mod encode;
pub mod keygen;
pub mod sign;

pub use encode::{secret_key_from_bytes, secret_key_to_bytes, SECRETKEY_BYTES};
pub use keygen::{protocols_keygen, SecretKey};
pub use sign::protocols_sign;
