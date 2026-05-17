//! SQIsign umbrella crate.
//!
//! Verify-only by default. The signing path is opt-in behind the `sign`
//! feature, so verify-only consumers (Katzenpost mix nodes and clients) pull
//! no signing dependencies. Not yet implemented.
#![forbid(unsafe_code)]

pub use sqisign_verify as verify;

#[cfg(feature = "sign")]
pub use sqisign_sign as sign;
