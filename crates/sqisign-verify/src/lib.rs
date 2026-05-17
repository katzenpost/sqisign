//! SQIsign `verification`: the top-level Verify protocol.
//!
//! Mirrors `vendor/the-sqisign/src/verification`. Ported in **Phase 2,
//! unit 9**. This crate must remain usable standalone, with `sqisign-sign`
//! and the full quaternion/id2iso paths excluded, so Katzenpost mix nodes
//! and clients pull no signing code. Not yet implemented.
#![forbid(unsafe_code)]
