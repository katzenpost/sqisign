//! Differential / invariant fuzz target for the NIST CTR-DRBG.
//!
//! Status: harness present, runner deferred (see FUZZING.md). Ready for
//! `cargo +nightly fuzz run ctr_drbg` on a fuzzing host.
//!
//! Asserts determinism and that an arbitrary draw is reproducible from the
//! same seed, including the continuation after a second call (the state
//! must evolve identically). Linking the reference's
//! `randombytes_ctrdrbg.c` for byte-equality against C is the documented
//! next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::{CtrDrbg, RngSource};

fuzz_target!(|data: &[u8]| {
    // Need at least 48 bytes of entropy; pad deterministically if short.
    let mut entropy = [0u8; 48];
    for (i, slot) in entropy.iter_mut().enumerate() {
        *slot = data.get(i).copied().unwrap_or(i as u8);
    }
    let n = 1 + (data.len() % 257);

    let mut a = CtrDrbg::new(&entropy, None);
    let mut b = CtrDrbg::new(&entropy, None);
    let mut xa = vec![0u8; n];
    let mut xb = vec![0u8; n];
    a.fill(&mut xa);
    b.fill(&mut xb);
    assert_eq!(xa, xb, "CTR-DRBG not deterministic");

    // The stream must continue identically across a further call.
    let mut ya = [0u8; 48];
    let mut yb = [0u8; 48];
    a.fill(&mut ya);
    b.fill(&mut yb);
    assert_eq!(ya, yb, "CTR-DRBG state diverged after a second draw");
});
