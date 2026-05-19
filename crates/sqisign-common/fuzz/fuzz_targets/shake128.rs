//! Differential / invariant fuzz target for `sqisign_common::shake128`.
//!
//! Status: harness present, runner deferred, exactly as for the `shake256`
//! target (see FUZZING.md). The target builds and is ready for
//! `cargo +nightly fuzz run shake128` on a fuzzing host.
//!
//! The plan's end goal is a true differential against the linked C
//! reference. As a first step it asserts the structural invariants the
//! reference relies on (determinism, XOF prefix stability) plus the
//! distinctness of SHAKE128 from SHAKE256, which catches an accidental
//! mis-wiring of one rate to the other without needing the C toolchain
//! linked in. Linking `fips202.c` here for byte-equality against C is the
//! documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::hash::{shake128_vec, shake256_vec};

fuzz_target!(|data: &[u8]| {
    // First byte chooses an output length in 0..=255; the rest is input.
    let (out_len, input) = match data.split_first() {
        Some((n, rest)) => (*n as usize, rest),
        None => (0, &[][..]),
    };

    let a = shake128_vec(input, out_len);
    let b = shake128_vec(input, out_len);
    assert_eq!(a, b, "shake128 is not deterministic");

    // XOF prefix stability: the shorter squeeze is a prefix of a longer one.
    let longer = shake128_vec(input, out_len + 17);
    assert_eq!(&longer[..out_len], &a[..], "XOF prefix not stable");

    // SHAKE128 and SHAKE256 must never coincide on a non-empty squeeze.
    if out_len > 0 {
        assert_ne!(
            a,
            shake256_vec(input, out_len),
            "shake128 collided with shake256 (rate mis-wiring?)"
        );
    }
});
