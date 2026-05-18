//! Differential / invariant fuzz target for `sqisign_common::shake256`.
//!
//! Status: harness present, runner deferred. `cargo-fuzz` + libFuzzer are
//! not installed in the current environment (as with valgrind for the
//! iai-callgrind gate); this is recorded as a deferred dependency in
//! FUZZING.md, not silently dropped. The target builds and is ready for
//! `cargo +nightly fuzz run shake256` on a fuzzing host.
//!
//! The plan's end goal for this target is a true differential against the
//! linked C reference. As a first step it asserts the structural
//! invariants the reference relies on (determinism, XOF prefix stability),
//! which already catches any divergence in those properties without
//! needing the C toolchain linked into the fuzz binary. Linking
//! `vendor/the-sqisign/src/common/generic/fips202.c` here for byte-equality
//! against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::hash::shake256_vec;

fuzz_target!(|data: &[u8]| {
    // First byte chooses an output length in 0..=255; the rest is input.
    let (out_len, input) = match data.split_first() {
        Some((n, rest)) => (*n as usize, rest),
        None => (0, &[][..]),
    };

    let a = shake256_vec(input, out_len);
    let b = shake256_vec(input, out_len);
    assert_eq!(a, b, "shake256 is not deterministic");

    // XOF prefix stability: the shorter squeeze is a prefix of a longer one.
    let longer = shake256_vec(input, out_len + 17);
    assert_eq!(&longer[..out_len], &a[..], "XOF prefix not stable");
});
