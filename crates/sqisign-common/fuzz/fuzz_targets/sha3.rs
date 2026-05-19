//! Differential / invariant fuzz target for the fixed-output SHA3 digests.
//!
//! Status: harness present, runner deferred (see FUZZING.md). Ready for
//! `cargo +nightly fuzz run sha3` on a fuzzing host.
//!
//! Asserts determinism, the intrinsic digest lengths, and mutual
//! distinctness of the three widths on the same input. Linking `fips202.c`
//! for byte-equality against the C path is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::{sha3_256, sha3_384, sha3_512};

fuzz_target!(|data: &[u8]| {
    let a = sha3_256(data);
    let b = sha3_384(data);
    let c = sha3_512(data);

    assert_eq!(a, sha3_256(data), "sha3_256 not deterministic");
    assert_eq!(b, sha3_384(data), "sha3_384 not deterministic");
    assert_eq!(c, sha3_512(data), "sha3_512 not deterministic");

    assert_eq!(a.len(), 32);
    assert_eq!(b.len(), 48);
    assert_eq!(c.len(), 64);

    assert_ne!(&a[..], &b[..32], "sha3_256/384 collided");
    assert_ne!(&a[..], &c[..32], "sha3_256/512 collided");
    assert_ne!(&b[..], &c[..48], "sha3_384/512 collided");
});
