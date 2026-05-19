//! Invariant fuzz target for `secure_clear`.
//!
//! Status: harness present, runner deferred (see FUZZING.md). Ready for
//! `cargo +nightly fuzz run secure_clear` on a fuzzing host.
//!
//! Asserts the total contract: clearing a chosen prefix zeroes exactly that
//! prefix and leaves the tail untouched, for an arbitrary buffer and split.
//! Linking the reference's `mem.c` for byte-equality against C is the
//! documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::secure_clear;

fuzz_target!(|data: &[u8]| {
    let (split_byte, buf_src) = match data.split_first() {
        Some((n, rest)) => (*n as usize, rest),
        None => (0, &[][..]),
    };
    let original = buf_src.to_vec();
    let n = if original.is_empty() {
        0
    } else {
        split_byte % (original.len() + 1)
    };

    let mut work = original.clone();
    secure_clear(&mut work[..n]);

    assert!(work[..n].iter().all(|&b| b == 0), "prefix not cleared");
    assert_eq!(&work[n..], &original[n..], "tail was modified");
});
