//! Differential / invariant fuzz target for the incremental SHAKE API.
//!
//! Status: harness present, runner deferred (see FUZZING.md). Ready for
//! `cargo +nightly fuzz run shake_inc` on a fuzzing host.
//!
//! The contract under fire is fragmentation-invisibility: an arbitrary
//! absorb chunking followed by an arbitrary squeeze chunking must equal the
//! one-shot result on the concatenated input, for both SHAKE256 and
//! SHAKE128. The first two bytes pick the absorb and squeeze chunk sizes;
//! the rest is the input. Linking `fips202.c` for byte-equality against the
//! C incremental path is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_common::hash::{shake128_vec, shake256_vec, Shake128Absorb, Shake256Absorb};

fuzz_target!(|data: &[u8]| {
    let (a_chunk, rest) = match data.split_first() {
        Some((n, r)) => (1usize + *n as usize, r),
        None => (1, &[][..]),
    };
    let (s_chunk, input) = match rest.split_first() {
        Some((n, r)) => (1usize + *n as usize, r),
        None => (1, &[][..]),
    };
    let out_len = 1 + (input.len() % 200);

    // SHAKE256 incremental, chunked both sides, vs one-shot.
    let mut a = Shake256Absorb::new();
    for c in input.chunks(a_chunk) {
        a.absorb(c);
    }
    let mut sq = a.finalize();
    let mut got = vec![0u8; out_len];
    let mut off = 0;
    while off < out_len {
        let n = s_chunk.min(out_len - off);
        sq.squeeze(&mut got[off..off + n]);
        off += n;
    }
    assert_eq!(got, shake256_vec(input, out_len), "shake256_inc != one-shot");

    // SHAKE128 incremental likewise.
    let mut a = Shake128Absorb::new();
    for c in input.chunks(a_chunk) {
        a.absorb(c);
    }
    let mut sq = a.finalize();
    let mut got128 = vec![0u8; out_len];
    let mut off = 0;
    while off < out_len {
        let n = s_chunk.min(out_len - off);
        sq.squeeze(&mut got128[off..off + n]);
        off += n;
    }
    assert_eq!(got128, shake128_vec(input, out_len), "shake128_inc != one-shot");
});
