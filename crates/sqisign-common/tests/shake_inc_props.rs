//! Property tests for the incremental SHAKE API.
//!
//! The defining contract of the incremental boundary is that fragmentation
//! is invisible: absorbing in any chunking and squeezing in any chunking
//! must equal the one-shot result on the concatenated input. These cover
//! that space around the C-derived vectors. The one-shot functions are the
//! oracle here, and they are themselves vector-proven against the reference.

use proptest::prelude::*;
use sqisign_common::hash::{shake128_vec, shake256_vec, Shake128Absorb, Shake256Absorb};

/// Turn a flat input and a list of chunk sizes into successive slices,
/// the last chunk taking whatever remains so the cover is exact.
fn chunked(input: &[u8], sizes: &[usize]) -> Vec<(usize, usize)> {
    let mut spans = Vec::new();
    let mut off = 0;
    for (i, &s) in sizes.iter().enumerate() {
        let n = if i + 1 == sizes.len() {
            input.len() - off
        } else {
            s.min(input.len() - off)
        };
        spans.push((off, n));
        off += n;
    }
    if off < input.len() {
        spans.push((off, input.len() - off));
    }
    spans
}

proptest! {
    #[test]
    fn shake256_inc_equals_one_shot(
        input in proptest::collection::vec(any::<u8>(), 0..600),
        a_sizes in proptest::collection::vec(0usize..200, 1..8),
        s_sizes in proptest::collection::vec(1usize..200, 1..8),
        out_len in 1usize..400,
    ) {
        let mut a = Shake256Absorb::new();
        for (o, n) in chunked(&input, &a_sizes) {
            a.absorb(&input[o..o + n]);
        }
        let mut sq = a.finalize();
        let mut got = vec![0u8; out_len];
        let spans = chunked(&got.clone(), &s_sizes);
        for (o, n) in spans {
            sq.squeeze(&mut got[o..o + n]);
        }
        prop_assert_eq!(got, shake256_vec(&input, out_len));
    }

    #[test]
    fn shake128_inc_equals_one_shot(
        input in proptest::collection::vec(any::<u8>(), 0..600),
        a_sizes in proptest::collection::vec(0usize..200, 1..8),
        s_sizes in proptest::collection::vec(1usize..200, 1..8),
        out_len in 1usize..400,
    ) {
        let mut a = Shake128Absorb::new();
        for (o, n) in chunked(&input, &a_sizes) {
            a.absorb(&input[o..o + n]);
        }
        let mut sq = a.finalize();
        let mut got = vec![0u8; out_len];
        let spans = chunked(&got.clone(), &s_sizes);
        for (o, n) in spans {
            sq.squeeze(&mut got[o..o + n]);
        }
        prop_assert_eq!(got, shake128_vec(&input, out_len));
    }

    // The two incremental sponges must not coincide (rate mis-wiring guard).
    #[test]
    fn inc_128_distinct_from_256(input in proptest::collection::vec(any::<u8>(), 0..256)) {
        let mut a = Shake256Absorb::new();
        a.absorb(&input);
        let mut s = a.finalize();
        let mut o256 = [0u8; 32];
        s.squeeze(&mut o256);

        let mut a = Shake128Absorb::new();
        a.absorb(&input);
        let mut s = a.finalize();
        let mut o128 = [0u8; 32];
        s.squeeze(&mut o128);

        prop_assert_ne!(o256, o128);
    }
}
