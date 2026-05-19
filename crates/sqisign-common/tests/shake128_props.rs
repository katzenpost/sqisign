//! Property tests for `sqisign_common::shake128`.
//!
//! The C-derived vectors prove equivalence at the points the reference
//! happened to traverse. These cover the space around them with invariants
//! the reference *relies on*. They deliberately do not re-derive expected
//! output via another SHAKE library: that would be tautological (our
//! implementation is sha3-backed). They assert structure instead. The same
//! invariants hold for SHAKE256; only the rate differs.

use proptest::prelude::*;
use sqisign_common::hash::{shake128, shake128_vec};

proptest! {
    // Same input, same output. A hash that is not a function is no hash.
    #[test]
    fn deterministic(input in proptest::collection::vec(any::<u8>(), 0..512),
                      out_len in 0usize..512) {
        prop_assert_eq!(
            shake128_vec(&input, out_len),
            shake128_vec(&input, out_len)
        );
    }

    // The defining property of an extendable-output function, and the one
    // the reference depends on when it squeezes varying lengths from the
    // same absorb: a shorter output is a prefix of a longer one.
    #[test]
    fn xof_prefix_stable(input in proptest::collection::vec(any::<u8>(), 0..256),
                         short in 0usize..128, extra in 0usize..256) {
        let long = shake128_vec(&input, short + extra);
        let pre = shake128_vec(&input, short);
        prop_assert_eq!(&long[..short], &pre[..]);
    }

    // The buffer-filling and allocating forms must never disagree, and the
    // output length is exactly what was asked for.
    #[test]
    fn forms_agree_and_length_exact(
        input in proptest::collection::vec(any::<u8>(), 0..256),
        out_len in 0usize..300,
    ) {
        let v = shake128_vec(&input, out_len);
        prop_assert_eq!(v.len(), out_len);
        let mut b = vec![0u8; out_len];
        shake128(&input, &mut b);
        prop_assert_eq!(v, b);
    }

    // Sensitivity: flipping one input bit should change the digest. Not a
    // security claim, just a smoke test that input is actually absorbed.
    #[test]
    fn input_is_absorbed(mut input in proptest::collection::vec(any::<u8>(), 1..128),
                          bit in 0usize..1024) {
        let before = shake128_vec(&input, 32);
        let idx = (bit / 8) % input.len();
        input[idx] ^= 1u8 << (bit % 8);
        prop_assert_ne!(before, shake128_vec(&input, 32));
    }

    // SHAKE128 and SHAKE256 are distinct functions: for the same input and
    // length their outputs must differ (different domain/rate). This guards
    // against an accidental mis-wiring of one to the other.
    #[test]
    fn distinct_from_shake256(input in proptest::collection::vec(any::<u8>(), 0..256)) {
        use sqisign_common::hash::shake256_vec;
        prop_assert_ne!(shake128_vec(&input, 32), shake256_vec(&input, 32));
    }
}
