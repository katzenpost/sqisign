//! Property tests for the fixed-output SHA3 digests.
//!
//! The C-derived vectors prove equivalence at the points the reference
//! traversed. These cover the space with structural invariants: the digest
//! is a deterministic function of its input, has the intrinsic length, is
//! sensitive to every input bit, and the three widths are mutually distinct
//! (a guard against accidentally wiring one to another).

use proptest::prelude::*;
use sqisign_common::{sha3_256, sha3_384, sha3_512};

proptest! {
    #[test]
    fn deterministic(input in proptest::collection::vec(any::<u8>(), 0..1024)) {
        prop_assert_eq!(sha3_256(&input), sha3_256(&input));
        prop_assert_eq!(sha3_384(&input), sha3_384(&input));
        prop_assert_eq!(sha3_512(&input), sha3_512(&input));
    }

    #[test]
    fn lengths_are_intrinsic(input in proptest::collection::vec(any::<u8>(), 0..256)) {
        prop_assert_eq!(sha3_256(&input).len(), 32);
        prop_assert_eq!(sha3_384(&input).len(), 48);
        prop_assert_eq!(sha3_512(&input).len(), 64);
    }

    #[test]
    fn input_is_absorbed(mut input in proptest::collection::vec(any::<u8>(), 1..128),
                          bit in 0usize..1024) {
        let before = sha3_512(&input);
        let idx = (bit / 8) % input.len();
        input[idx] ^= 1u8 << (bit % 8);
        prop_assert_ne!(before, sha3_512(&input));
    }

    // The three widths share a Keccak permutation but differ in capacity and
    // padding; their common prefix must not coincide for the same input.
    #[test]
    fn widths_are_distinct(input in proptest::collection::vec(any::<u8>(), 0..256)) {
        let a = sha3_256(&input);
        let b = sha3_384(&input);
        let c = sha3_512(&input);
        prop_assert_ne!(&a[..], &b[..32]);
        prop_assert_ne!(&a[..], &c[..32]);
        prop_assert_ne!(&b[..], &c[..48]);
    }
}
