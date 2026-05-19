//! Property tests for `mp_copy`: it is the identity into `b`, total and
//! length-preserving.

use proptest::prelude::*;
use sqisign_mp::mp_copy;

proptest! {
    #[test]
    fn copies_exactly(a in proptest::collection::vec(any::<u64>(), 0..256)) {
        let mut b = vec![0xa5a5_a5a5_a5a5_a5a5u64; a.len()];
        mp_copy(&mut b, &a);
        prop_assert_eq!(b, a);
    }

    // Destination prior contents are fully overwritten (no residue).
    #[test]
    fn overwrites_destination(a in proptest::collection::vec(any::<u64>(), 1..128),
                              fill in any::<u64>()) {
        let mut b = vec![fill; a.len()];
        mp_copy(&mut b, &a);
        prop_assert_eq!(&b, &a);
    }
}
