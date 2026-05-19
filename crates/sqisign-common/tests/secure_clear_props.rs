//! Property tests for `secure_clear`.
//!
//! The observable contract is narrow and total: every byte of the slice
//! passed becomes zero, and nothing outside it is touched. These cover that
//! across arbitrary buffers and split points. (The optimiser-resistance of
//! the wipe is a property of the audited `zeroize` crate, not something a
//! behavioural test can observe; it is not asserted here.)

use proptest::prelude::*;
use sqisign_common::secure_clear;

proptest! {
    #[test]
    fn clears_every_byte(mut buf in proptest::collection::vec(any::<u8>(), 0..1024)) {
        secure_clear(&mut buf);
        prop_assert!(buf.iter().all(|&b| b == 0));
    }

    // Clearing a prefix zeroes exactly that prefix and preserves the tail,
    // mirroring the reference's `size`-bounded wipe.
    #[test]
    fn clears_only_the_given_span(
        buf in proptest::collection::vec(any::<u8>(), 1..512),
        split in 0usize..512,
    ) {
        let n = split % (buf.len() + 1);
        let mut work = buf.clone();
        secure_clear(&mut work[..n]);
        prop_assert!(work[..n].iter().all(|&b| b == 0));
        prop_assert_eq!(&work[n..], &buf[n..]);
    }
}
