//! Property tests for `fp_set_one`.
//!
//! `fp_set_one` writes the Montgomery representative of `1`,
//! `[0x19, 0, 0, 0, 0x300000000000]`, the same bit pattern the
//! reference exposes as `extern const ONE` at lines 526..530 of
//! `fp_p5248_64.c`. The reference computes it as `nres(positional 1)`;
//! the port writes the constant directly. Both are bit-equal at the
//! `fp_t` boundary, which is the contract the differential test pins.
//!
//! Two properties suffice for a constant setter:
//!  1. **Produces the Montgomery `ONE` for any pre-fill.** For arbitrary
//!     five-limb destination, `fp_set_one(out)` leaves
//!     `out == MONTGOMERY_ONE`, bit for bit.
//!  2. **Idempotent.** Calling `fp_set_one` a second time on a buffer
//!     already overwritten by `fp_set_one` leaves it bit-exact equal to
//!     a single call.

use proptest::prelude::*;
use sqisign_gf::{fp_set_one, Fp, NWORDS_FIELD};

/// The Montgomery representative of `1`.
const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];

fn uniform5() -> impl Strategy<Value = Fp> {
    (
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
    )
        .prop_map(|(a, b, c, d, e)| [a, b, c, d, e])
}

proptest! {
    #[test]
    fn produces_montgomery_one(prefill in uniform5()) {
        let mut out: Fp = prefill;
        fp_set_one(&mut out);
        prop_assert_eq!(out, MONTGOMERY_ONE);
        // Defence in depth: every limb explicitly compared, in case Fp
        // were ever widened past NWORDS_FIELD.
        for i in 0..NWORDS_FIELD {
            prop_assert_eq!(out[i], MONTGOMERY_ONE[i]);
        }
    }

    #[test]
    fn idempotent(prefill in uniform5()) {
        let mut a: Fp = prefill;
        let mut b: Fp = prefill;
        fp_set_one(&mut a);
        fp_set_one(&mut b);
        fp_set_one(&mut b); // second call: must leave b unchanged
        prop_assert_eq!(a, b);
        prop_assert_eq!(a, MONTGOMERY_ONE);
    }
}
