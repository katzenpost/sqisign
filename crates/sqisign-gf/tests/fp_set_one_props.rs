//! Property test for `fp_set_one`: the value it writes is the
//! multiplicative identity in the Montgomery domain, whatever the
//! destination held before.

use proptest::prelude::*;
use sqisign_gf::{fp_is_equal, fp_mul, fp_set_one, Fp};

const MASK51: u64 = (1 << 51) - 1;

/// Canonical-range limbs: 51 bits each, top limb below `p4`.
fn element() -> impl Strategy<Value = Fp> {
    any::<[u64; 5]>().prop_map(|l| {
        [
            l[0] & MASK51,
            l[1] & MASK51,
            l[2] & MASK51,
            l[3] & MASK51,
            l[4] & ((1 << 44) - 1),
        ]
    })
}

proptest! {
    #[test]
    fn is_multiplicative_identity(a in element(), prefill in any::<[u64; 5]>()) {
        let mut one: Fp = prefill;
        fp_set_one(&mut one);
        let mut out: Fp = [0; 5];
        fp_mul(&mut out, &a, &one);
        prop_assert_eq!(fp_is_equal(&out, &a), 0xFFFF_FFFF);
    }
}
