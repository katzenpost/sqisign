//! Property tests for `mp_mul2`.
//!
//! `mp_mul2` is the reference's partial two-digit product: it omits the
//! `a1*b0` cross term, so `c == a*b - (a1*b0)*2^64`. That defect-or-by-
//! design identity is the invariant; it collapses to the true product
//! exactly when `a1 == 0` or `b0 == 0`.

use proptest::prelude::*;
use sqisign_mp::mp_mul2;

fn mul256(a: &[u64; 2], b: &[u64; 2]) -> [u64; 4] {
    let mut r = [0u128; 5];
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            let p = (ai as u128) * (bj as u128);
            r[i + j] += p & 0xffff_ffff_ffff_ffff;
            r[i + j + 1] += p >> 64;
        }
    }
    let mut out = [0u64; 4];
    let mut carry = 0u128;
    for (k, slot) in out.iter_mut().enumerate() {
        let v = r[k] + carry;
        *slot = v as u64;
        carry = v >> 64;
    }
    out
}

fn sub_at(mut w: [u64; 4], limb: usize, mut amount: u128) -> [u64; 4] {
    let mut k = limb;
    while amount != 0 && k < 4 {
        let cur = w[k] as u128;
        let s = amount & 0xffff_ffff_ffff_ffff;
        if cur >= s {
            w[k] = (cur - s) as u64;
            amount >>= 64;
        } else {
            w[k] = (cur + (1u128 << 64) - s) as u64;
            amount = (amount >> 64) + 1;
        }
        k += 1;
    }
    w
}

proptest! {
    // The defining identity, over the full 256-bit result.
    #[test]
    fn equals_product_minus_a1b0(a0 in any::<u64>(), a1 in any::<u64>(),
                                 b0 in any::<u64>(), b1 in any::<u64>()) {
        let a = [a0, a1];
        let b = [b0, b1];
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &a, &b);
        let expect = sub_at(mul256(&a, &b), 1, (a1 as u128) * (b0 as u128));
        prop_assert_eq!(c, expect);
    }

    // When a1 == 0 the omitted term is zero: the true full product.
    #[test]
    fn full_product_when_a1_zero(a0 in any::<u64>(), b0 in any::<u64>(), b1 in any::<u64>()) {
        let a = [a0, 0];
        let b = [b0, b1];
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &a, &b);
        prop_assert_eq!(c, mul256(&a, &b));
    }

    // When b0 == 0 likewise.
    #[test]
    fn full_product_when_b0_zero(a0 in any::<u64>(), a1 in any::<u64>(), b1 in any::<u64>()) {
        let a = [a0, a1];
        let b = [0, b1];
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &a, &b);
        prop_assert_eq!(c, mul256(&a, &b));
    }

    // Multiplying by zero is zero (both terms vanish).
    #[test]
    fn times_zero_is_zero(a0 in any::<u64>(), a1 in any::<u64>()) {
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &[a0, a1], &[0, 0]);
        prop_assert_eq!(c, [0u64; 4]);
    }
}
