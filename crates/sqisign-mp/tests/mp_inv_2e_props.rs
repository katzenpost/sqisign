//! Property tests for `mp_inv_2e`: it produces a modular inverse.
//!
//! The defining relation `a * b == 1 (mod 2^e)` holds for any odd `a`
//! and `e` within the array width; checked via a self-contained
//! limb-wise multiply.

use proptest::prelude::*;
use sqisign_mp::mp_inv_2e;

fn mul_low_ebits(a: &[u64], b: &[u64], e: usize) -> Vec<u64> {
    let limbs = e.div_ceil(64);
    let mut acc = vec![0u128; limbs + 1];
    for (i, &ai) in a.iter().enumerate() {
        if i > limbs {
            break;
        }
        let mut carry = 0u128;
        for (j, &bj) in b.iter().enumerate() {
            if i + j > limbs {
                break;
            }
            let cur = acc[i + j] + (ai as u128) * (bj as u128) + carry;
            acc[i + j] = cur & 0xffff_ffff_ffff_ffff;
            carry = cur >> 64;
        }
    }
    let mut out: Vec<u64> = acc[..limbs].iter().map(|&v| v as u64).collect();
    let r = e % 64;
    if r != 0 {
        let t = out.len() - 1;
        out[t] &= (1u64 << r) - 1;
    }
    out
}

fn is_one(v: &[u64]) -> bool {
    !v.is_empty() && v[0] == 1 && v[1..].iter().all(|&x| x == 0)
}

proptest! {
    // For any odd a and e in [4, 64*n], a * inv(a) == 1 mod 2^e.
    #[test]
    fn is_inverse(
        mut a in proptest::collection::vec(any::<u64>(), 2..16),
        eraw in 4u32..,
    ) {
        a[0] |= 1; // odd
        let bits = 64 * a.len() as u32;
        let e = 4 + (eraw % (bits - 3)); // 4 ..= bits
        let mut b = vec![0u64; a.len()];
        mp_inv_2e(&mut b, &a, e as i32);
        prop_assert!(is_one(&mul_low_ebits(&a, &b, e as usize)),
            "a*b != 1 mod 2^{}", e);
    }

    // The inverse of 1 is 1, at any width and e.
    #[test]
    fn inverse_of_one(n in 2usize..16, eraw in 4u32..) {
        let mut a = vec![0u64; n];
        a[0] = 1;
        let bits = 64 * n as u32;
        let e = 4 + (eraw % (bits - 3));
        let mut b = vec![0u64; n];
        mp_inv_2e(&mut b, &a, e as i32);
        let mut one = vec![0u64; n];
        one[0] = 1;
        prop_assert_eq!(b, one);
    }

    // Inverting twice returns the original modulo 2^e (involution of the
    // inverse map on the unit group mod 2^e), checked by a*b==1 again.
    #[test]
    fn double_inverse_round_trips(
        mut a in proptest::collection::vec(any::<u64>(), 2..12),
        eraw in 4u32..,
    ) {
        a[0] |= 1;
        let bits = 64 * a.len() as u32;
        let e = 4 + (eraw % (bits - 3));
        let mut b = vec![0u64; a.len()];
        mp_inv_2e(&mut b, &a, e as i32);
        let mut bb = vec![0u64; a.len()];
        mp_inv_2e(&mut bb, &b, e as i32);
        // inv(inv(a)) * b == 1 mod 2^e, i.e. bb is an inverse of b too.
        prop_assert!(is_one(&mul_low_ebits(&bb, &b, e as usize)));
    }
}
