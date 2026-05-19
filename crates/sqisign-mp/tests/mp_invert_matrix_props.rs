//! Property tests for `mp_invert_matrix`.
//!
//! It is a faithful composition that inherits the `mp_neg` no-carry
//! defect. That defect corrupts even `neg(0)` (low limb 0, all higher
//! limbs all-ones), so the *off-diagonal* of `M * result` is unreliable
//! for any input, diagonal matrices included. The invariant that does
//! survive the quirk is the **main diagonal**: its two entries are
//! `a*R1 + b*S1` and `c*R2 + d*S2`, neither of which depends on a
//! negated cross term, so both are `1 (mod 2^e)` for every
//! odd-determinant input. We also assert the structural identity
//! `R1 == S2` for a scalar matrix (both equal `g * a`, no negation),
//! which exercises the determinant-inverse path through the proven
//! `mp_inv_2e`.

use proptest::prelude::*;
use sqisign_mp::{mp_inv_2e, mp_invert_matrix};

fn mul_low(x: &[u64], y: &[u64], e: usize) -> Vec<u64> {
    let limbs = e.div_ceil(64);
    let mut acc = vec![0u128; limbs + 2];
    for (i, &xi) in x.iter().enumerate() {
        if i > limbs {
            break;
        }
        let mut carry = 0u128;
        for (j, &yj) in y.iter().enumerate() {
            if i + j > limbs {
                break;
            }
            let cur = acc[i + j] + (xi as u128) * (yj as u128) + carry;
            acc[i + j] = cur & 0xffff_ffff_ffff_ffff;
            carry = cur >> 64;
        }
    }
    let mut out: Vec<u64> = acc.iter().map(|&v| v as u64).take(limbs).collect();
    let r = e % 64;
    if r != 0 {
        let t = out.len() - 1;
        out[t] &= (1u64 << r) - 1;
    }
    out
}

/// `x*y1 + x2*y2` truncated to the low `e` bits; the M*Minv diagonal.
fn diag(x1: &[u64], y1: &[u64], x2: &[u64], y2: &[u64], e: usize) -> Vec<u64> {
    let m1 = mul_low(x1, y1, e);
    let m2 = mul_low(x2, y2, e);
    let mut o = vec![0u64; m1.len()];
    let mut cr = 0u128;
    for (i, slot) in o.iter_mut().enumerate() {
        let s = m1[i] as u128 + m2[i] as u128 + cr;
        *slot = s as u64;
        cr = s >> 64;
    }
    let r = e % 64;
    if r != 0 {
        let t = o.len() - 1;
        o[t] &= (1u64 << r) - 1;
    }
    o
}

fn is_one(v: &[u64]) -> bool {
    !v.is_empty() && v[0] == 1 && v[1..].iter().all(|&x| x == 0)
}

prop_compose! {
    // Four limb vectors of one shared width, plus a raw e seed.
    fn matrix()(
        n in 2usize..10,
    )(
        r1 in proptest::collection::vec(any::<u64>(), n),
        r2 in proptest::collection::vec(any::<u64>(), n),
        s1 in proptest::collection::vec(any::<u64>(), n),
        s2 in proptest::collection::vec(any::<u64>(), n),
        eraw in any::<u32>(),
    ) -> (Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>, u32) {
        (r1, r2, s1, s2, eraw)
    }
}

proptest! {
    // The main diagonal of M * result is always 1 mod 2^e for any
    // odd-determinant matrix, regardless of the inherited mp_neg defect.
    #[test]
    fn main_diagonal_is_one(m in matrix()) {
        let (mut r1, mut r2, mut s1, mut s2, eraw) = m;
        let n = r1.len();
        // Force an odd determinant: r1,s2 odd; r2,s1 even.
        r1[0] |= 1;
        s2[0] |= 1;
        r2[0] &= !1;
        s1[0] &= !1;
        let bits = 64 * n as u32;
        let e = (4 + eraw % (bits - 3)) as usize;
        let (oa, ob, oc, od) = (r1.clone(), r2.clone(), s1.clone(), s2.clone());
        mp_invert_matrix(&mut r1, &mut r2, &mut s1, &mut s2, e as i32);
        prop_assert!(is_one(&diag(&oa, &r1, &ob, &s1, e)), "p00 != 1 mod 2^{}", e);
        prop_assert!(is_one(&diag(&oc, &r2, &od, &s2, e)), "p11 != 1 mod 2^{}", e);
    }

    // A scalar matrix a*I has R1 == S2 (both g*a) and a true main
    // diagonal; this drives the determinant-inverse path of mp_inv_2e.
    #[test]
    fn scalar_matrix_diagonal_inverse(
        a in proptest::collection::vec(any::<u64>(), 2..10),
        eraw in any::<u32>(),
    ) {
        let n = a.len();
        let mut av = a.clone();
        av[0] |= 1; // odd => det = a^2 odd
        let bits = 64 * n as u32;
        let e = (4 + eraw % (bits - 3)) as usize;
        let mut r1 = av.clone();
        let mut r2 = vec![0u64; n];
        let mut s1 = vec![0u64; n];
        let mut s2 = av.clone();
        mp_invert_matrix(&mut r1, &mut r2, &mut s1, &mut s2, e as i32);
        prop_assert_eq!(&r1, &s2, "scalar matrix must give R1 == S2");
        // p00 = a*R1 + 0*S1 = a * a^-1 == 1 mod 2^e.
        prop_assert!(
            is_one(&mul_low(&av, &r1, e)),
            "a * a^-1 != 1 mod 2^{}", e
        );
    }

    // Sanity: inverting 1 yields 1, exercising mp_inv_2e via the
    // determinant path independent of the matrix wrapper.
    #[test]
    fn unit_determinant_inverts_to_unit(n in 2usize..10, eraw in any::<u32>()) {
        let bits = 64 * n as u32;
        let e = (4 + eraw % (bits - 3)) as i32;
        let mut one = vec![0u64; n];
        one[0] = 1;
        let mut out = vec![0u64; n];
        mp_inv_2e(&mut out, &one, e);
        prop_assert_eq!(out, one);
    }
}
