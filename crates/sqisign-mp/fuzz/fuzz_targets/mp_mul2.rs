//! Invariant fuzz target for `sqisign_mp::mp_mul2`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the reproduced partial-product identity
//! `c == a*b - (a1*b0)*2^64` for arbitrary two-digit operands, and that
//! it collapses to the full product when a1 or b0 is zero. Linking
//! `mp.c` for byte-equality vs C is the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
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

fuzz_target!(|data: &[u8]| {
    if data.len() < 32 {
        return;
    }
    let rd = |o: usize| {
        let mut w = [0u8; 8];
        w.copy_from_slice(&data[o..o + 8]);
        u64::from_le_bytes(w)
    };
    let a = [rd(0), rd(8)];
    let b = [rd(16), rd(24)];

    let mut c = [0u64; 4];
    mp_mul2(&mut c, &a, &b);

    let expect = sub_at(mul256(&a, &b), 1, (a[1] as u128) * (b[0] as u128));
    assert_eq!(c, expect, "mp_mul2 != a*b - a1*b0*2^64");
});
