//! Invariant fuzz target for `sqisign_mp::mp_inv_2e`: for any odd `a`
//! and `e` within the array width, `a * b == 1 (mod 2^e)`. Linking
//! `mp.c` for byte-equality vs C is the next increment. Runner deferred
//! (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
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

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 {
        return;
    }
    let mut a: Vec<u64> = data[4..]
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    if a.len() < 2 {
        return;
    }
    a[0] |= 1; // odd, as the reference requires
    let bits = 64 * a.len() as u32;
    let eraw = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    let e = 4 + (eraw % (bits - 3)); // 4 ..= bits

    let mut b = vec![0u64; a.len()];
    mp_inv_2e(&mut b, &a, e as i32);

    let prod = mul_low_ebits(&a, &b, e as usize);
    assert!(
        prod[0] == 1 && prod[1..].iter().all(|&x| x == 0),
        "a*b != 1 mod 2^{e}"
    );
});
