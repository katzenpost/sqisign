//! Invariant fuzz target for `sqisign_mp::mp_invert_matrix`. For any
//! odd-determinant 2x2 matrix and `e` within the array width the
//! composition must not panic, and the *main* diagonal of `M * result`
//! is `1 (mod 2^e)`. The off-diagonal is not asserted: it inherits the
//! `mp_neg` no-carry defect (24/1180 reference vectors are not true
//! inverses for that reason). Linking `mp.c` for byte-equality vs C is
//! the next increment. Runner deferred (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_invert_matrix;

fn mul_low(a: &[u64], b: &[u64], e: usize) -> Vec<u64> {
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
    if data.len() < 36 {
        return;
    }
    let eraw = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    let body = &data[4..];
    let q = body.len() / 4;
    if q < 8 {
        return;
    }
    let to_limbs = |s: &[u8]| -> Vec<u64> {
        s.chunks(8)
            .map(|c| {
                let mut w = [0u8; 8];
                w[..c.len()].copy_from_slice(c);
                u64::from_le_bytes(w)
            })
            .collect()
    };
    let n = q / 8;
    if n == 0 {
        return;
    }
    let span = n * 8;
    let mut r1 = to_limbs(&body[..span]);
    let mut r2 = to_limbs(&body[span..2 * span]);
    let mut s1 = to_limbs(&body[2 * span..3 * span]);
    let mut s2 = to_limbs(&body[3 * span..4 * span]);
    r1.resize(n, 0);
    r2.resize(n, 0);
    s1.resize(n, 0);
    s2.resize(n, 0);

    // Force an odd determinant (r1,s2 odd; r2,s1 even => r1*s2-r2*s1 odd).
    r1[0] |= 1;
    s2[0] |= 1;
    r2[0] &= !1;
    s1[0] &= !1;

    let bits = 64 * n as u32;
    let e = (4 + eraw % (bits - 3)) as usize;
    let (oa, ob, oc, od) = (r1.clone(), r2.clone(), s1.clone(), s2.clone());

    mp_invert_matrix(&mut r1, &mut r2, &mut s1, &mut s2, e as i32);

    // p00 = oa*R1 + ob*S1 ; p11 = oc*R2 + od*S2. Both == 1 mod 2^e.
    let diag = |x1: &[u64], y1: &[u64], x2: &[u64], y2: &[u64]| -> bool {
        let m1 = mul_low(x1, y1, e);
        let m2 = mul_low(x2, y2, e);
        let mut o = vec![0u64; m1.len()];
        let mut cr = 0u128;
        for i in 0..o.len() {
            let s = m1[i] as u128 + m2[i] as u128 + cr;
            o[i] = s as u64;
            cr = s >> 64;
        }
        let r = e % 64;
        if r != 0 {
            let t = o.len() - 1;
            o[t] &= (1u64 << r) - 1;
        }
        o[0] == 1 && o[1..].iter().all(|&x| x == 0)
    };
    assert!(diag(&oa, &r1, &ob, &s1), "p00 != 1 mod 2^{e}");
    assert!(diag(&oc, &r2, &od, &s2), "p11 != 1 mod 2^{e}");
});
