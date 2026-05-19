//! Differential test of the ported `mp_mul2` against the committed
//! C-derived vectors. `mp_mul2` is a *partial* product: the reference
//! omits the `a1*b0` cross term, computing `a*b - (a1*b0)*2^64`. The port
//! reproduces the reference's structure, so all 2296 vectors must match
//! bit-for-bit; the test also asserts the documented partial-product
//! identity holds for every vector, so a future upstream re-pin that
//! changed the behaviour would be noticed.

use sqisign_mp::mp_mul2;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/mp_mul2.json");

fn limbs(label: &str, bytes: &[u8]) -> Vec<u64> {
    assert_eq!(
        bytes.len() % 8,
        0,
        "{label} not a whole number of u64 limbs"
    );
    bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

/// 256-bit value of up to four little-endian limbs.
fn val(w: &[u64]) -> u128_pair {
    // a*b for two-limb inputs fits 256 bits; track as (lo128, hi128).
    let mut lo = 0u128;
    let mut hi = 0u128;
    for (i, &limb) in w.iter().enumerate() {
        if i < 2 {
            lo |= (limb as u128) << (64 * i);
        } else {
            hi |= (limb as u128) << (64 * (i - 2));
        }
    }
    u128_pair { lo, hi }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy)]
struct u128_pair {
    lo: u128,
    hi: u128,
}

#[test]
fn mp_mul2_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_mul2");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let expected = limbs("c", &decode("c", &v.outputs["c"]).expect("c hex"));
        assert_eq!(a.len(), 2);
        assert_eq!(b.len(), 2);
        assert_eq!(expected.len(), 4);

        let mut c = [0u64; 4];
        mp_mul2(&mut c, &a, &b);
        assert_eq!(
            c.to_vec(),
            expected,
            "vector {} diverged from the C reference",
            v.id
        );

        // Independently confirm the documented partial-product identity:
        // C == a*b - (a1*b0)*2^64, as a 256-bit check.
        let prod = mul256(&a, &b);
        let adjusted = sub_at(prod, 1, (a[1] as u128) * (b[0] as u128));
        let got = val(&c);
        assert!(
            got.lo == adjusted.0 && got.hi == adjusted.1,
            "vector {}: partial-product identity violated (upstream re-pin?)",
            v.id
        );
    }
}

fn mul256(a: &[u64], b: &[u64]) -> (u128, u128) {
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
        let vv = r[k] + carry;
        *slot = vv as u64;
        carry = vv >> 64;
    }
    (
        (out[0] as u128) | ((out[1] as u128) << 64),
        (out[2] as u128) | ((out[3] as u128) << 64),
    )
}

fn sub_at(v: (u128, u128), limb: usize, amount: u128) -> (u128, u128) {
    // Represent (lo,hi) as four 64-bit limbs, subtract `amount` (<=128
    // bits) starting at `limb`, repack.
    let mut w = [
        v.0 as u64,
        (v.0 >> 64) as u64,
        v.1 as u64,
        (v.1 >> 64) as u64,
    ];
    let mut amt = amount;
    let mut k = limb;
    while amt != 0 && k < 4 {
        let cur = w[k] as u128;
        let s = amt & 0xffff_ffff_ffff_ffff;
        if cur >= s {
            w[k] = (cur - s) as u64;
            amt >>= 64;
        } else {
            w[k] = (cur + (1u128 << 64) - s) as u64;
            amt = (amt >> 64) + 1;
        }
        k += 1;
    }
    (
        (w[0] as u128) | ((w[1] as u128) << 64),
        (w[2] as u128) | ((w[3] as u128) << 64),
    )
}
