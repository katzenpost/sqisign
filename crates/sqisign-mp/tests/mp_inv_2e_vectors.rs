//! Differential test of the ported `mp_inv_2e` against the committed
//! C-derived vectors. Beyond bit-equality, it independently confirms the
//! defining relation `a * b == 1 (mod 2^e)` for every vector via a
//! self-contained limb-wise multiply (no bignum dependency), so a future
//! upstream re-pin that broke the inverse would be caught.

use sqisign_mp::mp_inv_2e;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_inv_2e.json"
);

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

fn le_i32(b: &[u8]) -> i32 {
    assert_eq!(b.len(), 4, "e is an i32");
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

/// Low `e` bits of `a * b`, as a limb vector wide enough to hold them.
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
    // Mask the top partial limb to e bits.
    let r = e % 64;
    if r != 0 {
        let top = out.len() - 1;
        out[top] &= (1u64 << r) - 1;
    }
    out
}

#[test]
fn mp_inv_2e_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_inv_2e");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let e = le_i32(&decode("e", &v.inputs["e"]).expect("e hex"));
        let expected = limbs("b", &decode("b", &v.outputs["b"]).expect("b hex"));

        let mut b = vec![0u64; a.len()];
        mp_inv_2e(&mut b, &a, e);
        assert_eq!(
            b, expected,
            "vector {} diverged from the C reference (e={})",
            v.id, e
        );

        // a * b == 1 mod 2^e.
        let prod = mul_low_ebits(&a, &b, e as usize);
        let mut one = vec![0u64; prod.len()];
        one[0] = 1;
        assert_eq!(
            prod, one,
            "vector {}: a*b != 1 mod 2^{} (upstream re-pin?)",
            v.id, e
        );
    }
}
