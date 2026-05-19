//! Differential test of the ported `mp_invert_matrix` against the
//! committed C-derived vectors.
//!
//! Bit-equality to the recorded reference output is the oracle. The test
//! also classifies every vector by whether `M * result == I (mod 2^e)`
//! and asserts the count that are *not* true inverses is exactly 24 --
//! these are the inherited `mp_neg` no-carry defect propagating through
//! the two negations. Pinning the count means a future upstream re-pin
//! (e.g. once the mp_neg fix lands) is noticed rather than silently
//! re-baselined.

use sqisign_mp::mp_invert_matrix;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_invert_matrix.json"
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

/// `x * y` truncated to the low `e` bits, as a limb vector.
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
        if i + y.len() <= limbs {
            acc[i + y.len()] += carry;
        }
    }
    mask(
        &acc.iter()
            .map(|&v| v as u64)
            .take(limbs)
            .collect::<Vec<_>>(),
        e,
    )
}

fn add_low(x: &[u64], y: &[u64], e: usize) -> Vec<u64> {
    let n = x.len().max(y.len());
    let mut out = vec![0u64; n];
    let mut carry = 0u128;
    for (i, slot) in out.iter_mut().enumerate() {
        let xi = *x.get(i).unwrap_or(&0) as u128;
        let yi = *y.get(i).unwrap_or(&0) as u128;
        let s = xi + yi + carry;
        *slot = s as u64;
        carry = s >> 64;
    }
    mask(&out, e)
}

fn mask(v: &[u64], e: usize) -> Vec<u64> {
    let limbs = e.div_ceil(64);
    let mut out: Vec<u64> = v
        .iter()
        .copied()
        .chain(std::iter::repeat(0))
        .take(limbs)
        .collect();
    let r = e % 64;
    if r != 0 && !out.is_empty() {
        let t = out.len() - 1;
        out[t] &= (1u64 << r) - 1;
    }
    out
}

fn is_zero(v: &[u64]) -> bool {
    v.iter().all(|&x| x == 0)
}
fn is_one(v: &[u64]) -> bool {
    !v.is_empty() && v[0] == 1 && v[1..].iter().all(|&x| x == 0)
}

#[test]
fn mp_invert_matrix_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_invert_matrix");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut not_true_inverse = 0usize;
    for v in &file.vectors {
        let a = limbs("r1", &decode("r1", &v.inputs["r1"]).expect("r1"));
        let b = limbs("r2", &decode("r2", &v.inputs["r2"]).expect("r2"));
        let c = limbs("s1", &decode("s1", &v.inputs["s1"]).expect("s1"));
        let d = limbs("s2", &decode("s2", &v.inputs["s2"]).expect("s2"));
        let e = le_i32(&decode("e", &v.inputs["e"]).expect("e")) as usize;
        let er1 = limbs(
            "r1_out",
            &decode("r1_out", &v.outputs["r1_out"]).expect("r1o"),
        );
        let er2 = limbs(
            "r2_out",
            &decode("r2_out", &v.outputs["r2_out"]).expect("r2o"),
        );
        let es1 = limbs(
            "s1_out",
            &decode("s1_out", &v.outputs["s1_out"]).expect("s1o"),
        );
        let es2 = limbs(
            "s2_out",
            &decode("s2_out", &v.outputs["s2_out"]).expect("s2o"),
        );

        let (mut r1, mut r2, mut s1, mut s2) = (a.clone(), b.clone(), c.clone(), d.clone());
        mp_invert_matrix(&mut r1, &mut r2, &mut s1, &mut s2, e as i32);
        assert_eq!(r1, er1, "vector {} r1_out diverged", v.id);
        assert_eq!(r2, er2, "vector {} r2_out diverged", v.id);
        assert_eq!(s1, es1, "vector {} s1_out diverged", v.id);
        assert_eq!(s2, es2, "vector {} s2_out diverged", v.id);

        // M = ((a,b),(c,d)); Minv = ((r1,r2),(s1,s2)). Classify.
        let p00 = add_low(&mul_low(&a, &r1, e), &mul_low(&b, &s1, e), e);
        let p01 = add_low(&mul_low(&a, &r2, e), &mul_low(&b, &s2, e), e);
        let p10 = add_low(&mul_low(&c, &r1, e), &mul_low(&d, &s1, e), e);
        let p11 = add_low(&mul_low(&c, &r2, e), &mul_low(&d, &s2, e), e);
        let true_inverse = is_one(&p00) && is_zero(&p01) && is_zero(&p10) && is_one(&p11);
        if !true_inverse {
            not_true_inverse += 1;
        }
    }
    assert_eq!(
        not_true_inverse, 24,
        "expected exactly the 24 inherited-mp_neg-defect vectors; a \
         change here means upstream was re-pinned (regenerate + review)"
    );
}
