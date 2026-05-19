//! Differential test of the ported `mp_neg` against the committed
//! C-derived vectors. The reference adds the two's-complement `+1` to
//! limb 0 only, with no carry propagation, so it equals `-a` only when
//! `a[0] != 0`. The port reproduces that; all 1042 vectors must match
//! bit-for-bit, and the test independently re-derives the faithful model
//! and the `== -a iff a[0] != 0` characterisation so a future upstream
//! re-pin changing the behaviour is noticed.

use sqisign_mp::mp_neg;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/mp_neg.json");

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

#[test]
fn mp_neg_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_neg");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut quirk = 0;
    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = limbs("r", &decode("r", &v.outputs["r"]).expect("r hex"));

        let mut got = a.clone();
        mp_neg(&mut got);
        assert_eq!(got, expected, "vector {} diverged from the reference", v.id);

        // Faithful model: complement then +1 on limb 0 with no carry.
        let mut model: Vec<u64> = a.iter().map(|x| !x).collect();
        model[0] = model[0].wrapping_add(1);
        assert_eq!(got, model, "vector {}: faithful model mismatch", v.id);

        // Characterisation: equals true (-a) mod 2^(64n) iff a[0] != 0.
        let is_true_neg = is_negation(&a, &got);
        if a[0] != 0 {
            assert!(is_true_neg, "vector {}: a[0]!=0 must be true -a", v.id);
        } else if !is_true_neg {
            quirk += 1;
        }
    }
    assert!(
        quirk > 0,
        "expected some a[0]==0 quirk vectors; battery changed (upstream re-pin?)"
    );
}

/// True iff `r == (-a) mod 2^(64*len)`, i.e. `a + r == 0` with carry.
fn is_negation(a: &[u64], r: &[u64]) -> bool {
    let mut carry = 0u128;
    for i in 0..a.len() {
        let s = a[i] as u128 + r[i] as u128 + carry;
        if (s as u64) != 0 {
            return false;
        }
        carry = s >> 64;
    }
    true
}
