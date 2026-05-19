//! Differential test of the ported `swap_ct` against the committed
//! C-derived vectors. In-place on a and b; both mutated arrays are
//! recorded. The battery includes arbitrary options (not just 0 /
//! all-ones), pinning the per-bit conditional swap; independently
//! checked against `(a&!opt)|(b&opt)` and its mirror.

use sqisign_mp::swap_ct;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/swap_ct.json");

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

fn le_u64(b: &[u8]) -> u64 {
    assert_eq!(b.len(), 8, "option is a u64");
    u64::from_le_bytes(b.try_into().unwrap())
}

#[test]
fn swap_ct_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::swap_ct");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a0 = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b0 = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let option = le_u64(&decode("option", &v.inputs["option"]).expect("option hex"));
        let ea = limbs(
            "a_out",
            &decode("a_out", &v.outputs["a_out"]).expect("a_out hex"),
        );
        let eb = limbs(
            "b_out",
            &decode("b_out", &v.outputs["b_out"]).expect("b_out hex"),
        );

        let mut a = a0.clone();
        let mut b = b0.clone();
        swap_ct(&mut a, &mut b, option);
        assert_eq!(a, ea, "vector {}: a_out diverged", v.id);
        assert_eq!(b, eb, "vector {}: b_out diverged", v.id);

        let exp_a: Vec<u64> = (0..a0.len())
            .map(|i| (a0[i] & !option) | (b0[i] & option))
            .collect();
        let exp_b: Vec<u64> = (0..a0.len())
            .map(|i| (b0[i] & !option) | (a0[i] & option))
            .collect();
        assert_eq!(a, exp_a, "vector {}: not the per-bit swap (a)", v.id);
        assert_eq!(b, exp_b, "vector {}: not the per-bit swap (b)", v.id);
    }
}
