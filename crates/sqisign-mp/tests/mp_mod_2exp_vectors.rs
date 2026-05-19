//! Differential test of the ported `mp_mod_2exp` against the committed
//! C-derived vectors. The battery includes `e` straddling word
//! boundaries and `e` at/over the full width (the no-op path).

use sqisign_mp::mp_mod_2exp;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_mod_2exp.json"
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

fn le_u32(b: &[u8]) -> u32 {
    assert_eq!(b.len(), 4, "u32 field");
    u32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

#[test]
fn mp_mod_2exp_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_mod_2exp");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let mut a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let e = le_u32(&decode("e", &v.inputs["e"]).expect("e hex"));
        let expected = limbs("r", &decode("r", &v.outputs["r"]).expect("r hex"));

        mp_mod_2exp(&mut a, e);
        assert_eq!(
            a, expected,
            "vector {} diverged from the C reference (e={})",
            v.id, e
        );
    }
}
