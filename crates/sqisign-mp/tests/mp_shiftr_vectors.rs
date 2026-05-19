//! Differential test of the ported `mp_shiftr` against the committed
//! C-derived vectors. Both outputs are checked: the shifted array `r` and
//! the returned `bit_out` (the value's bit 0 before shifting).

use sqisign_mp::mp_shiftr;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_shiftr.json"
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

fn le_u64(b: &[u8]) -> u64 {
    assert_eq!(b.len(), 8, "u64 field");
    u64::from_le_bytes(b.try_into().unwrap())
}

#[test]
fn mp_shiftr_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_shiftr");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let mut x = limbs("x", &decode("x", &v.inputs["x"]).expect("x hex"));
        let shift = le_u32(&decode("shift", &v.inputs["shift"]).expect("shift hex"));
        let expected_r = limbs("r", &decode("r", &v.outputs["r"]).expect("r hex"));
        let expected_bit = le_u64(&decode("bit_out", &v.outputs["bit_out"]).expect("bit_out hex"));

        let bit = mp_shiftr(&mut x, shift);
        assert_eq!(
            x, expected_r,
            "vector {} array diverged (shift={})",
            v.id, shift
        );
        assert_eq!(
            bit, expected_bit,
            "vector {} returned bit diverged (shift={})",
            v.id, shift
        );
    }
}
