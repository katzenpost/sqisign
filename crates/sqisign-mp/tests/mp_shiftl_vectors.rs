//! Differential test of the ported `mp_shiftl` against the committed
//! C-derived vectors. The recorded `r` is the input shifted left and
//! truncated to `nwords`; the harness only exercises the reference's
//! defined shift domain (1..=63).

use sqisign_mp::mp_shiftl;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_shiftl.json"
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
fn mp_shiftl_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_shiftl");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let mut x = limbs("x", &decode("x", &v.inputs["x"]).expect("x hex"));
        let shift = le_u32(&decode("shift", &v.inputs["shift"]).expect("shift hex"));
        let nwords = le_u32(&decode("nwords", &v.inputs["nwords"]).expect("nwords hex")) as usize;
        let expected = limbs("r", &decode("r", &v.outputs["r"]).expect("r hex"));

        assert_eq!(x.len(), nwords, "vector {}: x vs nwords", v.id);

        mp_shiftl(&mut x, shift);
        assert_eq!(
            x, expected,
            "vector {} diverged from the C reference (nwords={}, shift={})",
            v.id, nwords, shift
        );
    }
}
