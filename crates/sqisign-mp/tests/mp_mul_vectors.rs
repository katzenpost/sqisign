//! Differential test of the ported `mp_mul` against the committed
//! C-derived vectors. This includes the 18 single-limb cases where the
//! reference double-counts column 0 (`2*(a*b) mod 2^64`): the port
//! faithfully reproduces that upstream defect, so every one of the 1048
//! vectors must match bit-for-bit, the buggy ones included. The plan's
//! rule stands: the C reference is the oracle, never silently deviated
//! from. (A correction PR has been opened upstream.)

use sqisign_mp::mp_mul;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/mp_mul.json");

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

fn le_u32(b: &[u8]) -> usize {
    assert_eq!(b.len(), 4, "nwords must be a u32");
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize
}

#[test]
fn mp_mul_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_mul");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut single_limb_doublings = 0;
    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let nwords = le_u32(&decode("nwords", &v.inputs["nwords"]).expect("nwords hex"));
        let expected = limbs("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut c = vec![0u64; nwords];
        mp_mul(&mut c, &a, &b);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference (nwords={})",
            v.id, nwords
        );

        // Tally the cases that exhibit the documented upstream defect, so
        // a future upstream fix (changing the vectors) makes this test
        // notice rather than silently pass a different oracle.
        if nwords == 1 {
            let true_low = (a[0] as u128 * b[0] as u128) as u64;
            if expected[0] == true_low.wrapping_mul(2) && expected[0] != true_low {
                single_limb_doublings += 1;
            }
        }
    }
    assert_eq!(
        single_limb_doublings, 18,
        "expected exactly the 18 known single-limb doubling vectors; a \
         change here means upstream was re-pinned (regenerate + review)"
    );
}
