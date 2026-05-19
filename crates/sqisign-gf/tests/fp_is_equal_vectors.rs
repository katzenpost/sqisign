//! Differential test of the ported `fp_is_equal` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation (two operands, each five little-endian 8-byte
//! limbs) and the raw returned `uint32_t` mask (`0xFFFFFFFF` for equal,
//! `0` for unequal), not a boolean: a port that returns `0x1` rather
//! than `0xFFFFFFFF` for the positive outcome would silently break the
//! downstream `fp_select` consumers (which AND the returned mask with
//! field limbs), so the test asserts *bit-equality* of the returned
//! `u32`.
//!
//! `fp_is_equal` is the second predicate boundary in the gf battery
//! (after `fp_is_zero`) and the first binary predicate. `modcmp`
//! `redc`s **both** operands first (canonicalising any redundant
//! representative of the same field element to bit-equal limb
//! vectors), then per limb applies the same `(x - 1) >> 51 & 1`
//! zero-detect trick `modis0` uses to `c[i] ^ d[i]` and AND-folds the
//! five resulting bits into `eq`. The vectors therefore must include
//! at least one pair `(a, b)` where `a != b` bit-wise but `redc(a) ==
//! redc(b)` (the non-canonical-equality case the redc-first structure
//! exists to exercise); the edge battery in
//! `tools/cdump/src/dump_main.c::emit_fp_predicate2_edges` pairs the
//! canonical zero with the radix-2^51 encoding of `p` (which reduces
//! to zero) for exactly this purpose.

use sqisign_gf::{fp_is_equal, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_is_equal.json"
);

fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be exactly {NWORDS_FIELD} u64 limbs"
    );
    let mut limbs = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        limbs[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    limbs
}

#[test]
fn fp_is_equal_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_is_equal");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut equal_records = 0usize;
    let mut unequal_records = 0usize;
    let mut non_canonical_equal_records = 0usize;

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = fp_from("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let expected_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(
            expected_bytes.len(),
            4,
            "result is a 4-byte little-endian u32 (vector {})",
            v.id
        );
        let expected = u32::from_le_bytes(expected_bytes.try_into().unwrap());

        let got = fp_is_equal(&a, &b);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp_is_equal boundary",
            v.id
        );

        match expected {
            0xFFFF_FFFF => {
                equal_records += 1;
                // A non-canonical-equality record is one whose recorded
                // a and b are bit-unequal yet the reference says they
                // represent the same field element (i.e. redc(a) ==
                // redc(b) via different non-canonical representatives).
                // The edge battery's (canonical zero, p as limbs) pair
                // is exactly such a record, in both orderings.
                if a != b {
                    non_canonical_equal_records += 1;
                }
            }
            0 => unequal_records += 1,
            other => panic!(
                "vector {}: recorded mask {other:#x} is neither 0 nor 0xFFFFFFFF",
                v.id
            ),
        }
    }

    // Count pin: equal + unequal == total, both non-zero. The 12x12
    // cross-product of the edge pattern set yields 12 diagonal equal
    // pairs plus 2 non-canonical-equal pairs (canonical zero with p
    // as limbs, in both orderings); the remaining 130 edge pairs and
    // the 1000-seed pseudo-random sweep are dominated by the unequal
    // outcome (two independent random 320-bit values have negligible
    // probability of landing on the same canonical representative).
    // Pinning the counts here guards against the edge patterns
    // silently dropping out of the emitter, and against a port that
    // returns the wrong mask width and accidentally counts on both
    // sides.
    assert_eq!(
        equal_records + unequal_records,
        file.vectors.len(),
        "every vector must be either equal or unequal (no third outcome)"
    );
    assert!(
        equal_records > 0,
        "expected at least one equal record (the diagonal of the edge cross-product)"
    );
    assert!(
        unequal_records > 0,
        "expected at least one unequal record (the off-diagonal of the edge cross-product, plus the sweep)"
    );
    assert!(
        non_canonical_equal_records > 0,
        "expected at least one non-canonical-equality record (a != b bit-wise but redc(a) == redc(b)); the edge battery's (canonical zero, p as limbs) pair is the witness"
    );
}

#[test]
fn fp_is_equal_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical inputs yields identical output. modcmp/redc
    // has no randomness, no global state, and no hidden inputs; this
    // assertion is cheap and catches any future state leak.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = fp_from("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let r1 = fp_is_equal(&a, &b);
        let r2 = fp_is_equal(&a, &b);
        assert_eq!(
            r1, r2,
            "fp_is_equal not deterministic on input from vector {}",
            v.id
        );
    }
}
