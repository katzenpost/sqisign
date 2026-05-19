//! Differential test of the ported `fp_is_zero` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation (five little-endian 8-byte limbs) and the raw
//! returned `uint32_t` mask (`0xFFFFFFFF` for zero, `0` for nonzero), not
//! a boolean: a port that returns `0x1` rather than `0xFFFFFFFF` for the
//! positive outcome would silently break the downstream `fp_select`
//! consumers (which AND the returned mask with field limbs), so the test
//! asserts *bit-equality* of the returned `u32`.
//!
//! `fp_is_zero` is the first predicate boundary in the gf battery.
//! `modis0` first `redc`s its argument (canonicalising any redundant
//! representative of `0 mod p` to the bit-exact `[0, 0, 0, 0, 0]`), then
//! OR-folds the canonical limbs and returns `0`/`1`; the `-(uint32_t)`
//! wrapper turns that into `0`/`0xFFFFFFFF`. The vectors therefore must
//! include at least one non-canonical encoding of the field zero (a
//! representative of `p`, `2p`, etc.) to exercise the Montgomery
//! reduction path, not just the canonical `[0, 0, 0, 0, 0]`; see the
//! edge battery in `tools/cdump/src/dump_main.c::emit_fp_predicate_edges`.

use sqisign_gf::{fp_is_zero, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_is_zero.json"
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
fn fp_is_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_is_zero");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut zero_count = 0usize;
    let mut nonzero_count = 0usize;

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(
            expected_bytes.len(),
            4,
            "result is a 4-byte little-endian u32 (vector {})",
            v.id
        );
        let expected = u32::from_le_bytes(expected_bytes.try_into().unwrap());

        let got = fp_is_zero(&a);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp_is_zero boundary",
            v.id
        );

        match expected {
            0xFFFF_FFFF => zero_count += 1,
            0 => nonzero_count += 1,
            other => panic!(
                "vector {}: recorded mask {other:#x} is neither 0 nor 0xFFFFFFFF",
                v.id
            ),
        }
    }

    // Pin the count of zero-detected and nonzero records. The edge
    // battery emits two zero-representative inputs (the canonical
    // `[0, 0, 0, 0, 0]` and the non-canonical `[mask51, mask51, mask51,
    // mask51, p4 - 1]`, the radix-2^51 encoding of `p` itself) and ten
    // nonzero edges; the 1000-seed pseudo-random sweep is all nonzero
    // (random 320-bit values have negligible probability of landing on
    // a representative of `0 mod p`). Pinning the counts here guards
    // against the edge patterns silently dropping out of the emitter,
    // and against a port that returns the wrong mask width and
    // accidentally counts on both sides.
    assert_eq!(
        zero_count, 2,
        "expected exactly two zero-detected vectors (canonical zero plus the non-canonical `p` encoding)"
    );
    assert_eq!(
        nonzero_count, 1010,
        "expected exactly 1010 nonzero vectors (ten edge nonzeros plus the 1000-seed sweep)"
    );
    assert_eq!(
        zero_count + nonzero_count,
        file.vectors.len(),
        "every vector must be either zero-detected or nonzero (no third outcome)"
    );
}

#[test]
fn fp_is_zero_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical output. modis0/redc has
    // no randomness, no global state, and no hidden inputs; this
    // assertion is cheap and catches any future state leak.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let r1 = fp_is_zero(&a);
        let r2 = fp_is_zero(&a);
        assert_eq!(
            r1, r2,
            "fp_is_zero not deterministic on input from vector {}",
            v.id
        );
    }
}
