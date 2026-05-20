//! Differential test of the ported `fp_is_square` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation (five little-endian 8-byte limbs) and the raw
//! returned `uint32_t` mask (`0xFFFFFFFF` for a quadratic residue or the
//! field zero, `0` for a non-residue), not a boolean: a port that returns
//! `0x1` rather than `0xFFFFFFFF` for the positive outcome would silently
//! break the downstream `fp_select` consumers (which AND the returned
//! mask with field limbs), so the test asserts *bit-equality* of the
//! returned `u32`.
//!
//! `fp_is_square` is the unary predicate `-(uint32_t)modqr(NULL, *a)`.
//! `modqr` evaluates the Euler criterion via the progenitor: square the
//! progenitor, multiply by `x`, test for unity via `modis1`. The
//! `modis0(x)` OR makes the field zero positive by convention. The
//! battery includes the canonical zero, the radix-2^51 encoding of `p`
//! (also zero), `MONTGOMERY_ONE` (positive), and arbitrary full-width
//! limbs (roughly half negative).

use sqisign_gf::{fp_is_square, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_is_square.json"
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
fn fp_is_square_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_is_square");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut square_count = 0usize;
    let mut nonsquare_count = 0usize;

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

        let got = fp_is_square(&a);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp_is_square boundary",
            v.id
        );

        match expected {
            0xFFFF_FFFF => square_count += 1,
            0 => nonsquare_count += 1,
            other => panic!(
                "vector {}: recorded mask {other:#x} is neither 0 nor 0xFFFFFFFF",
                v.id
            ),
        }
    }

    // Every vector must be either a square or a non-residue (no third
    // outcome); both classes are well-represented so the test isn't
    // accidentally exercising only one branch.
    assert_eq!(
        square_count + nonsquare_count,
        file.vectors.len(),
        "every vector must be either a square or a non-residue (no third outcome)"
    );
    assert!(
        square_count > 0,
        "expected at least one quadratic-residue record (edges include zero and the Montgomery one)"
    );
    assert!(
        nonsquare_count > 0,
        "expected at least one non-residue record (the random sweep covers both branches)"
    );
}

#[test]
fn fp_is_square_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical output. modqr is a
    // fixed addition chain over modpro/modsqr/modmul/modis1/modis0 with
    // no randomness, no global state, and no hidden inputs.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let r1 = fp_is_square(&a);
        let r2 = fp_is_square(&a);
        assert_eq!(
            r1, r2,
            "fp_is_square not deterministic on input from vector {}",
            v.id
        );
    }
}
