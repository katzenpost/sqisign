//! Differential test of the ported `fp_set_zero` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! `fp_set_zero` is the reference's `modzer`: a plain five-limb
//! zero-fill with no `prop`, no `2p` correction and no reduction. The
//! recorded output is therefore the canonical all-zero representative
//! on every record regardless of the destination's prior contents. The
//! record's "input" is the destination pre-fill the setter is asked to
//! overwrite: feeding the ported function the same pre-fill and
//! comparing the resulting limbs against the C-recorded output is the
//! only way to catch a no-op or partial-write port. The test asserts
//! both: bit-equality to the recorded reference output, *and* that
//! every recorded output is the all-zero limb vector, the latter pinned
//! as a `count == total` assertion so the zero-output property cannot
//! quietly drop out of the emitter or the port.

use sqisign_gf::{fp_set_zero, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_set_zero.json"
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
fn fp_set_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_set_zero");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    const ZERO: Fp = [0u64; NWORDS_FIELD];

    let total = file.vectors.len();
    let mut zero_hits = 0usize;

    for v in &file.vectors {
        let prefill = fp_from(
            "prefill",
            &decode("prefill", &v.inputs["prefill"]).expect("prefill hex"),
        );
        let expected = fp_from("out", &decode("out", &v.outputs["out"]).expect("out hex"));

        // Pre-fill the destination exactly as the C harness did; the
        // setter must overwrite it. A port that wrote fewer than five
        // limbs would leave the corresponding pre-fill bytes visible
        // and diverge from the recorded all-zero output.
        let mut x: Fp = prefill;
        fp_set_zero(&mut x);
        assert_eq!(
            x, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
        if expected == ZERO {
            zero_hits += 1;
        }
    }

    // Count pin: every record's output is the canonical all-zero
    // representative. modzer is bit-exact zero by construction; this
    // guards against the emitter silently emitting an unrelated op or
    // the port gaining unintended residual limbs.
    assert_eq!(
        zero_hits, total,
        "fp_set_zero must produce the bit-exact canonical zero on every vector"
    );
}
