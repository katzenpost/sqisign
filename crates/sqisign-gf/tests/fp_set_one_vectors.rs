//! Differential test of the ported `fp_set_one` against the committed
//! C-derived vectors.
//!
//! The reference's `fp_set_one` wraps `modone`, which writes positional
//! `1` then calls `nres(a, a)` to convert it to its Montgomery
//! representative. The on-the-wire output is therefore not
//! `[1, 0, 0, 0, 0]` but the Montgomery `ONE`,
//! `[0x19, 0, 0, 0, 0x300000000000]`, the same value the reference
//! exposes as `extern const ONE` at lines 526..530 of `fp_p5248_64.c`.
//! The differential boundary is the raw five-limb `fp_t`
//! representation: the port must reproduce the recorded reference
//! output bit-for-bit, and the destination's prior contents must not
//! leak into the result.
//!
//! Two assertions per vector: bit-equality to the recorded reference
//! output AND a `count == total` pin that every recorded output equals
//! the Montgomery `ONE` constant.

use sqisign_gf::{fp_set_one, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_set_one.json"
);

/// Montgomery representative of `1`; must match the reference's
/// `extern const ONE` at `fp_p5248_64.c:526..530`.
const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];

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
fn fp_set_one_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_set_one");
    let total = file.vectors.len();
    assert!(
        total >= 1000,
        "expected the full battery, found {total} vectors"
    );

    let mut unit_hits = 0usize;
    for v in &file.vectors {
        let prefill = fp_from(
            "prefill",
            &decode("prefill", &v.inputs["prefill"]).expect("hex"),
        );
        let expected = fp_from("out", &decode("out", &v.outputs["out"]).expect("hex"));

        // Pre-fill so a partial-write port diverges visibly from the
        // recorded reference output.
        let mut out: Fp = prefill;
        fp_set_one(&mut out);
        assert_eq!(
            out, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );

        if expected == MONTGOMERY_ONE {
            unit_hits += 1;
        }
    }

    assert_eq!(
        unit_hits, total,
        "fp_set_one must produce the bit-exact Montgomery ONE on every vector"
    );
}
