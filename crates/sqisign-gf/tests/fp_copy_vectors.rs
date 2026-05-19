//! Differential test of the ported `fp_copy` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! Unlike `fp_add`/`fp_sub`/`fp_neg`, `fp_copy` is the reference's
//! `modcpy`, a plain five-limb assignment with no `prop`, no `2p`
//! correction and no reduction. The recorded output is therefore exactly
//! the recorded input limb for limb, including for the non-canonical
//! patterns in the edge battery. The test asserts both: bit-equality to
//! the recorded reference output, *and* identity (output == input) on
//! every record, the latter pinned as a `count == total` assertion
//! (precedent `mp_copy_vectors.rs`) so the identity property cannot
//! quietly drop out of the emitter or the port.

use sqisign_gf::{fp_copy, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp_copy.json");

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
fn fp_copy_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_copy");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let total = file.vectors.len();
    let mut identity_hits = 0usize;

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        // A non-trivial pre-fill so a no-op or partial-write port would
        // diverge visibly, mirroring mp_copy_vectors.rs's 0xdead pre-fill.
        let mut c: Fp = [0xdead_beef_dead_beefu64; NWORDS_FIELD];
        fp_copy(&mut c, &a);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
        if c == a {
            identity_hits += 1;
        }
    }

    // Count pin: every record's output is the recorded input limb for
    // limb. modcpy is bit-exact identity by construction; this guards
    // against the emitter silently emitting an unrelated op or the port
    // gaining an unintended reduction.
    assert_eq!(
        identity_hits, total,
        "fp_copy must be the bit-exact identity on every vector"
    );
}
