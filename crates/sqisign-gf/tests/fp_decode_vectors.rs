//! Differential test of the ported `fp_decode` against the committed
//! C-derived vectors. The differential boundary is the 32-byte input
//! buffer the reference reads from `void *src`, plus the raw internal
//! `fp_t` (five little-endian 8-byte limbs) and the raw returned
//! `uint32_t` mask (`0xFFFFFFFF` on canonical in-range input, `0` on
//! out-of-range), recorded together so the per-call (input -> output +
//! result) shape is preserved.
//!
//! `fp_decode` is the modified `modimp` (see
//! `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:677..698`): fold the
//! 32 bytes into `d` in descending address order via `modshl(8) + add`
//! into limb 0, then `modfsb` returns 1 iff the decoded value is below
//! `p`, the reference negates that to a full-width mask `res`, runs
//! `nres` to convert to Montgomery form, ANDs `res` into every limb so
//! the out-of-range branch ends up zeroed. The returned `uint32_t` is
//! the low 32 bits of `res`.
//!
//! The edge battery partitions inputs into a canonical class (positive
//! result, non-zero d in general) and a non-canonical class (negative
//! result, d zeroed); the test counts both to ensure neither has
//! silently dropped out of the emitter and to pin the `& res` zeroing
//! on the negative branch.

use sqisign_gf::{fp_decode, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_decode.json"
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
fn fp_decode_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_decode");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut canonical_count = 0usize;
    let mut non_canonical_count = 0usize;

    for v in &file.vectors {
        let src_bytes = decode("src", &v.inputs["src"]).expect("src hex");
        assert_eq!(
            src_bytes.len(),
            32,
            "src is exactly 32 bytes (vector {})",
            v.id
        );
        let mut src = [0u8; 32];
        src.copy_from_slice(&src_bytes);

        let expected_d = fp_from("d", &decode("d", &v.outputs["d"]).expect("d hex"));
        let expected_result_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(
            expected_result_bytes.len(),
            4,
            "result is a 4-byte little-endian u32 (vector {})",
            v.id
        );
        let expected_result = u32::from_le_bytes(expected_result_bytes.try_into().unwrap());

        let mut d = [0u64; NWORDS_FIELD];
        let got_result = fp_decode(&mut d, &src);
        assert_eq!(
            got_result, expected_result,
            "vector {} diverged from the C reference at the fp_decode result mask",
            v.id
        );
        assert_eq!(
            d, expected_d,
            "vector {} diverged from the C reference at the fp_decode d output",
            v.id
        );

        match expected_result {
            0xFFFF_FFFF => canonical_count += 1,
            0 => {
                non_canonical_count += 1;
                // On the non-canonical branch the reference zeroes every
                // limb of d via `& res`; pin this here so a port that
                // forgets the mask doesn't sneak past.
                assert_eq!(
                    d, [0u64; NWORDS_FIELD],
                    "vector {} non-canonical: d must be zeroed on out-of-range input",
                    v.id
                );
            }
            other => panic!(
                "vector {}: recorded mask {other:#x} is neither 0 nor 0xFFFFFFFF",
                v.id
            ),
        }
    }

    assert_eq!(
        canonical_count + non_canonical_count,
        file.vectors.len(),
        "every vector must be either canonical or non-canonical (no third outcome)"
    );
    // The edge battery emits 12 canonical edges and 12 non-canonical
    // edges; the 100-record forced-non-canonical sub-sweep biases the
    // tail toward non-canonical; the 1000-record random sweep mixes
    // both classes (random 256-bit values exceed p with probability
    // roughly 1 - 5/2^8 ~ 0.98). Pinning that both classes are
    // well-represented (a non-zero count each) guards against the edge
    // patterns silently dropping out of the emitter.
    assert!(
        canonical_count >= 12,
        "expected at least the 12 canonical edges, got {canonical_count}"
    );
    assert!(
        non_canonical_count >= 12 + 100,
        "expected at least the 12 non-canonical edges plus the 100 forced-non-canonical sub-sweep, got {non_canonical_count}"
    );
}

#[test]
fn fp_decode_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical output. modshl + modfsb
    // + nres is a fixed chain with no randomness, no global state, and
    // no hidden inputs.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let src_bytes = decode("src", &v.inputs["src"]).expect("src hex");
        let mut src = [0u8; 32];
        src.copy_from_slice(&src_bytes);
        let mut d1 = [0u64; NWORDS_FIELD];
        let mut d2 = [0u64; NWORDS_FIELD];
        let r1 = fp_decode(&mut d1, &src);
        let r2 = fp_decode(&mut d2, &src);
        assert_eq!(
            (r1, d1),
            (r2, d2),
            "fp_decode not deterministic on input from vector {}",
            v.id
        );
    }
}
