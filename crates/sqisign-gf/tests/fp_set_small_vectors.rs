//! Differential test of the ported `fp_set_small` against the committed
//! C-derived vectors.
//!
//! The reference's `fp_set_small` is the thin wrapper
//! `modint((int)val, *x)`: it narrows the `digit_t` argument to `int`
//! (`int32_t` on every reasonable platform), writes the sign-extended
//! result at limb 0 with limbs 1..=4 zeroed, then calls `nres(a, a)`
//! to convert the positional residue to its Montgomery representative.
//! The on-the-wire output is therefore the Montgomery representative
//! of `(val as i32 as i64) mod p`, in the redundant radix-2^51
//! representation.
//!
//! The differential boundary is the raw five-limb `fp_t` representation
//! plus the recorded `val` (an 8-byte little-endian `u64`). The port
//! must reproduce the recorded reference output bit-for-bit, including
//! on records whose `val` exceeds `int32` range: those are exactly the
//! ones that pin the high-bits-ignored narrowing.
//!
//! Three assertions:
//!  1. Per-vector bit-equality at the `fp_t` boundary.
//!  2. `count > 0` for records exercising the high-bits-ignored
//!     narrowing (records whose `val` differs from
//!     `val as i32 as u64`). The C dump's seeded sweep alone produces
//!     ~1000 such records; the edge battery adds explicit witnesses.
//!  3. `val == 1` records produce the Montgomery `ONE` constant
//!     bit-for-bit, the independent algorithmic check that `nres`
//!     applied through the boundary to the positional `1` produces the
//!     same constant `fp_set_one` writes directly.

use sqisign_gf::{fp_set_small, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_set_small.json"
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

fn u64_from(label: &str, bytes: &[u8]) -> u64 {
    assert_eq!(bytes.len(), 8, "{label} must be exactly 8 bytes (one u64)");
    u64::from_le_bytes(bytes.try_into().unwrap())
}

#[test]
fn fp_set_small_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_set_small");
    let total = file.vectors.len();
    assert!(
        total >= 1000,
        "expected the full battery, found {total} vectors"
    );

    let mut high_bits_hits = 0usize;
    let mut one_hits = 0usize;

    for v in &file.vectors {
        let prefill = fp_from(
            "prefill",
            &decode("prefill", &v.inputs["prefill"]).expect("prefill hex"),
        );
        let val = u64_from("val", &decode("val", &v.inputs["val"]).expect("val hex"));
        let expected = fp_from("out", &decode("out", &v.outputs["out"]).expect("out hex"));

        // Pre-fill the destination exactly as the C harness did; the
        // setter must overwrite all five limbs. A port that wrote fewer
        // than five limbs would leave the corresponding pre-fill bytes
        // visible and diverge from the recorded output.
        let mut out: Fp = prefill;
        fp_set_small(&mut out, val);
        assert_eq!(
            out, expected,
            "vector {} (val=0x{val:016x}) diverged from the C reference at the fp_t boundary",
            v.id
        );

        // (2) Pin: at least one record exercises the high-bits-ignored
        // narrowing (val differs from `val as i32 as u64`).
        let narrowed_then_sx = (val as i32) as i64 as u64;
        if val != narrowed_then_sx {
            high_bits_hits += 1;
        }

        // (3) Pin: val == 1 records produce the Montgomery ONE
        // bit-pattern, the algorithmic confirmation that NRES_C is
        // correct through the boundary.
        if val == 1 {
            assert_eq!(
                expected, MONTGOMERY_ONE,
                "fp_set_small(_, 1) must produce the Montgomery ONE constant at vector {}",
                v.id
            );
            one_hits += 1;
        }
    }

    assert!(
        high_bits_hits > 0,
        "expected at least one record exercising the high-bits-ignored narrowing \
         (val != val as i32 as u64); found {high_bits_hits}"
    );
    assert!(
        one_hits > 0,
        "expected at least one val == 1 record to anchor the nres(1) == MONTGOMERY_ONE \
         oracle through the boundary; found {one_hits}"
    );
}
