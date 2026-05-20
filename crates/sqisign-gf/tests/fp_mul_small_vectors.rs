//! Differential test of the ported `fp_mul_small` against the committed
//! C-derived vectors.
//!
//! The reference's `fp_mul_small` is the thin wrapper
//! `modmli(*a, (int)val, *x)`: narrow the `uint32_t val` to `int`
//! (`int32_t` on every reasonable platform), build the Montgomery
//! representative of the narrowed-and-sign-extended integer via
//! `modint(b, t)` into a five-limb scratch, then run a single
//! `modmul(a, t, x)` cross product. The on-the-wire output is therefore
//! the Montgomery representative of `a * (val as i32 as i64) mod p`, in
//! the redundant radix-2^51 representation.
//!
//! The differential boundary is the raw five-limb `fp_t` representation
//! for both `a` and `out`, plus the recorded `val` (an 8-byte
//! little-endian `u64`, with only the low 32 bits observable through the
//! boundary by construction). The port must reproduce the recorded
//! reference output bit-for-bit, including on records whose `val`
//! exceeds `int32` range: those are exactly the ones that pin the
//! high-bits-ignored narrowing.
//!
//! Two assertions:
//!  1. Per-vector bit-equality at the `fp_t` boundary.
//!  2. `count > 0` for records exercising the high-bits-ignored
//!     narrowing (records whose `val` differs from
//!     `val as i32 as u64`). The C dump's seeded sweep alone produces
//!     ~1000 such records; the edge battery adds explicit witnesses.

use sqisign_gf::{fp_mul_small, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_mul_small.json"
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

fn u64_from(label: &str, bytes: &[u8]) -> u64 {
    assert_eq!(bytes.len(), 8, "{label} must be exactly 8 bytes (one u64)");
    u64::from_le_bytes(bytes.try_into().unwrap())
}

#[test]
fn fp_mul_small_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_mul_small");
    let total = file.vectors.len();
    assert!(
        total >= 1000,
        "expected the full battery, found {total} vectors"
    );

    let mut high_bits_hits = 0usize;

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let val = u64_from("val", &decode("val", &v.inputs["val"]).expect("val hex"));
        let expected = fp_from("out", &decode("out", &v.outputs["out"]).expect("out hex"));

        // The destination is overwritten in full by modmul's per-column
        // writes, so the starting buffer is irrelevant; a port that
        // forgot to write any limb would diverge against the recorded
        // output regardless of the pre-fill.
        let mut out: Fp = [0u64; NWORDS_FIELD];
        fp_mul_small(&mut out, &a, val as u32);
        assert_eq!(
            out, expected,
            "vector {} (val=0x{val:016x}) diverged from the C reference at the fp_t boundary",
            v.id
        );

        // (2) Pin: at least one record exercises the high-bits-ignored
        // narrowing (val differs from `val as i32 as u64`). The C
        // wrapper narrows `uint32_t val` to `(int)val` before calling
        // modmli, so any val with bit 31 set has the same observable
        // behaviour as the `i32 as u64` sign-extended image.
        let narrowed_then_sx = (val as i32) as i64 as u64;
        if val != narrowed_then_sx {
            high_bits_hits += 1;
        }
    }

    assert!(
        high_bits_hits > 0,
        "expected at least one record exercising the high-bits-ignored narrowing \
         (val != val as i32 as u64); found {high_bits_hits}"
    );
}
