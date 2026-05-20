//! Differential test of the ported `fp_inv` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! `fp_inv` is the *in-place* one-liner `modinv(*x, NULL, *x)` (see
//! `the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:628..632`). The
//! reference passes the same buffer as both input `x` and destination
//! `z`; the cdump harness snapshots the input into a separate `a` and
//! the in-place result into `c`, so the recorded record carries both
//! sides of the in-place mutation. The port wrapper snapshots `*x` into
//! a local to resolve the borrow conflict, then writes back; the
//! resulting `x` must be bit-equal to the reference's in-place
//! destination.
//!
//! `modinv` builds `x^-1 = x^(p-2)` via Fermat: progenitor `x^((p-3)/4)`,
//! two further squarings, multiply by `x`. On the field zero the chain
//! squares-and-multiplies down to the canonical zero (no division by
//! zero); the differential battery includes that edge as id 0.

use sqisign_gf::{fp_inv, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp_inv.json");

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
fn fp_inv_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_inv");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        // The port's wrapper is in-place; the recorded boundary is
        // (input -> output), so we feed the snapshotted input through
        // the in-place primitive and read the resulting buffer.
        let mut x = a;
        fp_inv(&mut x);
        assert_eq!(
            x, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
    }
}

#[test]
fn fp_inv_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical limbs. modinv is a
    // fixed addition chain over modpro/modnsqr/modmul with no
    // randomness, no global state, and no hidden inputs.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let mut x1 = a;
        let mut x2 = a;
        fp_inv(&mut x1);
        fp_inv(&mut x2);
        assert_eq!(
            x1, x2,
            "fp_inv not deterministic on input from vector {}",
            v.id
        );
    }
}
