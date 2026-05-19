//! Invariant fuzz target for `sqisign_gf::fp_sqr`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_sqr` on a fuzzing host.
//!
//! Asserts the two sound raw-limb invariants for arbitrary (possibly
//! non-canonical) five-limb inputs: bit-exact equivalence to
//! `fp_mul(a, a)` (the per-column set of partial products is the same
//! multiset under both orderings and `u128` accumulation is associative
//! and commutative, so the running `t` is bit-equal at every masking
//! point), and the structural carry-propagation invariant that limbs
//! 0..=3 are below 2^51 (limb 4 is left fully unmasked by design and is
//! not constrained: the reference's final `c[4] = (spint)t` is a 64-bit
//! truncation, no `& mask`). Linking the reference's `fp_p5248_64.c` for
//! byte-equality against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_mul, fp_sqr, Fp, NWORDS_FIELD};

fn fp(bytes: &[u8]) -> Fp {
    let mut n = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks(8).take(NWORDS_FIELD).enumerate() {
        let mut w = [0u8; 8];
        w[..chunk.len()].copy_from_slice(chunk);
        n[i] = u64::from_le_bytes(w);
    }
    n
}

fuzz_target!(|data: &[u8]| {
    if data.len() < NWORDS_FIELD * 8 {
        return;
    }
    let a = fp(&data[..NWORDS_FIELD * 8]);

    let mut sqr = [0u64; NWORDS_FIELD];
    let mut mul = [0u64; NWORDS_FIELD];
    fp_sqr(&mut sqr, &a);
    fp_mul(&mut mul, &a, &a);
    assert_eq!(sqr, mul, "fp_sqr(a) not bit-exact equal to fp_mul(a, a)");

    for (k, &limb) in sqr.iter().take(4).enumerate() {
        assert!(limb < (1u64 << 51), "limb {k} not reduced below 2^51");
    }
});
