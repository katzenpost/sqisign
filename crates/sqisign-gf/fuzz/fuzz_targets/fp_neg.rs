//! Invariant fuzz target for `sqisign_gf::fp_neg`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_neg` on a fuzzing host.
//!
//! `fp_neg` is the unary analogue of `fp_sub` (`0 - b[i]` limbwise), not
//! an involution on the redundant form, so no bit-exact
//! `fp_neg(fp_neg(a)) == a` law transfers. This asserts the two sound
//! raw-limb invariants for arbitrary (possibly non-canonical) five-limb
//! inputs: `fp_neg` of the all-zero limb vector is the all-zero limb
//! vector (modneg's limbwise `0 - 0` is exactly zero, prop signals no
//! borrow, no correction fires), and the structural carry-propagation
//! invariant that limbs 0..=3 are below 2^51 (limb 4 is left unmasked by
//! design and is not constrained). The zero check is driven from a fixed
//! all-zero input rather than the fuzzed bytes, because `fp_neg(a)` for
//! arbitrary `a` is *not* all-zero. Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_neg, Fp, NWORDS_FIELD};

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

    let mut out = [0u64; NWORDS_FIELD];
    fp_neg(&mut out, &a);
    for (k, &limb) in out.iter().take(4).enumerate() {
        assert!(limb < (1u64 << 51), "limb {k} not reduced below 2^51");
    }

    let mut zero = [0u64; NWORDS_FIELD];
    fp_neg(&mut zero, &[0u64; NWORDS_FIELD]);
    assert_eq!(
        zero,
        [0u64; NWORDS_FIELD],
        "fp_neg(0) not the canonical all-zero representative"
    );
});
