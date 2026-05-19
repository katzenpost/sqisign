//! Invariant fuzz target for `sqisign_gf::fp_sub`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_sub` on a fuzzing host.
//!
//! Subtraction is not commutative, so the bit-exact commutativity law
//! used for `fp_add` does not transfer. Instead this asserts the two
//! sound raw-limb invariants for arbitrary (possibly non-canonical)
//! five-limb inputs: `fp_sub(a, a)` is the canonical all-zero limb
//! vector (modsub's limbwise difference is exactly zero, prop signals no
//! borrow, no correction fires), and the structural carry-propagation
//! invariant that limbs 0..=3 are below 2^51 (limb 4 is left unmasked by
//! design and is not constrained). Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_sub, Fp, NWORDS_FIELD};

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
    if data.len() < 2 * NWORDS_FIELD * 8 {
        return;
    }
    let a = fp(&data[..NWORDS_FIELD * 8]);
    let b = fp(&data[NWORDS_FIELD * 8..2 * NWORDS_FIELD * 8]);

    let mut diff = [0u64; NWORDS_FIELD];
    fp_sub(&mut diff, &a, &b);
    for (k, &limb) in diff.iter().take(4).enumerate() {
        assert!(limb < (1u64 << 51), "limb {k} not reduced below 2^51");
    }

    let mut self_diff = [0u64; NWORDS_FIELD];
    fp_sub(&mut self_diff, &a, &a);
    assert_eq!(
        self_diff,
        [0u64; NWORDS_FIELD],
        "fp_sub(a, a) not the canonical all-zero representative"
    );
});
