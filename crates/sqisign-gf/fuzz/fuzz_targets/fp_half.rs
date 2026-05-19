//! Invariant fuzz target for `sqisign_gf::fp_half`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_half` on a fuzzing host.
//!
//! `fp_half` is literally `modmul(TWO_INV, a, out)`; the cross-oracle
//! `fp_half(a) == fp_mul(&TWO_INV, a)` is therefore a sanity tautology
//! at the source level (both call the same `modmul` with the same
//! operands), and any divergence on the fuzz target would indicate a
//! mis-transcribed `TWO_INV` constant or a corrupt build, not an
//! algorithmic defect. The structural carry-propagation invariant
//! (limbs 0..=3 below 2^51; limb 4 left fully unmasked by design)
//! inherits from `modmul`'s final column writes. Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_half, fp_mul, Fp, NWORDS_FIELD};

/// Montgomery representative of `2^-1 mod p`, transcribed verbatim from
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:532..536`. Kept in
/// sync with the internal `TWO_INV` constant in `sqisign_gf::lib`; the
/// fuzz crate cannot access private items so the constant is duplicated
/// here. Any drift between the two is exactly what the cross-oracle
/// below catches on the first input.
const TWO_INV: Fp = [
    0x0000_0000_0000_000c,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_4000_0000_0000,
];

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

    let mut half = [0u64; NWORDS_FIELD];
    let mut mul = [0u64; NWORDS_FIELD];
    fp_half(&mut half, &a);
    fp_mul(&mut mul, &TWO_INV, &a);
    assert_eq!(
        half, mul,
        "fp_half(a) not bit-exact equal to fp_mul(&TWO_INV, a)"
    );

    for (k, &limb) in half.iter().take(4).enumerate() {
        assert!(limb < (1u64 << 51), "limb {k} not reduced below 2^51");
    }
});
