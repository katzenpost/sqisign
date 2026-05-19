//! Invariant fuzz target for `sqisign_gf::fp_is_zero`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_is_zero` on a fuzzing host.
//!
//! Asserts the sound raw-mask invariants for arbitrary (possibly
//! non-canonical) five-limb inputs: the returned `u32` is either `0` or
//! `0xFFFFFFFF` (the C `-(uint32_t)int01` invariant; any other bit
//! pattern would indicate the port mishandled the negation cast chain
//! that turns modis0's `{0, 1}` result into the `{0, 0xFFFFFFFF}` mask
//! the rest of the codebase consumes), and `fp_is_zero` of the
//! canonical all-zero limb pattern is the all-ones mask. The all-zero
//! check is driven from a fixed all-zero input rather than the fuzzed
//! bytes, because `fp_is_zero(a)` for arbitrary fuzzed `a` is almost
//! always `0` (random 320-bit values have negligible probability of
//! landing on a representative of `0 mod p`). Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_is_zero, Fp, NWORDS_FIELD};

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

    let m = fp_is_zero(&a);
    assert!(
        m == 0 || m == 0xFFFF_FFFF,
        "fp_is_zero returned non-mask value {m:#x}"
    );

    let zero = fp_is_zero(&[0u64; NWORDS_FIELD]);
    assert_eq!(
        zero, 0xFFFF_FFFF,
        "fp_is_zero(canonical zero) not the all-ones mask"
    );
});
