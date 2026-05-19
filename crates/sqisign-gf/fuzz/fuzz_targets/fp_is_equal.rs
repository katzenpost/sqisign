//! Invariant fuzz target for `sqisign_gf::fp_is_equal`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_is_equal` on a fuzzing host.
//!
//! Asserts the sound raw-mask invariants for arbitrary (possibly
//! non-canonical) five-limb inputs: the returned `u32` is either `0` or
//! `0xFFFFFFFF` (the C `-(uint32_t)int01` invariant; any other bit
//! pattern would indicate the port mishandled the negation cast chain
//! that turns modcmp's `{0, 1}` result into the `{0, 0xFFFFFFFF}` mask
//! the rest of the codebase consumes), reflexivity `fp_is_equal(a, a)
//! == 0xFFFFFFFF` for arbitrary `a` (modcmp's redc is deterministic on
//! both operands, so the per-limb XORs are all zero and the AND-fold
//! is `1`), and symmetry `fp_is_equal(a, b) == fp_is_equal(b, a)`
//! bit-exact (XOR is symmetric, the AND-fold is commutative). Linking
//! the reference's `fp_p5248_64.c` for byte-equality against C is the
//! documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_is_equal, Fp, NWORDS_FIELD};

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

    let m_ab = fp_is_equal(&a, &b);
    assert!(
        m_ab == 0 || m_ab == 0xFFFF_FFFF,
        "fp_is_equal returned non-mask value {m_ab:#x}"
    );

    let m_ba = fp_is_equal(&b, &a);
    assert_eq!(m_ab, m_ba, "fp_is_equal not bit-exact symmetric");

    let m_aa = fp_is_equal(&a, &a);
    assert_eq!(
        m_aa, 0xFFFF_FFFF,
        "fp_is_equal(a, a) not the all-ones mask (reflexivity)"
    );
});
