//! Invariant fuzz target for `sqisign_gf::fp_copy`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_copy` on a fuzzing host.
//!
//! `fp_copy` is the reference's `modcpy`: a plain five-limb assignment.
//! The asserted invariant is therefore the only one that transfers from
//! the bit-exact form, namely that the output equals the input limb for
//! limb on arbitrary (possibly non-canonical) inputs. A non-trivial
//! pre-fill in the destination guards against a no-op or partial-write
//! port. Linking the reference's `fp_p5248_64.c` for byte-equality
//! against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_copy, Fp, NWORDS_FIELD};

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

    let mut out: Fp = [0xa5a5_a5a5_a5a5_a5a5u64; NWORDS_FIELD];
    fp_copy(&mut out, &a);
    assert_eq!(out, a, "fp_copy is not the bit-exact identity");
});
