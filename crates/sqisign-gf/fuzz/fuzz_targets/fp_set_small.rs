//! Invariant fuzz target for `sqisign_gf::fp_set_small`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_set_small` on a fuzzing host.
//!
//! `fp_set_small` is the reference's `modint((int)val, *x)`: narrow
//! `val` to `int32`, sign-extend to `u64` at limb 0, zero limbs 1..=4,
//! then `nres` to Montgomery. Three invariants survive the bit-exact
//! form into property space and are exercised by the harness:
//!
//! 1. `val == 0` produces the canonical all-zero limb vector,
//!    regardless of pre-fill (the Montgomery image of zero is zero).
//! 2. `val == 1` produces the Montgomery `ONE` constant
//!    `[0x19, 0, 0, 0, 0x300000000000]`, regardless of pre-fill (the
//!    same constant `fp_set_one` writes directly; this is the
//!    independent algorithmic confirmation that `NRES_C` is correct
//!    through the boundary).
//! 3. High-bits-ignored narrowing: for any `val`,
//!    `fp_set_small(out, val) == fp_set_small(out, val as i32 as u64)`
//!    bit-for-bit. The C wrapper drops the high 32 bits of `val`
//!    before calling `modint`, so the boundary is insensitive to them.
//!
//! The fuzz input is split into a five-limb destination pre-fill and a
//! single u64 `val`, both consumed verbatim. Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_set_one, fp_set_small, Fp, NWORDS_FIELD};

const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
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
    // 5 * 8 (prefill) + 8 (val) = 48 bytes minimum.
    if data.len() < NWORDS_FIELD * 8 + 8 {
        return;
    }
    let prefill = fp(&data[..NWORDS_FIELD * 8]);
    let val_bytes: [u8; 8] = data[NWORDS_FIELD * 8..NWORDS_FIELD * 8 + 8]
        .try_into()
        .unwrap();
    let val = u64::from_le_bytes(val_bytes);

    // (1) val == 0 -> canonical zero.
    let mut zero_out: Fp = prefill;
    fp_set_small(&mut zero_out, 0);
    assert_eq!(
        zero_out,
        [0u64; NWORDS_FIELD],
        "fp_set_small(_, 0) must produce the canonical all-zero limb vector"
    );

    // (2) val == 1 -> Montgomery ONE, bit-equal to fp_set_one.
    let mut one_small: Fp = prefill;
    fp_set_small(&mut one_small, 1);
    let mut one_set: Fp = prefill;
    fp_set_one(&mut one_set);
    assert_eq!(
        one_small, MONTGOMERY_ONE,
        "fp_set_small(_, 1) must produce the Montgomery ONE constant"
    );
    assert_eq!(
        one_small, one_set,
        "fp_set_small(_, 1) must equal fp_set_one(_) bit-exact"
    );

    // (3) High-bits-ignored narrowing.
    let narrowed = (val as i32) as u64;
    let mut full_out: Fp = prefill;
    fp_set_small(&mut full_out, val);
    let mut narrow_out: Fp = prefill;
    fp_set_small(&mut narrow_out, narrowed);
    assert_eq!(
        full_out, narrow_out,
        "fp_set_small must ignore the high 32 bits of val (only the int32 narrowing is observable)"
    );
});
