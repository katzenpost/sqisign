//! Invariant fuzz target for `sqisign_mp::mp_shiftl`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the multiply-by-2^s identity (shift by 1 == self-add) and that
//! the low `shift` bits clear, for an arbitrary limb count and a shift in
//! the reference's 1..=63 domain. Linking `mp.c` for byte-equality vs C
//! is the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::{mp_add, mp_shiftl};

fuzz_target!(|data: &[u8]| {
    if data.len() < 9 {
        return;
    }
    // First byte picks a shift in 1..=63; the rest is little-endian limbs.
    let shift = 1 + (data[0] as u32 % 63);
    let v: Vec<u64> = data[1..]
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    if v.is_empty() {
        return;
    }

    let mut x = v.clone();
    mp_shiftl(&mut x, shift);

    // Low `shift` bits of the least-significant limb must be zero.
    assert_eq!(x[0] & ((1u64 << shift) - 1), 0, "low bits not cleared");

    // Shift by 1 == self-add.
    if shift == 1 {
        let mut doubled = vec![0u64; v.len()];
        mp_add(&mut doubled, &v, &v);
        assert_eq!(x, doubled, "x<<1 != x+x");
    }
});
