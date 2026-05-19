//! Invariant fuzz target for `sqisign_mp::mp_mod_2exp`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the per-limb low-bit-mask identity, idempotence, and the
//! over-width no-op, for an arbitrary width and `e`. Linking `mp.c` for
//! byte-equality vs C is the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_mod_2exp;

fuzz_target!(|data: &[u8]| {
    if data.len() < 5 {
        return;
    }
    let e = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) % 8192;
    let v: Vec<u64> = data[4..]
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

    let full = 64 * v.len() as u32;
    let mut a = v.clone();
    mp_mod_2exp(&mut a, e);

    for (i, (&got, &was)) in a.iter().zip(v.iter()).enumerate() {
        let bit_lo = i as u32 * 64;
        let expect = if e >= full || bit_lo + 64 <= e {
            was
        } else if bit_lo >= e {
            0
        } else {
            was & ((1u64 << (e - bit_lo)) - 1)
        };
        assert_eq!(got, expect, "limb {i} wrong (e={e})");
    }

    // Idempotent.
    let once = a.clone();
    mp_mod_2exp(&mut a, e);
    assert_eq!(a, once, "mp_mod_2exp not idempotent");
});
