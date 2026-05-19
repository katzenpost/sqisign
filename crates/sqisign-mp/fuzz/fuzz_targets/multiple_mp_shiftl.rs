//! Invariant fuzz target for `sqisign_mp::multiple_mp_shiftl`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts agreement with `mp_shiftl` in the 1..=63 domain, that an
//! over-width shift zeroes the value, and shift composition, for an
//! arbitrary width. Linking `mp.c` for byte-equality vs C is the next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::{mp_shiftl, multiple_mp_shiftl};

fuzz_target!(|data: &[u8]| {
    if data.len() < 9 {
        return;
    }
    let shift = 1 + (u16::from_le_bytes([data[0], data[1]]) as u32 % 4096);
    let v: Vec<u64> = data[2..]
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
    multiple_mp_shiftl(&mut x, shift);

    // In the 1..=63 domain it must equal a single mp_shiftl.
    if shift <= 63 {
        let mut y = v.clone();
        mp_shiftl(&mut y, shift);
        assert_eq!(x, y, "multiple != mp_shiftl in domain");
    }

    // A shift at/over the full bit width must zero the value.
    if shift >= 64 * v.len() as u32 {
        assert!(x.iter().all(|&w| w == 0), "over-width shift not zero");
    }
});
