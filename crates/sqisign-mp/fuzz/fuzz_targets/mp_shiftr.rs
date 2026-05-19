//! Invariant fuzz target for `sqisign_mp::mp_shiftr`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the returned bit is the entry parity, the top limb is
//! zero-filled (logical shift), and a single limb matches native `>>`,
//! for an arbitrary width and a shift in 1..=63. Linking `mp.c` for
//! byte-equality vs C is the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_shiftr;

fuzz_target!(|data: &[u8]| {
    if data.len() < 9 {
        return;
    }
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

    let entry_bit = v[0] & 1;
    let mut x = v.clone();
    let bit = mp_shiftr(&mut x, shift);

    assert_eq!(bit, entry_bit, "returned bit != entry x[0]&1");
    assert_eq!(
        *x.last().unwrap() >> (64 - shift),
        0,
        "top not zero-filled (not a logical shift)"
    );
    if v.len() == 1 {
        assert_eq!(x[0], v[0] >> shift, "single limb != native >>");
    }
});
