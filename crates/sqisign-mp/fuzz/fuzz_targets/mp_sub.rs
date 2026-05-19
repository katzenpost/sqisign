//! Invariant fuzz target for `sqisign_mp::mp_sub`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the mutual-inverse and negation laws against `mp_add` and that
//! a single limb matches native wrapping subtraction, for an arbitrary
//! limb count. Linking `mp.c` for byte-equality vs C is the next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::{mp_add, mp_sub};

fn limbs(bytes: &[u8]) -> Vec<u64> {
    bytes
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    let half = data.len() / 2;
    let a = limbs(&data[..half]);
    let mut b = limbs(&data[half..half + half]);
    if a.is_empty() {
        return;
    }
    b.resize(a.len(), 0);

    // (a + b) - b == a.
    let mut s = vec![0u64; a.len()];
    mp_add(&mut s, &a, &b);
    let mut back = vec![0u64; a.len()];
    mp_sub(&mut back, &s, &b);
    assert_eq!(back, a, "(a + b) - b != a");

    // a - a == 0.
    let mut z = vec![0u64; a.len()];
    mp_sub(&mut z, &a, &a);
    assert!(z.iter().all(|&w| w == 0), "a - a != 0");

    if a.len() == 1 {
        let mut c = [0u64; 1];
        mp_sub(&mut c, &a, &b);
        assert_eq!(c[0], a[0].wrapping_sub(b[0]), "single limb != wrapping_sub");
    }
});
