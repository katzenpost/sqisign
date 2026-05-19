//! Invariant fuzz target for `sqisign_mp::mp_compare`: reflexive,
//! antisymmetric, and consistent with the top differing limb, for an
//! arbitrary width. Linking `mp.c` for byte-equality vs C is the next
//! increment. Runner deferred (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_compare;

fuzz_target!(|data: &[u8]| {
    let half = data.len() / 2;
    let a: Vec<u64> = data[..half]
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    if a.is_empty() {
        return;
    }
    let mut b: Vec<u64> = data[half..half + half]
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    b.resize(a.len(), 0);

    let ab = mp_compare(&a, &b);
    assert_eq!(ab, -mp_compare(&b, &a), "not antisymmetric");
    assert_eq!(mp_compare(&a, &a), 0, "not reflexive");

    let mut expect = 0i32;
    for i in (0..a.len()).rev() {
        if a[i] != b[i] {
            expect = if a[i] > b[i] { 1 } else { -1 };
            break;
        }
    }
    assert_eq!(ab, expect, "not consistent with top differing limb");
});
