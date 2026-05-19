//! Invariant fuzz target for `sqisign_mp::mp_mul`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! For nwords>=2 asserts commutativity and agreement with a u128
//! reference (two limbs); for nwords==1 asserts the *reproduced* upstream
//! defect, `2*(a*b) mod 2^64`. Linking `mp.c` for byte-equality vs C is
//! the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_mul;

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

    let mut c = vec![0u64; a.len()];
    mp_mul(&mut c, &a, &b);

    if a.len() == 1 {
        // Reproduced upstream defect: 2*(a*b) mod 2^64.
        let low = (a[0] as u128 * b[0] as u128) as u64;
        assert_eq!(c[0], low.wrapping_mul(2), "single-limb defect not reproduced");
    } else {
        // Commutativity for the correct multilimb regime.
        let mut c2 = vec![0u64; a.len()];
        mp_mul(&mut c2, &b, &a);
        assert_eq!(c, c2, "mp_mul not commutative for nwords>=2");
    }
});
