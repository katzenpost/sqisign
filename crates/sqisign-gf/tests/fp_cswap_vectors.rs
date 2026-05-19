//! Differential test of the ported `fp_cswap` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs per operand
//! (the reference's `digit_t = uint64_t` memory layout,
//! `NWORDS_FIELD == 5`), recorded for *both* operands twice: once
//! pre-call (`g_in`, `f_in`) and once post-call (`g_out`, `f_out`).
//!
//! `fp_cswap` is a branchless constant-time conditional swap whose
//! contract consults only the LSB of `ctl`: `ctl & 1 == 0` is a no-op,
//! `ctl & 1 == 1` swaps `a` and `b` limb for limb. The battery exercises
//! both endpoints on every `(g, f)` pair and additionally records two
//! non-canonical ctl values (`0xfffffffe` and `0xffffffff`) so the
//! LSB-only narrowing is pinned at the recorded boundary: `0xfffffffe`
//! must behave as `0`, `0xffffffff` must behave as `1`.
//!
//! Four assertions: (1) bit-equality of `g_out` AND `f_out` to the
//! recorded reference outputs on every vector; (2) the LSB partition
//! (count of `ctl & 1 == 0` records + count of `ctl & 1 == 1` records
//! == total), so no record was recorded outside the LSB binary split;
//! (3) both endpoint counts are strictly positive, so the battery
//! actually exercises both contract endpoints; (4) at least one record
//! at each of the two non-canonical ctl values (`0xfffffffe`,
//! `0xffffffff`) is present, pinning that the LSB-only contract is
//! covered.

use sqisign_gf::{fp_cswap, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_cswap.json"
);

fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be exactly {NWORDS_FIELD} u64 limbs"
    );
    let mut limbs = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        limbs[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    limbs
}

fn ctl_from(bytes: &[u8]) -> u32 {
    assert_eq!(bytes.len(), 4, "ctl must be exactly 4 bytes");
    u32::from_le_bytes(bytes.try_into().unwrap())
}

#[test]
fn fp_cswap_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_cswap");
    let total = file.vectors.len();
    assert!(
        total >= 1000,
        "expected the full battery, found {total} vectors"
    );

    let mut lsb_zero_records = 0usize;
    let mut lsb_one_records = 0usize;
    let mut nonlsb_zero_seen = false;
    let mut nonlsb_one_seen = false;

    for v in &file.vectors {
        let g_in = fp_from(
            "g_in",
            &decode("g_in", &v.inputs["g_in"]).expect("g_in hex"),
        );
        let f_in = fp_from(
            "f_in",
            &decode("f_in", &v.inputs["f_in"]).expect("f_in hex"),
        );
        let ctl = ctl_from(&decode("ctl", &v.inputs["ctl"]).expect("ctl hex"));
        let g_expected = fp_from(
            "g_out",
            &decode("g_out", &v.outputs["g_out"]).expect("g_out hex"),
        );
        let f_expected = fp_from(
            "f_out",
            &decode("f_out", &v.outputs["f_out"]).expect("f_out hex"),
        );

        // Replay: feed the recorded pre-call limbs back into the port
        // and compare both post-call limb vectors against the recorded
        // outputs. A port that wrote only g (or only f) would diverge
        // visibly on at least one of the two assertions.
        let mut g: Fp = g_in;
        let mut f: Fp = f_in;
        fp_cswap(&mut g, &mut f, ctl);
        assert_eq!(
            g, g_expected,
            "vector {} (ctl=0x{:08x}) diverged from the C reference on g_out",
            v.id, ctl
        );
        assert_eq!(
            f, f_expected,
            "vector {} (ctl=0x{:08x}) diverged from the C reference on f_out",
            v.id, ctl
        );

        if (ctl & 1) == 0 {
            lsb_zero_records += 1;
            if ctl != 0 {
                nonlsb_zero_seen = true;
            }
        } else {
            lsb_one_records += 1;
            if ctl != 1 {
                nonlsb_one_seen = true;
            }
        }
    }

    // LSB partition pin: every record falls into one of the two LSB
    // halves, and both halves are exercised.
    assert_eq!(
        lsb_zero_records + lsb_one_records,
        total,
        "every recorded ctl must partition by its LSB"
    );
    assert!(
        lsb_zero_records > 0 && lsb_one_records > 0,
        "both LSB endpoints must be exercised: lsb_zero={lsb_zero_records}, lsb_one={lsb_one_records}"
    );
    // LSB-only contract pin: the battery records at least one ctl with
    // LSB clear but higher bits set, and at least one with LSB set and
    // higher bits set. Without these the differential would not
    // distinguish "consume the LSB" from "consume the full mask".
    assert!(
        nonlsb_zero_seen,
        "battery must record at least one ctl with LSB clear and higher bits set"
    );
    assert!(
        nonlsb_one_seen,
        "battery must record at least one ctl with LSB set and higher bits set"
    );
}
