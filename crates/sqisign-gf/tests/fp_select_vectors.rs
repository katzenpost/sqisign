//! Differential test of the ported `fp_select` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! `fp_select` is a branchless constant-time conditional select whose
//! contract restricts `ctl` to the two endpoints `0x00000000` (select
//! `a0`) and `0xFFFFFFFF` (select `a1`); any other `ctl` value is
//! undefined per the reference. The battery therefore drives `ctl`
//! through *both* endpoints on every `(a0, a1)` pair. The record's
//! `prefill` field is the destination pre-fill the C harness fed into
//! `fp_select`; feeding the ported function the same pre-fill and
//! comparing the resulting limbs against the C-recorded output is the
//! only way to catch a no-op or partial-write port at the `fp_t`
//! boundary.
//!
//! Three assertions: (1) bit-equality to the recorded reference output
//! on every vector; (2) the count of `ctl == 0` records and the count
//! of `ctl == 0xFFFFFFFF` records sum to the total, pinning that no
//! record was recorded at a third undefined `ctl` value; (3) both
//! endpoint counts are strictly positive, so the battery actually
//! exercises both contract endpoints rather than silently degenerating
//! to one.

use sqisign_gf::{fp_select, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_select.json"
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
fn fp_select_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_select");
    let total = file.vectors.len();
    assert!(
        total >= 1000,
        "expected the full battery, found {total} vectors"
    );

    let mut ctl_zero_records = 0usize;
    let mut ctl_ones_records = 0usize;

    for v in &file.vectors {
        let a0 = fp_from("a0", &decode("a0", &v.inputs["a0"]).expect("a0 hex"));
        let a1 = fp_from("a1", &decode("a1", &v.inputs["a1"]).expect("a1 hex"));
        let ctl = ctl_from(&decode("ctl", &v.inputs["ctl"]).expect("ctl hex"));
        let prefill = fp_from(
            "prefill",
            &decode("prefill", &v.inputs["prefill"]).expect("prefill hex"),
        );
        let expected = fp_from("d", &decode("d", &v.outputs["d"]).expect("d hex"));

        // Pre-fill the destination exactly as the C harness did; a port
        // that quietly skipped any limb would leave the corresponding
        // pre-fill byte visible and diverge from the recorded output.
        let mut d: Fp = prefill;
        fp_select(&mut d, &a0, &a1, ctl);
        assert_eq!(
            d, expected,
            "vector {} (ctl=0x{:08x}) diverged from the C reference at the fp_t boundary",
            v.id, ctl
        );

        match ctl {
            0x00000000 => ctl_zero_records += 1,
            0xFFFFFFFF => ctl_ones_records += 1,
            other => panic!(
                "vector {} recorded an undefined ctl value 0x{:08x}; the reference contract \
                 restricts ctl to 0x00000000 or 0xFFFFFFFF",
                v.id, other
            ),
        }
    }

    // Count pin: every record falls into one of the two declared
    // endpoints, and both endpoints are exercised.
    assert_eq!(
        ctl_zero_records + ctl_ones_records,
        total,
        "every recorded ctl must be one of the two declared endpoints"
    );
    assert!(
        ctl_zero_records > 0 && ctl_ones_records > 0,
        "both contract endpoints must be exercised: ctl_zero={ctl_zero_records}, ctl_ones={ctl_ones_records}"
    );
}
