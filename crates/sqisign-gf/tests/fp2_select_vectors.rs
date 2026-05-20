//! Differential test of the ported `fp2_select` against the committed
//! C-derived vectors. Record shape:
//! (a0_re, a0_im, a1_re, a1_im, ctl, prefill_re, prefill_im) ->
//! (d_re, d_im).
//!
//! Same `ctl` contract as `fp_select`: only `0x00000000` (select `a0`)
//! and `0xFFFFFFFF` (select `a1`) are exercised; any other `ctl` is
//! undefined per the reference. The pre-fills exercise the
//! complete-overwrite property on both halves.

use sqisign_gf::{fp2_select, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_select.json"
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

#[test]
fn fp2_select_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_select");
    let total = file.vectors.len();
    assert!(total >= 1000);

    let mut ctl_zero = 0usize;
    let mut ctl_ones = 0usize;
    for v in &file.vectors {
        let a0_re = fp_from(
            "a0_re",
            &decode("a0_re", &v.inputs["a0_re"]).expect("a0_re hex"),
        );
        let a0_im = fp_from(
            "a0_im",
            &decode("a0_im", &v.inputs["a0_im"]).expect("a0_im hex"),
        );
        let a1_re = fp_from(
            "a1_re",
            &decode("a1_re", &v.inputs["a1_re"]).expect("a1_re hex"),
        );
        let a1_im = fp_from(
            "a1_im",
            &decode("a1_im", &v.inputs["a1_im"]).expect("a1_im hex"),
        );
        let pre_re = fp_from(
            "prefill_re",
            &decode("prefill_re", &v.inputs["prefill_re"]).expect("prefill_re hex"),
        );
        let pre_im = fp_from(
            "prefill_im",
            &decode("prefill_im", &v.inputs["prefill_im"]).expect("prefill_im hex"),
        );
        let ctl_bytes = decode("ctl", &v.inputs["ctl"]).expect("ctl hex");
        assert_eq!(ctl_bytes.len(), 4);
        let ctl = u32::from_le_bytes(ctl_bytes.try_into().unwrap());
        let exp_re = fp_from(
            "d_re",
            &decode("d_re", &v.outputs["d_re"]).expect("d_re hex"),
        );
        let exp_im = fp_from(
            "d_im",
            &decode("d_im", &v.outputs["d_im"]).expect("d_im hex"),
        );

        let a0 = Fp2 {
            re: a0_re,
            im: a0_im,
        };
        let a1 = Fp2 {
            re: a1_re,
            im: a1_im,
        };
        let mut d = Fp2 {
            re: pre_re,
            im: pre_im,
        };
        fp2_select(&mut d, &a0, &a1, ctl);
        assert_eq!(
            d.re, exp_re,
            "vector {} (ctl=0x{:08x}) diverged on d_re",
            v.id, ctl
        );
        assert_eq!(
            d.im, exp_im,
            "vector {} (ctl=0x{:08x}) diverged on d_im",
            v.id, ctl
        );

        match ctl {
            0x00000000 => ctl_zero += 1,
            0xFFFFFFFF => ctl_ones += 1,
            other => panic!(
                "vector {}: undefined ctl 0x{:08x} recorded; the reference contract \
                 restricts ctl to 0x00000000 or 0xFFFFFFFF",
                v.id, other
            ),
        }
    }
    assert_eq!(ctl_zero + ctl_ones, total);
    assert!(
        ctl_zero > 0 && ctl_ones > 0,
        "both endpoints must be exercised: ctl_zero={ctl_zero}, ctl_ones={ctl_ones}"
    );
}
