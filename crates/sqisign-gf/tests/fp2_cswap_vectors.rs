//! Differential test of the ported `fp2_cswap` against the committed
//! C-derived vectors. Record shape:
//! (g_in_re, g_in_im, f_in_re, f_in_im, ctl) ->
//! (g_out_re, g_out_im, f_out_re, f_out_im).
//!
//! Same `ctl & 1` LSB-only contract as `fp_cswap`. The battery
//! records four ctl values (0, 1, 0xfffffffe, 0xffffffff) so the
//! LSB-only narrowing is pinned at the recorded boundary.

use sqisign_gf::{fp2_cswap, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_cswap.json"
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
fn fp2_cswap_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_cswap");
    let total = file.vectors.len();
    assert!(total >= 1000);

    let mut lsb_zero = 0usize;
    let mut lsb_one = 0usize;
    let mut nonlsb_zero_seen = false;
    let mut nonlsb_one_seen = false;
    for v in &file.vectors {
        let g_in_re = fp_from(
            "g_in_re",
            &decode("g_in_re", &v.inputs["g_in_re"]).expect("g_in_re hex"),
        );
        let g_in_im = fp_from(
            "g_in_im",
            &decode("g_in_im", &v.inputs["g_in_im"]).expect("g_in_im hex"),
        );
        let f_in_re = fp_from(
            "f_in_re",
            &decode("f_in_re", &v.inputs["f_in_re"]).expect("f_in_re hex"),
        );
        let f_in_im = fp_from(
            "f_in_im",
            &decode("f_in_im", &v.inputs["f_in_im"]).expect("f_in_im hex"),
        );
        let ctl_bytes = decode("ctl", &v.inputs["ctl"]).expect("ctl hex");
        assert_eq!(ctl_bytes.len(), 4);
        let ctl = u32::from_le_bytes(ctl_bytes.try_into().unwrap());
        let g_exp_re = fp_from(
            "g_out_re",
            &decode("g_out_re", &v.outputs["g_out_re"]).expect("g_out_re hex"),
        );
        let g_exp_im = fp_from(
            "g_out_im",
            &decode("g_out_im", &v.outputs["g_out_im"]).expect("g_out_im hex"),
        );
        let f_exp_re = fp_from(
            "f_out_re",
            &decode("f_out_re", &v.outputs["f_out_re"]).expect("f_out_re hex"),
        );
        let f_exp_im = fp_from(
            "f_out_im",
            &decode("f_out_im", &v.outputs["f_out_im"]).expect("f_out_im hex"),
        );

        let mut g = Fp2 {
            re: g_in_re,
            im: g_in_im,
        };
        let mut f = Fp2 {
            re: f_in_re,
            im: f_in_im,
        };
        fp2_cswap(&mut g, &mut f, ctl);
        assert_eq!(g.re, g_exp_re, "vector {} ctl=0x{:08x} g_out_re", v.id, ctl);
        assert_eq!(g.im, g_exp_im, "vector {} ctl=0x{:08x} g_out_im", v.id, ctl);
        assert_eq!(f.re, f_exp_re, "vector {} ctl=0x{:08x} f_out_re", v.id, ctl);
        assert_eq!(f.im, f_exp_im, "vector {} ctl=0x{:08x} f_out_im", v.id, ctl);

        if (ctl & 1) == 0 {
            lsb_zero += 1;
            if ctl != 0 {
                nonlsb_zero_seen = true;
            }
        } else {
            lsb_one += 1;
            if ctl != 1 {
                nonlsb_one_seen = true;
            }
        }
    }
    assert_eq!(lsb_zero + lsb_one, total);
    assert!(lsb_zero > 0 && lsb_one > 0);
    assert!(
        nonlsb_zero_seen,
        "battery must record at least one ctl with LSB clear and higher bits set"
    );
    assert!(
        nonlsb_one_seen,
        "battery must record at least one ctl with LSB set and higher bits set"
    );
}
