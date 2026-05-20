//! Differential test of `ec_curve_to_basis_2f_to_hint`. Inputs: curve,
//! `f`. Outputs: basis, mutated curve, hint byte.
mod common;

use sqisign_ec::{ec_curve_to_basis_2f_to_hint, EcBasis};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_curve_to_basis_2f_to_hint.json"
);

#[test]
fn ec_curve_to_basis_2f_to_hint_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_curve_to_basis_2f_to_hint");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let mut e = common::ec_curve_from("e_in", &v.inputs);
        let f = common::i32_field("f", &v.inputs);
        let exp_b = common::ec_basis_from("b", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);
        let hint_bytes = decode("hint", &v.outputs["hint"]).expect("hint hex");
        assert_eq!(hint_bytes.len(), 1);
        let exp_hint = hint_bytes[0];

        let mut b = EcBasis::zero();
        let got_hint = ec_curve_to_basis_2f_to_hint(&mut b, &mut e, f);
        assert_eq!(b, exp_b, "vector {} basis diverged", v.id);
        assert_eq!(e, exp_e, "vector {} curve diverged", v.id);
        assert_eq!(got_hint, exp_hint, "vector {} hint diverged", v.id);
    }
}
