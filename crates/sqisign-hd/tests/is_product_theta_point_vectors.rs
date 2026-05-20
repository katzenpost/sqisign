//! Differential test of `is_product_theta_point`.
mod common;

use sqisign_hd::is_product_theta_point;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/is_product_theta_point.json"
);

#[test]
fn is_product_theta_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::is_product_theta_point");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a = common::theta_point_from("a", &v.inputs);
        let exp = common::u32_field("result", &v.outputs);
        let got = is_product_theta_point(&a);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
