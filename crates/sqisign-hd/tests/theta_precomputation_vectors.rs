//! Differential test of `theta_precomputation`.
mod common;

use sqisign_hd::theta_precomputation;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/theta_precomputation.json"
);

#[test]
fn theta_precomputation_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::theta_precomputation");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let mut a = common::theta_structure_from("a", &v.inputs);
        let exp = common::theta_structure_from("c", &v.outputs);
        // The C harness explicitly forces precomputation=false before
        // invoking the boundary, so the Rust input has that field set
        // accordingly too; mirror it here.
        a.precomputation = false;
        theta_precomputation(&mut a);
        assert_eq!(a, exp, "vector {} diverged", v.id);
    }
}
