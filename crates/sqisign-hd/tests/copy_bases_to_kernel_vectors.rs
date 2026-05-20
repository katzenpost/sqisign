//! Differential test of `copy_bases_to_kernel`.
mod common;

use sqisign_hd::{copy_bases_to_kernel, ThetaKernelCouplePoints};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/copy_bases_to_kernel.json"
);

#[test]
fn copy_bases_to_kernel_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::copy_bases_to_kernel");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let b1 = common::ec_basis_from("b1", &v.inputs);
        let b2 = common::ec_basis_from("b2", &v.inputs);
        let exp = common::theta_kernel_couple_points_from("ker", &v.outputs);
        let mut got = ThetaKernelCouplePoints::zero();
        copy_bases_to_kernel(&mut got, &b1, &b2);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
