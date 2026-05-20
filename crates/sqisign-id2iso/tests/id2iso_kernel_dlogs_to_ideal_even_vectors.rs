//! Differential test of `id2iso_kernel_dlogs_to_ideal_even`.

mod common;

use sqisign_id2iso::id2iso_kernel_dlogs_to_ideal_even;
use sqisign_precomp::QUATALG_PINFTY;
use sqisign_quaternion::QuatLeftIdeal;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/id2iso_kernel_dlogs_to_ideal_even.json"
);

#[test]
fn id2iso_kernel_dlogs_to_ideal_even_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_id2iso::id2iso_kernel_dlogs_to_ideal_even"
    );
    assert!(!f.vectors.is_empty());
    let alg = &*QUATALG_PINFTY;
    for v in &f.vectors {
        let vec2 = common::read_vec2("vec2", &v.inputs);
        let fparam = common::read_i32("f", &v.inputs);
        let exp = common::read_lideal("l", &v.outputs);

        let mut l = QuatLeftIdeal::new();
        id2iso_kernel_dlogs_to_ideal_even(&mut l, &vec2, fparam, alg);
        assert!(common::lideal_eq(&l, &exp), "vector {} diverged", v.id);
    }
}
