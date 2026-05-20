//! Differential test of `id2iso_ideal_to_kernel_dlogs_even`.

mod common;

use sqisign_id2iso::id2iso_ideal_to_kernel_dlogs_even;
use sqisign_precomp::QUATALG_PINFTY;
use sqisign_quaternion::ibz_vec_2_new;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/id2iso_ideal_to_kernel_dlogs_even.json"
);

#[test]
fn id2iso_ideal_to_kernel_dlogs_even_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_id2iso::id2iso_ideal_to_kernel_dlogs_even"
    );
    assert!(!f.vectors.is_empty());
    let alg = &*QUATALG_PINFTY;
    for v in &f.vectors {
        let l = common::read_lideal("l", &v.inputs);
        let exp = common::read_vec2("vec", &v.outputs);

        let mut out = ibz_vec_2_new();
        id2iso_ideal_to_kernel_dlogs_even(&mut out, &l, alg);
        assert!(common::vec2_eq(&out, &exp), "vector {} diverged", v.id);
    }
}
