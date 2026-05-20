//! Differential test of `quat_change_to_O0_basis`.
#![allow(non_snake_case)]
mod common;
use common::{read_ibz, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_new, quat_change_to_O0_basis, QuatAlgElem};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_change_to_O0_basis.json"
);

fn read_elem(prefix: &str, inputs: &std::collections::BTreeMap<String, String>) -> QuatAlgElem {
    let denom_key = format!("{prefix}_denom");
    QuatAlgElem {
        denom: read_ibz(&denom_key, inputs),
        coord: [
            read_ibz(&format!("{prefix}_c0"), inputs),
            read_ibz(&format!("{prefix}_c1"), inputs),
            read_ibz(&format!("{prefix}_c2"), inputs),
            read_ibz(&format!("{prefix}_c3"), inputs),
        ],
    }
}

#[test]
fn quat_change_to_O0_basis_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_change_to_O0_basis");
    for v in &f.vectors {
        let e = read_elem("e", &v.inputs);
        let exp = read_vec4("v", &v.outputs);
        let mut vv = ibz_vec_4_new();
        quat_change_to_O0_basis(&mut vv, &e);
        assert!(vec4_eq(&vv, &exp), "vector {}", v.id);
    }
}
