//! Differential test of `quat_lattice_contains`.
mod common;
use common::{read_i32, read_ibz, read_lattice, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_new, quat_lattice_contains, QuatAlgElem};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_contains.json"
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
fn quat_lattice_contains_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_contains");
    for v in &f.vectors {
        let l = read_lattice("l", &v.inputs);
        let x = read_elem("x", &v.inputs);
        let exp_coord = read_vec4("coord", &v.outputs);
        let exp_ok = read_i32("ok", &v.outputs);
        let mut coord = ibz_vec_4_new();
        let ok = quat_lattice_contains(Some(&mut coord), &l, &x);
        assert_eq!(ok & 0xff, exp_ok & 0xff, "vector {}: ok", v.id);
        if ok != 0 {
            assert!(vec4_eq(&coord, &exp_coord), "vector {}: coord", v.id);
        }
    }
}
