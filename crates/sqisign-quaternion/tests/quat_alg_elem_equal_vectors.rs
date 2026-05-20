//! Differential test of `quat_alg_elem_equal`.
mod common;
use common::{read_i32, read_ibz};
use sqisign_quaternion::{quat_alg_elem_equal, QuatAlgElem};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_elem_equal.json"
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
fn quat_alg_elem_equal_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_elem_equal");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let a = read_elem("a", &v.inputs);
        let b = read_elem("b", &v.inputs);
        let exp = read_i32("res", &v.outputs);
        let r = quat_alg_elem_equal(&a, &b);
        assert_eq!(r, exp, "vector {}", v.id);
    }
}
