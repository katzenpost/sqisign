//! Differential test of `quat_alg_sub`.
mod common;
use common::{ibz_eq, read_ibz, vec4_eq};
use sqisign_quaternion::{quat_alg_sub, QuatAlgElem};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_sub.json"
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
fn quat_alg_sub_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_sub");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let a = read_elem("a", &v.inputs);
        let b = read_elem("b", &v.inputs);
        let exp = read_elem("r", &v.outputs);
        let mut r = QuatAlgElem::new();
        quat_alg_sub(&mut r, &a, &b);
        assert!(ibz_eq(&r.denom, &exp.denom), "vector {}: denom", v.id);
        assert!(vec4_eq(&r.coord, &exp.coord), "vector {}: coord", v.id);
    }
}
