//! Differential test of `quat_alg_elem_mul_by_scalar`.
use sqisign_quaternion::{quat_alg_elem_mul_by_scalar, Ibz, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_elem_mul_by_scalar.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_elem(prefix: &str, m: &std::collections::BTreeMap<String, String>) -> QuatAlgElem {
    QuatAlgElem {
        denom: read_ibz("denom", &m[&format!("{prefix}_denom")]),
        coord: [
            read_ibz("c0", &m[&format!("{prefix}_c0")]),
            read_ibz("c1", &m[&format!("{prefix}_c1")]),
            read_ibz("c2", &m[&format!("{prefix}_c2")]),
            read_ibz("c3", &m[&format!("{prefix}_c3")]),
        ],
    }
}

#[test]
fn quat_alg_elem_mul_by_scalar_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::quat_alg_elem_mul_by_scalar"
    );
    assert!(f.vectors.len() >= 50);
    for v in &f.vectors {
        let scalar = read_ibz("scalar", &v.inputs["scalar"]);
        let e = read_elem("e", &v.inputs);
        let exp = read_elem("r", &v.outputs);
        let mut r = QuatAlgElem::new();
        quat_alg_elem_mul_by_scalar(&mut r, &scalar, &e);
        assert_eq!(r.denom.0, exp.denom.0, "vector {}: denom", v.id);
        for k in 0..4 {
            assert_eq!(r.coord[k].0, exp.coord[k].0, "vector {}: coord {}", v.id, k);
        }
    }
}
