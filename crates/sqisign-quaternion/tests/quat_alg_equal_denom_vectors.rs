//! Differential test of `quat_alg_equal_denom`.
use sqisign_quaternion::{quat_alg_equal_denom, Ibz, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_equal_denom.json"
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
fn quat_alg_equal_denom_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_equal_denom");
    assert!(f.vectors.len() >= 50);
    for v in &f.vectors {
        let a = read_elem("a", &v.inputs);
        let b = read_elem("b", &v.inputs);
        let exp_a = read_elem("ra", &v.outputs);
        let exp_b = read_elem("rb", &v.outputs);
        let mut ra = QuatAlgElem::new();
        let mut rb = QuatAlgElem::new();
        quat_alg_equal_denom(&mut ra, &mut rb, &a, &b);
        assert_eq!(ra.denom.0, exp_a.denom.0, "vector {}: ra denom", v.id);
        assert_eq!(rb.denom.0, exp_b.denom.0, "vector {}: rb denom", v.id);
        for k in 0..4 {
            assert_eq!(
                ra.coord[k].0, exp_a.coord[k].0,
                "vector {}: ra coord {}",
                v.id, k
            );
            assert_eq!(
                rb.coord[k].0, exp_b.coord[k].0,
                "vector {}: rb coord {}",
                v.id, k
            );
        }
    }
}
