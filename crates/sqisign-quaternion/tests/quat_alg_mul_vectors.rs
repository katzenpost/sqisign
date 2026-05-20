//! Differential test of `quat_alg_mul`.
use sqisign_quaternion::{quat_alg_mul, Ibz, QuatAlg, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_mul.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_elem(prefix: &str, inputs: &std::collections::BTreeMap<String, String>) -> QuatAlgElem {
    let denom_key = format!("{prefix}_denom");
    QuatAlgElem {
        denom: read_ibz(&denom_key, &inputs[&denom_key]),
        coord: [
            read_ibz(&format!("{prefix}_c0"), &inputs[&format!("{prefix}_c0")]),
            read_ibz(&format!("{prefix}_c1"), &inputs[&format!("{prefix}_c1")]),
            read_ibz(&format!("{prefix}_c2"), &inputs[&format!("{prefix}_c2")]),
            read_ibz(&format!("{prefix}_c3"), &inputs[&format!("{prefix}_c3")]),
        ],
    }
}

#[test]
fn quat_alg_mul_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_mul");
    assert!(f.vectors.len() >= 80);
    for v in &f.vectors {
        let a = read_elem("a", &v.inputs);
        let b = read_elem("b", &v.inputs);
        let p = read_ibz("p", &v.inputs["p"]);
        let exp = read_elem("r", &v.outputs);
        let alg = QuatAlg::init_set(&p);
        let mut r = QuatAlgElem::new();
        quat_alg_mul(&mut r, &a, &b, &alg);
        assert_eq!(r.denom.0, exp.denom.0, "vector {}: denom", v.id);
        for k in 0..4 {
            assert_eq!(r.coord[k].0, exp.coord[k].0, "vector {}: coord {}", v.id, k);
        }
    }
}
