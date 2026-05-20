//! Differential test of `quat_alg_norm`.
use sqisign_quaternion::{quat_alg_norm, Ibz, QuatAlg, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_norm.json"
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
fn quat_alg_norm_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_norm");
    assert!(f.vectors.len() >= 100);
    for v in &f.vectors {
        let a = read_elem("a", &v.inputs);
        let p = read_ibz("p", &v.inputs["p"]);
        let exp_num = read_ibz("num", &v.outputs["num"]);
        let exp_den = read_ibz("den", &v.outputs["den"]);
        let alg = QuatAlg::init_set(&p);
        let mut num = Ibz::zero();
        let mut den = Ibz::zero();
        quat_alg_norm(&mut num, &mut den, &a, &alg);
        assert_eq!(num.0, exp_num.0, "vector {}: num", v.id);
        assert_eq!(den.0, exp_den.0, "vector {}: den", v.id);
    }
}
