//! Differential test of `quat_alg_elem_is_zero`.
mod common;
use common::{read_i32, read_ibz};
use sqisign_quaternion::{quat_alg_elem_is_zero, QuatAlgElem};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_elem_is_zero.json"
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
fn quat_alg_elem_is_zero_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_elem_is_zero");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let x = read_elem("x", &v.inputs);
        let exp = read_i32("res", &v.outputs);
        // The C `res` is an `int` cast to bytes; only the low byte is
        // significant for this 0/1 predicate. Compare on the low byte to
        // tolerate sign-extension differences across platforms.
        let r = quat_alg_elem_is_zero(&x);
        assert_eq!(r & 0xff, exp & 0xff, "vector {}", v.id);
    }
}
