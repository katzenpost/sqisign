//! Differential test of `ibz_2x2_mul_mod`.
mod common;
use common::{ibz_eq, read_ibz};
use sqisign_quaternion::{ibz_2x2_mul_mod, ibz_mat_2x2_new, IbzMat2x2};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_2x2_mul_mod.json"
);

fn read_mat2x2(prefix: &str, inputs: &std::collections::BTreeMap<String, String>) -> IbzMat2x2 {
    let mut m = ibz_mat_2x2_new();
    m[0][0] = read_ibz(&format!("{prefix}_00"), inputs);
    m[0][1] = read_ibz(&format!("{prefix}_01"), inputs);
    m[1][0] = read_ibz(&format!("{prefix}_10"), inputs);
    m[1][1] = read_ibz(&format!("{prefix}_11"), inputs);
    m
}

fn mat2_eq(a: &IbzMat2x2, b: &IbzMat2x2) -> bool {
    (0..2).all(|i| (0..2).all(|j| ibz_eq(&a[i][j], &b[i][j])))
}

#[test]
fn ibz_2x2_mul_mod_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_2x2_mul_mod");
    for v in &f.vectors {
        let a = read_mat2x2("a", &v.inputs);
        let b = read_mat2x2("b", &v.inputs);
        let m = read_ibz("m", &v.inputs);
        let exp = read_mat2x2("r", &v.outputs);
        let mut r = ibz_mat_2x2_new();
        ibz_2x2_mul_mod(&mut r, &a, &b, &m);
        assert!(mat2_eq(&r, &exp), "vector {}", v.id);
    }
}
