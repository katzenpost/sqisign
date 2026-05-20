//! Differential test of `ibz_mat_2x2_inv_mod`.
mod common;
use common::{ibz_eq, read_i32, read_ibz};
use sqisign_quaternion::{ibz_mat_2x2_inv_mod, ibz_mat_2x2_new, IbzMat2x2};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_2x2_inv_mod.json"
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
fn ibz_mat_2x2_inv_mod_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_2x2_inv_mod");
    for v in &f.vectors {
        let m = read_mat2x2("m", &v.inputs);
        let mod_ = read_ibz("mod", &v.inputs);
        let exp_inv = read_mat2x2("inv", &v.outputs);
        let exp_ok = read_i32("ok", &v.outputs);
        let mut inv = ibz_mat_2x2_new();
        let ok = ibz_mat_2x2_inv_mod(&mut inv, &m, &mod_);
        assert_eq!(ok & 0xff, exp_ok & 0xff, "vector {}: ok", v.id);
        assert!(mat2_eq(&inv, &exp_inv), "vector {}: inv", v.id);
    }
}
