//! Differential test of `quat_lideal_create_principal`.
mod common;
use common::{ibz_eq, mat4x4_eq, read_ibz, read_lattice};
use sqisign_quaternion::{
    quat_lattice_O0_set, quat_lideal_create_principal, QuatAlg, QuatAlgElem, QuatLattice,
    QuatLeftIdeal,
};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lideal_create_principal.json"
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
fn quat_lideal_create_principal_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::quat_lideal_create_principal"
    );
    for v in &f.vectors {
        let x = read_elem("x", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let exp_norm = read_ibz("L_norm", &v.outputs);
        let exp_lat = read_lattice("L_lat", &v.outputs);
        let alg = QuatAlg::init_set(&p);
        let mut o0 = QuatLattice::new();
        quat_lattice_O0_set(&mut o0);
        let mut l = QuatLeftIdeal::new();
        quat_lideal_create_principal(&mut l, &x, &o0, &alg);
        assert!(ibz_eq(&l.norm, &exp_norm), "vector {}: norm", v.id);
        assert!(
            ibz_eq(&l.lattice.denom, &exp_lat.denom),
            "vector {}: lat.denom",
            v.id
        );
        assert!(
            mat4x4_eq(&l.lattice.basis, &exp_lat.basis),
            "vector {}: lat.basis",
            v.id
        );
    }
}
