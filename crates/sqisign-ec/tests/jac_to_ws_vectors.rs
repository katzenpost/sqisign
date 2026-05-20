//! Differential test of `jac_to_ws`.
mod common;

use sqisign_ec::{jac_to_ws, JacPoint};
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_to_ws.json"
);

#[test]
fn jac_to_ws_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_to_ws");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp_q = common::jac_point_from("c", &v.outputs);
        let exp_t = common::fp2_from("t", &v.outputs);
        let exp_ao3 = common::fp2_from("ao3", &v.outputs);
        let mut q = JacPoint::zero();
        let mut t = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        let mut ao3 = t;
        jac_to_ws(&mut q, &mut t, &mut ao3, &p, &curve);
        assert_eq!(q, exp_q, "vector {} Q diverged", v.id);
        assert_eq!(t, exp_t, "vector {} t diverged", v.id);
        assert_eq!(ao3, exp_ao3, "vector {} ao3 diverged", v.id);
    }
}
