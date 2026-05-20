//! Differential test of `DBLW` (modified Jacobian Weierstrass doubling).
mod common;

use sqisign_ec::{dblw, JacPoint};
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/DBLW.json");

#[test]
fn dblw_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::DBLW");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let t = common::fp2_from("t", &v.inputs);
        let exp_q = common::jac_point_from("c", &v.outputs);
        let exp_u = common::fp2_from("u", &v.outputs);
        let mut q = JacPoint::zero();
        let mut u = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        dblw(&mut q, &mut u, &p, &t);
        assert_eq!(q, exp_q, "vector {} Q diverged", v.id);
        assert_eq!(u, exp_u, "vector {} u diverged", v.id);
    }
}
