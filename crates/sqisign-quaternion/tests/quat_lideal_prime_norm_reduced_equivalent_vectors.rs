//! Differential test of `quat_lideal_prime_norm_reduced_equivalent`
//! against C-derived vectors.
//!
//! Each vector seeds a fresh KAT-only [`CtrDrbg`] with a recorded 48-byte
//! entropy block (no personalization), reads the input generator `x` (a
//! `QuatAlgElem`) and the algebra's prime `p`, and reconstructs the
//! input left ideal via `quat_lideal_create_principal(&mut L, &x, &O0, &alg)`
//! (the same O0 and algebra the C harness used). Then calls
//! `quat_lideal_prime_norm_reduced_equivalent(&mut drbg, &mut L, &alg,
//! primality_num_iter, equiv_bound_coeff, &O0)`. The recorded `ok` flag,
//! the norm of the produced left ideal, and the canonical bytes of the
//! ideal's lattice (denominator + 16 basis entries) must all match the C
//! reference byte-for-byte.

mod common;

use common::{ibz_eq, mat4x4_eq, read_ibz};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{
    quat_lattice_O0_set, quat_lideal_create_principal,
    quat_lideal_prime_norm_reduced_equivalent, Ibz, QuatAlg, QuatAlgElem, QuatLattice,
    QuatLeftIdeal,
};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lideal_prime_norm_reduced_equivalent.json"
);

fn fixed48(v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "entropy must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

fn read_i32_le(name: &str, m: &std::collections::BTreeMap<String, String>) -> i32 {
    let b = decode(name, &m[name]).expect("hex decode");
    assert_eq!(b.len(), 4, "{name}: expected 4-byte i32");
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

fn read_ok(name: &str, m: &std::collections::BTreeMap<String, String>) -> i32 {
    let b = decode(name, &m[name]).expect("hex decode");
    assert_eq!(b.len(), 1, "{name}: expected 1-byte ok flag");
    b[0] as i8 as i32
}

fn read_quat_elem(prefix: &str, m: &std::collections::BTreeMap<String, String>) -> QuatAlgElem {
    let mut e = QuatAlgElem::new();
    e.denom = read_ibz(&format!("{prefix}_denom"), m);
    e.coord[0] = read_ibz(&format!("{prefix}_c0"), m);
    e.coord[1] = read_ibz(&format!("{prefix}_c1"), m);
    e.coord[2] = read_ibz(&format!("{prefix}_c2"), m);
    e.coord[3] = read_ibz(&format!("{prefix}_c3"), m);
    e
}

#[test]
fn quat_lideal_prime_norm_reduced_equivalent_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::quat_lideal_prime_norm_reduced_equivalent"
    );
    assert!(!f.vectors.is_empty(), "vector battery empty");

    for v in &f.vectors {
        let entropy = fixed48(&decode("entropy", &v.inputs["entropy"]).expect("entropy hex"));
        let x = read_quat_elem("x", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let primality_num_iter = read_i32_le("primality_num_iter", &v.inputs);
        let equiv_bound_coeff = read_i32_le("equiv_bound_coeff", &v.inputs);

        let ok_exp = read_ok("ok", &v.outputs);
        let l_norm_exp = read_ibz("lideal_norm", &v.outputs);
        let lat_denom_exp = read_ibz("lideal_lat_denom", &v.outputs);
        let lat_basis_exp = common::read_mat4x4("lideal_lat_basis", &v.outputs);

        // Rebuild O0 and the algebra from the recorded prime.
        let mut o0 = QuatLattice::new();
        quat_lattice_O0_set(&mut o0);
        let alg = QuatAlg::init_set(&p);

        // Rebuild the input left ideal via the same call the C harness used.
        let mut lideal = QuatLeftIdeal::new();
        quat_lideal_create_principal(&mut lideal, &x, &o0, &alg);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let ok_got = quat_lideal_prime_norm_reduced_equivalent(
            &mut drbg,
            &mut lideal,
            &alg,
            primality_num_iter,
            equiv_bound_coeff,
            &o0,
        );

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        assert!(
            ibz_eq(&lideal.norm, &l_norm_exp),
            "vector {}: lideal.norm diverged (got {}, want {})",
            v.id,
            lideal.norm.0,
            l_norm_exp.0
        );
        assert!(
            ibz_eq(&lideal.lattice.denom, &lat_denom_exp),
            "vector {}: lideal.lattice.denom diverged",
            v.id
        );
        assert!(
            mat4x4_eq(&lideal.lattice.basis, &lat_basis_exp),
            "vector {}: lideal.lattice.basis diverged",
            v.id
        );

        // Silence dead-code warnings for the local helper imports when
        // the test consumes only a subset of `Ibz`'s API.
        let _ = Ibz::zero();
    }
}
