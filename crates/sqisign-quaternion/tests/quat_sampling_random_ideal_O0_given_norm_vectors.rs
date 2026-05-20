//! Differential test of `quat_sampling_random_ideal_O0_given_norm` against
//! C-derived vectors.
//!
//! Each vector seeds a fresh KAT-only `CtrDrbg` with a recorded 48-byte
//! entropy block (no personalization), constructs the standard extremal
//! maximal order via `quat_lattice_O0_set_extremal`, builds a `QuatAlg`
//! from the lvl1 prime in `vectors/precomp/QUATALG_PINFTY.json`, and
//! calls `quat_sampling_random_ideal_O0_given_norm(&mut drbg, &mut lideal,
//! &norm, is_prime, &params, prime_cofactor)`. The recorded `ok` flag, the
//! norm of the produced left ideal, and the canonical bytes of the
//! ideal's lattice (denominator plus 16 basis entries) must all match
//! the C reference byte-for-byte.
//!
//! `prime_cofactor` is conveyed by the C emitter as an `Ibz`; for the
//! prime path the emitter writes a zero, which the Rust test interprets
//! as `None`.

mod common;

use common::{ibz_eq, mat4x4_eq, read_ibz};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{
    ibz_is_zero, quat_lattice_O0_set_extremal, quat_sampling_random_ideal_O0_given_norm, Ibz,
    QuatAlg, QuatLeftIdeal, QuatPExtremalMaximalOrder, QuatRepresentIntegerParams,
};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_sampling_random_ideal_O0_given_norm.json"
);

const QUATALG_PINFTY_VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/precomp/QUATALG_PINFTY.json"
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

fn load_pinfty_prime() -> Ibz {
    let f = load(QUATALG_PINFTY_VECTORS).expect("QUATALG_PINFTY.json load");
    let r = f
        .vectors
        .first()
        .expect("QUATALG_PINFTY.json: at least one record");
    read_ibz("p", &r.outputs)
}

#[test]
#[allow(non_snake_case)]
fn quat_sampling_random_ideal_O0_given_norm_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::quat_sampling_random_ideal_O0_given_norm"
    );
    assert!(!f.vectors.is_empty(), "vector battery empty");

    let p = load_pinfty_prime();
    let alg = QuatAlg::init_set(&p);
    let mut order = QuatPExtremalMaximalOrder::new();
    quat_lattice_O0_set_extremal(&mut order);

    let params = QuatRepresentIntegerParams {
        primality_test_iterations: 25,
        order: &order,
        algebra: &alg,
    };

    for v in &f.vectors {
        let entropy = fixed48(&decode("entropy", &v.inputs["entropy"]).expect("entropy hex"));
        let norm = read_ibz("norm", &v.inputs);
        let is_prime = read_i32_le("is_prime", &v.inputs);
        let cof_raw = read_ibz("prime_cofactor", &v.inputs);
        let cof_ref = if ibz_is_zero(&cof_raw) != 0 {
            None
        } else {
            Some(&cof_raw)
        };

        let ok_exp = read_ok("ok", &v.outputs);
        let l_norm_exp = read_ibz("lideal_norm", &v.outputs);
        let lat_denom_exp = read_ibz("lideal_lat_denom", &v.outputs);
        let lat_basis_exp = common::read_mat4x4("lideal_lat_basis", &v.outputs);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut lideal = QuatLeftIdeal::new();
        let ok_got = quat_sampling_random_ideal_O0_given_norm(
            &mut drbg, &mut lideal, &norm, is_prime, &params, cof_ref,
        );

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        assert!(
            ibz_eq(&lideal.norm, &l_norm_exp),
            "vector {}: lideal.norm",
            v.id
        );
        assert!(
            ibz_eq(&lideal.lattice.denom, &lat_denom_exp),
            "vector {}: lideal.lattice.denom",
            v.id
        );
        assert!(
            mat4x4_eq(&lideal.lattice.basis, &lat_basis_exp),
            "vector {}: lideal.lattice.basis",
            v.id
        );
    }
}
