//! Differential test of `quat_represent_integer` against C-derived vectors.
//!
//! Each vector seeds a fresh KAT-only `CtrDrbg` with a recorded 48-byte
//! entropy block (no personalization), constructs the standard extremal
//! maximal order `MAXORD_O0` via `quat_lattice_O0_set_extremal`, builds a
//! `QuatAlg` from the lvl1 prime `QUATALG_PINFTY.p`, and calls
//! `quat_represent_integer(&mut drbg, &mut gamma, &n_gamma, non_diag, &params)`.
//! The recorded `ok` flag and the four canonical-bytes coordinates of
//! `gamma` (plus its denominator) must match the C reference byte-for-byte.

mod common;

use common::read_ibz;
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{
    quat_lattice_O0_set_extremal, quat_represent_integer, Ibz, QuatAlg, QuatAlgElem,
    QuatPExtremalMaximalOrder, QuatRepresentIntegerParams,
};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_represent_integer.json"
);

const QUATALG_PINFTY_VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/precomp/QUATALG_PINFTY.json"
);

/// Load the lvl1 quaternion-algebra ramification prime directly from the
/// committed precomp JSON. We avoid pulling in `sqisign-precomp` as a
/// dev-dependency because that crate already depends on `sqisign-quaternion`
/// (introducing a cycle); the JSON file is the single source of truth on
/// both sides anyway.
fn load_pinfty_prime() -> Ibz {
    let f = load(QUATALG_PINFTY_VECTORS).expect("QUATALG_PINFTY.json load");
    let r = f
        .vectors
        .first()
        .expect("QUATALG_PINFTY.json: at least one record");
    read_ibz("p", &r.outputs)
}

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

fn ibz_eq(a: &Ibz, b: &Ibz) -> bool {
    a.0 == b.0
}

#[test]
fn quat_represent_integer_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_quaternion::quat_represent_integer");
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
        let n_gamma = read_ibz("n_gamma", &v.inputs);
        let non_diag = read_i32_le("non_diag", &v.inputs);

        let ok_exp = read_ok("ok", &v.outputs);
        let gamma_denom_exp = read_ibz("gamma_denom", &v.outputs);
        let gamma_c0_exp = read_ibz("gamma_c0", &v.outputs);
        let gamma_c1_exp = read_ibz("gamma_c1", &v.outputs);
        let gamma_c2_exp = read_ibz("gamma_c2", &v.outputs);
        let gamma_c3_exp = read_ibz("gamma_c3", &v.outputs);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut gamma = QuatAlgElem::new();
        let ok_got = quat_represent_integer(&mut drbg, &mut gamma, &n_gamma, non_diag, &params);

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        if ok_exp != 0 {
            assert!(
                ibz_eq(&gamma.denom, &gamma_denom_exp),
                "vector {}: gamma.denom",
                v.id
            );
            assert!(
                ibz_eq(&gamma.coord[0], &gamma_c0_exp),
                "vector {}: gamma.coord[0]",
                v.id
            );
            assert!(
                ibz_eq(&gamma.coord[1], &gamma_c1_exp),
                "vector {}: gamma.coord[1]",
                v.id
            );
            assert!(
                ibz_eq(&gamma.coord[2], &gamma_c2_exp),
                "vector {}: gamma.coord[2]",
                v.id
            );
            assert!(
                ibz_eq(&gamma.coord[3], &gamma_c3_exp),
                "vector {}: gamma.coord[3]",
                v.id
            );
        }
    }
}
