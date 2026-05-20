//! Loader: parses the canonical-bytes JSON vectors committed under
//! `vectors/precomp/` into the typed precomputed constants exported by
//! the crate root.
//!
//! Each JSON is embedded at compile time via `include_str!` and parsed
//! once per process via the `LazyLock` machinery in `lib.rs`; no
//! filesystem access happens at runtime. The differential-test side of
//! the integrity gate lives in `tests/precomp_vectors.rs`.

use std::collections::BTreeMap;

use sqisign_ec::{EcBasis, EcCurve, EcPoint};
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_quaternion::dim2::IbzMat2x2;
use sqisign_quaternion::dim4::{ibz_mat_4x4_new, IbzMat4x4};
use sqisign_quaternion::{Ibz, QuatAlgElem, QuatLattice, QuatLeftIdeal};

use crate::{CurveWithEndomorphismRing, QuatPExtremalMaximalOrder};

const TWO_TO_SECURITY_BITS_JSON: &str =
    include_str!("../../../vectors/precomp/TWO_TO_SECURITY_BITS.json");
const TORSION_PLUS_2POWER_JSON: &str =
    include_str!("../../../vectors/precomp/TORSION_PLUS_2POWER.json");
const SEC_DEGREE_JSON: &str = include_str!("../../../vectors/precomp/SEC_DEGREE.json");
const COM_DEGREE_JSON: &str = include_str!("../../../vectors/precomp/COM_DEGREE.json");
const QUAT_PRIME_COFACTOR_JSON: &str =
    include_str!("../../../vectors/precomp/QUAT_prime_cofactor.json");
const QUATALG_PINFTY_JSON: &str = include_str!("../../../vectors/precomp/QUATALG_PINFTY.json");
const EXTREMAL_ORDERS_JSON: &str = include_str!("../../../vectors/precomp/EXTREMAL_ORDERS.json");
const CONNECTING_IDEALS_JSON: &str =
    include_str!("../../../vectors/precomp/CONNECTING_IDEALS.json");
const CONJUGATING_ELEMENTS_JSON: &str =
    include_str!("../../../vectors/precomp/CONJUGATING_ELEMENTS.json");
const CURVES_WITH_ENDOMORPHISMS_JSON: &str =
    include_str!("../../../vectors/precomp/CURVES_WITH_ENDOMORPHISMS.json");

/// One record's worth of recorded fields. We keep the canonical-bytes
/// hex strings here for direct lookup; conversion to `Ibz` or limb
/// arrays happens at the call site.
struct Record {
    outputs: BTreeMap<String, String>,
}

fn parse_records(json: &str, boundary_expected: &str, expected_count: usize) -> Vec<Record> {
    let value: serde_json::Value =
        serde_json::from_str(json).expect("precomp JSON must parse at compile-time-known content");
    let boundary = value["boundary"]
        .as_str()
        .expect("precomp JSON: missing boundary string");
    assert_eq!(
        boundary, boundary_expected,
        "precomp JSON boundary mismatch: expected {boundary_expected}, found {boundary}"
    );
    let arr = value["vectors"]
        .as_array()
        .expect("precomp JSON: vectors must be an array");
    assert_eq!(
        arr.len(),
        expected_count,
        "precomp JSON {boundary}: expected {expected_count} records, found {}",
        arr.len()
    );
    let mut out = Vec::with_capacity(arr.len());
    for (i, v) in arr.iter().enumerate() {
        let id = v["id"].as_u64().expect("precomp JSON: missing id");
        assert_eq!(
            id as usize, i,
            "precomp JSON {boundary}: non-monotonic ids (got {id} at index {i})"
        );
        let outputs = v["outputs"]
            .as_object()
            .expect("precomp JSON: missing outputs");
        let mut m = BTreeMap::new();
        for (k, val) in outputs {
            m.insert(
                k.clone(),
                val.as_str()
                    .expect("precomp JSON: output value must be a hex string")
                    .to_string(),
            );
        }
        out.push(Record { outputs: m });
    }
    out
}

fn decode_hex(s: &str) -> Vec<u8> {
    let body = s.strip_prefix("0x").expect("hex must be 0x-prefixed");
    hex::decode(body).expect("hex decode")
}

fn ibz_from_field(rec: &Record, name: &str) -> Ibz {
    let hex = rec
        .outputs
        .get(name)
        .unwrap_or_else(|| panic!("precomp: missing output {name}"));
    Ibz::from_canonical_bytes(&decode_hex(hex))
        .unwrap_or_else(|e| panic!("precomp: ibz canonical decode failed for {name}: {e}"))
}

fn u32_from_field(rec: &Record, name: &str) -> u32 {
    let bytes = decode_hex(
        rec.outputs
            .get(name)
            .unwrap_or_else(|| panic!("precomp: missing output {name}")),
    );
    assert_eq!(
        bytes.len(),
        4,
        "precomp: u32 field {name} must be 4 bytes, got {}",
        bytes.len()
    );
    u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

fn fp_from_field(rec: &Record, name: &str) -> [u64; NWORDS_FIELD] {
    let bytes = decode_hex(
        rec.outputs
            .get(name)
            .unwrap_or_else(|| panic!("precomp: missing output {name}")),
    );
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "precomp: fp_t field {name} must be {} bytes, got {}",
        NWORDS_FIELD * 8,
        bytes.len()
    );
    let mut out = [0u64; NWORDS_FIELD];
    for i in 0..NWORDS_FIELD {
        let mut a = [0u8; 8];
        a.copy_from_slice(&bytes[8 * i..8 * (i + 1)]);
        out[i] = u64::from_le_bytes(a);
    }
    out
}

fn fp2_from_fields(rec: &Record, prefix: &str) -> Fp2 {
    Fp2 {
        re: fp_from_field(rec, &format!("{prefix}_re")),
        im: fp_from_field(rec, &format!("{prefix}_im")),
    }
}

fn ec_point_from_fields(rec: &Record, prefix: &str) -> EcPoint {
    EcPoint {
        x: fp2_from_fields(rec, &format!("{prefix}_x")),
        z: fp2_from_fields(rec, &format!("{prefix}_z")),
    }
}

fn ec_basis_from_fields(rec: &Record, prefix: &str) -> EcBasis {
    EcBasis {
        P: ec_point_from_fields(rec, &format!("{prefix}_P")),
        Q: ec_point_from_fields(rec, &format!("{prefix}_Q")),
        PmQ: ec_point_from_fields(rec, &format!("{prefix}_PmQ")),
    }
}

fn ec_curve_from_fields(rec: &Record, prefix: &str) -> EcCurve {
    let is_a24 = u32_from_field(rec, &format!("{prefix}_is_A24"));
    EcCurve {
        A: fp2_from_fields(rec, &format!("{prefix}_A")),
        C: fp2_from_fields(rec, &format!("{prefix}_C")),
        A24: ec_point_from_fields(rec, &format!("{prefix}_A24")),
        is_A24_computed_and_normalized: is_a24 != 0,
    }
}

fn ibz_mat_2x2_from_fields(rec: &Record, prefix: &str) -> IbzMat2x2 {
    [
        [
            ibz_from_field(rec, &format!("{prefix}_00")),
            ibz_from_field(rec, &format!("{prefix}_01")),
        ],
        [
            ibz_from_field(rec, &format!("{prefix}_10")),
            ibz_from_field(rec, &format!("{prefix}_11")),
        ],
    ]
}

fn ibz_mat_4x4_from_fields(rec: &Record, prefix: &str) -> IbzMat4x4 {
    let mut out = ibz_mat_4x4_new();
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = ibz_from_field(rec, &format!("{prefix}_{i}_{j}"));
        }
    }
    out
}

fn quat_lattice_from_fields(rec: &Record, prefix: &str) -> QuatLattice {
    QuatLattice {
        denom: ibz_from_field(rec, &format!("{prefix}_denom")),
        basis: ibz_mat_4x4_from_fields(rec, &format!("{prefix}_basis")),
    }
}

fn quat_alg_elem_from_fields(rec: &Record, denom_field: &str, coord_prefix: &str) -> QuatAlgElem {
    QuatAlgElem {
        denom: ibz_from_field(rec, denom_field),
        coord: [
            ibz_from_field(rec, &format!("{coord_prefix}_0")),
            ibz_from_field(rec, &format!("{coord_prefix}_1")),
            ibz_from_field(rec, &format!("{coord_prefix}_2")),
            ibz_from_field(rec, &format!("{coord_prefix}_3")),
        ],
    }
}

pub(crate) fn load_single_ibz(boundary_short: &str) -> Ibz {
    let json = match boundary_short {
        "TWO_TO_SECURITY_BITS" => TWO_TO_SECURITY_BITS_JSON,
        "TORSION_PLUS_2POWER" => TORSION_PLUS_2POWER_JSON,
        "SEC_DEGREE" => SEC_DEGREE_JSON,
        "COM_DEGREE" => COM_DEGREE_JSON,
        "QUAT_prime_cofactor" => QUAT_PRIME_COFACTOR_JSON,
        other => panic!("load_single_ibz: unknown boundary {other}"),
    };
    let recs = parse_records(json, &format!("sqisign_precomp::{boundary_short}"), 1);
    ibz_from_field(&recs[0], "value")
}

pub(crate) fn load_named_ibz(
    boundary_short: &str,
    extra_records: &[&str],
    record_index: usize,
    field_name: &str,
) -> Ibz {
    let json = match boundary_short {
        "QUATALG_PINFTY" => QUATALG_PINFTY_JSON,
        other => panic!("load_named_ibz: unknown boundary {other}"),
    };
    let _ = extra_records;
    let recs = parse_records(json, &format!("sqisign_precomp::{boundary_short}"), 1);
    ibz_from_field(&recs[record_index], field_name)
}

pub(crate) fn load_extremal_orders() -> [QuatPExtremalMaximalOrder; 7] {
    let recs = parse_records(EXTREMAL_ORDERS_JSON, "sqisign_precomp::EXTREMAL_ORDERS", 7);
    // Build them one by one to avoid requiring Default on the array slot.
    let mut iter = recs.into_iter().map(|rec| QuatPExtremalMaximalOrder {
        order: quat_lattice_from_fields(&rec, "order"),
        z: quat_alg_elem_from_fields(&rec, "z_denom", "z"),
        t: quat_alg_elem_from_fields(&rec, "t_denom", "t"),
        q: u32_from_field(&rec, "q"),
    });
    [
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
    ]
}

pub(crate) fn load_connecting_ideals() -> [QuatLeftIdeal; 7] {
    let recs = parse_records(
        CONNECTING_IDEALS_JSON,
        "sqisign_precomp::CONNECTING_IDEALS",
        7,
    );
    let mut iter = recs.into_iter().map(|rec| QuatLeftIdeal {
        norm: ibz_from_field(&rec, "ideal_norm"),
        lattice: quat_lattice_from_fields(&rec, "ideal_lat"),
    });
    [
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
    ]
}

pub(crate) fn load_conjugating_elements() -> [QuatAlgElem; 7] {
    let recs = parse_records(
        CONJUGATING_ELEMENTS_JSON,
        "sqisign_precomp::CONJUGATING_ELEMENTS",
        7,
    );
    let mut iter = recs
        .into_iter()
        .map(|rec| quat_alg_elem_from_fields(&rec, "denom", "c"));
    [
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
    ]
}

pub(crate) fn load_curves_with_endomorphisms() -> [CurveWithEndomorphismRing; 7] {
    let recs = parse_records(
        CURVES_WITH_ENDOMORPHISMS_JSON,
        "sqisign_precomp::CURVES_WITH_ENDOMORPHISMS",
        7,
    );
    let mut iter = recs.into_iter().map(|rec| CurveWithEndomorphismRing {
        curve: ec_curve_from_fields(&rec, "curve"),
        basis_even: ec_basis_from_fields(&rec, "basis"),
        action_i: ibz_mat_2x2_from_fields(&rec, "action_i"),
        action_j: ibz_mat_2x2_from_fields(&rec, "action_j"),
        action_k: ibz_mat_2x2_from_fields(&rec, "action_k"),
        action_gen2: ibz_mat_2x2_from_fields(&rec, "action_gen2"),
        action_gen3: ibz_mat_2x2_from_fields(&rec, "action_gen3"),
        action_gen4: ibz_mat_2x2_from_fields(&rec, "action_gen4"),
    });
    [
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
        iter.next().unwrap(),
    ]
}
