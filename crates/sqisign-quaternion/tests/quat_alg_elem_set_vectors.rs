//! Differential test of `quat_alg_elem_set`.
use sqisign_quaternion::{quat_alg_elem_set, Ibz, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_elem_set.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_i32(l: &str, h: &str) -> i32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

fn read_elem(prefix: &str, m: &std::collections::BTreeMap<String, String>) -> QuatAlgElem {
    QuatAlgElem {
        denom: read_ibz("denom", &m[&format!("{prefix}_denom")]),
        coord: [
            read_ibz("c0", &m[&format!("{prefix}_c0")]),
            read_ibz("c1", &m[&format!("{prefix}_c1")]),
            read_ibz("c2", &m[&format!("{prefix}_c2")]),
            read_ibz("c3", &m[&format!("{prefix}_c3")]),
        ],
    }
}

#[test]
fn quat_alg_elem_set_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_elem_set");
    assert!(f.vectors.len() >= 80);
    for v in &f.vectors {
        let denom = read_le_i32("denom", &v.inputs["denom"]);
        let coords_bytes = decode("coords", &v.inputs["coords"]).unwrap();
        assert_eq!(coords_bytes.len(), 16);
        let cs: [i32; 4] = [
            i32::from_le_bytes(coords_bytes[0..4].try_into().unwrap()),
            i32::from_le_bytes(coords_bytes[4..8].try_into().unwrap()),
            i32::from_le_bytes(coords_bytes[8..12].try_into().unwrap()),
            i32::from_le_bytes(coords_bytes[12..16].try_into().unwrap()),
        ];
        let exp = read_elem("r", &v.outputs);
        let mut r = QuatAlgElem::new();
        quat_alg_elem_set(&mut r, denom, cs[0], cs[1], cs[2], cs[3]);
        assert_eq!(r.denom.0, exp.denom.0, "vector {}: denom", v.id);
        for k in 0..4 {
            assert_eq!(r.coord[k].0, exp.coord[k].0, "vector {}: coord {}", v.id, k);
        }
    }
}
