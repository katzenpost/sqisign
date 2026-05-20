//! Shared decode helpers for the hd differential vector tests.
//!
//! The hd record schema is built on top of the ec/gf schemas: every
//! `theta_point` field is four [`Fp2`]s; every `theta_couple_point` is
//! two [`EcPoint`]s; every `theta_couple_curve` is two [`EcCurve`]s; and
//! [`ThetaStructure`] adds the precomputation flag plus eight precomputed
//! [`Fp2`] factors. The dump emits each leaf [`Fp2`] as one
//! `<prefix>_re` / `<prefix>_im` pair of `digits` fields per
//! `put_fp2_named`.

#![allow(dead_code)]

use std::collections::BTreeMap;

use sqisign_ec::{EcBasis, EcCurve, EcPoint, JacPoint};
use sqisign_gf::{Fp, Fp2, NWORDS_FIELD};
use sqisign_hd::{
    ThetaCoupleCurve, ThetaCoupleJacPoint, ThetaCouplePoint, ThetaKernelCouplePoints, ThetaPoint,
    ThetaStructure,
};
use sqisign_vectors::decode;

fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be exactly {NWORDS_FIELD} u64 limbs"
    );
    let mut limbs = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        limbs[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    limbs
}

pub fn fp2_from(prefix: &str, fields: &BTreeMap<String, String>) -> Fp2 {
    let re_key = format!("{prefix}_re");
    let im_key = format!("{prefix}_im");
    let re_bytes = decode(&re_key, &fields[&re_key]).expect("re hex");
    let im_bytes = decode(&im_key, &fields[&im_key]).expect("im hex");
    Fp2 {
        re: fp_from(&re_key, &re_bytes),
        im: fp_from(&im_key, &im_bytes),
    }
}

pub fn ec_point_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcPoint {
    EcPoint {
        x: fp2_from(&format!("{prefix}_x"), fields),
        z: fp2_from(&format!("{prefix}_z"), fields),
    }
}

pub fn jac_point_from(prefix: &str, fields: &BTreeMap<String, String>) -> JacPoint {
    JacPoint {
        x: fp2_from(&format!("{prefix}_x"), fields),
        y: fp2_from(&format!("{prefix}_y"), fields),
        z: fp2_from(&format!("{prefix}_z"), fields),
    }
}

pub fn ec_curve_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcCurve {
    let is_a24_key = format!("{prefix}_is_A24");
    let is_a24_bytes = decode(&is_a24_key, &fields[&is_a24_key]).expect("is_A24 hex");
    assert_eq!(is_a24_bytes.len(), 4);
    let is_a24 = u32::from_le_bytes(is_a24_bytes.as_slice().try_into().unwrap()) != 0;
    EcCurve {
        A: fp2_from(&format!("{prefix}_A"), fields),
        C: fp2_from(&format!("{prefix}_C"), fields),
        A24: ec_point_from(&format!("{prefix}_A24"), fields),
        is_A24_computed_and_normalized: is_a24,
    }
}

pub fn ec_basis_from(prefix: &str, fields: &BTreeMap<String, String>) -> EcBasis {
    EcBasis {
        P: ec_point_from(&format!("{prefix}_P"), fields),
        Q: ec_point_from(&format!("{prefix}_Q"), fields),
        PmQ: ec_point_from(&format!("{prefix}_PmQ"), fields),
    }
}

pub fn theta_point_from(prefix: &str, fields: &BTreeMap<String, String>) -> ThetaPoint {
    ThetaPoint {
        x: fp2_from(&format!("{prefix}_x"), fields),
        y: fp2_from(&format!("{prefix}_y"), fields),
        z: fp2_from(&format!("{prefix}_z"), fields),
        t: fp2_from(&format!("{prefix}_t"), fields),
    }
}

pub fn theta_couple_point_from(
    prefix: &str,
    fields: &BTreeMap<String, String>,
) -> ThetaCouplePoint {
    ThetaCouplePoint {
        p1: ec_point_from(&format!("{prefix}_P1"), fields),
        p2: ec_point_from(&format!("{prefix}_P2"), fields),
    }
}

pub fn theta_couple_jac_point_from(
    prefix: &str,
    fields: &BTreeMap<String, String>,
) -> ThetaCoupleJacPoint {
    ThetaCoupleJacPoint {
        p1: jac_point_from(&format!("{prefix}_P1"), fields),
        p2: jac_point_from(&format!("{prefix}_P2"), fields),
    }
}

pub fn theta_couple_curve_from(
    prefix: &str,
    fields: &BTreeMap<String, String>,
) -> ThetaCoupleCurve {
    ThetaCoupleCurve {
        e1: ec_curve_from(&format!("{prefix}_E1"), fields),
        e2: ec_curve_from(&format!("{prefix}_E2"), fields),
    }
}

pub fn theta_kernel_couple_points_from(
    prefix: &str,
    fields: &BTreeMap<String, String>,
) -> ThetaKernelCouplePoints {
    ThetaKernelCouplePoints {
        t1: theta_couple_point_from(&format!("{prefix}_T1"), fields),
        t2: theta_couple_point_from(&format!("{prefix}_T2"), fields),
        t1m2: theta_couple_point_from(&format!("{prefix}_T1m2"), fields),
    }
}

pub fn theta_structure_from(prefix: &str, fields: &BTreeMap<String, String>) -> ThetaStructure {
    let precomp_key = format!("{prefix}_precomp");
    let precomp_bytes = decode(&precomp_key, &fields[&precomp_key]).expect("precomp hex");
    assert_eq!(precomp_bytes.len(), 4);
    let precomp = u32::from_le_bytes(precomp_bytes.as_slice().try_into().unwrap()) != 0;
    ThetaStructure {
        null_point: theta_point_from(&format!("{prefix}_null"), fields),
        precomputation: precomp,
        xyz_big_0: fp2_from(&format!("{prefix}_XYZ0"), fields),
        yzt_big_0: fp2_from(&format!("{prefix}_YZT0"), fields),
        xzt_big_0: fp2_from(&format!("{prefix}_XZT0"), fields),
        xyt_big_0: fp2_from(&format!("{prefix}_XYT0"), fields),
        xyz_0: fp2_from(&format!("{prefix}_xyz0"), fields),
        yzt_0: fp2_from(&format!("{prefix}_yzt0"), fields),
        xzt_0: fp2_from(&format!("{prefix}_xzt0"), fields),
        xyt_0: fp2_from(&format!("{prefix}_xyt0"), fields),
    }
}

pub fn u32_field(key: &str, fields: &BTreeMap<String, String>) -> u32 {
    let bytes = decode(key, &fields[key]).expect("u32 hex");
    assert_eq!(bytes.len(), 4);
    u32::from_le_bytes(bytes.as_slice().try_into().unwrap())
}
