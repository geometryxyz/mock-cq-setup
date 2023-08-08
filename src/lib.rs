mod powers;
mod utils;
mod lagrange;

pub use lagrange::lagrange_commitments; 
pub use powers::compute_g_powers;
pub use utils::deserialize_points;

#[cfg(feature = "serialize")]
#[test]
fn srs_roundtrip() {
    use ark_bn254::{G1Projective, Fr, G2Projective};
    use lagrange::lagrange_commitments;
    use powers::compute_g_powers;

    let k = 10; 
    let n = 1 << k; 
    let tau = Fr::from(100 as u64);

    let _ = lagrange_commitments::<G1Projective>(tau, n, Some("lagrange.bin"));
    let _ = compute_g_powers::<G1Projective>(tau, n as usize, Some("srs_g1.bin"));
    let _ = compute_g_powers::<G2Projective>(tau, n as usize, Some("srs_g2.bin"));
}

#[cfg(feature = "serialize")]
#[test]
fn serialize_srs_test() {
    use ark_bn254::{G1Projective, Fr};
    use lagrange::lagrange_commitments;
    let k = 3; 
    let n = 1 << k; 
    let tau = Fr::from(100 as u64);

    let l_basis = lagrange_commitments::<G1Projective>(tau, n, Some("lagrange.bin"));
}

