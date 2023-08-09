mod common;
mod fk;
mod lagrange;
mod powers;
mod quotients;
mod utils;

pub use lagrange::lagrange_commitments;
pub use powers::compute_g_powers;
pub use utils::deserialize_points;

// pub struct CommonPreprocessedInput<E: PairingEngine> {
//     pub(crate) zv_2: E::G2Affine,
//     pub(crate) t_2: E::G2Affine,
// }

// pub struct Index<E: PairingEngine> {
//     pub(crate) common: CommonPreprocessedInput<E>,
//     pub(crate) qs: Vec<E::G1Affine>,
//     pub(crate) ls: Vec<E::G1Affine>,
//     pub(crate) ls_at_0: Vec<E::G1Affine>,
// }

/*

pub struct StaticTableValues<E: MultiMillerLoop> {
    size: usize,
    // Mapping from value to its index in the table
    value_index_mapping: BTreeMap<E::Scalar, usize>,
    // quotient commitments
    qs: Vec<E::G1>,
}
*/

#[cfg(feature = "serialize")]
#[test]
fn srs_roundtrip() {
    use ark_bn254::{Fr, G1Projective, G2Projective};
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
    use ark_bn254::{Fr, G1Projective};
    use lagrange::lagrange_commitments;
    let k = 3;
    let n = 1 << k;
    let tau = Fr::from(100 as u64);

    let l_basis = lagrange_commitments::<G1Projective>(tau, n, Some("lagrange.bin"));
}
