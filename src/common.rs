use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, One};
use std::ops::Mul;

pub struct CommonPreprocessedInput<E: Pairing> {
    pub(crate) zv_2: E::G2Affine,
    pub(crate) t_2: E::G2Affine,
    pub(crate) x_b0_bound: E::G2Affine,
    pub(crate) srs_g1_len: usize,
}

impl<E: Pairing> CommonPreprocessedInput<E> {
    fn compute(
        tau: E::ScalarField,
        table_coeffs: &[E::ScalarField],
        srs_g2_len: usize,
        srs_g1_len: usize,
        circuit_domain: usize,
    ) -> Self {
        let g2 = E::G2::generator();

        // zv_2 = x^n - 1
        // We always use the max power of g2
        let zv = tau.pow(&[(srs_g2_len - 1) as u64]) - E::ScalarField::one();
        let zv_2: E::G2Affine = g2.mul(zv).into();

        // TODO: compute powers of tau once and serialize that array
        let powers_of_tau: Vec<E::ScalarField> =
            std::iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(srs_g1_len)
                .collect();
        assert_eq!(powers_of_tau.len(), table_coeffs.len());
        let table_at_tau: E::ScalarField = table_coeffs
            .iter()
            .zip(powers_of_tau.iter())
            .map(|(&t_i, tau_pow_i)| t_i * tau_pow_i)
            .sum();

        let t_2: E::G2Affine = g2.mul(table_at_tau).into();

        let b0_bound_index = srs_g1_len - 1 - (circuit_domain - 2);
        let x_b0_bound: E::G2Affine = g2.mul(tau.pow(&[b0_bound_index as u64])).into();

        Self {
            zv_2,
            t_2,
            x_b0_bound,
            srs_g1_len,
        }
    }
}
