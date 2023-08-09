use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, One};
use ark_std::cfg_iter;
use std::ops::Mul;

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::utils::is_pow_2;

pub struct CommonPreprocessedInput<E: Pairing> {
    pub(crate) zv_2: E::G2Affine,
    pub(crate) t_2: E::G2Affine,
    pub(crate) x_b0_bound: E::G2Affine,
    pub(crate) srs_g1_len: usize,
}

impl<E: Pairing> CommonPreprocessedInput<E> {
    fn compute(
        powers_of_tau: &[E::ScalarField],
        table_coeffs: &[E::ScalarField],
        srs_g1_len: usize,
        circuit_domain: usize,
    ) -> Self {
        assert_eq!(powers_of_tau.len(), table_coeffs.len());
        assert!(is_pow_2(srs_g1_len));

        let tau = powers_of_tau[1];
        let g2 = E::G2::generator();

        // zv_2 = x^n - 1
        // force that g2_srs is always 1 more longer than g1_srs
        let zv = tau.pow(&[(srs_g1_len + 1 - 1) as u64]) - E::ScalarField::one();
        let zv_2: E::G2Affine = g2.mul(zv).into();

        let table_at_tau: E::ScalarField = cfg_iter!(table_coeffs)
            .zip(cfg_iter!(powers_of_tau))
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
