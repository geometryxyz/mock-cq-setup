use crate::fk::UpperToeplitz;
use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use std::ops::Mul;

pub fn compute_qs<E: Pairing>(
    t: &DensePolynomial<E::ScalarField>,
    domain: &GeneralEvaluationDomain<E::ScalarField>,
    tau_powers: &[E::ScalarField],
) -> Vec<E::G1Affine> {
    /*
        - N (table size) is always pow2
        - Toeplitz multiplication will happen in 2 * N, so appending zero commitments on hs is not needed
    */

    let toeplitz = UpperToeplitz::from_poly(t);

    let mut tau_rev = tau_powers.to_vec();
    tau_rev.reverse();

    let hs: Vec<E::ScalarField> = toeplitz.mul_by_vec(&tau_powers);
    assert_eq!(hs.len(), 2 * domain.size());

    let ks: Vec<_> = domain.fft(&hs[..domain.size()]);

    let n_inv = domain.size_as_field_element().inverse().unwrap();
    let normalized_roots = domain.elements().map(|g_i| g_i * n_inv);

    let gen = E::G1::generator();
    let qs_at_tau: Vec<E::G1> = ks
        .iter()
        .zip(normalized_roots)
        .map(|(&ki, normalizer_i)| gen.mul(ki * normalizer_i))
        .collect();

    E::G1::normalize_batch(&qs_at_tau)
}
