use crate::fk::UpperToeplitz;
use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use std::ops::Mul;

pub fn compute_qs<E: Pairing>(
    t: &DensePolynomial<E::ScalarField>,
    domain: &GeneralEvaluationDomain<E::ScalarField>,
    tau_powers: &[E::ScalarField],
    path: Option<&str>
) -> Option<Vec<E::G1Affine>> {
    /*
        - N (table size) is always pow2
        - Toeplitz multiplication will happen in 2 * N, so appending zero commitments on hs is not needed
    */

    let toeplitz = UpperToeplitz::from_poly(t);

    let mut tau_rev = tau_powers.to_vec();
    tau_rev.reverse();

    let hs: Vec<E::ScalarField> = toeplitz.mul_by_vec(&tau_rev);
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

    #[cfg(feature = "serialize")]
    {
        let path = path.expect("Path not provided");
        let qs_affine = G::normalize_batch(&qs_at_tau);
        let data = serialize_points(&qs_affine);
        write_points(path, &data);
        None
    }

    #[cfg(not(feature = "serialize"))]
    Some(E::G1::normalize_batch(&qs_at_tau))
}

#[cfg(not(feature = "serialize"))]
#[cfg(test)]
mod powers_test {
    use ark_ec::Group;
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_ff::{Field, One};
    use ark_bn254::{Fr, G1Projective, Bn254};
    use ark_std::{test_rng, UniformRand};
    use std::ops::Mul;

    // cargo test test_qs
    #[test]
    fn test_qs() {
        let k = 5;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let n_inv = domain.size_as_field_element().inverse().unwrap();
        let normalized_roots: Vec<_> = domain.elements().map(|g_i| g_i * n_inv).collect();

        let mut rng = test_rng();

        let t_evals: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let t_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&t_evals));

        let g: G1Projective = G1Projective::generator();
        let q_commitments: Vec<G1Projective> = (0..n).map(|i| {
            let x_minus_w_i = DensePolynomial::from_coefficients_slice(&[-domain.element(i), Fr::one()]);

            let mut t_poly_i = t_poly.clone();
            t_poly_i[0] -= t_evals[i];

            let q = &t_poly_i / &x_minus_w_i;
            assert_eq!(t_poly_i, &q * &x_minus_w_i);

            g.mul(q.evaluate(&tau) * normalized_roots[i])
        }).collect();

        let tau_powers: Vec<Fr> =
        std::iter::successors(Some(Fr::one()), |p| Some(*p * tau))
            .take(n)
            .collect();

        let qs = super::compute_qs::<Bn254>(&t_poly, &domain, &tau_powers, None);
        assert_eq!(qs.unwrap(), q_commitments);
    }
}