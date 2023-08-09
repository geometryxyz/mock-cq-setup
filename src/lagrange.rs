#[cfg(feature = "parallel")]
use crate::utils::parallelize;
use ark_ec::CurveGroup;
use ark_ff::{FftField, Field, One};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

pub fn lagrange_commitments<G: CurveGroup>(tau: G::ScalarField, n: u64) -> Vec<G::Affine> {
    let mut g_lagrange_projective = vec![G::zero(); n as usize];
    let w = G::ScalarField::get_root_of_unity(n).unwrap();

    /*
               w^i * zh(X)
        Li =  -------------
              N * (X - w^i)
    */

    let gen = G::generator();
    let zh = tau.pow(&[n]) - G::ScalarField::one();
    let const_multiplier = zh * (G::ScalarField::from(n).inverse().unwrap());

    #[cfg(not(feature = "parallel"))]
    for (i, li) in g_lagrange_projective.iter_mut().enumerate() {
        let w_pow_i = w.pow(&[i as u64]);
        *li = gen.mul(w_pow_i * const_multiplier * (tau - w_pow_i).inverse().unwrap());
    }

    #[cfg(feature = "parallel")]
    parallelize(&mut g_lagrange_projective, |g, start| {
        for (idx, g) in g.iter_mut().enumerate() {
            let offset = start + idx;
            let w_pow_i = w.pow(&[offset as u64]);
            *g = gen.mul(w_pow_i * const_multiplier * (tau - w_pow_i).inverse().unwrap());
        }
    });

    G::normalize_batch(&g_lagrange_projective)
}

pub fn lagrange_openings_commitments_at_zero<G: CurveGroup>(
    tau: G::ScalarField,
    n: usize,
) -> Vec<G::Affine> {
    assert!(crate::utils::is_pow_2(n));

    /*
               w^i * zh(X)
        Li =  -------------
              N * (X - w^i)

        Li(0) =   -w^i
                ---------- = 1/N
                N * -w^i
    */

    let domain = GeneralEvaluationDomain::<G::ScalarField>::new(n).unwrap();

    let lagrange_at_tau = domain.evaluate_all_lagrange_coefficients(tau);

    let li_at_zero = G::ScalarField::from(n as u64).inverse().unwrap();
    let w: G::ScalarField = domain.element(1);
    let gen = G::generator();
    // X^{N-1} evaluated at tau
    let x_to_n_minus_one = tau.pow(&[n as u64 - 1]);

    let mut lagrange_openings_at_zero = vec![G::zero(); n as usize];

    // w^(N - i)*L_i(tau) - 1/N * X^(N-1)
    #[cfg(not(feature = "parallel"))]
    for (i, li) in lagrange_openings_at_zero.iter_mut().enumerate() {
        let w_inv_pow_i = w.pow(&[(n - i) as u64]);
        *li = gen.mul(w_inv_pow_i * lagrange_at_tau[i] - li_at_zero * x_to_n_minus_one);
    }

    #[cfg(feature = "parallel")]
    parallelize(&mut lagrange_openings_at_zero, |g, start| {
        for (idx, g) in g.iter_mut().enumerate() {
            let offset = start + idx;
            let w_inv_pow_i = w.pow(&[(n - offset) as u64]);
            *g = gen.mul(w_inv_pow_i * lagrange_at_tau[idx] - li_at_zero * x_to_n_minus_one);
        }
    });

    G::normalize_batch(&lagrange_openings_at_zero)
}

#[cfg(test)]
mod test_lagrange {
    use std::ops::Mul;

    use ark_bn254::{Fr, G1Affine, G1Projective};
    use ark_ec::{CurveGroup, Group};
    use ark_ff::{FftField, Field, One, Zero};
    use ark_poly::{
        univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain,
        Polynomial,
    };

    // given x coords construct Li polynomials
    pub fn construct_lagrange_basis<F: FftField>(
        evaluation_domain: &[F],
    ) -> Vec<DensePolynomial<F>> {
        let mut bases = Vec::with_capacity(evaluation_domain.len());
        for i in 0..evaluation_domain.len() {
            let mut l_i = DensePolynomial::from_coefficients_slice(&[F::one()]);
            let x_i = evaluation_domain[i];
            for (j, _) in evaluation_domain.iter().enumerate() {
                if j != i {
                    let xi_minus_xj_inv = (x_i - evaluation_domain[j]).inverse().unwrap();
                    l_i = &l_i
                        * &DensePolynomial::from_coefficients_slice(&[
                            -evaluation_domain[j] * xi_minus_xj_inv,
                            xi_minus_xj_inv,
                        ]);
                }
            }

            bases.push(l_i);
        }

        bases
    }

    // cargo test --features=parallel compute_lagrange_commitments
    #[test]
    fn compute_lagrange_commitments() {
        let k = 3;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let l_evals = domain.evaluate_all_lagrange_coefficients(tau);

        let g: ark_ec::short_weierstrass::Projective<ark_bn254::g1::Config> =
            G1Projective::generator();
        let l_coms: Vec<G1Projective> = l_evals.iter().map(|li| g.mul(li)).collect();
        let l_coms: Vec<G1Affine> = G1Projective::normalize_batch(&l_coms);

        let l_basis_coms = super::lagrange_commitments::<G1Projective>(tau, n as u64);
        assert_eq!(l_basis_coms, l_coms);
    }

    // cargo test --features=parallel compute_lagrange_opening_commitments_at_zero
    #[test]
    fn compute_lagrange_opening_commitments_at_zero() {
        use ark_bn254::{Fr, G1Projective};
        let k = 3;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let roots: Vec<_> = domain.elements().collect();
        let l_basis = construct_lagrange_basis::<Fr>(&roots);

        let l_at_zero = Fr::from(n as u64).inverse().unwrap();
        let x_poly = DensePolynomial::from_coefficients_slice(&[Fr::zero(), Fr::one()]);

        let g = G1Projective::generator();
        let q_commitments: Vec<G1Projective> = l_basis
            .iter()
            .map(|li| {
                let mut li_x_minus_li_0 = li.clone();
                li_x_minus_li_0[0] -= l_at_zero;

                let qi = &li_x_minus_li_0 / &x_poly;
                g.mul(qi.evaluate(&tau))
            })
            .collect();

        let lagrange_openings_commitments =
            super::lagrange_openings_commitments_at_zero::<G1Projective>(tau, n);
        assert_eq!(lagrange_openings_commitments, q_commitments);
    }
}
