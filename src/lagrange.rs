#[cfg(feature = "parallel")]
use crate::utils::parallelize;
#[cfg(feature = "serialize")]
use crate::utils::{serialize_points, write_points};
use ark_ec::CurveGroup;
use ark_ff::{FftField, Field, One, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

pub fn lagrange_commitments<G: CurveGroup>(
    tau: G::ScalarField,
    n: u64,
    path: Option<&str>,
) -> Option<Vec<G::Affine>> {
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

    #[cfg(feature = "serialize")]
    {
        let path = path.expect("Path not provided");
        let lagrange_affine = G::normalize_batch(&g_lagrange_projective);
        let data = serialize_points(&lagrange_affine);
        write_points(path, &data);
        None
    }

    #[cfg(not(feature = "serialize"))]
    Some(G::normalize_batch(&g_lagrange_projective))
}


pub fn lagrange_commitments_at_zero<G: CurveGroup>(
    tau: G::ScalarField,
    n: usize,
    path: Option<&str>,
) -> Option<Vec<G::Affine>> {
    fn is_pow_2(x: usize) -> bool {
        (x & (x - 1)) == 0
    }
    assert!(is_pow_2(n));

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

    let mut lagrange_at_zero = vec![G::zero(); n as usize];

    // w^(N - i)*L_i(tau) - 1/N * X^(N-1)
    #[cfg(not(feature = "parallel"))]
    for (i, li) in lagrange_at_zero.iter_mut().enumerate() {
        let w_inv_pow_i = w.pow(&[(n - i) as u64]);
        *li = gen.mul(w_inv_pow_i * lagrange_at_tau[i] - li_at_zero * x_to_n_minus_one);
    }

    #[cfg(feature = "parallel")]
    parallelize(&mut lagrange_at_zero, |g, start| {
        for (idx, g) in g.iter_mut().enumerate() {
            let offset = start + idx;
            let w_inv_pow_i = w.pow(&[(n - offset) as u64]);
            *g = gen.mul(w_inv_pow_i * lagrange_at_tau[idx] - li_at_zero * x_to_n_minus_one);
        }
    });

    #[cfg(feature = "serialize")]
    {
        let path = path.expect("Path not provided");
        let lagrange_at_zero_affine = G::normalize_batch(&lagrange_at_zero);
        let data = serialize_points(&lagrange_at_zero_affine);
        write_points(path, &data);
        None
    }

    #[cfg(not(feature = "serialize"))]
    Some(G::normalize_batch(&lagrange_at_zero))
}

#[cfg(not(feature = "serialize"))]
#[cfg(test)]
mod test_lagrange {
    use std::ops::Mul;

    use ark_bn254::{Fr, G1Affine};
    use ark_ec::{CurveGroup, Group};
    use ark_ff::{batch_inversion, Field};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

    // cargo test --features=parallel compute_lagrange_commitments
    #[test]
    fn compute_lagrange_commitments() {
        use ark_bn254::{Fr, G1Projective};
        let k = 3;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let l_evals = domain.evaluate_all_lagrange_coefficients(tau);

        let g: ark_ec::short_weierstrass::Projective<ark_bn254::g1::Config> =
            G1Projective::generator();
        let l_coms: Vec<G1Projective> = l_evals.iter().map(|li| g.mul(li)).collect();
        let l_coms: Vec<G1Affine> = G1Projective::normalize_batch(&l_coms);

        let l_basis_coms = super::lagrange_commitments::<G1Projective>(tau, n as u64, None);
        assert_eq!(l_basis_coms.unwrap(), l_coms);
    }

    #[test]
    fn test_root_inverses() {
        let n = 16;
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let w = domain.element(1);
        let mut roots_inv: Vec<Fr> = domain.elements().collect();
        batch_inversion(&mut roots_inv);

        for (i, &w_inv_pow_i) in roots_inv.iter().enumerate() {
            let w_inv_true = w.pow(&[(n - i) as u64]);

            assert_eq!(w_inv_true, w_inv_pow_i);
        }
    }
}
