use ark_ec::CurveGroup;
use ark_ff::{Field, FftField, One};
#[cfg(feature = "parallel")]
use crate::utils::parallelize;
#[cfg(feature = "serialize")]
use crate::utils::{write_points, serialize_points};

pub fn lagrange_commitments<G: CurveGroup>(tau: G::ScalarField, n: u64, path: Option<&str>) -> Option<Vec<G::Affine>> {
    let mut g_lagrange_projective = vec![G::zero(); n as usize];
    let w = G::ScalarField::get_root_of_unity(n).unwrap();

    /*
        Li =  w^i * zh(X)
              ------------
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

#[cfg(not(feature = "serialize"))]
#[cfg(test)]
mod test_lagrange {
    use std::ops::Mul;

    use ark_bn254::G1Affine;
    use ark_ec::{Group, CurveGroup};
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};

    // cargo test --all --features=parallel -- --nocapture
    #[test]
    fn compute_lagrange_commitments() {
        use ark_bn254::{G1Projective, Fr};
        let k = 3; 
        let n = 1 << k; 
        let tau = Fr::from(100 as u64);

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let l_evals = domain.evaluate_all_lagrange_coefficients(tau);

        let g: ark_ec::short_weierstrass::Projective<ark_bn254::g1::Config> = G1Projective::generator();
        let l_coms: Vec<G1Projective> = l_evals.iter().map(|li| g.mul(li)).collect();
        let l_coms: Vec<G1Affine> = G1Projective::normalize_batch(&l_coms);

        let l_basis_coms = super::lagrange_commitments::<G1Projective>(tau, n as u64, None);
        assert_eq!(l_basis_coms.unwrap(), l_coms);
    }
}