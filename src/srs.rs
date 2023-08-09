#[cfg(feature = "parallel")]
use crate::utils::parallelize;
use ark_ec::CurveGroup;
#[cfg(feature = "parallel")]
use rayon::{self, prelude::*};

/// [tau^0]G, ..., [tau^{n-1}]G
pub fn compute_g_powers<G: CurveGroup>(
    tau: G::ScalarField,
    n: usize,
) -> Vec<G::Affine> {
    let mut g_srs = vec![G::zero(); n - 1];

    #[cfg(not(feature = "parallel"))]
    let g_srs: Vec<G> = std::iter::once(G::generator())
        .chain(g_srs.iter().scan(G::generator(), |state, _| {
            *state *= &tau;
            Some(*state)
        }))
        .collect();

    #[cfg(feature = "parallel")]
    {
        use ark_ff::Field;
        use ark_ff::Zero;
        parallelize(&mut g_srs, |g, start| {
            let mut current_g: G = G::generator();
            current_g = current_g.mul(tau.pow(&[start as u64]));
            for g in g.iter_mut() {
                *g = current_g;
                current_g *= tau;
            }
        });
    }

    G::normalize_batch(&g_srs)
}

#[cfg(test)]
mod powers_test {
    use ark_ec::{pairing::Pairing, AffineRepr};
    use ark_ff::One;
    use std::ops::Mul;

    fn sanity_srs<E: Pairing>(n: usize, tau: E::ScalarField) -> Vec<E::G1Affine> {
        let powers_of_tau: Vec<E::ScalarField> =
            std::iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(n)
                .collect();

        let g1_gen = E::G1Affine::generator();

        let srs_g1: Vec<E::G1Affine> = powers_of_tau
            .iter()
            .map(|tp| g1_gen.mul(tp).into())
            .collect();

        srs_g1
    }

    // cargo test --features=parallel test_srs
    #[test]
    fn test_srs() {
        use ark_bn254::{Bn254, Fr, G1Projective};
        let k = 1;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let srs_sanity = sanity_srs::<Bn254>(n, tau);
        let g1_srs = super::compute_g_powers::<G1Projective>(tau, n);

        assert_eq!(srs_sanity, g1_srs);
    }
}
