#[cfg(feature = "parallel")]
use crate::utils::parallelize;
#[cfg(feature = "serialize")]
use crate::utils::{serialize_points, write_points};
use ark_ec::CurveGroup;
#[cfg(feature = "parallel")]
use rayon::{self, prelude::*};

pub fn compute_g_powers<G: CurveGroup>(
    tau: G::ScalarField,
    n: usize,
    path: Option<&str>,
) -> Option<Vec<G::Affine>> {
    let mut g_srs = vec![G::zero(); n];

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
        g_srs.push(G::zero());
        parallelize(&mut g_srs, |g, start| {
            let mut current_g: G = G::generator();
            current_g = current_g.mul(tau.pow(&[start as u64]));
            for g in g.iter_mut() {
                *g = current_g;
                current_g *= tau;
            }
        });
    }

    #[cfg(feature = "serialize")]
    {
        let path = path.expect("Path not provided");
        let g_affine = G::normalize_batch(&g_srs);
        let data = serialize_points(&g_affine);
        write_points(path, &data);
        None
    }

    #[cfg(not(feature = "serialize"))]
    Some(G::normalize_batch(&g_srs))
}

#[cfg(not(feature = "serialize"))]
#[cfg(test)]
mod powers_test {
    use ark_ec::{pairing::Pairing, AffineRepr};
    use ark_ff::One;
    use std::ops::Mul;

    fn sanity_srs<E: Pairing>(n: usize, tau: E::ScalarField) -> Vec<E::G1Affine> {
        let powers_of_tau: Vec<E::ScalarField> =
            std::iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(n + 1)
                .collect();

        let g1_gen = E::G1Affine::generator();

        let srs_g1: Vec<E::G1Affine> = powers_of_tau
            .iter()
            .map(|tp| g1_gen.mul(tp).into())
            .collect();

        srs_g1
    }

    // cargo test --all --features=parallel -- --nocapture
    #[test]
    fn test_srs() {
        use ark_bn254::{Bn254, Fr, G1Projective};
        let k = 1;
        let n = 1 << k;
        let tau = Fr::from(100 as u64);

        let srs_sanity = sanity_srs::<Bn254>(n, tau);
        let g1_srs = super::compute_g_powers::<G1Projective>(tau, n, None).unwrap();

        assert_eq!(srs_sanity, g1_srs);
    }
}
