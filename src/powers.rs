use ark_ff::FftField;
#[cfg(feature = "parallel")]
use crate::utils::parallelize;
#[cfg(feature = "serialize")]
use crate::utils::{serialize_points, write_points};
#[cfg(feature = "parallel")]
use rayon::{self, prelude::*};


pub fn compute_tau_powers<F: FftField>(tau: F, n: usize, path: Option<&str>) -> Option<Vec<F>> {
    let mut t_pows = vec![F::zero(); n];

    #[cfg(not(feature = "parallel"))]
    let t_pows: Vec<F> = std::iter::once(F::one())
        .chain(t_pows.iter().scan(F::one(), |state, _| {
            *state *= &tau;
            Some(*state)
        }))
        .collect();

    #[cfg(feature = "parallel")]
    {
        use ark_ff::Field;
        use ark_ff::Zero;
        parallelize(&mut t_pows, |tau_chunk, start| {
            let mut current_tau: G = tau.pow(&[start as u64]);
            for tau_i in tau_chunk.iter_mut() {
                *tau_i = current_tau;
                current_tau *= tau;
            }
        });
    }

    #[cfg(feature = "serialize")]
    {
        let path = path.expect("Path not provided");
        let data = serialize_points(&t_pows);
        write_points(path, &data);
        None
    }

    #[cfg(not(feature = "serialize"))]
    Some(t_pows)
}

#[cfg(not(feature = "serialize"))]
#[cfg(test)]
mod powers_test {
    use ark_ff::One;

    #[test]
    fn test_srs() {
        use ark_bn254::Fr;
        let k = 5;
        let tau = Fr::from(2 as u64);

        let tau_pows = super::compute_tau_powers::<Fr>(tau, k, None).unwrap();
        let tau_successors: Vec<Fr> =
        std::iter::successors(Some(Fr::one()), |p| Some(*p * tau))
            .take(k + 1)
            .collect();

        assert_eq!(tau_pows, tau_successors);
    }
}