#[cfg(feature = "parallel")]
use crate::utils::parallelize;
use ark_ff::FftField;
#[cfg(feature = "parallel")]
use rayon::{self, prelude::*};

/// [tau^0, ..., tau^{n-1}]
pub fn compute_tau_powers<F: FftField>(tau: F, n: usize) -> Vec<F> {
    let mut t_pows = vec![F::zero(); n - 1];

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
            let mut current_tau: F = tau.pow(&[start as u64]);
            for tau_i in tau_chunk.iter_mut() {
                *tau_i = current_tau;
                current_tau *= tau;
            }
        });
    }

    t_pows
}

#[cfg(test)]
mod powers_test {
    use ark_ff::One;

    // cargo test --features=parallel test_tau_pows
    #[test]
    fn test_tau_pows() {
        use ark_bn254::Fr;
        let n = 5;
        let tau = Fr::from(2 as u64);

        let tau_pows = super::compute_tau_powers::<Fr>(tau, n);
        let tau_successors: Vec<Fr> = std::iter::successors(Some(Fr::one()), |p| Some(*p * tau))
            .take(n)
            .collect();

        assert_eq!(tau_pows, tau_successors);
    }
}
