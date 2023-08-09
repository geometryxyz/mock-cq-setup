use ark_poly::DenseUVPolynomial;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain};
use ark_ff::FftField;
use rand::SeedableRng;
use rand::rngs::StdRng;

/// Generates random table from seed
/// If seed is not provided, table is just simple sequence: [1..n]
pub fn gen_table<F: FftField>(k: usize, seed_string: Option<&str>) -> DensePolynomial<F> {
    let n = 1 << k;
    let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

    let t_evals: Vec<F> = match seed_string {
        Some(seed_string) => {
            let seed: [u8; 32] = {
                let mut seed_bytes = [0; 32];
                let string_bytes: &[u8] = seed_string.as_bytes();
                let take = std::cmp::min(32, string_bytes.len());
                seed_bytes[0..take].copy_from_slice(&string_bytes[..take]);
                seed_bytes
            };
            let mut rng = StdRng::from_seed(seed);
            (0..n).map(|_| F::rand(&mut rng)).collect()
        },
        None => {
            (0..n).map(|i| F::from((i + 1) as u64)).collect()
        },
    };

    DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&t_evals))
}