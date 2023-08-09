use ark_ec::pairing::Pairing;
use ark_ff::FftField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain,
};
use mock_cq_setup::{
    compute_g_powers, compute_qs, compute_tau_powers, deserialize_vec, lagrange_commitments,
    lagrange_openings_commitments_at_zero,
};
use std::env;
use std::time::Instant;

// Given N, runs the setup
fn run<E: Pairing>(n: usize, tau: E::ScalarField, table_path: &str) {
    let now = Instant::now();
    let powers_of_tau = compute_tau_powers(tau, n);
    let elapsed_time = now.elapsed();
    println!(
        "Running compute_tau_powers() took {} seconds.",
        elapsed_time.as_secs()
    );

    let now = Instant::now();
    let _ = lagrange_commitments::<E::G1>(tau, n as u64);
    let elapsed_time = now.elapsed();
    println!(
        "Running lagrange_commitments() took {} seconds.",
        elapsed_time.as_secs()
    );

    let now = Instant::now();
    let _ = lagrange_openings_commitments_at_zero::<E::G1>(tau, n);
    let elapsed_time = now.elapsed();
    println!(
        "Running lagrange_openings_commitments_at_zero() took {} seconds.",
        elapsed_time.as_secs()
    );

    let now = Instant::now();
    let _ = compute_g_powers::<E::G1>(tau, n);
    let elapsed_time = now.elapsed();
    println!(
        "Running compute_g_powers G1() took {} seconds.",
        elapsed_time.as_secs()
    );

    let now = Instant::now();
    let _ = compute_g_powers::<E::G2>(tau, n + 1);
    let elapsed_time = now.elapsed();
    println!(
        "Running compute_g_powers G2() took {} seconds.",
        elapsed_time.as_secs()
    );

    let now = Instant::now();
    let t = read_table::<E::ScalarField>(table_path);
    let domain = GeneralEvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let _ = compute_qs::<E>(&t, &domain, &powers_of_tau);
    let elapsed_time = now.elapsed();
    println!(
        "Running compute_qs took {} seconds.",
        elapsed_time.as_secs()
    );
}

fn read_table<F: FftField>(path: &str) -> DensePolynomial<F> {
    let table_coeffs = deserialize_vec::<F>(path);
    DensePolynomial::from_coefficients_vec(table_coeffs)
}

fn parse_args<F: FftField>() -> (u64, F, String) {
    let args: Vec<String> = env::args().collect();

    let to_u64 = |arg: &String| -> u64 {
        match arg.parse::<u64>() {
            Ok(arg) => arg,
            Err(_) => {
                panic!("Failed to parse argument as u64.");
            }
        }
    };

    let k = to_u64(&args[1]);
    let tau = to_u64(&args[2]);

    let n = 1 << k;
    let tau = F::from(tau);

    (n, tau, args[3].clone())
}

use ark_bn254::{Bn254, Fr};
// cargo run --bin run_setup {k} {tau} -- --features=parallel
fn main() {
    let (n, tau, table_path) = parse_args::<Fr>();
    run::<Bn254>(n as usize, tau, table_path.as_str());
}
