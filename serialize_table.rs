use std::env;
use ark_bn254::Fr;
use mock_cq_setup::{gen_table, serialize_vec, write_bytes};

// cargo run --bin serialize_table {k} {seed} {path}
fn main() {
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

    // there is a table seed
    if args.len() == 4 {
        let table_coeffs = gen_table::<Fr>(k as usize, Some(args[2].as_str()));
        let data = serialize_vec(&table_coeffs.coeffs);
        write_bytes(&args[3], &data);
    } else if args.len() == 3 {
        let table_coeffs = gen_table::<Fr>(k as usize, None);
        let data = serialize_vec(&table_coeffs.coeffs);
        write_bytes(&args[2], &data);
    } else {
        panic!("Not enough arguments")
    };
}
