use ark_ec::AffineRepr;
use ark_ff::FftField;
use ark_serialize::{CanonicalDeserialize, Read};
use ark_std::log2;
use std::fs::File;

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}

pub fn next_pow2(n: usize) -> usize {
    let two: u32 = 2;
    let a: u32 = log2(n);

    if two.pow(a - 1) == n as u32 {
        return n;
    }

    two.pow(a).try_into().unwrap()
}

#[cfg(feature = "parallel")]
pub fn parallelize<T: Send, F: Fn(&mut [T], usize) + Send + Sync + Clone>(v: &mut [T], f: F) {
    let n = v.len();
    let num_threads = rayon::current_num_threads();
    let mut chunk = (n as usize) / num_threads;
    if chunk < num_threads {
        chunk = n as usize;
    }

    rayon::scope(|scope| {
        for (chunk_num, v) in v.chunks_mut(chunk).enumerate() {
            let f = f.clone();
            scope.spawn(move |_| {
                let start = chunk_num * chunk;
                f(v, start);
            });
        }
    });
}

pub fn write_bytes(f_name: &str, data: &[u8]) {
    use std::io::Write;

    let mut file = File::create(f_name).unwrap();
    file.write_all(&data).unwrap();
}

use ark_serialize::CanonicalSerialize;
pub fn serialize_vec<T: CanonicalSerialize>(x: &[T]) -> Vec<u8> {
    let mut data = Vec::<u8>::new();
    x.serialize_compressed(&mut data).unwrap();
    data
}

pub fn deserialize_vec<T: CanonicalDeserialize>(path: &str) -> Vec<T> {
    let mut file = File::open(path).expect("File could not be opened.");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    Vec::<T>::deserialize_compressed(buffer.as_slice()).unwrap()
}
