use ark_ec::AffineRepr;
use ark_serialize::{CanonicalDeserialize, Read};
use std::fs::File;

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

#[cfg(feature = "serialize")]
pub fn write_points(f_name: &str, data: &[u8]) {
    use std::io::Write;

    let mut file = File::create(f_name).unwrap();
    file.write_all(&data).unwrap();
}

#[cfg(feature = "serialize")]
use ark_serialize::CanonicalSerialize;
#[cfg(feature = "serialize")]
pub fn serialize_points<G: AffineRepr>(x: &[G]) -> Vec<u8> {
    let mut data = Vec::<u8>::new();
    x.serialize_compressed(&mut data).unwrap();
    data
}

pub fn deserialize_points<G: AffineRepr>(path: &str) -> Vec<G> {
    let mut file = File::open(path).expect("File could not be opened.");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    Vec::<G>::deserialize_compressed(buffer.as_slice()).unwrap()
}
