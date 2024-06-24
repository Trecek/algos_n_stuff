
/* 
Implement a custom bit vector type to avoid the overhead of Vec<u8>.
Use a more compact representation for neighbors, such as storing only the changed index and new value.
Parallelize the neighbor generation using Rayon or another parallel processing library.
*/
use fxhash::FxHashSet;
use std::simd::num::SimdUint;
use std::simd::Simd;

const A: u8 = 0b011;
const C: u8 = 0b110;
const G: u8 = 0b101;
const T: u8 = 0b000;

pub fn encode_dna(sequence: &str) -> Vec<u8> {
    sequence
        .bytes()
        .map(|b| match b {
            b'A' => A,
            b'C' => C,
            b'G' => G,
            b'T' => T,
            _ => panic!("Invalid DNA base"),
        })
        .collect()
}

// We could probably have something like (conceptually) a read ahead interator
// It could push or queue characters onto a buffer while the heading function just
// continously processes from buffer.
// Indexes could be used to keep track of the original word
// We could have the buffer be pinned memory for gpu accelleration
pub fn hamming_distance_simd(a: &[u8], b: &[u8]) -> usize {
    let mut distance = 0;
    for (chunk_a, chunk_b) in a.chunks(32).zip(b.chunks(32)) {
        let va = Simd::<u8, 32>::from_slice(chunk_a);
        let vb = Simd::<u8, 32>::from_slice(chunk_b);
        let vxor = va ^ vb;
        distance += vxor.cast::<u32>().reduce_sum() as usize;
    }
    distance / 2
}

pub fn neighbors_simd(sequence: &[u8]) -> FxHashSet<Vec<u8>> {
    let mut neighbors = FxHashSet::default();
    // Assuming `u8x32` was a placeholder for `Simd<u8, 32>`
    let masks = [
        Simd::<u8, 32>::splat(0b001),
        Simd::<u8, 32>::splat(0b010),
        Simd::<u8, 32>::splat(0b100),
    ];

    for i in 0..sequence.len() {
        let chunk_start = i - (i % 32);
        let chunk_end = (chunk_start + 32).min(sequence.len());
        // Ensure the slice is the correct length for SIMD operations
        if chunk_end - chunk_start == 32 {
            let chunk = Simd::from_slice(&sequence[chunk_start..chunk_end]);

            for &mask in &masks {
                let mut neighbor = sequence.to_vec();
                let mut modified_chunk = chunk;
                modified_chunk[i % 32] ^= mask[i % 32];
                // Convert the modified SIMD vector to an array
                let modified_array: [u8; 32] = modified_chunk.to_array();
                // Ensure we don't exceed the bounds of the neighbor vector
                let copy_len = modified_array.len().min(neighbor.len() - chunk_start);
                // Copy the modified array back to the slice
                neighbor[chunk_start..chunk_start + copy_len].copy_from_slice(&modified_array[..copy_len]);
                neighbors.insert(neighbor);
            }
        }
    }

    neighbors
}

pub fn decode_dna(encoded: &[u8]) -> String {
    encoded
        .iter()
        .map(|&b| match b {
            A => 'A',
            C => 'C',
            G => 'G',
            T => 'T',
            _ => panic!("Invalid encoded DNA base"),
        })
        .collect()
}
