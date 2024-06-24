use std::simd::u64x2;
use rayon::prelude::*;
/*
We are bit packing 3bit encoded DNA bases
Our words our fed in ascii, when we re-encode, we pack those bits into u64s
We are packing bits as they appear regardless of word boundaries
A single u64 could contain bits from multiple words
We use SIMD to XOR all the data, this will give us a result of all base to base comparisons
Our hamming score is the sum of 1s in the XOR result
Since our word length is fixed, we can extract a words edit distance from the XOR result by index
*/
const BITS_PER_BASE: usize = 3;

// Constants for DNA base encoding (as u8)
const A: u8 = 0b011;
const C: u8 = 0b110;
const G: u8 = 0b101;
const T: u8 = 0b000;

fn encode_base(base: u8) -> u8 {
    match base {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        _ => panic!("Invalid DNA base"),
    }
}

const K1: u64 = 0x5555555555555555; // binary: 0101...
const K2: u64 = 0x3333333333333333; // binary: 00110011..
const K4: u64 = 0x0f0f0f0f0f0f0f0f; // binary: 4 zeros, 4 ones ...

// This worked faster than count_ones on my machines so far
fn popcount64(mut x: u64) -> u32 {
    x = x - ((x >> 1) & K1);
    x = (x & K2) + ((x >> 2) & K2);
    x = (x + (x >> 4)) & K4;
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x + (x >> 32);
    (x & 0x7f) as u32
}

pub struct CompactDNA {
    data: Vec<u64>,
    word_length: usize,
}

impl CompactDNA {
    fn new(sequences: &[Vec<u8>]) -> Self {
        // TODO: implement checks for fixed word length
        let word_length = sequences[0].len();
        let sequence_count = sequences.len();
        let bits_per_sequence = word_length * BITS_PER_BASE;
        let u64_per_sequence = (bits_per_sequence + 63) / 64;
        let total_u64s = u64_per_sequence * sequence_count;

        let mut compact_dna = CompactDNA {
            data: vec![0; total_u64s],
            word_length,
            //sequence_count,
        };
        compact_dna.push_sequences(sequences);
        compact_dna
    }

    // Pack base bits into u64s, a words bits can span multiple u64s, we use indexing to keep track
    // interleave word pairs so they are consecutive in memory, increase cache hit rate
    fn push_sequences(&mut self, sequences: &[Vec<u8>]) {
        let mut current_u64_index = 0;
        let mut current_bit_index = 0;

        for sequence in sequences {
            for &base in sequence.iter().take(self.word_length) {
                let encoded = encode_base(base) as u64;
                self.data[current_u64_index] |= encoded << current_bit_index;

                current_bit_index += BITS_PER_BASE;
                if current_bit_index >= 64 {
                    current_u64_index += 1;
                    current_bit_index = 0;
                }
            }
            
            // Move to the next u64 for the next sequence
            if current_bit_index > 0 {
                current_u64_index += 1;
                current_bit_index = 0;
            }
        }
    }

    fn calculate_hamming_distance(&self, index_a: usize, index_b: usize) -> usize {
        let bits_per_sequence = self.word_length * BITS_PER_BASE;
        let u64_per_sequence = (bits_per_sequence + 63) / 64;
        let start_a = index_a * u64_per_sequence;
        let start_b = index_b * u64_per_sequence;

        // Determine the number of full u64x2 operations and any remaining u64
        let full_simd_ops = u64_per_sequence / 2;
        let has_remainder = u64_per_sequence % 2 != 0;

        // Perform SIMD XOR and popcount
        let mut total_diff = 0u32;

        // Process full u64x2 chunks
        for i in 0..full_simd_ops {
            let a = u64x2::from_array([self.data[start_a + 2*i], self.data[start_a + 2*i + 1]]);
            let b = u64x2::from_array([self.data[start_b + 2*i], self.data[start_b + 2*i + 1]]);
            let xor = a ^ b;

            total_diff += popcount64(xor[0]) + popcount64(xor[1]);
        }

        // Handle the remainder u64 if it exists
        if has_remainder {
            let last_index = 2 * full_simd_ops;
            let last_xor = self.data[start_a + last_index] ^ self.data[start_b + last_index];
            total_diff += popcount64(last_xor);
        }

        // Handle the last u64 if it's not fully used
        if bits_per_sequence % 64 != 0 {
            let mask = (1u64 << (bits_per_sequence % 64)) - 1;
            let last_xor = self.data[start_a + u64_per_sequence - 1] ^ self.data[start_b + u64_per_sequence - 1];
            total_diff -= popcount64(last_xor) - popcount64(last_xor & mask);
        }

        total_diff as usize / 2
    }

}

// placeholder struct
pub struct BitHamProcessor {
}

impl BitHamProcessor {
    pub fn new() -> Self {
        BitHamProcessor { }
    }

    pub fn process_sequences(&self, sequences: &[Vec<u8>]) -> Vec<usize> {
        let compact_dna = CompactDNA::new(sequences);

        // Generate pairs to be compared for edit distance calculation
        (0..sequences.len())
            .into_par_iter()
            .flat_map(|i| {
                ((i + 1)..sequences.len()).into_par_iter().map(move |j| (i, j))
            })
            .map(|(i, j)| compact_dna.calculate_hamming_distance(i, j))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_stream_processor() {
        let sequences = vec![
            b"ATCGATCGATCGATCGATCGA".to_vec(),
            b"TAGCTAGCTAGCTAGCTAGCT".to_vec(),
            b"GCATGCATGCATGCATGCATG".to_vec(),
        ];

        println!("Input sequences:");
        for (i, seq) in sequences.iter().enumerate() {
            println!("Sequence {}: {:?}", i, String::from_utf8_lossy(seq));
        }

        let processor = BitHamProcessor::new();

        let results = processor.process_sequences(&sequences);

        println!("Results: {:?}", results);

        assert_eq!(results.len(), 3, "Expected 3 pairwise comparisons");

        for (i, &distance) in results.iter().enumerate() {
            println!("Pair {}: Distance = {}", i, distance);
            assert_eq!(
                distance, 21,
                "Distance should be 21 for completely different sequences"
            );
        }
    }

    #[test]
    fn test_dna_stream_processor_diff() {
        let sequences = vec![
            b"ATCGATCGATCGATCGATCGA".to_vec(),
            b"TAGCTAGCTAGCTAGCTAGGA".to_vec(),
            b"GCATGCATGCATGCATGCATG".to_vec(),
            b"GCCGATTACGTACGTACGTAC".to_vec(),
        ];

        println!("Input sequences:");
        for (i, seq) in sequences.iter().enumerate() {
            println!("Sequence {}: {:?}", i, String::from_utf8_lossy(seq));
        }

        let processor = BitHamProcessor::new();

        let results = processor.process_sequences(&sequences);

        println!("Results: {:?}", results);

        assert_eq!(results.len(), 6, "Expected 6 pairwise comparisons");

        // Expected distances for each pair
        let expected_distances = [19, 21, 17, 21, 21, 19];

        for (i, &expected_distance) in expected_distances.iter().enumerate() {
            assert_eq!(
                results[i], expected_distance,
                "Pair {}: Expected distance = {}, got = {}",
                i, expected_distance, results[i]
            );
        }
    }
}
