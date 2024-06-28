use safe_arch::*;

use crate::algos::common::*;
use once_cell::sync::OnceCell;
use std::simd::u64x4;
use std::sync::Arc;

const BITS_PER_BASE: usize = 3;
const BASES_PER_U64: usize = 21; // 63 bits for bases, 1 bit unused

/*
We are bit packing 3bit encoded DNA bases
Our words our fed in ascii, when we re-encode, we pack those bits into u64s
A single u64 could contain bits from multiple words
We use SIMD to XOR all the data, this will give us a result of all base2base comparisons
Since our word length is fixed, we can easily extract a words edit distance from the XOR result by index
Our hamming score is the sum of 1s in the XOR result (sum'ed on word range)
*/

pub struct CompactDNA {
    packed_data: Box<[u64x4]>,
    word_length: usize,
}

impl CompactDNA {
    fn new(sequences: &[Vec<u8>]) -> Self {
        let word_length = sequences[0].len();
        let num_u64s = (word_length + BASES_PER_U64 - 1) / BASES_PER_U64;
        let num_u64x4 = (num_u64s + 3) / 4; // Round up to nearest multiple of 4
        let mut packed_data = vec![u64x4::splat(0); sequences.len() * num_u64x4];

        for (i, seq) in sequences.iter().enumerate() {
            let packed_sequence = &mut packed_data[i * num_u64x4..(i + 1) * num_u64x4];
            for (j, &base) in seq.iter().enumerate() {
                let shift = (j % BASES_PER_U64) * BITS_PER_BASE;
                let u64_idx = j / BASES_PER_U64;
                let simd_idx = u64_idx / 4;
                let simd_offset = u64_idx % 4;
                let mut arr = packed_sequence[simd_idx].to_array();
                arr[simd_offset] |= (encode_dna(base) as u64) << shift;
                packed_sequence[simd_idx] = u64x4::from_array(arr);
            }
        }

        CompactDNA {
            packed_data: packed_data.into_boxed_slice(),
            word_length,
        }
    }

    fn calculate_hamming_distance(&self) -> Vec<usize> {
        let num_words = self.packed_data.len();
        let num_pairs = num_words * (num_words - 1) / 2;
        let mut results = vec![0; num_pairs];

        if is_x86_feature_detected!("avx2") {
            self.calculate_hamming_distance_avx2(&mut results);
        } else {
            self.calculate_hamming_distance_fallback(&mut results);
        }

        results
    }

    #[inline(always)]
    fn calculate_hamming_distance_avx2(&self, results: &mut [usize]) {
        let num_words = self.packed_data.len();
        let word_length = self.word_length;
        let used_bits = word_length * BITS_PER_BASE;
        let full_u64_count = used_bits / 64;
        let remaining_bits = used_bits % 64;

        for i in 0..num_words {
            for j in (i + 1)..num_words {
                let mut total_diff = 0;

                // Process full u64s
                for k in 0..full_u64_count {
                    let xor_result = self.packed_data[i].to_array()[k] ^ self.packed_data[j].to_array()[k];
                    total_diff += xor_result.count_ones();
                }

                // Process remaining bits
                if remaining_bits > 0 {
                    let mask = (1u64 << remaining_bits) - 1;
                    let xor_result = self.packed_data[i].to_array()[full_u64_count] ^ self.packed_data[j].to_array()[full_u64_count];
                    let masked_xor = xor_result & mask;
                    total_diff += masked_xor.count_ones();
                }

                let pair_diff = total_diff as usize / 2;

                let pair_index = i * (num_words - 1) - (i * (i + 1) / 2) + j - 1;
                results[pair_index] = pair_diff;
            }
        }
    }
    #[inline(always)]
    fn calculate_hamming_distance_fallback(&self, results: &mut [usize]) {
        let num_words = self.packed_data.len();

        for i in 0..num_words {
            for j in (i + 1)..num_words {
                if j + 2 < num_words {
                    prefetch_t2(&self.packed_data[j + 2]);
                }

                let xor_result = self.packed_data[i] ^ self.packed_data[j];
                let pair_diff = xor_result.to_array().iter()
                    .map(|&v| v.count_ones() as usize) // Use count_ones here
                    .sum::<usize>() / 2;

                let pair_index = i * (num_words - 1) - (i * (i + 1) / 2) + j - 1;
                results[pair_index] = pair_diff;
            }
        }
    }
}

pub struct BitHamProcessor {
    compact_dna: Arc<OnceCell<CompactDNA>>,
}

impl BitHamProcessor {
    pub fn new() -> Self {
        BitHamProcessor {
            compact_dna: Arc::new(OnceCell::new()),
        }
    }

    pub fn initialize(&self, sequences: &[Vec<u8>]) {
        self.compact_dna.get_or_init(|| CompactDNA::new(sequences));
    }

    pub fn process_sequences(&self) -> Vec<usize> {
        let compact_dna = Arc::clone(&self.compact_dna);
        compact_dna.get().unwrap().calculate_hamming_distance()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_ham_small() {
        let sequences = vec![
            b"ATCG".to_vec(),
            b"TAGC".to_vec(),
            b"ATCC".to_vec(),
            b"TACC".to_vec(),
            b"GTCA".to_vec(),
        ];

        println!("Input sequences:");
        for (i, seq) in sequences.iter().enumerate() {
            println!("Sequence {}: {:?}", i, String::from_utf8_lossy(seq));
        }

        let processor = BitHamProcessor::new();
        processor.initialize(&sequences);
        let results = processor.process_sequences();

        println!("Results: {:?}", results);

        assert_eq!(results.len(), 10, "Expected 10 pairwise comparisons");

        let expected_distances = [4, 1, 3, 2, 3, 1, 4, 2, 2, 3];
        println!("Expected distances: {:?}", expected_distances);
        for (i, &distance) in results.iter().enumerate() {
            assert_eq!(
                distance, expected_distances[i],
                "Pair {}: Expected distance = {}, got = {}",
                i, expected_distances[i], distance
            );
        }
    }

    #[test]
    fn test_bit_ham_long() {
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
        processor.initialize(&sequences);
        let results = processor.process_sequences();

        println!("Results: {:?}", results);

        assert_eq!(results.len(), 6, "Expected 6 pairwise comparisons");

        let expected_distances = [19, 21, 17, 21, 21, 19];
        for (i, &distance) in results.iter().enumerate() {
            assert_eq!(
                distance, expected_distances[i],
                "Pair {}: Expected distance = {}, got = {}",
                i, expected_distances[i], distance
            );
        }
    }
}
