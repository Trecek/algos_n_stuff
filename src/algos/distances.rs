
use serde::{Deserialize, Serialize};
use std::simd::cmp::SimdOrd;
use std::simd::num::SimdUint;
use std::simd::*;
use std::ops::BitAnd;
use std::simd::num::SimdInt;

// This is the size of the bit vector, corresponding to ascii characters
pub const PEQ_SIZE: usize = 256;

pub trait Distance<T: ?Sized> {
    fn distance(&self, a: &T, b: &T) -> isize;
    // Input must be a byte slice for SIMD myers variant algorithms
    fn find_distance(&self, t: &[u8], p: &[u8]) -> isize;
}

pub struct SequenceLevenshteinDistanceSimd;

impl SequenceLevenshteinDistanceSimd {
    pub fn new() -> Self {
        SequenceLevenshteinDistanceSimd
    }

    /*
    TODO: use array size 4 for DNA, use something like "ord(nucleotide) % 4" instead of reencoding
    Each element in simd peq is a size 8 element aray of u16 values
    These values corresond the the Ascii value of the character
    The simd peq array is size 256, one for each ascii character
    We have 8 elements for every element in the 256 array
    In my use cases I window sequence data, so this allows me to process 8 windows at a time
        for each bitwise operation
    peq: [
    u16x8 { ... },  // For ASCII character 0
    u16x8 { ... },  // For ASCII character 1
    ...
    u16x8 { ... },  // For ASCII character 255
    ]
    */
    #[inline(always)]
    fn sequence_levenshtein_simd(&self, read: &[u8], barcode: &[u8]) -> Vec<(usize, usize)> {
        // TODO: make SIMD vector width configurable
        // TODO: make simd type configurable (u8, u16, u32, u64, etc.)
        const SIMD_WIDTH: usize = 8; // Assuming 8 windows processed simultaneously
        let mut matches = vec![];

        let windows: Vec<_> = read.windows(SIMD_WIDTH).collect();
        //println!("windows: {:?}", windows);
        let mut min_last_col = i16x8::splat(barcode.len() as i16);
        let mut score = i16x8::splat(barcode.len() as i16);
        let mut peq = [0u16; PEQ_SIZE];

        // Fill bit vectors for each character in barcode
        // Same peq used for all windows (this would need changed if we feed windows from different sequences)
        for i in 0..barcode.len() {
            peq[barcode[i] as usize] |= 1 << i;
        }
        //println!("peq: {:?}", peq);
        let mut pv = u16x8::splat(!0);
        let mut mv = u16x8::splat(0);
        let hb = u16x8::splat(1 << (barcode.len() - 1) as u16);

        for j in 0..barcode.len() {
            let mut eq_values = [0; SIMD_WIDTH];
            for i in 0..SIMD_WIDTH {
                let window_char = windows[i][j] as usize;
                eq_values[i] = peq[window_char];
            }
            let eq = u16x8::from(eq_values);

            // Our data is in a SIMD vector, we don't need to do anything different
            //     for bitwise operations. The compiler will handle it
            let xv = eq | mv;
            let xh = (((eq & pv) + pv) ^ pv) | eq;
            let ph = mv | !(xh | pv);
            let mh = pv & xh;

            let ph_mask: i16x8 = ph.bitand(&hb).cast();
            let mh_mask: i16x8 = mh.bitand(&hb).cast();
            // This is hacky way to deal with bools which are poorly supported in simd
            let ph_add = ph_mask.signum();
            let mh_sub = mh_mask.signum();

            score = score + ph_add;
            score = score - mh_sub;

            let ph = ph << 1 | u16x8::splat(1);
            let mh = mh << 1;
            pv = mh | !(xv | ph);
            mv = ph & xv;

            min_last_col = min_last_col.simd_min(score);
        }

        let min_last_col_array = min_last_col.to_array();

        for i in 0..SIMD_WIDTH {
            if min_last_col_array[i] <= 1 {
                let window_index = i;
                matches.push((window_index, window_index + barcode.len() - 1));
            }
        }

        matches
    }
}

impl<T: AsRef<[u8]> + ?Sized> Distance<T> for SequenceLevenshteinDistanceSimd {
    #[inline(always)]
    fn distance(&self, a: &T, b: &T) -> isize {
        let read = a.as_ref();
        let barcode = b.as_ref();
        let matches = self.sequence_levenshtein_simd(read, barcode);
        if matches.is_empty() {
            read.len() as isize
        } else {
            matches[0].0 as isize
        }
    }

    #[inline(always)]
    fn find_distance(&self, read: &[u8], barcode: &[u8]) -> isize {
        let matches = self.sequence_levenshtein_simd(read, barcode);
        if matches.is_empty() {
            read.len() as isize
        } else {
            matches[0].0 as isize
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceLevenshteinDistance;

impl SequenceLevenshteinDistance {
    pub fn new() -> Self {
        SequenceLevenshteinDistance
    }

    // This is sequence levenshtein distance modified myers algorithm
    #[inline(always)]
    fn sequence_levenshtein(&self, t: &[u8], n: usize, p: &[u8], m: usize) -> isize {
        let mut min_last_col = n as i16;
        let mut score = n as i16;
        let mut peq = [0u16; PEQ_SIZE];
        // Fill bit vector
        for i in 0..n {
            peq[t[i] as usize] |= 1u16 << i;
        }
        let mut pv = !0u16;
        let mut mv = 0u16;
        let hb = 1u16 << (n - 1);
        for j in 0..m {
            let eq = peq[p[j] as usize];
            let xv = eq | mv;
            let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;
            let ph = mv | !(xh | pv);
            let mh = pv & xh;
            if ph & hb != 0 {
                score += 1;
            }
            if mh & hb != 0 {
                score -= 1;
            }
            let ph = ph << 1 | 1;
            let mh = mh << 1;
            pv = mh | !(xv | ph);
            mv = ph & xv;
            if score < min_last_col {
                min_last_col = score;
            }
        }
        min_last_col as isize
    }
}

impl<T: AsRef<[u8]> + ?Sized> Distance<T> for SequenceLevenshteinDistance {
    #[inline(always)]
    fn distance(&self, a: &T, b: &T) -> isize {
        let t = a.as_ref();
        let p = b.as_ref();

        let n = t.len();
        let m = p.len();

        // Instead of calculating twice, this would be easy to adapt with SIMD
        // We can collect the windows of each sequence and process them in parallel
        let score_t = self.sequence_levenshtein(t, n, p, m);
        let score_p = self.sequence_levenshtein(p, m, t, n);
        std::cmp::min(score_t, score_p)
    }

    #[inline(always)]
    fn find_distance(&self, t: &[u8], p: &[u8]) -> isize {
        let n = t.len();
        let m = p.len();

        let score_t = self.sequence_levenshtein(t, n, p, m);
        let score_p = self.sequence_levenshtein(p, m, t, n);
        std::cmp::min(score_t, score_p)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HammingDistance;

impl HammingDistance {
    pub fn new() -> Self {
        HammingDistance
    }
}

impl<T: AsRef<str> + ?Sized> Distance<T> for HammingDistance {
    fn distance(&self, a: &T, b: &T) -> isize {
        let a = a.as_ref();
        let b = b.as_ref();
        a.chars().zip(b.chars()).filter(|(a, b)| a != b).count() as isize
    }

    fn find_distance(&self, a: &[u8], b: &[u8]) -> isize {
        let a = std::str::from_utf8(a).unwrap();
        let b = std::str::from_utf8(b).unwrap();
        a.chars().zip(b.chars()).filter(|(a, b)| a != b).count() as isize
    }
}

// Wagner-Fischer algorithm
// This is what you want if the sub-string you want to identify is large (>64 characters)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceLevenshteinDistanceWagner;

impl SequenceLevenshteinDistanceWagner {
    pub fn new() -> Self {
        SequenceLevenshteinDistanceWagner
    }

    // Its not needless, using clippy suggestion WILL BREAK IT
    #[allow(clippy::needless_range_loop)]
    fn wagner_distance(&self, s1: &[u8], s2: &[u8]) -> isize {
        let len1 = s1.len();
        let len2 = s2.len();
        let mut current_row = vec![0; len2 + 1];
        let mut previous_row: Vec<usize> = vec![0; len2 + 1];

        for j in 1..=len2 {
            current_row[j] = j;
        }

        for i in 1..=len1 {
            std::mem::swap(&mut previous_row, &mut current_row);
            current_row[0] = i;
            // We can't use enumerate in this case
            // b/c we need to mutate current_row inside the loop
            for j in 1..=len2 {
                if s1[i - 1] == s2[j - 1] {
                    current_row[j] = previous_row[j - 1];
                } else {
                    current_row[j] = std::cmp::min(
                        std::cmp::min(previous_row[j] + 1, current_row[j - 1] + 1),
                        previous_row[j - 1] + 1,
                    );
                }
            }
        }

        // Find the minimum value in the last row and the last column
        let min_last_row = *current_row.iter().min().unwrap();
        let min_last_col = previous_row[len2];

        std::cmp::min(min_last_row as isize, min_last_col as isize)
    }
}

impl<T: AsRef<[u8]> + ?Sized> Distance<T> for SequenceLevenshteinDistanceWagner {
    fn distance(&self, a: &T, b: &T) -> isize {
        let s1 = a.as_ref();
        let s2 = b.as_ref();
        self.wagner_distance(s1, s2)
    }

    fn find_distance(&self, s1: &[u8], s2: &[u8]) -> isize {
        self.wagner_distance(s1, s2)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LevenshteinDistance;

// This is from the bk tree repo. I kept around for cross checking.
impl LevenshteinDistance {
    pub fn new() -> Self {
        LevenshteinDistance
    }

    fn levenshtein_distance(&self, a: &str, b: &str) -> isize {
        if a == b {
            return 0;
        }

        let a_len = a.chars().count();
        let b_len = b.chars().count();

        if a_len == 0 {
            return b_len as isize;
        }

        if b_len == 0 {
            return a_len as isize;
        }

        let mut res = 0;
        let mut cache: Vec<usize> = (1..).take(a_len).collect();
        let mut a_dist;
        let mut b_dist;

        for (ib, cb) in b.chars().enumerate() {
            res = ib;
            a_dist = ib;
            for (ia, ca) in a.chars().enumerate() {
                b_dist = if ca == cb { a_dist } else { a_dist + 1 };
                a_dist = cache[ia];

                res = if a_dist > res {
                    if b_dist > res {
                        res + 1
                    } else {
                        b_dist
                    }
                } else if b_dist > a_dist {
                    a_dist + 1
                } else {
                    b_dist
                };

                cache[ia] = res;
            }
        }

        res as isize
    }
}

impl<T: AsRef<str> + ?Sized + Clone> Distance<T> for LevenshteinDistance {
    fn distance(&self, a: &T, b: &T) -> isize {
        let a = a.as_ref();
        let b = b.as_ref();
        self.levenshtein_distance(a, b)
    }

    fn find_distance(&self, a: &[u8], b: &[u8]) -> isize {
        let a = std::str::from_utf8(a).unwrap();
        let b = std::str::from_utf8(b).unwrap();
        self.levenshtein_distance(a, b)
    }
}

#[cfg(test)]
mod tests {
    // TODO: Re-add example from papers as unit tests
    use super::*;

    #[test]
    fn test_sequence_levenshtein_simd() {
        let dist = SequenceLevenshteinDistanceSimd::new();

        // Test case
        let read = b"ACGTACGTGGGGGGG";
        let barcode = b"ACGTACGT";
        let expected_matches = vec![(0, 7), (1, 8)];
        let matches = dist.sequence_levenshtein_simd(read, barcode);
        assert_eq!(matches, expected_matches);
    }

    #[test]
    fn test_sequence_levenshtein_normal() {
        let dist = SequenceLevenshteinDistance::new();

        // Test case
        let window1 = b"ACGTACGT";
        let barcode = b"ACGTACGT";
        let expected_distance = 0;
        let matches = dist.sequence_levenshtein(window1, 8, barcode, 8);
        assert_eq!(matches, expected_distance);
    }
}
