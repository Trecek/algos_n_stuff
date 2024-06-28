// Constants for DNA base encoding (as u8)
const A: u8 = 0b011;
const C: u8 = 0b110;
const G: u8 = 0b101;
const T: u8 = 0b000;

pub fn encode_dna(base: u8) -> u8 {
    match base {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        _ => panic!("Invalid DNA base"),
    }
}

pub fn decode_dna(base: u8) -> u8 {
    match base {
        A => b'A',
        C => b'C',
        G => b'G',
        T => b'T',
        _ => panic!("Invalid DNA base"),
    }
}

/* 
#[inline(always)]
fn popcount_u64x4(v: u64x4) -> u32 {
    //let masked_v = v & u64x4::splat(mask);
    let m256i_vec = m256i::from(v.to_array());
    //println!("Masked v: {:?}", masked_v);
    //println!("m256i_vec: {:?}", m256i_vec);
    let mut total_popcount = 0;

    for i in 0..4 {
        let word = match i {
            0 => safe_arch::extract_i64_from_m256i::<0>(m256i_vec),
            1 => safe_arch::extract_i64_from_m256i::<1>(m256i_vec),
            2 => safe_arch::extract_i64_from_m256i::<2>(m256i_vec),
            _ => safe_arch::extract_i64_from_m256i::<3>(m256i_vec),
        };
        let popcount = safe_arch::population_count_i64(word) as u32;
        total_popcount += popcount;
    }

    total_popcount
}
*/
