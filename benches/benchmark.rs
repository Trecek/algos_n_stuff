// benches/benchmark.rs
use algos_n_stuff::algos::distances::{Distance, HammingDistanceSimd};
use algos_n_stuff::algos::seq_gen::*;
use algos_n_stuff::algos::stream_edit_calc::BitHamProcessor;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fxhash::FxHashSet;
use rand::Rng;

fn benchmark_encode_dna(c: &mut Criterion) {
    let sequence = "ATCGATCGATCGATCGATCG".repeat(1000); // 20000 bases

    c.bench_function("encode_dna", |b| {
        b.iter(|| {
            let encoded = encode_dna(black_box(&sequence));
            black_box(encoded)
        })
    });
}

fn benchmark_neighbors_simd(c: &mut Criterion) {
    let sequence = "ATCGATCGATCGATCGATCG".repeat(100); // 2000 bases
    let encoded = encode_dna(&sequence);

    c.bench_function("neighbors_simd", |b| {
        b.iter(|| {
            let result: FxHashSet<Vec<u8>> = neighbors_simd(black_box(&encoded));
            black_box(result)
        })
    });
}

fn benchmark_hamming_distance_simd(c: &mut Criterion) {
    let sequence1 = "ATCGATCGATCGATCGATCG".repeat(100); // 2000 bases
    let sequence2 = "ATCGATCGATCGATCGATCG".repeat(99) + "ATCGATCGATCGATCGATCT"; // 2000 bases with 1 difference at the end

    let hamming = HammingDistanceSimd::new();

    c.bench_function("hamming_distance_simd", |b| {
        b.iter(|| {
            let distance = hamming.distance(
                black_box(sequence1.as_bytes()),
                black_box(sequence2.as_bytes()),
            );
            black_box(distance)
        })
    });
}

fn generate_random_dna_sequence(length: usize) -> Vec<u8> {
    let bases = b"ATCG";
    let mut rng = rand::thread_rng();
    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}

fn benchmark_dna_stream_processor(c: &mut Criterion) {
    let sequence_length = 10; // Each sequence is 10 bases long
    let num_sequences = 1536;
    println!(
        "Generating {} random DNA sequences of length {}",
        num_sequences, sequence_length
    );
    let sequences: Vec<Vec<u8>> = (0..num_sequences)
        .map(|_| generate_random_dna_sequence(sequence_length))
        .collect();

    // Calculate the number of unique pairs
    let num_pairs = (num_sequences * (num_sequences - 1)) / 2;
    println!("Number of pairs: {}", num_pairs);

    // Calculate the total number of bases
    let total_bases = num_pairs * sequence_length;
    println!("Total number of bases: {}", total_bases);

    c.bench_function("dna_stream_processor_full_pipeline", |b| {
        b.iter(|| {
            let processor = BitHamProcessor::new();
            let results = processor.process_sequences(black_box(&sequences));
            //let results = processor.get_results();
            black_box(results)
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10).measurement_time(std::time::Duration::from_secs(15));
    targets = benchmark_dna_stream_processor
}
criterion_main!(benches);
