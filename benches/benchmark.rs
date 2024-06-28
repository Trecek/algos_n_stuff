// benches/benchmark.rs
use algos_n_stuff::algos::bit_packed_ham::BitHamProcessor;
//use algos_n_stuff::algos::common::*;
//use algos_n_stuff::algos::distances::{Distance, HammingDistanceSimd};
//use algos_n_stuff::algos::seq_gen::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;
use std::fs::File;
use std::io::Write;
use serde_json::json;

#[allow(dead_code)]
fn benchmark_bit_ham_flame(c: &mut Criterion) {
    let sequence_length = 10; // Each sequence is 10 bases long
    let num_sequences = 8192;
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
    let processor = BitHamProcessor::new();
    processor.initialize(black_box(&sequences));

    c.bench_function("process_sequences", |b| {
        b.iter(|| {
            processor.process_sequences();
        });
    });
}

fn generate_random_dna_sequence(length: usize) -> Vec<u8> {
    let bases = b"ATCG";
    let mut rng = rand::thread_rng();
    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}


fn benchmark_bit_ham_process_sequences(c: &mut Criterion) {
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
    
    let processor = BitHamProcessor::new();
    processor.initialize(black_box(&sequences));

    let mut bases_per_microsecond = 0.0;

    let mut benchmark = c.benchmark_group("BitHamProcessor");
    benchmark.sample_size(10);
    benchmark.measurement_time(std::time::Duration::from_secs(60));

    // Warm-up phase (not measured)
    for _ in 0..100 {
        black_box(processor.process_sequences());
    }

    benchmark.bench_function("process_sequences", |b| {
        b.iter_custom(|iters| {
            let start = std::time::Instant::now();
            for _ in 0..iters {
                black_box(processor.process_sequences());
            }
            let duration = start.elapsed();

            // Calculate bases per microsecond
            let total_bases_processed = total_bases as f64 * iters as f64;
            let total_microseconds = duration.as_micros() as f64;
            bases_per_microsecond = total_bases_processed / total_microseconds;

            duration
        });
    });

    benchmark.finish();

    // Calculate sequences per microsecond
    let sequences_per_microsecond = bases_per_microsecond / sequence_length as f64;

    // Write benchmark results
    let data = json!({
        "num_sequences": num_sequences,
        "sequence_length": sequence_length,
        "num_pairs": num_pairs,
        "total_bases": total_bases,
        "bases_per_microsecond": bases_per_microsecond,
        "sequences_per_microsecond": sequences_per_microsecond,
        "algorithm_names": ["BitHamProcessor"],
        "processing_speeds": [bases_per_microsecond]
    });

    let mut file = File::create("benches/benchmark_results.json").unwrap();
    file.write_all(data.to_string().as_bytes()).unwrap();
}

criterion_group!(benches, benchmark_bit_ham_process_sequences);
criterion_main!(benches);

